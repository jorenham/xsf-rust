//! Generalized hypergeometric functions pFq of real arguments.
//!
//! The internal functions are generic over the numbers of parameters. They take the lower
//! parameters together with the 1 of `n!` as `b`, so that the terms of the series are
//! `t_n = z^n prod_i (a_i)_n / prod_j (b_j)_n`.

use core::f64::consts::{PI, TAU};
use core::iter::{from_fn, once, successors, zip};
use core::ops::{Add, ControlFlow};

use super::{Scaled, add_exp, add_exp_in, fold_exp, ln_dd, ln_gamma, mul_gamma_ratio};
use crate::utils::{EPSILON_SQ, is_nonpos_int, neg_one_pow, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{
    DoubleDouble as Dd, FRAC_2_PI_DD, LN_PI_DD, LN_SQRT_2PI_DD, frexp_f64, ldexp_f64, two_prod,
    two_sum,
};
use crate::xsf::cephes::hyp2f1::HYP2F1_LOST;
use crate::xsf::cephes::{gamma, gammasgn, lgam, rgamma};
use crate::xsf::{cospi, sinpi, xsf_hypu_digits, xsf_hypu_min_digits};

mod conformal;
mod recurrence;
mod residue;

use residue::Residues;

/// Maximum number of terms of `pfq_series`
const PFQ_MAXITER: u32 = 10_000;

/// Maximum number of terms of `pfq_series` in double-double arithmetic, fewer as a term takes
/// about 20 times as long
const PFQ_DD_MAXITER: u32 = PFQ_MAXITER / 4;

/// Maximum number of terms of `pfq_series_long`
///
/// Of 199 `None` of 3F2 on 20000 points, it gives 13, and twice as many terms would give 4
/// more, in up to twice the time (measured).
const PFQ_LONG_MAXITER: u32 = 2 * PFQ_MAXITER;

/// Largest estimated relative error of an accurate result, after which no other method is tried
const PFQ_TOL: f64 = 1e-13;

/// Digits of U beyond those of `PFQ_TOL` that `chgu` must estimate for the Borel sum of 2F0 to be
/// taken before the other methods, as these are often more accurate (measured)
const CHGU_DIGITS_MARGIN: f64 = 1.0;

/// The smallest estimated relative error above `PFQ_TOL`
const PFQ_INACCURATE: f64 = f64::from_bits(PFQ_TOL.to_bits() + 1);

/// A result and an estimate of its relative error, which is infinite if there is none
#[derive(Clone, Copy, Debug)]
struct Estimate {
    y: f64,
    err: f64,
}

impl Estimate {
    /// `y` with the relative error `err`, where NaN is no estimate
    fn new(y: f64, err: f64) -> Self {
        let err = if err.is_nan() { f64::INFINITY } else { err };
        Self { y, err }
    }

    /// A result that is exact, or of a function that is NaN where it is inaccurate
    const fn exact(y: f64) -> Self {
        Self { y, err: 0.0 }
    }

    /// `y` with the absolute error `err`, where a `y` that isn't finite or a NaN `err` has no
    /// estimate
    fn absolute((y, err): (f64, f64)) -> Self {
        // (an exact `y = 0` as well)
        let err = if !y.is_finite() || err.is_nan() {
            f64::INFINITY
        } else if err == 0.0 {
            0.0
        } else if err <= PFQ_TOL * y.abs() {
            // (accurate if `err <= PFQ_TOL |y|`, as `err / |y|` may round past `PFQ_TOL`)
            (err / y.abs()).min(PFQ_TOL)
        } else {
            (err / y.abs()).max(PFQ_INACCURATE)
        };
        Self::new(y, err)
    }

    /// `y` with the relative error `err` before it is rounded to `f64`, which may overflow, plus
    /// the spacing of a subnormal `y` relative to it
    fn rounded(y: f64, err: f64) -> Self {
        let normal = y.is_normal() || y.is_infinite();
        let spacing = f64::from_bits(1) / y.abs();
        Self::new(y, if normal { err } else { err + spacing })
    }

    /// `y` as `f64`, with its relative error plus the size of `other` relative to it
    fn scaled(y: Scaled, other: Option<Scaled>) -> Self {
        let y = fold_exp(y);
        let err = other.map(fold_exp).map_or(0.0, |v| {
            let size = v.m.abs() * (1.0 + v.err) / y.m.abs();
            ldexp_f64(size, v.e.saturating_sub(y.e))
        });
        Self::rounded(ldexp_f64(y.m, y.e), y.err + err)
    }

    /// `scaled`, or `None` if `y` overflows at only one end of the interval of its error
    fn scaled_certain(y: Scaled) -> Option<Self> {
        let y = fold_exp(y);
        let sum = Self::rounded(ldexp_f64(y.m, y.e), y.err);
        // (fused, as a rounded `1 + err` can hide an end above the threshold; an end is NaN for
        // the infinite `m` of an exponent that `fold_exp` clamps)
        let overflows = |err: f64| !ldexp_f64(y.m.mul_add(err, y.m), y.e).is_finite();
        (overflows(-sum.err) == overflows(sum.err)).then_some(sum)
    }

    /// Whether its estimated error is at most `PFQ_TOL`
    fn is_accurate(self) -> bool {
        self.err <= PFQ_TOL
    }

    /// The result, if it is accurate
    fn accurate(self) -> Option<f64> {
        self.is_accurate().then_some(self.y)
    }

    /// The result, unless its estimated error says that no digit may be correct
    fn trusted(self) -> Option<f64> {
        // (relative to the exact result, which y can miss by `err |y|`: `err / (1 - err)`)
        (self.err < HYP2F1_LOST * (1.0 - self.err)).then_some(self.y)
    }

    /// The trusted result, or NaN if there is none
    fn or_nan(y: Option<Self>) -> f64 {
        y.and_then(Self::trusted).unwrap_or(f64::NAN)
    }
}

/// The choice between the results of methods that are tried in order
trait OrBetter {
    /// This if it is accurate, else the one of it and `other()` with the smaller estimated
    /// error, or this one if they are equal
    fn or_better(self, other: impl FnOnce() -> Self) -> Self;
}

impl OrBetter for Option<Estimate> {
    fn or_better(self, other: impl FnOnce() -> Self) -> Self {
        match (self, self.is_some_and(Estimate::is_accurate)) {
            (Some(_), true) => self,
            (Some(y), false) => other().filter(|v| v.err < y.err).or(self),
            (None, _) => other(),
        }
    }
}

/// Most parameters, with the 1 of `n!`, for which the rounding errors of the methods in `f64`
/// are measured
const PFQ_MEASURED_PARAMS: f64 = 6.0;

/// Largest size of terms that add nothing to a sum, even if there are `PFQ_MAXITER` of them
const PFQ_NEGLIGIBLE: f64 = 1e-150;

/// The sizes of the parameters for which `pfq_products` is exact enough for all `n`
const PFQ_MODERATE: core::ops::RangeInclusive<f64> = 1e-30..=1e30;

/// Most parameters in `a` or in `b` for which `pfq_products` is exact
const PFQ_MAX_UNCHECKED: usize = 4;

/// Below the logarithm of the smallest positive normal number
const PFQ_LN_UNDERFLOW: f64 = -709.0;

/// Largest distance of a parameter to a non-positive integer `-m` at which its rounding error
/// is counted for the whole tail of the series from `n = m` on, as errors of that size are seen
/// up to a distance of 0.01
const PFQ_NEAR_INTEGER: f64 = 0.0625;

/// The `|x|` below which `ln_gamma` is by `lgam`, so that `mul_gammas` shifts `x` to Stirling's
/// formula
const PFQ_GAMMA_MAX: f64 = 13.0;

/// Maximum order of `pfq_weniger`
const PFQ_WENIGER_MAXITER: usize = 40;

/// Number of the last steps of `pfq_weniger` in its estimate where its noise stopped it
/// (measured)
const PFQ_WENIGER_STEPS: usize = 4;

/// Factor by which the error of the `f64` Weniger sum of 2F0 can exceed its estimate (up to 4 in
/// sweeps of the `a_i` and `z` where it has digits)
const HYP2F0_WENIGER_F64_MARGIN: f64 = 4.0;

/// Largest `z` for `pfq_weniger`, above which the series needs less time (measured)
const PFQ_WENIGER_MAX_Z: f64 = -0.6;

/// Smallest `z` for which `pfq_weniger` is tried before `pfq_continuation`, as it is more often
/// accurate there (measured)
const PFQ_WENIGER_FIRST_Z: f64 = -1.5;

/// From this `|z|` on, `hyp2f2` tries the asymptotic expansion before the series, as for
/// `z < 0` the terms of the series then reach about `e^|z| > 10^13`, and the expansion
/// is often accurate and cheaper
const HYP2F2_ASYMPTOTIC_FIRST: f64 = 30.0;

/// From this `-z` on, `hyp1f2` tries the asymptotic expansion before the series, as the terms of
/// the series then reach about `e^(2 sqrt(-z)) > 10^13`, and the expansion is often accurate
/// and cheaper (for `z > 0`, the series is more accurate)
const HYP1F2_ASYMPTOTIC_FIRST: f64 = 225.0;

/// If the terms of the series of pFq for `p <= q` and `z < 0` reach about this power of
/// `e` times the result, `hyppfq` tries the asymptotic expansion before the series, as for
/// `HYP2F2_ASYMPTOTIC_FIRST` and `HYP1F2_ASYMPTOTIC_FIRST`
const PFQ_ASYMPTOTIC_FIRST: f64 = 30.0;

/// The arithmetic in which a series is summed, with the analogs of the `f64` constants
trait SeriesReal: Real + Add<f64, Output = Self> {
    /// Relative rounding error of an operation
    const EPSILON: f64;
    /// Smallest size of a result that is that accurate: it did not underflow
    const MIN_POSITIVE: f64;
    /// The mantissa, of size in `[1/2, 1)`, and the exponent
    fn frexp(self) -> (Self, i32);
    /// `self` times `2^e`
    fn ldexp(self, e: i32) -> Self;
}

impl SeriesReal for f64 {
    const EPSILON: f64 = f64::EPSILON;
    const MIN_POSITIVE: f64 = f64::MIN_POSITIVE;

    fn frexp(self) -> (Self, i32) {
        frexp_f64(self)
    }

    fn ldexp(self, e: i32) -> Self {
        ldexp_f64(self, e)
    }
}

impl SeriesReal for Dd {
    /// Not 2^-105, as its operations are not correctly rounded to 106 bits
    const EPSILON: f64 = pow2(-100);
    /// Below it, its low part is not normal
    const MIN_POSITIVE: f64 = pow2(106 - 1022);

    fn frexp(self) -> (Self, i32) {
        Dd::frexp(self)
    }

    fn ldexp(self, e: i32) -> Self {
        Dd::ldexp(self, e)
    }
}

/// The smallest `n` with `v + n = 0` for a `v` in `params`, or infinity if there is none
fn first_zero(params: &[f64]) -> f64 {
    let zeros = params.iter().filter(|&&v| is_nonpos_int(v));
    zeros.map(|v| -v).fold(f64::INFINITY, f64::min)
}

/// The degree of the polynomial if an `a_i` is a non-positive integer, else infinity, or `None`
/// at a pole, where a `(b_j)_n` vanishes before
fn pfq_degree(a: &[f64], b: &[f64]) -> Option<f64> {
    let degree = first_zero(a);
    (first_zero(b) >= degree).then_some(degree)
}

/// Whether `v` is a pole of the Gamma function
fn is_pole(v: &Dd) -> bool {
    v.hi <= 0.0 && v.is_integer()
}

/// The `n` from which on all `a_i + n` and `b_j + n` are positive
fn pfq_settled(a: &[f64], b: &[f64]) -> f64 {
    a.iter().chain(b).fold(0.0, |n, &v| n.max(-v))
}

/// The factor of the measured rounding errors for more than `PFQ_MEASURED_PARAMS` parameters,
/// as each of them adds its own to a term
fn pfq_rounding_scale(a: &[f64], b: &[f64]) -> f64 {
    ((a.len() + b.len()) as f64 / PFQ_MEASURED_PARAMS).max(1.0)
}

/// `v` in descending order, as `pfq_ratio_bound` is tightest if `a` and `b` are sorted alike
fn descending<P: AsMut<[f64]>>(mut v: P) -> P {
    v.as_mut().sort_unstable_by(|x, y| y.total_cmp(x));
    v
}

/// Bound on `|t_{m+1} / t_m|` for all `m >= n`, if all `a_i + n` and `b_j + n` are positive,
/// except that for a longer `a` the unpaired `a_i + m` are taken at `m = n` (and they grow)
///
/// Each `a_i` is paired with `b_i`: `(a_i + m) / (b_i + m)` decreases to 1, or stays below 1.
/// The other `1 / (b_j + m)` decrease. This is tightest if `a` and `b` are sorted alike.
fn pfq_ratio_bound(a: &[f64], b: &[f64], z: f64, n: f64) -> f64 {
    let paired = a.iter().zip(b).map(|(a, b)| ((a + n) / (b + n)).max(1.0));
    let unpaired_a = a.iter().skip(b.len()).map(|a| a + n);
    let unpaired_b = b.iter().skip(a.len()).map(|b| b + n);
    z.abs() * paired.product::<f64>() * unpaired_a.product::<f64>() / unpaired_b.product::<f64>()
}

/// `pfq_ratio_bound`, but infinite unless the ratios `|t_{m+1} / t_m|` from `m = n` on stay
/// at least 1 once they reach it, for a longer `a`
///
/// They are at most 1 up to the `m0` where `pfq_ratio_bound` is 1 with each unpaired `a_i + n`
/// replaced by `x + m0`, for the largest unpaired `a_i = x`. From `m0` on they grow if the
/// derivative of their log is positive: the `k` unpaired `a_i + m` add at least `k / (x + m)`,
/// and a pair with `a_i > b_i` subtracts `(a_i - b_i) / ((a_i + m) (b_i + m))`. Times `x + m`,
/// that is at most its value at `m0`, with `(x + m) / (a_i + m)` at least 1, as its factors
/// then decrease or stay below 1.
fn pfq_regrowth_bound(a: &[f64], b: &[f64], z: f64, n: f64) -> f64 {
    let bound = pfq_ratio_bound(a, b, z, n);
    let Some(unpaired) = a.get(b.len()..).filter(|v| !v.is_empty()) else {
        return bound;
    };
    let (k, x) = (unpaired.len() as f64, unpaired[0]);
    let paired = zip(a, b).map(|(a, b)| ((a + n) / (b + n)).max(1.0));
    let m0 = (paired.product::<f64>() * z.abs()).powf(-k.recip()) - x;
    let pulls = zip(a, b).filter(|(a, b)| a > b);
    let pull = pulls.map(|(a, b)| (a - b) / (b + m0) * ((x + m0) / (a + m0)).max(1.0));
    if m0 >= n && pull.sum::<f64>() <= k {
        bound
    } else {
        f64::INFINITY
    }
}

/// The products of the `a_i + n` and of the `b_j + n`, whose quotient is `t_{n+1} / (z t_n)`
///
/// They and their quotient are normal for up to `PFQ_MAX_UNCHECKED` parameters with sizes in
/// `PFQ_MODERATE`, as each nonzero `|a_i + n|` is then in `[2^-100, 2^101]`.
fn pfq_products<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> [T; 2] {
    let shifted = |p: T, &v: &f64| p * (T::from_f64(v) + n);
    let one = T::from_f64(1.0);
    [a.iter().fold(one, shifted), b.iter().fold(one, shifted)]
}

/// Whether `v` is normal in the arithmetic of `T`
fn is_normal<T: SeriesReal>(v: T) -> bool {
    (T::MIN_POSITIVE..=f64::MAX).contains(&v.to_f64().abs())
}

/// `t_{n+1} / (z t_n)` for any parameters and `n`, as a mantissa and an `i64` exponent, so that
/// no partial product over- or underflows
fn pfq_frexp_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> (T, i64) {
    let product = |params: &[f64]| {
        let times = |(m, e): (T, i64), &v: &f64| {
            let (vm, ve) = (T::from_f64(v) + n).frexp();
            let (m, me) = (m * vm).frexp();
            (m, e + i64::from(ve) + i64::from(me))
        };
        params.iter().fold((T::from_f64(1.0), 0), times)
    };
    let ((num, ne), (den, de)) = (product(a), product(b));
    (num / den, ne - de)
}

/// `m` times `2^e`, where beyond `i32` it over- or underflows anyway
fn ldexp_i64<T: SeriesReal>(m: T, e: i64) -> T {
    m.ldexp(e.clamp(i32::MIN.into(), i32::MAX.into()) as i32)
}

/// `t_{n+1} / (z t_n)` for any parameters and `n`, which is zero if it underflows, and infinite
/// if it overflows or would be inexact, as it is not normal
#[cold]
fn pfq_checked_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> T {
    let (ratio, e) = pfq_frexp_ratio::<T>(a, b, n);
    let ratio = ldexp_i64(ratio, e);
    let ln_abs = |params: &[f64]| params.iter().map(|v| (v + n).abs().ln()).sum::<f64>();
    if is_normal(ratio) {
        ratio
    } else if ln_abs(a) - ln_abs(b) < PFQ_LN_UNDERFLOW {
        T::from_f64(0.0)
    } else {
        T::from_f64(f64::INFINITY)
    }
}

/// Bound on the size of the term after `term` if it underflows, or if its ratio does
#[cold]
#[inline(never)]
fn pfq_underflow_bound<T: SeriesReal>(term: T, z: f64) -> T {
    T::from_f64(T::MIN_POSITIVE * term.to_f64().abs().max(1.0) * z.abs().max(1.0))
}

/// Whether `pfq_products` is not exact enough for these parameters
fn pfq_extreme(a: &[f64], b: &[f64]) -> bool {
    let many = a.len().max(b.len()) > PFQ_MAX_UNCHECKED;
    many || a.iter().chain(b).any(|v| !PFQ_MODERATE.contains(&v.abs()))
}

/// `t_{n+1} / (z t_n)`, which is `checked` if the parameters are extreme
#[inline]
fn pfq_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, checked: bool) -> T {
    if checked {
        pfq_checked_ratio(a, b, n)
    } else {
        let [num, den]: [T; 2] = pfq_products(a, b, n);
        num / den
    }
}

/// The term `t_{n+1}` after `term = t_n`, with a ratio that is `checked` if the parameters are
/// extreme
#[inline]
fn pfq_next<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, z: f64, term: T, checked: bool) -> T {
    let ratio = pfq_ratio::<T>(a, b, n, checked) * T::from_f64(z);
    if is_normal(ratio) {
        term * ratio
    } else {
        pfq_rescaled_next(a, b, n, z, term, term * ratio)
    }
}

/// `pfq_next` from the mantissas and exponents of its factors, as `t_{n+1} / t_n` can over- or
/// underflow where `t_{n+1}` is normal; or `next`, as computed, if `t_{n+1}` is not normal either
#[cold]
#[inline(never)]
fn pfq_rescaled_next<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, z: f64, term: T, next: T) -> T {
    let (ratio, e) = pfq_frexp_ratio::<T>(a, b, n);
    let ((zm, ze), (tm, te)) = (frexp_f64(z), term.frexp());
    let rescaled = ldexp_i64(tm * (ratio * T::from_f64(zm)), e + i64::from(ze + te));
    if is_normal(rescaled) { rescaled } else { next }
}

/// The terms `t_1, t_2, ...` of the series after `t_0 = 1`, in the arithmetic of `T`
///
/// A term that underflows is a bound on it. The terms end before one that overflows or is NaN,
/// or that is such a bound and not negligible.
fn pfq_terms<T: SeriesReal>(a: &[f64], b: &[f64], z: f64) -> impl Iterator<Item = T> {
    let checked = pfq_extreme(a, b);
    pfq_terms_of(first_zero(a), z, move |n, term| pfq_next(a, b, n, z, term, checked))
}

/// `pfq_terms` with `next(n, t_n) = t_{n+1}`, for the degree of a polynomial, else infinity
fn pfq_terms_of<T: SeriesReal>(
    degree: f64,
    z: f64,
    next: impl Fn(f64, T) -> T,
) -> impl Iterator<Item = T> {
    // the largest term, which is small once the terms are only upper bounds
    let mut limit = f64::MAX;
    let indices = (0_u32..).map(f64::from);
    indices.scan(T::from_f64(1.0), move |term, n| {
        let next = next(n, *term);
        // (the terms after the last one of a polynomial are zero; tested last, as that is faster)
        *term = if next.to_f64().abs() < T::MIN_POSITIVE && n < degree {
            // the term or its ratio underflows, after which the terms can become large again:
            // follow a bound on them
            limit = PFQ_NEGLIGIBLE;
            pfq_underflow_bound(*term, z)
        } else {
            next
        };
        (term.to_f64().abs() <= limit).then_some(*term)
    })
}

/// The series `sum_n t_n` (DLMF 16.2.1) and an estimate of its absolute error
///
/// The estimate is infinite if the series did not converge within `PFQ_MAXITER` terms, or if
/// its terms overflow, or underflow before they become large again. If it converged, but its
/// terms cancel, then it is summed again in double-double arithmetic.
pub(crate) fn pfq_series(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> (f64, f64) {
    let (y, err, _) = pfq_series_tails(a, b, z);
    (y, err)
}

/// The series in double-double with up to `PFQ_LONG_MAXITER` terms, if it is accurate, as near
/// `z = 1` it can need more than `pfq_series`, and an `f64` sum can drop its small terms
fn pfq_series_long(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> Option<(f64, f64)> {
    let (y, err, ..) = pfq_sum::<Dd>(a, b, z, PFQ_LONG_MAXITER);
    let y = (f64::from(y), err);
    Estimate::absolute(y).is_accurate().then_some(y)
}

/// Largest partial sum, in units of `unit`, whose square is added to the sums of squares
const PFQ_SQUARE_MAX: f64 = 1e150;

/// The factor 2^-400 by which the units of those sums shrink once it is exceeded
const PFQ_UNIT: f64 = pow2(-400);

/// Fewest terms that `pfq_flat_series` sums again, as fewer than about 900 roundings of at most
/// 2^-53 each stay below `PFQ_TOL`
const PFQ_FLAT_MIN_TERMS: f64 = 512.0;

/// Replay cutoff of `pfq_flat_series` on `|y - 1|`, which keeps the replay to sums of small
/// terms, much of which `f64` can drop below half an ulp of 1
const PFQ_FLAT_CUTOFF: f64 = pow2(-20);

/// `pfq_series`, and the root of the sum of its squared tails `s - s_n`, which is at least the
/// largest of them
fn pfq_series_tails(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> (f64, f64, f64) {
    let (y, err, tails, terms) = pfq_sum::<f64>(a, b, z, PFQ_MAXITER);
    // (not in double-double, whose `EPSILON` leaves room for many more parameters)
    let err = err * pfq_rounding_scale(a.as_ref(), b.as_ref());
    if terms >= PFQ_FLAT_MIN_TERMS && err <= PFQ_TOL * y.abs() {
        pfq_flat_series(a, b, z, (y, err, tails), terms)
    } else if err <= PFQ_TOL * y.abs() || err == f64::INFINITY {
        (y, err, tails)
    } else {
        let (y, err, tails, _) = pfq_sum::<Dd>(a, b, z, PFQ_DD_MAXITER);
        (f64::from(y), err, tails)
    }
}

/// `pfq_series_tails` of the accurate `f64` sum `series` of `terms` terms after `t_0 = 1`, which
/// for the nonterminating `_{q+1}F_q` with `0 < z < 1` are summed again without it if they are
/// small, and taken plus 1 if that is accurate and `series` is not within `PFQ_TOL` of it
///
/// Added to a sum near 1, `f64` rounds away their parts below half an ulp, which the estimate of
/// `pfq_sum` misses.
/// Each term has at most `pfq_ratio_roundings` per ratio before it, and each addition one, so
/// their sum bounds the error.
#[cold]
#[inline(never)]
fn pfq_flat_series(
    a: &impl AsRef<[f64]>,
    b: &impl AsRef<[f64]>,
    z: f64,
    series: (f64, f64, f64),
    terms: f64,
) -> (f64, f64, f64) {
    let ((a, b), (y, _, tails)) = ((a.as_ref(), b.as_ref()), series);
    let gate = a.len() == b.len() && first_zero(a).is_infinite() && 0.0 < z && z < 1.0;
    // (|y - 1| is at most about the sum of the terms)
    if !gate || (y - 1.0).abs() > PFQ_FLAT_CUTOFF {
        return series;
    }
    // rounded up, as a bound that rounds to below 1 can miss terms that grow
    let ratio_roundings = pfq_ratio_roundings(a, b);
    let rho = pfq_ratio_bound(a, b, z, terms) * (1.0 + ratio_roundings * f64::EPSILON);
    if rho >= 1.0 {
        return series;
    }
    // (a term that underflows is a bound of about `f64::MIN_POSITIVE`, which adds nothing)
    let replay = pfq_terms::<f64>(a, b, z).take(terms as usize);
    let (delta, mass, last) = replay.fold((0.0, 0.0, 0.0), |(d, m, _), t| (d + t, m + t.abs(), t));
    let roundings = terms * (ratio_roundings + 0.5) * f64::EPSILON * mass;
    let err = roundings + last.abs() * rho / (1.0 - rho);
    let flat = 1.0 + delta;
    let flat_err = err + f64::EPSILON / 2.0 * flat;
    if flat_err <= PFQ_TOL * flat && (delta - (y - 1.0)).abs() + err > PFQ_TOL * y {
        (flat, flat_err, tails)
    } else {
        series
    }
}

/// Maximum number of roundings of a ratio `t_{k+1} / t_k`: each `a_i + k` and `b_j + k`, the
/// product, the quotient, and the factors `z` and `t_k`
fn pfq_ratio_roundings(a: &[f64], b: &[f64]) -> f64 {
    2.0 * (a.len() + b.len()) as f64 + 3.0
}

/// `pfq_series_tails` in the arithmetic of `T`, with at most `maxiter` terms, and with the sum
/// in that arithmetic, and the number of terms after `t_0` that it added
///
/// The rounding error of each ratio `t_{n+1} / t_n` affects the terms after it, which sum to
/// `s - s_n` for the partial sum `s_n` and the sum `s`. These errors are added in quadrature,
/// like those of the partial sums. A part of them adds up instead, as `a + n` rounds alike
/// for all `n` between two powers of 2. The constants are measured.
fn pfq_sum<T: SeriesReal>(
    a: &impl AsRef<[f64]>,
    b: &impl AsRef<[f64]>,
    z: f64,
    maxiter: u32,
) -> (T, f64, f64, f64) {
    let terms = pfq_terms::<T>(a.as_ref(), b.as_ref(), z);
    pfq_sum_terms(a, b, z, maxiter, terms, f64::EPSILON)
}

/// `pfq_sum` of the given terms `t_1, t_2, ...`, up to a bound on the rest of them of `tol` times
/// the sum. The terms can be those of a more exact `b`.
#[inline]
fn pfq_sum_terms<T: SeriesReal>(
    a: &impl AsRef<[f64]>,
    b: &impl AsRef<[f64]>,
    z: f64,
    maxiter: u32,
    mut terms: impl Iterator<Item = T>,
    tol: f64,
) -> (T, f64, f64, f64) {
    // generic over the storage of the parameters, as the sum is faster if it is compiled for
    // the numbers of parameters of its caller
    let (a, b) = (a.as_ref(), b.as_ref());
    let settled = pfq_settled(a, b);
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let mut sum = T::from_f64(1.0);
    // the partial sum as `f64`
    let mut s = sum.to_f64();
    // the sums of s_n and s_n^2, in units of `unit` so that the squares don't overflow
    let (mut sum1, mut sum2, mut unit) = (0.0, 0.0_f64, 1.0);
    // the sum of |t_k|, in units of `unit`
    let mut abs1 = 1.0;
    let ratio_roundings = pfq_ratio_roundings(a, b);
    // the rounding error, and the root of the sum of (s - s_n)^2
    let estimates = |n: f64, sum: f64, sum1: f64, sum2: f64, unit: f64, abs1: f64| {
        // accurate enough for an estimate, relative to s so that its parts don't overflow
        let root2 = sum2.sqrt();
        let tails = if sum == 0.0 {
            root2
        } else {
            let sum = sum * unit;
            ((n + 1.0) - 2.0 * (sum1 / sum) + sum2 / sum / sum)
                .max(0.0)
                .sqrt()
                * sum.abs()
        };
        let mut rounding = T::EPSILON * ((4.0 + 0.25 * (n + 1.0).sqrt()) * tails + root2);
        if sum2 > f64::MAX * unit * unit {
            // The measured estimate was only accepted where the squares don't overflow. Beyond,
            // t_k has at most k such roundings that all go the same way, so the sum has at most
            // sum_k k |t_k| <= (n + 1) sum_k |t_k| of them.
            rounding += T::EPSILON * (n + 1.0) * ratio_roundings * abs1;
        }
        (rounding / unit, tails / unit)
    };
    for n in (0..maxiter).map(f64::from) {
        while s.abs() * unit > PFQ_SQUARE_MAX {
            // the squares would overflow, so the sums are rescaled (exactly, by a power of 2)
            (sum1, sum2, unit) = (sum1 * PFQ_UNIT, sum2 * PFQ_UNIT * PFQ_UNIT, unit * PFQ_UNIT);
            abs1 *= PFQ_UNIT;
        }
        sum1 += s * unit;
        sum2 += s * unit * (s * unit);
        if n == degree {
            // the series terminated
            let (rounding, tails) = estimates(n, s, sum1, sum2, unit, abs1);
            return (sum, rounding, tails, n);
        }
        // t_(n+1)
        let Some(term) = terms.next() else {
            return (sum, f64::INFINITY, f64::INFINITY, n);
        };
        let t = term.to_f64();
        abs1 += t.abs() * unit;
        sum = sum + term;
        s = sum.to_f64();
        if n >= settled && t.abs() <= tol * s.abs() {
            let rho = pfq_ratio_bound(a, b, z, n + 1.0);
            // geometric bound on the remaining terms, if the series converges (p <= q + 1)
            let tail = t.abs() * rho / (1.0 - rho);
            if a.len() <= b.len() && rho < 1.0 && tail <= tol * s.abs() {
                let (rounding, tails) = estimates(n, s, sum1, sum2, unit, abs1);
                return (sum, rounding + tail, tails, n + 1.0);
            }
        }
    }
    (sum, f64::INFINITY, f64::INFINITY, f64::from(maxiter))
}

/// Weniger's delta transformation of the series, and an estimate of its absolute error
///
/// This is a factorial Levin-type sequence transformation \[1\], which is used for
/// pFq in HypergeometricFunctions.jl \[2\]. For `z < 0` it converges much faster than the
/// series, and also for `z < -1`, where the series diverges. With the term `t_{n+1}` as the
/// remainder estimate of the partial sum `s_n`, its order `k` is
///
/// ```text
/// delta_k = (sum_{j=0}^k c_j s_j / t_{j+1}) / (sum_{j=0}^k c_j / t_{j+1}),
/// c_j = (-1)^j binom(k, j) (j + 1)_{k-1},
/// ```
///
/// by the three-term recursion of its numerator and denominator. It is applied from the `n`
/// where the terms start to alternate, so that the weights of the `s_j` are positive. The same
/// recursion for `|s_j|` then gives their weighted mean, which scales the rounding errors.
///
/// It is computed in the arithmetic of `T`. The estimate is infinite if the orders up to
/// `PFQ_WENIGER_MAXITER` did not converge, or if a term overflows, or underflows before the terms
/// become large again. It is `None` if the terms only alternate after more than `PFQ_MAXITER`
/// terms, or if a term before that overflows or underflows.
///
/// 1. E. J. Weniger, "Nonlinear sequence transformations for the acceleration of convergence
///    and the summation of divergent series", Comput. Phys. Rep. 10 (1989), 189-371
/// 2. R. M. Slevinsky, "Fast and stable rational approximation of generalized hypergeometric
///    functions", Numer. Algor. 98 (2025), 587-624, <https://doi.org/10.1007/s11075-024-01808-w>
fn pfq_weniger<T: SeriesReal>(
    a: &impl AsRef<[f64]>,
    b: &impl AsRef<[f64]>,
    z: f64,
) -> Option<(f64, f64)> {
    let (a, b) = (a.as_ref(), b.as_ref());
    let settled = pfq_settled(a, b).ceil();
    if settled > f64::from(PFQ_MAXITER) {
        return None;
    }
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    // the transformation starts at this n
    let start = settled.min(degree);
    let checked = pfq_extreme(a, b);
    // the rounding errors of the partial sums and of the terms (the factor is measured)
    let eps = 4.0 * T::EPSILON * pfq_rounding_scale(a, b);
    let one = T::from_f64(1.0);
    let (mut term, mut sum) = (one, one);
    // the largest |s_n| before `start`
    let mut scale = 1.0_f64;
    for n in (0..start as u32).map(f64::from) {
        term = pfq_next(a, b, n, z, term, checked);
        if !is_normal(term) {
            // overflow or underflow
            return None;
        }
        sum = sum + term;
        scale = scale.max(sum.to_f64().abs());
    }
    if start == degree {
        // the series terminated
        return Some((sum.to_f64(), eps * scale));
    }
    let mut rows = [[T::from_f64(0.0); 3]; PFQ_WENIGER_MAXITER];
    let mut delta = sum;
    // the last steps, the newest first
    let mut steps = [f64::INFINITY; PFQ_WENIGER_STEPS];
    for k in 0..PFQ_WENIGER_MAXITER {
        let n = start + k as f64;
        term = pfq_next(a, b, n, z, term, checked);
        if !is_normal(term) {
            // overflow, or underflow, after which the terms only matter if they can increase
            let negligible = term.to_f64().abs() < 1.0
                && a.len() <= b.len()
                && pfq_ratio_bound(a, b, z, n) <= 1.0;
            let noise = eps * scale;
            return Some((sum.to_f64(), if negligible { noise } else { f64::INFINITY }));
        }
        // the numerator, the denominator, and the numerator for |s_j|
        rows[k] = [sum / term, one / term, sum.abs() / term];
        for j in (0..k).rev() {
            // the order of `rows[j]`, which starts at s_j
            let (m, i) = (j as f64, (k - j - 1) as f64);
            let c = if i == 0.0 {
                one
            } else {
                let (s, t) = (m + i, m + 2.0 * i);
                T::from_f64((s + 1.0) * s) / T::from_f64((t + 1.0) * t)
            };
            rows[j] = core::array::from_fn(|l| rows[j + 1][l] - c * rows[j][l]);
        }
        let [num, den, abs] = rows[0];
        let next = num / den;
        let noise = eps * ((abs / den).to_f64().abs() + next.to_f64().abs() + scale);
        // the steps need not be below the rounding of the `f64` result, which only double-double
        // gets under (a NaN noise stays NaN)
        let rounding = f64::EPSILON * next.to_f64().abs();
        let tol = if noise < rounding { rounding } else { noise };
        steps.copy_within(..PFQ_WENIGER_STEPS - 1, 1);
        steps[0] = (next - delta).to_f64().abs();
        delta = next;
        // delta is a weighted mean of the last one and of the one of that order from s_1
        // (NaN for k = 0). For a large |z| the weight of the latter is small, and so are the
        // steps, also while the two still differ, so it has to agree as well.
        let other = (rows[1][0] / rows[1][1] - delta).to_f64().abs();
        if steps[0] <= tol && steps[1] <= tol && other <= tol {
            // The steps also count, as they can be smaller than the error once the noise grows
            // past them (e.g. by 4x per order for z = -16). Then they can also dip far below the
            // error where it changes sign, so the largest of the last few (that it took) counts.
            let counted = if noise < rounding { 2 } else { steps.len() }.min(k + 1);
            let err = noise + steps[..counted].iter().fold(other, |m, &s| m.max(s));
            return Some((delta.to_f64(), err));
        }
        sum = sum + term;
    }
    Some((delta.to_f64(), f64::INFINITY))
}

/// `pfq_weniger`, also with the rounding errors (of relative size `eps`) that the terms `t_n`
/// accumulated up to `settled` added, as `a + n` can round alike
fn weniger_drift(
    (y, err): (f64, f64),
    terms: impl Iterator<Item = f64>,
    settled: f64,
    eps: f64,
) -> Estimate {
    // (each scaled, as their sum can overflow)
    let drift = terms
        .zip(0..=settled as u32)
        .map(|(t, n)| 4.0 * eps * f64::from(n) * t.abs());
    Estimate::absolute((y, err + drift.sum::<f64>()))
}

/// `psi(x) e`: the relative change of `Gamma(x)` by the low part `e` of its argument, to first
/// order
fn gamma_low((x, e): (f64, f64)) -> f64 {
    if e == 0.0 { 0.0 } else { crate::digamma(x) * e }
}

/// `Gamma(u) / Gamma(w)` and an estimate of its relative error, for `u` and `w` as `two_sum`s
///
/// The low part `e` of an argument `x + e` is its rounding error, which matters near a pole.
/// It is applied to first order: `Gamma(x + e) ~ Gamma(x) (1 + psi(x) e)`, where the
/// square of `psi(x) e` estimates the error of that.
///
/// The ratio is zero only if `w` is a pole, and `None` if it underflows or overflows.
pub(crate) fn gamma_ratio(u: (f64, f64), w: (f64, f64)) -> Option<(f64, f64)> {
    // `gamma` and `rgamma` are within 2.3 ulp for 99% of |x| <= 170 (measured)
    let (mut ratio, mut err) = (gamma(u.0) * rgamma(w.0), 3.0 * f64::EPSILON);
    if !ratio.is_normal() {
        if is_nonpos_int(w.0) && w.1 == 0.0 {
            return Some((0.0, 0.0));
        }
        // overflow or underflow
        let (ln_u, ln_w) = (lgam(u.0), lgam(w.0));
        ratio = gammasgn(u.0) * gammasgn(w.0) * (ln_u - ln_w).exp();
        err += f64::EPSILON * (ln_u.abs() + ln_w.abs());
    }
    ratio.is_normal().then(|| {
        let (low_u, low_w) = (gamma_low(u), gamma_low(w));
        let err = err + low_u * low_u + low_w * low_w;
        (ratio * (1.0 + low_u - low_w), err)
    })
}

/// The parameters `a` or `b` of a method for pFq: an array, or a `Vec` if their number
/// is only known at run time
pub(crate) trait Params: Clone + AsRef<[f64]> + AsMut<[f64]> {}

impl<P: Clone + AsRef<[f64]> + AsMut<[f64]>> Params for P {}

/// `v` without `v[i]`
fn skip(v: &impl AsRef<[f64]>, i: usize) -> impl Iterator<Item = f64> + Clone {
    let v = v.as_ref();
    (0..v.len()).filter(move |&k| k != i).map(move |k| v[k])
}

/// The analytic continuation to `z < -1` (DLMF 16.8.8) and an estimate of its absolute error, or
/// `None` if a factor of a coefficient overflows or underflows
///
/// This is for `p = q + 1`, where `b` contains the 1 of `n!`. It is a sum over `j` of series in
/// `1 / z` with the parameters `a_j + 1 - b` and `a_j + 1 - a`, or its limit where `a_i` differ
/// by integers, which is `None` if it underflows or overflows.
pub(crate) fn pfq_continuation<P: Params>(a: &P, b: &P, z: f64) -> Option<(f64, f64)> {
    let one = b.as_ref().iter().position(|&v| v == 1.0)?;
    let mut diffs = zip(a.as_ref(), 0..).flat_map(|(aj, j)| skip(a, j).map(move |ak| ak - aj));
    if diffs.any(|d| d.fract() == 0.0) {
        // the terms of DLMF 16.8.8 are infinite: their limit is the sum of the terms of `pfq_alg`
        let b: Vec<f64> = skip(b, one).chain([1.0]).collect();
        let mut terms = vec![None; a.as_ref().len()];
        pfq_alg(a, &b, z, Alg::Value, &mut terms)?;
        let (y, err) = add_exp(terms.into_iter().flatten())?.to_f64();
        return y.is_normal().then_some((y, err * y.abs()));
    }
    let (mut sum, mut err) = (0.0, 0.0);
    for (j, &aj) in a.as_ref().iter().enumerate() {
        // Gamma(a_k - a_j) / Gamma(a_k) and Gamma(b_k) / Gamma(b_k - a_j) for each k != j
        let ratios = skip(a, j).map(|ak| gamma_ratio(two_sum(ak, -aj), (ak, 0.0)));
        let ratios = ratios.chain(skip(b, one).map(|bk| gamma_ratio((bk, 0.0), two_sum(bk, -aj))));
        // (-z)^(-a_j) is the square of this root, which stays normal for larger a_j
        let root = Some((-z).powf(-0.5 * aj)).filter(|r| r.is_normal())?;
        // the coefficient as a mantissa and a power of 2, as products of its factors can
        // underflow or overflow while it does not; their errors are added in quadrature
        let (mut c, mut exp, mut c_err2) = (1.0, 0, 0.0);
        for factor in [Some((root, f64::EPSILON)); 2].into_iter().chain(ratios) {
            let (factor, factor_err) = factor?;
            let (mantissa, e) = frexp_f64(factor);
            (c, exp) = (c * mantissa, exp + e);
            c_err2 += factor_err * factor_err;
        }
        if c == 0.0 {
            // a b_k - a_j is a non-positive integer
            continue;
        }
        let c = ldexp_f64(c, exp);
        // a_j + 1 - b and a_j + 1 - a, and the relative effect of their rounding errors on the
        // terms of their series, which is large near a non-positive integer -m: there it is on
        // the terms from n = m on, and not on those before them
        let (mut rounding, mut near) = (0.0, 0.0);
        let mut shift = |v: f64| {
            let (d, e1) = two_sum(aj, -v);
            let (x, e2) = two_sum(d, 1.0);
            if e1 + e2 != 0.0 {
                let distance = (x - x.round().min(0.0)).abs();
                rounding += (e1 + e2).abs() / distance;
                if distance < PFQ_NEAR_INTEGER {
                    near += (e1 + e2).abs() / distance;
                }
            }
            x
        };
        let mut shifted = |v: &P| {
            let mut v = v.clone();
            v.as_mut().iter_mut().for_each(|x| *x = shift(*x));
            descending(v)
        };
        let (upper, lower) = (shifted(b), shifted(a));
        let (f, f_err, tails) = pfq_series_tails(&upper, &lower, 1.0 / z);
        sum += c * f;
        // those terms sum to a tail of the series; the last part is for a subnormal c
        err += c.abs() * (f_err + (c_err2.sqrt() + rounding) * f.abs() + near * tails)
            + f64::MIN_POSITIVE * f.abs();
    }
    Some((sum, err))
}

/// A method for `p = q + 1`, which returns the result and an estimate of its absolute error, or
/// `None`
type Method<P> = fn(&P, &P, f64) -> Option<(f64, f64)>;

/// `_{q+1}F_q` for `z < 1`, or for a polynomial, by the first method for `z` that is accurate,
/// else by the most accurate one
///
/// `b` contains the 1 of `n!`. Unless it is a polynomial, the last methods are `pfq_conformal`,
/// which is the slowest, and for `|z| < 1` `pfq_recurrence`, `pfq_series_long` and
/// `pfq_recurrence_lowered`, of which the first and last only apply to 3F2.
fn pfq_accurate<P: Params>(a: &P, b: &P, z: f64) -> Option<Estimate> {
    let series: Method<P> = |a, b, z| Some(pfq_series(a, b, z));
    let polynomial = first_zero(a.as_ref()).is_finite();
    let methods: &[Method<P>] = if polynomial || z > PFQ_WENIGER_MAX_Z {
        &[series]
    } else if z > -1.0 {
        &[pfq_weniger::<f64>, series]
    } else if z == -1.0 {
        &[pfq_weniger::<f64>]
    } else if z >= PFQ_WENIGER_FIRST_Z {
        &[pfq_weniger::<f64>, pfq_continuation]
    } else {
        &[pfq_continuation, pfq_weniger::<f64>]
    };
    let last: Option<Method<P>> = (!polynomial).then_some(conformal::pfq_conformal);
    let inside: &[Method<P>] = if !polynomial && z.abs() < 1.0 {
        &[recurrence::pfq_recurrence, pfq_series_long, recurrence::pfq_recurrence_lowered]
    } else {
        &[]
    };
    let methods = methods.iter().chain(&last).chain(inside);
    methods.fold(None, |y, method| y.or_better(|| method(a, b, z).map(Estimate::absolute)))
}

/// `pfq_series` after sorting the parameters, which tightens `pfq_ratio_bound`
fn pfq_series_sorted(a: impl Params, b: impl Params, z: f64) -> Estimate {
    Estimate::absolute(pfq_series(&descending(a), &descending(b), z))
}

/// The sum of `v` in double-double
fn sum_dd(v: &[f64]) -> Dd {
    v.iter().fold(Dd::ZERO, |s, &v| s + v)
}

/// Whether the series of `_{q+1}F_q` diverges at `z = 1`, where `sum_j b_j <= sum_i a_i`
/// for `b` without the 1 of `n!`, or `None` if that is uncertain
///
/// The sums are in double-double, as their rounding can change the sign of their difference.
fn diverges_at_one(a: &[f64], b: &[f64]) -> Option<bool> {
    let excess = (sum_dd(b) - sum_dd(a)).hi;
    let size = a.iter().chain(b).map(|v| v.abs()).sum::<f64>();
    // certain if no partial sum of the excess rounds
    let mut terms = a.iter().map(|v| -v).chain(b.iter().copied());
    let exact = terms.try_fold(0.0, |s, v| {
        let (s, e) = two_sum(s, v);
        (e == 0.0).then_some(s)
    });
    // a bound on the rounding errors of the additions, each at most `EPSILON_SQ * size / 2`
    let bound = (a.len() + b.len()).max(8) as f64 * EPSILON_SQ * size;
    let uncertain = exact.is_none() && excess.abs() <= bound;
    (!uncertain).then_some(excess <= 0.0)
}

/// The sign of `prod Gamma(a_i) / prod Gamma(b_j)`, and so of the late terms of the series for
/// `z > 0`, as `(v)_n = Gamma(v + n) / Gamma(v)`
fn late_sign(a: &[f64], b: &[f64]) -> f64 {
    a.iter().chain(b).map(|&v| gammasgn(v)).product()
}

/// `_{q+1}F_q` for finite arguments, where the last of `b` is the 1 of `n!`
fn pfq_real<P: Params>(a: P, b: P, z: f64) -> Option<Estimate> {
    if z == 0.0 {
        return Some(Estimate::exact(1.0));
    }
    let Some(degree) = pfq_degree(a.as_ref(), b.as_ref()) else {
        return Some(Estimate::exact(f64::INFINITY));
    };
    if degree.is_infinite() && z > 1.0 {
        // complex
        Some(Estimate::exact(f64::NAN))
    } else if degree.is_infinite() && z == 1.0 {
        // +-inf where the series diverges, and no method where it converges
        let (a, b) = (a.as_ref(), b.as_ref());
        let diverges = diverges_at_one(a, b.split_last()?.1)?;
        diverges.then(|| Estimate::exact(f64::INFINITY * late_sign(a, b)))
    } else {
        pfq_accurate(&descending(a), &descending(b), z)
    }
}

/// Generalized hypergeometric function $_3F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// For $\abs z < 1$,
///
/// $$
/// \hyp 3 2 {a_1\enspace a_2\enspace a_3} {b_1\enspace b_2} z
/// = \sum\_{n=0}^\infty
///     {\rpow{a_1}{n} \ \rpow{a_2}{n} \ \rpow{a_3}{n} \over \rpow{b_1}{n} \ \rpow{b_2}{n}}
///     {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1], and
/// the analytic continuation for $z \le -1$ [^2]. A non-positive integer $a_i = -n$ makes it a
/// polynomial of degree at most $n$, defined for all $z$. An $a_i$ equal to a $b_j$ cancels it,
/// also if it is a non-positive integer.
///
/// It is:
///
/// - $+\infty$ if a $b_j$ is a non-positive integer $-m$ and $z \ne 0$, unless it is a
///   polynomial of degree $\lt m$;
/// - $\pm\infty$ at $z = 1$ for finite parameters if $b_1 + b_2 \le a_1 + a_2 + a_3$, the sign
///   of its limit from below, unless it is a polynomial or that can't be determined (then
///   NaN);
/// - NaN for $z > 1$, where it is complex, unless it is a polynomial or cancels to
///   $(1 - z)^{-a_i}$ with an integer $a_i$ or to a $_2F_1$ that is a rational function of $z$,
///   and if an argument is NaN;
/// - for an infinite argument, the limit as in [`hyppfq`].
///
/// # Returns NaN
///
/// Where its internal error estimate says that none of its digits are correct:
///
/// - for many $1 - 10^{-5} \lesssim z \le 1$, also where an $a_i$ cancels a $b_j$;
/// - for many $z \lesssim -10^{30}$;
/// - if the terms cancel too much (e.g. for large parameters), overflow, or underflow;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.8.E8>
///
/// # See also
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
#[must_use]
pub fn hyp3f2(a1: f64, a2: f64, a3: f64, b1: f64, b2: f64, z: f64) -> f64 {
    Estimate::or_nan(hyp3f2_estimate(a1, a2, a3, b1, b2, z))
}

/// `hyp3f2` and an estimate of its relative error
fn hyp3f2_estimate(a1: f64, a2: f64, a3: f64, b1: f64, b2: f64, z: f64) -> Option<Estimate> {
    if ![a1, a2, a3, b1, b2, z].into_iter().all(f64::is_finite) {
        return pfq_estimate(&[a1, a2, a3], &[b1, b2], z);
    }
    let (a, b) = ([a1, a2, a3], [b1, b2, 1.0]);
    // a common parameter cancels: 2F1, or 1F0 if that has one as well
    let mut pairs = (0..3).flat_map(|i| (0..2).map(move |j| (i, j)));
    let Some((i, j)) = pairs.find(|&(i, j)| a[i] == b[j]) else {
        return pfq_real(a, b, z);
    };
    let (a, b) = ([a[(i + 1) % 3], a[(i + 2) % 3]], [b[1 - j], 1.0]);
    match a.iter().position(|&v| v == b[0]) {
        // the closed form of `pfq_estimate` is also real for z > 1 if `a` is an integer
        Some(i) => pfq_estimate(&[a[1 - i]], &[], z),
        None if z == 1.0 => gauss_sum(a[0], a[1], b[0]).or_better(|| pfq_real(a, b, z)),
        // `pfq_real` is NaN there, as complex, but this 2F1 can be real
        None if z > 1.0 && first_zero(&a).is_infinite() => pfq_estimate(&a, &b[..1], z),
        None => pfq_real(a, b, z),
    }
}

/// Gauss's sum `2F1(a, b; c; 1)` (DLMF 15.4.20), if it converges and is finite
fn gauss_sum(a: f64, b: f64, c: f64) -> Option<Estimate> {
    let c = Dd::new(c);
    let excess = c - a - b;
    if !(excess.hi > 0.0) {
        return None;
    }
    let y = mul_gamma_ratio([c, excess], [c - a, c - b], || Scaled::new(1.0, 0, 0.0))?;
    // its error bound can't tell if it overflows
    Some(Estimate::scaled(y, None)).filter(|sum| sum.y.is_finite())
}

/// A term of a series with a bound on its size
trait Term: Copy {
    /// The arithmetic in which the series is summed
    type Value: SeriesReal;
    fn value(self) -> Self::Value;
    fn bound(self) -> f64;
}

impl<V: SeriesReal> Term for V {
    type Value = Self;

    fn value(self) -> V {
        self
    }

    fn bound(self) -> f64 {
        self.to_f64().abs()
    }
}

/// A value and a bound on its size
impl Term for (f64, f64) {
    type Value = f64;

    fn value(self) -> f64 {
        self.0
    }

    fn bound(self) -> f64 {
        self.1
    }
}

/// The sum of an asymptotic series in `1/z` up to about its smallest term from term `settled`
/// on, and an estimate of its absolute error
///
/// From term `settled` on, all parameters plus `n` must be positive. Before it, the terms may
/// vanish or grow. From term `settled + |z|` on, the term ratios are about `n / |z|`, so that
/// growing terms keep growing. The size of a term `t_n` is the larger of the bounds on `t_n` and
/// on `t_{n+1}`, as a `t_n` of a recurrence of order 2 can vanish while `t_{n+1}` does not.
///
/// The error is the smallest size times `sqrt(2 pi |z|)`, about the sum of the terms near
/// it, plus rounding. For terms of one sign, the remainder can be half of that, as for
/// `sum_n n! / z^(n+1)` with `z > 0`, whose term ratios are also about `n / z`. The sum is the
/// one whose error is smallest, as the rounding grows with `n`. It stops where the remainder is
/// negligible if `ratio_bound(n) < 1`, so that no later term is larger before they grow for good.
///
/// Before term `settled + |z|`, the terms can grow again after a small one, and then shrink to
/// smaller ones, e.g. after a tiny `a_m` or a small positive `1 + a_m - a_i + n` of a term of
/// the algebraic part of pFq. So the error of a sum is taken to be at least its distance to any
/// later sum, less the error of that one (a heuristic, as that error is an estimate too), and
/// infinite if `PFQ_MAXITER` ends the terms while the last is larger than the smallest.
fn asymptotic_sum<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    z: f64,
    ratio_bound: impl Fn(f64) -> f64,
) -> (T::Value, f64) {
    let eps = T::Value::EPSILON;
    let (mut sum, mut abs_sum, mut prev) = (T::Value::from_f64(0.0), 0.0, f64::INFINITY);
    // the smallest size, and the sum before a term with the smallest error
    let (mut min, mut best) = (f64::INFINITY, (sum, f64::INFINITY));
    let remainder = (TAU * z.abs()).sqrt().max(1.0);
    let mut terms = terms.peekable();
    for n in (0..PFQ_MAXITER).map(f64::from) {
        let Some(term) = terms.next() else { break };
        let (t, bound) = (term.value(), term.bound());
        let size = bound.max(terms.peek().map_or(0.0, |next| next.bound()));
        if !t.to_f64().is_finite() || n >= settled + z.abs() && bound >= prev {
            break;
        } else if n >= settled {
            // a later sum can disqualify the current best even when its term is larger
            let err = remainder * size + n * eps * abs_sum;
            best.1 = best.1.max((sum - best.0).to_f64().abs() - err);
            if err < best.1 {
                best = (sum, err);
            }
            if size < min {
                min = size;
                if remainder * min <= eps * sum.to_f64().abs() && ratio_bound(n) < 1.0 {
                    break;
                }
            }
        }
        sum = sum + t;
        abs_sum += bound;
        prev = bound;
        // (at the cap, terms that grew after the smallest may still shrink to smaller ones)
        if n + 1.0 == f64::from(PFQ_MAXITER) && size > min {
            best.1 = f64::INFINITY;
        }
    }
    best
}

/// The sum of a convergent series, in the arithmetic of its terms, and an estimate of its
/// absolute error, which is infinite if it did not converge within `PFQ_MAXITER` terms
///
/// From term `settled` on, the bound on each term is at most `ratio_bound(n)` times the one
/// before it, and `ratio_bound` does not increase.
fn convergent_sum<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    ratio_bound: impl Fn(f64) -> f64,
) -> (T::Value, f64) {
    let eps = T::Value::EPSILON;
    // the partial sums, until the tail is negligible
    let sums = terms
        .zip((0..PFQ_MAXITER).map(f64::from))
        .take_while(|(term, _)| term.value().to_f64().is_finite())
        .try_fold((T::Value::from_f64(0.0), 0.0), |(sum, abs_sum), (term, n)| {
            let (sum, abs_sum) = (sum + term.value(), abs_sum + term.bound());
            let rho = if n >= settled { ratio_bound(n) } else { 1.0 };
            // geometric bound on the remaining terms
            let tail = term.bound() * rho / (1.0 - rho);
            if rho < 1.0 && tail <= eps * sum.to_f64().abs() {
                ControlFlow::Break((sum, tail + n * eps * abs_sum))
            } else {
                ControlFlow::Continue((sum, abs_sum))
            }
        });
    match sums {
        ControlFlow::Break(sum_err) => sum_err,
        ControlFlow::Continue((sum, _)) => (sum, f64::INFINITY),
    }
}

/// The sum of the series in `1 / z` of a term of the algebraic part, where the term ratios are
/// about `n^kappa / |z|`: by `asymptotic_sum`, or by `convergent_sum` for `kappa = 0`, which
/// also stands for `p > q + 1`
fn pfq_alg_sum<T: Term<Value = f64>>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    kappa: usize,
    z: f64,
    ratio_bound: impl Fn(f64) -> f64,
) -> (f64, f64) {
    let az = z.abs();
    let scale = match kappa {
        0 => return convergent_sum(terms, settled, ratio_bound),
        1 => az,
        2 => az.sqrt(),
        _ => az.powf(1.0 / kappa as f64),
    };
    asymptotic_sum(terms, settled, scale, ratio_bound)
}

/// The terms `t_0 = 1, t_1, ...` of a recurrence of order 2, where `next(k, t1, t2)` is `t_k`
/// for `t_{k-1}` and `t_{k-2}`, with `t_{-1} = 0`
fn recurrence_terms(
    next: impl Fn(f64, f64, f64) -> f64 + Clone,
) -> impl Iterator<Item = f64> + Clone {
    // k, t_k, and t_(k-1)
    let terms = successors(Some((0.0, 1.0, 0.0)), move |&(k, t1, t2)| {
        Some((k + 1.0, next(k + 1.0, t1, t2), t1))
    });
    terms.map(|(_, t, _)| t)
}

/// The coefficients `B_2k / (2k (2k - 1))` of Stirling's series (DLMF 5.11.1) in double-double,
/// for `k = 16, ..., 1`
#[rustfmt::skip]
const STIRLING_DD: [Dd; 16] = [
    Dd { hi: -15_238_221.539_407_415, lo: -8.767_745_224_906_25e-10 },
    Dd { hi: 691_472.268_851_313, lo: 2.558_529_630_515_8e-11 },
    Dd { hi: -36_108.771_253_724_99, lo: 5.897_583_353_514_365e-13 },
    Dd { hi: 2_193.103_333_333_333_5, lo: -1.333_925_562_600_294_8e-13 },
    Dd { hi: -156.848_284_626_002_03, lo: 9.391_823_141_715_389e-15 },
    Dd { hi: 13.402_864_044_168_393, lo: -6.154_114_101_993_966e-16 },
    Dd { hi: -1.392_432_216_905_901_1, lo: 1.583_705_698_923_030_3e-17 },
    Dd { hi: 0.179_644_372_368_830_57, lo: -6.401_600_482_710_946e-19 },
    Dd { hi: -0.029_550_653_594_771_242, lo: 4.861_760_957_508_855e-19 },
    Dd { hi: 0.006_410_256_410_256_41, lo: 2.224_004_456_380_521_7e-19 },
    Dd { hi: -0.001_917_526_917_526_917_6, lo: 1.067_570_277_687_247_5e-19 },
    Dd { hi: 0.000_841_750_841_750_841_7, lo: 3.687_017_488_923_769_4e-20 },
    Dd { hi: -0.000_595_238_095_238_095_3, lo: 5.369_382_187_547_26e-20 },
    Dd { hi: 0.000_793_650_793_650_793_7, lo: 6.883_823_317_368_282e-22 },
    Dd { hi: -0.002_777_777_777_777_778, lo: 1.060_108_790_874_715_4e-19 },
    Dd { hi: 0.083_333_333_333_333_33, lo: 4.625_929_269_271_485e-18 },
];

/// The term `k = 17` of Stirling's series at `PFQ_GAMMA_MAX`, which bounds the sum of the terms
/// after those of `STIRLING_DD` from there on
const STIRLING_DD_REST: f64 = 6.7e-29;

/// The `g < 0` down to which `mul_gammas_in` shifts `g` to `ln_gamma_stirling`, in about `|g|`
/// steps
const PFQ_GAMMA_DD_MIN: f64 = -1024.0;

/// `ln Gamma(g)` for `g >= PFQ_GAMMA_MAX` by Stirling's series (DLMF 5.11.1) in double-double,
/// its sign, and a bound on its absolute error
fn ln_gamma_stirling(g: Dd) -> (Dd, f64, f64) {
    let w = Dd::ONE / (g * g);
    let corr = STIRLING_DD.iter().fold(Dd::ZERO, |s, &c| s * w + c) / g;
    let v = (g - 0.5) * ln_dd(g) - g + LN_SQRT_2PI_DD + corr;
    (v, 1.0, STIRLING_DD_REST + Dd::EPSILON * v.hi.abs())
}

/// `y` times `v`, with the low part of the mantissa in the logarithm
fn mul_dd(mut y: Scaled, v: Dd) -> Scaled {
    let p = v * y.m;
    y.m = p.hi;
    y.mul_exp(Dd::new(p.lo / p.hi))
}

/// `y` times `Gamma(num_0) ... / (Gamma(den_0) ...)`, with the error bounds of the factors
/// added, or `None` at a pole (or if the result is far out of range)
///
/// Each factor is by `ln_gamma`, which for `|g| < PFQ_GAMMA_MAX` is at `g + n >= PFQ_GAMMA_MAX`
/// instead, divided by the rising factorial `(g)_n` (DLMF 5.2.5). There it is Stirling's formula
/// in double-double, whose error is smaller than that of `lgam` at `g`. The `(g)_n` are multiplied
/// in double-double, and rounded once.
fn mul_gammas(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    y: Scaled,
) -> Option<Scaled> {
    mul_gammas_in(num, den, y, false)
}

/// `mul_gammas` if not `dd`. Otherwise each factor is by `ln_gamma_stirling` at `g + n`, also for
/// `PFQ_GAMMA_DD_MIN < g <= -PFQ_GAMMA_MAX`, where `ln_gamma` reflects `g` in `f64`, and the
/// `(g)_n` are not rounded, so that the factors are as exact as double-double arithmetic.
fn mul_gammas_in(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    mut y: Scaled,
    dd: bool,
) -> Option<Scaled> {
    // the product of the `(g)_n^-power`, as a mantissa and an exponent
    // (`rising_e` is an `f64`, which sums exactly, as many tiny factors can take it out of the
    // range of `i32`)
    let (mut rising, mut rising_e) = (Dd::ONE, 0.0);
    let gammas = num.into_iter().map(|g| (g, 1.0));
    for (mut g, power) in gammas.chain(den.into_iter().map(|g| (g, -1.0))) {
        let min = if dd { PFQ_GAMMA_DD_MIN } else { -PFQ_GAMMA_MAX };
        if min < g.hi && g.hi < PFQ_GAMMA_MAX {
            // (not the ceiling, which can give `g + n < PFQ_GAMMA_MAX` as the difference rounds)
            let n = (PFQ_GAMMA_MAX - g.hi).floor() + 1.0;
            let ks = (0_u32..).map(f64::from).take_while(|&k| k < n);
            // `(g)_n` as a product of mantissas and an exponent, so that it is normal also for
            // a tiny factor (the product is at least `2^-n`)
            let (r, r_e) = ks.fold((Dd::ONE, 0), |(p, e), k| {
                let (f, j) = (g + k).frexp();
                (p * f, e + j)
            });
            // (else a factor is 0: `g` is a pole, where `ln_gamma` is `None`)
            if r != 0.0 {
                let (q, s) = if power > 0.0 {
                    (rising / r, -r_e)
                } else {
                    (rising * r, r_e)
                };
                let (m, e) = q.frexp();
                (rising, rising_e, g) = (m, rising_e + f64::from(e + s), g + n);
                if dd {
                    // (the rounding of each factor of `(g)_n`)
                    y.err += n * Dd::EPSILON;
                }
            }
        }
        let (v, sign, v_err) = if dd && g.hi >= PFQ_GAMMA_MAX {
            ln_gamma_stirling(g)
        } else {
            ln_gamma(g)?
        };
        (y.m, y.err) = (sign * y.m, y.err + v_err);
        y = y.mul_exp(v * power);
    }
    if rising != 1.0 {
        // (one eps bounds the rounding of `rising` to `rising.hi` and of the product of the
        // mantissas, which can't overflow or underflow; `None` for an exponent beyond the range
        // that `fold_exp` takes, which needs about a million tiny Gamma arguments)
        let ((m, e), (r, k)) = (frexp_f64(y.m), frexp_f64(rising.hi));
        let e = f64::from(y.e + e + k) + rising_e;
        y.e = (e.abs() <= 1e9).then_some(e as i32)?;
        if dd {
            (y.m, y.err) = (m, y.err + Dd::EPSILON);
            y = mul_dd(y, rising.ldexp(-k));
        } else {
            (y.m, y.err) = (m * r, y.err + f64::EPSILON);
        }
    }
    Some(y)
}

/// `exp(l) sum` times `Gamma(num_0) ... / (Gamma(den_0) ...)`, where `err` is the absolute
/// error of `sum`
fn gamma_series(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    (sum, err): (f64, f64),
    l: Dd,
) -> Option<Scaled> {
    mul_gammas(num, den, Scaled::new(sum, 0, err / sum.abs())?.with_l(l))
}

/// The exponentially large part of the asymptotic expansion of 2F2 (DLMF 16.11.1, 16.11.7),
/// with `|z|^nu` instead of `z^nu`:
///
/// ```text
/// Gamma(b_1) Gamma(b_2) / (Gamma(a_1) Gamma(a_2)) e^z |z|^nu sum_k c_k z^-k,
/// nu = a_1 + a_2 - b_1 - b_2
/// ```
///
/// By the differential equation (DLMF 16.8.3), `c_0 = 1` and
/// `k c_k = A_1 c_{k-1} + A_0 c_{k-2}`, with `s = nu + 1 - k` and
///
/// - `A_1 = 2 s^2 + (2 b_1 + 2 b_2 - a_1 - a_2 - 1) s + b_1 b_2 - a_1 a_2`,
/// - `A_0 = (s + 1) (s + b_1) (s + b_2)`.
fn hyp2f2_exp(a: [f64; 2], b: [f64; 2], z: f64) -> Option<Scaled> {
    let nu = Dd::new(a[0]) + a[1] - b[0] - b[1];
    let (sum_a, sum_b, dprod) = (a[0] + a[1], b[0] + b[1], b[0] * b[1] - a[0] * a[1]);
    // the terms t_k = c_k z^-k
    let terms = recurrence_terms(|k, t1, t2| {
        let s = nu.hi + 1.0 - k;
        let a1 = (2.0 * s + 2.0 * sum_b - sum_a - 1.0) * s + dprod;
        let a0 = (s + 1.0) * (s + b[0]) * (s + b[1]);
        (a1 * t1 + a0 / z * t2) / (k * z)
    });
    // (the factors of A_0 are negative from here on)
    let settled = (nu.hi + 2.0 + b[0].max(b[1]).max(0.0)).max(0.0);
    // (without small divisors, the terms don't grow again before `settled + |z|`)
    let series = asymptotic_sum(terms, settled, z, |_| 0.0);
    let l = nu * ln_dd(z.abs()) + z;
    gamma_series(b.map(Dd::new), a.map(Dd::new), series, l)
}

/// `(x, k, factor)` for the rounded parameter `x = 1 + am - v` of a series, where `-k` is the
/// non-positive integer nearest to it: times `factor`, the terms after `t_k` have the factor
/// `x + k` of the exact parameter if it is an `upper` one, and `1 / (x + k)` if it is a lower one
///
/// The rounding error of `x` is large relative to a small `x + k`. An inexact `x` at `-k` would
/// end the series if upper, or divide by zero if lower, so that it is a few ulp from `-k` instead.
fn shifted(am: f64, v: f64, upper: bool) -> (f64, f64, f64) {
    let (s, e1) = two_sum(1.0, am);
    let (nearest, e2) = two_sum(s, -v);
    let (k, e) = (-nearest.round().min(0.0), e1 + e2);
    let offset = nearest + k;
    let lands = offset == 0.0 && e != 0.0;
    // (more than an ulp of -k, also for k = 0)
    let x = nearest - if lands { (k + 1.0) * f64::EPSILON } else { 0.0 };
    let (rounded, exact) = (x + k, offset + e);
    // (an upper parameter at -k is exact here: the series ends)
    let upper_factor = if rounded == 0.0 { 1.0 } else { exact / rounded };
    (x, k, if upper { upper_factor } else { rounded / exact })
}

/// The `k` and factors of the parameters of a series from `shifted`
struct Shifts<P> {
    ks: P,
    factors: P,
}

impl<P: Params> Shifts<P> {
    fn iter(&self) -> impl Iterator<Item = (&f64, &f64)> {
        zip(self.ks.as_ref(), self.factors.as_ref())
    }
}

/// `shifted(am, v, upper)` for each `v` of `params`, but with `value` as the parameter at `own`
fn shifted_params<P: Params>(
    params: &P,
    am: f64,
    upper: bool,
    (own, value): (usize, f64),
) -> (P, Shifts<P>) {
    let (mut xs, mut ks, mut factors) = (params.clone(), params.clone(), params.clone());
    for (i, &v) in params.as_ref().iter().enumerate() {
        (xs.as_mut()[i], ks.as_mut()[i], factors.as_mut()[i]) = shifted(am, v, upper);
    }
    (xs.as_mut()[own], ks.as_mut()[own], factors.as_mut()[own]) = (value, 0.0, 1.0);
    (xs, Shifts { ks, factors })
}

/// The series of the term of `a_m` of the algebraic part (DLMF 16.11.2), as its sum and
/// absolute error, where `b` contains the 1 of `n!` as its last parameter
///
/// With `sizes`, it is the sum of the sizes of its terms.
fn pfq_alg_series<A: Params, B: Params>(
    a: &A,
    b: &B,
    m: usize,
    z: f64,
    sizes: bool,
) -> Option<(f64, f64)> {
    let am = a.as_ref()[m];
    let last = b.as_ref().len().checked_sub(1)?;
    // a_m instead of the 1 of n!, and 1 + a_m - b_j
    let (mut upper, up_shifts) = shifted_params(b, am, true, (last, am));
    // the 1 of n!, and 1 + a_m - a_i
    let (mut lower, low_shifts) = shifted_params(a, am, false, (m, 1.0));
    let (up, low) = (&*descending(upper.as_mut()), &*descending(lower.as_mut()));
    let bound = |n| pfq_regrowth_bound(up, low, 1.0 / z, n);
    // kappa = q - p + 1, or 0 for p > q + 1, where the series converges
    let (p, q) = (a.as_ref().len(), last);
    let kappa = (q + 1).saturating_sub(p);
    // their k and factors, and the largest k with the factor of the terms after it
    let shifts = || up_shifts.iter().chain(low_shifts.iter());
    let (k_max, all) = shifts().fold((0.0_f64, 1.0), |(top, all), (&k, &f)| (k.max(top), all * f));
    let series = once(1.0)
        // in (-1)^(p+q+1) / z
        .chain(pfq_terms(up, low, neg_one_pow(p + q + 1) / z))
        .zip((0_u32..).map(f64::from))
        .map(|(t, n)| {
            // (one product is faster)
            if n > k_max {
                t * all
            } else {
                shifts().fold(t, |t, (&k, &f)| if n > k { t * f } else { t })
            }
        })
        .map(|t| if sizes { t.abs() } else { t });
    Some(pfq_alg_sum(series, pfq_settled(up, low), kappa, z, bound))
}

/// `pfq_alg_series` for `p = q` or `p > q + 1` in double-double, with the exact parameters, and
/// an infinite error if no term bounds the rest
///
/// A product of the factors of a term ratio that isn't normal is less exact, and ends the terms.
fn pfq_alg_series_dd(a: &[f64], b: &[f64], m: usize, z: f64) -> (Dd, f64) {
    let [mut up, mut low] = residue::series_params(a, &b[..b.len() - 1], m);
    // (so that the factors of a term ratio pair by size, with the smallest of `up` unpaired)
    for v in [&mut up, &mut low] {
        v.sort_unstable_by(|x, y| y.hi.total_cmp(&x.hi));
    }
    let poles = up.iter().filter(|x| is_pole(x));
    let degree = poles.fold(f64::INFINITY, |n, x| n.min(-x.hi));
    let ratio = |n: f64| {
        let normal = |p: Dd| if is_normal(p) { p } else { Dd::NAN };
        let prod = |v: &[Dd]| v.iter().fold(Dd::ONE, |p, &x| normal(p * (x + n)));
        // in (-1)^(p+q+1) / z
        if n == degree {
            Dd::ZERO
        } else {
            prod(&up) / (prod(&low) * (neg_one_pow(a.len() + b.len()) * z))
        }
    };
    let terms = once(Dd::ONE).chain(pfq_terms_of(degree, 1.0, |n, t| t * ratio(n)));
    if degree.is_finite() {
        // (all of its terms)
        return asymptotic_sum(terms, degree + 1.0, z, |_| 0.0);
    } else if up.len() < low.len() {
        // (it converges for p > q + 1)
        //
        // The bound is on the factors that are shifted in double-double, as `x.hi + n` is not
        // the rounded `x + n` where that nearly vanishes, from where all of them are positive.
        let shifted = |v: &[Dd], n: f64| v.iter().map(|&x| (x + n).hi).collect::<Vec<_>>();
        let bound = |n| pfq_ratio_bound(&shifted(&up, n), &shifted(&low, n), 1.0 / z, 0.0);
        let positive_from = |x: &Dd| (-*x).floor().hi + 1.0;
        let settled = up.iter().chain(&low).map(positive_from).fold(0.0, f64::max);
        return convergent_sum(terms, settled, bound);
    }
    // Where its factors are positive, each paired one of a term ratio is monotone with the limit
    // 1, and the unpaired one grows. So for `m >= n` the ratio is at most a `w(m)` that is linear
    // in `m`, with `w(n) = above`, and at least `w(m) / spread`.
    //
    // The product of `w(n)` to `w(m)` is at most 1 while `w(n) + w(m) <= 2`: no term up to there
    // is larger than term `n`. And from `w(m) >= spread` on, the terms grow for good. The sum
    // only ends at an `n` with `above + spread <= 2` (above their rounding), and so not at a
    // term that a later one exceeds.
    let ends = |n: f64| {
        let shift = |x: &Dd| (*x + n).hi;
        let unpaired = shift(&up[low.len()]) / z.abs();
        let quotients = zip(&up, &low).map(|(u, l)| shift(u) / shift(l));
        let [above, spread] = quotients.fold([unpaired, 1.0], |[above, spread], f| {
            [above * f.max(1.0), spread * f.max(f.recip())]
        });
        let positive = up.iter().chain(&low).all(|x| shift(x) > 0.0);
        positive && unpaired.is_normal() && (above + spread) * (1.0 + f64::EPSILON.sqrt()) <= 2.0
    };
    // (the unpaired factor alone is above 1 after that)
    let ns = (0..PFQ_MAXITER).map(f64::from);
    let mut ns = ns.take_while(|&n| (up[low.len()] + n).hi <= z.abs());
    let Some(start) = ns.find(|&n| ends(n)) else {
        return (Dd::ONE, f64::INFINITY);
    };
    let terms = zip(terms, (0_u32..).map(f64::from))
        .take_while(|&(_, n)| n < start || ends(n))
        .map(|(term, _)| term);
    asymptotic_sum(terms, start, z, |_| 0.0)
}

/// The term of `pfq_alg` of the class of `a_m`: the `a_i` that differ from `a_m` by integers
///
/// The terms of these `a_i` are infinite, unlike the limit of their sum, which is the sum of
/// the `Residues` at the poles that they share. It is the `term` of the first of the smallest
/// of them, with the residue `(-1)^k / k!` of the Gamma function in place of each `Gamma(-k)`
/// of its factor. The `term` stays `None` for the other `a_i`, and if the class has no residues.
///
/// The residues are complex for `z > 0`, where it is the sum of the bounds on their sizes, as
/// with `sizes`.
///
/// For a polynomial, it is `None` for the class of its non-positive integer `a_i`, which the
/// `Residues` exclude. The terms of the other classes vanish, with the `1 / Gamma(a_i)` of their
/// factor.
fn pfq_class_term(
    a: &[f64],
    b: &[f64],
    m: usize,
    z: f64,
    sizes: bool,
    term: &mut Option<Scaled>,
) -> Option<()> {
    let (am, lower) = (a[m], &b[..b.len() - 1]);
    let others = skip(&a, m).map(Dd::new);
    let diffs = others.clone().map(|ai| ai - am);
    let before = |(d, i): (Dd, usize)| d.is_integer() && (d.hi < 0.0 || d.hi == 0.0 && i < m);
    if diffs.clone().zip(0..).any(before) {
        return Some(());
    }
    if first_zero(a).is_finite() {
        return (!is_nonpos_int(am)).then_some(());
    }
    let residues = Residues::of_class(a, lower, m, z);
    let [up, low] = residues.params().map(descending);
    // (their order is final after the last zero factor, and their bounds are valid once all
    // parameters plus n exceed 1)
    let settled = pfq_settled(&up, &low).floor() + 2.0;
    let bound = |n| pfq_ratio_bound(&up, &low, 1.0 / z, n) * Residues::growth_bound(&up, &low, n);
    let sizes = sizes || z > 0.0;
    let series = residues.map(|(t, size)| (if sizes { size } else { t }, size));
    let (sum, err) = pfq_alg_sum(series, settled, b.len().saturating_sub(a.len()), z, bound);
    if (sum, err) != (0.0, 0.0) {
        let bs = lower.iter().copied().map(Dd::new);
        let args = bs.clone().map(|bj| bj - am);
        let poles = args.clone().filter(is_pole);
        let sign = neg_one_pow(poles.clone().filter(|g| g.hi % 2.0 != 0.0).count());
        // (Gamma(a_i - a_m) has the residue 1 for a_i = a_m)
        let num = bs
            .chain(diffs.filter(|d| d.hi != 0.0))
            .chain(poles.map(|g| -g + 1.0));
        let den = others.chain(args.filter(|g| !is_pole(g)));
        let l = -ln_dd(z.abs()) * am;
        let mut y = gamma_series(num, den, (sign * sum, err), l)?;
        if sizes {
            y.m = y.m.abs();
        }
        *term = Some(y);
    }
    Some(())
}

/// What `pfq_alg` gives for each term
#[derive(Clone, Copy, PartialEq, Eq)]
enum Alg {
    /// Its value
    Value,
    /// A bound on its size
    Sizes,
    /// Its value for `p = q` or `p > q + 1`, with its series and its factor in double-double
    Dd,
}

/// The terms of the algebraic part of the asymptotic expansion of pFq for `p <= q` (DLMF
/// 16.11.2), with `|z|` instead of `-z`, or `None` if a factor of a term overflows or underflows,
/// or if a non-positive integer `a_i` differs from another `a_i` by an integer
///
/// For `p > q + 1` and `z < 0`, their series converge, and they sum to `pfq_borel`.
///
/// `b` contains the 1 of `n!` as its last parameter. The term of `a_m` is `|z|^-a_m` times
/// a `_{q+1}F_{p-1}` series in `(-1)^(p+q+1) / z`, with the parameters `a_m` and
/// `1 + a_m - b_j`, and `1 + a_m - a_i`. It stays `None` if it vanishes, where
/// `1 / Gamma(b_j - a_m)` does. The `a_i` that differ by integers have one term, by
/// `pfq_class_term`, which for `z > 0` is a bound on its size.
///
/// With `Alg::Sizes`, each term is the sum of the sizes of the terms of its series, times the
/// size of its factor. Their sum bounds the algebraic part at any phase of `z` and of its terms.
///
/// With `Alg::Dd`, the series and the factor of each term are in double-double, except for the
/// `a_i` that differ by integers.
fn pfq_alg<A: Params, B: Params>(
    a: &A,
    b: &B,
    z: f64,
    mode: Alg,
    terms: &mut [Option<Scaled>],
) -> Option<()> {
    let sizes = mode == Alg::Sizes;
    let (_, lower) = b.as_ref().split_last()?;
    let ln_z = ln_dd(z.abs());
    let bs = lower.iter().copied().map(Dd::new);
    for (m, &am) in a.as_ref().iter().enumerate() {
        let others = skip(a, m).map(Dd::new);
        let diffs = others.clone().map(|ai| ai - am);
        if diffs.clone().any(Dd::is_integer) {
            pfq_class_term(a.as_ref(), b.as_ref(), m, z, sizes, &mut terms[m])?;
            continue;
        }
        let num = bs.clone().chain(diffs);
        let den = others.chain(bs.clone().map(|bj| bj - am));
        if den.clone().any(|g| is_pole(&g)) {
            continue;
        }
        let l = -ln_z * am;
        let mut y = if mode == Alg::Dd {
            let (sum, err) = pfq_alg_series_dd(a.as_ref(), b.as_ref(), m, z);
            let y = Scaled::new(1.0, 0, err / sum.hi.abs())?.with_l(l);
            mul_gammas_in(num, den, mul_dd(y, sum), true)?
        } else {
            gamma_series(num, den, pfq_alg_series(a, b, m, z, sizes)?, l)?
        };
        if sizes {
            y.m = y.m.abs();
        }
        terms[m] = Some(y);
    }
    Some(())
}

/// 2F2 by its asymptotic expansion (DLMF 16.11.7)
///
/// Only the exponentially large part is used for `z > 0`, and only the algebraic part for
/// `z < 0`. The size of the other part is added to the error.
fn hyp2f2_asymptotic(a: [f64; 2], b: [f64; 2], z: f64) -> Option<Estimate> {
    // the algebraic part
    let alg = |a: [f64; 2]| {
        let mut terms = [None; 2];
        if a[1] == b[1] {
            // they cancel
            pfq_alg(&[a[0]], &[b[0], 1.0], z, Alg::Value, &mut terms)?;
        } else {
            pfq_alg(&a, &[b[0], b[1], 1.0], z, Alg::Value, &mut terms)?;
        }
        // for `z > 0`, the terms lack the phase of (-z)^(-a_m), whose real part is cos(pi a_m)
        // (the term of a class of `a_i` that differ by integers is a size, so it keeps its sign)
        let class = a[1] != b[1] && (Dd::new(a[0]) - Dd::new(a[1])).is_integer();
        let phased = z > 0.0 && !class;
        Some(terms.into_iter().zip(a).filter_map(move |(t, a)| {
            if !phased {
                return t;
            }
            // (a zero phase drops the term)
            let phase = cospi(a);
            t.filter(|_| phase != 0.0).map(|y| y.mul_f64(phase))
        }))
    };
    let (y, other) = if z > 0.0 {
        let alg = alg(a)?;
        // (no terms have no size, but terms that cancel leave their uncertainty)
        let other = match alg.clone().next() {
            None => None,
            Some(_) => Some(add_exp(alg)?),
        };
        (hyp2f2_exp(a, b, z)?, other)
    } else {
        (add_exp(alg(a)?)?, Some(hyp2f2_exp(a, b, z)?))
    };
    Some(Estimate::scaled(y, other))
}

/// Generalized hypergeometric function $_2F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// $$
/// \hyp 2 2 {a_1\enspace a_2} {b_1\enspace b_2} z
/// = \sum\_{n=0}^\infty
///     {\rpow{a_1}{n} \ \rpow{a_2}{n} \over \rpow{b_1}{n} \ \rpow{b_2}{n}}
///     {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1]. A
/// non-positive integer $a_i = -n$ makes it a polynomial of degree at most $n$. An $a_i$ equal to a
/// $b_j$ cancels it, also if it is a non-positive integer.
///
/// It is:
///
/// - $+\infty$ if a $b_j$ is a non-positive integer $-m$ and $z \ne 0$, unless it is a
///   polynomial of degree $\lt m$;
/// - $\pm\infty$ where it detects that it overflows, for large $\abs z$, and otherwise NaN;
/// - NaN if an argument is NaN;
/// - for an infinite argument, the limit as in [`hyppfq`].
///
/// # Returns NaN
///
/// Where its internal error estimate says that none of its digits are correct:
///
/// - for some large parameters, e.g. in random samples, for about 16% of $-10^3 < z < 0$ and 9% of
///   $0 < z < 10^3$ if all $\abs{a_i}, \abs{b_j} \le 100$;
/// - near a zero, or if the result underflows;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp3f2`]: Generalized hypergeometric function,
///   $\hyp 3 2 {a_1\enspace a_2\enspace a_3} {b_1\enspace b_2} z$
#[must_use]
pub fn hyp2f2(a1: f64, a2: f64, b1: f64, b2: f64, z: f64) -> f64 {
    Estimate::or_nan(hyp2f2_estimate(a1, a2, b1, b2, z))
}

/// `hyp2f2` and an estimate of its relative error
fn hyp2f2_estimate(a1: f64, a2: f64, b1: f64, b2: f64, z: f64) -> Option<Estimate> {
    if ![a1, a2, b1, b2, z].into_iter().all(f64::is_finite) {
        return pfq_estimate(&[a1, a2], &[b1, b2], z);
    } else if z == 0.0 {
        return Some(Estimate::exact(1.0));
    }
    let (mut a, mut b) = ([a1, a2], [b1, b2]);
    // a common parameter cancels, also at a pole; it becomes a_2 = b_2 = 1, which cancels like
    // the 1 of n! in the series, and which `hyp2f2_asymptotic` recognizes
    let pairs = [(0, 0), (0, 1), (1, 0), (1, 1)];
    if let Some((i, j)) = pairs.into_iter().find(|&(i, j)| a[i] == b[j]) {
        (a, b) = ([a[1 - i], 1.0], [b[1 - j], 1.0]);
        if a[0] == b[0] {
            // (with the rounding of `exp`)
            return Some(Estimate::rounded(z.exp(), f64::EPSILON));
        }
    }
    let Some(degree) = pfq_degree(&a, &b) else {
        return Some(Estimate::exact(f64::INFINITY));
    };
    let series = || Some(pfq_series_sorted(a, [b[0], b[1], 1.0], z));
    let asymptotic = || hyp2f2_asymptotic(a, b, z);
    // last, as it is the slowest
    let asymptotic_dd = || pfq_asymptotic_dd(&a, &[b[0], b[1], 1.0], z);
    if degree.is_finite() {
        series()
    } else if z.abs() < HYP2F2_ASYMPTOTIC_FIRST {
        series().or_better(asymptotic).or_better(asymptotic_dd)
    } else {
        asymptotic().or_better(series).or_better(asymptotic_dd)
    }
}

/// `v^(1 / k)` in double-double precision, for `v > 0`
fn root_dd(v: impl Into<Dd>, k: usize) -> Dd {
    let v = v.into();
    let s = v.hi.powf(1.0 / k as f64);
    let pow = |n: usize| (0..n).fold(Dd::ONE, |p, _| p * s);
    match k {
        1 => v,
        2 => v.sqrt(),
        // (one Newton step)
        _ => Dd::new(s) + (v - pow(k)) / (k as f64 * pow(k - 1).hi),
    }
}

/// Each `step`-th term of `iter`, from its first one on, with alternating signs
fn alternating(iter: impl Iterator<Item = f64>, step: usize) -> impl Iterator<Item = f64> {
    let signs = [1.0, -1.0].into_iter().cycle();
    iter.step_by(step).zip(signs).map(|(t, sign)| sign * t)
}

/// The exponential part of the asymptotic expansion of 1F2 (DLMF 16.11.1, 16.11.8), with
/// `x = 2 sqrt(|z|)` and `nu = a - b_1 - b_2 + 1/2`:
///
/// ```text
/// Gamma(b_1) Gamma(b_2) / Gamma(a) |z|^(nu/2) / sqrt(pi) S,
/// S = 1/2 e^x sum_k c_k x^-k                    (z > 0)
/// S = sum_k c_k x^-k cos(x + pi (nu - k) / 2)   (z < 0)
/// ```
///
/// For `z < 0`, it's the sum of the expansions with `e^(+-i x)`. By the differential
/// equation (DLMF 16.8.3) in `x`, `c_0 = 1` and `2 k c_k = A_1 c_{k-1} + A_0 c_{k-2}`, with
/// `s = nu + 1 - k`, `beta_j = 2 b_j - 2`, and
///
/// - `A_1 = (2 s + beta_1 + beta_2 + 1) (s + 1) + (s + beta_1) (s + beta_2)`,
/// - `A_0 = (s + 1) (s + 1 + beta_1) (s + 1 + beta_2)`.
fn hyp1f2_exp(a: f64, b: [f64; 2], z: f64) -> Option<Scaled> {
    let nu = Dd::new(a) - b[0] - b[1] + 0.5;
    let sqrt_z = Dd::new(z.abs()).sqrt();
    let (x, beta) = (2.0 * sqrt_z.hi, b.map(|bj| 2.0 * bj - 2.0));
    // the terms t_k = c_k x^-k
    let terms = recurrence_terms(|k, t1, t2| {
        let s = nu.hi + 1.0 - k;
        let a1 = (2.0 * s + beta[0] + beta[1] + 1.0) * (s + 1.0) + (s + beta[0]) * (s + beta[1]);
        let a0 = (s + 1.0) * (s + 1.0 + beta[0]) * (s + 1.0 + beta[1]);
        (a1 * t1 + a0 / x * t2) / (2.0 * k * x)
    });
    // (the factors of A_0 are negative from here on)
    let settled = (nu.hi + 2.0 + beta[0].max(beta[1]).max(0.0)).max(0.0);
    let l = nu * 0.5 * ln_dd(z.abs()) - LN_PI_DD * 0.5;
    // (the term ratios are about k / x, and without small divisors, the terms don't grow again
    // before `settled + x`)
    let y = if z > 0.0 {
        let (sum, err) = asymptotic_sum(terms, settled, x, |_| 0.0);
        Scaled::new(sum, -1, err / sum.abs())?.with_l(l + sqrt_z * 2.0)
    } else {
        // the terms of even and odd k, with the cosines as +-cos(theta) and +-sin(theta)
        let [(even, even_err), (odd, odd_err)] = [0, 1].map(|parity| {
            let terms = alternating(terms.clone().skip(parity), 2);
            asymptotic_sum(terms, settled / 2.0, x / 2.0, |_| 0.0)
        });
        // theta / pi = x / pi + nu / 2 in double-double, reduced to [-1, 1], and a bound of the
        // error of theta
        let t = sqrt_z * FRAC_2_PI_DD + nu * 0.5;
        let t = t - 2.0 * (0.5 * t.hi).round();
        let dtheta = 4.0 * EPSILON_SQ * (x + PI * nu.hi.abs());
        // (to first order in t.lo)
        let (cos, sin) = (cospi(t.hi), sinpi(t.hi));
        let (cos, sin) = (cos - PI * t.lo * sin, sin + PI * t.lo * cos);
        let m = even * cos + odd * sin;
        // (and the roundings of `cospi`, `sinpi`, and the products and sum)
        let err = even_err * cos.abs()
            + odd_err * sin.abs()
            + dtheta * (even.abs() + odd.abs())
            + 3.0 * f64::EPSILON * ((even * cos).abs() + (odd * sin).abs());
        Scaled::new(m, 0, err / m.abs())?.with_l(l)
    };
    mul_gammas(b.map(Dd::new), [Dd::new(a)], y)
}

/// 1F2 by its asymptotic expansion (DLMF 16.11.8)
///
/// For `z > 0`, only the exponentially large part is used, and the size of the algebraic part
/// is added to the error. The terms with `e^-x` are smaller than the truncation error.
fn hyp1f2_asymptotic(a: f64, b: [f64; 2], z: f64) -> Option<Estimate> {
    let exp = hyp1f2_exp(a, b, z)?;
    // the algebraic part, `None` if it vanishes
    let mut alg = [None];
    pfq_alg(&[a], &[b[0], b[1], 1.0], z, Alg::Value, &mut alg)?;
    if z > 0.0 {
        Some(Estimate::scaled(exp, alg[0]))
    } else {
        Some(Estimate::scaled(add_exp(alg.into_iter().flatten().chain([exp]))?, None))
    }
}

/// Generalized hypergeometric function $_1F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// $$
/// \hyp 1 2 a {b_1\enspace b_2} z
/// = \sum\_{n=0}^\infty {\rpow a n \over \rpow{b_1}{n} \ \rpow{b_2}{n}} {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1]. A
/// non-positive integer $a = -n$ makes it a polynomial of degree $n$. An $a$ equal to a $b_j$
/// cancels it, also if it is a non-positive integer.
///
/// It is:
///
/// - $+\infty$ if a $b_j$ is a non-positive integer $-m$ and $z \ne 0$, unless it is a
///   polynomial of degree $\lt m$;
/// - $\pm\infty$ where it detects that it overflows, for large $\abs z$, and otherwise NaN;
/// - NaN if an argument is NaN;
/// - for an infinite argument, the limit as in [`hyppfq`].
///
/// # Returns NaN
///
/// Where its internal error estimate says that none of its digits are correct:
///
/// - for about half of $z < -10^{40}$;
/// - for some large parameters, e.g. in random samples, for about 13% of $-10^5 < z < -10^3$ and
///   7% of $0 < z < 10^5$ if all $\abs a, \abs{b_j} \le 100$;
/// - near a zero, or if the result underflows;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function, $_0F_1(b; z)$
/// - [`hyp0f2`]: Generalized hypergeometric function,
///   $_0F_2\left[b_1\enspace b_2\middle\| z\right]$
/// - [`hyp2f2`]: Generalized hypergeometric function,
///   $\hyp 2 2 {a_1\enspace a_2} {b_1\enspace b_2} z$
#[must_use]
pub fn hyp1f2(a: f64, b1: f64, b2: f64, z: f64) -> f64 {
    Estimate::or_nan(hyp1f2_estimate(a, b1, b2, z))
}

/// `hyp1f2` and an estimate of its relative error
fn hyp1f2_estimate(a: f64, b1: f64, b2: f64, z: f64) -> Option<Estimate> {
    if ![a, b1, b2, z].into_iter().all(f64::is_finite) {
        return pfq_estimate(&[a], &[b1, b2], z);
    } else if z == 0.0 {
        return Some(Estimate::exact(1.0));
    }
    let (mut a, mut b) = (a, [b1, b2]);
    // a common parameter cancels, also at a pole; it becomes a = b_2 = 1, which cancels like the
    // 1 of n! in the series
    if let Some(j) = b.iter().position(|&bj| bj == a) {
        (a, b) = (1.0, [b[1 - j], 1.0]);
    }
    let Some(degree) = pfq_degree(&[a], &b) else {
        return Some(Estimate::exact(f64::INFINITY));
    };
    let series = || Some(pfq_series_sorted([a], [b[0], b[1], 1.0], z));
    let asymptotic = || hyp1f2_asymptotic(a, b, z);
    if degree.is_finite() {
        series()
    } else if z > -HYP1F2_ASYMPTOTIC_FIRST {
        series().or_better(asymptotic)
    } else {
        asymptotic().or_better(series)
    }
}

/// Generalized hypergeometric function $_0F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// $$
/// _0F_2\left[b_1\enspace b_2\middle\| z\right]
/// = \sum\_{n=0}^\infty {1 \over \rpow{b_1}{n} \ \rpow{b_2}{n}} {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1].
///
/// It is:
///
/// - $+\infty$ if a $b_j$ is a non-positive integer and $z \ne 0$;
/// - $\pm\infty$ where it detects that it overflows, for large $\abs z$, and otherwise NaN;
/// - NaN if an argument is NaN;
/// - for an infinite argument, the limit as in [`hyppfq`].
///
/// # Returns NaN
///
/// Where its internal error estimate says that none of its digits are correct:
///
/// - for most $z < -10^{40}$;
/// - for some large parameters, e.g. in random samples, for about 7% of $-10^8 < z < 0$ if
///   $\abs{b_1}, \abs{b_2} \le 100$;
/// - near a zero;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function, $_0F_1(b; z)$
/// - [`hyp1f2`]: Generalized hypergeometric function, $\hyp 1 2 a {b_1\enspace b_2} z$
#[must_use]
pub fn hyp0f2(b1: f64, b2: f64, z: f64) -> f64 {
    Estimate::or_nan(pfq_estimate(&[], &[b1, b2], z))
}

/// Generalized hypergeometric function $_2F_0$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// $$
/// _2F_0\left[a_1\enspace a_2\middle\| z\right]
/// = \sum\_{n=0}^\infty \rpow{a_1}{n} \ \rpow{a_2}{n} {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1], and for
/// $z < 0$, where the series diverges, its Borel sum, as in mpmath [^2]. A non-positive integer
/// $a_i = -n$ makes it a polynomial of degree at most $n$, defined for all $z$.
///
/// It is:
///
/// - $\pm\infty$ where it detects that it overflows, for large $\abs z$, and otherwise NaN;
/// - NaN for $z > 0$, where the series diverges, unless it is a polynomial, and if an argument is
///   NaN;
/// - for an infinite argument, the limit as in [`hyppfq`].
///
/// # Returns NaN
///
/// Where its internal error estimate says that none of its digits are correct:
///
/// - if the result underflows, e.g. in random samples, for about 2% of $-10^6 < z < 0$ if
///   $\abs{a_1}, \abs{a_2} \le 100$, and for more of $z < -10^{12}$, or overflows undetected;
/// - near a zero, also an exact one of a polynomial;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.6.E21>
///
/// # See also
/// - [`hyp1f0`](crate::hyp1f0): Hypergeometric function, $_1F_0\left[a\middle\| z\right]$
/// - [`hypu`](crate::hypu): Confluent hypergeometric function, $U(a_1,a_2,x)$
#[must_use]
pub fn hyp2f0(a1: f64, a2: f64, z: f64) -> f64 {
    Estimate::or_nan(pfq_estimate(&[a1, a2], &[], z))
}

/// The terms `t_k = c_k x^-k` of the exponential part of the asymptotic expansion of
/// pFq for `p <= q` (DLMF 16.11.1), where `b` contains the 1 of `n!`
///
/// With `D = x d/dx` and `kappa = q - p + 1`, the differential equation (DLMF 16.8.3) is
/// `(prod_j (D + kappa (b_j - 1)) - x^kappa prod_i (D + kappa a_i)) w = 0`, and each factor
/// `D + g` maps `e^x x^s` to `e^x (x^(s+1) + (s + g) x^s)`. So `c_0 = 1` and
/// `kappa k c_k = sum_{m=1}^q d_{m+1}(nu - k + m) c_{k-m}`, where `d_j(s)` is the
/// coefficient of `e^x x^(s+q+1-j)` in the image of `e^x x^s`.
fn pfq_exp_terms(a: &[f64], b: &[f64], nu: f64, x: f64) -> impl Iterator<Item = f64> + Clone {
    let (q, kappa) = (b.len() - 1, (b.len() - a.len()) as f64);
    let left: Vec<f64> = b.iter().map(|bj| kappa * (bj - 1.0)).collect();
    let right: Vec<f64> = a.iter().map(|ai| kappa * ai).collect();
    // the coefficients of e^x x^(s+n-j) in prod_(i<n) (D + g_i) e^x x^s, for j <= q + 1
    let coefs = move |g: &[f64], s: f64| {
        let mut c = vec![0.0; q + 2];
        c[0] = 1.0;
        for (i, &gi) in g.iter().enumerate() {
            for j in (1..=i + 1).rev() {
                c[j] += c[j - 1] * (s + (i + 1 - j) as f64 + gi);
            }
        }
        c
    };
    // t_k, and the sums so far of t_(k+1), ..., t_(k+q)
    let (mut k, mut t, mut sums) = (0.0, 1.0, vec![0.0; q]);
    from_fn(move || {
        let (term, l, r) = (t, coefs(&left, nu - k), coefs(&right, nu - k));
        let mut scale = 1.0;
        for (m, sum) in sums.iter_mut().enumerate() {
            scale /= x;
            *sum += (l[m + 2] - r[m + 2]) * scale * t;
        }
        k += 1.0;
        t = sums[0] / (kappa * k);
        sums.rotate_left(1);
        sums[q - 1] = 0.0;
        Some(term)
    })
}

/// `cos(pi / k)` and `sin(pi / k)` in double-double precision, for `k >= 2`
fn cos_sin_pi_div(k: usize) -> (Dd, Dd) {
    let phi = 2.0 / (FRAC_2_PI_DD * k as f64);
    // n, and the terms n - 1 of their Taylor series, of the degrees 2 n - 2 and 2 n - 1, which
    // are below 2^-106 from n = 21 on
    let terms = successors(Some((1.0, Dd::ONE, phi)), |&(n, c, _)| {
        let next = -(c * phi * phi) / ((2.0 * n - 1.0) * (2.0 * n));
        Some((n + 1.0, next, next * phi / (2.0 * n + 1.0)))
    });
    let zero = (Dd::ZERO, Dd::ZERO);
    terms
        .take(20)
        .fold(zero, |(cos, sin), (_, c, s)| (cos + c, sin + s))
}

/// The sum of the two exponential parts of `pfq_exp` for `z < 0` and `p < q`, which have
/// `x e^(+-pi i / kappa)` as `x` (DLMF 16.11.8, 16.11.9), relative to
/// `e^(x cos(pi / kappa))`:
///
/// ```text
/// 2 sum_k t_k cos(theta - pi k / kappa),
/// theta = x sin(pi / kappa) + pi nu / kappa
/// ```
///
/// with its absolute error, and `cos(pi / kappa)`, where `terms` are the `t_k = c_k x^-k`
/// and `root` is `|z|^(1 / kappa)`.
fn pfq_exp_cos(
    terms: &(impl Iterator<Item = f64> + Clone),
    settled: f64,
    kappa: usize,
    root: Dd,
    nu: Dd,
) -> ((f64, f64), Dd) {
    let kf = kappa as f64;
    let x = kf * root.hi;
    let (cos, sin) = cos_sin_pi_div(kappa);
    // theta / pi in double-double, reduced to [-1, 1], and a bound of the error of theta
    let turns = root * kf * sin * (FRAC_2_PI_DD * 0.5) + nu / kf;
    let t = turns - 2.0 * (0.5 * turns.hi).round();
    let dtheta = 4.0 * EPSILON_SQ * (x + PI * nu.hi.abs());
    // the terms k = r mod kappa have the cosines +-cos(theta - pi r / kappa)
    let sums = (0..kappa).map(|r| {
        let sub = alternating(terms.clone().skip(r), kappa);
        let (s, s_err) = asymptotic_sum(sub, settled / kf, x / kf, |_| 0.0);
        let u = t.hi - r as f64 / kf;
        // (to first order in t.lo)
        let c = cospi(u) - PI * t.lo * sinpi(u);
        // (and the roundings of `u`, `cospi`, `sinpi`, the product, and the sum)
        let rounding = (dtheta + 2.0 * f64::EPSILON) * s.abs() + 3.0 * f64::EPSILON * (s * c).abs();
        [s * c, s_err * c.abs() + rounding, s.abs()]
    });
    let [sum, err, size] = sums.fold([0.0; 3], |u, v| [u[0] + v[0], u[1] + v[1], u[2] + v[2]]);
    // the other exponentials of DLMF 16.11.9, with x e^(+-3 pi i / kappa), relative to these
    let next = 2.0 * (x * (cospi(3.0 / kf) - cos.hi)).exp();
    let others = if kappa == 2 { 0.0 } else { next };
    ((2.0 * sum, 2.0 * (err + others * size)), cos)
}

/// The exponential part of the asymptotic expansion of pFq for `p <= q` (DLMF 16.11.1),
/// with `|z|^nu` instead of `z^nu`, where `b` contains the 1 of `n!`:
///
/// ```text
/// prod_j Gamma(b_j) / prod_i Gamma(a_i)
///     (2 pi)^((1 - kappa) / 2) kappa^(-nu - 1/2) e^x |x|^nu sum_k c_k x^-k,
/// nu = sum_i a_i - sum_j b_j + (kappa - 1) / 2
/// ```
///
/// with `kappa = q - p + 1` and `x = kappa z^(1 / kappa)`, which is real for `z > 0` or
/// `p = q`. Otherwise, it is the sum of the two parts with
/// `x = kappa |z|^(1 / kappa) e^(+-pi i / kappa)`.
fn pfq_exp(a: &[f64], b: &[f64], z: f64) -> Option<Scaled> {
    let (_, lower) = b.split_last()?;
    let kappa = b.len() - a.len();
    let kf = kappa as f64;
    let nu = sum_dd(a) - sum_dd(lower) + 0.5 * (kf - 1.0);
    let real = z > 0.0 || kappa == 1;
    let root = root_dd(z.abs(), kappa) * if real { z.signum() } else { 1.0 };
    let x = kf * root.hi;
    // (the factors of the last coefficient are negative from here on)
    let beta = lower.iter().fold(0.0_f64, |m, bj| m.max(kf * (bj - 1.0)));
    let settled = (nu.hi + lower.len() as f64 + beta).max(0.0);
    let terms = pfq_exp_terms(a, b, nu.hi, x);
    let (series, growth) = if real {
        let (sum, err) = asymptotic_sum(terms, settled, x, |_| 0.0);
        // the other exponentials of DLMF 16.11.8 and 16.11.9, with x e^(+-2 pi i / kappa),
        // relative to this one
        let next = 2.0 * (x * (cospi(2.0 / kf) - 1.0)).exp();
        let others = if kappa == 1 { 0.0 } else { next };
        ((sum, err + others * sum.abs()), Dd::ONE)
    } else {
        pfq_exp_cos(&terms, settled, kappa, root, nu)
    };
    // (2 pi)^((1 - kappa) / 2) kappa^(-1/2) |z|^(nu / kappa) e^x
    let ln_power = nu / kf * ln_dd(z.abs()) - ln_dd(kf) * 0.5;
    let l = ln_power - LN_SQRT_2PI_DD * (kf - 1.0) + root * kf * growth;
    let num = lower.iter().copied().map(Dd::new);
    let den = a.iter().copied().map(Dd::new);
    gamma_series(num, den, series, l)
}

/// pFq for `p <= q` by its asymptotic expansion (DLMF 16.11.7 to 16.11.9), where `b` contains
/// the 1 of `n!`
///
/// Only the exponentially large part is used for `z > 0`, and only the algebraic part for
/// `z < 0` and `p = q`. The size of the other part is added to the error, for `z > 0` as the
/// sizes of the terms of the algebraic part. For `z < 0` and `p < q`, it is the sum of both
/// parts.
fn pfq_asymptotic<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<Estimate> {
    let (a_ref, b_ref) = (a.as_ref(), b.as_ref());
    let kappa = b_ref.len() - a_ref.len();
    let mut terms = vec![None; a_ref.len()];
    let mode = if z > 0.0 { Alg::Sizes } else { Alg::Value };
    pfq_alg(a, b, z, mode, &mut terms)?;
    let alg = terms.into_iter().flatten();
    let exp = pfq_exp(a_ref, b_ref, z)?;
    if z > 0.0 {
        Some(Estimate::scaled(exp, add_exp(alg)))
    } else if kappa == 1 {
        Some(Estimate::scaled(add_exp(alg)?, Some(exp)))
    } else {
        // (DLMF 16.11.9 for p <= q - 2 has the exponential part only; with the algebraic part,
        // it is theorem 3 of Paris, https://arxiv.org/abs/1708.04824)
        Some(Estimate::scaled(add_exp(alg.chain([exp]))?, None))
    }
}

/// `add_exp` in double-double
///
/// The mantissas of the terms are scaled to `[1/2, 1)` first: a term whose `exp` underflows to
/// zero is then below `2^-1021` of the largest one, which is within the error bound of that one.
fn add_exp_dd(terms: impl IntoIterator<Item = Scaled, IntoIter: Clone>) -> Option<Scaled> {
    let terms = terms.into_iter().map(|mut y| {
        let (m, k) = frexp_f64(y.m);
        (y.m, y.e) = (m, y.e + k);
        y
    });
    add_exp_in(terms, Dd::exp, Dd::EPSILON)
}

/// The sum of the terms of `pfq_alg` in double-double, for `p = q` or `p > q + 1`
fn pfq_alg_dd<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<Scaled> {
    let mut terms = vec![None; a.as_ref().len()];
    pfq_alg(a, b, z, Alg::Dd, &mut terms)?;
    add_exp_dd(terms.into_iter().flatten())
}

/// pFq for `p > q + 1` and `z < 0` by `pfq_alg_dd`, or `None` if it may or may not overflow
fn pfq_borel_dd<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<Estimate> {
    Estimate::scaled_certain(pfq_alg_dd(a, b, z)?)
}

/// pFq for `p = q` and `z < 0` by the algebraic part of `pfq_asymptotic` in double-double, where
/// its terms cancel, e.g. for nearly equal `a_i`, or the terms of their series do
fn pfq_asymptotic_dd<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<Estimate> {
    if z > 0.0 || a.as_ref().len() + 1 != b.as_ref().len() {
        return None;
    }
    let alg = pfq_alg_dd(a, b, z)?;
    Some(Estimate::scaled(alg, Some(pfq_exp(a.as_ref(), b.as_ref(), z)?)))
}

/// The limit at `z = +-inf` for finite parameters if cheap, with `b` without the 1 of `n!`
fn pfq_limit(a: &[f64], b: &[f64], z: f64, degree: f64) -> Option<f64> {
    if z.is_finite() {
        None
    } else if degree == 0.0 {
        Some(1.0)
    } else if degree.is_finite() {
        // the sign of the last term: the parity of its negative factors
        let odd = |v: &f64| (-v).ceil().clamp(0.0, degree) % 2.0;
        let parity = a.iter().chain(b).chain([&z]).map(odd).sum::<f64>() % 2.0;
        Some(f64::INFINITY * (1.0 - 2.0 * parity))
    } else if z > 0.0 && a.len() <= b.len() {
        // the sign of the exponentially large term (DLMF 16.11)
        Some(f64::INFINITY * late_sign(a, b))
    } else if (b.len()..=b.len() + 1).contains(&a.len()) && a.iter().all(|&v| v > 0.0) {
        // all terms of DLMF 16.11 vanish
        (z < 0.0).then_some(0.0)
    } else {
        None
    }
}

/// pFq for `p > q + 1` and `z < 0`, where the series diverges: its Borel sum, where `b` contains
/// the 1 of `n!`
///
/// That is the Meijer G function of DLMF 16.5.1, whose residues at the poles of the
/// `Gamma(a_i + s)` sum to `pfq_alg`. For small `|z|`, where those cancel, it is the series up
/// to its smallest term, as for 2F0, or else its Weniger transformation, last in double-double,
/// which is the slowest. The residues come first, as the steps of Weniger can stall where some
/// terms nearly vanish, e.g. for an `a_i` near a non-positive integer. Where none of these is
/// accurate, the residues are in double-double, as for 2F0.
fn pfq_borel<P: Params>(a: &P, b: &P, z: f64) -> Option<Estimate> {
    let (a_ref, b_ref) = (a.as_ref(), b.as_ref());
    // the term ratios are about n^d |z|, so that the terms are smallest near n = |z|^(-1/d)
    let d = (a_ref.len() - b_ref.len()) as f64;
    // (capped as for 2F0)
    let scale = z.abs().powf(-1.0 / d).min(f64::EPSILON.recip());
    let terms = || once(1.0).chain(pfq_terms(a_ref, b_ref, z));
    let settled = pfq_settled(a_ref, b_ref).ceil();
    let series = || Some(Estimate::absolute(asymptotic_sum(terms(), settled, scale, |_| 0.0)));
    // also with the rounding errors of the terms, as for 2F0
    let weniger = |weniger: (f64, f64), eps: f64| weniger_drift(weniger, terms(), settled, eps);
    let weniger_f64 = || Some(weniger(pfq_weniger::<f64>(a, b, z)?, f64::EPSILON));
    // as for 2F0, where the rounding errors in f64 dominate
    let weniger_dd = || Some(weniger(pfq_weniger::<Dd>(a, b, z)?, Dd::EPSILON));
    let residues = || {
        let mut terms = vec![None; a_ref.len()];
        pfq_alg(a, b, z, Alg::Value, &mut terms)?;
        Some(Estimate::scaled(add_exp(terms.into_iter().flatten())?, None))
    };
    series()
        .or_better(residues)
        .or_better(weniger_f64)
        .or_better(weniger_dd)
        .or_better(|| pfq_borel_dd(a, b, z))
}

/// The decays `exp(-m)` of `v_k` from where `hyp2f0_miller` recurs down: below `f64::EPSILON`, and
/// further, as the difference of the two results estimates the truncation error
const MILLER_DECAY: [f64; 2] = [40.0, 60.0];

/// 2F0 for `z < 0` by Miller's algorithm (DLMF 13.29.5 to 13.29.7), and an estimate of its
/// absolute error, or `None` if it would take more than `PFQ_DD_MAXITER` steps
///
/// With `x = -1/z`, the `v_k = (a1)_k (a2)_k / k! U(a1 + k, a1 - a2 + 1, x)` are the minimal
/// solution of
///
/// `(a1 + k - 1) (a2 + k - 1) v_(k-1) = k ((a1 + a2 + 2 k - 1 + x) v_k - (k + 1) v_(k+1))`
///
/// and sum to `x^(-a1)`, so that 2F0 is `v_0 / sum_k v_k` (DLMF 13.6.21): without a gamma
/// function, a power, or the rounded `x` and `a1 - a2 + 1`. It is in double-double, with a
/// first-order estimate of its rounding errors, which the recurrence amplifies for `a1, a2 < 0`.
/// The estimate is not a bound: it propagates the errors of the steps with one choice of signs.
///
/// As the sum cancels unless `a1, a2 > 0`, it is that of `a1 + shift` and `a2 + shift`, and the
/// `f_k = 2F0(a1 + k, a2 + k;; z)` then go down from `k = shift` by (DLMF 13.3.7)
///
/// `f_(k-1) = (1 - (a1 + a2 + 2 k - 1) z) f_k - (a1 + k) (a2 + k) z^2 f_(k+1)`
///
/// with the rounding errors of its steps and of `f_shift` and `f_(shift+1)`, each times the
/// sensitivity of `f_0` to it.
fn hyp2f0_miller(a1: f64, a2: f64, z: f64, shift: f64) -> Option<(f64, f64)> {
    let x = Dd::new(-1.0) / z;
    let sum = x + a1 + a2 + (2.0 * shift - 1.0);
    // v_0 / sum_k v_k from v_n = 1 and v_(n+1) = 0, and its relative rounding error
    let ratio = |n: f64| {
        let (mut v, mut e) = ([Dd::ONE, Dd::ZERO], [0.0; 2]);
        let (mut total, mut e_total) = (v[0], 0.0);
        for k in (1..=n as u32).rev().map(f64::from) {
            let (p, q) = ((sum + 2.0 * k) * v[0], v[1] * (k + 1.0));
            let den = (Dd::new(a1) + (shift + k - 1.0)) * (Dd::new(a2) + (shift + k - 1.0)) / k;
            // the rounding error of this step, with the sign of the errors so far
            let rounding = Dd::EPSILON * (p.hi.abs() + q.hi.abs());
            let propagated = (sum.hi + 2.0 * k) * e[0] - (k + 1.0) * e[1];
            e = [(propagated + rounding.copysign(propagated)) / den.hi, e[0]];
            v = [(p - q) / den, v[0]];
            total = total + v[0];
            e_total += e[0].abs();
        }
        // f_shift and f_(shift+1) times sum_k v_k (DLMF 13.6.21)
        let den = (Dd::new(a1) + shift) * (Dd::new(a2) + shift) * -z;
        let mut f = [v[0], v[1] / den];
        // and their rounding errors
        let e_f1 = (e[1] / den.hi).abs() + Dd::EPSILON * f[1].hi.abs();
        let mut steps = Vec::with_capacity(shift as usize);
        for k in (1..=shift as u32).rev().map(f64::from) {
            let c = (Dd::new(a1) + a2 + (2.0 * k - 1.0)) * -z + 1.0;
            let d = (Dd::new(a1) + k) * (Dd::new(a2) + k) * z * z;
            let (p, q) = (c * f[0], d * f[1]);
            steps.push((c.hi, d.hi, Dd::EPSILON * (p.hi.abs() + q.hi.abs())));
            f = [p - q, f[0]];
        }
        // the sensitivities of f_0 to f_(k-1) and f_k, up from k = 1
        let (mut w, mut e_f) = ([1.0, 0.0], 0.0);
        for (c, d, rounding) in steps.into_iter().rev() {
            e_f += (w[0] * rounding).abs();
            w = [w[0] * c + w[1], -w[0] * d];
        }
        e_f += (w[0] * e[0]).abs() + (w[1] * e_f1).abs();
        ((f[0] / total).hi, (e_f / f[0].hi).abs() + (e_total / total.hi).abs())
    };
    // A heuristic: v_k is about k^((a1 + a2) / 2 - 1) K_(a1 - a2)(2 sqrt(k x)) for k >> a1, a2
    // (DLMF 13.8.11), which decays by exp(-m) once 2 sqrt(k x) is |a1 - a2| + m, but for the power
    // of k (left to the truncation error), from where a1 + k and a2 + k are positive (and over at
    // least m steps)
    let settled = (-a1 - shift).max(-a2 - shift).max(0.0);
    let steps = |m: f64| settled + m + ((a1 - a2).abs() + m).powi(2) / (4.0 * x.hi);
    let [short, long] = MILLER_DECAY.map(steps);
    let ((y_short, _), (y, err)) =
        (long + shift <= f64::from(PFQ_DD_MAXITER)).then(|| (ratio(short), ratio(long)))?;
    // (and the rounding of the quotient and of `y`; a subnormal `y` of the shifted run has lost
    // digits that its estimate misses)
    (shift == 0.0 || y.is_normal())
        .then(|| (y, (err + f64::EPSILON) * y.abs() + (y - y_short).abs()))
}

/// Generalized hypergeometric function $_pF_q$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// For $p \le q$, and for $\abs z < 1$ if $p = q + 1$,
///
/// $$
/// \hyp p q {a_1\enspace \cdots\enspace a_p} {b_1\enspace \cdots\enspace b_q} z
/// = \sum\_{n=0}^\infty
///     {\rpow{a_1}{n} \cdots \rpow{a_p}{n} \over \rpow{b_1}{n} \cdots \rpow{b_q}{n}}
///     {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1]. For
/// $p = q + 1$ and $z \le -1$, it is the analytic continuation, and for $p > q + 1$ and $z < 0$,
/// where the series diverges, its Borel sum, as in mpmath [^2] [^3]. A non-positive integer
/// $a_i = -n$ makes it a polynomial of degree at most $n$, defined for all $z$.
///
/// Each $a_i$ equal to a $b_j$ cancels it, also if it is a non-positive integer or infinite. These
/// $p$ and $q$ also have their own function:
///
/// | $p$ | $q = 0$                   | $q = 1$                   | $q = 2$    |
/// | --- | ------------------------- | ------------------------- | ---------- |
/// | 0   | [`hyp0f0`](crate::hyp0f0) | [`hyp0f1`](crate::hyp0f1) | [`hyp0f2`] |
/// | 1   | [`hyp1f0`](crate::hyp1f0) | [`hyp1f1`](crate::hyp1f1) | [`hyp1f2`] |
/// | 2   | [`hyp2f0`]                | [`hyp2f1`](crate::hyp2f1) | [`hyp2f2`] |
/// | 3   |                           |                           | [`hyp3f2`] |
///
/// It is:
///
/// - 1 at $z = 0$ unless an argument is NaN, and usually for a $b_j = +\infty$ and a finite $z$
///   if $p \le q + 1$ or it is a polynomial, unless it is complex;
/// - $+\infty$ if a $b_j$ is a non-positive integer $-m$ and $z \ne 0$, unless it is a
///   polynomial of degree $\lt m$;
/// - $+\infty$ at $z = 1$ for $_1F_0$ if $a_1 > 0$, and $\pm\infty$, the sign of its limit from
///   below, for $p = q + 1 \ge 3$ after cancellations and finite parameters if
///   $\sum_j b_j \le \sum_i a_i$, unless it is a polynomial or that can't be determined (then
///   `None`);
/// - infinite at $z = \pm\infty$ if it is a polynomial of degree $\ge 1$, and at $z = +\infty$
///   for $p \le q$ if it isn't a polynomial;
/// - 0 at $z = -\infty$ for $q \le p \le q + 1$ if all $a_i > 0$;
/// - NaN for $p = q + 1$ and $z > 1$, where it is complex, unless it is a polynomial,
///   $(1 - z)^{-a_1}$ with an integer $a_1$, or a $_2F_1$ that is a rational function of $z$,
///   and if an argument is NaN.
///
/// # Returns `None`
///
/// - for $p > q + 1$ and $z > 0$, unless it is a polynomial;
/// - for $p > q + 1$ and $z < 0$ where its Borel sum can't be computed accurately, e.g. for some
///   $\abs z \lesssim 1$ and $\abs{a_i} \gtrsim 10$, or $a_i$ that are nearly equal or differ
///   by integers, and where it underflows;
/// - for the other $p \le q + 1$ where it can't be computed accurately, as for [`hyp3f2`] if
///   $p = q + 1$, for [`hyp2f2`] if $p = q$, and for [`hyp1f2`] if $p < q$;
/// - for the other infinite arguments, unless $p \le 2$ and $q \le 1$ where the function in the
///   table has a limit there.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.6.E21>
///
/// [^3]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.18.E1>
#[must_use]
pub fn hyppfq(a: &[f64], b: &[f64], z: f64) -> Option<f64> {
    pfq_estimate(a, b, z)?.accurate()
}

/// The signature of `xsf_hypu_digits`, U with the digits that `chgu` estimates
type Hypu = fn(f64, f64, f64) -> (f64, i32, i32);

/// `hyppfq` and an estimate of its relative error, or `None` where it has no method
fn pfq_estimate(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let exact = |y| Some(Estimate::exact(y));
    if a.iter().chain(b).chain([&z]).any(|v| v.is_nan()) {
        return exact(f64::NAN);
    } else if z == 0.0 {
        return exact(1.0);
    }
    // each b_j cancels an equal a_i, if any
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    let mut cancel = |bj: &f64| a.iter().position(|ai| ai == bj).map(|i| a.remove(i));
    b.retain(|bj| cancel(bj).is_none());
    let Some(degree) = pfq_degree(&a, &b) else {
        return exact(f64::INFINITY);
    };
    let complex = a.len() == b.len() + 1 && z > 1.0 && degree.is_infinite();
    if a.iter().chain(&b).any(|v| v.is_infinite()) {
        // the terms n >= 1 vanish for b_j = +inf
        let finite = a.iter().chain([&z]).all(|v| v.is_finite());
        let converges = a.len() <= b.len() + 1 || degree.is_finite();
        let one = finite && converges && !b.contains(&f64::NEG_INFINITY) && !complex;
        return one.then_some(Estimate::exact(1.0));
    } else if let Some(y) = pfq_limit(&a, &b, z, degree) {
        return exact(y);
    }
    // a polynomial is short enough for double-double
    let polynomial = degree.is_finite().then(|| {
        let (y, err, ..) = pfq_sum::<Dd>(&a, &[&b[..], &[1.0]].concat(), z, PFQ_DD_MAXITER);
        Estimate::absolute((f64::from(y), err))
    });
    if polynomial.is_some_and(Estimate::is_accurate) {
        return polynomial;
    }
    // the xsf functions are NaN where they are inaccurate
    let xsf = |y: f64| exact(y).filter(|_| !y.is_nan());
    let y = match (a.as_slice(), b.as_slice()) {
        (&[], &[]) => exact(crate::hyp0f0(z)),
        // +inf where it diverges, and NaN where it is complex
        (&[a], &[]) if z == 1.0 && a > 0.0 => exact(f64::INFINITY),
        (&[a], &[]) if !complex || a.fract() == 0.0 => exact(crate::hyp1f0(a, z)),
        // real only where it is rational for the exact parameters: (1 - z)^(c - a - b) times a
        // polynomial, for an integer b and c - a (DLMF 15.8.1), which `hyp2f1` can round to
        (&[a1, a2], &[b]) if complex => {
            let rational = |a: f64, other: f64| is_pole(&(Dd::new(b) - a)) && other.fract() == 0.0;
            let y = crate::hyp2f1(a1, a2, b, z);
            let real = y.is_finite() && (rational(a1, a2) || rational(a2, a1));
            exact(if real { y } else { f64::NAN })
        },
        _ if complex => exact(f64::NAN),
        (&[], &[b]) => xsf(crate::hyp0f1(b, z)),
        (&[a], &[b]) => xsf(crate::hyp1f1(a, b, z)),
        (&[a1, a2], &[b]) => xsf(crate::hyp2f1(a1, a2, b, z)),
        _ if z.is_infinite() => None,
        // the Borel sum of the divergent series (DLMF 13.6.21), as in mpmath, if `chgu` estimates
        // that U has the digits of PFQ_TOL; otherwise, or if a factor isn't normal (their product
        // could be an over- times an underflow), the series up to its smallest term, as for tiny
        // |z|
        (&[a1, a2], &[]) if z < 0.0 => {
            // the smaller |a_i| first, as x^a1 then overflows less often
            let (a1, a2) = if a1.abs() <= a2.abs() {
                (a1, a2)
            } else {
                (a2, a1)
            };
            let x = -1.0 / z;
            // the rounding error of x (as x z is near -1, x z + 1 is exact)
            let (xz, e) = two_prod(x, z);
            let dx = -((xz + 1.0) + e) / z;
            let borel = |a1: f64, a2: f64, hypu: Hypu, margin: f64| {
                const T: f64 = pow2(20);
                const POLE_SPACINGS: f64 = 8.0;
                // with the digits that chgu estimates
                let y = |b: f64, x: f64| {
                    let (p, (u, _, digits)) = (x.powf(a1), hypu(a1, b, x));
                    let digits = f64::from(digits) - margin;
                    (p.is_normal() && u.is_normal()).then_some((p, u, digits))
                };
                let (d, e1) = two_sum(a1, -a2);
                let (b, e2) = two_sum(d, 1.0);
                // chgu takes U as a polynomial where its `a1 - b + 1` (that should be a2) rounds to
                // a non-positive integer
                let (s, f1) = two_sum(a1, -b);
                let (a2_rounded, f2) = two_sum(s, 1.0);
                if is_nonpos_int(a2_rounded) && a2_rounded != a2 {
                    return None;
                }
                // The rounding errors of b and x change y by about (y_t - y) / T, with b and x
                // moved by T times them (far enough to stand out from the error of U, and near
                // enough to be linear). That of chgu's `a1 - b + 1` adds to that of b, as chgu
                // then returns about U at `b = a1 - a2_rounded + 1`. They may take half of
                // PFQ_TOL, as y can depend much more on b and x than on a_i and z, e.g. for an a2
                // near a nonpositive integer. A move that would round away (e.g. that of b for a
                // tiny a1) is a spacing instead, a larger step.
                let moved = |v: f64, d: f64| {
                    let v_t = v + T * d;
                    let spacing = (f64::EPSILON * v.abs()).max(f64::MIN_POSITIVE);
                    if v_t == v && d != 0.0 {
                        v + spacing.copysign(d)
                    } else {
                        v_t
                    }
                };
                let db = ((e1 + e2).abs() + (f1 + f2).abs()).copysign(e1 + e2);
                let ((p, u, digits), (p_t, u_t, digits_t)) =
                    (y(b, x)?, y(moved(b, db), moved(x, dx))?);
                // (relative to y, also where it overflows)
                let rounding = ((p_t / p) * u_t - u).abs() / u.abs() / T;
                let (y, digits) = (p * u, digits.min(digits_t));
                // (with the rounding doubled, as y can curve within the move, e.g. for an a1 near a
                // non-positive integer)
                let err = 10.0_f64.powf(-digits) + 2.0 * rounding;
                // (if accurate, U has the digits of `PFQ_TOL` plus the margin, and the rounding
                // moves y by at most `PFQ_TOL / 2`; an infinite y overflows by more than `err`)
                let accurate = digits >= -PFQ_TOL.log10()
                    && if y.is_infinite() {
                        (p / (1.0 + err) * u).is_infinite()
                    } else {
                        rounding <= PFQ_TOL / 2.0
                    };
                let err = if accurate {
                    PFQ_TOL
                } else {
                    // Within a few spacings of a non-positive integer n (also 0), the rounding of
                    // chgu's a2 can move a cancellation in U, so only an accurate y is kept there
                    // (unless U is the polynomial of a2 = n).
                    let n = a2_rounded.round().min(0.0);
                    let spacing = f64::EPSILON * a1.abs().max(b.abs()).max(1.0);
                    let near =
                        (a2 - n).abs().min((a2_rounded - n).abs()) <= POLE_SPACINGS * spacing;
                    if near && (a2 != n || a2_rounded != n) {
                        return None;
                    }
                    // (above `PFQ_TOL`, or NaN for no estimate, as for an infinite y)
                    if y.is_infinite() {
                        f64::NAN
                    } else {
                        err.max(PFQ_INACCURATE)
                    }
                };
                // (with the spacing of a subnormal y)
                Some(Estimate::rounded(y, err))
            };
            let terms = (0..).map(f64::from).scan(1.0, |t, n| {
                let t_n = *t;
                *t *= (a1 + n) * (a2 + n) * z / (n + 1.0);
                Some(t_n)
            });
            let settled = (-a1).max(-a2).max(0.0).ceil();
            // (with x capped, as sqrt(2 pi x) in the error estimate could overflow; it stays above
            // 1e8, while past `settled` the terms alternate, so the first omitted one bounds the
            // remainder)
            let scale = x.min(f64::EPSILON.recip());
            let series = || {
                let sum = asymptotic_sum(terms.clone(), settled, scale, |_| 0.0);
                Some(Estimate::absolute(sum))
            };
            // Try Weniger when both U evaluations and the truncated series miss the tolerance
            let weniger = |weniger: (f64, f64), eps: f64| {
                Some(weniger_drift(weniger, terms.clone(), settled, eps))
            };
            // then with a1 and a2 swapped (DLMF 13.2.40), as chgu may take another method, e.g.
            // chgubi instead of chguit for an integer a1 - a2; then Weniger in double-double, as
            // the rounding errors in f64 can dominate, where the terms before `settled` cancel
            // or the partial sums grow (z < -1); then Miller's algorithm, which is slower, but
            // holds where U lacks a digit by the estimate of chgu; then U for a >= 1 where chgu
            // needs the digits of PFQ_TOL of its first methods, so that it goes on to chguit (not
            // for a < 1, where its other methods are slow and seldom gain the digit); then
            // Miller's algorithm from where all a_i + k > 0, as its sum cancels before; then U
            // with the digits of PFQ_TOL; and last the residues of `pfq_borel` in double-double:
            // the two series of U (DLMF 13.2.42), which lose digits in chgu where they cancel or
            // near a pole of a Gamma function
            let hypu_tol = |a, b, x| {
                xsf_hypu_min_digits(a, b, x, (CHGU_DIGITS_MARGIN - PFQ_TOL.log10()) as i32)
            };
            let miller = |shift| hyp2f0_miller(a1, a2, z, shift).map(Estimate::absolute);
            borel(a1, a2, xsf_hypu_digits, CHGU_DIGITS_MARGIN)
                .or_better(series)
                .or_better(|| borel(a2, a1, xsf_hypu_digits, CHGU_DIGITS_MARGIN))
                .or_better(|| {
                    let (y, err) = pfq_weniger::<f64>(&[a1, a2], &[1.0], z)?;
                    weniger((y, HYP2F0_WENIGER_F64_MARGIN * err), f64::EPSILON)
                })
                .or_better(|| weniger(pfq_weniger::<Dd>(&[a1, a2], &[1.0], z)?, Dd::EPSILON))
                .or_better(|| miller(0.0))
                .or_better(|| (a1 >= 1.0).then(|| borel(a1, a2, hypu_tol, CHGU_DIGITS_MARGIN))?)
                .or_better(|| (a2 >= 1.0).then(|| borel(a2, a1, hypu_tol, CHGU_DIGITS_MARGIN))?)
                .or_better(|| {
                    // (past a shifted a_i = 0)
                    let shift = (-a1).max(-a2).floor() + 1.0;
                    (shift > 0.0).then(|| miller(shift))?
                })
                .or_better(|| borel(a1, a2, xsf_hypu_digits, 0.0))
                .or_better(|| borel(a2, a1, xsf_hypu_digits, 0.0))
                .or_better(|| pfq_borel_dd(&[a1, a2], &[1.0], z))
        },
        (&[a], &[b1, b2]) => hyp1f2_estimate(a, b1, b2, z),
        (&[a1, a2], &[b1, b2]) => hyp2f2_estimate(a1, a2, b1, b2, z),
        (&[a1, a2, a3], &[b1, b2]) => hyp3f2_estimate(a1, a2, a3, b1, b2, z),
        _ if a.len() == b.len() + 1 => {
            b.push(1.0);
            pfq_real(a, b, z)
        },
        _ if a.len() <= b.len() => {
            b.push(1.0);
            let series = || Some(pfq_series_sorted(a.clone(), b.clone(), z));
            let asymptotic = || pfq_asymptotic(&a, &b, z);
            // the terms of the series for z < 0 reach about e^x, and the result
            // e^(x cos(pi / kappa)), or a power of z for p = q
            let kappa = (b.len() - a.len()) as f64;
            let x = kappa * z.abs().powf(1.0 / kappa);
            let cancellation = x * (1.0 - cospi(1.0 / kappa).max(0.0));
            // last, as it is the slowest
            let asymptotic_dd = || pfq_asymptotic_dd(&a, &b, z);
            if degree.is_finite() {
                series()
            } else if z > 0.0 || cancellation < PFQ_ASYMPTOTIC_FIRST {
                series().or_better(asymptotic).or_better(asymptotic_dd)
            } else {
                asymptotic().or_better(series).or_better(asymptotic_dd)
            }
        },
        _ if z < 0.0 => {
            b.push(1.0);
            pfq_borel(&a, &b, z)
        },
        _ => None,
    };
    polynomial.or_better(|| y)
}

/// Generalized hypergeometric function $_pF_q$ as `hyp![a_1, ..., a_p; b_1, ..., b_q; z]`
///
/// Shorthand for [`hyppfq`], where either list may be empty.
///
/// # Examples
///
/// ```
/// use xsf::{hyp, hyppfq};
///
/// let y = hyp![1.0, 2.0; 3.0; -0.5];
/// assert_eq!(y, hyppfq(&[1.0, 2.0], &[3.0], -0.5));
/// // empty lists are allowed
/// assert_eq!(hyp![; ; 0.0], Some(1.0));
/// ```
#[macro_export]
macro_rules! hyp {
    ($($a:expr),* $(,)? ; $($b:expr),* $(,)? ; $z:expr $(,)?) => {
        $crate::hyppfq(&[$($a),*], &[$($b),*], $z)
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_asymptotic_sum_cap_after_hump() {
        // small terms, a hump before term `settled + |z|`, then terms that shrink too slowly to
        // end before the cap
        let term = |n: u32| {
            let t = match n {
                0..5 => 10.0_f64.powi(-i32::try_from(n).unwrap()),
                5..10 => 10.0_f64.powi(i32::try_from(n).unwrap() - 7),
                _ => 1e3 * 0.9999_f64.powi(i32::try_from(n).unwrap() - 10),
            };
            (t, t)
        };
        let (_, err) = asymptotic_sum((0..).map(term), 0.0, 20.0, |_| 0.0);
        assert_eq!(err, f64::INFINITY);
    }

    #[test]
    fn test_estimate_absolute_at_tol() {
        // `err / y` rounds to the float above `PFQ_TOL`
        let y = 0.012_571_498_510_638_38;
        assert!(Estimate::absolute((y, PFQ_TOL * y)).is_accurate());
    }

    #[test]
    fn test_pfq_sum_rescales_repeatedly() {
        // the partial sum jumps from 1 to 1e300, whose square needs two rescalings
        let (y, err, tails, _) =
            pfq_sum_terms::<f64>(&[-1.0], &[1.0], 1.0, 2, [1e300].into_iter(), 1e-13);
        assert_eq!(y, 1e300);
        assert!(err.is_finite() && tails.is_finite());
    }

    #[test]
    fn test_hyp2f0_miller_err() {
        // the estimate covers the rounding of the result, here 2.5e-25 (Arb)
        let (y, err) = hyp2f0_miller(-36.5, -51.8, -0.026, 0.0).unwrap();
        let expected = -6.120_387_741_382_293e-9;
        assert!(err >= f64::EPSILON * y.abs() && (y - expected).abs() <= err);
        // and that of the recurrence in a1 and a2, here 3.6e-13 (Arb)
        let (y, err) = hyp2f0_miller(
            -50.999_999_999_699_924,
            -104.536_535_990_336_39,
            -0.056_613_881_920_597_09,
            105.0,
        )
        .unwrap();
        let expected = 2.794_502_416_895_517_4e22;
        assert!((y - expected).abs() <= err);
    }

    #[test]
    fn test_root_dd() {
        // of a double-double with a low part
        let v = Dd::new(10.0) / 3.0;
        for k in 1..=5 {
            let root = root_dd(v, k);
            let pow = (0..k).fold(Dd::ONE, |p, _| p * root);
            assert!((pow - v).hi.abs() <= pow2(-100) * v.hi, "k = {k}");
        }
    }

    #[test]
    fn test_mul_gammas() {
        // Gamma(g) Gamma(180.5) / (Gamma(178.25) Gamma(1.5)) for g = -3 + 2^-40 + lo, where the
        // large arguments are by `ln_gamma`, and where (lo 2^40)^2 is negligible or not (mpmath)
        let one = Scaled::new(1.0, 0, 0.0).unwrap();
        let (num, den) = ([Dd::new(180.5)], [178.25, 1.5].map(Dd::new));
        let check = |y: Option<Scaled>, expected: f64, tol: f64| {
            let (y, err) = y.unwrap().to_f64();
            assert!(err <= tol && (y / expected - 1.0).abs() <= err, "{y} != {expected}");
        };
        for (e, expected, tol) in
            [(-75, -2.419_553_175_200_6e16, 1e-14), (-60, -2.419_550_867_807_498_4e16, 1e-12)]
        {
            let g = Dd::new(-3.0 + pow2(-40)) + pow2(e);
            check(mul_gammas(num.into_iter().chain([g]), den, one), expected, tol);
        }
        // equal factors, whose errors are equal too
        let g = Dd::new(-8.000_005_264_759_276);
        check(mul_gammas([g; 3], [], one), -104.541_370_932_890_13, 1e-14);
        // the mantissa of `y` times `Gamma(g)` overflows
        let (g, big) = (Dd::new(-3.0 + pow2(-40)), Scaled::new(1e300, 0, 0.0));
        check(mul_gammas([g], [g], big.unwrap()), 1e300, 1e-14);
        // `PFQ_GAMMA_MAX - g` rounds to an integer, and `(g)_n` is below the normal range
        check(
            mul_gammas([Dd::new(-8.000_000_000_000_002)], [], one),
            -13_962_052_416.203_12,
            1e-15,
        );
        check(mul_gammas([Dd::new(1e-290)], [], one), 9.999_999_999_999_999e289, 1e-15);
        check(mul_gammas([], [Dd::new(-1e-290)], one), -1e-290, 1e-15);
        // a pole
        assert!(mul_gammas([Dd::new(-2.0)], den, one).is_none());
        assert!(mul_gammas(num, [Dd::new(-2.0)], one).is_none());
        // in double-double, below `-PFQ_GAMMA_MAX`: Gamma(-30.5) / Gamma(-28.5)
        let y = mul_gammas_in([Dd::new(-30.5)], [Dd::new(-28.5)], one, true);
        check(y, 1.0 / (30.5 * 29.5), 1e-15);
    }

    #[test]
    fn test_pfq_borel_dd() {
        // the Borel sum for p > q + 1, whose series are in -1 / z (2F0; Arb) or in 1 / z (3F0;
        // mpmath meijerg)
        for (a, z, expected) in [
            (vec![0.3, 1.7], -2.5, 0.651_478_922_494_888_7),
            (vec![0.3, 1.7, 0.6], -2.5, 0.790_728_400_803_319_9),
        ] {
            let y = pfq_borel_dd(&a, &[1.0], z).unwrap();
            assert!(y.err < 8.0 * f64::EPSILON, "{y:?}");
            assert!((y.y / expected - 1.0).abs() <= y.err, "{y:?}");
        }
    }

    #[test]
    fn test_scaled_certain() {
        // m (1 + err) 2^1024 is above the overflow threshold, unlike with `1 + err` rounded,
        // for the err = 3.25390625 eps after `fold_exp`
        let (m, err) = (1.0 - 3.5 * f64::EPSILON, 1.253_906_25 * f64::EPSILON);
        assert!(Estimate::scaled_certain(Scaled::new(m, 1024, err).unwrap()).is_none());
        // (finite at both ends)
        let y = Estimate::scaled_certain(Scaled::new(m, 1024, 0.0).unwrap()).unwrap();
        assert!(y.y < f64::MAX, "{y:?}");
        // an exponent that `fold_exp` clamps overflows at both ends
        let l = Dd::new(2.5e12);
        let y = Estimate::scaled_certain(Scaled::new(-1.0, 0, 0.0).unwrap().with_l(l)).unwrap();
        assert!(y.y == f64::NEG_INFINITY && y.is_accurate(), "{y:?}");
    }

    #[test]
    fn test_add_exp_dd() {
        // 1 + 2^1023 e^-710, where `exp` underflows before it is times the mantissa (mpmath)
        let one = Scaled::new(1.0, 0, 0.0).unwrap();
        let big = Scaled::new(pow2(1023), 0, 0.0)
            .unwrap()
            .with_l(Dd::new(-710.0));
        let (y, err) = add_exp_dd([one, big]).unwrap().to_f64();
        assert!(err <= 1e-15 && (y / 1.402_349_450_878_746_4 - 1.0).abs() <= err, "{y}, {err}");
    }

    #[test]
    fn test_ln_gamma_stirling() {
        // where its error is the largest: Gamma(13) = 12!
        let (v, sign, err) = ln_gamma_stirling(Dd::new(PFQ_GAMMA_MAX));
        assert!(sign == 1.0 && err <= 1e-28);
        assert!((v - ln_dd(479_001_600.0)).hi.abs() <= err);
    }

    #[test]
    fn test_pfq_alg_series_dd_growing() {
        // the terms of the series of a_1 grow from the fourth one on, long before `n = |z|`
        let (a, b) = ([pow2(-70), -19.75, -29.625], [1.0, -219.75, -229.625, 1.0]);
        assert!(pfq_alg_series_dd(&a, &b, 0, -100.0).1.is_infinite());
        // the low parts of the parameters make the term after the one where they settle larger
        let (s, u) = (pow2(-49), pow2(-52));
        let (a, b) = ([s - 16.0, 1.0 + 3.0 * u], [1.0 + u, 2.0 * s - 16.0, 1.0]);
        assert!(pfq_alg_series_dd(&a, &b, 0, -1.1 * s).1.is_infinite());
        // the bound of the unpaired factor underflows, and the terms grow to ~1e72 from the fifth
        let c = (1..=10).map(|i| f64::from(i) / 16.0);
        let a: Vec<_> = once(f64::from_bits(1)).chain(c.clone()).collect();
        let b: Vec<_> = once(3.0)
            .chain(c.map(|c| c - pow2(49)))
            .chain([1.0])
            .collect();
        assert!(pfq_alg_series_dd(&a, &b, 0, -2.0).1.is_infinite());
    }

    #[test]
    fn test_pfq_alg_series_dd_polynomial() {
        // of degree 0, which is below `|z|`
        let (sum, err) = pfq_alg_series_dd(&[0.5, 0.75], &[1.5, 0.25, 1.0], 0, -0.5);
        assert!(sum.hi == 1.0 && err <= Dd::EPSILON);
    }

    #[test]
    fn test_pfq_asymptotic_integer_difference() {
        // the bound on the algebraic part is too large for equal a_i, and near them the
        // expansion is accurate or `None` (Arb reference)
        let c = 1.0 + pow2(-10);
        let (b, z) = ([c, c, 11.25, 1.0], 21.25);
        let accurate = |a: &[f64; 3]| pfq_asymptotic(a, &b, z).and_then(Estimate::accurate);
        assert_eq!(accurate(&[1.0, 1.0, 10.5]), None);
        let near = accurate(&[1.0, 1.0 + pow2(-20), 10.5]);
        let expected = 739_413_030.005_544_4;
        assert!(near.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected));
    }

    #[test]
    fn test_pfq_asymptotic_nearly_cancelled() {
        // equal a_i within sqrt(eps) of b_2 and b_3, where the algebraic part is about 1 and
        // changes much with a_i: accurate or `None` (Arb reference)
        let t = -pow2(-40);
        let (root, root_2) = (f64::EPSILON.sqrt(), (2.0 * f64::EPSILON).sqrt());
        let y = pfq_asymptotic(&[t; 3], &[t - 1.0, t + root, t + root_2, 1.0], 30.0)
            .and_then(Estimate::accurate);
        let expected = -816_436.169_307_593_7_f64;
        assert!(y.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected.abs()));
    }
}
