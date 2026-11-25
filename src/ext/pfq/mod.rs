//! Generalized hypergeometric functions pFq of real arguments.
//!
//! The internal functions take `a` and `b` directly (without a synthetic `1` for `n!`), and
//! supply the `(1)_n = n!` factor in the denominator where series terms or ratio bounds are
//! evaluated.

use core::cell::OnceCell;
use core::f64::consts::{FRAC_2_PI, FRAC_PI_2};
use core::iter::zip;
use core::ops::Add;

use super::quadpack::quad;
use super::{Scaled, add_exp, add_exp_in, fold_exp, ln_dd, ln_gamma};
use crate::utils::{EPSILON_SQ, is_nonpos_int, neg_one_pow, next_up_nonneg, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{
    DoubleDouble as Dd, FRAC_PI_2_DD, frexp_f64, ldexp_f64, two_sum,
};
use crate::xsf::cephes::gammasgn;
use crate::xsf::cephes::hyp2f1::HYP2F1_LOST;
use crate::xsf::cospi;

mod asymptotic;
mod borel;
mod conformal;
mod gamma;
mod hyp3f2;
mod hyp3f2_endpoint;
pub use hyp3f2::hyp3f2;
mod recurrence;
mod residue;
mod series;

/// Maximum number of terms of `pfq_series`
const PFQ_MAXITER: u32 = 10_000;

/// Maximum number of terms of `pfq_series` in double-double arithmetic, fewer as a term takes
/// about 20 times as long
const PFQ_DD_MAXITER: u32 = PFQ_MAXITER / 4;

/// Default term limit for `pfq_series_long`
///
/// Of 199 `None` of 3F2 on 20000 points, it gives 13, and twice as many terms would give 4
/// more, in up to twice the time (measured).
const PFQ_LONG_MAXITER: u32 = 2 * PFQ_MAXITER;

/// Largest estimated relative error of an accurate result, after which no other method is tried
const PFQ_TOL: f64 = 1e-13;

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
        match self {
            Some(y) if y.is_accurate() => self,
            Some(y) => other().filter(|v| v.err < y.err).or(self),
            None => other(),
        }
    }
}

/// The sizes of the parameters for which `pfq_products` is exact enough for all `n`
const PFQ_MODERATE: core::ops::RangeInclusive<f64> = 1e-30..=1e30;

/// Most effective factors in `a` or in `b` (including the 1 of `n!`) for which `pfq_products` is
/// exact
const PFQ_MAX_UNCHECKED: usize = 4;

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
    /// Adds `term` to `(sum, comp)` with Neumaier compensation in `f64` (no-op `comp` in `Dd`)
    fn add_compensated(sum: Self, comp: Self, term: Self) -> (Self, Self);
    /// The compensated partial sum and its offset from `t_0 = 1`, as `f64`
    fn sum_and_delta(sum: Self, comp: Self) -> (f64, f64);
    /// Upper bound on `|self|` as `f64`
    fn abs_up(self) -> f64;
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

    #[inline]
    fn add_compensated(sum: Self, comp: Self, term: Self) -> (Self, Self) {
        let next = sum + term;
        if next.is_finite() && term.is_finite() {
            let delta = if sum.abs() >= term.abs() {
                (sum - next) + term
            } else {
                (term - next) + sum
            };
            let next_comp = comp + delta;
            (
                next,
                if next_comp.is_finite() {
                    next_comp
                } else {
                    comp
                },
            )
        } else {
            (next, 0.0)
        }
    }

    #[inline]
    fn sum_and_delta(sum: Self, comp: Self) -> (f64, f64) {
        (sum + comp, (sum - 1.0) + comp)
    }

    #[inline]
    fn abs_up(self) -> f64 {
        self.abs()
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

    #[inline]
    fn add_compensated(sum: Self, _comp: Self, term: Self) -> (Self, Self) {
        (sum + term, Self::ZERO)
    }

    #[inline]
    fn sum_and_delta(sum: Self, _comp: Self) -> (f64, f64) {
        (sum.hi, (sum + (-1.0)).hi)
    }

    #[inline]
    fn abs_up(self) -> f64 {
        dd_hi_up(self.abs())
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

/// The smallest non-negative integer `n` with `v + n = 0` for a `v` in `params`, or infinity
fn first_zero_dd(params: &[Dd]) -> f64 {
    params
        .iter()
        .filter(|x| is_pole(x))
        .map(|x| -x.hi)
        .fold(f64::INFINITY, f64::min)
}

/// The smallest integer `n` with `x + n > 0`
fn positive_from(x: &Dd) -> f64 {
    (-*x).floor().hi + 1.0
}

/// Upper bound on a non-negative `Dd` as `f64`
fn dd_hi_up(x: Dd) -> f64 {
    if x.lo > 0.0 && x.hi >= 0.0 {
        next_up_nonneg(x.hi)
    } else {
        x.hi
    }
}

/// Lower bound on a positive `Dd` as `f64`
fn dd_hi_down(x: Dd) -> f64 {
    if x.lo < 0.0 && x.hi > 0.0 && x.hi.is_finite() {
        f64::from_bits(x.hi.to_bits() - 1)
    } else {
        x.hi
    }
}

/// Sorts `v` in descending order, comparing low words on ties
fn sort_descending_dd(v: &mut [Dd]) {
    v.sort_unstable_by(|x, y| y.hi.total_cmp(&x.hi).then_with(|| y.lo.total_cmp(&x.lo)));
}

/// The `n` from which on all `a_i + n` and `b_j + n` are positive
fn pfq_settled(a: &[f64], b: &[f64]) -> f64 {
    a.iter().chain(b).fold(0.0, |n, &v| n.max(-v))
}

/// Sorts `v` in descending order, as `pfq_ratio_bound` is tightest if `a` and `b` are sorted alike
fn sort_descending(v: &mut [f64]) {
    v.sort_unstable_by(|x, y| y.total_cmp(x));
}

/// `v` in descending order
fn descending<P: AsMut<[f64]>>(mut v: P) -> P {
    sort_descending(v.as_mut());
    v
}

/// Upper bound on `1 / |z|`
fn recip_up(z: f64) -> f64 {
    let az = z.abs();
    if az.is_infinite() {
        0.0
    } else {
        next_up_nonneg(1.0 / az)
    }
}

/// Small stack-backed scratch buffer (up to 8 elements before spilling to `Vec`)
struct Scratch<T> {
    buf: [T; 8],
    vec: Vec<T>,
}

impl<T: Copy> Scratch<T> {
    fn new(len: usize, fill: T) -> Self {
        Self {
            buf: [fill; 8],
            vec: if len > 8 { vec![fill; len] } else { Vec::new() },
        }
    }

    fn fill(&mut self, len: usize, iter: impl IntoIterator<Item = T>) -> &mut [T] {
        let dst = if len <= self.buf.len() {
            &mut self.buf[..len]
        } else {
            self.vec.resize(len, self.buf[0]);
            &mut self.vec[..]
        };
        for (d, x) in zip(&mut *dst, iter) {
            *d = x;
        }
        dst
    }
}

/// Whether `v` is normal in the arithmetic of `T`
fn is_normal<T: SeriesReal>(v: T) -> bool {
    (T::MIN_POSITIVE..=f64::MAX).contains(&v.to_f64().abs())
}

/// `m` times `2^e`, where beyond `i32` it over- or underflows anyway
fn ldexp_i64<T: SeriesReal>(m: T, e: i64) -> T {
    m.ldexp(e.clamp(i32::MIN.into(), i32::MAX.into()) as i32)
}

/// `v` without `v[i]`
fn skip(v: &[f64], i: usize) -> impl Iterator<Item = f64> + Clone + '_ {
    (0..v.len()).filter(move |&k| k != i).map(move |k| v[k])
}

/// A method for `p = q + 1`, returning an `Estimate` or `None`
type Method = fn(&[f64], &[f64], f64) -> Option<Estimate>;

/// `_{q+1}F_q` for `z < 1`, or for a polynomial, by the first method for `z` that is accurate,
/// else by the most accurate one
///
/// Unless it is a polynomial, the last methods are `pfq_conformal`, which is the slowest, and for
/// `|z| < 1` `pfq_recurrence`, `pfq_series_long`, and `pfq_recurrence_extended` (3F2 only for the
/// recurrences), and then `pfq_overflow`.
fn pfq_accurate(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let series: Method = |a, b, z| Some(series::pfq_series(a, b, z));
    let weniger: Method = |a, b, z| series::pfq_weniger::<f64>(a, b, z).map(Estimate::absolute);
    let polynomial = first_zero(a).is_finite();
    let methods: &[Method] = if polynomial {
        &[series, hyp3f2::legendre_square]
    } else if z > PFQ_WENIGER_MAX_Z {
        &[series]
    } else if z > -1.0 {
        &[weniger, series]
    } else if z == -1.0 {
        &[weniger]
    } else if z >= PFQ_WENIGER_FIRST_Z {
        &[weniger, asymptotic::pfq_continuation]
    } else {
        &[asymptotic::pfq_continuation, weniger]
    };
    let last: Option<Method> = (!polynomial).then_some(conformal::pfq_conformal);
    let inside: &[Method] = if !polynomial && z.abs() < 1.0 {
        &[recurrence::pfq_recurrence, series::pfq_series_long, recurrence::pfq_recurrence_extended]
    } else {
        &[]
    };
    let overflow: Option<Method> = (!polynomial).then_some(series::pfq_overflow);
    let methods = methods.iter().chain(&last).chain(inside).chain(&overflow);
    methods.fold(None, |y, method| y.or_better(|| method(a, b, z)))
}

/// The sum of `v` in double-double
fn sum_dd(v: &[f64]) -> Dd {
    v.iter().fold(Dd::ZERO, |s, &v| s + v)
}

/// Whether the series of `_{q+1}F_q` diverges at `z = 1`, where `sum_j b_j <= sum_i a_i`, or
/// `None` if that is uncertain
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

/// `_{q+1}F_q` for finite arguments (sorts `a` and `b` in place when evaluating)
fn pfq_real(a: &mut [f64], b: &mut [f64], z: f64) -> Option<Estimate> {
    if z == 0.0 {
        return Some(Estimate::exact(1.0));
    }
    let Some(degree) = pfq_degree(a, b) else {
        return Some(Estimate::exact(f64::INFINITY));
    };
    if degree.is_infinite() && z > 1.0 {
        // complex
        Some(Estimate::exact(f64::NAN))
    } else if degree.is_infinite() && z == 1.0 {
        if let Some(value) = hyp3f2::unit_upper_at_one(a, b) {
            return Some(value);
        }
        // +-inf where the series diverges, and no result where it converges without a method
        if diverges_at_one(a, b)? {
            return Some(Estimate::exact(f64::INFINITY * late_sign(a, b)));
        }
        let (&[a1, a2, a3], &[b1, b2]) = (&*a, &*b) else {
            return None;
        };
        hyp3f2_endpoint::estimate([a1, a2, a3], [b1, b2])
    } else {
        sort_descending(a);
        sort_descending(b);
        let result = pfq_accurate(a, b, z);
        if z == 1.0 && degree <= f64::from(PFQ_DD_MAXITER) {
            result.or_better(|| hyp3f2::saalschutz(a, b, degree))
        } else {
            result
        }
    }
}

/// `pFq` for `p <= q` by series and asymptotic expansion, trying `series` first when `series_first`
fn pfq_entire(a: &[f64], b: &[f64], z: f64, series_first: bool) -> Option<Estimate> {
    pfq_entire_with_exp(a, b, z, series_first, None)
}

fn pfq_entire_with_exp(
    a: &[f64],
    b: &[f64],
    z: f64,
    series_first: bool,
    exp_cache: Option<&OnceCell<Option<Scaled>>>,
) -> Option<Estimate> {
    let Some(degree) = pfq_degree(a, b) else {
        return Some(Estimate::exact(f64::INFINITY));
    };
    let series = || Some(series::pfq_series_sorted(a, b, z));
    if degree.is_finite() {
        return series();
    }
    let asymptotic = || asymptotic::pfq_asymptotic(a, b, z, exp_cache);
    let result = if series_first {
        series().or_better(asymptotic)
    } else {
        asymptotic().or_better(series)
    };
    // last, as it is the slowest (and applies to `p = q`, `z <= 0` only)
    result.or_better(|| asymptotic::pfq_asymptotic_dd(a, b, z))
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
/// - [`hyp3f2`](crate::hyp3f2): Generalized hypergeometric function,
///   $\hyp 3 2 {a_1\enspace a_2\enspace a_3} {b_1\enspace b_2} z$
#[must_use]
pub fn hyp2f2(a1: f64, a2: f64, b1: f64, b2: f64, z: f64) -> f64 {
    Estimate::or_nan(hyp2f2_estimate(a1, a2, b1, b2, z))
}

/// `hyp2f2` and an estimate of its relative error
fn hyp2f2_estimate(a1: f64, a2: f64, b1: f64, b2: f64, z: f64) -> Option<Estimate> {
    if ![a1, a2, b1, b2, z].into_iter().all(f64::is_finite) {
        pfq_estimate(&[a1, a2], &[b1, b2], z)
    } else if z == 0.0 {
        Some(Estimate::exact(1.0))
    } else {
        let (a, b) = ([a1, a2], [b1, b2]);
        let series_first = z.abs() < HYP2F2_ASYMPTOTIC_FIRST;
        // a common parameter cancels, also at a pole; stay on `pfq_series`/`pfq_asymptotic`
        // rather than `hyp1f1` to avoid its 500-term polynomial cap and Kummer exp(-750) underflow
        let pairs = [(0, 0), (0, 1), (1, 0), (1, 1)];
        match pairs.into_iter().find(|&(i, j)| a[i] == b[j]) {
            Some((i, j)) => {
                let (a0, b0) = (a[1 - i], b[1 - j]);
                if a0 == b0 {
                    // (with the rounding of `exp`)
                    Some(Estimate::rounded(z.exp(), f64::EPSILON))
                } else {
                    pfq_entire(&[a0], &[b0], z, series_first)
                }
            },
            None => pfq_entire(&a, &b, z, series_first),
        }
    }
}

/// `v^(1 / k)` in double-double precision, for `v > 0`
fn root_dd(v: impl Into<Dd>, k: usize) -> Dd {
    let v = v.into();
    if k == 1 || v.hi == 0.0 || !v.hi.is_finite() {
        return v;
    }
    let (m, e) = v.frexp();
    let ki = i32::try_from(k).unwrap_or(i32::MAX);
    // scale `v` by `2^(-q * k)` so `u` stays normal: `u in [0.5, 2^(k - 1))` for `k <= 512`, and
    // symmetric `u in [2^-513, 2^512]` for `513..=1024` (no subnormal loss, `pow(k)` overflow, or
    // `|ln v|` amplification in Newton's step)
    let (u, q) = if k <= 512 {
        let q = e.div_euclid(ki);
        (m.ldexp(e - q * ki), q)
    } else if k <= 1024 {
        let q = (e + ki / 2).div_euclid(ki);
        (m.ldexp(e - q * ki), q)
    } else {
        (v, 0)
    };
    let r = if k == 2 {
        u.sqrt()
    } else {
        let s = u.hi.powf(1.0 / k as f64);
        let prev = (1..k).fold(Dd::ONE, |p, _| p * s);
        // (one Newton step)
        Dd::new(s) + (u - prev * s) / (k as f64 * prev.hi)
    };
    r.ldexp(q)
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
        pfq_estimate(&[a], &[b1, b2], z)
    } else if z == 0.0 {
        Some(Estimate::exact(1.0))
    } else {
        let b = [b1, b2];
        let series_first = z > -HYP1F2_ASYMPTOTIC_FIRST;
        // a common parameter cancels, also at a pole; keep on `pfq_series`/`pfq_asymptotic`
        // rather than `hyp0f1` so large negative `z` stays accurate
        if let Some(j) = b.iter().position(|&bj| bj == a) {
            pfq_entire(&[], &[b[1 - j]], z, series_first)
        } else {
            let exp = OnceCell::new();
            pfq_entire_with_exp(&[a], &b, z, series_first, Some(&exp))
                .or_better(|| {
                    // DLMF 16.12.1, for the exact parameters in either order.
                    let (b, low) = two_sum(a, 0.5);
                    let pair = (b, 2.0 * a);
                    if !series_first
                        && a > 0.0
                        && low == 0.0
                        && (pair == (b1, b2) || pair == (b2, b1))
                    {
                        let root = pfq_entire(&[], &[b], z * 0.25, false)?;
                        let y = root.y * root.y;
                        // Propagate the relative error of the factor and the rounding of its square.
                        let err = (2.0 + root.err) * root.err + f64::EPSILON;
                        (err <= PFQ_TOL && y.is_finite() && y.mul_add(-err, y) > f64::MIN_POSITIVE)
                            .then(|| Estimate::new(y, err))
                    } else {
                        None
                    }
                })
                .or_better(|| {
                    asymptotic::hyp1f2_asymptotic_dd(a, &b, z, &exp)
                        .filter(|estimate| estimate.is_accurate())
                })
        }
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

/// `add_exp` in double-double
fn add_exp_dd(terms: impl IntoIterator<Item = Scaled, IntoIter: Clone>) -> Option<Scaled> {
    add_exp_in(terms, Dd::exp, Dd::EPSILON)
}

/// The limit at `z = +-inf` for finite parameters if cheap
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
/// | 3   |                           |                           | [`hyp3f2`](crate::hyp3f2) |
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
/// - $\pm\infty$ for $p = q + 1 \ge 3$ and $0 < z < 1$ where it detects that it overflows;
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
/// - for the other $p \le q + 1$ where it can't be computed accurately, as for [`hyp3f2`](crate::hyp3f2) if
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
pub fn hyppfq(a: impl AsRef<[f64]>, b: impl AsRef<[f64]>, z: f64) -> Option<f64> {
    pfq_estimate(a.as_ref(), b.as_ref(), z)?.accurate()
}

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
    let polynomial = degree
        .is_finite()
        .then(|| Estimate::from(series::pfq_sum::<Dd>(&a, &b, z, PFQ_DD_MAXITER)));
    if polynomial.is_some_and(Estimate::is_accurate) {
        return polynomial;
    }
    if z == 1.0 && degree.is_finite() {
        let arrays = a
            .as_slice()
            .try_into()
            .ok()
            .zip(b.as_slice().try_into().ok());
        if let Some(estimate) = arrays.and_then(|(a, b)| hyp4f3_racah(a, b)) {
            return Some(estimate);
        }
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
        (&[a0], &[b0]) => xsf(crate::hyp1f1(a0, b0, z)).or_better(|| {
            if a0 > 0.0 && b0 > a0 && z.is_finite() && -z >= HYP2F2_ASYMPTOTIC_FIRST {
                // For b-a >= 1, the beta integral gives M(a,b,-t) <= (b/t)^a
                // (DLMF 13.4.1 and Jensen's inequality for the gamma ratio).
                // Use only exponents to skip values below half the smallest subnormal.
                let (_, b_exp) = frexp_f64(b0);
                let (_, t_exp) = frexp_f64(-z);
                if b0 - a0 > 1.0 && a0 * f64::from(t_exp - b_exp - 1) > 1075.0 {
                    return None;
                }
                // As in 2F2 with a canceled parameter, try the large-argument expansion.
                asymptotic::pfq_asymptotic(&a, &b, z, None)
                    .filter(|v| v.y.is_finite() && v.is_accurate())
            } else {
                None
            }
        }),
        (&[a1, a2], &[b]) => xsf(crate::hyp2f1(a1, a2, b, z)),
        _ if z.is_infinite() => None,
        (&[a1, a2], &[]) if z < 0.0 => borel::hyp2f0_borel(a1, a2, z),
        (&[a], &[b1, b2]) => hyp1f2_estimate(a, b1, b2, z),
        (&[a1, a2], &[b1, b2]) => hyp2f2_estimate(a1, a2, b1, b2, z),
        (&[a1, a2, a3], &[b1, b2]) => hyp3f2::estimate(a1, a2, a3, b1, b2, z),
        // Near one the elliptic integral is faster than summing this series.
        (&[0.5, 0.5, 0.5, 0.5], &[1.0, 1.0, 1.0]) if (0.9375..=1.0).contains(&z) => {
            Some(hyp4f3_half_integral(z)).or_better(|| pfq_real(&mut a, &mut b, z))
        },
        _ if a.len() == b.len() + 1 => pfq_real(&mut a, &mut b, z),
        _ if a.len() <= b.len() => {
            // the terms of the series for z < 0 reach about e^x, and the result
            // e^(x cos(pi / kappa)), or a power of z for p = q
            let kappa = (b.len() + 1 - a.len()) as f64;
            let x = kappa * z.abs().powf(1.0 / kappa);
            let cancellation = x * (1.0 - cospi(1.0 / kappa).max(0.0));
            pfq_entire(&a, &b, z, z > 0.0 || cancellation < PFQ_ASYMPTOTIC_FIRST).or_better(|| {
                match (a.as_slice(), b.as_slice()) {
                    (&[a1, a2], &[b1, b2, b3]) => hyp2f3_bessel_product([a1, a2], [b1, b2, b3], z),
                    (&[], &[_, _, _]) if [0.25, 0.5, 0.75].iter().all(|v| b.contains(v)) => {
                        hyp0f3_quarter_parameters(z)
                    },
                    _ => None,
                }
            })
        },
        _ if z < 0.0 => {
            borel::pfq_borel(&a, &b, z).or_better(|| borel::hyp3f0_near_coincident(&a, &b, z))
        },
        _ => None,
    };
    polynomial.or_better(|| y)
}

// Use the coefficient identity DLMF 18.26.20 when the terminating Racah series loses accuracy.
#[cold]
fn hyp4f3_racah(mut a: [f64; 4], mut b: [f64; 3]) -> Option<Estimate> {
    a.sort_by(f64::total_cmp);
    b.sort_by(f64::total_cmp);
    let [neg_n, neg_y, upper_a, upper_b] = a;
    let [neg_size, lower_a, lower_b] = b;
    let (n, y, size) = (-neg_n, -neg_y, -neg_size);
    if !([n, y, size].into_iter().all(|x| x.fract() == 0.0)
        && 1.0 <= y
        && n < size
        && size <= f64::from(PFQ_DD_MAXITER)
        && upper_a > 0.0
        && lower_a > 0.0)
    {
        return None;
    }

    // Reject any inexact parameter arithmetic. In particular, rounding a
    // nearly balanced input into the Racah family would change the function.
    let exact_add = |x, y| {
        let (sum, err) = two_sum(x, y);
        (err == 0.0 && sum.is_finite()).then_some(sum)
    };
    let evaluate = |an, ay, beta_delta_plus_one, gamma_plus_one| {
        let beta = exact_add(an, size - n)?;
        let gamma = exact_add(gamma_plus_one, -1.0)?;
        let delta = exact_add(exact_add(beta_delta_plus_one, -beta)?, -1.0)?;
        let g1 = exact_add(y, gamma_plus_one)?;
        if !(delta > -1.0 && beta > exact_add(size, gamma)? && ay == exact_add(g1, delta)?) {
            return None;
        }
        // These index differences are exact integers in the admitted range.
        let gap = size - n - y;
        let start = (n + y - size).max(0.0);
        let remaining = n - start;
        let f = [-y, exact_add(exact_add(beta, -gamma)?, -y)?];
        let g = [y - size, g1];
        let c = exact_add(-size, -delta)?;
        // Normalize by the first convolution term. Reversing g's
        // coefficient recurrence gives a terminating hypergeometric series.
        let upper =
            [start - y, exact_add(f[1], start)?, -remaining, exact_add(1.0 - remaining, -c)?];
        let lower = [
            exact_add(beta_delta_plus_one, start)?,
            1.0 + gap.abs(),
            exact_add(1.0 - remaining, -g1)?,
        ];
        let sum = series::pfq_sum::<Dd>(&upper, &lower, 1.0, PFQ_DD_MAXITER);

        // Keep full double-double precision and use a symmetric exponent range.
        let min = <Dd as SeriesReal>::MIN_POSITIVE;
        let safe = min..=min.recip();
        let scaled_rising = |a: &[f64], b: &[f64], count: f64| {
            let checked = series::pfq_extreme(a, b);
            let product = (0..count as u32).try_fold(Dd::ONE, |value, j| {
                let value = value * series::pfq_ratio::<Dd>(a, b, f64::from(j), checked);
                safe.contains(&value.hi.abs()).then_some(value)
            })?;
            Some(product.frexp())
        };
        if !sum.value.hi.is_normal() {
            return None;
        }
        // At most 2n fixed-arity ratios plus the final scaled product.
        let rel = sum.error / sum.value.hi.abs()
            + 32.0 * <Dd as SeriesReal>::EPSILON * (n + 1.0)
            + f64::EPSILON;
        let (s, se) = sum.value.frexp();
        let (f, fe) = scaled_rising(&f, &[beta_delta_plus_one], start)?;
        let (g, ge) = scaled_rising(&g, &[c], remaining)?;
        let (scale, e) = scaled_rising(&[-size, gamma_plus_one], &[c], n)?;
        let m = (s * f * g / scale).hi;
        let e = se + fe + ge - e;
        let est = Estimate::scaled_certain(Scaled::new(m, e, rel)?)?;
        est.is_accurate().then_some(est)
    };
    [(upper_a, upper_b), (upper_b, upper_a)]
        .into_iter()
        .find_map(|(an, ay)| {
            evaluate(an, ay, lower_a, lower_b).or_else(|| evaluate(an, ay, lower_b, lower_a))
        })
}

// Integrate the three-dimensional hyper-BCC elliptic square over the fourth angle
// (Guttmann, J. Phys. A 43 (2010) 305205, sections 1 and 2.1).
#[cold]
fn hyp4f3_half_integral(z: f64) -> Estimate {
    let gap = 1.0 - z;
    let f = |theta: f64| {
        let s = theta.sin();
        crate::ellipkm1(0.5 * (1.0 + (gap + z * s * s).sqrt())).powi(2)
    };
    let (value, error) = if z == 1.0 {
        quad(f, 0.0, FRAC_PI_2, 0.0, PFQ_TOL / 4.0)
    } else {
        // Resolve the layer at theta = 0 with theta = sqrt(1-z) sinh(u).
        let d = gap.sqrt();
        let hi = (FRAC_PI_2 / d).asinh();
        quad(|u| f(d * u.sinh()) * (d * u.cosh()), 0.0, hi, 0.0, PFQ_TOL / 4.0)
    };
    let norm = FRAC_2_PI.powi(3);
    // Allow for the transformed abscissae, Jacobian, elliptic evaluation and normalization
    // separately from the quadrature estimate.
    Estimate::absolute((norm * value, norm * (error + 100.0 * f64::EPSILON * value)))
}

// The defining series (DLMF 16.2.1) for these parameters is sum (256z)^k/(4k)!.
// For z < 0, grouping the four exponential series gives cosh(t) cos(t), t^4 = -64z.
#[cold]
fn hyp0f3_quarter_parameters(z: f64) -> Option<Estimate> {
    // For 1 <= t <= 1024, reduce to the nearest zero of cos using double-double π/2.
    if !(-pow2(34)..=-1.0 / 64.0).contains(&z) {
        return None;
    }
    let t = Dd::new(-z).sqrt().mul_pwr2(8.0).sqrt();
    let zero = (0.5 * (t.hi / FRAC_PI_2_DD.hi - 1.0)).round();
    let reduced = t - (2.0 * zero + 1.0) * FRAC_PI_2_DD;
    let cosine = neg_one_pow(zero + 1.0) * reduced.hi.sin();
    // Bound the absolute phase error from the two roots and the reduction,
    // then allow for rounding sin, the cosh factor and their product.
    let err = 8.0 * EPSILON_SQ * t.hi / cosine.abs() + 8.0 * f64::EPSILON;
    if err > PFQ_TOL {
        return None;
    }
    let amplitude = cosine * (1.0 + (-2.0 * t.hi).exp());
    Estimate::scaled_certain(Scaled::new(amplitude, -1, err)?.with_l(t))
}

// DLMF 16.12.1 gives 0F1(u; z/4) 0F1(v; z/4) for this family.
#[cold]
fn hyp2f3_bessel_product(a: [f64; 2], b: [f64; 3], z: f64) -> Option<Estimate> {
    let w = z * 0.25;
    if z >= 0.0 || !w.is_normal() {
        return None;
    }
    let [sum, v, u] = descending(b);
    let (u_minus_one, rem) = two_sum(u, -1.0);
    if u < 1.0
        || rem != 0.0
        || two_sum(v, u_minus_one) != (sum, 0.0)
        || a[0].min(a[1]) * 2.0 != sum
        || two_sum(sum, 1.0) != (a[0].max(a[1]) * 2.0, 0.0)
    {
        return None;
    }
    let series_first = z.abs().sqrt() < PFQ_ASYMPTOTIC_FIRST;
    let factor = |v| pfq_entire(&[], &[v], w, series_first).filter(|e| e.is_accurate());
    let (left, right) = (factor(u)?, factor(v)?);
    let y = left.y * right.y;
    let err = left.err + right.err + left.err * right.err + f64::EPSILON;
    let magnitude = y.abs();
    (err <= PFQ_TOL
        && y.is_normal()
        && magnitude.mul_add(err, magnitude).is_finite()
        && magnitude.mul_add(-err, magnitude) >= f64::MIN_POSITIVE)
        .then(|| Estimate::new(y, err))
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
/// assert_eq!(y, hyppfq([1.0, 2.0], [3.0], -0.5));
/// // empty lists are allowed
/// assert_eq!(hyp![; ; 0.0], Some(1.0));
/// ```
#[macro_export]
macro_rules! hyp {
    ($($a:expr),* $(,)? ; $($b:expr),* $(,)? ; $z:expr $(,)?) => {
        $crate::hyppfq([$($a),*], [$($b),*], $z)
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_absolute_at_tol() {
        // `err / y` rounds to the float above `PFQ_TOL`
        let y = 0.012_571_498_510_638_38;
        assert!(Estimate::absolute((y, PFQ_TOL * y)).is_accurate());
    }

    #[test]
    fn test_hyp4f3_racah() {
        // Noncanceling interior modes from exact rational sums.
        let a = [-18.0, 19.75, -23.0, 24.75];
        let b = [-24.0, 27.25, 1.25];
        let expected = -8.744_675_708_675_688e-6;
        let value = hyppfq(a, b, 1.0).unwrap();
        assert!((value / expected - 1.0).abs() < 1e-14);
        assert_eq!(hyppfq([a[2], a[0], a[3], a[1]], [b[2], b[0], b[1]], 1.0), Some(value));

        let value = hyppfq([-20.0, 21.75, -22.0, 23.75], b, 1.0).unwrap();
        assert!((value / 2.514_243_963_935_837_5e-6 - 1.0).abs() < 1e-14);

        // A true denominator pole occurs before termination.
        assert_eq!(hyppfq(a, [-17.0, b[1], b[2]], 1.0), Some(f64::INFINITY));
        // Private fallback decisions leave room for better generic answers.
        assert!(hyp4f3_racah(a, [-18.0, b[1], b[2]]).is_none());
        assert!(hyp4f3_racah(a, [b[0], b[1].next_up(), b[2]]).is_none());
        assert!(hyp4f3_racah([-1.0, 8.0, -2.0, 3.0], [-4.0, 12.0, 1.0]).is_none());
    }

    #[test]
    fn test_root_dd() {
        use crate::utils::EPSILON_SQ;
        // of a double-double with a low part
        let v = Dd::new(10.0) / 3.0;
        for k in 1..=5 {
            let root = root_dd(v, k);
            let pow = (0..k).fold(Dd::ONE, |p, _| p * root);
            assert!((pow - v).hi.abs() <= pow2(-100) * v.hi, "k = {k}");
        }
        // subnormal input for k = 2 stays full double-double precision
        let sub = 3.0 * f64::from_bits(1);
        let r_sub = root_dd(sub, 2);
        let expected_sub = Dd::new(0.75).sqrt().ldexp(-536);
        assert!(
            (r_sub - expected_sub).hi.abs() <= 4.0 * EPSILON_SQ * expected_sub.hi,
            "{r_sub:?} vs {expected_sub:?}"
        );
        // k = 4 near f64::MAX does not overflow `pow(4)` in Newton's step
        let r_max = root_dd(f64::MAX, 4);
        assert!(r_max.hi.is_finite());
        let (rm, re) = r_max.frexp();
        let pow4 = (0..4).fold(Dd::ONE, |p, _| p * rm);
        let (vm, ve) = Dd::new(f64::MAX).frexp();
        let diff = (pow4 - vm.ldexp(ve - 4 * re)).hi.abs();
        assert!(diff <= 20.0 * EPSILON_SQ * vm.hi, "{diff:e}");
        // huge k = 3 input stays within `pfq_exp_x_rel`
        let big = 3.0 * pow2(900);
        let r_big = root_dd(big, 3);
        let (bm, be) = r_big.frexp();
        let pow3 = bm * bm * bm;
        let (zm, ze) = Dd::new(big).frexp();
        let diff3 = (pow3 - zm.ldexp(ze - 3 * be)).hi.abs();
        assert!(diff3 <= 3.0 * asymptotic::pfq_exp_x_rel(3) * zm.hi, "{diff3:e}");
        // symmetric scaling keeps 513 <= k <= 1024 accurate on subnormals and near f64::MAX
        for (k, val) in [(513, f64::from_bits(1)), (1024, f64::from_bits(1)), (1024, f64::MAX)] {
            let root = root_dd(val, k);
            assert!(root.hi.is_finite());
            let (mut pm, mut pe) = (Dd::ONE, 0_i64);
            let (rm, re) = root.frexp();
            for _ in 0..k {
                let (m, e) = (pm * rm).frexp();
                pm = m;
                pe += i64::from(re) + i64::from(e);
            }
            let (vm, ve) = Dd::new(val).frexp();
            let rel = (pm - vm.ldexp((i64::from(ve) - pe) as i32)).hi.abs() / vm.hi;
            assert!(
                rel <= k as f64 * asymptotic::pfq_exp_x_rel(k),
                "k={k} val={val:e} rel={rel:e}"
            );
        }
        assert!(asymptotic::pfq_exp_x_rel(1025).is_infinite());
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
    fn test_add_exp_subnormal_with_zero_err() {
        use crate::xsf::cephes::dd_real::LOG2;
        // value 1 represented with subnormal mantissa 2^-1074 and exp(1074 * ln 2), plus 0 +- 0.5
        let sub = Scaled::new(f64::from_bits(1), 0, 0.0)
            .unwrap()
            .with_l(LOG2 * 1074.0);
        let zero_err = Scaled::new(0.5, 0, 0.0).unwrap().into_zero_err();
        let (y, err) = add_exp([sub, zero_err]).unwrap().to_f64();
        assert!((y - 1.0).abs() <= 1e-14 && (err - 0.5).abs() <= 1e-13, "{y} (err {err})");
        // tiny nonzero term whose value underflows to 0 relative to `lmax` (`2^-1075`) still
        // adds its absolute uncertainty `0.5 * 2^-1075 * 1.7e308 ~ 2.10e-16` (`~4.20e-16` rel)
        let half = Scaled::new(0.5, 0, 0.0).unwrap();
        let (_, base_err) = add_exp([half]).unwrap().to_f64();
        let tiny_unc = Scaled::new(0.5, -1075, 1.7e308).unwrap();
        let (y, err) = add_exp([half, tiny_unc]).unwrap().to_f64();
        assert!(
            (y - 0.5).abs() <= 1e-15 && err >= base_err + 4e-16,
            "{y} (err {err:e}, base {base_err:e})"
        );
        // single term with huge relative error preserves the full `Dd` exponent low word
        let huge_err = Scaled::new(0.5, 0, 1e308).unwrap();
        let (y, err) = add_exp([huge_err]).unwrap().to_f64();
        assert!((y - 0.5).abs() <= 1e-15 && err >= 1e308, "{y} (err {err:e})");
    }
}
