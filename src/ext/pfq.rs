//! Generalized hypergeometric functions pFq of real arguments.
//!
//! The internal functions take `a` and `b` directly (without a synthetic `1` for `n!`), and
//! supply the `(1)_n = n!` factor in the denominator where series terms or ratio bounds are
//! evaluated.

use core::cell::OnceCell;
use core::f64::consts::{PI, TAU};
use core::iter::{from_fn, once, successors, zip};
use core::ops::{Add, ControlFlow};

use super::{Scaled, add_exp, add_exp_in, fold_exp, ln_dd, ln_gamma, mul_gamma_ratio};
use crate::utils::{EPSILON_SQ, is_nonpos_int, neg_one_pow, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{
    DoubleDouble as Dd, FRAC_2_PI_DD, LN_SQRT_2PI_DD, frexp_f64, ldexp_f64, two_prod, two_sum,
};
use crate::xsf::cephes::hyp2f1::HYP2F1_LOST;
use crate::xsf::cephes::{gamma, gammasgn, lgam, rgamma};
use crate::xsf::{cospi, sinpi, xsf_hypu_digits};

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

/// Most effective factors in `a` or in `b` (including the 1 of `n!`) for which `pfq_products` is
/// exact
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

/// The factor of the measured rounding errors for more than `PFQ_MEASURED_PARAMS` parameters,
/// as each of them adds its own to a term
fn pfq_rounding_scale(a: &[f64], b: &[f64]) -> f64 {
    ((a.len() + b.len() + 1) as f64 / PFQ_MEASURED_PARAMS).max(1.0)
}

/// `v` in descending order, as `pfq_ratio_bound` is tightest if `a` and `b` are sorted alike
fn descending<P: AsMut<[f64]>>(mut v: P) -> P {
    v.as_mut().sort_unstable_by(|x, y| y.total_cmp(x));
    v
}

/// Smallest `f64` above a finite non-negative `x`, or `x` if non-finite
fn next_up_nonneg(x: f64) -> f64 {
    if x.is_finite() {
        f64::from_bits((x + 0.0).to_bits() + 1)
    } else {
        x
    }
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

/// `b` with the `1` of `n!` inserted in descending order (for a descending `b`)
fn with_one(b: &[f64]) -> impl Iterator<Item = f64> + Clone + '_ {
    let k = b.partition_point(|&v| v > 1.0);
    b[..k]
        .iter()
        .copied()
        .chain(once(1.0))
        .chain(b[k..].iter().copied())
}

/// Bound on `|t_{m+1} / t_m|` for all `m >= n`, if all `a_i + n` and `b_j + n` are positive,
/// except that for a longer `a` the unpaired `a_i + m` are taken at `m = n` (and they grow)
///
/// `b` is the complete denominator (including the `1` of `n!`), of length `b_len`. Each `a_i` is
/// paired with `b_i`: `(a_i + m) / (b_i + m)` decreases to 1, or stays below 1. The other
/// `1 / (b_j + m)` decrease. This is tightest if `a` and `b` are sorted alike.
///
/// The factors are multiplied as mantissas and binary exponents so intermediate products cannot
/// overflow to `inf` and divide the bound down to `0`, and any inexact product or quotient is
/// rounded up (including subnormal conversion) so the result never falls below the exact bound
/// for the given `a`, `b`, `z`, and `n`.
fn pfq_ratio_bound_iter(
    a: &[f64],
    b: impl Iterator<Item = f64> + Clone,
    b_len: usize,
    z: f64,
    n: f64,
) -> f64 {
    let az = z.abs();
    let scale = |(m, e, r): (f64, i64, usize), f: f64, divide: bool| {
        let ok = f.is_finite() && if divide { f > 0.0 } else { f >= 0.0 };
        ok.then(|| {
            let (fm, fe) = frexp_f64(f);
            let (m, me) = frexp_f64(if divide { m / fm } else { m * fm });
            let e = if divide {
                e - i64::from(fe)
            } else {
                e + i64::from(fe)
            };
            (m, e + i64::from(me), r + 2)
        })
    };
    let acc = (az > 0.0 && az.is_finite()).then(|| {
        let (zm, ze) = frexp_f64(az);
        let init = (zm, i64::from(ze), 0_usize);
        let acc = zip(a, b.clone()).try_fold(init, |acc, (&ai, bj)| {
            if ai > bj {
                scale(scale(acc, ai + n, false)?, bj + n, true)
            } else {
                (ai + n >= 0.0 && bj + n > 0.0 && (bj + n).is_finite()).then_some(acc)
            }
        })?;
        let acc = a
            .iter()
            .skip(b_len)
            .try_fold(acc, |acc, &ai| scale(acc, ai + n, false))?;
        b.skip(a.len())
            .try_fold(acc, |acc, bj| scale(acc, bj + n, true))
    });
    if az == 0.0 {
        0.0
    } else if let Some((m, e, roundings)) = acc.flatten() {
        if roundings == 0 {
            az
        } else if m == 0.0 {
            0.0
        } else {
            let slack = 1.0 + (roundings as f64 + 2.0) * f64::EPSILON;
            let bound = ldexp_i64(m * slack, e);
            if !bound.is_finite() {
                f64::INFINITY
            } else if bound < f64::MIN_POSITIVE {
                next_up_nonneg(bound)
            } else {
                bound
            }
        }
    } else {
        f64::INFINITY
    }
}

/// `pfq_ratio_bound` for a complete descending denominator `b` that already includes the `n + 1`
/// factorial factor (e.g. from `Residues` or already-shifted `pfq_alg_series_dd`)
fn pfq_ratio_bound_full(a: &[f64], b: &[f64], z: f64, n: f64) -> f64 {
    pfq_ratio_bound_iter(a, b.iter().copied(), b.len(), z, n)
}

/// Bound on `|t_{m+1} / t_m|` for all `m >= n` (if all `a_i + n` and `b_j + n` are positive,
/// except that for `a.len() > b.len() + 1` the unpaired `a_i + m` are taken at `m = n` and grow),
/// for descending `a` and `b` without the `1` of `n!`, which is inserted in descending order
fn pfq_ratio_bound(a: &[f64], b: &[f64], z: f64, n: f64) -> f64 {
    pfq_ratio_bound_iter(a, with_one(b), b.len() + 1, z, n)
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

/// Certifies that for `up_s.len() > low_s.len()`, term ratios stay at least 1 once they reach it
///
/// The ratios are at most 1 up to `m0 = n + delta`, where `base` bounds `|z|` times the paired
/// ratios and `delta = base^(-1/k) - x_s` for the largest unpaired `x_s = up_s[low_s.len()]`. From
/// `m0` on they grow if the derivative of their log is positive: the `k` unpaired factors add at
/// least `k / (x_s + m - n)`, and each pair with `u > l` subtracts
/// `(u - l) / ((u + m - n) (l + m - n))`. Times `x_s + m - n`, that is at most
/// `(u - l) / (l + delta)` times `max(1, (x_s + delta) / (u + delta))` (for descending `up_s`,
/// `(x_s + m - n) / (u + m - n)` is at most 1 and increases toward 1).
fn regrowth_certified(
    up_s: &[f64],
    low_s: &[f64],
    unpaired: &[f64],
    z: f64,
    bound: f64,
) -> Option<()> {
    bound.is_finite().then_some(())?;
    let (k, x_s) = (unpaired.len() as f64, unpaired[0]);
    // conservative upper bound on `base`, giving a conservative lower bound on `delta`
    let base = pfq_ratio_bound_full(&up_s[..low_s.len()], low_s, z, 0.0);
    base.is_normal().then_some(())?;
    // certify `base * root^k <= 1` via an upward-rounded scaled product
    let fits = |r: f64| {
        (r > 0.0 && r.is_finite()).then(|| {
            let ((bm, be), (rm, re)) = (frexp_f64(base), frexp_f64(r));
            let (m, e) = (0..unpaired.len()).fold((bm, i64::from(be)), |(m, e), _| {
                let (nm, ne) = frexp_f64(m * rm);
                (nm, e + i64::from(re) + i64::from(ne))
            });
            let slack = 1.0 + (2 * unpaired.len() + 2) as f64 * f64::EPSILON;
            ldexp_i64(m * slack, e) <= 1.0
        }) == Some(true)
    };
    let r0 = base.powf(-k.recip()) * (1.0 - 128.0 * f64::EPSILON);
    let root = [r0, r0 * (1.0 - 1024.0 * f64::EPSILON)]
        .into_iter()
        .find(|&r| fits(r))?;
    let delta = root - x_s;
    let delta = if delta > 0.0 {
        f64::from_bits(delta.to_bits() - 1)
    } else {
        delta
    };
    let pulls = zip(up_s, low_s).filter(|&(u, l)| u > l);
    let pull = pulls.map(|(&u, &l)| (u - l) / (l + delta) * ((x_s + delta) / (u + delta)).max(1.0));
    let pull_slack = 1.0 + (4.0 * low_s.len() as f64 + 4.0) * f64::EPSILON;
    (delta >= 0.0 && pull.sum::<f64>() * pull_slack <= k).then_some(())
}

/// Reusable scratch buffers for shifted `Dd` ratio/regrowth bounds
struct BoundScratch {
    up: Scratch<f64>,
    low: Scratch<f64>,
}

impl BoundScratch {
    fn new(up_len: usize, low_len: usize) -> Self {
        Self {
            up: Scratch::new(up_len, 0.0),
            low: Scratch::new(low_len, 0.0),
        }
    }

    fn shift<'a>(&'a mut self, up: &[Dd], low: &[Dd], n: f64) -> (&'a [f64], &'a [f64]) {
        let up_s = self.up.fill(up.len(), up.iter().map(|&x| dd_hi_up(x + n)));
        let low_s = self
            .low
            .fill(low.len(), low.iter().map(|&x| dd_hi_down(x + n)));
        (up_s, low_s)
    }

    /// `pfq_ratio_bound_full` on `up + n` (rounded up) and `low + n` (rounded down) in `Dd`
    fn ratio_bound(&mut self, up: &[Dd], low: &[Dd], z: f64, n: f64) -> f64 {
        let (up_s, low_s) = self.shift(up, low, n);
        pfq_ratio_bound_full(up_s, low_s, z, 0.0)
    }

    /// `ratio_bound`, or infinite unless ratios for `up.len() > low.len()` stay at least 1 once
    /// they reach it
    fn regrowth_bound(&mut self, up: &[Dd], low: &[Dd], z: f64, n: f64) -> f64 {
        let (up_s, low_s) = self.shift(up, low, n);
        let bound = pfq_ratio_bound_full(up_s, low_s, z, 0.0);
        match up_s.get(low_s.len()..).filter(|v| !v.is_empty()) {
            None => bound,
            Some(unpaired) => regrowth_certified(up_s, low_s, unpaired, z, bound)
                .map_or(f64::INFINITY, |()| bound),
        }
    }
}

/// The products of the `a_i + n` and of `(1 + n)` times the `b_j + n`, whose quotient is
/// `t_{n+1} / (z t_n)`
///
/// They and their quotient are normal for up to `PFQ_MAX_UNCHECKED` factors with sizes in
/// `PFQ_MODERATE`, as each nonzero `|a_i + n|` is then in `[2^-100, 2^101]`.
fn pfq_products<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> [T; 2] {
    let shifted = |p: T, &v: &f64| p * (T::from_f64(v) + n);
    let one = T::from_f64(1.0);
    [a.iter().fold(one, shifted), b.iter().fold(one + n, shifted)]
}

/// Whether `v` is normal in the arithmetic of `T`
fn is_normal<T: SeriesReal>(v: T) -> bool {
    (T::MIN_POSITIVE..=f64::MAX).contains(&v.to_f64().abs())
}

/// `t_{n+1} / (z t_n)` for any parameters and `n`, as a mantissa and an `i64` exponent, so that
/// no partial product over- or underflows
fn pfq_frexp_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> (T, i64) {
    let product = |params: &[f64], init: (T, i64)| {
        let times = |(m, e): (T, i64), &v: &f64| {
            let (vm, ve) = (T::from_f64(v) + n).frexp();
            let (m, me) = (m * vm).frexp();
            (m, e + i64::from(ve) + i64::from(me))
        };
        params.iter().fold(init, times)
    };
    let (fm, fe) = (T::from_f64(1.0) + n).frexp();
    let ((num, ne), (den, de)) =
        (product(a, (T::from_f64(1.0), 0)), product(b, (fm, i64::from(fe))));
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
    } else if ln_abs(a) - ln_abs(b) - (1.0 + n).ln() < PFQ_LN_UNDERFLOW {
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
    let many = a.len().max(b.len() + 1) > PFQ_MAX_UNCHECKED;
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

/// `t_{n+1} / t_n` as `Ok(ratio)` if the products, their quotient, and the product with `z` are
/// normal, else `Err(ratio)` with the already-computed `q * z`
///
/// With `checked = false`, the parameters must not be `pfq_extreme` so that partial products of
/// `pfq_products` stay normal.
#[inline]
fn pfq_normal_ratio<T: SeriesReal>(
    a: &[f64],
    b: &[f64],
    n: f64,
    z: f64,
    checked: bool,
) -> Result<T, T> {
    let (q, products_normal) = if checked {
        (pfq_checked_ratio(a, b, n), true)
    } else {
        let [num, den]: [T; 2] = pfq_products(a, b, n);
        (num / den, is_normal(num) && is_normal(den))
    };
    let ratio = q * T::from_f64(z);
    if products_normal && is_normal(q) && is_normal(ratio) {
        Ok(ratio)
    } else {
        Err(ratio)
    }
}

/// The term `t_{n+1}` after `term = t_n`, with a ratio that is `checked` if the parameters are
/// extreme
#[inline]
fn pfq_next<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, z: f64, term: T, checked: bool) -> T {
    match pfq_normal_ratio::<T>(a, b, n, z, checked) {
        Ok(ratio) => term * ratio,
        Err(ratio) => pfq_rescaled_next(a, b, n, z, term).unwrap_or_else(|| term * ratio),
    }
}

/// `pfq_next` without its fallback to the rounded ratio, which can have more rounding errors, so
/// `None` where `t_{n+1}` is not normal
fn pfq_normal_next(a: &[f64], b: &[f64], n: f64, z: f64, term: f64, checked: bool) -> Option<f64> {
    match pfq_normal_ratio::<f64>(a, b, n, z, checked) {
        Ok(ratio) => Some(term * ratio).filter(|next| next.is_normal()),
        Err(_) => pfq_rescaled_next(a, b, n, z, term),
    }
}

/// `pfq_next` from the mantissas and exponents of its factors, as an intermediate product or
/// quotient can over- or underflow where `t_{n+1}` is normal; `None` if `t_{n+1}` is not normal
/// either
#[cold]
#[inline(never)]
fn pfq_rescaled_next<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, z: f64, term: T) -> Option<T> {
    let (ratio, e) = pfq_frexp_ratio::<T>(a, b, n);
    let ((zm, ze), (tm, te)) = (frexp_f64(z), term.frexp());
    let rescaled = ldexp_i64(tm * (ratio * T::from_f64(zm)), e + i64::from(ze + te));
    is_normal(rescaled).then_some(rescaled)
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
fn pfq_series(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> (f64, f64) {
    let (y, err, _) = pfq_series_tails(a, b, z);
    (y, err)
}

/// The series in double-double with up to `PFQ_LONG_MAXITER` terms, if it is accurate, as near
/// `z = 1` it can need more than `pfq_series`
fn pfq_series_long(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> Option<(f64, f64)> {
    let (y, err, _) = pfq_sum::<Dd>(a, b, z, PFQ_LONG_MAXITER);
    let y = (f64::from(y), err);
    Estimate::absolute(y).is_accurate().then_some(y)
}

/// Largest partial sum, in units of `unit`, whose square is added to the sums of squares
const PFQ_SQUARE_MAX: f64 = 1e150;

/// The factor 2^-400 by which the units of those sums shrink once it is exceeded
const PFQ_UNIT: f64 = pow2(-400);

/// `pfq_series`, and the root of the sum of its squared tails `s - s_n`, which is at least the
/// largest of them
fn pfq_series_tails(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> (f64, f64, f64) {
    let (y, err, tails) = pfq_sum::<f64>(a, b, z, PFQ_MAXITER);
    // (not in double-double, whose `EPSILON` leaves room for many more parameters)
    let err = err * pfq_rounding_scale(a.as_ref(), b.as_ref());
    if err <= PFQ_TOL * y.abs() || err == f64::INFINITY {
        (y, err, tails)
    } else {
        let (y, err, tails) = pfq_sum::<Dd>(a, b, z, PFQ_DD_MAXITER);
        (f64::from(y), err, tails)
    }
}

/// Maximum number of roundings of a ratio `t_{k+1} / t_k`: each `a_i + k`, `b_j + k`, and
/// `1 + k`, the product, the quotient, and the factors `z` and `t_k`
fn pfq_ratio_roundings(a: &[f64], b: &[f64]) -> f64 {
    2.0 * (a.len() + b.len() + 1) as f64 + 3.0
}

/// `pfq_series_tails` in the arithmetic of `T`, with at most `maxiter` terms, and with the sum
/// in that arithmetic
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
) -> (T, f64, f64) {
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
) -> (T, f64, f64) {
    // generic over the storage of the parameters, as the sum is faster if it is compiled for
    // the numbers of parameters of its caller
    let (a, b) = (a.as_ref(), b.as_ref());
    let settled = pfq_settled(a, b);
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let (mut sum, mut comp) = (T::from_f64(1.0), T::from_f64(0.0));
    // the compensated partial sum `s_n` and its offset `d_n = s_n - 1` as `f64`
    let (mut s, mut d) = T::sum_and_delta(sum, comp);
    // the sums of `d_n`, `d_n^2`, and `s_n^2`, in units of `unit` so that squares don't overflow;
    // shifting by `t_0 = 1` keeps `s - s_n = d - d_n` accurate near a flat sum
    let (mut dsum1, mut dsum2, mut sum2, mut unit) = (0.0, 0.0_f64, 0.0_f64, 1.0);
    // the sum of |t_k|, in units of `unit`
    let mut abs1 = 1.0;
    let ratio_roundings = pfq_ratio_roundings(a, b);
    // the rounding error, and the root of the sum of (s - s_n)^2 = (d - d_n)^2
    let estimates = |n: f64, d: f64, dsum1: f64, dsum2: f64, sum2: f64, unit: f64, abs1: f64| {
        // scaled by max(|d * unit|, sqrt(dsum2)) so that neither tiny `d` nor huge sums overflow
        let root2 = sum2.sqrt();
        let (du, droot2) = (d * unit, dsum2.sqrt());
        let scale = du.abs().max(droot2);
        let tails = if scale == 0.0 {
            0.0
        } else {
            let (r, r2) = (du / scale, droot2 / scale);
            ((n + 1.0) * (r * r) - 2.0 * r * (dsum1 / scale) + r2 * r2)
                .max(0.0)
                .sqrt()
                * scale
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
    let fail = ControlFlow::Break(None);
    let res = (0..maxiter).map(f64::from).try_for_each(|n| {
        while s.abs().max(d.abs()) * unit > PFQ_SQUARE_MAX {
            // the squares would overflow, so the sums are rescaled (exactly, by a power of 2)
            let u2 = PFQ_UNIT * PFQ_UNIT;
            (dsum1, dsum2, sum2, unit) = (dsum1 * PFQ_UNIT, dsum2 * u2, sum2 * u2, unit * PFQ_UNIT);
            abs1 *= PFQ_UNIT;
        }
        let (su, du) = (s * unit, d * unit);
        dsum1 += du;
        dsum2 += du * du;
        sum2 += su * su;
        if n == degree {
            // the series terminated
            ControlFlow::Break(Some(estimates(n, d, dsum1, dsum2, sum2, unit, abs1)))
        } else if let Some(term) = terms.next() {
            let t = term.to_f64();
            abs1 += t.abs() * unit;
            (sum, comp) = T::add_compensated(sum, comp, term);
            (s, d) = T::sum_and_delta(sum, comp);
            if !s.is_finite() {
                fail
            } else if n >= settled && t.abs() <= tol * s.abs() {
                let rho = pfq_ratio_bound(a, b, z, n + 1.0);
                // geometric bound on the remaining terms, if the series converges (p <= q + 1)
                let tail = t.abs() * rho / (1.0 - rho);
                if a.len() <= b.len() + 1 && rho < 1.0 && tail <= tol * s.abs() {
                    let (rounding, tails) = estimates(n, d, dsum1, dsum2, sum2, unit, abs1);
                    ControlFlow::Break(Some((rounding + tail, tails)))
                } else {
                    ControlFlow::Continue(())
                }
            } else {
                ControlFlow::Continue(())
            }
        } else {
            fail
        }
    });
    let (err, tails) = res
        .break_value()
        .flatten()
        .unwrap_or((f64::INFINITY, f64::INFINITY));
    (sum + comp, err, tails)
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
                && a.len() <= b.len() + 1
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
fn gamma_ratio(u: (f64, f64), w: (f64, f64)) -> Option<(f64, f64)> {
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
trait Params: Clone + AsRef<[f64]> + AsMut<[f64]> {}

impl<P: Clone + AsRef<[f64]> + AsMut<[f64]>> Params for P {}

/// `v` without `v[i]`
fn skip(v: &impl AsRef<[f64]>, i: usize) -> impl Iterator<Item = f64> + Clone {
    let v = v.as_ref();
    (0..v.len()).filter(move |&k| k != i).map(move |k| v[k])
}

/// The analytic continuation to `z < -1` (DLMF 16.8.8) and an estimate of its absolute error, or
/// `None` if a factor of a coefficient overflows or underflows
///
/// This is for `p = q + 1`. It is a sum over `j` of series in `1 / z` with the parameters
/// `a_j`, `a_j + 1 - b`, and `a_j + 1 - a` (`k != j`), or its limit where `a_i` differ by
/// integers, which is `None` if it underflows or overflows.
fn pfq_continuation<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<(f64, f64)> {
    let (a_ref, b_ref) = (a.as_ref(), b.as_ref());
    let mut diffs = zip(a_ref, 0..).flat_map(|(aj, j)| skip(&a_ref, j).map(move |ak| ak - aj));
    if diffs.any(|d| d.fract() == 0.0) {
        // the terms of DLMF 16.8.8 are infinite: their limit is the sum of the terms of `pfq_alg`
        let (y, err) = add_exp(pfq_alg(a_ref, b_ref, z, Alg::Value)?)?.to_f64();
        return y.is_normal().then_some((y, err * y.abs()));
    }
    let (mut sum, mut err) = (0.0, 0.0);
    for (j, &aj) in a_ref.iter().enumerate() {
        // Gamma(a_k - a_j) / Gamma(a_k) and Gamma(b_k) / Gamma(b_k - a_j) for each k != j
        let ratios = skip(&a_ref, j).map(|ak| gamma_ratio(two_sum(ak, -aj), (ak, 0.0)));
        let ratios = ratios.chain(
            b_ref
                .iter()
                .map(|&bk| gamma_ratio((bk, 0.0), two_sum(bk, -aj))),
        );
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
        // a_j, a_j + 1 - b, and a_j + 1 - a (k != j), and the relative effect of their rounding
        // errors on the terms of their series, which is large near a non-positive integer -m:
        // there it is on the terms from n = m on, and not on those before them
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
        let mut upper = a.clone();
        upper.as_mut()[0] = aj;
        for (dst, &bk) in zip(&mut upper.as_mut()[1..], b_ref) {
            *dst = shift(bk);
        }
        let mut lower = b.clone();
        for (dst, ak) in zip(lower.as_mut(), skip(&a_ref, j)) {
            *dst = shift(ak);
        }
        let (upper, lower) = (descending(upper), descending(lower));
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
type Method<A, B> = fn(&A, &B, f64) -> Option<(f64, f64)>;

/// `_{q+1}F_q` for `z < 1`, or for a polynomial, by the first method for `z` that is accurate,
/// else by the most accurate one
///
/// Unless it is a polynomial, the last methods are `pfq_conformal`, which is the slowest, and for
/// `|z| < 1` `pfq_recurrence`, `pfq_series_long`, `pfq_recurrence_lowered`,
/// `pfq_recurrence_balanced` and `pfq_recurrence_nested`, of which all but `pfq_series_long` only
/// apply to 3F2, and then `pfq_overflow`.
fn pfq_accurate<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<Estimate> {
    let series: Method<A, B> = |a, b, z| Some(pfq_series(a, b, z));
    let polynomial = first_zero(a.as_ref()).is_finite();
    let methods: &[Method<A, B>] = if polynomial || z > PFQ_WENIGER_MAX_Z {
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
    let last: Option<Method<A, B>> = (!polynomial).then_some(conformal::pfq_conformal);
    let inside: &[Method<A, B>] = if !polynomial && z.abs() < 1.0 {
        &[
            recurrence::pfq_recurrence,
            pfq_series_long,
            recurrence::pfq_recurrence_lowered,
            recurrence::pfq_recurrence_balanced,
            recurrence::pfq_recurrence_nested,
        ]
    } else {
        &[]
    };
    let methods = methods.iter().chain(&last).chain(inside);
    let y = methods.fold(None, |y, method| y.or_better(|| method(a, b, z).map(Estimate::absolute)));
    if polynomial {
        y
    } else {
        y.or_better(|| pfq_overflow(a.as_ref(), b.as_ref(), z))
    }
}

/// Most terms that `pfq_overflow` follows
///
/// Of its 57 results without it on 1657 points where the other methods are not accurate, 55 are
/// within it, and 54 within half of it (measured).
const PFQ_OVERFLOW_MAXITER: u32 = 4 * PFQ_MAXITER;

/// `+-inf` for `0 < z < 1` where the series shows that `_{q+1}F_q` overflows, else `None`
///
/// From the first `n >= pfq_settled` on, the terms have the sign of `t_n`. Once the sum of their
/// sizes is at least `2^1025`, and that of the terms before at most `2^1023`, the result has that
/// sign and a size of at least `2^1024`. Each term has at most `pfq_ratio_roundings` roundings per
/// step before it, and each sum one per term.
fn pfq_overflow(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    if z <= 0.0 {
        return None;
    }
    let (settled, checked) = (pfq_settled(a, b), pfq_extreme(a, b));
    let roundings = pfq_ratio_roundings(a, b) + 1.0;
    let next = |n, term| pfq_normal_next(a, b, n, z, term, checked);
    let indices = (0..PFQ_OVERFLOW_MAXITER).map(f64::from);

    // the sum of the sizes of the terms before `settled`, and the first term after them
    let (before, mut t) = indices
        .clone()
        .take_while(|&n| n < settled)
        .try_fold((0.0, 1.0), |(sum, t): (f64, f64), n| Some((sum + t.abs(), next(n, t)?)))?;

    // the sum of the sizes of the terms from `settled` on, and `t`, both in units of `2^e`
    let mut after = 0.0;
    let mut e = 0;
    for n in indices.skip_while(|&n| n < settled) {
        // the sum in [1/2, 1), exactly, and t with it, which is exact unless it underflows
        let (m, k) = frexp_f64(after + t.abs());
        after = m;
        e += k;
        t = ldexp_f64(t, -k);
        if !t.is_normal() {
            return None;
        }
        // (a bound of the relative rounding errors, with a factor 2 to spare)
        let err = (n + 1.0) * roundings * f64::EPSILON;
        // (the lower bound `>= 2^1025` by its exponent, as it can be beyond `f64`)
        if before * (1.0 + err) <= pow2(1023) && frexp_f64(after * (1.0 - err)).1 + e > 1025 {
            return Some(Estimate::exact(f64::INFINITY.copysign(t)));
        }
        let next = next(n, t)?;
        if next.abs() < t.abs() {
            // the rest is at most a geometric series
            let rho = pfq_ratio_bound(a, b, z, n + 1.0);
            if rho < 1.0 && before + ldexp_f64(after + next.abs() / (1.0 - rho), e) < pow2(1023) {
                return None;
            }
        }
        t = next;
    }
    None
}

/// `pfq_series` after sorting the parameters, which tightens `pfq_ratio_bound`
fn pfq_series_sorted(a: impl Params, b: impl Params, z: f64) -> Estimate {
    Estimate::absolute(pfq_series(&descending(a), &descending(b), z))
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

/// `_{q+1}F_q` for finite arguments
fn pfq_real<A: Params, B: Params>(a: A, b: B, z: f64) -> Option<Estimate> {
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
        let diverges = diverges_at_one(a, b)?;
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
/// - $\pm\infty$ where it detects that it overflows, for $0 < z < 1$;
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
    let (a, b) = ([a1, a2, a3], [b1, b2]);
    // a common parameter cancels: 2F1, or 1F0 if that has one as well
    let mut pairs = (0..3).flat_map(|i| (0..2).map(move |j| (i, j)));
    let Some((i, j)) = pairs.find(|&(i, j)| a[i] == b[j]) else {
        return pfq_real(a, b, z);
    };
    let (a, b) = ([a[(i + 1) % 3], a[(i + 2) % 3]], [b[1 - j]]);
    match a.iter().position(|&v| v == b[0]) {
        // the closed form of `pfq_estimate` is also real for z > 1 if `a` is an integer
        Some(i) => pfq_estimate(&[a[1 - i]], &[], z),
        None if z == 1.0 => gauss_sum(a[0], a[1], b[0]).or_better(|| pfq_real(a, b, z)),
        // `pfq_real` is NaN there, as complex, but this 2F1 can be real
        None if z > 1.0 && first_zero(&a).is_infinite() => pfq_estimate(&a, &b, z),
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

/// Incremental state for `asymptotic_sum` and each residue class of `pfq_exp_cos`
#[derive(Clone, Copy)]
struct AsymptoticSum<V> {
    settled: f64,
    z_abs: f64,
    remainder: f64,
    sum: V,
    abs_sum: f64,
    prev: f64,
    min: f64,
    best: (V, f64),
    done: bool,
}

impl<V: SeriesReal> AsymptoticSum<V> {
    fn new(settled: f64, z: f64) -> Self {
        let zero = V::from_f64(0.0);
        let z_abs = z.abs();
        Self {
            settled,
            z_abs,
            remainder: (TAU * z_abs).sqrt().max(1.0),
            sum: zero,
            abs_sum: 0.0,
            prev: f64::INFINITY,
            min: f64::INFINITY,
            best: (zero, f64::INFINITY),
            done: false,
        }
    }

    /// Advances by term `n` (with `next_bound` from the lookahead term); breaks when done
    fn step(
        &mut self,
        n: f64,
        t: V,
        bound: f64,
        next_bound: f64,
        mut ratio_bound: impl FnMut(f64) -> f64,
    ) -> ControlFlow<()> {
        let size = bound.max(next_bound);
        let stop = !t.to_f64().is_finite() || n >= self.settled + self.z_abs && bound >= self.prev;
        let converged = !stop && n >= self.settled && {
            // a later sum can disqualify the current best even when its term is larger
            let err = self.remainder * size + n * V::EPSILON * self.abs_sum;
            self.best.1 = self
                .best
                .1
                .max((self.sum - self.best.0).to_f64().abs() - err);
            if err < self.best.1 {
                self.best = (self.sum, err);
            }
            if size < self.min {
                self.min = size;
                self.remainder * self.min <= V::EPSILON * self.sum.to_f64().abs()
                    && ratio_bound(n) < 1.0
            } else {
                false
            }
        };
        if stop || converged {
            self.done = true;
            ControlFlow::Break(())
        } else {
            self.sum = self.sum + t;
            self.abs_sum += bound;
            self.prev = bound;
            // (at the cap, terms that grew after the smallest may still shrink to smaller ones)
            if n + 1.0 == f64::from(PFQ_MAXITER) && size > self.min {
                self.best.1 = f64::INFINITY;
            }
            ControlFlow::Continue(())
        }
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
    mut ratio_bound: impl FnMut(f64) -> f64,
) -> (T::Value, f64) {
    let mut state = AsymptoticSum::<T::Value>::new(settled, z);
    let mut terms = terms.peekable();
    let _ = (0..PFQ_MAXITER)
        .map(f64::from)
        .try_for_each(|n| match terms.next() {
            Some(term) => {
                let next = terms.peek().map_or(0.0, |next| next.bound());
                state.step(n, term.value(), term.bound(), next, &mut ratio_bound)
            },
            None => ControlFlow::Break(()),
        });
    state.best
}

/// The sum of a convergent series, in the arithmetic of its terms, and an estimate of its
/// absolute error, which is infinite if it did not converge within `PFQ_MAXITER` terms
///
/// From term `settled` on, the bound on each term is at most `ratio_bound(n)` times the one
/// before it, and `ratio_bound` does not increase.
fn convergent_sum<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    mut ratio_bound: impl FnMut(f64) -> f64,
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

/// Sum of a terminating series of known `degree`, with its rounding error estimate
///
/// `step_roundings` bounds the rounding errors accumulated per recurrence step when generating
/// each term (beyond the summation rounding itself).
fn polynomial_sum<T: SeriesReal>(
    mut terms: impl Iterator<Item = T>,
    degree: f64,
    step_roundings: f64,
) -> (T, f64) {
    let zero = T::from_f64(0.0);
    let step = |(sum, err_sum): (T, f64), k: u32| match terms.next() {
        Some(t) if t.to_f64().is_finite() => {
            let next = sum + t;
            if next.to_f64().is_finite() {
                let weight = ((degree + 1.0) + f64::from(k) * step_roundings) * T::EPSILON;
                ControlFlow::Continue((next, err_sum + weight * t.to_f64().abs()))
            } else {
                ControlFlow::Break(next)
            }
        },
        _ => ControlFlow::Break(sum),
    };
    if degree >= f64::from(PFQ_MAXITER) {
        (zero, f64::INFINITY)
    } else {
        match (0..=degree as u32).try_fold((zero, 0.0), step) {
            ControlFlow::Continue((sum, err_sum)) => (sum, err_sum),
            ControlFlow::Break(sum) => (sum, f64::INFINITY),
        }
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
    ratio_bound: impl FnMut(f64) -> f64,
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

/// Whether `up` and `low` (complete `Dd` denominator) need `frexp` rescaling in `pfq_next_dd`
fn pfq_extreme_dd(up: &[Dd], low: &[Dd]) -> bool {
    let many = up.len().max(low.len()) > PFQ_MAX_UNCHECKED;
    let extreme = |x: &Dd| {
        !PFQ_MODERATE.contains(&x.hi.abs()) || (x.lo != 0.0 && !PFQ_MODERATE.contains(&x.lo.abs()))
    };
    many || up.iter().chain(low).any(extreme)
}

/// `t_{n+1}` after `term = t_n` for complete `Dd` parameters `up` and `low`, preserving `Dd` low
/// parts through `(x + n).frexp()` when intermediate products are extreme or non-normal
#[cold]
#[inline(never)]
fn pfq_rescaled_next_dd(up: &[Dd], low: &[Dd], n: f64, w: f64, term: f64) -> f64 {
    let product = |params: &[Dd]| {
        params.iter().fold((1.0_f64, 0_i64), |(m, e), &x| {
            let (vm, ve) = (x + n).frexp();
            let (m, me) = frexp_f64(m * vm.hi);
            (m, e + i64::from(ve) + i64::from(me))
        })
    };
    let ((num, ne), (den, de)) = (product(up), product(low));
    let ((wm, we), (tm, te)) = (frexp_f64(w), frexp_f64(term));
    ldexp_i64(tm * ((num / den) * wm), (ne - de) + i64::from(we + te))
}

/// `t_{n+1}` after `term = t_n` for complete `Dd` parameters `up` and `low`
#[inline]
fn pfq_next_dd(up: &[Dd], low: &[Dd], n: f64, w: f64, term: f64, checked: bool) -> f64 {
    let normal = |v: f64| is_normal(v).then_some(v);
    let fast = (!checked)
        .then(|| {
            let num = normal(up.iter().fold(1.0_f64, |p, &x| p * (x + n).hi))?;
            let den = normal(low.iter().fold(1.0_f64, |p, &x| p * (x + n).hi))?;
            let ratio = normal(normal(num / den)? * w)?;
            Some(term * ratio)
        })
        .flatten();
    fast.unwrap_or_else(|| pfq_rescaled_next_dd(up, low, n, w, term))
}

/// The series of the term of `a_m` of the algebraic part (DLMF 16.11.2), as its sum and
/// absolute error
///
/// The shifted parameters `a_m`, `1 + a_m - b_j`, and `1 + a_m - a_i` (including `1` for `i = m`)
/// are kept in double-double so cancellations in `x + n`, pole checks, settlement, and bounds do
/// not lose their low parts. With `sizes`, it is the sum of the sizes of its terms.
fn pfq_alg_series(a: &[f64], b: &[f64], m: usize, z: f64, sizes: bool) -> (f64, f64) {
    let am = Dd::new(a[m]);
    let (mut up_s, mut low_s) =
        (Scratch::new(b.len() + 1, Dd::ZERO), Scratch::new(a.len(), Dd::ZERO));
    let up = up_s.fill(b.len() + 1, once(am).chain(b.iter().map(|&bj| (am - bj) + 1.0)));
    let low = low_s.fill(a.len(), a.iter().map(|&ai| (am - ai) + 1.0));
    sort_descending_dd(up);
    sort_descending_dd(low);
    let degree = first_zero_dd(up);
    let (p, q) = (a.len(), b.len());
    let w = neg_one_pow(p + q + 1) / z;
    let checked = pfq_extreme_dd(up, low);
    let series = once(1.0)
        .chain(pfq_terms_of(degree, w, |n, term| pfq_next_dd(up, low, n, w, term, checked)))
        .map(|t| if sizes { t.abs() } else { t });
    if degree.is_finite() {
        // each step rounds `(x + n).hi` and multiplies each factor, plus division and scaling
        polynomial_sum(series, degree, 2.0 * (up.len() + low.len()) as f64 + 2.0)
    } else {
        let settled = up
            .iter()
            .chain(&*low)
            .map(positive_from)
            .fold(0.0, f64::max);
        let kappa = (q + 1).saturating_sub(p);
        let mut scratch = BoundScratch::new(up.len(), low.len());
        let bound = |n| scratch.regrowth_bound(up, low, recip_up(z), n);
        pfq_alg_sum(series, settled, kappa, z, bound)
    }
}

/// `pfq_alg_series` for `p = q` or `p > q + 1` in double-double, with the exact parameters, and
/// an infinite error if no term bounds the rest
///
/// A product of the factors of a term ratio that isn't normal is less exact, and ends the terms.
fn pfq_alg_series_dd(a: &[f64], b: &[f64], m: usize, z: f64) -> (Dd, f64) {
    let [mut up, mut low] = residue::series_params(a, b, m);
    // (so that the factors of a term ratio pair by size, with the smallest of `up` unpaired)
    for v in [&mut up, &mut low] {
        sort_descending_dd(v);
    }
    let degree = first_zero_dd(&up);
    let ratio = |n: f64| {
        let normal = |p: Dd| if is_normal(p) { p } else { Dd::NAN };
        let prod = |v: &[Dd]| v.iter().fold(Dd::ONE, |p, &x| normal(p * (x + n)));
        // in (-1)^(p+q+1) / z
        if n == degree {
            Dd::ZERO
        } else {
            prod(&up) / normal(prod(&low) * (neg_one_pow(a.len() + b.len() + 1) * z))
        }
    };
    let terms = once(Dd::ONE).chain(pfq_terms_of(degree, 1.0, |n, t| t * ratio(n)));
    if degree.is_finite() {
        // (all of its terms, with per-step factor-product roundings in double-double)
        return polynomial_sum(terms, degree, (up.len() + low.len() + 2) as f64);
    } else if up.len() < low.len() {
        // (it converges for p > q + 1)
        let mut scratch = BoundScratch::new(up.len(), low.len());
        let bound = |n: f64| scratch.ratio_bound(&up, &low, recip_up(z), n);
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
    terms: &mut Vec<Scaled>,
) -> Option<()> {
    let am = a[m];
    let others = skip(&a, m).map(Dd::new);
    let diffs = others.clone().map(|ai| ai - am);
    let before = |(d, i): (Dd, usize)| d.is_integer() && (d.hi < 0.0 || d.hi == 0.0 && i < m);
    if diffs.clone().zip(0..).any(before) {
        return Some(());
    }
    if first_zero(a).is_finite() {
        return (!is_nonpos_int(am)).then_some(());
    }
    let residues = Residues::of_class(a, b, m, z);
    let [up, low] = residues.params();
    // (their order is final after the last zero factor, and their bounds are valid once all
    // parameters plus n exceed 1)
    let settled = up
        .iter()
        .chain(&low)
        .map(|&x| positive_from(&(x - 1.0)))
        .fold(0.0, f64::max);
    let mut scratch = BoundScratch::new(up.len(), low.len());
    let bound =
        |n| scratch.ratio_bound(&up, &low, recip_up(z), n) * Residues::growth_bound(&up, &low, n);
    let sizes = sizes || z > 0.0;
    let series = residues.map(|(t, size)| (if sizes { size } else { t }, size));
    let kappa = (b.len() + 1).saturating_sub(a.len());
    let (sum, err) = pfq_alg_sum(series, settled, kappa, z, bound);
    if (sum, err) != (0.0, 0.0) {
        let bs = b.iter().copied().map(Dd::new);
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
        terms.push(y);
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

/// The non-vanishing terms of the algebraic part of the asymptotic expansion of pFq for `p <= q`
/// (DLMF 16.11.2), with `|z|` instead of `-z`, or `None` if a factor of a term overflows or
/// underflows, or if a non-positive integer `a_i` differs from another `a_i` by an integer
///
/// For `p > q + 1` and `z < 0`, their series converge, and they sum to `pfq_borel`.
///
/// The term of `a_m` is `|z|^-a_m` times a `_{q+1}F_{p-1}` series in `(-1)^(p+q+1) / z`, with the
/// parameters `a_m` and `1 + a_m - b_j`, and `1 + a_m - a_i` (`i != m`). It is omitted if it
/// vanishes, where `1 / Gamma(b_j - a_m)` does. The `a_i` that differ by integers have one term,
/// by `pfq_class_term`, which for `z > 0` is a bound on its size.
///
/// With `Alg::Sizes`, each term is the sum of the sizes of the terms of its series, times the
/// size of its factor. Their sum bounds the algebraic part at any phase of `z` and of its terms.
///
/// With `Alg::Dd`, the series and the factor of each term are in double-double, except for the
/// `a_i` that differ by integers.
fn pfq_alg(a: &[f64], b: &[f64], z: f64, mode: Alg) -> Option<Vec<Scaled>> {
    let sizes = mode == Alg::Sizes;
    let ln_z = ln_dd(z.abs());
    let bs = b.iter().copied().map(Dd::new);
    let mut terms = Vec::with_capacity(a.len());
    for (m, &am) in a.iter().enumerate() {
        let others = skip(&a, m).map(Dd::new);
        let diffs = others.clone().map(|ai| ai - am);
        if diffs.clone().any(Dd::is_integer) {
            pfq_class_term(a, b, m, z, sizes, &mut terms)?;
        } else {
            let num = bs.clone().chain(diffs);
            let den = others.chain(bs.clone().map(|bj| bj - am));
            if !den.clone().any(|g| is_pole(&g)) {
                let l = -ln_z * am;
                let mut y = if mode == Alg::Dd {
                    let (sum, err) = pfq_alg_series_dd(a, b, m, z);
                    let y = Scaled::new(1.0, 0, err / sum.hi.abs())?.with_l(l);
                    mul_gammas_in(num, den, mul_dd(y, sum), true)?
                } else {
                    gamma_series(num, den, pfq_alg_series(a, b, m, z, sizes), l)?
                };
                if sizes {
                    y.m = y.m.abs();
                }
                terms.push(y);
            }
        }
    }
    Some(terms)
}

/// `pFq` for `p <= q` by series and asymptotic expansion, trying `series` first when `series_first`
fn pfq_entire<A: Params, B: Params>(a: &A, b: &B, z: f64, series_first: bool) -> Option<Estimate> {
    let Some(degree) = pfq_degree(a.as_ref(), b.as_ref()) else {
        return Some(Estimate::exact(f64::INFINITY));
    };
    let series = || Some(pfq_series_sorted(a.clone(), b.clone(), z));
    let asymptotic = || pfq_asymptotic(a.as_ref(), b.as_ref(), z);
    // last, as it is the slowest (and applies to `p = q`, `z <= 0` only)
    let asymptotic_dd = || pfq_asymptotic_dd(a.as_ref(), b.as_ref(), z);
    if degree.is_finite() {
        series()
    } else if series_first {
        series().or_better(asymptotic).or_better(asymptotic_dd)
    } else {
        asymptotic().or_better(series).or_better(asymptotic_dd)
    }
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
    let s = v.hi.powf(1.0 / k as f64);
    let pow = |n: usize| (0..n).fold(Dd::ONE, |p, _| p * s);
    match k {
        1 => v,
        2 => v.sqrt(),
        // (one Newton step)
        _ => Dd::new(s) + (v - pow(k)) / (k as f64 * pow(k - 1).hi),
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
        pfq_estimate(&[a], &[b1, b2], z)
    } else if z == 0.0 {
        Some(Estimate::exact(1.0))
    } else {
        let b = [b1, b2];
        let series_first = z > -HYP1F2_ASYMPTOTIC_FIRST;
        // a common parameter cancels, also at a pole; keep on `pfq_series`/`pfq_asymptotic`
        // rather than `hyp0f1` so large negative `z` stays accurate
        match b.iter().position(|&bj| bj == a) {
            Some(j) => pfq_entire(&[], &[b[1 - j]], z, series_first),
            None => pfq_entire(&[a], &b, z, series_first),
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

/// The terms `t_k = c_k x^-k` of the exponential part of the asymptotic expansion of
/// pFq for `p <= q` (DLMF 16.11.1)
///
/// With `D = x d/dx` and `kappa = q - p + 1`, the differential equation (DLMF 16.8.3) is
/// `(D prod_j (D + kappa (b_j - 1)) - x^kappa prod_i (D + kappa a_i)) w = 0`, and each factor
/// `D + g` maps `e^x x^s` to `e^x (x^(s+1) + (s + g) x^s)`. So `c_0 = 1` and
/// `kappa k c_k = sum_{m=1}^q d_{m+1}(nu - k + m) c_{k-m}`, where `d_j(s)` is the
/// coefficient of `e^x x^(s+q+1-j)` in the image of `e^x x^s`.
fn pfq_exp_terms(a: &[f64], b: &[f64], nu: f64, x: f64) -> impl Iterator<Item = f64> {
    let (q, kappa) = (b.len(), (b.len() + 1 - a.len()) as f64);
    let left: Vec<f64> = b
        .iter()
        .map(|bj| kappa * (bj - 1.0))
        .chain(once(0.0))
        .collect();
    let right: Vec<f64> = a.iter().map(|ai| kappa * ai).collect();
    // the coefficients of e^x x^(s+n-j) in prod_(i<n) (D + g_i) e^x x^s, for j <= q + 1
    let fill_coefs = |g: &[f64], s: f64, c: &mut [f64]| {
        c[0] = 1.0;
        c[1..].fill(0.0);
        for (i, &gi) in g.iter().enumerate() {
            for j in (1..=i + 1).rev() {
                c[j] += c[j - 1] * (s + (i + 1 - j) as f64 + gi);
            }
        }
    };
    // t_k, and the sums so far of t_(k+1), ..., t_(k+q)
    let (mut k, mut t, mut sums) = (0.0, 1.0, vec![0.0; q]);
    let (mut l, mut r) = (vec![0.0; q + 2], vec![0.0; q + 2]);
    from_fn(move || {
        let term = t;
        fill_coefs(&left, nu - k, &mut l);
        fill_coefs(&right, nu - k, &mut r);
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
    if k == 2 {
        (Dd::ZERO, Dd::ONE)
    } else {
        let phi = 2.0 / (FRAC_2_PI_DD * k as f64);
        // n, and the terms n - 1 of their Taylor series, of the degrees 2 n - 2 and 2 n - 1, which
        // are below 2^-106 from n = 21 on
        let terms = successors(Some((1.0, Dd::ONE, phi)), |&(n, c, _)| {
            let next = -(c * phi * phi) / ((2.0 * n - 1.0) * (2.0 * n));
            Some((n + 1.0, next, next * phi / (2.0 * n + 1.0)))
        });
        terms
            .take(20)
            .fold((Dd::ZERO, Dd::ZERO), |(cos, sin), (_, c, s)| (cos + c, sin + s))
    }
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
    terms: impl Iterator<Item = f64>,
    settled: f64,
    kappa: usize,
    root: Dd,
    nu: Dd,
) -> ((f64, f64), Dd) {
    let mut terms = terms.fuse();
    let kf = kappa as f64;
    let x = kf * root.hi;
    let (cos, sin) = cos_sin_pi_div(kappa);
    // theta / pi in double-double, reduced to [-1, 1], and a bound of the error of theta
    let turns = root * kf * sin * (FRAC_2_PI_DD * 0.5) + nu / kf;
    let t = turns - 2.0 * (0.5 * turns.hi).round();
    let dtheta = 4.0 * EPSILON_SQ * (x + PI * nu.hi.abs());
    // single pass over `terms`: each residue class `k = r mod kappa` tracks its own
    // `AsymptoticSum` state, with a `kappa`-term lookahead window for `t_{k + kappa}`
    let mut classes = vec![AsymptoticSum::<f64>::new(settled / kf, x / kf); kappa];
    let mut win: Vec<Option<f64>> = (0..kappa).map(|_| terms.next()).collect();
    let mut active = kappa;
    let _ = (0..PFQ_MAXITER).try_for_each(|j| {
        let (nf, sign) = (f64::from(j), if j % 2 == 0 { 1.0 } else { -1.0 });
        for (cls, slot) in zip(&mut classes, &mut win) {
            let raw = core::mem::replace(slot, terms.next());
            if !cls.done {
                let next = slot.map_or(0.0, f64::abs);
                match raw {
                    None => {
                        cls.done = true;
                        active -= 1;
                    },
                    Some(raw)
                        if cls
                            .step(nf, sign * raw, raw.abs(), next, |_| 0.0)
                            .is_break() =>
                    {
                        active -= 1;
                    },
                    Some(_) => {},
                }
            }
        }
        if active == 0 {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    });
    let (cos0, sin0) = if kappa == 2 {
        (cospi(t.hi), sinpi(t.hi))
    } else {
        (0.0, 0.0)
    };
    let [sum, err, size] =
        classes
            .iter()
            .enumerate()
            .fold([0.0; 3], |[sum, err, size], (r, cls)| {
                let (s, s_err) = cls.best;
                let (c, u_err) = if kappa == 2 {
                    let c = if r == 0 {
                        cos0 - PI * t.lo * sin0
                    } else {
                        sin0 + PI * t.lo * cos0
                    };
                    (c, 0.0)
                } else {
                    let u = t.hi - r as f64 / kf;
                    (cospi(u) - PI * t.lo * sinpi(u), 2.0 * f64::EPSILON)
                };
                let rounding = (dtheta + u_err) * s.abs() + 3.0 * f64::EPSILON * (s * c).abs();
                [sum + s * c, err + s_err * c.abs() + rounding, size + s.abs()]
            });
    // the other exponentials of DLMF 16.11.9, with x e^(+-3 pi i / kappa), relative to these
    let next = 2.0 * (x * (cospi(3.0 / kf) - cos.hi)).exp();
    let others = if kappa == 2 { 0.0 } else { next };
    ((2.0 * sum, 2.0 * (err + others * size)), cos)
}

/// The exponential part of the asymptotic expansion of pFq for `p <= q` (DLMF 16.11.1),
/// with `|z|^nu` instead of `z^nu`:
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
    let kappa = b.len() + 1 - a.len();
    let kf = kappa as f64;
    let nu = sum_dd(a) - sum_dd(b) + 0.5 * (kf - 1.0);
    let real = z > 0.0 || kappa == 1;
    let root = root_dd(z.abs(), kappa) * if real { z.signum() } else { 1.0 };
    let x = kf * root.hi;
    // (the factors of the last coefficient are negative from here on)
    let beta = b.iter().fold(0.0_f64, |m, bj| m.max(kf * (bj - 1.0)));
    let settled = (nu.hi + b.len() as f64 + beta).max(0.0);
    let terms = pfq_exp_terms(a, b, nu.hi, x);
    let (series, growth) = if real {
        let (sum, err) = asymptotic_sum(terms, settled, x, |_| 0.0);
        // the other exponentials of DLMF 16.11.8 and 16.11.9, with x e^(+-2 pi i / kappa),
        // relative to this one
        let next = 2.0 * (x * (cospi(2.0 / kf) - 1.0)).exp();
        let others = if kappa == 1 { 0.0 } else { next };
        ((sum, err + others * sum.abs()), Dd::ONE)
    } else {
        pfq_exp_cos(terms, settled, kappa, root, nu)
    };
    // (2 pi)^((1 - kappa) / 2) kappa^(-1/2) |z|^(nu / kappa) e^x
    let ln_power = nu / kf * ln_dd(z.abs()) - ln_dd(kf) * 0.5;
    let l = ln_power - LN_SQRT_2PI_DD * (kf - 1.0) + root * kf * growth;
    let num = b.iter().copied().map(Dd::new);
    let den = a.iter().copied().map(Dd::new);
    gamma_series(num, den, series, l)
}

/// pFq for `p <= q` by its asymptotic expansion (DLMF 16.11.7 to 16.11.9)
///
/// Only the exponentially large part is used for `z > 0`, and only the algebraic part for
/// `z < 0` and `p = q`. The size of the other part is added to the error, for `z > 0` as the
/// sizes of the terms of the algebraic part. For `z < 0` and `p < q`, it is the sum of both
/// parts.
fn pfq_asymptotic(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let kappa = b.len() + 1 - a.len();
    let alg = pfq_alg(a, b, z, if z > 0.0 { Alg::Sizes } else { Alg::Value })?;
    let exp = pfq_exp(a, b, z)?;
    if z > 0.0 {
        // distinguish an empty algebraic part (zero error) from a failed `add_exp`
        let other = if alg.is_empty() {
            None
        } else {
            Some(add_exp(alg)?)
        };
        Some(Estimate::scaled(exp, other))
    } else if kappa == 1 {
        Some(Estimate::scaled(add_exp(alg)?, Some(exp)))
    } else {
        // (DLMF 16.11.9 for p <= q - 2 has the exponential part only; with the algebraic part,
        // it is theorem 3 of Paris, https://arxiv.org/abs/1708.04824)
        Some(Estimate::scaled(add_exp(alg.into_iter().chain([exp]))?, None))
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
fn pfq_alg_dd(a: &[f64], b: &[f64], z: f64) -> Option<Scaled> {
    add_exp_dd(pfq_alg(a, b, z, Alg::Dd)?)
}

/// pFq for `p > q + 1` and `z < 0` by `pfq_alg_dd`, or `None` if it may or may not overflow
fn pfq_borel_dd(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    Estimate::scaled_certain(pfq_alg_dd(a, b, z)?)
}

/// pFq for `p = q` and `z < 0` by the algebraic part of `pfq_asymptotic` in double-double, where
/// its terms cancel, e.g. for nearly equal `a_i`, or the terms of their series do
fn pfq_asymptotic_dd(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    if z > 0.0 || a.len() != b.len() {
        return None;
    }
    let alg = pfq_alg_dd(a, b, z)?;
    Some(Estimate::scaled(alg, Some(pfq_exp(a, b, z)?)))
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

/// pFq for `p > q + 1` and `z < 0`, where the series diverges: its Borel sum
///
/// That is the Meijer G function of DLMF 16.5.1, whose residues at the poles of the
/// `Gamma(a_i + s)` sum to `pfq_alg`. For small `|z|`, where those cancel, it is the series up
/// to its smallest term, as for 2F0, or else its Weniger transformation in double-double. The
/// residues come first, as the steps of Weniger can stall where some terms nearly vanish, e.g.
/// for an `a_i` near a non-positive integer. Where none of these is accurate, the residues are
/// in double-double, as for 2F0.
fn pfq_borel(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    // the term ratios are about n^d |z|, so that the terms are smallest near n = |z|^(-1/d)
    let d = (a.len() - b.len() - 1) as f64;
    // (capped as for 2F0)
    let scale = z.abs().powf(-1.0 / d).min(f64::EPSILON.recip());
    let terms = || once(1.0).chain(pfq_terms(a, b, z));
    let settled = pfq_settled(a, b).ceil();
    let series = || Some(Estimate::absolute(asymptotic_sum(terms(), settled, scale, |_| 0.0)));
    // also with the rounding errors of the terms, as for 2F0
    let weniger_dd =
        || Some(weniger_drift(pfq_weniger::<Dd>(&a, &b, z)?, terms(), settled, Dd::EPSILON));
    let residues = || Some(Estimate::scaled(add_exp(pfq_alg(a, b, z, Alg::Value)?)?, None));
    series()
        .or_better(residues)
        .or_better(weniger_dd)
        .or_better(|| pfq_borel_dd(a, b, z))
}

/// Cached evaluation of the Borel sum `x^a1 U(a1, a1 - a2 + 1, x)` for `2F0(a1, a2;; -1/x)`,
/// with its finite-difference rounding sensitivity and pole proximity flag
struct RawBorel {
    p: f64,
    u: f64,
    digits: f64,
    rounding: f64,
    near_pole: bool,
}

impl RawBorel {
    fn eval(a1: f64, a2: f64, x: f64, dx: f64) -> Option<Self> {
        const T: f64 = pow2(20);
        const POLE_SPACINGS: f64 = 8.0;
        let (d, e1) = two_sum(a1, -a2);
        let (b, e2) = two_sum(d, 1.0);
        // chgu takes U as a polynomial where its `a1 - b + 1` (that should be a2) rounds to
        // a non-positive integer
        let (s, f1) = two_sum(a1, -b);
        let (a2_rounded, f2) = two_sum(s, 1.0);
        if is_nonpos_int(a2_rounded) && a2_rounded != a2 {
            return None;
        }
        let (p, (u, _, digits)) = (x.powf(a1), xsf_hypu_digits(a1, b, x));
        if !(p.is_normal() && u.is_normal()) {
            return None;
        }
        // The rounding errors of b and x change y by about (y_t - y) / T, with b and x moved
        // by T times them (far enough to stand out from the error of U, and near enough to be
        // linear). That of chgu's `a1 - b + 1` adds to that of b, as chgu then returns about U
        // at `b = a1 - a2_rounded + 1`. A move that would round away is a spacing instead.
        let moved = |v: f64, dv: f64| {
            let v_t = v + T * dv;
            let spacing = (f64::EPSILON * v.abs()).max(f64::MIN_POSITIVE);
            if v_t == v && dv != 0.0 {
                v + spacing.copysign(dv)
            } else {
                v_t
            }
        };
        let db = ((e1 + e2).abs() + (f1 + f2).abs()).copysign(e1 + e2);
        let (b_t, x_t) = (moved(b, db), moved(x, dx));
        let (digits, rounding) = if b_t == b && x_t == x {
            (f64::from(digits), 0.0)
        } else {
            let (p_t, (u_t, _, digits_t)) = (x_t.powf(a1), xsf_hypu_digits(a1, b_t, x_t));
            if !(p_t.is_normal() && u_t.is_normal()) {
                return None;
            }
            // relative to y, also where it overflows
            let rounding = ((p_t / p) * u_t - u).abs() / u.abs() / T;
            (f64::from(digits.min(digits_t)), rounding)
        };
        // Within a few spacings of a non-positive integer n (also 0), the rounding of chgu's
        // a2 can move a cancellation in U, so only an accurate y is kept there (unless U is
        // the polynomial of a2 = n).
        let n = a2_rounded.round().min(0.0);
        let spacing = f64::EPSILON * a1.abs().max(b.abs()).max(1.0);
        let near = (a2 - n).abs().min((a2_rounded - n).abs()) <= POLE_SPACINGS * spacing;
        let near_pole = near && (a2 != n || a2_rounded != n);
        Some(Self {
            p,
            u,
            digits,
            rounding,
            near_pole,
        })
    }

    fn estimate(&self, margin: f64) -> Option<Estimate> {
        let y = self.p * self.u;
        let digits = self.digits - margin;
        // with the rounding doubled, as y can curve within the move, e.g. for an a1 near a
        // non-positive integer
        let err = 10.0_f64.powf(-digits) + 2.0 * self.rounding;
        // if accurate, U has the digits of `PFQ_TOL` plus the margin, and the rounding moves
        // y by at most `PFQ_TOL / 2`; an infinite y overflows by more than `err`
        let accurate = digits >= -PFQ_TOL.log10()
            && if y.is_infinite() {
                (self.p / (1.0 + err) * self.u).is_infinite()
            } else {
                self.rounding <= PFQ_TOL / 2.0
            };
        let err = if accurate {
            PFQ_TOL
        } else if self.near_pole {
            return None;
        } else if y.is_infinite() {
            f64::NAN
        } else {
            err.max(PFQ_INACCURATE)
        };
        // with the spacing of a subnormal y
        Some(Estimate::rounded(y, err))
    }
}

/// Borel sum of `2F0(a1, a2;; z)` for `z < 0` (DLMF 13.6.21)
///
/// Tries `x^a1 U(a1, a1 - a2 + 1, x)` with `CHGU_DIGITS_MARGIN`, the truncated asymptotic
/// series, the swapped `U` orientation (DLMF 13.2.40), double-double Weniger, double-double
/// residues (`pfq_borel_dd`), and finally the cached `U` evaluations at `margin = 0.0`.
fn hyp2f0_borel(a1: f64, a2: f64, z: f64) -> Option<Estimate> {
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
    let raw1 = RawBorel::eval(a1, a2, x, dx);
    let raw2 = OnceCell::new();
    let eval2 = || {
        (a1 != a2)
            .then(|| raw2.get_or_init(|| RawBorel::eval(a2, a1, x, dx)).as_ref())
            .flatten()
    };
    let a = [a1, a2];
    let terms = || once(1.0).chain(pfq_terms(&a, &[], z));
    // keep the 2-factor asymptotic terms split in `(mantissa, exponent)` so an underflowing or
    // subnormal term (e.g. tiny `a1 * a2 * z`) stays a conservative bound rather than ending the
    // series before `t_1`
    let (zm, ze) = frexp_f64(z);
    let series_terms = || {
        (0_u32..)
            .map(f64::from)
            .scan((1.0_f64, 0_i64), move |(m, e), n| {
                let t = ldexp_i64(*m, *e);
                let bound = if *m == 0.0 {
                    0.0
                } else if t.is_finite() {
                    t.abs().max(f64::MIN_POSITIVE)
                } else {
                    f64::INFINITY
                };
                let ((m1, e1), (m2, e2)) = (frexp_f64(a1 + n), frexp_f64(a2 + n));
                let (nm, ne) = frexp_f64(*m * m1 * m2 * zm / (n + 1.0));
                *m = nm;
                *e += i64::from(e1) + i64::from(e2) + i64::from(ze) + i64::from(ne);
                Some((t, bound))
            })
    };
    let settled = (-a1).max(-a2).max(0.0).ceil();
    // with x capped, as sqrt(2 pi x) in the error estimate could overflow; it stays above
    // 1e8, while past `settled` the terms alternate, so the first omitted one bounds the
    // remainder
    let scale = x.min(f64::EPSILON.recip());
    let series =
        || Some(Estimate::absolute(asymptotic_sum(series_terms(), settled, scale, |_| 0.0)));
    let weniger_dd =
        || Some(weniger_drift(pfq_weniger::<Dd>(&a, &[], z)?, terms(), settled, Dd::EPSILON));
    raw1.as_ref()
        .and_then(|r| r.estimate(CHGU_DIGITS_MARGIN))
        .or_better(series)
        .or_better(|| eval2()?.estimate(CHGU_DIGITS_MARGIN))
        .or_better(weniger_dd)
        .or_better(|| pfq_borel_dd(&a, &[], z))
        .or_better(|| raw1.as_ref()?.estimate(0.0))
        .or_better(|| eval2()?.estimate(0.0))
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
    let polynomial = degree.is_finite().then(|| {
        let (y, err, _) = pfq_sum::<Dd>(&a, &b, z, PFQ_DD_MAXITER);
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
        (&[a1, a2], &[]) if z < 0.0 => hyp2f0_borel(a1, a2, z),
        (&[a], &[b1, b2]) => hyp1f2_estimate(a, b1, b2, z),
        (&[a1, a2], &[b1, b2]) => hyp2f2_estimate(a1, a2, b1, b2, z),
        (&[a1, a2, a3], &[b1, b2]) => hyp3f2_estimate(a1, a2, a3, b1, b2, z),
        _ if a.len() == b.len() + 1 => pfq_real(a, b, z),
        _ if a.len() <= b.len() => {
            // the terms of the series for z < 0 reach about e^x, and the result
            // e^(x cos(pi / kappa)), or a power of z for p = q
            let kappa = (b.len() + 1 - a.len()) as f64;
            let x = kappa * z.abs().powf(1.0 / kappa);
            let cancellation = x * (1.0 - cospi(1.0 / kappa).max(0.0));
            pfq_entire(&a, &b, z, z > 0.0 || cancellation < PFQ_ASYMPTOTIC_FIRST)
        },
        _ if z < 0.0 => pfq_borel(&a, &b, z),
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
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-1.0], &[], 1.0, 2, [1e300].into_iter(), 1e-13);
        assert_eq!(y, 1e300);
        assert!(err.is_finite() && tails.is_finite());
    }

    #[test]
    fn test_pfq_ratio_bound_safe() {
        // huge denominator factors with huge `z`: exact ratio is 2^-7, not 0 from `inf` division
        let rho = pfq_ratio_bound(&[], &[pow2(520), pow2(510)], pow2(1023), 0.0);
        assert!(rho >= pow2(-7) && rho < 1.0, "{rho}");
        let rho = pfq_ratio_bound_full(&[], &[pow2(520), pow2(510)], pow2(1023), 0.0);
        assert!(rho >= pow2(-7) && rho < 1.0, "{rho}");
        // full bound at n > 0 does not insert an extra (1 + n) factor (exact 1/3, not 1/6)
        let rho = pfq_ratio_bound_full(&[], &[2.0], 1.0, 1.0);
        assert!((1.0 / 3.0..0.34).contains(&rho), "{rho}");
        // exact ratio exceeds 1 by ~4.4e-17 while the unpadded f64 product rounds below 1
        let b = 5.0 * pow2(78);
        let a = b + 3.0 * pow2(28);
        let z = -(1.0 + 7.0 * f64::EPSILON);
        let rho = pfq_ratio_bound(&[a, a, a, 1.0], &[b, b, b], 1.0 / z, 9999.0);
        assert!(rho >= 1.0, "{rho}");
        // exact bound of 1 when every a_i <= b_i (including implicit 1) at |z| = 1
        assert_eq!(pfq_ratio_bound(&[1.0, 1.0, 1.0], &[1e300, 2.0], -1.0, 1.0), 1.0);
        // implicit 1 is inserted in descending order: [0.75] -> [1.0, 0.75] pairs with [1.5, 0.5]
        let rho = pfq_ratio_bound(&[1.5, 0.5], &[0.75], -1.0, 1.0);
        assert!((1.25..1.3).contains(&rho), "{rho}");
        // subnormal scaling rounds upward, never down or to zero
        let u = f64::from_bits(1);
        assert!(pfq_ratio_bound(&[], &[4.0], 5.0 * u, 0.0) >= 2.0 * u);
        assert!(pfq_ratio_bound(&[], &[4.0], u, 0.0) >= u);
        // reciprocal bound rounds upward even on the exact-z fast path
        assert!(pfq_ratio_bound(&[1.0], &[], recip_up(-3.0), 0.0) > 1.0 / 3.0);
        assert!(pfq_ratio_bound_full(&[1.0], &[2.0], recip_up(-3.0), 0.0) > 1.0 / 3.0);
        // before settled (a_i + n < 0), the bound is infinite
        assert_eq!(pfq_ratio_bound(&[-2.0], &[], 0.5, 0.0), f64::INFINITY);
    }

    #[test]
    fn test_pfq_sum_flat_compensated() {
        // 1000 terms of 1e-17 each: preserved by Neumaier compensation, with nonzero tail norm
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-1000.0], &[], 1.0, 1001, core::iter::repeat(1e-17), 1e-13);
        assert!((y - (1.0 + 1e-14)).abs() <= 2.0 * f64::EPSILON, "{y}");
        assert!(tails > 1e-13 && err.is_finite(), "{tails}, {err}");
        // cancelling partial sums leaving a tiny final offset: tail norm stays finite near 1
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-3.0], &[], 1.0, 4, [1.0, -1.0, 1e-200].into_iter(), 1e-13);
        assert_eq!(y, 1.0);
        assert!(err.is_finite() && (tails - 1.0).abs() <= 4.0 * f64::EPSILON, "{err}, {tails}");
        let (y, err, tails) = pfq_sum_terms::<Dd>(
            &[-3.0],
            &[],
            1.0,
            4,
            [1.0, -1.0, 1e-200].into_iter().map(Dd::new),
            1e-13,
        );
        assert_eq!(y.hi, 1.0);
        assert!(err.is_finite() && (tails - 1.0).abs() <= 4.0 * f64::EPSILON, "{err}, {tails}");
        // non-finite partial sum bails out cleanly with infinite error
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-2.0], &[], 1.0, 3, [f64::MAX, f64::MAX].into_iter(), 1e-13);
        assert!(y.is_infinite() && err == f64::INFINITY && tails == f64::INFINITY);
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
            let y = pfq_borel_dd(&a, &[], z).unwrap();
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
        let (a, b) = ([pow2(-70), -19.75, -29.625], [1.0, -219.75, -229.625]);
        assert!(pfq_alg_series_dd(&a, &b, 0, -100.0).1.is_infinite());
        // the low parts of the parameters make the term after the one where they settle larger
        let (s, u) = (pow2(-49), pow2(-52));
        let (a, b) = ([s - 16.0, 1.0 + 3.0 * u], [1.0 + u, 2.0 * s - 16.0]);
        assert!(pfq_alg_series_dd(&a, &b, 0, -1.1 * s).1.is_infinite());
        // the bound of the unpaired factor underflows, and the terms grow to ~1e72 from the fifth
        let c = (1..=10).map(|i| f64::from(i) / 16.0);
        let a: Vec<_> = once(f64::from_bits(1)).chain(c.clone()).collect();
        let b: Vec<_> = once(3.0).chain(c.map(|c| c - pow2(49))).collect();
        assert!(pfq_alg_series_dd(&a, &b, 0, -2.0).1.is_infinite());
    }

    #[test]
    fn test_pfq_alg_series_dd_polynomial() {
        // of degree 0, which is below `|z|`
        let (sum, err) = pfq_alg_series_dd(&[0.5, 0.75], &[1.5, 0.25], 0, -0.5);
        assert!(sum.hi == 1.0 && err <= Dd::EPSILON);
    }

    #[test]
    fn test_pfq_asymptotic_integer_difference() {
        // the bound on the algebraic part is too large for equal a_i, and near them the
        // expansion is accurate or `None` (Arb reference)
        let c = 1.0 + pow2(-10);
        let (b, z) = ([c, c, 11.25], 21.25);
        let accurate = |a: &[f64; 3]| pfq_asymptotic(a, &b, z).and_then(Estimate::accurate);
        assert_eq!(accurate(&[1.0, 1.0, 10.5]), None);
        let near = accurate(&[1.0, 1.0 + pow2(-20), 10.5]);
        let expected = 739_413_030.005_544_4;
        assert!(near.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected));
    }

    #[test]
    fn test_cos_sin_pi_div_two() {
        let (cos, sin) = cos_sin_pi_div(2);
        assert_eq!((cos.hi, cos.lo, sin.hi, sin.lo), (0.0, 0.0, 1.0, 0.0));
    }

    #[test]
    fn test_pfq_next_rescaled_subnormal_quotient() {
        // extreme parameters whose ratio underflows while `ratio * z` is normal (checked mode)
        let next = pfq_next::<f64>(&[1e-160, 1e-160], &[], 0.0, 1e200, 1.0, true);
        assert!((next / 1e-120 - 1.0).abs() <= 4.0 * f64::EPSILON, "{next}");
        // smallest subnormal term is not double-rounded to zero by the rescaled path
        let a = [f64::from_bits(0x0010_0000_0000_0001)];
        let z = f64::from_bits(0x3c9f_ffff_ffff_ffff);
        assert_eq!(pfq_next::<f64>(&a, &[], 0.0, z, 1.0, true), f64::from_bits(1));
    }

    #[test]
    fn test_pfq_asymptotic_nearly_cancelled() {
        // equal a_i within sqrt(eps) of b_2 and b_3, where the algebraic part is about 1 and
        // changes much with a_i: accurate or `None` (Arb reference)
        let t = -pow2(-40);
        let (root, root_2) = (f64::EPSILON.sqrt(), (2.0 * f64::EPSILON).sqrt());
        let y = pfq_asymptotic(&[t; 3], &[t - 1.0, t + root, t + root_2], 30.0)
            .and_then(Estimate::accurate);
        let expected = -816_436.169_307_593_7_f64;
        assert!(y.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected.abs()));
    }

    #[test]
    fn test_pfq_alg_dd_params() {
        // exact terminating parameter (-2) beside one with identical hi word and nonzero lo
        // (-2 - 2^-52): degree is 2, not infinite or 0
        let [up, _] = residue::series_params(&[-2.0, 0.5], &[1.0 + pow2(-52)], 0);
        assert_eq!(up[0].hi, -2.0);
        assert_eq!(up[1].hi, -2.0);
        assert_ne!(up[1].lo, 0.0);
        assert_eq!(first_zero_dd(&up), 2.0);
        let (sum, err) = pfq_alg_series(&[-2.0, 0.5], &[1.0 + pow2(-52)], 0, -10.0, false);
        // 1 + (-2)(-2 - 2^-52)/(1 * -1.5) * (-0.1) + ...
        let expected = 1.0 + 4.0 / 15.0 + 4.0 / 150.0;
        assert!((sum - expected).abs() <= 1e-13 && err <= 1e-13, "{sum}, {err}");
        // parameter with hi = -3 and nonzero lo = 2^-60 does not falsely terminate
        let [up_nonpole, _] = residue::series_params(&[pow2(-60), 0.5], &[4.0], 0);
        assert_eq!(up_nonpole[1].hi, -3.0);
        assert_eq!(up_nonpole[1].lo, pow2(-60));
        assert!(first_zero_dd(&up_nonpole).is_infinite());
        // |a_m| >> 1 with v ~ a_m preserves +1 in (a_m - v) + 1
        let big = pow2(60);
        let [up_big, _] = residue::series_params(&[big, 0.5], &[big - 256.0], 0);
        assert_eq!(up_big[1], Dd::new(257.0));
        // genuine b_j = 1.0 with a_j = 2^-60 keeps exact 2^-60 in up (not rounded to 0)
        let tiny = pow2(-60);
        let [up_one, _] = residue::series_params(&[tiny, 0.5], &[1.0], 0);
        assert_eq!(up_one[0], Dd::new(tiny));
        assert_eq!(up_one[1], Dd::new(tiny));
        let (sum, err) = pfq_alg_series(&[tiny, 0.5], &[1.0], 0, -10.0, false);
        assert!(
            (sum - 1.0).abs() <= 4.0 * f64::EPSILON && err <= 4.0 * f64::EPSILON,
            "{sum}, {err}"
        );
        // exact termination at huge |z| does not overflow `TAU * |z|`
        let (sum, err) = pfq_alg_series(&[-1.0], &[1.0, 8.0], 0, -f64::MAX, false);
        assert_eq!(sum, 1.0);
        assert!(err.is_finite() && err <= 4.0 * f64::EPSILON, "{err}");
        // regrowth guard rounds conservatively on threshold pull = 1 + eps
        let mut scratch = BoundScratch::new(2, 1);
        let rho = scratch.regrowth_bound(
            &[Dd::new(3.0), Dd::new(3.0)],
            &[Dd::new(1.0)],
            next_up_nonneg(1.0 / 12.0),
            0.0,
        );
        assert_eq!(rho, f64::INFINITY);
        // six-factor threshold where `1/5` rounding in `powf` would exceed a small fixed margin
        let mut scratch = BoundScratch::new(6, 1);
        let rho = scratch.regrowth_bound(
            &[Dd::new(8.590_444_340_720_38e49); 6],
            &[Dd::new(1.0)],
            next_up_nonneg(1e-300),
            0.0,
        );
        assert_eq!(rho, f64::INFINITY);
        // repeated-pole shifted-factor case: early f64 rounding drops the combined bound below 1,
        // while preserving Dd low parts keeps it strictly above 1
        let am = 1.0 + pow2(-41);
        let a = [am, am, 8191.5 + pow2(-40)];
        let b = [8191.5, 1.0];
        let z = -5.002_442_002_406_743;
        let residues = Residues::of_class(&a, &b, 0, z);
        let [up_r, low_r] = residues.params();
        let mut scratch = BoundScratch::new(up_r.len(), low_r.len());
        let bound = scratch.ratio_bound(&up_r, &low_r, recip_up(z), 8191.0)
            * Residues::growth_bound(&up_r, &low_r, 8191.0);
        assert!(bound > 1.0, "{bound}");
        // terminating algebraic series with many parameters accounts for term-generation rounding
        let mut b_many = vec![2.5];
        b_many.resize(51, pow2(-53));
        let (sum, err) = pfq_alg_series(&[0.5], &b_many, 0, pow2(20), false);
        let expected = 305.041_624_171_279_54;
        assert!(
            err.is_finite() && (sum - expected).abs() <= err,
            "{sum} vs {expected} (err {err})"
        );
        // huge finite degree-1 term does not overflow error accumulation before scaling by eps
        let (sum_big, err_big) = pfq_alg_series(&[0.5, 1.0], &[2.5, pow2(1023)], 0, -2.0, false);
        assert!(sum_big.is_finite() && err_big.is_finite(), "{sum_big}, {err_big}");
        // subnormal scaled denominator in double-double rejects the lossy term
        let (_, err_sub) =
            pfq_alg_series_dd(&[0.5, 0.75], &[2.5, 1.5 + pow2(-52)], 0, -3.0 * f64::from_bits(1));
        assert!(err_sub.is_infinite(), "{err_sub}");
    }

    #[test]
    fn test_raw_borel_exact_and_margin() {
        // exact dyadic parameters (db == 0, dx == 0) skip perturbed evaluation with zero rounding
        let raw = RawBorel::eval(0.5, 1.5, 4.0, 0.0).unwrap();
        assert_eq!(raw.rounding, 0.0);
        assert!(raw.estimate(CHGU_DIGITS_MARGIN).unwrap().is_accurate());
        // 13-digit cached evaluation is inaccurate (still trusted) at margin 1, accurate at margin 0
        let raw13 = RawBorel {
            p: 1.0,
            u: 2.0,
            digits: 13.0,
            rounding: 0.0,
            near_pole: false,
        };
        let strict = raw13.estimate(CHGU_DIGITS_MARGIN).unwrap();
        assert!(!strict.is_accurate() && strict.trusted().is_some());
        assert!(raw13.estimate(0.0).unwrap().is_accurate());
        // near-pole gate rejects inaccurate estimates while keeping accurate ones
        let raw_pole = RawBorel {
            near_pole: true,
            ..raw13
        };
        assert!(raw_pole.estimate(CHGU_DIGITS_MARGIN).is_none());
        assert!(raw_pole.estimate(0.0).unwrap().is_accurate());
        // overflow gate: overflows comfortably vs barely past f64::MAX
        let raw_over = RawBorel {
            p: 1e200,
            u: 1e200,
            digits: 13.0,
            rounding: 0.0,
            near_pole: false,
        };
        assert!(
            raw_over
                .estimate(CHGU_DIGITS_MARGIN)
                .unwrap()
                .trusted()
                .is_none()
        );
        assert_eq!(raw_over.estimate(0.0).unwrap().accurate(), Some(f64::INFINITY));
        let raw_edge = RawBorel {
            p: f64::MAX,
            u: 1.0 + 1e-14,
            digits: 13.0,
            rounding: 0.0,
            near_pole: false,
        };
        assert!(raw_edge.estimate(0.0).unwrap().trusted().is_none());
    }

    #[test]
    fn test_pfq_exp_cos_single_pass() {
        // kappa = 2: class 0 has a zero followed by a nonzero `k + kappa` lookahead term,
        // while class 1 decays and stops independently; nu = 0 gives phase theta = 200 rad
        // and cos(pi / 2) = 0; verify the underlying iterator is never polled past its first None
        let mut inner = [1.0, 0.5, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0].into_iter();
        let mut exhausted = false;
        let seq = from_fn(move || {
            assert!(!exhausted, "polled past first None");
            let next = inner.next();
            if next.is_none() {
                exhausted = true;
            }
            next
        });
        let ((sum, err), cos) = pfq_exp_cos(seq, 0.0, 2, Dd::new(100.0), Dd::ZERO);
        assert_eq!(cos, Dd::ZERO);
        // class 0 exact alternating sum is 1.0 - 0.0 + 0.25 = 1.25; must not stop at the 0.0
        let expected = 2.0 * 1.25 * 200.0_f64.cos() + 2.0 * 0.5 * 200.0_f64.sin();
        assert!((sum - expected).abs() <= err && err <= 1e-10, "{sum} vs {expected} (err {err})");
    }
}
