//! Generalized hypergeometric functions pFq of real arguments.
//!
//! The internal functions are generic over the numbers of parameters. They take the lower
//! parameters together with the 1 of `n!` as `b`, so that the terms of the series are
//! `t_n = z^n prod_i (a_i)_n / prod_j (b_j)_n`.

use core::f64::consts::{PI, TAU};
use core::iter::{from_fn, once, successors, zip};
use core::ops::Add;

use super::{
    Scaled, add_exp, fold_exp, ln_dd, ln_gamma, mul_gamma_ratio, pole_distance, trigamma_lo,
};
use crate::utils::{EPSILON_SQ, is_nonpos_int, neg_one_pow, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{
    DoubleDouble as Dd, FRAC_2_PI_DD, LN_PI_DD, LN_SQRT_2PI_DD, frexp_f64, ldexp_f64, two_sum,
};
use crate::xsf::cephes::{gamma, gammasgn, lgam, rgamma};
use crate::xsf::{cospi, sinpi, xsf_hypu_digits};

mod conformal;
mod residue;

use residue::Residues;

/// Maximum number of terms of `pfq_series`
const PFQ_MAXITER: u32 = 10_000;

/// Maximum number of terms of `pfq_series` in double-double arithmetic, fewer as a term takes
/// about 20 times as long
const PFQ_DD_MAXITER: u32 = PFQ_MAXITER / 4;

/// Largest estimated relative error of a result
const PFQ_TOL: f64 = 1e-13;

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

/// The `|x|` below which `ln_gamma` is by `lgam`, so that `gamma` and `rgamma` can be more accurate
const PFQ_GAMMA_MAX: f64 = 13.0;

/// Bound of the relative error of `gamma` and `rgamma` for `|x| < PFQ_GAMMA_MAX`, with the
/// `1.5 eps` of the products of `gamma_pow` and `mul_gammas`: they are within `4.4 eps` at 1.3e7
/// points (measured)
const PFQ_GAMMA_ERR: f64 = 7.0 * f64::EPSILON;

/// Maximum order of `pfq_weniger`
const PFQ_WENIGER_MAXITER: usize = 40;

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

/// The products of the `a_i + n` and of the `b_j + n`, whose quotient is `t_{n+1} / (z t_n)`
///
/// They and their quotient are normal for up to `PFQ_MAX_UNCHECKED` parameters with sizes in
/// `PFQ_MODERATE`, as each nonzero `|a_i + n|` is then in `[2^-100, 2^101]`.
fn pfq_products<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> [T; 2] {
    let shifted = |p: T, &v: &f64| p * (T::from_f64(v) + n);
    let one = T::from_f64(1.0);
    [a.iter().fold(one, shifted), b.iter().fold(one, shifted)]
}

/// `t_{n+1} / (z t_n)` for any parameters and `n`, which is zero if it underflows, and infinite
/// if it overflows or would be inexact, as it is not normal
#[cold]
fn pfq_checked_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> T {
    // normal in the arithmetic of `T`
    let is_normal = |v: T| (T::MIN_POSITIVE..=f64::MAX).contains(&v.to_f64().abs());
    // as mantissa and `i64` exponent, so that no partial product over- or underflows
    let product = |params: &[f64]| {
        let times = |(m, e): (T, i64), &v: &f64| {
            let (vm, ve) = (T::from_f64(v) + n).frexp();
            let (m, me) = (m * vm).frexp();
            (m, e + i64::from(ve) + i64::from(me))
        };
        params.iter().fold((T::from_f64(1.0), 0), times)
    };
    let ((num, ne), (den, de)) = (product(a), product(b));
    // beyond `i32`, it over- or underflows anyway
    let ratio = (num / den).ldexp((ne - de).clamp(i32::MIN.into(), i32::MAX.into()) as i32);
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

/// The terms `t_1, t_2, ...` of the series after `t_0 = 1`, in the arithmetic of `T`
///
/// A term that underflows is a bound on it. The terms end before one that overflows or is NaN,
/// or that is such a bound and not negligible.
fn pfq_terms<T: SeriesReal>(a: &[f64], b: &[f64], z: f64) -> impl Iterator<Item = T> {
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let checked = pfq_extreme(a, b);
    // the largest term, which is small once the terms are only upper bounds
    let mut limit = f64::MAX;
    let indices = (0_u32..).map(f64::from);
    indices.scan(T::from_f64(1.0), move |term, n| {
        let ratio: T = pfq_ratio(a, b, n, checked);
        let next = *term * (ratio * T::from_f64(z));
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

/// `pfq_series`, and the root of the sum of its squared tails `s - s_n`, which is at least the
/// largest of them
fn pfq_series_tails(a: &impl AsRef<[f64]>, b: &impl AsRef<[f64]>, z: f64) -> (f64, f64, f64) {
    let (y, err, tails) = pfq_sum::<f64>(a, b, z, PFQ_MAXITER);
    // (not in double-double, whose `EPSILON` leaves room for many more parameters)
    let err = err * pfq_rounding_scale(a.as_ref(), b.as_ref());
    if err <= PFQ_TOL * y.abs() || err == f64::INFINITY {
        (y, err, tails)
    } else {
        pfq_sum::<Dd>(a, b, z, PFQ_DD_MAXITER)
    }
}

/// `pfq_series_tails` in the arithmetic of `T`, with at most `maxiter` terms
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
) -> (f64, f64, f64) {
    // generic over the storage of the parameters, as the sum is faster if it is compiled for
    // the numbers of parameters of its caller
    let (a, b) = (a.as_ref(), b.as_ref());
    let settled = pfq_settled(a, b);
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let mut sum = T::from_f64(1.0);
    // the partial sum as `f64`
    let mut s = sum.to_f64();
    // the sums of s_n and s_n^2
    let (mut sum1, mut sum2) = (0.0, 0.0_f64);
    // the rounding error, and the root of the sum of (s - s_n)^2
    let estimates = |n: f64, sum: f64, sum1: f64, sum2: f64| {
        // accurate enough for an estimate, relative to s so that its parts don't overflow
        let tails = if sum == 0.0 {
            sum2.sqrt()
        } else {
            ((n + 1.0) - 2.0 * (sum1 / sum) + sum2 / sum / sum)
                .max(0.0)
                .sqrt()
                * sum.abs()
        };
        let rounding = T::EPSILON * ((4.0 + 0.25 * (n + 1.0).sqrt()) * tails + sum2.sqrt());
        (rounding, tails)
    };
    let mut terms = pfq_terms::<T>(a, b, z);
    for n in (0..maxiter).map(f64::from) {
        sum1 += s;
        sum2 += s * s;
        if n == degree {
            // the series terminated
            let (rounding, tails) = estimates(n, s, sum1, sum2);
            return (s, rounding, tails);
        }
        // t_(n+1)
        let Some(term) = terms.next() else { break };
        let t = term.to_f64();
        sum = sum + term;
        s = sum.to_f64();
        if n >= settled && t.abs() <= f64::EPSILON * s.abs() {
            let rho = pfq_ratio_bound(a, b, z, n + 1.0);
            // geometric bound on the remaining terms, if the series converges (p <= q + 1)
            let tail = t.abs() * rho / (1.0 - rho);
            if a.len() <= b.len() && rho < 1.0 && tail <= f64::EPSILON * s.abs() {
                let (rounding, tails) = estimates(n, s, sum1, sum2);
                return (s, rounding + tail, tails);
            }
        }
    }
    (s, f64::INFINITY, f64::INFINITY)
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
/// The estimate is infinite if the orders up to `PFQ_WENIGER_MAXITER` did not converge, or if a
/// term overflows, or underflows before the terms become large again. It is `None` if the terms
/// only alternate after more than `PFQ_MAXITER` terms, or if a term before that overflows or
/// underflows.
///
/// 1. E. J. Weniger, "Nonlinear sequence transformations for the acceleration of convergence
///    and the summation of divergent series", Comput. Phys. Rep. 10 (1989), 189-371
/// 2. R. M. Slevinsky, "Fast and stable rational approximation of generalized hypergeometric
///    functions", Numer. Algor. 98 (2025), 587-624, <https://doi.org/10.1007/s11075-024-01808-w>
pub(crate) fn pfq_weniger(
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
    let ratio = |n| -> f64 { pfq_ratio(a, b, n, checked) };
    // the rounding errors of the partial sums and of the terms (the factor is measured)
    let eps = 4.0 * f64::EPSILON * pfq_rounding_scale(a, b);
    let (mut term, mut sum) = (1.0_f64, 1.0_f64);
    // the largest |s_n| before `start`
    let mut scale = 1.0_f64;
    for n in (0..start as u32).map(f64::from) {
        term *= ratio(n) * z;
        if !term.is_normal() {
            // overflow or underflow
            return None;
        }
        sum += term;
        scale = scale.max(sum.abs());
    }
    if start == degree {
        // the series terminated
        return Some((sum, eps * scale));
    }
    let mut rows = [[0.0; 3]; PFQ_WENIGER_MAXITER];
    let (mut delta, mut step) = (sum, f64::INFINITY);
    for k in 0..PFQ_WENIGER_MAXITER {
        let n = start + k as f64;
        term *= ratio(n) * z;
        if !term.is_normal() {
            // overflow, or underflow, after which the terms only matter if they can increase
            let negligible =
                term.abs() < 1.0 && a.len() <= b.len() && pfq_ratio_bound(a, b, z, n) <= 1.0;
            let noise = eps * scale;
            return Some((sum, if negligible { noise } else { f64::INFINITY }));
        }
        // the numerator, the denominator, and the numerator for |s_j|
        rows[k] = [sum / term, 1.0 / term, sum.abs() / term];
        for j in (0..k).rev() {
            // the order of `rows[j]`, which starts at s_j
            let (m, i) = (j as f64, (k - j - 1) as f64);
            let c = if i == 0.0 {
                1.0
            } else {
                let (s, t) = (m + i, m + 2.0 * i);
                (s + 1.0) * s / ((t + 1.0) * t)
            };
            rows[j] = core::array::from_fn(|l| rows[j + 1][l] - c * rows[j][l]);
        }
        let [num, den, abs] = rows[0];
        let noise = eps * ((abs / den).abs() + (num / den).abs() + scale);
        let last = step;
        step = (num / den - delta).abs();
        delta = num / den;
        // delta is a weighted mean of the last one and of the one of that order from s_1
        // (NaN for k = 0). For a large |z| the weight of the latter is small, and so are the
        // steps, also while the two still differ, so it has to agree as well.
        let other = rows[1][0] / rows[1][1];
        if step <= noise && last <= noise && (other - delta).abs() <= noise {
            return Some((delta, noise));
        }
        sum += term;
    }
    Some((delta, f64::INFINITY))
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
        pfq_alg(a, &b, z, false, &mut terms)?;
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

/// `_{q+1}F_q` for `z < 1`, or for a polynomial, by the first method for `z` that is accurate
///
/// `b` contains the 1 of `n!`. Unless it is a polynomial, the last method is `pfq_conformal`,
/// which is the slowest.
pub(crate) fn pfq_accurate<P: Params>(a: &P, b: &P, z: f64) -> Option<f64> {
    let series: Method<P> = |a, b, z| Some(pfq_series(a, b, z));
    let polynomial = first_zero(a.as_ref()).is_finite();
    let methods: &[Method<P>] = if polynomial || z > PFQ_WENIGER_MAX_Z {
        &[series]
    } else if z > -1.0 {
        &[pfq_weniger, series]
    } else if z == -1.0 {
        &[pfq_weniger]
    } else if z >= PFQ_WENIGER_FIRST_Z {
        &[pfq_weniger, pfq_continuation]
    } else {
        &[pfq_continuation, pfq_weniger]
    };
    let last: Option<Method<P>> = (!polynomial).then_some(conformal::pfq_conformal);
    let mut methods = methods.iter().chain(&last);
    methods.find_map(|method| accurate(method(a, b, z)?))
}

/// The result `y`, if it is finite and its estimated absolute error is at most `PFQ_TOL |y|`
fn accurate((y, err): (f64, f64)) -> Option<f64> {
    (y.is_finite() && err <= PFQ_TOL * y.abs()).then_some(y)
}

/// `pfq_series` if it is `accurate`, after sorting the parameters, which tightens
/// `pfq_ratio_bound`
fn pfq_series_sorted(a: impl Params, b: impl Params, z: f64) -> Option<f64> {
    accurate(pfq_series(&descending(a), &descending(b), z))
}

/// The sum of `v` in double-double
fn sum_dd(v: &[f64]) -> Dd {
    v.iter().fold(Dd::new(0.0), |s, &v| s + v)
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

/// `_{q+1}F_q` for finite arguments, where the last of `b` is the 1 of `n!`
fn pfq_real<P: Params>(a: P, b: P, z: f64) -> Option<f64> {
    if z == 0.0 {
        return Some(1.0);
    }
    let Some(degree) = pfq_degree(a.as_ref(), b.as_ref()) else {
        return Some(f64::INFINITY);
    };
    if degree.is_infinite() && z > 1.0 {
        // complex
        Some(f64::NAN)
    } else if degree.is_infinite() && z == 1.0 {
        // +inf where the series diverges, and no method where it converges
        diverges_at_one(a.as_ref(), b.as_ref().split_last()?.1)?.then_some(f64::INFINITY)
    } else {
        pfq_accurate(&descending(a), &descending(b), z)
    }
}

/// Generalized hypergeometric function $_3F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Notes
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
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)). For
/// $z < -1$, it is the analytic continuation [^1] [^2]. A non-positive integer $a_i = -n$ makes it
/// a polynomial of degree at most $n$, defined for all $z$.
///
/// As in [`hyp2f1`](crate::hyp2f1), an $a_i$ equal to a $b_j$ cancels it, also if it is a
/// non-positive integer, so that the result is the $_2F_1$ of the other parameters. Where that
/// converges at $z = 1$, it is Gauss's sum [^3]
///
/// $$
/// \hyp 2 1 {a_1\enspace a_2} b 1
/// = {\Gamma(b) \ \Gamma(b - a_1 - a_2) \over \Gamma(b - a_1) \ \Gamma(b - a_2)}
/// .
/// $$
///
/// It is $+\infty$:
///
/// - if a $b_j$ is a non-positive integer $-m$, unless it is a polynomial of degree $\lt m$;
/// - at $z = 1$ if $b_1 + b_2 \le a_1 + a_2 + a_3$, unless it is a polynomial.
///
/// It is NaN for $z > 1$, where it is complex, unless it is a polynomial, and if an argument is
/// NaN. For an infinite argument, it is the limit as in [`hyppfq`].
///
/// # Returns `None`
///
/// If the relative error could exceed about $10^{-13}$:
///
/// - at $z = 1$, unless a parameter cancels so that Gauss's sum applies, and for
///   $1 - 10^{-5} \lesssim z < 1$;
/// - if the terms cancel too much (e.g. large parameters), overflow or underflow;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.8.E8>
/// [^3]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/15.4.E20>
///
/// # See also
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
#[must_use]
pub fn hyp3f2(a1: f64, a2: f64, a3: f64, b1: f64, b2: f64, z: f64) -> Option<f64> {
    if ![a1, a2, a3, b1, b2, z].into_iter().all(f64::is_finite) {
        return hyppfq(&[a1, a2, a3], &[b1, b2], z);
    }
    let (a, b) = ([a1, a2, a3], [b1, b2, 1.0]);
    // a common parameter cancels: 2F1, or 1F0 if that has one as well
    let mut pairs = (0..3).flat_map(|i| (0..2).map(move |j| (i, j)));
    let Some((i, j)) = pairs.find(|&(i, j)| a[i] == b[j]) else {
        return pfq_real(a, b, z);
    };
    let (a, b) = ([a[(i + 1) % 3], a[(i + 2) % 3]], [b[1 - j], 1.0]);
    match a.iter().position(|&v| v == b[0]) {
        Some(i) => pfq_real([a[1 - i]], [1.0], z),
        None if z == 1.0 => gauss_sum(a[0], a[1], b[0]).or_else(|| pfq_real(a, b, z)),
        None => pfq_real(a, b, z),
    }
}

/// Gauss's sum `2F1(a, b; c; 1)` (DLMF 15.4.20), if it converges and its estimated relative
/// error is at most `PFQ_TOL`
fn gauss_sum(a: f64, b: f64, c: f64) -> Option<f64> {
    let c = Dd::new(c);
    let excess = c - a - b;
    if !(excess.hi > 0.0) {
        return None;
    }
    let y = mul_gamma_ratio([c, excess], [c - a, c - b], || Scaled::new(1.0, 0, 0.0))?;
    // its error bound can't tell if it overflows
    checked_f64(y, None).filter(|y| y.is_finite())
}

/// A term of a series with a bound on its size
trait Term: Copy {
    fn value(self) -> f64;
    fn bound(self) -> f64;
}

impl Term for f64 {
    fn value(self) -> f64 {
        self
    }

    fn bound(self) -> f64 {
        self.abs()
    }
}

/// A value and a bound on its size
impl Term for (f64, f64) {
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
/// negligible if `ratio_bound(n) < 1`, so that the next terms are smaller.
fn asymptotic_sum<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    z: f64,
    ratio_bound: impl Fn(f64) -> f64,
) -> (f64, f64) {
    let (mut sum, mut abs_sum, mut prev) = (0.0_f64, 0.0, f64::INFINITY);
    // the smallest size, and the sum before a term with the smallest error
    let (mut min, mut best) = (f64::INFINITY, (sum, f64::INFINITY));
    let remainder = (TAU * z.abs()).sqrt().max(1.0);
    let mut terms = terms.peekable();
    for n in (0..PFQ_MAXITER).map(f64::from) {
        let Some(term) = terms.next() else { break };
        let (t, bound) = (term.value(), term.bound());
        let size = bound.max(terms.peek().map_or(0.0, |next| next.bound()));
        if !t.is_finite() || n >= settled + z.abs() && bound >= prev {
            break;
        } else if n >= settled && size < min {
            min = size;
            let err = remainder * min + n * f64::EPSILON * abs_sum;
            if err < best.1 {
                best = (sum, err);
            }
            if remainder * min <= f64::EPSILON * sum.abs() && ratio_bound(n) < 1.0 {
                break;
            }
        }
        sum += t;
        abs_sum += bound;
        prev = bound;
    }
    best
}

/// The sum of a convergent series and an estimate of its absolute error, which is infinite if
/// it did not converge within `PFQ_MAXITER` terms
///
/// From term `settled` on, the bound on each term is at most `ratio_bound(n)` times the one
/// before it, and `ratio_bound` does not increase.
fn convergent_sum<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    ratio_bound: impl Fn(f64) -> f64,
) -> (f64, f64) {
    let (mut sum, mut abs_sum) = (0.0_f64, 0.0);
    for (term, n) in terms.zip((0..PFQ_MAXITER).map(f64::from)) {
        if !term.value().is_finite() {
            break;
        }
        sum += term.value();
        abs_sum += term.bound();
        let rho = if n >= settled { ratio_bound(n) } else { 1.0 };
        // geometric bound on the remaining terms
        let tail = term.bound() * rho / (1.0 - rho);
        if rho < 1.0 && tail <= f64::EPSILON * sum.abs() {
            return (sum, tail + n * f64::EPSILON * abs_sum);
        }
    }
    (sum, f64::INFINITY)
}

/// The sum of the series in `1 / z` of a term of the algebraic part, where the term ratios are
/// about `n^kappa / |z|`: by `asymptotic_sum`, or by `convergent_sum` for `kappa = 0`
fn pfq_alg_sum<T: Term>(
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

/// `Gamma(g)^power` for `power = +-1` and a bound of its relative error, if `|g| < PFQ_GAMMA_MAX`
/// and that bound is below `ln_err`, the one of `ln_gamma`
///
/// It is to first order in `g.lo`: `Gamma(x + e) = Gamma(x) exp(u + t)` with `u = psi(x) e`,
/// where `|t|` is at most the second-order term of `ln_gamma`.
fn gamma_pow(g: Dd, power: f64, ln_err: f64) -> Option<(f64, f64)> {
    if ln_err <= PFQ_GAMMA_ERR || g.hi.abs() >= PFQ_GAMMA_MAX {
        return None;
    }
    let f: fn(f64) -> f64 = if power > 0.0 { gamma } else { rgamma };
    let t = trigamma_lo(g, pole_distance(g)?) * g.lo.abs() / 2.0;
    let u = power * gamma_low((g.hi, g.lo));
    // |exp(u + t) - 1 - u| <= w^2 + |t| if w = |u| + |t| <= 1 (else `err >= ln_err`)
    let (v, w) = (f(g.hi) * (1.0 + u), u.abs() + t);
    let err = PFQ_GAMMA_ERR + w * w + t;
    (v.is_normal() && err < ln_err).then_some((v, err))
}

/// `y` times `Gamma(num_0) ... / (Gamma(den_0) ...)`, with the error bounds of the factors
/// added, or `None` at a pole
///
/// Each factor is by `ln_gamma`, whose error for `|g| < PFQ_GAMMA_MAX` is that of a logarithm
/// and grows with `|ln|Gamma(g)||`, or by `gamma_pow` if its error bound is smaller.
fn mul_gammas(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    mut y: Scaled,
) -> Option<Scaled> {
    let gammas = num.into_iter().map(|g| (g, 1.0));
    for (g, power) in gammas.chain(den.into_iter().map(|g| (g, -1.0))) {
        let (v, sign, v_err) = ln_gamma(g)?;
        if let Some((v, v_err)) = gamma_pow(g, power, v_err) {
            // (the product of the mantissas can't overflow)
            let ((m, e), (v, k)) = (frexp_f64(y.m), frexp_f64(v));
            (y.m, y.e, y.err) = (m * v, y.e + e + k, y.err + v_err);
        } else {
            (y.m, y.err) = (sign * y.m, y.err + v_err);
            y = y.mul_exp(v * power);
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
    let l = nu * ln_dd(Dd::new(z.abs())) + z;
    gamma_series(b.map(Dd::new), a.map(Dd::new), series, l)
}

/// `(x, k, factor)` for the rounded parameter `x = 1 + am - v` of a series, where `-k` is the
/// non-positive integer nearest to it: times `factor`, the terms after `t_k` have the factor
/// `x + k` of the exact parameter if it is an `upper` one, and `1 / (x + k)` if it is a lower one
///
/// The rounding error of `x` is large relative to a small `x + k`. An inexact upper `x` at `-k`
/// would end the series, so that it is a few ulp from `-k` instead.
fn shifted(am: f64, v: f64, upper: bool) -> (f64, f64, f64) {
    let (s, e1) = two_sum(1.0, am);
    let (nearest, e2) = two_sum(s, -v);
    let (k, e) = (-nearest.round().min(0.0), e1 + e2);
    let offset = nearest + k;
    let ends = upper && offset == 0.0 && e != 0.0;
    // (more than an ulp of -k, also for k = 0)
    let x = nearest - if ends { (k + 1.0) * f64::EPSILON } else { 0.0 };
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
    let bound = |n| pfq_ratio_bound(up, low, 1.0 / z, n);
    // kappa = q - p + 1, as `b` contains the 1 of n!
    let kappa = b.as_ref().len() - a.as_ref().len();
    // their k and factors, and the largest k with the factor of the terms after it
    let shifts = || up_shifts.iter().chain(low_shifts.iter());
    let (k_max, all) = shifts().fold((0.0_f64, 1.0), |(top, all), (&k, &f)| (k.max(top), all * f));
    let series = once(1.0)
        // in (-1)^(p+q+1) / z, where (-1)^(p+q+1) = (-1)^kappa
        .chain(pfq_terms(up, low, neg_one_pow(kappa) / z))
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

/// The term of `pfq_alg` of the class of `a_m`: the `a_i` that differ from `a_m` by integers
///
/// The terms of these `a_i` are infinite, unlike the limit of their sum, which is the sum of
/// the `Residues` at the poles that they share. It is the `term` of the first of the smallest
/// of them, with the residue `(-1)^k / k!` of the Gamma function in place of each `Gamma(-k)`
/// of its factor. The `term` stays `None` for the other `a_i`, and if the class has no residues.
///
/// The residues are complex for `z > 0`, where it is the sum of the bounds on their sizes, as
/// with `sizes`.
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
    let residues = Residues::of_class(a, lower, m, z);
    let [up, low] = residues.params().map(descending);
    // (their order is final after the last zero factor, and their bounds are valid once all
    // parameters plus n exceed 1)
    let settled = pfq_settled(&up, &low).floor() + 2.0;
    let bound = |n| pfq_ratio_bound(&up, &low, 1.0 / z, n) * Residues::growth_bound(&up, &low, n);
    let sizes = sizes || z > 0.0;
    let series = residues.map(|(t, size)| (if sizes { size } else { t }, size));
    let (sum, err) = pfq_alg_sum(series, settled, b.len() - a.len(), z, bound);
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
        let l = -ln_dd(Dd::new(z.abs())) * am;
        let mut y = gamma_series(num, den, (sign * sum, err), l)?;
        if sizes {
            y.m = y.m.abs();
        }
        *term = Some(y);
    }
    Some(())
}

/// The terms of the algebraic part of the asymptotic expansion of pFq for `p <= q` (DLMF
/// 16.11.2), with `|z|` instead of `-z`, or `None` if a factor of a term overflows or underflows
///
/// `b` contains the 1 of `n!` as its last parameter. The term of `a_m` is `|z|^-a_m` times
/// a `_{q+1}F_{p-1}` series in `(-1)^(p+q+1) / z`, with the parameters `a_m` and
/// `1 + a_m - b_j`, and `1 + a_m - a_i`. It stays `None` if it vanishes, where
/// `1 / Gamma(b_j - a_m)` does. The `a_i` that differ by integers have one term, by
/// `pfq_class_term`, which for `z > 0` is a bound on its size.
///
/// With `sizes`, each term is the sum of the sizes of the terms of its series instead, times the
/// size of its factor. Their sum bounds the algebraic part at any phase of `z` and of its terms.
fn pfq_alg<A: Params, B: Params>(
    a: &A,
    b: &B,
    z: f64,
    sizes: bool,
    terms: &mut [Option<Scaled>],
) -> Option<()> {
    let (_, lower) = b.as_ref().split_last()?;
    let ln_z = ln_dd(Dd::new(z.abs()));
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
        let sum = pfq_alg_series(a, b, m, z, sizes)?;
        let mut y = gamma_series(num, den, sum, -ln_z * am)?;
        if sizes {
            y.m = y.m.abs();
        }
        terms[m] = Some(y);
    }
    Some(())
}

/// 2F2 by its asymptotic expansion (DLMF 16.11.7), if its estimated relative error is at
/// most `PFQ_TOL`
///
/// Only the exponentially large part is used for `z > 0`, and only the algebraic part for
/// `z < 0`. The size of the other part is added to the error.
fn hyp2f2_asymptotic(a: [f64; 2], b: [f64; 2], z: f64) -> Option<f64> {
    // the algebraic part
    let alg = |a: [f64; 2]| {
        let mut terms = [None; 2];
        if a[1] == b[1] {
            // they cancel
            pfq_alg(&[a[0]], &[b[0], 1.0], z, false, &mut terms)?;
        } else {
            pfq_alg(&a, &[b[0], b[1], 1.0], z, false, &mut terms)?;
        }
        Some(terms.into_iter().flatten())
    };
    let (y, other) = if z > 0.0 {
        (hyp2f2_exp(a, b, z)?, add_exp(alg(a)?))
    } else {
        (add_exp(alg(a)?)?, Some(hyp2f2_exp(a, b, z)?))
    };
    checked_f64(y, other)
}

/// `y` as `f64`, if its relative error plus the size of `other` relative to it is at most
/// `PFQ_TOL`, and if it does not underflow
fn checked_f64(y: Scaled, other: Option<Scaled>) -> Option<f64> {
    let y = fold_exp(y);
    let err = other.map(fold_exp).map_or(0.0, |v| {
        let size = v.m.abs() * (1.0 + v.err) / y.m.abs();
        ldexp_f64(size, v.e.saturating_sub(y.e))
    });
    // the result may overflow, but a subnormal one is not that accurate
    let accurate = (y.err + err <= PFQ_TOL).then(|| ldexp_f64(y.m, y.e));
    accurate.filter(|r| r.is_normal() || r.is_infinite())
}

/// Generalized hypergeometric function $_2F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Notes
///
/// For all $z$,
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
/// non-positive integer $a_i = -n$ makes it a polynomial of degree at most $n$.
///
/// As in [`hyp3f2`], an $a_i$ equal to a $b_j$ cancels it, also if it is a non-positive integer, so
/// that the result is the $_1F_1$ of the other parameters (or $e^z$), and it is $+\infty$ if a
/// $b_j$ is a non-positive integer $-m$, unless it is a polynomial of degree $\lt m$. It is NaN if
/// an argument is NaN, and overflows to $\pm\infty$ for large $z > 0$. For an infinite
/// argument, it is the limit as in [`hyppfq`].
///
/// # Returns `None`
///
/// If the relative error could exceed about $10^{-13}$:
///
/// - for some $z < 0$, e.g. $-50 < z < -25$ if all $\abs{a_i}, \abs{b_j} \le 3$, and wider for
///   larger parameters;
/// - for some large parameters, also for $z > 0$;
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
pub fn hyp2f2(a1: f64, a2: f64, b1: f64, b2: f64, z: f64) -> Option<f64> {
    if ![a1, a2, b1, b2, z].into_iter().all(f64::is_finite) {
        return hyppfq(&[a1, a2], &[b1, b2], z);
    } else if z == 0.0 {
        return Some(1.0);
    }
    let (mut a, mut b) = ([a1, a2], [b1, b2]);
    // a common parameter cancels, also at a pole; it becomes a_2 = b_2 = 1, which cancels like
    // the 1 of n! in the series, and which `hyp2f2_asymptotic` recognizes
    let pairs = [(0, 0), (0, 1), (1, 0), (1, 1)];
    if let Some((i, j)) = pairs.into_iter().find(|&(i, j)| a[i] == b[j]) {
        (a, b) = ([a[1 - i], 1.0], [b[1 - j], 1.0]);
        if a[0] == b[0] {
            // (not if it underflows, as in `hyp2f2_asymptotic`)
            return Some(z.exp()).filter(|r| r.is_normal() || r.is_infinite());
        }
    }
    let Some(degree) = pfq_degree(&a, &b) else {
        return Some(f64::INFINITY);
    };
    let series = || pfq_series_sorted(a, [b[0], b[1], 1.0], z);
    let asymptotic = || hyp2f2_asymptotic(a, b, z);
    if degree.is_finite() {
        series()
    } else if z.abs() < HYP2F2_ASYMPTOTIC_FIRST {
        series().or_else(asymptotic)
    } else {
        asymptotic().or_else(series)
    }
}

/// `v^(1 / k)` in double-double precision, for `v > 0`
fn root_dd(v: Dd, k: usize) -> Dd {
    let s = v.hi.powf(1.0 / k as f64);
    let pow = |n: usize| (0..n).fold(Dd::new(1.0), |p, _| p * s);
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
    let l = nu * 0.5 * ln_dd(Dd::new(z.abs())) - LN_PI_DD * 0.5;
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

/// 1F2 by its asymptotic expansion (DLMF 16.11.8), if its estimated relative error is at
/// most `PFQ_TOL`
///
/// For `z > 0`, only the exponentially large part is used, and the size of the algebraic part
/// is added to the error. The terms with `e^-x` are smaller than the truncation error.
fn hyp1f2_asymptotic(a: f64, b: [f64; 2], z: f64) -> Option<f64> {
    let exp = hyp1f2_exp(a, b, z)?;
    // the algebraic part, `None` if it vanishes
    let mut alg = [None];
    pfq_alg(&[a], &[b[0], b[1], 1.0], z, false, &mut alg)?;
    if z > 0.0 {
        checked_f64(exp, alg[0])
    } else {
        checked_f64(add_exp(alg.into_iter().flatten().chain([exp]))?, None)
    }
}

/// Generalized hypergeometric function $_1F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Notes
///
/// For all $z$,
///
/// $$
/// \hyp 1 2 a {b_1\enspace b_2} z
/// = \sum\_{n=0}^\infty {\rpow a n \over \rpow{b_1}{n} \ \rpow{b_2}{n}} {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1]. A
/// non-positive integer $a = -n$ makes it a polynomial of degree $n$.
///
/// As in [`hyp2f2`], an $a$ equal to a $b_j$ cancels it, also if it is a non-positive integer, so
/// that the result is the $_0F_1$ of the other $b_j$, and it is $+\infty$ if a $b_j$ is a
/// non-positive integer $-m$, unless it is a polynomial of degree $\lt m$. It is NaN if an argument
/// is NaN, and overflows to $\pm\infty$ for large $z > 0$. For an infinite argument, it is the
/// limit as in [`hyppfq`].
///
/// # Returns `None`
///
/// If the relative error could exceed about $10^{-13}$:
///
/// - for some $z < 0$: rarely if all $\abs a, \abs{b_j} \le 10$, but e.g. for about half of
///   $-10^5 < z < -10^3$ if they are up to 100;
/// - for most $z \lesssim -10^{32}$;
/// - near a zero, or if the result underflows;
/// - for some infinite arguments, as for [`hyppfq`].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function, $_0F_1(b; z)$
/// - [`hyp2f2`]: Generalized hypergeometric function,
///   $\hyp 2 2 {a_1\enspace a_2} {b_1\enspace b_2} z$
#[must_use]
pub fn hyp1f2(a: f64, b1: f64, b2: f64, z: f64) -> Option<f64> {
    if ![a, b1, b2, z].into_iter().all(f64::is_finite) {
        return hyppfq(&[a], &[b1, b2], z);
    } else if z == 0.0 {
        return Some(1.0);
    }
    let (mut a, mut b) = (a, [b1, b2]);
    // a common parameter cancels, also at a pole; it becomes a = b_2 = 1, which cancels like the
    // 1 of n! in the series
    if let Some(j) = b.iter().position(|&bj| bj == a) {
        (a, b) = (1.0, [b[1 - j], 1.0]);
    }
    let Some(degree) = pfq_degree(&[a], &b) else {
        return Some(f64::INFINITY);
    };
    let series = || pfq_series_sorted([a], [b[0], b[1], 1.0], z);
    let asymptotic = || hyp1f2_asymptotic(a, b, z);
    if degree.is_finite() {
        series()
    } else if z > -HYP1F2_ASYMPTOTIC_FIRST {
        series().or_else(asymptotic)
    } else {
        asymptotic().or_else(series)
    }
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
    let terms = successors(Some((1.0, Dd::new(1.0), phi)), |&(n, c, _)| {
        let next = -(c * phi * phi) / ((2.0 * n - 1.0) * (2.0 * n));
        Some((n + 1.0, next, next * phi / (2.0 * n + 1.0)))
    });
    let zero = (Dd::new(0.0), Dd::new(0.0));
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
    let root = root_dd(Dd::new(z.abs()), kappa) * if real { z.signum() } else { 1.0 };
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
        ((sum, err + others * sum.abs()), Dd::new(1.0))
    } else {
        pfq_exp_cos(&terms, settled, kappa, root, nu)
    };
    // (2 pi)^((1 - kappa) / 2) kappa^(-1/2) |z|^(nu / kappa) e^x
    let ln_power = nu / kf * ln_dd(Dd::new(z.abs())) - ln_dd(Dd::new(kf)) * 0.5;
    let l = ln_power - LN_SQRT_2PI_DD * (kf - 1.0) + root * kf * growth;
    let num = lower.iter().copied().map(Dd::new);
    let den = a.iter().copied().map(Dd::new);
    gamma_series(num, den, series, l)
}

/// pFq for `p <= q` by its asymptotic expansion (DLMF 16.11.7 to 16.11.9), if its
/// estimated relative error is at most `PFQ_TOL`, where `b` contains the 1 of `n!`
///
/// Only the exponentially large part is used for `z > 0`, and only the algebraic part for
/// `z < 0` and `p = q`. The size of the other part is added to the error, for `z > 0` as the
/// sizes of the terms of the algebraic part. For `z < 0` and `p < q`, it is the sum of both
/// parts.
fn pfq_asymptotic<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<f64> {
    let (a_ref, b_ref) = (a.as_ref(), b.as_ref());
    let kappa = b_ref.len() - a_ref.len();
    let mut terms = vec![None; a_ref.len()];
    pfq_alg(a, b, z, z > 0.0, &mut terms)?;
    let alg = terms.into_iter().flatten();
    let exp = pfq_exp(a_ref, b_ref, z)?;
    if z > 0.0 {
        checked_f64(exp, add_exp(alg))
    } else if kappa == 1 {
        checked_f64(add_exp(alg)?, Some(exp))
    } else {
        // (DLMF 16.11.9 for p <= q - 2 has the exponential part only; with the algebraic part,
        // it is theorem 3 of Paris, https://arxiv.org/abs/1708.04824)
        checked_f64(add_exp(alg.chain([exp]))?, None)
    }
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
        Some(f64::INFINITY * a.iter().chain(b).map(|&v| gammasgn(v)).product::<f64>())
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
/// # Notes
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
/// $p = q + 1$ and $z \le -1$, it is the analytic continuation. A non-positive integer
/// $a_i = -n$ makes it a polynomial of degree at most $n$.
///
/// Each $a_i$ equal to a $b_j$ cancels it, also if it is a non-positive integer or infinite. The
/// supported $p$ and $q$ that remain are:
///
/// | $p$ | $q = 0$                   | $q = 1$                   | $q = 2$    |
/// | --- | ------------------------- | ------------------------- | ---------- |
/// | 0   | [`hyp0f0`](crate::hyp0f0) | [`hyp0f1`](crate::hyp0f1) |            |
/// | 1   | [`hyp1f0`](crate::hyp1f0) | [`hyp1f1`](crate::hyp1f1) | [`hyp1f2`] |
/// | 2   |                           | [`hyp2f1`](crate::hyp2f1) | [`hyp2f2`] |
/// | 3   |                           |                           | [`hyp3f2`] |
///
/// and each other $p \le q + 1$.
///
/// It is:
///
/// - 1 at $z = 0$, and for a $b_j = +\infty$ and a finite $z$ if $p \le q + 1$ or it is a
///   polynomial, unless it is complex;
/// - $+\infty$ if a $b_j$ is a non-positive integer $-m$, unless it is a polynomial of degree
///   $\lt m$, for $_1F_0$ at $z = 1$ if $a_1 > 0$, and for $p = q + 1 \ge 3$ at $z = 1$ if
///   $\sum_j b_j \le \sum_i a_i$, unless it is a polynomial;
/// - infinite at $z = \pm\infty$ if it is a polynomial of degree $\ge 1$, with the sign of its
///   last term, and otherwise at $z = +\infty$ for $p \le q$, with the sign of
///   $\prod_j \Gamma(b_j) / \prod_i \Gamma(a_i)$;
/// - 0 at $z = -\infty$ for $q \le p \le q + 1$ if all $a_i > 0$;
/// - NaN for $p = q + 1$ and $z > 1$, where it is complex, unless it is a polynomial or
///   $_1F_0 = (1 - z)^{-a_1}$ with an integer $a_1$, and if an argument is NaN.
///
/// For $p = 2$, $q = 0$ and $z \lt 0$, where the series diverges, it is its Borel sum
/// $(-1/z)^{a_1} U(a_1, a_1 - a_2 + 1, -1/z)$ ([`hypu`](crate::hypu)) [^2], as in mpmath.
///
/// # Returns `None`
///
/// - for $p > q + 1$, unless it is a polynomial or $_2F_0$ with $z \lt 0$;
/// - for $_2F_0$ where all its methods fail their accuracy checks: $U$ in both parameter orders
///   (or where $(-1/z)^{a_1}$ or $U$ isn't normal), its series up to its smallest term, and its
///   Weniger transformation, e.g. for some $z \lt -1$ or $a_1, a_2 \ll 0$;
/// - where it can't be computed accurately: for $p = q + 1 \ge 3$ as for [`hyp3f2`], for the
///   other $p = q$ as for [`hyp2f2`], and for the other $p < q$ as for [`hyp1f2`];
/// - for the other infinite arguments, unless $p \le 2$ and $q \le 1$ where the function in the
///   table has a limit.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.6.E21>
#[must_use]
pub fn hyppfq(a: &[f64], b: &[f64], z: f64) -> Option<f64> {
    if a.iter().chain(b).chain([&z]).any(|v| v.is_nan()) {
        return Some(f64::NAN);
    } else if z == 0.0 {
        return Some(1.0);
    }
    // each b_j cancels an equal a_i, if any
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    let mut cancel = |bj: &f64| a.iter().position(|ai| ai == bj).map(|i| a.remove(i));
    b.retain(|bj| cancel(bj).is_none());
    let Some(degree) = pfq_degree(&a, &b) else {
        return Some(f64::INFINITY);
    };
    let complex = a.len() == b.len() + 1 && z > 1.0 && degree.is_infinite();
    if a.iter().chain(&b).any(|v| v.is_infinite()) {
        // the terms n >= 1 vanish for b_j = +inf
        let finite = a.iter().chain([&z]).all(|v| v.is_finite());
        let converges = a.len() <= b.len() + 1 || degree.is_finite();
        let one = finite && converges && !b.contains(&f64::NEG_INFINITY) && !complex;
        return one.then_some(1.0);
    } else if let Some(y) = pfq_limit(&a, &b, z, degree) {
        return Some(y);
    }
    // a polynomial is short enough for double-double
    let polynomial = || {
        let (y, err, _) = pfq_sum::<Dd>(&a, &[&b[..], &[1.0]].concat(), z, PFQ_DD_MAXITER);
        accurate((y, err))
    };
    if let Some(y) = degree.is_finite().then(polynomial).flatten() {
        return Some(y);
    }
    // the xsf functions are NaN where they are inaccurate
    let xsf = |y: f64| Some(y).filter(|y| !y.is_nan());
    match (a.as_slice(), b.as_slice()) {
        (&[], &[]) => Some(crate::hyp0f0(z)),
        // +inf where it diverges, and NaN where it is complex
        (&[a], &[]) if z == 1.0 && a > 0.0 => Some(f64::INFINITY),
        (&[a], &[]) if !complex || a.fract() == 0.0 => Some(crate::hyp1f0(a, z)),
        _ if complex => Some(f64::NAN),
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
            let borel = |a1: f64, a2: f64| {
                let (p, (u, digits)) = (x.powf(a1), xsf_hypu_digits(a1, a1 - a2 + 1.0, x));
                let enough = f64::from(digits) >= -PFQ_TOL.log10();
                (p.is_normal() && u.is_normal() && enough).then_some(p * u)
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
            let series = || accurate(asymptotic_sum(terms.clone(), settled, scale, |_| 0.0));
            // Try Weniger when both U evaluations and the truncated series miss the tolerance,
            // with the rounding errors of the ratios up to `settled`, which its estimate misses
            // (4 eps each, as in `pfq_weniger`, added linearly as `a + n` can round alike)
            let weniger = || {
                let (y, err) = pfq_weniger(&[a1, a2], &[1.0], z)?;
                let drift = terms
                    .clone()
                    .zip(0..=settled as u32)
                    .map(|(t, n)| f64::from(n) * t.abs());
                accurate((y, err + 4.0 * f64::EPSILON * drift.sum::<f64>()))
            };
            // then with a1 and a2 swapped (DLMF 13.2.40), as chgu may take another method, e.g.
            // chgubi instead of chguit for an integer a1 - a2
            borel(a1, a2)
                .or_else(series)
                .or_else(|| borel(a2, a1))
                .or_else(weniger)
        },
        (&[a], &[b1, b2]) => hyp1f2(a, b1, b2, z),
        (&[a1, a2], &[b1, b2]) => hyp2f2(a1, a2, b1, b2, z),
        (&[a1, a2, a3], &[b1, b2]) => hyp3f2(a1, a2, a3, b1, b2, z),
        _ if a.len() == b.len() + 1 => {
            b.push(1.0);
            pfq_real(a, b, z)
        },
        _ if a.len() <= b.len() => {
            b.push(1.0);
            let series = || pfq_series_sorted(a.clone(), b.clone(), z);
            let asymptotic = || pfq_asymptotic(&a, &b, z);
            // the terms of the series for z < 0 reach about e^x, and the result
            // e^(x cos(pi / kappa)), or a power of z for p = q
            let kappa = (b.len() - a.len()) as f64;
            let x = kappa * z.abs().powf(1.0 / kappa);
            let cancellation = x * (1.0 - cospi(1.0 / kappa).max(0.0));
            if degree.is_finite() {
                series()
            } else if z > 0.0 || cancellation < PFQ_ASYMPTOTIC_FIRST {
                series().or_else(asymptotic)
            } else {
                asymptotic().or_else(series)
            }
        },
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_root_dd() {
        // of a double-double with a low part
        let v = Dd::new(10.0) / 3.0;
        for k in 1..=5 {
            let root = root_dd(v, k);
            let pow = (0..k).fold(Dd::new(1.0), |p, _| p * root);
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
        // a pole
        assert!(mul_gammas([Dd::new(-2.0)], den, one).is_none());
        assert!(mul_gammas(num, [Dd::new(-2.0)], one).is_none());
    }

    #[test]
    fn test_pfq_asymptotic_integer_difference() {
        // the bound on the algebraic part is too large for equal a_i, and near them the
        // expansion is accurate or `None` (Arb reference)
        let c = 1.0 + pow2(-10);
        let (b, z) = ([c, c, 11.25, 1.0], 21.25);
        assert_eq!(pfq_asymptotic(&[1.0, 1.0, 10.5], &b, z), None);
        let near = pfq_asymptotic(&[1.0, 1.0 + pow2(-20), 10.5], &b, z);
        let expected = 739_413_030.005_544_4;
        assert!(near.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected));
    }

    #[test]
    fn test_pfq_asymptotic_nearly_cancelled() {
        // equal a_i within sqrt(eps) of b_2 and b_3, where the algebraic part is about 1 and
        // changes much with a_i: accurate or `None` (Arb reference)
        let t = -pow2(-40);
        let (root, root_2) = (f64::EPSILON.sqrt(), (2.0 * f64::EPSILON).sqrt());
        let y = pfq_asymptotic(&[t; 3], &[t - 1.0, t + root, t + root_2, 1.0], 30.0);
        let expected = -816_436.169_307_593_7_f64;
        assert!(y.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected.abs()));
    }
}
