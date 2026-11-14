//! Generalized hypergeometric functions ${}_pF_q$ of real arguments.
//!
//! The internal functions are generic over the numbers of parameters. They take the lower
//! parameters together with the 1 of $n!$ as `b`, so that the terms of the series are
//! $t_n = z^n \prod_i (a_i)_n / \prod_j (b_j)_n$.

use core::ops::Add;

use crate::utils::{is_nonpos_int, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{DoubleDouble, frexp_f64, ldexp_f64, two_sum};
use crate::xsf::cephes::{gamma, gammasgn, lgam, rgamma};

/// Maximum number of terms of `pfq_series`
const PFQ_MAXITER: u32 = 10_000;

/// Maximum number of terms of `pfq_series` in double-double arithmetic, fewer as a term takes
/// about 20 times as long
const PFQ_DD_MAXITER: u32 = PFQ_MAXITER / 4;

/// Largest estimated relative error of a result
const PFQ_TOL: f64 = 1e-13;

/// Largest size of terms that add nothing to a sum, even if there are `PFQ_MAXITER` of them
const PFQ_NEGLIGIBLE: f64 = 1e-150;

/// The sizes of the parameters for which `pfq_products` is exact enough for all $n$
const PFQ_MODERATE: core::ops::RangeInclusive<f64> = 1e-30..=1e30;

/// Below the logarithm of the smallest positive normal number
const PFQ_LN_UNDERFLOW: f64 = -709.0;

/// Largest distance of a parameter to a non-positive integer $-m$ at which its rounding error
/// is counted for the whole tail of the series from $n = m$ on, as errors of that size are seen
/// up to a distance of 0.01
const PFQ_NEAR_INTEGER: f64 = 0.0625;

/// Maximum order of `pfq_weniger`
const PFQ_WENIGER_MAXITER: usize = 40;

/// Largest `z` for `pfq_weniger`, above which the series needs less time (measured)
const PFQ_WENIGER_MAX_Z: f64 = -0.6;

/// Smallest `z` for which `pfq_weniger` is tried before `pfq_continuation`, as it is more often
/// accurate there (measured)
const PFQ_WENIGER_FIRST_Z: f64 = -1.5;

/// The arithmetic in which a series is summed, with the analogs of the `f64` constants
trait SeriesReal: Real + Add<f64, Output = Self> {
    /// Relative rounding error of an operation
    const EPSILON: f64;
    /// Smallest size of a result that is that accurate: it did not underflow
    const MIN_POSITIVE: f64;
}

impl SeriesReal for f64 {
    const EPSILON: f64 = f64::EPSILON;
    const MIN_POSITIVE: f64 = f64::MIN_POSITIVE;
}

impl SeriesReal for DoubleDouble {
    /// Not 2^-105, as its operations are not correctly rounded to 106 bits
    const EPSILON: f64 = pow2(-100);
    /// Below it, its low part is not normal
    const MIN_POSITIVE: f64 = pow2(106 - 1022);
}

/// The smallest $n$ with $v + n = 0$ for a `v` in `params`, or infinity if there is none
fn first_zero(params: &[f64]) -> f64 {
    let zeros = params.iter().filter(|&&v| is_nonpos_int(v));
    zeros.map(|v| -v).fold(f64::INFINITY, f64::min)
}

/// The $n$ from which on all $a_i + n$ and $b_j + n$ are positive
fn pfq_settled(a: &[f64], b: &[f64]) -> f64 {
    a.iter().chain(b).fold(0.0, |n, &v| n.max(-v))
}

/// `v` in descending order, as `pfq_ratio_bound` is tightest if `a` and `b` are sorted alike
fn descending<const N: usize>(mut v: [f64; N]) -> [f64; N] {
    v.sort_unstable_by(|x, y| y.total_cmp(x));
    v
}

/// Bound on $|t_{m+1} / t_m|$ for all $m \ge n$, if all $a_i + n$ and $b_j + n$ are positive
///
/// Each $a_i$ is paired with $b_i$: $(a_i + m) / (b_i + m)$ decreases to 1, or stays below 1.
/// The other $1 / (b_j + m)$ decrease. This is tightest if `a` and `b` are sorted alike.
fn pfq_ratio_bound<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
    n: f64,
) -> f64 {
    if P > Q {
        return f64::INFINITY;
    }
    let paired = a.iter().zip(b).map(|(a, b)| ((a + n) / (b + n)).max(1.0));
    let unpaired = b[P..].iter().map(|b| b + n);
    z.abs() * paired.product::<f64>() / unpaired.product::<f64>()
}

/// The products of the $a_i + n$ and of the $b_j + n$, whose quotient is $t_{n+1} / (z t_n)$
///
/// For up to four parameters each of a size in `PFQ_MODERATE`, an $a_i + n$ that is not zero is
/// at least $2^{-100}$ and at most $2^{101}$, so that both products and their quotient are normal.
fn pfq_products<T: SeriesReal, const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    n: f64,
) -> [T; 2] {
    let shifted = |p: T, &v: &f64| p * (T::from_f64(v) + n);
    let one = T::from_f64(1.0);
    [a.iter().fold(one, shifted), b.iter().fold(one, shifted)]
}

/// $t_{n+1} / (z t_n)$ for any parameters and $n$, which is zero if it underflows, and infinite
/// if it overflows or would be inexact, as a product in it is not normal
#[cold]
fn pfq_checked_ratio<T: SeriesReal, const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    n: f64,
) -> T {
    // normal in the arithmetic of `T`
    let is_normal = |v: T| (T::MIN_POSITIVE..=f64::MAX).contains(&v.to_f64().abs());
    let [num, den]: [T; 2] = pfq_products(a, b, n);
    let ln_abs = |params: &[f64]| params.iter().map(|v| (v + n).abs().ln()).sum::<f64>();
    if is_normal(num) && is_normal(den) && is_normal(num / den) {
        num / den
    } else if ln_abs(a) - ln_abs(b) < PFQ_LN_UNDERFLOW {
        T::from_f64(0.0)
    } else {
        // the products of the quotients of the paired parameters are more often normal
        let factor = |v: Option<&f64>| v.map_or(T::from_f64(1.0), |&v| T::from_f64(v) + n);
        let mut quotients = (0..P.max(Q)).map(|i| factor(a.get(i)) / factor(b.get(i)));
        let normal = |p: T, q: T| Some(p * q).filter(|&r| is_normal(q) && is_normal(r));
        let ratio = quotients.try_fold(T::from_f64(1.0), normal);
        ratio.unwrap_or(T::from_f64(f64::INFINITY))
    }
}

/// Bound on the size of the term after `term` if it underflows, or if its ratio does
#[cold]
#[inline(never)]
fn pfq_underflow_bound<T: SeriesReal>(term: T, z: f64) -> T {
    T::from_f64(T::MIN_POSITIVE * term.to_f64().abs().max(1.0) * z.abs().max(1.0))
}

/// Whether `pfq_products` is not exact enough for these parameters
fn pfq_extreme<const P: usize, const Q: usize>(a: &[f64; P], b: &[f64; Q]) -> bool {
    a.iter().chain(b).any(|v| !PFQ_MODERATE.contains(&v.abs()))
}

/// $t_{n+1} / (z t_n)$, which is `checked` if the parameters are extreme
#[inline]
fn pfq_ratio<T: SeriesReal, const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    n: f64,
    checked: bool,
) -> T {
    if checked {
        pfq_checked_ratio(a, b, n)
    } else {
        let [num, den]: [T; 2] = pfq_products(a, b, n);
        num / den
    }
}

/// The terms $t_1, t_2, \dots$ of the series after $t_0 = 1$, in the arithmetic of `T`
///
/// A term that underflows is a bound on it. The terms end before one that overflows or is NaN,
/// or that is such a bound and not negligible.
fn pfq_terms<T: SeriesReal, const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
) -> impl Iterator<Item = T> {
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

/// The series $\sum_n t_n$ (DLMF 16.2.1) and an estimate of its absolute error
///
/// The estimate is infinite if the series did not converge within `PFQ_MAXITER` terms, or if
/// its terms overflow, or underflow before they become large again. If it converged, but its
/// terms cancel, then it is summed again in double-double arithmetic.
pub(crate) fn pfq_series<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
) -> (f64, f64) {
    let (y, err, _) = pfq_series_tails(a, b, z);
    (y, err)
}

/// `pfq_series`, and the root of the sum of its squared tails $s - s_n$, which is at least the
/// largest of them
fn pfq_series_tails<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
) -> (f64, f64, f64) {
    let (y, err, tails) = pfq_sum::<f64, P, Q>(a, b, z, PFQ_MAXITER);
    if err <= PFQ_TOL * y.abs() || err == f64::INFINITY {
        (y, err, tails)
    } else {
        pfq_sum::<DoubleDouble, P, Q>(a, b, z, PFQ_DD_MAXITER)
    }
}

/// `pfq_series_tails` in the arithmetic of `T`, with at most `maxiter` terms
///
/// The rounding error of each ratio $t_{n+1} / t_n$ affects the terms after it, which sum to
/// $s - s_n$ for the partial sum $s_n$ and the sum $s$. These errors are added in quadrature,
/// like those of the partial sums. A part of them adds up instead, as `a + n` rounds alike
/// for all `n` between two powers of 2. The constants are measured.
fn pfq_sum<T: SeriesReal, const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
    maxiter: u32,
) -> (f64, f64, f64) {
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
    let mut terms = pfq_terms::<T, P, Q>(a, b, z);
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
            // geometric bound on the remaining terms
            let tail = t.abs() * rho / (1.0 - rho);
            if rho < 1.0 && tail <= f64::EPSILON * s.abs() {
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
/// ${}_pF_q$ in HypergeometricFunctions.jl \[2\]. For $z < 0$ it converges much faster than the
/// series, and also for $z < -1$, where the series diverges. With the term $t_{n+1}$ as the
/// remainder estimate of the partial sum $s_n$, its order $k$ is
///
/// $$
/// \delta_k = {\sum_{j=0}^k c_j s_j / t_{j+1} \over \sum_{j=0}^k c_j / t_{j+1}} , \quad
/// c_j = (-1)^j \binom k j (j + 1)_{k-1} ,
/// $$
///
/// by the three-term recursion of its numerator and denominator. It is applied from the $n$
/// where the terms start to alternate, so that the weights of the $s_j$ are positive. The same
/// recursion for $|s_j|$ then gives their weighted mean, which scales the rounding errors.
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
pub(crate) fn pfq_weniger<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
) -> Option<(f64, f64)> {
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
        return Some((sum, 4.0 * f64::EPSILON * scale));
    }
    let mut rows = [[0.0; 3]; PFQ_WENIGER_MAXITER];
    let (mut delta, mut step) = (sum, f64::INFINITY);
    for k in 0..PFQ_WENIGER_MAXITER {
        let n = start + k as f64;
        term *= ratio(n) * z;
        if !term.is_normal() {
            // overflow, or underflow, after which the terms only matter if they can increase
            let negligible = term.abs() < 1.0 && pfq_ratio_bound(a, b, z, n) <= 1.0;
            let noise = 4.0 * f64::EPSILON * scale;
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
        // the rounding errors of the partial sums and of the terms (the factor is measured)
        let noise = 4.0 * f64::EPSILON * ((abs / den).abs() + (num / den).abs() + scale);
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

/// $\Gamma(u) / \Gamma(w)$ and an estimate of its relative error, for `u` and `w` as `two_sum`s
///
/// The low part `e` of an argument `x + e` is its rounding error, which matters near a pole.
/// It is applied to first order: $\Gamma(x + e) \approx \Gamma(x) (1 + \psi(x) e)$, where the
/// square of $\psi(x) e$ estimates the error of that.
///
/// The ratio is zero only if `w` is a pole, and `None` if it underflows or overflows.
pub(crate) fn gamma_ratio(u: (f64, f64), w: (f64, f64)) -> Option<(f64, f64)> {
    let low = |(x, e): (f64, f64)| if e == 0.0 { 0.0 } else { crate::digamma(x) * e };
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
        let (low_u, low_w) = (low(u), low(w));
        let err = err + low_u * low_u + low_w * low_w;
        (ratio * (1.0 + low_u - low_w), err)
    })
}

/// The analytic continuation to $z < -1$ (DLMF 16.8.8) and an estimate of its absolute error, or
/// `None` if a factor of a coefficient overflows or underflows
///
/// This is for $p = q + 1$, where `b` contains the 1 of $n!$, and no two $a_i$ may differ by an
/// integer. It is a sum over $j$ of series in $1 / z$ with the parameters $a_j + 1 - b$ and
/// $a_j + 1 - a$.
pub(crate) fn pfq_continuation<const N: usize>(
    a: &[f64; N],
    b: &[f64; N],
    z: f64,
) -> Option<(f64, f64)> {
    let skip = |v: [f64; N], i| (0..N).filter(move |&k| k != i).map(move |k| v[k]);
    let one = b.iter().position(|&v| v == 1.0)?;
    let (mut sum, mut err) = (0.0, 0.0);
    for (j, &aj) in a.iter().enumerate() {
        if skip(*a, j).any(|ak| (ak - aj).fract() == 0.0) {
            return None;
        }
        // Gamma(a_k - a_j) / Gamma(a_k) and Gamma(b_k) / Gamma(b_k - a_j) for each k != j
        let ratios = skip(*a, j).map(|ak| gamma_ratio(two_sum(ak, -aj), (ak, 0.0)));
        let ratios = ratios.chain(skip(*b, one).map(|bk| gamma_ratio((bk, 0.0), two_sum(bk, -aj))));
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
        let (upper, lower) = (descending(b.map(&mut shift)), descending(a.map(&mut shift)));
        let (f, f_err, tails) = pfq_series_tails(&upper, &lower, 1.0 / z);
        sum += c * f;
        // those terms sum to a tail of the series; the last part is for a subnormal c
        err += c.abs() * (f_err + (c_err2.sqrt() + rounding) * f.abs() + near * tails)
            + f64::MIN_POSITIVE * f.abs();
    }
    Some((sum, err))
}

/// A method for $p = q + 1$, which returns the result and an estimate of its absolute error, or
/// `None`
type Method<const N: usize> = fn(&[f64; N], &[f64; N], f64) -> Option<(f64, f64)>;

/// ${}_{q+1}F_q$ for $z < 1$, or for a polynomial, by the first method for `z` that is accurate
///
/// `b` contains the 1 of $n!$.
pub(crate) fn pfq_accurate<const N: usize>(a: &[f64; N], b: &[f64; N], z: f64) -> Option<f64> {
    let series: Method<N> = |a, b, z| Some(pfq_series(a, b, z));
    let methods: &[Method<N>] = if first_zero(a).is_finite() || z > PFQ_WENIGER_MAX_Z {
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
    let mut results = methods.iter().filter_map(|method| method(a, b, z));
    let accurate = results.find(|&(y, err)| y.is_finite() && err <= PFQ_TOL * y.abs());
    accurate.map(|(y, _)| y)
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
/// non-positive integer, so that the result is that of `hyp2f1` for the other parameters, and
/// it is $+\infty$:
///
/// - if a $b_j$ is a non-positive integer $-m$, unless it is a polynomial of degree $\lt m$;
/// - at $z = 1$ if $b_1 + b_2 \le a_1 + a_2 + a_3$, unless it is a polynomial.
///
/// It is NaN for $z > 1$, where it is complex, unless it is a polynomial, and if an argument is
/// NaN.
///
/// # Returns `None`
///
/// If the relative error could exceed about $10^{-13}$:
///
/// - at $z = 1$, and for $0.996 \lesssim z < 1$;
/// - for $z < -1$ if two $a_i$ differ by an integer (nearly always for $z < -3$);
/// - if the terms cancel too much (e.g. large parameters), overflow or underflow;
/// - if an argument is infinite.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.8.E8>
///
/// # See also
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
#[must_use]
pub fn hyp3f2(a1: f64, a2: f64, a3: f64, b1: f64, b2: f64, z: f64) -> Option<f64> {
    let args = [a1, a2, a3, b1, b2, z];
    if args.into_iter().any(f64::is_nan) {
        return Some(f64::NAN);
    } else if args.into_iter().any(f64::is_infinite) {
        return None;
    }
    let (a, b) = ([a1, a2, a3], [b1, b2, 1.0]);
    let mut pairs = (0..3).flat_map(|i| (0..2).map(move |j| (i, j)));
    if let Some((i, j)) = pairs.find(|&(i, j)| a[i] == b[j]) {
        // a common parameter cancels; NaN marks an unreliable result of `hyp2f1`
        let y = crate::hyp2f1(a[(i + 1) % 3], a[(i + 2) % 3], b[1 - j], z);
        return Some(y).filter(|y| !y.is_nan());
    }
    if z == 0.0 {
        return Some(1.0);
    }
    // the degree of a polynomial, else infinity
    let degree = first_zero(&a);
    if first_zero(&b) < degree {
        return Some(f64::INFINITY);
    }
    if degree.is_infinite() {
        if z > 1.0 {
            // complex
            return Some(f64::NAN);
        } else if z == 1.0 && b1 + b2 <= a1 + a2 + a3 {
            // the series diverges
            return Some(f64::INFINITY);
        } else if z == 1.0 {
            return None;
        }
    }
    let (a, b) = (descending(a), descending(b));
    pfq_accurate(&a, &b, z)
}
