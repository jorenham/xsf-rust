//! Generalized hypergeometric functions ${}_pF_q$ of real arguments.
//!
//! The internal functions are generic over the numbers of parameters. They take the lower
//! parameters together with the 1 of $n!$ as `b`, so that the terms of the series are
//! $t_n = z^n \prod_i (a_i)_n / \prod_j (b_j)_n$.

use core::f64::consts::TAU;
use core::iter::{once, successors};
use core::ops::Add;

use crate::utils::{EPSILON_SQ, is_nonpos_int, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, frexp_f64, ldexp_f64, two_sum};
use crate::xsf::cephes::{gamma, gammasgn, lgam, rgamma};
use crate::xsf::{Scaled, add_exp, fold_exp, ln_dd, mul_gamma_ratio};

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

/// From this $\abs z$ on, `hyp2f2` tries the asymptotic expansion before the series, as for
/// $z < 0$ the terms of the series then reach about $e^{\abs z} > 10^{13}$, and the expansion
/// is often accurate and cheaper
const HYP2F2_ASYMPTOTIC_FIRST: f64 = 30.0;

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

impl SeriesReal for Dd {
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

/// The degree of the polynomial if an $a_i$ is a non-positive integer, else infinity, or `None`
/// at a pole, where a $(b_j)_n$ vanishes before
fn pfq_degree(a: &[f64], b: &[f64]) -> Option<f64> {
    let degree = first_zero(a);
    (first_zero(b) >= degree).then_some(degree)
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

/// Bound on $|t_{m+1} / t_m|$ for all $m \ge n$, if all $a_i + n$ and $b_j + n$ are positive,
/// except that for $P > Q$ the unpaired $a_i + m$ are taken at $m = n$ (and they grow)
///
/// Each $a_i$ is paired with $b_i$: $(a_i + m) / (b_i + m)$ decreases to 1, or stays below 1.
/// The other $1 / (b_j + m)$ decrease. This is tightest if `a` and `b` are sorted alike.
fn pfq_ratio_bound<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
    n: f64,
) -> f64 {
    let paired = a.iter().zip(b).map(|(a, b)| ((a + n) / (b + n)).max(1.0));
    let unpaired_a = a.iter().skip(Q).map(|a| a + n);
    let unpaired_b = b.iter().skip(P).map(|b| b + n);
    z.abs() * paired.product::<f64>() * unpaired_a.product::<f64>() / unpaired_b.product::<f64>()
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
        pfq_sum::<Dd, P, Q>(a, b, z, PFQ_DD_MAXITER)
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
            // geometric bound on the remaining terms, if the series converges (P <= Q)
            let tail = t.abs() * rho / (1.0 - rho);
            if P <= Q && rho < 1.0 && tail <= f64::EPSILON * s.abs() {
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
            let negligible = term.abs() < 1.0 && P <= Q && pfq_ratio_bound(a, b, z, n) <= 1.0;
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
    methods.iter().find_map(|method| accurate(method(a, b, z)?))
}

/// The result `y`, if it is finite and its estimated absolute error is at most `PFQ_TOL |y|`
fn accurate((y, err): (f64, f64)) -> Option<f64> {
    (y.is_finite() && err <= PFQ_TOL * y.abs()).then_some(y)
}

/// `pfq_series` if it is `accurate`, after sorting the parameters, which tightens
/// `pfq_ratio_bound`
fn pfq_series_sorted<const P: usize, const Q: usize>(
    a: [f64; P],
    b: [f64; Q],
    z: f64,
) -> Option<f64> {
    accurate(pfq_series(&descending(a), &descending(b), z))
}

/// Whether the series of ${}_{q+1}F_q$ diverges at $z = 1$, where $\sum_j b_j \le \sum_i a_i$
/// for `b` without the 1 of $n!$, or `None` if that is uncertain
///
/// The sums are in double-double, as their rounding can change the sign of their difference.
fn diverges_at_one(a: &[f64], b: &[f64]) -> Option<bool> {
    let sum = |v: &[f64]| v.iter().fold(Dd::new(0.0), |s, &v| s + v);
    let excess = (sum(b) - sum(a)).hi;
    let size = a.iter().chain(b).map(|v| v.abs()).sum::<f64>();
    // certain if no partial sum of the excess rounds
    let mut terms = a.iter().map(|v| -v).chain(b.iter().copied());
    let exact = terms.try_fold(0.0, |s, v| {
        let (s, e) = two_sum(s, v);
        (e == 0.0).then_some(s)
    });
    let uncertain = exact.is_none() && excess.abs() <= 8.0 * EPSILON_SQ * size;
    (!uncertain).then_some(excess <= 0.0)
}

/// ${}_{q+1}F_q$ for finite arguments, where the last of `b` is the 1 of $n!$
fn pfq_real<const N: usize>(a: [f64; N], b: [f64; N], z: f64) -> Option<f64> {
    if z == 0.0 {
        return Some(1.0);
    }
    let Some(degree) = pfq_degree(&a, &b) else {
        return Some(f64::INFINITY);
    };
    if degree.is_infinite() && z > 1.0 {
        // complex
        Some(f64::NAN)
    } else if degree.is_infinite() && z == 1.0 {
        // +inf where the series diverges, and no method where it converges
        diverges_at_one(&a, &b[..N - 1])?.then_some(f64::INFINITY)
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
/// non-positive integer, so that the result is the $_2F_1$ of the other parameters, and
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
    // a common parameter cancels: 2F1, or 1F0 if that has one as well
    let mut pairs = (0..3).flat_map(|i| (0..2).map(move |j| (i, j)));
    let Some((i, j)) = pairs.find(|&(i, j)| a[i] == b[j]) else {
        return pfq_real(a, b, z);
    };
    let (a, b) = ([a[(i + 1) % 3], a[(i + 2) % 3]], [b[1 - j], 1.0]);
    match a.iter().position(|&v| v == b[0]) {
        Some(i) => pfq_real([a[1 - i]], [1.0], z),
        None => pfq_real(a, b, z),
    }
}

/// The sum of an asymptotic series in $1/z$ up to its smallest term from term `settled` on, and
/// an estimate of its absolute error
///
/// From term `settled` on, all parameters plus $n$ must be positive. Before it, the terms may
/// vanish or grow. From term `settled + |z|` on, the term ratios are about $n / \abs z$, so that
/// growing terms keep growing. The size of a term $t_n$ is $\max(\abs{t_n}, \abs{t_{n+1}})$, as
/// a $t_n$ of a recurrence of order 2 can vanish while $t_{n+1}$ does not. It stops at a term of
/// a negligible size if `ratio_bound(n) < 1`, so that the next terms are smaller.
///
/// The error is the smallest size times $\sqrt{2 \pi \abs z}$, about the sum of the terms near
/// it, plus rounding. For terms of one sign, the remainder can be half of that, as for
/// $\sum_n n! / z^{n+1}$ with $z > 0$, whose term ratios are also about $n / z$.
fn asymptotic_sum(
    terms: impl Iterator<Item = f64>,
    settled: f64,
    z: f64,
    ratio_bound: impl Fn(f64) -> f64,
) -> (f64, f64) {
    let (mut sum, mut abs_sum, mut prev) = (0.0_f64, 0.0, f64::INFINITY);
    // the smallest size, the sum before its term, and its error
    let (mut min, mut best) = (f64::INFINITY, (sum, f64::INFINITY));
    let remainder = (TAU * z.abs()).sqrt().max(1.0);
    let mut terms = terms.peekable();
    for n in (0..PFQ_MAXITER).map(f64::from) {
        let Some(t) = terms.next() else { break };
        let size = t.abs().max(terms.peek().map_or(0.0, |next| next.abs()));
        if !t.is_finite() || n >= settled + z.abs() && t.abs() >= prev {
            break;
        } else if n >= settled && size < min {
            min = size;
            best = (sum, remainder * min + n * f64::EPSILON * abs_sum);
            if min <= f64::EPSILON * sum.abs() && ratio_bound(n) < 1.0 {
                break;
            }
        }
        sum += t;
        abs_sum += t.abs();
        prev = t.abs();
    }
    best
}

/// The terms $t_0 = 1, t_1, \dots$ of a recurrence of order 2, where `next(k, t1, t2)` is $t_k$
/// for $t_{k-1}$ and $t_{k-2}$, with $t_{-1} = 0$
fn recurrence_terms(next: impl Fn(f64, f64, f64) -> f64) -> impl Iterator<Item = f64> {
    // k, t_k, and t_(k-1)
    let terms = successors(Some((0.0, 1.0, 0.0)), move |&(k, t1, t2)| {
        Some((k + 1.0, next(k + 1.0, t1, t2), t1))
    });
    terms.map(|(_, t, _)| t)
}

/// `exp(l) sum` times `Gamma(num_0) ... / (Gamma(den_0) ...)`, where `err` is the absolute
/// error of `sum`
fn gamma_series(num: &[Dd], den: &[Dd], (sum, err): (f64, f64), l: Dd) -> Option<Scaled> {
    let y = Scaled::new(sum, 0, err / sum.abs())?.with_l(l);
    mul_gamma_ratio(num, den, || Some(y))
}

/// The exponentially large part of the asymptotic expansion of $_2F_2$ (DLMF 16.11.1, 16.11.7),
/// with $\abs z^\nu$ instead of $z^\nu$:
///
/// $$
/// {\Gamma(b_1) \Gamma(b_2) \over \Gamma(a_1) \Gamma(a_2)} e^z \abs z^\nu \sum_k c_k z^{-k},
/// \qquad \nu = a_1 + a_2 - b_1 - b_2
/// $$
///
/// By the differential equation (DLMF 16.8.3), $c_0 = 1$ and
/// $k c_k = A_1 c_{k-1} + A_0 c_{k-2}$, with $s = \nu + 1 - k$ and
///
/// - $A_1 = 2 s^2 + (2 b_1 + 2 b_2 - a_1 - a_2 - 1) s + b_1 b_2 - a_1 a_2$,
/// - $A_0 = (s + 1) (s + b_1) (s + b_2)$.
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
    gamma_series(&b.map(Dd::new), &a.map(Dd::new), series, l)
}

/// The algebraic part of the asymptotic expansion of $_2F_2$ (DLMF 16.11.2, 16.11.7), with
/// $\abs z$ instead of $-z$, as its 2 terms, or `None` if $a_1 - a_2$ is an integer (unless
/// $a_2 = b_2$)
///
/// The term of $a_m$ is $\abs z^{-a_m}$ times a $_3F_1$ series in $-1/z$, with the parameters
/// $a_m$, $1 + a_m - b_1$, $1 + a_m - b_2$, and $1 + a_m - a_n$. It vanishes (as `None`) if
/// $1 / \Gamma(b_j - a_m)$ does. If $a_2 = b_2$, they cancel, and only the term of $a_1$ is left.
fn hyp2f2_alg(a: [f64; 2], b: [f64; 2], z: f64) -> Option<[Option<Scaled>; 2]> {
    let is_pole = |v: &Dd| v.hi <= 0.0 && v.is_integer();
    let cancels = a[1] == b[1];
    if !cancels && (Dd::new(a[1]) - a[0]).is_integer() {
        return None;
    }
    let (b, ln_z) = (b.map(Dd::new), ln_dd(Dd::new(z.abs())));
    let mut terms = [None; 2];
    let count = if cancels { 1 } else { 2 };
    for m in 0..count {
        let (am, an) = (a[m], a[1 - m]);
        let (bm, an_am) = (b.map(|bj| bj - am), Dd::new(an) - am);
        let (num, den): (&[Dd], &[Dd]) = if cancels {
            // Gamma(b_2) / Gamma(a_2) and Gamma(a_2 - a_1) / Gamma(b_2 - a_1) cancel
            (&[b[0]], &[bm[0]])
        } else {
            (&[b[0], b[1], an_am], &[Dd::new(an), bm[0], bm[1]])
        };
        if den.iter().any(is_pole) {
            continue;
        }
        let up = descending([am, 1.0 - bm[0].hi, 1.0 - bm[1].hi]);
        let low = descending([1.0 + am - an, 1.0]);
        let bound = |n| pfq_ratio_bound(&up, &low, 1.0 / z, n);
        let series = once(1.0).chain(pfq_terms(&up, &low, -1.0 / z));
        let sum = asymptotic_sum(series, pfq_settled(&up, &low), z, bound);
        terms[m] = Some(gamma_series(num, den, sum, -ln_z * am)?);
    }
    Some(terms)
}

/// $_2F_2$ by its asymptotic expansion (DLMF 16.11.7), if its estimated relative error is at
/// most `PFQ_TOL`
///
/// Only the exponentially large part is used for $z > 0$, and only the algebraic part for
/// $z < 0$. The size of the other part is added to the error.
fn hyp2f2_asymptotic(a: [f64; 2], b: [f64; 2], z: f64) -> Option<f64> {
    // the algebraic part, `Some(None)` if it vanishes
    let alg = |a| Some(add_exp(hyp2f2_alg(a, b, z)?.into_iter().flatten()));
    let (y, other) = if z > 0.0 {
        // (its size is continuous in a_1, also at an integer a_1 - a_2)
        let nudged = [a[0] + f64::EPSILON.sqrt() * a[0].abs().max(1.0), a[1]];
        (hyp2f2_exp(a, b, z)?, alg(a).or_else(|| alg(nudged))?)
    } else {
        (alg(a)??, Some(hyp2f2_exp(a, b, z)?))
    };
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
/// an argument is NaN, and overflows to $\pm\infty$ for large $z > 0$.
///
/// # Returns `None`
///
/// If the relative error could exceed about $10^{-13}$:
///
/// - for some $z < 0$, e.g. $-50 < z < -25$ if all $\abs{a_i}, \abs{b_j} \le 3$, and wider for
///   larger parameters;
/// - for $z < 0$ beyond that if $a_1 - a_2$ is an integer (unless one of them cancels);
/// - for some large parameters, also for $z > 0$;
/// - near a zero, or if the result underflows;
/// - if an argument is infinite.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp3f2`]: Generalized hypergeometric function,
///   $\hyp 3 2 {a_1\enspace a_2\enspace a_3} {b_1\enspace b_2} z$
#[must_use]
pub fn hyp2f2(a1: f64, a2: f64, b1: f64, b2: f64, z: f64) -> Option<f64> {
    let args = [a1, a2, b1, b2, z];
    if args.into_iter().any(f64::is_nan) {
        return Some(f64::NAN);
    } else if args.into_iter().any(f64::is_infinite) {
        return None;
    } else if z == 0.0 {
        return Some(1.0);
    }
    let (mut a, mut b) = ([a1, a2], [b1, b2]);
    // a common parameter cancels, also at a pole; it becomes a_2 = b_2 = 1, which cancels like
    // the 1 of n! in the series, and which `hyp2f2_alg` recognizes
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
/// is NaN.
///
/// # Returns `None`
///
/// If the relative error could exceed about $10^{-13}$:
///
/// - for some $z < 0$, e.g. $z \lesssim -400$ if all $\abs a, \abs{b_j} \le 3$;
/// - near a zero;
/// - if the result overflows;
/// - if an argument is infinite.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function, $_0F_1(b; z)$
/// - [`hyp2f2`]: Generalized hypergeometric function,
///   $\hyp 2 2 {a_1\enspace a_2} {b_1\enspace b_2} z$
#[must_use]
pub fn hyp1f2(a: f64, b1: f64, b2: f64, z: f64) -> Option<f64> {
    let args = [a, b1, b2, z];
    if args.into_iter().any(f64::is_nan) {
        return Some(f64::NAN);
    } else if args.into_iter().any(f64::is_infinite) {
        return None;
    } else if z == 0.0 {
        return Some(1.0);
    }
    let (mut a, mut b) = (a, [b1, b2]);
    // a common parameter cancels, also at a pole; it becomes a = b_2 = 1, which cancels like the
    // 1 of n! in the series
    if let Some(j) = b.iter().position(|&bj| bj == a) {
        (a, b) = (1.0, [b[1 - j], 1.0]);
    }
    if pfq_degree(&[a], &b).is_none() {
        return Some(f64::INFINITY);
    }
    pfq_series_sorted([a], [b[0], b[1], 1.0], z)
}
