//! Generalized hypergeometric functions ${}_pF_q$ of real arguments.
//!
//! The internal functions are generic over the numbers of parameters. They take the lower
//! parameters together with the 1 of $n!$ as `b`, so that the terms of the series are
//! $t_n = z^n \prod_i (a_i)_n / \prod_j (b_j)_n$.

use crate::utils::is_nonpos_int;

/// Maximum number of terms of `pfq_series`
const PFQ_MAXITER: u32 = 10_000;

/// Largest estimated relative error of a result
const PFQ_TOL: f64 = 1e-13;

/// Largest size of terms that add nothing to a sum, even if there are `PFQ_MAXITER` of them
const PFQ_NEGLIGIBLE: f64 = 1e-150;

/// The sizes of the parameters for which `pfq_products` is exact enough for all $n$
const PFQ_MODERATE: core::ops::RangeInclusive<f64> = 1e-30..=1e30;

/// Below the logarithm of the smallest positive normal number
const PFQ_LN_UNDERFLOW: f64 = -709.0;

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
/// For up to five parameters each of a size in `PFQ_MODERATE`, an $a_i + n$ that is not zero is
/// at least $2^{-100}$ and at most $2^{101}$, so that both products and their quotient are normal.
fn pfq_products<const P: usize, const Q: usize>(a: &[f64; P], b: &[f64; Q], n: f64) -> [f64; 2] {
    let shifted = |p: f64, v: &f64| p * (v + n);
    [a.iter().fold(1.0, shifted), b.iter().fold(1.0, shifted)]
}

/// $t_{n+1} / (z t_n)$ for any parameters and $n$, which is zero if it underflows, and infinite
/// if it overflows or would be inexact, as a product in it is not normal
#[cold]
fn pfq_checked_ratio<const P: usize, const Q: usize>(a: &[f64; P], b: &[f64; Q], n: f64) -> f64 {
    let [num, den] = pfq_products(a, b, n);
    let ln_abs = |params: &[f64]| params.iter().map(|v| (v + n).abs().ln()).sum::<f64>();
    if num.is_normal() && den.is_normal() && (num / den).is_normal() {
        num / den
    } else if ln_abs(a) - ln_abs(b) < PFQ_LN_UNDERFLOW {
        0.0
    } else {
        // the products of the quotients of the paired parameters are more often normal
        let factor = |v: Option<&f64>| v.map_or(1.0, |v| v + n);
        let mut quotients = (0..P.max(Q)).map(|i| factor(a.get(i)) / factor(b.get(i)));
        let normal = |p: f64, q: f64| Some(p * q).filter(|r| q.is_normal() && r.is_normal());
        quotients.try_fold(1.0, normal).unwrap_or(f64::INFINITY)
    }
}

/// Bound on the size of the term after `term` if it underflows, or if its ratio does
#[cold]
#[inline(never)]
fn pfq_underflow_bound(term: f64, z: f64) -> f64 {
    f64::MIN_POSITIVE * term.abs().max(1.0) * z.abs().max(1.0)
}

/// The terms $t_1, t_2, \dots$ of the series after $t_0 = 1$
///
/// A term that underflows is a bound on it. The terms end before one that overflows or is NaN,
/// or that is such a bound and not negligible.
fn pfq_terms<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
) -> impl Iterator<Item = f64> {
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let checked = a.iter().chain(b).any(|v| !PFQ_MODERATE.contains(&v.abs()));
    // the largest term, which is small once the terms are only upper bounds
    let mut limit = f64::MAX;
    let indices = (0_u32..).map(f64::from);
    indices.scan(1.0_f64, move |term, n| {
        let ratio = if checked {
            pfq_checked_ratio(a, b, n)
        } else {
            let [num, den] = pfq_products(a, b, n);
            num / den
        };
        let next = *term * (ratio * z);
        // (the terms after the last one of a polynomial are zero; tested last, as that is faster)
        *term = if next.abs() < f64::MIN_POSITIVE && n < degree {
            // the term or its ratio underflows, after which the terms can become large again:
            // follow a bound on them
            limit = PFQ_NEGLIGIBLE;
            pfq_underflow_bound(*term, z)
        } else {
            next
        };
        (term.abs() <= limit).then_some(*term)
    })
}

/// The series $\sum_n t_n$ (DLMF 16.2.1) and an estimate of its absolute error
///
/// The estimate is infinite if the series did not converge within `PFQ_MAXITER` terms, or if
/// its terms overflow, or underflow before they become large again.
///
/// The rounding error of each ratio $t_{n+1} / t_n$ affects the terms after it, which sum to
/// $s - s_n$ for the partial sum $s_n$ and the sum $s$. These errors are added in quadrature,
/// like those of the partial sums. A part of them adds up instead, as `a + n` rounds alike
/// for all `n` between two powers of 2. The constants are measured.
pub(crate) fn pfq_series<const P: usize, const Q: usize>(
    a: &[f64; P],
    b: &[f64; Q],
    z: f64,
) -> (f64, f64) {
    let settled = pfq_settled(a, b);
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let mut sum = 1.0_f64;
    // the sums of s_n and s_n^2
    let (mut sum1, mut sum2) = (0.0, 0.0_f64);
    let rounding = |n: f64, sum: f64, sum1: f64, sum2: f64| {
        // the root of the sum of (s - s_n)^2, which is accurate enough for an estimate, relative
        // to s so that its parts don't overflow
        let tails = if sum == 0.0 {
            sum2.sqrt()
        } else {
            ((n + 1.0) - 2.0 * (sum1 / sum) + sum2 / sum / sum)
                .max(0.0)
                .sqrt()
                * sum.abs()
        };
        f64::EPSILON * ((4.0 + 0.25 * (n + 1.0).sqrt()) * tails + sum2.sqrt())
    };
    let mut terms = pfq_terms(a, b, z);
    for n in (0..PFQ_MAXITER).map(f64::from) {
        sum1 += sum;
        sum2 += sum * sum;
        if n == degree {
            // the series terminated
            return (sum, rounding(n, sum, sum1, sum2));
        }
        // t_(n+1)
        let Some(term) = terms.next() else { break };
        sum += term;
        if n >= settled && term.abs() <= f64::EPSILON * sum.abs() {
            let rho = pfq_ratio_bound(a, b, z, n + 1.0);
            // geometric bound on the remaining terms
            let tail = term.abs() * rho / (1.0 - rho);
            if rho < 1.0 && tail <= f64::EPSILON * sum.abs() {
                return (sum, rounding(n, sum, sum1, sum2) + tail);
            }
        }
    }
    (sum, f64::INFINITY)
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
/// $z < -1$, it is the analytic continuation [^1]. A non-positive integer $a_i = -n$ makes it
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
/// - for $\abs z \ge 1$, unless it is a polynomial;
/// - for $\abs z \gtrsim 0.996$;
/// - if the terms cancel too much (e.g. large parameters, near a zero), overflow or underflow;
/// - if an argument is infinite.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
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
        } else if z.abs() >= 1.0 {
            return None;
        }
    }
    let (a, b) = (descending(a), descending(b));
    let (y, err) = pfq_series(&a, &b, z);
    (y.is_finite() && err <= PFQ_TOL * y.abs()).then_some(y)
}
