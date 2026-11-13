//! Translated from `scipy/special/_orthogonal.py` (scipy/scipy@7d279ba)
//!
//! - <https://github.com/scipy/scipy/blob/7d279ba/scipy/special/_orthogonal.py>
//!
//! The nodes are in ascending order, and the weights sum to the integral of the weight function.

use core::f64::consts::{LN_2, PI};

use super::orthogonal_eval::{eval_gegenbauer, eval_jacobi, eval_legendre};
use super::tridiag::eigvals;
use crate::utils::MAXLOG;
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, LOG2, ldexp_f64, two_sum};
use crate::xsf::cephes::incbet_lnpow;
use Nodes::{Plain, Shifted, Symmetric};

/// Corresponds to `_gen_roots_and_weights`: Gauss quadrature by the Golub-Welsch algorithm
///
/// The monic orthogonal polynomials satisfy `p_{k+1}(x) = (x - a_k) p_k(x) - b_k^2 p_{k-1}(x)`,
/// so the roots of `f(n, x)` are the eigenvalues of the tridiagonal matrix with diagonal `a_k` and
/// off-diagonal `b_k`. They are refined by a Newton step with the derivative `df(n, x)`.
///
/// The weights sum to `mu0`, and are proportional to `1 / (sigma(x) df(n, x)^2)` at the roots,
/// where `sigma f'' + tau f' + lambda f = 0` (DLMF Table 18.8.1). For a root `x + dx` with the
/// rounding error `dx`, `sigma(x, dx)` is `sigma(x + dx)`, and `df(n, x + dx)` is
/// `df(n, x) (1 - tau(x) dx / sigma)` to first order, since `f(n, x + dx) = 0`.
fn gen_roots_and_weights(
    n: usize,
    mu0: f64,
    an: impl Fn(f64) -> f64,
    bn: impl Fn(f64) -> f64,
    f: impl Fn(i32, f64) -> f64,
    df: impl Fn(i32, f64) -> f64,
    sigma: impl Fn(f64, f64) -> f64,
    tau: impl Fn(f64) -> f64,
    nodes: Nodes,
) -> (Vec<f64>, Vec<f64>) {
    let d: Vec<f64> = (0..n).map(|k| an(k as f64)).collect();
    let e: Vec<f64> = (1..n).map(|k| bn(k as f64)).collect();
    let mut x = eigvals(&d, &e);

    // improve roots by Newton's method
    //
    // Unlike in SciPy, which takes one step, at most 4 steps are taken while
    // `step^2 > eps sigma^2`, i.e. while the error of a step, about `step^2 / sigma`, matters for
    // `sigma` and thus for the weights. This only happens for a root very close to an endpoint (a
    // zero of `sigma`), e.g. within 1e-14 of -1 for alpha and beta near -1.
    //
    // Unlike in SciPy, a root and its weight where the polynomials overflow (e.g. `L_n` for
    // n >= 400) are the eigenvalue and 0, instead of NaN for all.
    let n = i32::try_from(n).expect("n is too large");
    // the rounding errors of the refined roots
    let mut dx = vec![0.0; x.len()];
    for (xi, dxi) in x.iter_mut().zip(&mut dx) {
        for _ in 0..4 {
            let step = f(n, *xi) / df(n, *xi);
            if !step.is_finite() {
                break;
            }
            (*xi, *dxi) = two_sum(*xi, -step);
            let s = sigma(*xi, *dxi);
            if step * step <= f64::EPSILON * s * s {
                break;
            }
        }
    }

    // Unlike in SciPy, the weights are computed from `df` at the refined roots, instead of
    // `1 / (f(n - 1, x) df(n, x))` with `df` at the eigenvalues, which is more sensitive to the
    // rounding of `x` (about n^2 times near x = +-1).
    //
    // They are `mu0 / sum(1 / r^2) / r^2`, with `r = sqrt(sigma) df` scaled to `r >= 1`, so that
    // only the weights that are tiny relative to `mu0 / n` underflow. A non-finite `r` is an
    // overflow, with weight 0.
    let r: Vec<f64> = x
        .iter()
        .zip(&dx)
        .map(|(&xi, &dxi)| {
            let s = sigma(xi, dxi);
            s.sqrt() * df(n, xi) * (1.0 - tau(xi) * dxi / s)
        })
        .collect();
    // A zero `r` is an underflow, with weight NaN, and the others 0 (relatively tiny), instead of
    // NaN for all.
    let r_min = r
        .iter()
        .map(|ri| ri.abs())
        .filter(|&ri| ri > 0.0)
        .fold(f64::INFINITY, f64::min);
    let r: Vec<f64> = r
        .iter()
        .map(|&ri| {
            if ri.is_finite() {
                ri / r_min
            } else {
                f64::INFINITY
            }
        })
        .collect();
    let scale = mu0 / r.iter().map(|ri| 1.0 / (ri * ri)).sum::<f64>();
    let mut w: Vec<f64> = r.iter().map(|ri| scale / ri / ri).collect();

    match nodes {
        Plain => {},
        Symmetric => {
            w = w
                .iter()
                .zip(w.iter().rev())
                .map(|(a, b)| (a + b) / 2.0)
                .collect();
            x = x
                .iter()
                .zip(x.iter().rev())
                .map(|(a, b)| (a - b) / 2.0)
                .collect();
        },
        Shifted => {
            x = x
                .iter()
                .zip(&dx)
                .map(|(xi, dxi)| ((xi + 1.0) + dxi) / 2.0)
                .collect();
        },
    }
    (x, w)
}

/// What [`gen_roots_and_weights`] does with the roots
#[derive(Clone, Copy)]
enum Nodes {
    Plain,
    /// Makes the roots antisymmetric and the weights symmetric
    Symmetric,
    /// Shifts the roots `x` to `[0, 1]` as `(x + 1) / 2`, but with their rounding errors, so that
    /// the ones near 0 have a small relative error
    Shifted,
}

/// `(m, e)` with `m 2^e = exp(l)`, and `e >= 0` only as large as needed for a finite `m`
fn exp_scaled(l: Dd) -> (f64, f64) {
    let e = ((l.hi - MAXLOG) / LN_2).ceil().max(0.0);
    let r = l - LOG2 * e;
    (r.hi.exp() * (1.0 + r.lo), e)
}

/// The rule with the weights multiplied by `2^e`, for weights computed with `mu0 / 2^e`
fn ldexp_weights((x, w): (Vec<f64>, Vec<f64>), e: f64) -> (Vec<f64>, Vec<f64>) {
    (x, w.into_iter().map(|wi| ldexp_f64(wi, e as i32)).collect())
}

/// `1 - (x + dx)^2` for `|dx| <= ulp(x)`, with a small relative error near `x + dx = +-1`
fn one_minus_x2(x: f64, dx: f64) -> f64 {
    ((1.0 - x) - dx) * ((1.0 + x) + dx)
}

/// Panics if `n` is zero, like the `ValueError` in SciPy
fn check_n(n: usize) {
    assert!(n >= 1, "n must be a positive integer");
}

/// Gauss-Jacobi quadrature
///
/// Corresponds to [`scipy.special.roots_jacobi`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_jacobi.html
///
/// The nodes are the roots of the Jacobi polynomial $P_n^{(\alpha,\beta)}$, and the rule is exact
/// for polynomials of degree $2n - 1$ with the weight function $(1 - x)^\alpha (1 + x)^\beta$ on
/// $[-1, 1]$ (A&S 22.2.1).
///
/// # Panics
/// - if `n` is zero
/// - if $\alpha \le -1$ or $\beta \le -1$
///
/// # See also
/// - [`eval_jacobi`]: Evaluate Jacobi polynomials, $P_n^{(\alpha,\beta)}$
#[must_use]
pub fn roots_jacobi(n: usize, alpha: f64, beta: f64) -> (Vec<f64>, Vec<f64>) {
    check_jacobi(n, alpha, beta);
    if alpha == 0.0 && beta == 0.0 {
        return roots_legendre(n);
    }
    // Unlike in SciPy, the rule for alpha == beta < -1/2 isn't the Gauss-Gegenbauer one, which is
    // less accurate near -1 (see `roots_gegenbauer`).
    if alpha == beta && alpha >= -0.5 {
        return roots_gegenbauer(n, alpha + 0.5);
    }

    let (a, b) = (alpha, beta);
    // Unlike in SciPy, the weights sum to `m` and are then scaled by `2^e`, with
    // `m 2^e = mu0 = 2^s B(p, q)` for `s = a + b + 1`, `p = a + 1` and `q = b + 1`, so that they
    // can be finite where `mu0` overflows. Where `B(p, q)` underflows, `ln(mu0)` is computed by
    // `incbet_lnpow` (for `p, q >= 13`), instead of `s ln(2) + betaln(p, q)`, which cancels.
    let (s, p, q) = (a + b + 1.0, a + 1.0, b + 1.0);
    let beta_pq = crate::beta(p, q);
    let (m, e) = if beta_pq.is_normal() || p.min(q) < 13.0 {
        // `2^(s - e) <= 2^1000`, and also `2^(s - e) B(p, q) <= 2^1000` if `B(p, q) > 1`
        let e = (s - 1000.0 + beta_pq.log2().max(0.0)).ceil().max(0.0);
        ((s - e).exp2() * beta_pq, e)
    } else {
        // `ln(2^-(p + q) / B(p, q))` is symmetric, but `incbet_lnpow` needs `p >= q`
        exp_scaled(Dd::new(-incbet_lnpow(p.max(q), p.min(q), 0.5, 0.5) - LN_2))
    };
    let nodes = if a == b { Symmetric } else { Plain };
    ldexp_weights(jacobi(n, a, b, m, nodes), e)
}

/// Panics if `n` is zero, or if `alpha <= -1` or `beta <= -1`, like the `ValueError` in SciPy
fn check_jacobi(n: usize, alpha: f64, beta: f64) {
    check_n(n);
    assert!(
        alpha > -1.0 && beta > -1.0,
        "alpha and beta must be greater than -1"
    );
}

/// The Gauss-Jacobi rule, with the weights summing to `mu0`
///
/// Unlike in SciPy, `2 + a + b` is `(a + 1) + (b + 1)` and `b^2 - a^2` is `(b - a) (b + a)`, which
/// don't cancel for `a, b` near -1.
fn jacobi(n: usize, a: f64, b: f64, mu0: f64, nodes: Nodes) -> (Vec<f64>, Vec<f64>) {
    let ab2 = (a + 1.0) + (b + 1.0);
    let an = |k: f64| {
        if k == 0.0 {
            (b - a) / ab2
        } else {
            (b - a) * (b + a) / ((2.0 * k - 2.0 + ab2) * (2.0 * k + ab2))
        }
    };
    let bn = |k: f64| {
        let r = if k == 1.0 {
            1.0
        } else {
            (k * (k - 2.0 + ab2) / (2.0 * k - 3.0 + ab2)).sqrt()
        };
        let s = ((k - 1.0 + (a + 1.0)) * (k - 1.0 + (b + 1.0)) / (2.0 * k - 1.0 + ab2)).sqrt();
        2.0 / (2.0 * k - 2.0 + ab2) * s * r
    };
    let f = |n, x| eval_jacobi(n, a, b, x);
    let df = |n, x| 0.5 * (f64::from(n) - 1.0 + ab2) * eval_jacobi(n - 1, a + 1.0, b + 1.0, x);
    gen_roots_and_weights(
        n,
        mu0,
        an,
        bn,
        f,
        df,
        one_minus_x2,
        |x| b - a - ab2 * x,
        nodes,
    )
}

/// Gauss-Jacobi (shifted) quadrature
///
/// Has no corresponding function in `scipy.special`, see [`roots_gengegenbauer`] instead.
///
/// The nodes are the roots of the shifted Jacobi polynomial $P_n^{*(\alpha,\beta)}$, and the rule
/// is exact for polynomials of degree $2n - 1$ with the weight function
/// $(1 - x)^\alpha x^\beta$ on $[0, 1]$.
///
/// # Panics
/// - if `n` is zero
/// - if $\alpha \le -1$ or $\beta \le -1$
///
/// # See also
/// - [`eval_sh_jacobi`](crate::eval_sh_jacobi): Evaluate shifted Jacobi polynomials,
///   $P_n^{*(\alpha,\beta)}$
/// - [`roots_jacobi`]: Gauss-Jacobi quadrature
#[must_use]
pub fn roots_sh_jacobi(n: usize, alpha: f64, beta: f64) -> (Vec<f64>, Vec<f64>) {
    check_jacobi(n, alpha, beta);
    // the weights sum to `B(alpha + 1, beta + 1)`, without the factor `2^(alpha + beta + 1)` of
    // `roots_jacobi`, which can overflow
    let mu0 = crate::beta(alpha + 1.0, beta + 1.0);
    jacobi(n, alpha, beta, mu0, Shifted)
}

/// Gauss-Jacobi quadrature for the generalized Gegenbauer polynomials
///
/// Corresponds to [`scipy.special.roots_sh_jacobi`][docs], which, despite its name, isn't
/// [`roots_sh_jacobi`] (see [`eval_gengegenbauer`](crate::eval_gengegenbauer)).
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_jacobi.html
///
/// The nodes are the roots of $G_n^{(p,q)}$, and the rule is exact for polynomials of degree
/// $2n - 1$ with the weight function $(1 - x)^{p-q} x^{q-1}$ on $[0, 1]$ (A&S 22.2.2).
///
/// For small $q$, the nodes and weights have relative errors up to about $2^{-54} / q$, from the
/// rounding of $q - 1$ (also in SciPy).
///
/// # Panics
/// - if `n` is zero
/// - if $p - q \le -1$ or $q \le 0$, or if $q - 1$ rounds to $-1$ (for $q < 2^{-54}$)
///
/// # See also
/// - [`eval_gengegenbauer`](crate::eval_gengegenbauer): Evaluate generalized Gegenbauer
///   polynomials, $G_n^{(p,q)}$
#[must_use]
pub fn roots_gengegenbauer(n: usize, p: f64, q: f64) -> (Vec<f64>, Vec<f64>) {
    assert!(
        p - q > -1.0 && q > 0.0,
        "(p - q) must be greater than -1, and q must be greater than 0"
    );
    roots_sh_jacobi(n, p - q, q - 1.0)
}

/// Gauss-Gegenbauer quadrature
///
/// Corresponds to [`scipy.special.roots_gegenbauer`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_gegenbauer.html
///
/// The nodes are the roots of the Gegenbauer polynomial $C_n^{(\alpha)}$, and the rule is exact
/// for polynomials of degree $2n - 1$ with the weight function $(1 - x^2)^{\alpha - 1/2}$ on
/// $[-1, 1]$ (A&S 22.2.3). For $\alpha = 0$, these are the roots of $T_n$, see
/// [`roots_chebyshev_t`].
///
/// The weights are NaN if $C_n^{(\alpha)}$ overflows at all nodes, e.g. for $n = 150$ and
/// $\alpha = 10^4$ (also in SciPy). Near $\alpha = -1/2$, they lose accuracy (also in SciPy): for
/// odd $n$, the relative error is about $n \cdot 10^{-16} / (\alpha + 1/2)$, and within about
/// $10^{-15}$ of $-1/2$, some can be NaN. [`roots_jacobi`] with $\alpha - 1/2$ for both
/// parameters doesn't have this problem.
///
/// # Panics
/// - if `n` is zero
/// - if $\alpha \le -1/2$
///
/// # See also
/// - [`eval_gegenbauer`]: Evaluate Gegenbauer polynomials, $C_n^{(\alpha)}$
#[must_use]
pub fn roots_gegenbauer(n: usize, alpha: f64) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    assert!(alpha > -0.5, "alpha must be greater than -0.5");
    // Unlike in SciPy, a subnormal alpha is treated as 0, where `C_n` at the central root would
    // underflow.
    if alpha.abs() < f64::MIN_POSITIVE {
        // C(n,0,x) == 0 uniformly, however, as alpha->0, C(n,alpha,x)->T(n,x)
        return roots_chebyshev_t(n);
    }

    // Unlike in SciPy, `B(1/2, alpha + 1/2)` doesn't overflow, so it needs no approximation for
    // large alpha.
    let mu0 = crate::beta(0.5, alpha + 0.5);
    let an = |_| 0.0;
    let bn = |k: f64| {
        if k == 1.0 {
            // the general case, without the cancelling factors `2 alpha` for small alpha
            (1.0 / (2.0 * (1.0 + alpha))).sqrt()
        } else {
            (k * (k + 2.0 * alpha - 1.0) / (4.0 * (k + alpha) * (k + alpha - 1.0))).sqrt()
        }
    };
    let f = |n, x| eval_gegenbauer(n, alpha, x);
    let df = |n: i32, x: f64| {
        let c = f64::from(n);
        (-c * x * f(n, x) + ((c - 1.0) + 2.0 * alpha) * f(n - 1, x)) / one_minus_x2(x, 0.0)
    };
    gen_roots_and_weights(
        n,
        mu0,
        an,
        bn,
        f,
        df,
        one_minus_x2,
        |x| -(2.0 * alpha + 1.0) * x,
        Symmetric,
    )
}

/// Gauss-Legendre quadrature
///
/// Corresponds to [`scipy.special.roots_legendre`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_legendre.html
///
/// The nodes are the roots of the Legendre polynomial $P_n$, and the rule is exact for
/// polynomials of degree $2n - 1$ on $[-1, 1]$ (A&S 22.2.10).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_legendre`]: Evaluate Legendre polynomials, $P_n$
#[must_use]
pub fn roots_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    let an = |_| 0.0;
    let bn = |k: f64| k * (1.0 / (4.0 * k * k - 1.0)).sqrt();
    let f = |n, x| eval_legendre(n, x);
    let df = |n: i32, x: f64| {
        let c = f64::from(n);
        (-c * x * f(n, x) + c * f(n - 1, x)) / one_minus_x2(x, 0.0)
    };
    gen_roots_and_weights(n, 2.0, an, bn, f, df, one_minus_x2, |x| -2.0 * x, Symmetric)
}

/// Gauss-Legendre (shifted) quadrature
///
/// Corresponds to [`scipy.special.roots_sh_legendre`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_legendre.html
///
/// The nodes are the roots of the shifted Legendre polynomial $P_n^*$, and the rule is exact for
/// polynomials of degree $2n - 1$ on $[0, 1]$ (A&S 22.2.11).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_sh_legendre`](crate::eval_sh_legendre): Evaluate shifted Legendre polynomials, $P_n^*$
#[must_use]
pub fn roots_sh_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    // Unlike in SciPy, the nodes near 0 aren't `(x + 1) / 2` with the nodes `x` of
    // `roots_legendre`, which cancels.
    roots_sh_jacobi(n, 0.0, 0.0)
}

/// Gauss-Chebyshev (first kind) quadrature
///
/// Corresponds to [`scipy.special.roots_chebyt`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebyt.html
///
/// The nodes are the roots of the Chebyshev polynomial $T_n$, and the rule is exact for
/// polynomials of degree $2n - 1$ with the weight function $1 / \sqrt{1 - x^2}$ on $[-1, 1]$
/// (A&S 22.2.4).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_chebyshev_t`](crate::eval_chebyshev_t): Evaluate Chebyshev polynomials of the first
///   kind, $T_n$
#[doc(alias = "roots_chebyt")]
#[must_use]
pub fn roots_chebyshev_t(n: usize) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    let m = n as f64;
    // `sin(pi k / (2n))` for `k = 1 - n, 3 - n, ..., n - 1`
    let x = (0..n).map(|i| crate::sinpi((2.0 * i as f64 + 1.0 - m) / (2.0 * m)));
    (x.collect(), vec![PI / m; n])
}
