//! Translated from `scipy/special/_orthogonal.py` (scipy/scipy@7d279ba)
//!
//! - <https://github.com/scipy/scipy/blob/7d279ba/scipy/special/_orthogonal.py>
//!
//! The nodes are in ascending order, and the weights sum to the integral of the weight function.

use core::f64::consts::{LN_2, PI, SQRT_2};

use super::orthogonal_eval::{
    eval_gegenbauer, eval_genlaguerre, eval_hermite_h, eval_hermite_he, eval_jacobi, eval_legendre,
};
use super::tridiag::eigvals;
use crate::utils::{MAXLOG, SQRT_PI, SQRT_TAU};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, LOG2, ldexp_f64, two_sum};
use crate::xsf::cephes::incbet_lnpow;
use crate::xsf::ln_gamma;
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

/// The rule for `f(x / s)` from the one for `f(x)`
fn scale((x, w): (Vec<f64>, Vec<f64>), s: f64) -> (Vec<f64>, Vec<f64>) {
    let x = x.into_iter().map(|xi| xi * s).collect();
    let w = w.into_iter().map(|wi| wi * s).collect();
    (x, w)
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

/// Gauss-Chebyshev (second kind) quadrature
///
/// Corresponds to [`scipy.special.roots_chebyu`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebyu.html
///
/// The nodes are the roots of the Chebyshev polynomial $U_n$, and the rule is exact for
/// polynomials of degree $2n - 1$ with the weight function $\sqrt{1 - x^2}$ on $[-1, 1]$
/// (A&S 22.2.5).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_chebyshev_u`](crate::eval_chebyshev_u): Evaluate Chebyshev polynomials of the second
///   kind, $U_n$
#[doc(alias = "roots_chebyu")]
#[must_use]
pub fn roots_chebyshev_u(n: usize) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    let m = n as f64 + 1.0;
    // Unlike in SciPy, `cos(t)` and `sin(t)` with `t = pi k / (n + 1)` are `sinpi` of arguments in
    // `[-1/2, 1/2]`, as in `roots_chebyshev_t`, which is more accurate near `t = 0` and `t = pi`.
    let rule = (1..=n).rev().map(|k| {
        let sin = crate::sinpi(k.min(n + 1 - k) as f64 / m);
        (
            crate::sinpi((m - 2.0 * k as f64) / (2.0 * m)),
            PI * sin * sin / m,
        )
    });
    rule.unzip()
}

/// Gauss-Chebyshev (first kind) quadrature on $[-2, 2]$
///
/// Corresponds to [`scipy.special.roots_chebyc`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebyc.html
///
/// The nodes are the roots of $C_n(x) = 2 T_n(x / 2)$, and the rule is exact for polynomials of
/// degree $2n - 1$ with the weight function $1 / \sqrt{1 - (x/2)^2}$ on $[-2, 2]$ (A&S 22.2.6).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_chebyshev_c`](crate::eval_chebyshev_c): Evaluate $C_n$
#[doc(alias = "roots_chebyc")]
#[must_use]
pub fn roots_chebyshev_c(n: usize) -> (Vec<f64>, Vec<f64>) {
    scale(roots_chebyshev_t(n), 2.0)
}

/// Gauss-Chebyshev (second kind) quadrature on $[-2, 2]$
///
/// Corresponds to [`scipy.special.roots_chebys`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebys.html
///
/// The nodes are the roots of $S_n(x) = U_n(x / 2)$, and the rule is exact for polynomials of
/// degree $2n - 1$ with the weight function $\sqrt{1 - (x/2)^2}$ on $[-2, 2]$ (A&S 22.2.7).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_chebyshev_s`](crate::eval_chebyshev_s): Evaluate $S_n$
#[doc(alias = "roots_chebys")]
#[must_use]
pub fn roots_chebyshev_s(n: usize) -> (Vec<f64>, Vec<f64>) {
    scale(roots_chebyshev_u(n), 2.0)
}

/// Gauss-Chebyshev (first kind, shifted) quadrature
///
/// Corresponds to [`scipy.special.roots_sh_chebyt`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_chebyt.html
///
/// The nodes are the roots of $T_n^*(x) = T_n(2x - 1)$, and the rule is exact for polynomials of
/// degree $2n - 1$ with the weight function $1 / \sqrt{x - x^2}$ on $[0, 1]$ (A&S 22.2.8).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_sh_chebyshev_t`](crate::eval_sh_chebyshev_t): Evaluate $T_n^*$
#[doc(alias = "roots_sh_chebyt")]
#[must_use]
pub fn roots_sh_chebyshev_t(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (_, w) = roots_chebyshev_t(n);
    let m = n as f64;
    // Unlike in SciPy, the nodes are `sin(t / 2 + pi / 4)^2` instead of `(1 + sin(t)) / 2`, with
    // `sin(t)` the nodes of `roots_chebyshev_t`, which cancels near 0.
    let x = (0..n).map(|i| crate::sinpi((2.0 * i as f64 + 1.0) / (4.0 * m)).powi(2));
    (x.collect(), w)
}

/// Gauss-Chebyshev (second kind, shifted) quadrature
///
/// Corresponds to [`scipy.special.roots_sh_chebyu`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_chebyu.html
///
/// The nodes are the roots of $U_n^*(x) = U_n(2x - 1)$, and the rule is exact for polynomials of
/// degree $2n - 1$ with the weight function $\sqrt{x - x^2}$ on $[0, 1]$ (A&S 22.2.9).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_sh_chebyshev_u`](crate::eval_sh_chebyshev_u): Evaluate $U_n^*$
#[doc(alias = "roots_sh_chebyu")]
#[must_use]
pub fn roots_sh_chebyshev_u(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (_, w) = roots_chebyshev_u(n);
    let m = n as f64 + 1.0;
    // Unlike in SciPy, the nodes are `cos(t / 2)^2` instead of `(1 + cos(t)) / 2`, with `cos(t)`
    // the nodes of `roots_chebyshev_u`, which cancels near 0.
    let x = (1..=n).map(|k| crate::sinpi(k as f64 / (2.0 * m)).powi(2));
    // the weights sum to `B(3/2, 3/2) = pi/8` instead of `pi/2`
    (x.collect(), w.into_iter().map(|wi| wi / 4.0).collect())
}

/// Gauss-generalized Laguerre quadrature
///
/// Corresponds to [`scipy.special.roots_genlaguerre`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_genlaguerre.html
///
/// The nodes are the roots of the generalized Laguerre polynomial $L_n^{(\alpha)}$, and the rule
/// is exact for polynomials of degree $2n - 1$ with the weight function $x^\alpha e^{-x}$ on
/// $[0, \infty)$ (A&S 22.3.9).
///
/// # Panics
/// - if `n` is zero
/// - if $\alpha \le -1$
///
/// # See also
/// - [`eval_genlaguerre`]: Evaluate generalized Laguerre polynomials, $L_n^{(\alpha)}$
#[must_use]
pub fn roots_genlaguerre(n: usize, alpha: f64) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    assert!(alpha > -1.0, "alpha must be greater than -1");

    let mu0 = crate::gamma(alpha + 1.0);
    if n == 1 {
        return (vec![alpha + 1.0], vec![mu0]);
    }

    let an = |k: f64| 2.0 * k + alpha + 1.0;
    let bn = |k: f64| -(k * (k + alpha)).sqrt();
    let f = |n, x| eval_genlaguerre(n, alpha, x);
    // Unlike in SciPy, the derivative is `-L_{n-1}^(alpha+1)` (DLMF 18.9.23), which doesn't cancel
    // for alpha near -1.
    let df = |n: i32, x: f64| -eval_genlaguerre(n - 1, alpha + 1.0, x);

    // Unlike in SciPy, the weights can be finite where `mu0` overflows (alpha > 170.6), with
    // `ln(mu0)` in double-double.
    let (mu0, e) = if mu0.is_finite() {
        (mu0, 0.0)
    } else {
        ln_gamma(Dd::new(alpha) + 1.0).map_or((f64::NAN, 0.0), |(l, _, _)| exp_scaled(l))
    };
    let rule = gen_roots_and_weights(
        n,
        mu0,
        an,
        bn,
        f,
        df,
        |x, dx| x + dx,
        |x| alpha + 1.0 - x,
        Plain,
    );
    ldexp_weights(rule, e)
}

/// Gauss-Laguerre quadrature
///
/// Corresponds to [`scipy.special.roots_laguerre`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_laguerre.html
///
/// The nodes are the roots of the Laguerre polynomial $L_n$, and the rule is exact for
/// polynomials of degree $2n - 1$ with the weight function $e^{-x}$ on $[0, \infty)$
/// (A&S 22.2.13).
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_laguerre`](crate::eval_laguerre): Evaluate Laguerre polynomials, $L_n$
#[must_use]
pub fn roots_laguerre(n: usize) -> (Vec<f64>, Vec<f64>) {
    roots_genlaguerre(n, 0.0)
}

/// The largest `n` for which the Hermite rules use [`gen_roots_and_weights`] instead of
/// [`roots_hermite_asy`], like in SciPy
const HERMITE_GW_MAX: usize = 150;

/// Gauss-Hermite (physicist's) quadrature
///
/// Corresponds to [`scipy.special.roots_hermite`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_hermite.html
///
/// The nodes are the roots of the Hermite polynomial $H_n$, and the rule is exact for polynomials
/// of degree $2n - 1$ with the weight function $e^{-x^2}$ on $(-\infty, \infty)$ (A&S 22.2.14).
/// For $n > 150$, they are computed by the asymptotic method of Townsend, Trogdon, and Olver
/// [^TTO].
///
/// [^TTO]: A. Townsend, T. Trogdon, and S. Olver, "Fast computation of Gauss quadrature nodes and
///   weights on the whole real line", IMA J. Numer. Anal. 36 (2016),
///   <https://arxiv.org/abs/1410.5286>
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_hermite_h`]: Evaluate physicist's Hermite polynomials, $H_n$
#[doc(alias = "roots_hermite")]
#[must_use]
pub fn roots_hermite_h(n: usize) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    if n > HERMITE_GW_MAX {
        return roots_hermite_asy(n);
    }
    let an = |_| 0.0;
    let bn = |k: f64| (k / 2.0).sqrt();
    let f = |n: i32, x| eval_hermite_h(n as u32, x);
    let df = |n: i32, x| 2.0 * f64::from(n) * f(n - 1, x);
    let tau = |x: f64| -2.0 * x;
    gen_roots_and_weights(n, SQRT_PI, an, bn, f, df, |_, _| 1.0, tau, Symmetric)
}

/// Gauss-Hermite (statistician's) quadrature
///
/// Corresponds to [`scipy.special.roots_hermitenorm`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_hermitenorm.html
///
/// The nodes are the roots of the Hermite polynomial $He_n$, and the rule is exact for
/// polynomials of degree $2n - 1$ with the weight function $e^{-x^2/2}$ on $(-\infty, \infty)$
/// (A&S 22.2.15). For $n > 150$, they are scaled from [`roots_hermite_h`].
///
/// # Panics
/// - if `n` is zero
///
/// # See also
/// - [`eval_hermite_he`]: Evaluate probabilist's Hermite polynomials, $He_n$
#[doc(alias = "roots_hermitenorm")]
#[must_use]
pub fn roots_hermite_he(n: usize) -> (Vec<f64>, Vec<f64>) {
    check_n(n);
    if n > HERMITE_GW_MAX {
        return scale(roots_hermite_asy(n), SQRT_2);
    }
    let an = |_| 0.0;
    let bn = |k: f64| k.sqrt();
    let f = |n: i32, x| eval_hermite_he(n as u32, x);
    let df = |n: i32, x| f64::from(n) * f(n - 1, x);
    gen_roots_and_weights(n, SQRT_TAU, an, bn, f, df, |_, _| 1.0, |x| -x, Symmetric)
}

/// Corresponds to `_compute_tauk`: the root of `tau - sin(tau) = c` of Lemma 3.1 in [^TTO], by
/// 5 Newton steps
///
/// [^TTO]: see [`roots_hermite_h`]
fn compute_tauk(n: usize, k: f64) -> f64 {
    let a = (n % 2) as f64 - 0.5;
    let h = (n / 2) as f64;
    let c = (4.0 * h - 4.0 * k + 3.0) * PI / (4.0 * h + 2.0 * a + 2.0);
    let mut xi = 0.5 * PI;
    for _ in 0..5 {
        xi -= (xi - xi.sin() - c) / (1.0 - xi.cos());
    }
    xi
}

/// Corresponds to `_initial_nodes_a`: Tricomi's approximation of the square of the `k`-th
/// positive root of `H_n` (Lemma 3.1 in [^TTO]), except near `sqrt(2n + 1)`
///
/// [^TTO]: see [`roots_hermite_h`]
fn initial_nodes_a(n: usize, k: f64) -> f64 {
    let tauk = compute_tauk(n, k);
    let sigk = (0.5 * tauk).cos().powi(2);
    let a = (n % 2) as f64 - 0.5;
    let nu = 4.0 * (n / 2) as f64 + 2.0 * a + 2.0;
    // Initial approximation of Hermite roots (square)
    nu * sigk - 1.0 / (3.0 * nu) * (5.0 / (4.0 * (1.0 - sigk).powi(2)) - 1.0 / (1.0 - sigk) - 0.25)
}

/// Corresponds to `_initial_nodes_b`: Gatteschi's approximation of the square of the `k`-th
/// positive root of `H_n` from the top (Lemma 3.2 in [^TTO]), given the `k`-th zero `ak` of Ai
///
/// [^TTO]: see [`roots_hermite_h`]
fn initial_nodes_b(n: usize, ak: f64) -> f64 {
    let a = (n % 2) as f64 - 0.5;
    let nu = 4.0 * (n / 2) as f64 + 2.0 * a + 2.0;
    // Initial approximation of Hermite roots (square)
    nu + 2.0_f64.powf(2.0 / 3.0) * ak * nu.powf(1.0 / 3.0)
        + 1.0 / 5.0 * 2.0_f64.powf(4.0 / 3.0) * ak.powi(2) * nu.powf(-1.0 / 3.0)
        + (9.0 / 140.0 - 12.0 / 175.0 * ak.powi(3)) * nu.powi(-1)
        + (16.0 / 1575.0 * ak + 92.0 / 7875.0 * ak.powi(4))
            * 2.0_f64.powf(2.0 / 3.0)
            * nu.powf(-5.0 / 3.0)
        - (15152.0 / 3031875.0 * ak.powi(5) + 1088.0 / 121275.0 * ak.powi(2))
            * 2.0_f64.powf(1.0 / 3.0)
            * nu.powf(-7.0 / 3.0)
}

/// Corresponds to `_initial_nodes`: approximations of the non-negative roots of `H_n`, by
/// [`initial_nodes_a`] up to the turnover point, and [`initial_nodes_b`] above it
fn initial_nodes(n: usize) -> Vec<f64> {
    // Turnover point
    // linear polynomial fit to error of 10, 25, 40, ..., 1000 point rules
    let fit = 0.49082003 * n as f64 - 4.37859653;
    let turnover = fit.round_ties_even() as usize;
    let h = n / 2;
    // the indices `k` of `initial_nodes_b`, from `kmax` down to 1
    let kmax = h.saturating_sub(turnover + 1);
    let ak = crate::xsf::airy_ai_zeros_vec(kmax);

    // Central node is always zero
    let zero = (n % 2 == 1).then_some(0.0);
    let xa = (1..=h.min(turnover + 1)).map(|k| initial_nodes_a(n, k as f64));
    let xb = ak.into_iter().rev().map(|ak| initial_nodes_b(n, ak));
    zero.into_iter()
        .chain(xa.chain(xb).map(f64::sqrt))
        .collect()
}

/// Corresponds to `_pbcf`: the parabolic cylinder function `U(a, x)` (DLMF 12.10.35) and its
/// derivative (DLMF 12.10.36), with `a = -n - 1/2` and `x = sqrt(2n + 1) cos(theta)`, by its
/// asymptotic expansion in Airy functions (Sections 3.2 and 3.3 of [^TTO], with one more term)
///
/// [^TTO]: see [`roots_hermite_h`]
#[rustfmt::skip]
fn pbcf(n: usize, theta: f64) -> (f64, f64) {
    let (st, ct) = theta.sin_cos();
    // https://dlmf.nist.gov/12.10#vii
    let mu = 2.0 * n as f64 + 1.0;
    // https://dlmf.nist.gov/12.10#E23
    let eta = 0.5 * theta - 0.5 * st * ct;
    // https://dlmf.nist.gov/12.10#E39
    let zeta = -(3.0 * eta / 2.0).powf(2.0 / 3.0);
    // https://dlmf.nist.gov/12.10#E40
    let phi = (-zeta / st.powi(2)).powf(0.25);
    // Coefficients
    // https://dlmf.nist.gov/12.10#E43
    let a0 = 1.0;
    let a1 = 0.10416666666666666667;
    let a2 = 0.08355034722222222222;
    let a3 = 0.12822657455632716049;
    let a4 = 0.29184902646414046425;
    let a5 = 0.88162726744375765242;
    let b0 = 1.0;
    let b1 = -0.14583333333333333333;
    let b2 = -0.09874131944444444444;
    let b3 = -0.14331205391589506173;
    let b4 = -0.31722720267841354810;
    let b5 = -0.94242914795712024914;
    // Polynomials
    // https://dlmf.nist.gov/12.10#E9
    // https://dlmf.nist.gov/12.10#E10
    let ctp = |k: i32| ct.powi(k);
    let u0 = 1.0;
    let u1 = (1.0 * ctp(3) - 6.0 * ct) / 24.0;
    let u2 = (-9.0 * ctp(4) + 249.0 * ctp(2) + 145.0) / 1152.0;
    let u3 = (-4042.0 * ctp(9) + 18189.0 * ctp(7) - 28287.0 * ctp(5)
        - 151995.0 * ctp(3) - 259290.0 * ct) / 414720.0;
    let u4 = (72756.0 * ctp(10) - 321339.0 * ctp(8) - 154982.0 * ctp(6)
        + 50938215.0 * ctp(4) + 122602962.0 * ctp(2) + 12773113.0) / 39813120.0;
    let u5 = (82393456.0 * ctp(15) - 617950920.0 * ctp(13) + 1994971575.0 * ctp(11)
        - 3630137104.0 * ctp(9) + 4433574213.0 * ctp(7) - 37370295816.0 * ctp(5)
        - 119582875013.0 * ctp(3) - 34009066266.0 * ct) / 6688604160.0;
    let v0 = 1.0;
    let v1 = (1.0 * ctp(3) + 6.0 * ct) / 24.0;
    let v2 = (15.0 * ctp(4) - 327.0 * ctp(2) - 143.0) / 1152.0;
    let v3 = (-4042.0 * ctp(9) + 18189.0 * ctp(7) - 36387.0 * ctp(5)
        + 238425.0 * ctp(3) + 259290.0 * ct) / 414720.0;
    let v4 = (-121260.0 * ctp(10) + 551733.0 * ctp(8) - 151958.0 * ctp(6)
        - 57484425.0 * ctp(4) - 132752238.0 * ctp(2) - 12118727.0) / 39813120.0;
    let v5 = (82393456.0 * ctp(15) - 617950920.0 * ctp(13) + 2025529095.0 * ctp(11)
        - 3750839308.0 * ctp(9) + 3832454253.0 * ctp(7) + 35213253348.0 * ctp(5)
        + 130919230435.0 * ctp(3) + 34009066266.0 * ct) / 6688604160.0;
    // Airy Evaluation (Bi and Bip unused)
    let [ai, aip, _, _] = crate::airy(mu.powf(4.0 / 6.0) * zeta);
    // Prefactor for U
    let p = 2.0 * SQRT_PI * mu.powf(1.0 / 6.0) * phi;
    // Terms for U
    // https://dlmf.nist.gov/12.10#E42
    let phip = |k: i32| phi.powi(6 * k);
    let a_0 = b0 * u0;
    let a_1 = (b2 * u0 + phip(1) * b1 * u1 + phip(2) * b0 * u2) / zeta.powi(3);
    let a_2 = (b4 * u0 + phip(1) * b3 * u1 + phip(2) * b2 * u2 + phip(3) * b1 * u3
        + phip(4) * b0 * u4) / zeta.powi(6);
    let b_0 = -(a1 * u0 + phip(1) * a0 * u1) / zeta.powi(2);
    let b_1 = -(a3 * u0 + phip(1) * a2 * u1 + phip(2) * a1 * u2 + phip(3) * a0 * u3)
        / zeta.powi(5);
    let b_2 = -(a5 * u0 + phip(1) * a4 * u1 + phip(2) * a3 * u2 + phip(3) * a2 * u3
        + phip(4) * a1 * u4 + phip(5) * a0 * u5) / zeta.powi(8);
    // U
    // https://dlmf.nist.gov/12.10#E35
    let u = p * (ai * (a_0 + a_1 / mu.powi(2) + a_2 / mu.powi(4))
        + aip * (b_0 + b_1 / mu.powi(2) + b_2 / mu.powi(4)) / mu.powf(8.0 / 6.0));
    // Prefactor for derivative of U
    let pd = SQRT_TAU * mu.powf(2.0 / 6.0) / phi;
    // Terms for derivative of U
    // https://dlmf.nist.gov/12.10#E46
    let c_0 = -(b1 * v0 + phip(1) * b0 * v1) / zeta;
    let c_1 = -(b3 * v0 + phip(1) * b2 * v1 + phip(2) * b1 * v2 + phip(3) * b0 * v3)
        / zeta.powi(4);
    let c_2 = -(b5 * v0 + phip(1) * b4 * v1 + phip(2) * b3 * v2 + phip(3) * b2 * v3
        + phip(4) * b1 * v4 + phip(5) * b0 * v5) / zeta.powi(7);
    let d_0 = a0 * v0;
    let d_1 = (a2 * v0 + phip(1) * a1 * v1 + phip(2) * a0 * v2) / zeta.powi(3);
    let d_2 = (a4 * v0 + phip(1) * a3 * v1 + phip(2) * a2 * v2 + phip(3) * a1 * v3
        + phip(4) * a0 * v4) / zeta.powi(6);
    // Derivative of U
    // https://dlmf.nist.gov/12.10#E36
    let ud = pd * (ai * (c_0 + c_1 / mu.powi(2) + c_2 / mu.powi(4)) / mu.powf(4.0 / 6.0)
        + aip * (d_0 + d_1 / mu.powi(2) + d_2 / mu.powi(4)));
    (u, ud)
}

/// Corresponds to `_newton`: at most 5 Newton steps in `theta` (see [`pbcf`]) from the
/// approximate non-negative roots `x` of `H_n`, and the weights of those roots
fn newton(n: usize, x: &[f64]) -> (Vec<f64>, Vec<f64>) {
    // Variable transformation
    let mu = (2.0 * n as f64 + 1.0).sqrt();
    let mut theta: Vec<f64> = x.iter().map(|xi| (xi / mu).acos()).collect();
    let mut ud = Vec::new();
    // Newton iteration
    for _ in 0..5 {
        let (u, d): (Vec<f64>, Vec<f64>) = theta.iter().map(|&t| pbcf(n, t)).unzip();
        let mut max = 0.0_f64;
        for ((t, ui), udi) in theta.iter_mut().zip(u).zip(&d) {
            let dtheta = ui / (SQRT_2 * mu * t.sin() * udi);
            *t += dtheta;
            max = max.max(dtheta.abs());
        }
        ud = d;
        if max < 1e-14 {
            break;
        }
    }
    // Undo variable transformation
    let mut x: Vec<f64> = theta.iter().map(|t| mu * t.cos()).collect();
    // Central node is always zero
    if n % 2 == 1 {
        x[0] = 0.0;
    }
    // Compute weights
    let w: Vec<f64> = x
        .iter()
        .zip(ud)
        .map(|(xi, udi)| (-xi * xi).exp() / (2.0 * udi * udi))
        .collect();
    (x, w)
}

/// Corresponds to `_roots_hermite_asy`: [`roots_hermite_h`] for large `n`, in linear time
fn roots_hermite_asy(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (nodes, weights) = newton(n, &initial_nodes(n));
    // Combine with negative parts (without the central zero of odd n)
    let skip = n % 2;
    let x = nodes
        .iter()
        .skip(skip)
        .rev()
        .map(|x| -x)
        .chain(nodes.iter().copied());
    let w: Vec<f64> = weights
        .iter()
        .skip(skip)
        .rev()
        .chain(&weights)
        .copied()
        .collect();
    // Scale weights
    let scale = SQRT_PI / w.iter().sum::<f64>();
    (x.collect(), w.into_iter().map(|wi| wi * scale).collect())
}
