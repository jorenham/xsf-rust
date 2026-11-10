//! Translated from `xsf/orthogonal_eval.h` and `xsf/cpu/orthogonal_eval.h` (scipy/xsf@e1800e8).
//!
//! - <https://github.com/scipy/xsf/blob/e1800e8/include/xsf/orthogonal_eval.h>
//! - <https://github.com/scipy/xsf/blob/e1800e8/include/xsf/cpu/orthogonal_eval.h>

use core::f64::consts::SQRT_2;
use core::ops::Range;
use num_traits::Zero;

use crate::utils::{MAXGAM, SQRT_EPS, SQRT_PI};
use crate::xsf::cephes::{HYP2F1_EPS, HYP2F1_ETHRESH};

/// The degrees `n` with an exact `n - 1`
const EXACT_N_MINUS_1: core::ops::RangeInclusive<f64> = 1.0..=2.0 / f64::EPSILON;

mod sealed {
    pub trait OrthoPolyArg:
        crate::xsf::Hyp2F1Arg
        + crate::xsf::Hyp1F1Arg
        + num_traits::NumOps<f64>
        + num_traits::NumOps
        + Clone
    {
        const ZERO: Self;

        fn is_finite(&self) -> bool;
        fn is_nan(&self) -> bool;

        /// `Some(x)` if `Self` is `f64`
        fn as_f64(&self) -> Option<f64>;

        /// `2F1(-n, b; c; (1 - x) / 2)`, where `b` is `n + lam` (rounded)
        #[inline]
        fn hyp2f1_sh(self, n: f64, b: f64, c: f64, _lam: f64) -> Self {
            (self.clone() * -0.5 + 0.5).hyp2f1(-n, b, c)
        }
    }

    impl OrthoPolyArg for f64 {
        const ZERO: Self = 0.0;

        #[inline]
        fn is_finite(&self) -> bool {
            (*self).is_finite()
        }

        #[inline]
        fn is_nan(&self) -> bool {
            (*self).is_nan()
        }

        #[inline]
        fn as_f64(&self) -> Option<f64> {
            Some(*self)
        }

        #[inline]
        fn hyp2f1_sh(self, n: f64, b: f64, c: f64, lam: f64) -> Self {
            // Unlike in xsf, use the recurrence in the degree if -1 < x < 1, since the
            // hypergeometric series cancels catastrophically for large non-integer degrees.
            //
            // E.g. `eval_legendre(100.5, 0.3)` was 0.963 instead of 0.01165.
            let nu = n.max(-b); // the 2F1 is symmetric in n and -n - lam
            let y = crate::hyp2f1(-n, b, c, self * -0.5 + 0.5);
            super::jacobi_hyp2f1_recurrence(self, nu, lam, c, y).unwrap_or(y)
        }
    }

    impl OrthoPolyArg for num_complex::Complex<f64> {
        const ZERO: Self = Self::ZERO;

        #[inline]
        fn is_finite(&self) -> bool {
            (*self).is_finite()
        }

        #[inline]
        fn is_nan(&self) -> bool {
            (*self).is_nan()
        }

        #[inline]
        fn as_f64(&self) -> Option<f64> {
            None
        }
    }
}

/// `G_n = 2F1(-n, n + lam; c; (1 - x) / 2)` by forward recurrence in `n` from `n - floor(n)`, or
/// `None` if that isn't more accurate than `series` (`G_n` itself).
///
/// - With `alpha = c - 1` and `beta = lam - c`, this is the Jacobi recurrence (DLMF 18.9.2) of
///   `binom(n + alpha, n) G_n`, also for non-integer `n` (as a contiguous relation of 2F1).
/// - It needs `alpha, beta > -1` (positive denominators), `n > 2` (at least one step), and
///   `n < 1e5` (to bound its cost).
/// - `series` is estimated to be less accurate if it differs by more than relative errors `tol`
///   in the starting values explain: `HYP2F1_ETHRESH`, or even 1 if the recurrence is unstable
///   (the errors grow by more than `1 / SQRT_EPS`, leaving less than half the digits).
fn jacobi_hyp2f1_recurrence(x: f64, n: f64, lam: f64, c: f64, series: f64) -> Option<f64> {
    let (alpha, beta) = (c - 1.0, lam - c);
    if !(x.abs() < 1.0 && n > 2.0 && n < 1e5 && n.fract() != 0.0 && c > 0.0 && beta > -1.0) {
        return None;
    }
    let ab = lam - 1.0;
    // `n - j` is exact for integer `0 <= j <= n`
    let m = n.floor();
    let recur = |mut g0: f64, mut g1: f64| {
        let mut j = m - 2.0;
        while j >= 0.0 {
            let v = n - j;
            let t = 2.0 * v + ab;
            let g2 = ((t - 1.0) * (t * (t - 2.0) * x + (alpha - beta) * ab) * g1
                - 2.0 * (v - 1.0) * (v + beta - 1.0) * t * g0)
                / (2.0 * (v + ab) * (t - 2.0) * (v + alpha));
            (g0, g1) = (g1, g2);
            j -= 1.0;
        }
        g1
    };

    // `lam` is passed separately, so that `c - a - b` is the same as in the caller's 2F1
    let g = |v: f64| crate::hyp2f1(-v, v + lam, c, x * -0.5 + 0.5);
    let (g0, g1) = (g(n - m), g(n - m + 1.0));
    let y = recur(g0, g1);
    // the solutions with starting values (1, 0) and (0, 1)
    let growth = (recur(1.0, 0.0) * g0).abs() + (recur(0.0, 1.0) * g1).abs();
    let stable = growth * SQRT_EPS <= y.abs();
    let tol = if stable { HYP2F1_ETHRESH } else { 1.0 };
    (!series.is_finite() || (y - series).abs() > tol * growth).then_some(y)
}

///////////////////////////////
// Jacobi
// Legendre
//

pub trait JacobiArg<N>: sealed::OrthoPolyArg {
    fn eval_jacobi(self, n: N, alpha: f64, beta: f64) -> Self;
    fn eval_legendre(self, n: N) -> Self;
}

impl<Z: sealed::OrthoPolyArg> JacobiArg<f64> for Z {
    /// Corresponds to `xsf::detail::eval_jacobi`
    #[inline]
    fn eval_jacobi(self, n: f64, a: f64, b: f64) -> Self {
        if a == -1.0 && b.abs() == 1.0 {
            if n == 0.0 {
                return Self::ZERO + 1.0;
            } else if n == 1.0 {
                return (self - 1.0) * (0.5 * (1.0 + b));
            } else if n > 1.0 {
                let p = self.clone().eval_jacobi(n - 1.0, 1.0, b);
                return (self - 1.0) * ((n + b) / (2.0 * n)) * p;
            }
        }
        // Unlike in xsf (scipy/xsf#457), if `|alpha + 1| < HYP2F1_EPS` (where the real `hyp2f1`
        // loses precision), alpha is shifted up by a relation that follows from DLMF 18.5.7 term by
        // term.
        // E.g. `eval_jacobi(7.5, -1 + 1e-15, 2.0, -0.6)` was NaN instead of 1.337.
        if (a + 1.0).abs() < HYP2F1_EPS && EXACT_N_MINUS_1.contains(&n) {
            let c = a + 1.0;
            let p = self.clone().eval_jacobi(n, c, b) * c;
            let q = self.clone().eval_jacobi(n - 1.0, c + 1.0, b) * (self - 1.0) * (0.5 * (n + b));
            return (p + q) / (n + c);
        }
        let d = jacobi_binom(n, a);
        self.hyp2f1_sh(n, n + a + b + 1.0, a + 1.0, a + b + 1.0) * d
    }

    /// Corresponds to `xsf::detail::eval_legendre`
    #[inline]
    fn eval_legendre(self, n: f64) -> Self {
        if self.is_nan() {
            return self;
        }
        self.hyp2f1_sh(n, n + 1.0, 1.0, 1.0)
    }
}

impl JacobiArg<i32> for f64 {
    /// Corresponds to `xsf::detail::eval_jacobi_l`
    #[inline]
    fn eval_jacobi(self, n: i32, alpha: f64, beta: f64) -> Self {
        let x = self;
        if n < 0 {
            x.eval_jacobi(f64::from(n), alpha, beta)
        } else if n == 0 {
            1.0
        } else if n == 1 {
            0.5 * (2.0 * (alpha + 1.0) + (alpha + beta + 2.0) * (x - 1.0))
        } else if alpha == -1.0 {
            // Unlike in xsf, for every beta (not only |beta| = 1), since the recurrence below
            // divides by alpha + 1. This follows from DLMF 18.5.7 term by term.
            let n = f64::from(n);
            ((n + beta) / (2.0 * n)) * (x - 1.0) * x.eval_jacobi(n - 1.0, 1.0, beta)
        } else {
            let mut d = (alpha + beta + 2.0) * (x - 1.0) / (2.0 * (alpha + 1.0));
            let mut p = d + 1.0;
            for kk in 0..n - 1 {
                let k = f64::from(kk) + 1.0;
                let t = 2.0 * k + alpha + beta;
                d = ((t * (t + 1.0) * (t + 2.0)) * (x - 1.0) * p
                    + 2.0 * k * (k + beta) * (t + 2.0) * d)
                    / (2.0 * (k + alpha + 1.0) * (k + alpha + beta + 1.0) * t);
                p += d;
            }
            jacobi_binom(f64::from(n), alpha) * p
        }
    }

    /// Corresponds to `xsf::detail::eval_legendre_l`
    // Integer-degree P_n(x) by forward recurrence on the values p_k = P_k(x) and their
    // differences d_k = p_k - p_{k-1}:
    //   d_{k+1} = ((2k + 1) (x - 1) p_k + k d_k) / (k + 1),  p_{k+1} = p_k + d_{k+1},
    // starting from p_1 = x, d_1 = x - 1. Near x = 0 an explicit power series is used instead.
    #[inline]
    fn eval_legendre(self, n: i32) -> Self {
        let x = self;
        if x.is_nan() {
            return x;
        }

        // P_{-n-1}(x) = P_n(x); add before negating to avoid signed overflow.
        let n = if n < 0 { -(n + 1) } else { n };

        if n == 0 {
            1.0
        } else if n == 1 {
            x
        } else if x.abs() < 1e-5 {
            // Power series rather than recurrence due to loss of precision
            // http://functions.wolfram.com/Polynomials/LegendreP/02/
            let a = n / 2;
            let mut d = if a % 2 == 0 { 1.0 } else { -1.0 };
            if n % 2 == 0 {
                d *= -2.0 / crate::beta(f64::from(a) + 1.0, -0.5);
            } else {
                d *= 2.0 * x / crate::beta(f64::from(a) + 1.0, 0.5);
            }

            let mut p = 0.0;
            for kk in 0..=a {
                let (n, a, kk) = (f64::from(n), f64::from(a), f64::from(kk));
                p += d;
                d *= -2.0 * x * x * (a - kk) * (2.0 * n + 1.0 - 2.0 * a + 2.0 * kk)
                    / ((n + 1.0 - 2.0 * a + 2.0 * kk) * (n + 2.0 - 2.0 * a + 2.0 * kk));
                if d.abs() < 1e-20 * p.abs() {
                    // converged
                    break;
                }
            }
            p
        } else {
            let mut d = x - 1.0;
            let mut p = x;
            for kk in 0..n - 1 {
                let k = f64::from(kk) + 1.0;
                d = ((2.0 * k + 1.0) / (k + 1.0)) * (x - 1.0) * p + (k / (k + 1.0)) * d;
                p += d;
            }
            p
        }
    }
}

/// Evaluate Jacobi polynomial $P_n^{(\alpha, \beta)}$ at a point.
///
/// Corresponds to [`scipy.special.eval_jacobi`][jaceval], and `xsf::eval_jacobi`.
///
/// [jaceval]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_jacobi.html
///
/// # Definition
///
/// The Jacobi polynomials can be defined via the Gauss hypergeometric function $_2F_1$ as
///
/// $$
/// P_n^{(\alpha, \beta)}(z) =
///   {n + \alpha \choose n} \\,
///   \hyp{2}{1}{-n,\enspace n+1+\alpha+\beta}{1+\alpha}{-{z-1 \over 2}}
/// $$
///
/// When $n$ is an integer the result is a polynomial of degree $n$.
/// See Abramowitz & Stegun 22.5.42 [^AS] or DLMF 18.5.7 [^DLMF] for details.
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
///   Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E7>
///
/// # See also
/// - [`eval_legendre`]: Evaluate Legendre polynomials, $P_n$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $_2F_1$
#[inline]
pub fn eval_jacobi<N, Z: JacobiArg<N>>(n: N, alpha: f64, beta: f64, z: Z) -> Z {
    z.eval_jacobi(n, alpha, beta)
}

/// Evaluate Legendre polynomial $P_n$ at a point.
///
/// Corresponds to [`scipy.special.eval_legendre`][legeval], and `xsf::eval_legendre`.
///
/// [legeval]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_legendre.html
///
/// # Definition
///
/// The Legendre polynomials can be defined via the Gauss hypergeometric function $_2F_1$ as
///
/// $$
/// P_n(z) = \hyp{2}{1}{-n,\enspace n+1}{1}{-{z-1 \over 2}}
/// $$
///
/// When $n$ is an integer the result is a polynomial of degree $n$.
/// See Abramowitz & Stegun 22.5.42 [^AS] or DLMF 18.5.7 [^DLMF] for details.
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
///   Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E7>
///
/// # See also
/// - [`eval_jacobi`]: Evaluate Jacobi polynomials, $P_n^{(\alpha, \beta)}$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $_2F_1$
#[inline]
pub fn eval_legendre<N, Z: JacobiArg<N>>(n: N, z: Z) -> Z {
    z.eval_legendre(n)
}

/// Evaluate shifted Jacobi polynomial $P_n^{*(\alpha,\beta)}$ at a point.
///
/// # Definition
///
/// $$ P_n^{*(\alpha,\beta)}(z) = P_n^{(\alpha,\beta)}(2z - 1) $$
///
/// # Note
///
/// This differs from [`scipy.special.eval_sh_jacobi`][docs], which evaluates the generalized
/// Gegenbauer polynomial $G_n^{(p,q)}$ instead, i.e. [`eval_gengegenbauer`] (scipy/scipy#18988).
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_jacobi.html
///
/// # See also
/// - [`eval_jacobi`]: Evaluate Jacobi polynomials, $P_n^{(\alpha, \beta)}$
/// - [`eval_gengegenbauer`]: Evaluate generalized Gegenbauer polynomials, $G_n^{(p,q)}$
#[inline]
pub fn eval_sh_jacobi<N, Z: JacobiArg<N>>(n: N, alpha: f64, beta: f64, z: Z) -> Z {
    (z * 2.0 - 1.0).eval_jacobi(n, alpha, beta)
}

/// Evaluate generalized Gegenbauer polynomial $G_n^{(p,q)}$ at a point.
///
/// Corresponds to [`scipy.special.eval_sh_jacobi`][docs] (and `xsf::eval_sh_jacobi`), which,
/// despite its name, isn't the shifted Jacobi polynomial [`eval_sh_jacobi`], see scipy/scipy#18988.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_jacobi.html
///
/// # Definition
///
/// $$
/// G_n^{(p,q)}(z) = {2n + p - 1 \choose n}^{-1} P_n^{(p-q,\,q-1)}(2z - 1)
/// $$
///
/// See <https://mathworld.wolfram.com/JacobiPolynomial.html>, eq. 1.
///
/// # See also
/// - [`eval_jacobi`]: Evaluate Jacobi polynomials, $P_n^{(\alpha, \beta)}$
/// - [`eval_sh_jacobi`]: Evaluate shifted Jacobi polynomials, $P_n^{*(\alpha,\beta)}$
#[inline]
pub fn eval_gengegenbauer<N: Copy + Into<f64>, Z: JacobiArg<N>>(n: N, p: f64, q: f64, z: Z) -> Z {
    let nf: f64 = n.into();
    eval_sh_jacobi(n, p - q, q - 1.0, z) / crate::binom(2.0 * nf + p - 1.0, nf)
}

/// Evaluate shifted Legendre polynomial $P_n^*$ at a point.
///
/// Corresponds to [`scipy.special.eval_sh_legendre`][docs], and `xsf::eval_sh_legendre`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_legendre.html
///
/// # Definition
///
/// $$ P_n^*(z) = P_n(2z - 1) $$
///
/// # See also
/// - [`eval_legendre`]: Evaluate Legendre polynomials, $P_n$
#[inline]
pub fn eval_sh_legendre<N, Z: JacobiArg<N>>(n: N, z: Z) -> Z {
    (z * 2.0 - 1.0).eval_legendre(n)
}

///////////////////////////////
// Gegenbauer (Ultraspherical)
//

pub trait GegenbauerArg<N>: sealed::OrthoPolyArg {
    fn eval_gegenbauer(self, n: N, alpha: f64) -> Self;
}

impl<Z: sealed::OrthoPolyArg + From<f64>> GegenbauerArg<f64> for Z {
    /// Corresponds to `xsf::detail::eval_gegenbauer`
    #[inline]
    fn eval_gegenbauer(self, n: f64, alpha: f64) -> Self {
        // If n is an integer, use more stable `eval_gegenbauer_l`
        if let Some(x) = self.as_f64() {
            if (0.0..1e10).contains(&n) && n.fract() == 0.0 {
                return Self::from(eval_gegenbauer_l(n as i64, alpha, x));
            }
        }

        // Unlike in xsf, complex x with n = alpha = 0 gives 0, like `eval_gegenbauer_l`, instead
        // of the NaN of the Gamma ratio below (inf / inf).
        if n == 0.0 && alpha == 0.0 && !self.is_nan() {
            return Self::ZERO;
        }

        let mut d =
            crate::gamma(n + 2.0 * alpha) / crate::gamma(1.0 + n) / crate::gamma(2.0 * alpha);
        if !d.is_finite() || (1.0 + n).max(2.0 * alpha) >= MAXGAM {
            // Unlike in xsf (scipy/scipy#25818), fall back to the equivalent binomial if a Gamma
            // function overflows.
            d = crate::binom(n + 2.0 * alpha - 1.0, n);
        }
        self.hyp2f1_sh(n, n + 2.0 * alpha, alpha + 0.5, 2.0 * alpha) * d
    }
}

impl GegenbauerArg<i32> for f64 {
    #[inline]
    fn eval_gegenbauer(self, n: i32, alpha: f64) -> Self {
        eval_gegenbauer_l(i64::from(n), alpha, self)
    }
}

/// Corresponds to `xsf::detail::eval_gegenbauer_l`
// Integer-degree C_n^(alpha)(x). For |x| >= 1e-5, forward recurrence on the normalized values
// p_k = C_k^(alpha)(x) / binom(k + 2 alpha - 1, k) and their differences d_k = p_k - p_{k-1}:
//   d_{k+1} = 2 (k + alpha) / (k + 2 alpha) (x - 1) p_k + k / (k + 2 alpha) d_k,  p_{k+1} = p_k + d_{k+1},
// starting from p_1 = x, d_1 = x - 1. For |x| < 1e-5, the explicit power series in x is used instead.
fn eval_gegenbauer_l(n: i64, alpha: f64, x: f64) -> f64 {
    if alpha.is_nan() || x.is_nan() {
        return f64::NAN;
    }

    if alpha == 0.0 {
        return 0.0;
    }
    if n < 0 {
        return 0.0;
    } else if n == 0 {
        return 1.0;
    } else if n == 1 {
        return 2.0 * alpha * x;
    }

    // Symmetry property: C_n^alpha(-x) = (-1)^n C_n^alpha(x) for integer n
    let sign = if x < 0.0 && n % 2 != 0 { -1.0 } else { 1.0 };
    let x = x.abs();

    if x < 1e-5 {
        // Power series rather than recurrence due to loss of precision
        // http://functions.wolfram.com/Polynomials/GegenbauerC3/02/
        let a = n / 2;

        let mut d = if a % 2 == 0 { 1.0 } else { -1.0 };
        d /= crate::beta(alpha, 1.0 + a as f64);
        if n == 2 * a {
            d /= a as f64 + alpha;
        } else {
            d *= 2.0 * x;
        }

        let mut p = 0.0;
        for kk in 0..=a {
            let (n, a, kk) = (n as f64, a as f64, kk as f64);
            p += d;
            d *= -4.0 * x * x * (a - kk) * (alpha - a + kk + n)
                / ((n + 1.0 - 2.0 * a + 2.0 * kk) * (n + 2.0 - 2.0 * a + 2.0 * kk));
            if d.abs() < 1e-20 * p.abs() {
                // converged
                break;
            }
        }
        return sign * p;
    }

    let mut d = x - 1.0;
    let mut p = x;
    for kk in 0..n - 1 {
        let k = kk as f64 + 1.0;
        d = (2.0 * (k + alpha) / (k + 2.0 * alpha)) * (x - 1.0) * p + (k / (k + 2.0 * alpha)) * d;
        p += d;
    }

    let n = n as f64;
    if (alpha / n).abs() < 1e-8 {
        // avoid loss of precision
        return sign * 2.0 * alpha / n * p;
    }
    sign * crate::binom(n + 2.0 * alpha - 1.0, n) * p
}

/// Evaluate Gegenbauer polynomial $C_n^{(\alpha)}$ at a point.
///
/// Corresponds to [`scipy.special.eval_gegenbauer`][docs], and `xsf::eval_gegenbauer`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_gegenbauer.html
///
/// # Definition
///
/// The Gegenbauer polynomials can be defined via the Gauss hypergeometric function $_2F_1$ as
///
/// $$
/// C_n^{(\alpha)}(z) =
///     {n+2\alpha-1 \choose n}\\,
///     \hyp{2}{1}{-n,\enspace n+2\alpha}{{1 \over 2}+\alpha}{-{z-1 \over 2}}
/// $$
///
/// When $n$ is an integer the result is a polynomial of degree $n$.
/// See Abramowitz & Stegun 22.5.46 [^AS] or DLMF 18.5.9 [^DLMF] for details.
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
///   Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E9>
///
/// # See also
/// - [`eval_jacobi`]: Evaluate Jacobi polynomials, $P_n^{(\alpha, \beta)}$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $_2F_1$
#[doc(alias = "eval_ultraspherical")]
#[inline]
pub fn eval_gegenbauer<N, Z: GegenbauerArg<N>>(n: N, alpha: f64, z: Z) -> Z {
    z.eval_gegenbauer(n, alpha)
}

///////////////////////////////
// Chebyshev T (first kind)
// Chebyshev U (second kind)
//

pub trait ChebyshevArg<N>: sealed::OrthoPolyArg {
    fn eval_chebyshev_t(self, n: N) -> Self;
    fn eval_chebyshev_u(self, n: N) -> Self;
}

impl<Z: sealed::OrthoPolyArg> ChebyshevArg<f64> for Z {
    #[inline]
    fn eval_chebyshev_t(self, n: f64) -> Self {
        self.hyp2f1_sh(n, n, 0.5, 0.0)
    }

    #[inline]
    fn eval_chebyshev_u(self, n: f64) -> Self {
        self.hyp2f1_sh(n, n + 2.0, 1.5, 2.0) * (n + 1.0)
    }
}

impl ChebyshevArg<i32> for f64 {
    /// Corresponds to `xsf::detail::eval_chebyt_l`
    // Use the Chebyshev T recurrence directly, see Mason & Handscomb, Chebyshev Polynomials, CRC Press (2003).
    #[inline]
    fn eval_chebyshev_t(self, k: i32) -> Self {
        if self.is_nan() {
            return self;
        }

        // symmetry: T_{-k}(x) = T_k(x)
        let k = k.unsigned_abs();

        let mut b2 = 0.0;
        let mut b1 = -1.0;
        let mut b0 = 0.0;
        let x = 2.0 * self;
        for _ in 0..=k {
            b2 = b1;
            b1 = b0;
            b0 = x * b1 - b2;
        }
        (b0 - b2) / 2.0
    }

    #[inline]
    fn eval_chebyshev_u(self, n: i32) -> Self {
        if self.is_nan() {
            self
        } else if n < -1 {
            -chebyshev_recurrence(self, 0.0, -1.0, -1..-(n + 2))
        } else {
            chebyshev_recurrence(self, 0.0, -1.0, -1..n)
        }
    }
}

#[inline]
fn chebyshev_recurrence(x: f64, p1: f64, p2: f64, range: Range<i32>) -> f64 {
    let x2 = x * 2.0;
    range.fold((p1, p2), |(p1, p2), _| (x2 * p1 - p2, p1)).0
}

/// Evaluate Chebyshev polynomial of the first kind $T_n$ at a point.
///
/// This is a pure rust implementation equivalent to [`scipy.special.eval_chebyt`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebyt.html
///
/// # Definition
///
/// The Chebyshev polynomials of the first kind can be defined via Gauss' hypergeometric function,
/// $_2F_1$, as
///
/// $$ T_n(z) = \hyp{2}{1}{-n,\enspace n}{1 \over 2}{-{z-1 \over 2}} $$
///
/// When $n$ is an integer the result is an orthogonal polynomial of degree $\abs n$.
/// See Abramowitz & Stegun (1972) eq. 22.5.47 [^AS] for details.
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
///   Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
///
/// # See also
/// - [`eval_chebyshev_u`]: Evaluate Chebyshev polynomials of the second kind, $U_n$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $_2F_1$
#[doc(alias = "eval_chebyt")]
#[inline]
pub fn eval_chebyshev_t<N, Z: ChebyshevArg<N>>(n: N, z: Z) -> Z {
    z.eval_chebyshev_t(n)
}

/// Evaluate Chebyshev polynomial of the second kind $U_n$ at a point.
///
/// This is a pure rust implementation equivalent to [`scipy.special.eval_chebyu`][docs].
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebyu.html
///
/// # Definition
///
/// The Chebyshev polynomials of the second kind can be defined via Gauss' hypergeometric function,
/// $_2F_1$, as
///
/// $$ U_n(z) = (n+1)\\ \hyp{2}{1}{-n,\enspace 2+n}{3 \over 2}{-{z-1 \over 2}} $$
///
/// When $n$ is an integer the result is an orthogonal polynomial of degree $\abs n$.
/// See Abramowitz & Stegun (1972) eq. 22.5.48 [^AS] for details.
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
///   Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
///
/// # See also
/// - [`eval_chebyshev_t`]: Evaluate Chebyshev polynomials of the first kind, $T_n$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $_2F_1$
#[doc(alias = "eval_chebyu")]
#[inline]
pub fn eval_chebyshev_u<N, Z: ChebyshevArg<N>>(n: N, z: Z) -> Z {
    z.eval_chebyshev_u(n)
}

/// Evaluate Chebyshev polynomial of the first kind $C_n$ on $[-2, 2]$ at a point.
///
/// Corresponds to [`scipy.special.eval_chebyc`][docs], and `xsf::eval_chebyc`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebyc.html
///
/// # Definition
///
/// $$ C_n(z) = 2 T_n(z / 2) $$
///
/// # See also
/// - [`eval_chebyshev_t`]: Evaluate Chebyshev polynomials of the first kind, $T_n$
#[doc(alias = "eval_chebyc")]
#[inline]
pub fn eval_chebyshev_c<N, Z: ChebyshevArg<N>>(n: N, z: Z) -> Z {
    (z * 0.5).eval_chebyshev_t(n) * 2.0
}

/// Evaluate Chebyshev polynomial of the second kind $S_n$ on $[-2, 2]$ at a point.
///
/// Corresponds to [`scipy.special.eval_chebys`][docs], and `xsf::eval_chebys`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebys.html
///
/// # Definition
///
/// $$ S_n(z) = U_n(z / 2) $$
///
/// # See also
/// - [`eval_chebyshev_u`]: Evaluate Chebyshev polynomials of the second kind, $U_n$
#[doc(alias = "eval_chebys")]
#[inline]
pub fn eval_chebyshev_s<N, Z: ChebyshevArg<N>>(n: N, z: Z) -> Z {
    (z * 0.5).eval_chebyshev_u(n)
}

/// Evaluate shifted Chebyshev polynomial of the first kind $T_n^*$ at a point.
///
/// Corresponds to [`scipy.special.eval_sh_chebyt`][docs], and `xsf::eval_sh_chebyt`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_chebyt.html
///
/// # Definition
///
/// $$ T_n^*(z) = T_n(2z - 1) $$
///
/// # See also
/// - [`eval_chebyshev_t`]: Evaluate Chebyshev polynomials of the first kind, $T_n$
#[doc(alias = "eval_sh_chebyt")]
#[inline]
pub fn eval_sh_chebyshev_t<N, Z: ChebyshevArg<N>>(n: N, z: Z) -> Z {
    (z * 2.0 - 1.0).eval_chebyshev_t(n)
}

/// Evaluate shifted Chebyshev polynomial of the second kind $U_n^*$ at a point.
///
/// Corresponds to [`scipy.special.eval_sh_chebyu`][docs], and `xsf::eval_sh_chebyu`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_chebyu.html
///
/// # Definition
///
/// $$ U_n^*(z) = U_n(2z - 1) $$
///
/// # See also
/// - [`eval_chebyshev_u`]: Evaluate Chebyshev polynomials of the second kind, $U_n$
#[doc(alias = "eval_sh_chebyu")]
#[inline]
pub fn eval_sh_chebyshev_u<N, Z: ChebyshevArg<N>>(n: N, z: Z) -> Z {
    (z * 2.0 - 1.0).eval_chebyshev_u(n)
}

///////////////////////////////
// Generalized Laguerre
// Laguerre
//

pub trait LaguerreArg<N>: sealed::OrthoPolyArg {
    fn eval_genlaguerre(self, n: N, a: f64) -> Self;
}

/// `binom(n + alpha, n)`, for `eval_jacobi`
///
/// Unlike in xsf (scipy/xsf#457), this is `binom(n + alpha, n - 1) * b / n` with `b = alpha + 1` if
/// `alpha < 0` and `1 <= n <= 2^53` (so `n - 1` is exact), since rounding `n + alpha` loses the
/// digits of a small `b`. `binom_alpha` is more accurate, but near `alpha = -2, -3, ...` it exposes
/// the cancellation of the `eval_jacobi` recurrence for `beta = +-1`.
fn jacobi_binom(n: f64, alpha: f64) -> f64 {
    if alpha < 0.0 && EXACT_N_MINUS_1.contains(&n) {
        crate::binom(n + alpha, n - 1.0) * (alpha + 1.0) / n
    } else {
        crate::binom(n + alpha, n)
    }
}

// binom(n + alpha, n) = 1 / (n B(alpha + 1, n)) (DLMF 5.12.1), keeping the digits of alpha + 1 for
// alpha near -1
fn binom_alpha(n: f64, alpha: f64) -> f64 {
    if alpha < -0.5 && n >= 1.0 {
        return 1.0 / n / crate::beta(alpha + 1.0, n);
    }
    crate::binom(n + alpha, n)
}

impl<Z: sealed::OrthoPolyArg> LaguerreArg<f64> for Z {
    /// Corresponds to `xsf::detail::eval_genlaguerre` in `xsf/cpu/orthogonal_eval.h`
    #[inline]
    fn eval_genlaguerre(self, n: f64, alpha: f64) -> Z {
        if alpha <= -1.0 {
            // polynomial defined only for alpha > -1
            return Self::ZERO + f64::NAN;
        }

        let d = binom_alpha(n, alpha);
        self.hyp1f1(-n, alpha + 1.0) * d
    }
}

impl LaguerreArg<i32> for f64 {
    /// Corresponds to `xsf::detail::eval_genlaguerre_l`
    // Integer-degree L_n^(alpha)(x) by forward recurrence on the normalized values
    // p_k = L_k^(alpha)(x) / binom(k + alpha, k) and their differences d_k = p_k - p_{k-1}:
    //   d_{k+1} = (k d_k - x p_k) / (k + alpha + 1),  p_{k+1} = p_k + d_{k+1},
    // starting from p_0 = 1, d_0 = 0.
    #[inline]
    fn eval_genlaguerre(self, n: i32, alpha: f64) -> Self {
        let x = self;
        if alpha <= -1.0 || alpha.is_nan() || x.is_nan() {
            // polynomial defined only for alpha > -1
            return f64::NAN;
        }

        if n < 0 {
            return 0.0;
        } else if n == 0 {
            return 1.0;
        } else if n == 1 {
            return -x + alpha + 1.0;
        }

        let mut d = -x / (alpha + 1.0);
        let mut p = d + 1.0;
        for kk in 0..n - 1 {
            let k = f64::from(kk) + 1.0;
            d = -x / (k + alpha + 1.0) * p + (k / (k + alpha + 1.0)) * d;
            p += d;
        }
        binom_alpha(f64::from(n), alpha) * p
    }
}

/// Evaluate generalized Laguerre polynomial $L_n^{(\alpha)}$ at a point.
///
/// Corresponds to [`scipy.special.eval_genlaguerre`][glag], and `xsf::eval_genlaguerre`.
///
/// [glag]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_genlaguerre.html
///
/// # Definition
///
/// The generalized Laguerre polynomials can be defined via the confluent hypergeometric function
/// $_1F_1$ as
///
/// $$ L_n^{(\alpha)}(z) = {n+\alpha \choose n} \\, \hyp 1 1 {-n} {1+\alpha} z $$
///
/// When $n$ is an integer the result is a polynomial of degree $n$ [^AS] [^DLMF].
/// The Laguerre polynomials, $L_n = L_n^{(0)}$, are the special case where $\alpha = 0$.
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
/// Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E12>
///
/// # See also
/// - [`eval_laguerre`]: Evaluate Laguerre polynomials, $L_n = L_n^{(0)}$
/// - [`hyp1f1`](crate::hyp1f1): Confluent hypergeometric function, $_1F_1$
#[inline]
pub fn eval_genlaguerre<N, Z: LaguerreArg<N>>(n: N, alpha: f64, z: Z) -> Z {
    z.eval_genlaguerre(n, alpha)
}

/// Evaluate Laguerre polynomial $L_n$ at a point
///
/// Corresponds to [`scipy.special.eval_laguerre`][lag], and `xsf::eval_laguerre`.
///
/// [lag]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_laguerre.html
///
/// # Definition
///
/// The Laguerre polynomials can be defined via the confluent hypergeometric function $_1F_1$ as
///
/// $$ L_n(z) = \hyp 1 1 {-n} 1 z $$
///
/// See 22.5.16 and 22.5.54 in [^AS] for details.
///
/// When $n$ is an integer the result is a polynomial of degree $n$ [^DLMF].
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
/// Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E12>
///
/// # See also
/// - [`eval_genlaguerre`]: Evaluate generalized Laguerre polynomials, $L_n^{(\alpha)}$
/// - [`hyp1f1`](crate::hyp1f1): Confluent hypergeometric function, $_1F_1$
#[inline]
pub fn eval_laguerre<N, Z: LaguerreArg<N>>(n: N, z: Z) -> Z {
    z.eval_genlaguerre(n, 0.0)
}

///////////////////////////////
// Hermite H  (physicists')
// Hermite He (probabilists')
//

/// `H_n(x)` (`scale = 1`) or `He_n(x)` (`scale = 1/2`), for `x > 0` and non-integer `n > 0`:
///
/// - start at `v = n - ceil(n) - 3` and `v + 1`, with the integral of DLMF 12.5.1
///   `y_v(x) = 1/Gamma(-v) int_0^inf t^(-v-1) exp(-scale t^2 - 2 scale x t) dt`
/// - evaluate it by the trapezoidal rule in `u = ln(t / t0)`, with `t0` the peak of the integrand
/// - recur forward to `n` with `y_{v+1} = 2 scale (x y_v - v y_{v-1})` (DLMF 12.8.1)
fn hermite_rec(x: f64, n: f64, scale: f64) -> f64 {
    const H: f64 = 0.1;
    // the terms for u < 0 only shrink by exp(-k H) < 0.75 per step, so their tail is < 4x the last
    const TOL: f64 = f64::EPSILON / 16.0;

    let m = n.ceil() + 3.0;
    let k = m - n;
    let sx = scale * x;
    let t0 = k / (sx + (sx * sx + 2.0 * scale * k).sqrt());
    let (a, b) = (scale * t0 * t0, 2.0 * sx * t0);

    // the integrands for v and v + 1, relative to those at t0
    let integrands = |u: f64| {
        let e = u.exp();
        let g = (k * u - a * (e * e - 1.0) - b * (e - 1.0)).exp();
        (g, g / e)
    };

    // their sums over u = 0, H, 2H, ... and u = -H, -2H, ..., each until the terms are negligible
    let mut s = (0.0, 0.0);
    for (i0, du) in [(0, H), (1, -H)] {
        for i in i0.. {
            let (g0, g1) = integrands(f64::from(i) * du);
            s = (s.0 + g0, s.1 + g1);
            if g0 < TOL * s.0 && g1 < TOL * s.1 {
                break;
            }
        }
    }

    // y_v and y_{v+1} for v = -k, and then the recurrence up to v = n
    let c = H * t0.powf(k) * (-a - b).exp() / crate::gamma(k);
    let y = (c * s.0, c * s.1 * (k - 1.0) / t0);
    let (_, yn) = (1..m as u32).rev().fold(y, |(y0, y1), j| {
        let v = n - f64::from(j);
        (y1, 2.0 * scale * (x * y1 - v * y0))
    });
    yn
}

#[inline]
fn eval_hermite_impl(x: f64, n: f64, scale: f64) -> f64 {
    if x.is_nan() || n.is_nan() {
        f64::NAN
    } else if n == 0.0 {
        1.0
    } else if x.is_zero() {
        if n < 0.0 || n % 2.0 == 0.0 {
            SQRT_PI * (scale * n).exp2() / crate::gamma(0.5 - 0.5 * n)
        } else {
            0.0
        }
    } else if x.is_sign_negative() && n.is_finite() && n.fract() != 0.0 {
        // For non-integer n and x < 0, the negative-base `powf` of the `hypu` formula below is NaN,
        // so use the two Kummer M terms of U(-n/2, 1/2, x^2) (DLMF 13.2.42), with |x| -> x:
        //
        // H_n(x) = 2^n sqrt(pi) [M(-n/2, 1/2, x^2) / Gamma((1-n)/2)
        //                        - 2x M((1-n)/2, 3/2, x^2) / Gamma(-n/2)]
        //
        // For x > 0 these terms cancel (H_n(x) ~ (2x)^n, while both terms grow like e^(x^2)),
        // which is why `hypu` is used there.
        //
        // E.g. `eval_hermite_h(2.5, -5.0)` was NaN, not -5.77e8.
        let z = scale * x * x;
        let m0 = crate::hyp1f1(-0.5 * n, 0.5, z) * crate::rgamma(0.5 - 0.5 * n);
        let m1 = crate::hyp1f1(0.5 - 0.5 * n, 1.5, z) * crate::rgamma(-0.5 * n);
        SQRT_PI * (scale * n).exp2() * (m0 - 2.0 * scale.sqrt() * x * m1)
    } else {
        // For non-integer n > 0 and x > 0 (which xsf and SciPy don't support), use `hermite_rec`,
        // since `hypu` below has relative errors of up to 1e-6 there.
        //
        // E.g. `eval_hermite_h(55.5, 3.5)` was 2.862819e47 instead of 2.862818e47.
        if (0.0..1e4).contains(&x) && (0.0..1e3).contains(&n) && n.fract() != 0.0 {
            let h = hermite_rec(x, n, scale);
            if h.is_finite() {
                return h;
            }
        }
        let c = if x.is_sign_positive() {
            (scale * n).exp2()
        } else {
            (-2.0 * scale.sqrt()).powf(n)
        };
        c * crate::hypu(-0.5 * n, 0.5, scale * x * x)
    }
}

pub trait HermiteArg<N>: sealed::OrthoPolyArg {
    fn eval_hermite_he(self, n: N) -> Self;
    fn eval_hermite_h(self, n: N) -> Self;
}

impl HermiteArg<f64> for f64 {
    #[inline]
    fn eval_hermite_he(self, n: f64) -> f64 {
        eval_hermite_impl(self, n, 0.5)
    }

    #[inline]
    fn eval_hermite_h(self, n: f64) -> f64 {
        eval_hermite_impl(self, n, 1.0)
    }
}

impl HermiteArg<u32> for f64 {
    /// Corresponds to `xsf::detail::eval_hermitenorm`
    #[inline]
    fn eval_hermite_he(self, n: u32) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            let (mut y1, mut y2) = (1.0, 0.0);
            for k in (1..=n).rev() {
                (y1, y2) = (self * y1 - f64::from(k) * y2, y1);
            }
            y1
        }
    }

    /// Corresponds to `xsf::detail::eval_hermite`
    #[inline]
    fn eval_hermite_h(self, n: u32) -> f64 {
        (f64::from(n) / 2.0).exp2() * (SQRT_2 * self).eval_hermite_he(n)
    }
}

/// Evaluate probabilists' (normalized) Hermite polynomial $He_n(x)$ at a point
///
/// Corresponds to [`scipy.special.eval_hermitenorm`][docs], and `xsf::eval_hermitenorm`, but
/// also supports non-integer $n$.
///
/// `n` can be either `u32` or `f64`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermitenorm.html
///
/// # Definition
///
/// For non-negative integer $n$,
///
/// $$
/// He_n(x) =
///     x^n \\
///     {}_2F_0\left[-{n \over 2}, -{n-1 \over 2} \middle| -{2 \over x^2} \right]
/// $$
///
/// See 22.11.8 in [^AS] for details.
///
/// When `n: u32` the result is a polynomial of degree $n$ [^DLMF].
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
/// Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E13>
///
/// # See also
/// - [`eval_hermite_h`]: Evaluate physicists' Hermite polynomials, $H_n$
///
#[doc(alias = "eval_hermitenorm")]
#[inline]
pub fn eval_hermite_he<N>(n: N, x: f64) -> f64
where
    f64: HermiteArg<N>,
{
    x.eval_hermite_he(n)
}

/// Evaluate physicists' Hermite polynomial $H_n(x)$ at a point
///
/// Corresponds to [`scipy.special.eval_hermite`][docs], and `xsf::eval_hermite`, but also supports
/// non-integer $n$.
///
/// `n` can be either `u32` or `f64`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermite.html
///
/// # Definition
///
/// For non-negative integer $n$,
///
/// $$
/// H_n(x) =
///     (2x)^n \\
///     {}_2F_0\left[-{n \over 2}, -{n-1 \over 2} \middle| -{1 \over x^2} \right]
/// $$
///
/// See 22.11.8 in [^AS] for details.
///
/// When `n: u32` the result is a polynomial of degree $n$ [^DLMF].
///
/// [^AS]: Milton Abramowitz and Irene A. Stegun, eds. Handbook of Mathematical Functions with
/// Formulas, Graphs, and Mathematical Tables. New York: Dover, 1972.
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.5.E13>
///
/// # See also
/// - [`eval_hermite_he`]: Evaluate probabilists' (normalized) Hermite polynomials, $He_n$
///
#[doc(alias = "eval_hermite")]
#[inline]
pub fn eval_hermite_h<N>(n: N, x: f64) -> f64
where
    f64: HermiteArg<N>,
{
    x.eval_hermite_h(n)
}
