//! Translated from `xsf/orthogonal_eval.h` and `xsf/cpu/orthogonal_eval.h` in xsf main, which in
//! turn were translated from `scipy/special/orthogonal_eval.pxd`
//!
//! - <https://github.com/scipy/xsf/blob/b6cb425/include/xsf/orthogonal_eval.h>
//! - <https://github.com/scipy/xsf/blob/b6cb425/include/xsf/cpu/orthogonal_eval.h>

use core::f64::consts::SQRT_2;
use core::ops::Range;
use num_traits::Zero;

use crate::utils::{MAXGAM, SQRT_PI};

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

        #[inline]
        fn hyp2f1_sh(self, a1: f64, a2: f64, b: f64) -> Self {
            (self.clone() * -0.5 + 0.5).hyp2f1(a1, a2, b)
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
    #[allow(clippy::float_cmp, clippy::manual_midpoint)]
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
        let d = crate::binom(n + a, n);
        self.hyp2f1_sh(-n, n + a + b + 1.0, a + 1.0) * d
    }

    /// Corresponds to `xsf::detail::eval_legendre`
    #[inline]
    fn eval_legendre(self, n: f64) -> Self {
        if self.is_nan() {
            return self;
        }
        self.hyp2f1_sh(-n, n + 1.0, 1.0)
    }
}

impl JacobiArg<i32> for f64 {
    /// Corresponds to `xsf::detail::eval_jacobi_l`
    #[inline]
    #[allow(clippy::float_cmp, clippy::manual_midpoint)]
    fn eval_jacobi(self, n: i32, alpha: f64, beta: f64) -> Self {
        let x = self;
        if n < 0 {
            x.eval_jacobi(f64::from(n), alpha, beta)
        } else if n == 0 {
            1.0
        } else if n == 1 {
            0.5 * (2.0 * (alpha + 1.0) + (alpha + beta + 2.0) * (x - 1.0))
        } else if alpha == -1.0 && beta.abs() == 1.0 {
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
            let n = f64::from(n);
            crate::binom(n + alpha, n) * p
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
/// This is a translation of the [`scipy.special.eval_jacobi`][jaceval] Cython implementation
/// into Rust.
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
/// This is a translation of the [`scipy.special.eval_legendre`][legeval] Cython implementation
/// into Rust.
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
/// Corresponds to [`scipy.special.eval_sh_jacobi`][docs] (and `xsf::eval_sh_jacobi` in xsf
/// main), which, despite its name, isn't the shifted Jacobi polynomial [`eval_sh_jacobi`], see
/// scipy/scipy#18988.
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
/// Corresponds to [`scipy.special.eval_sh_legendre`][docs], and `xsf::eval_sh_legendre` in xsf
/// main.
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
    #[allow(clippy::cast_possible_truncation)]
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
            // function overflows, e.g. `eval_gegenbauer(0.5, 100.0, 1.0)` was NaN instead of 15.95.
            d = crate::binom(n + 2.0 * alpha - 1.0, n);
        }
        self.hyp2f1_sh(-n, n + 2.0 * alpha, alpha + 0.5) * d
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
#[allow(clippy::cast_precision_loss, clippy::float_cmp)]
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
/// This is a translation of the [`scipy.special.eval_gegenbauer`][docs] Cython implementation
/// into Rust.
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
/// See Abramowitz & Stegun 22.5.46 [^AS] or DLMF 18.5.7 [^DLMF] for details.
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
        self.hyp2f1_sh(-n, n, 0.5)
    }

    #[inline]
    fn eval_chebyshev_u(self, n: f64) -> Self {
        self.hyp2f1_sh(-n, n + 2.0, 1.5) * (n + 1.0)
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
            -chebyshev_recurrence(self, 0.0, -1.0, -1..(-n - 2))
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
/// Corresponds to [`scipy.special.eval_chebyc`][docs], and `xsf::eval_chebyc` in xsf main.
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
/// Corresponds to [`scipy.special.eval_chebys`][docs], and `xsf::eval_chebys` in xsf main.
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
/// Corresponds to [`scipy.special.eval_sh_chebyt`][docs], and `xsf::eval_sh_chebyt` in xsf main.
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
/// Corresponds to [`scipy.special.eval_sh_chebyu`][docs], and `xsf::eval_sh_chebyu` in xsf main.
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

impl<Z: sealed::OrthoPolyArg> LaguerreArg<f64> for Z {
    /// Corresponds to `xsf::detail::eval_genlaguerre` in `xsf/cpu/orthogonal_eval.h`
    #[inline]
    fn eval_genlaguerre(self, n: f64, alpha: f64) -> Z {
        if alpha <= -1.0 {
            // polynomial defined only for alpha > -1
            return Self::ZERO + f64::NAN;
        }

        let d = crate::binom(n + alpha, n);
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
        let n = f64::from(n);
        crate::binom(n + alpha, n) * p
    }
}

/// Evaluate generalized Laguerre polynomial $L_n^{(\alpha)}$ at a point.
///
/// This is a translation of the [`scipy.special.eval_genlaguerre`][glag] Cython implementation
/// into Rust.
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
/// This is a translation of the [`scipy.special.eval_laguerre`][lag] Cython implementation
/// into Rust.
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
        // E.g. `eval_hermite_h(2.5298681626028885, -5.382735923248713)` was NaN, not -2.26e10.
        let z = scale * x * x;
        let m0 = crate::hyp1f1(-0.5 * n, 0.5, z) * crate::rgamma(0.5 - 0.5 * n);
        let m1 = crate::hyp1f1(0.5 - 0.5 * n, 1.5, z) * crate::rgamma(-0.5 * n);
        SQRT_PI * (scale * n).exp2() * (m0 - 2.0 * scale.sqrt() * x * m1)
    } else {
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
    /// Corresponds to `eval_hermitenorm` in `scipy/special/orthogonal_eval.pxd`
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

    /// Corresponds to `eval_hermite` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_hermite_h(self, n: u32) -> f64 {
        (f64::from(n) / 2.0).exp2() * (SQRT_2 * self).eval_hermite_he(n)
    }
}

/// Evaluate probabilists' (normalized) Hermite polynomial $He_n(x)$ at a point
///
/// This is a translation of the [`scipy.special.eval_hermitenorm`][docs] Cython implementation
/// into Rust, with additional support for non-integer $n$.
///
/// `n` can be either `u32` or `f64`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermitenorm.html
///
/// # Definition
///
/// $$
/// He_n(x) =
///     2^{n \over 2} x^n \\
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
/// This is a translation of the [`scipy.special.eval_hermite`][docs] Cython implementation
/// into Rust, with additional support for non-integer $n$.
///
/// `n` can be either `u32` or `f64`.
///
/// [docs]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermite.html
///
/// # Definition
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

///////////////////////////////
// Tests
//

#[cfg(test)]
mod tests {
    use crate::{
        eval_chebyshev_t, eval_chebyshev_u, eval_gegenbauer, eval_genlaguerre, eval_hermite_h,
        eval_hermite_he, eval_jacobi, eval_laguerre, eval_legendre, np_assert_allclose,
    };
    use num_complex::c64;
    use num_traits::ToPrimitive;

    #[test]
    fn test_eval_jacobi_eq_legendre() {
        // manual test cases
        let xs = [-5.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 5.0];
        let ns = [0, 1, 2, 3, 4, 8, 15];
        for n in ns {
            let p0 = xs.map(|x| crate::legendre_p(n, x));

            let pi = xs.map(|x| eval_jacobi(n, 0.0, 0.0, x));
            np_assert_allclose!(pi, p0, rtol = 1e-13);

            let pf = xs.map(|x| eval_jacobi(f64::from(n), 0.0, 0.0, x));
            np_assert_allclose!(pf, p0, rtol = 1e-12, atol = f64::EPSILON);

            let pc = xs.map(|x| eval_jacobi(f64::from(n), 0.0, 0.0, c64(x, 0.0)).re);
            // the xsf hyp2f1 implementation for complex numbers isn't very accurate
            np_assert_allclose!(pc, p0, rtol = 1e-7, atol = 1e-8);
        }
    }

    #[test]
    fn test_eval_jacobi_scipy() {
        // test data obtained from `scipy.special.eval_jacobi`
        let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
        let ns = [0, 1, 2, 3];

        let expect_1_0 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-7.0, -1.0, 0.35, 0.5, 0.65, 2.0, 8.0],
            [57.0, 1.0, -0.575, -0.5, -0.375, 3.0, 67.0],
            [-491.0, -1.0, -0.173_125, -0.375, -0.539_375, 4.0, 584.0],
        ];
        let expect_0_1 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-8.0, -2.0, -0.65, -0.5, -0.35, 1.0, 7.0],
            [67.0, 3.0, -0.375, -0.5, -0.575, 1.0, 57.0],
            [-584.0, -4.0, 0.539_375, 0.375, 0.173_125, 1.0, 491.0],
        ];
        let expect_1_1 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-10.0, -2.0, -0.2, 0.0, 0.2, 2.0, 10.0],
            [93.0, 3.0, -0.7125, -0.75, -0.7125, 3.0, 93.0],
            [-860.0, -4.0, 0.293, 0.0, -0.293, 4.0, 860.0],
        ];
        let expect_h_h = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-7.5, -1.5, -0.15, 0.0, 0.15, 1.5, 7.5],
            [61.875, 1.875, -0.6, -0.625, -0.6, 1.875, 61.875],
            [
                -535.9375, -2.1875, 0.214_375, 0.0, -0.214_375, 2.1875, 535.9375,
            ],
        ];
        let expect_h_nh = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-4.5, -0.5, 0.4, 0.5, 0.6, 1.5, 5.5],
            [33.375, 0.375, -0.435, -0.375, -0.285, 1.875, 40.875],
            [
                -275.3125, -0.3125, -0.1775, -0.3125, -0.4225, 2.1875, 337.1875,
            ],
        ];
        let expect_nh_nh = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-2.5, -0.5, -0.05, 0.0, 0.05, 0.5, 2.5],
            [18.375, 0.375, -0.3675, -0.375, -0.3675, 0.375, 18.375],
            [-151.5625, -0.3125, 0.0925, 0.0, -0.0925, 0.3125, 151.5625],
        ];
        let expect = [
            (1.0, 0.0, expect_1_0),
            (0.0, 1.0, expect_0_1),
            (1.0, 1.0, expect_1_1),
            (0.5, 0.5, expect_h_h),
            (0.5, -0.5, expect_h_nh),
            (-0.5, -0.5, expect_nh_nh),
        ];

        for (alpha, beta, expected) in expect {
            for (i, &n) in ns.iter().enumerate() {
                let p_i32 = xs.map(|x| eval_jacobi(n, alpha, beta, x));
                np_assert_allclose!(p_i32, expected[i], rtol = 1e-14, atol = f64::EPSILON);

                let p_f64 = xs.map(|x| eval_jacobi(f64::from(n), alpha, beta, x));
                np_assert_allclose!(p_f64, expected[i], rtol = 1e-14, atol = f64::EPSILON);
            }
        }
    }

    #[test]
    fn test_eval_legendre_eq_legendre_f64() {
        // manual test cases
        let xs = [-5.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 5.0];
        let ns = [0, 1, 2, 3, 4, 8, 15];
        for n in ns {
            let p0 = xs.map(|x| crate::legendre_p(n, x));

            let pi = xs.map(|x| eval_legendre(n, x));
            np_assert_allclose!(pi, p0, rtol = 1e-13);

            let pf = xs.map(|x| eval_legendre(f64::from(n), x));
            np_assert_allclose!(pf, p0, rtol = 1e-12, atol = f64::EPSILON);
        }
    }

    #[test]
    fn test_eval_gegenbauer_eq_self() {
        // manual test cases
        let xs = [
            0.0,
            1.0,
            -1.0,
            1e-6,
            -1e-6,
            2.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        let zs = [c64(0.5, 0.0), c64(-0.5, 0.0), c64(1.0, 0.0), c64(-1.0, 0.0)];
        let ns = [0, 1, 2, 3, 4, 11];
        let alphas = [0.0, 0.5, 1.0, 2.0];
        for &n in &ns {
            for &alpha in &alphas {
                let c_i32 = xs.map(|x| eval_gegenbauer(n, alpha, x));
                let c_f64 = xs.map(|x| eval_gegenbauer(f64::from(n), alpha, x));
                np_assert_allclose!(c_i32, c_f64, rtol = 1.5e-14, atol = 1.5e-14);

                let c_f64 = zs.map(|z| eval_gegenbauer(f64::from(n), alpha, z.re));
                let c_c64 = zs.map(|z| eval_gegenbauer(f64::from(n), alpha, z).re);
                np_assert_allclose!(c_c64, c_f64, rtol = 1e-8, atol = 1e-12);
            }
        }
    }

    #[test]
    fn test_eval_genlaguerre_scipy() {
        // test data obtained from `scipy.special.eval_genlaguerre`
        let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
        let ns = [0, 1, 2, 3];

        let expect_0 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [6.0, 2.0, 1.1, 1.0, 0.9, 0.0, -4.0],
            [23.5, 3.5, 1.205, 1.0, 0.805, -0.5, 3.5],
            [
                74.333_333_333_333_33,
                5.666_666_666_666_667,
                1.315_166_666_666_666_7,
                1.0,
                0.714_833_333_333_333_3,
                -0.666_666_666_666_666_7,
                2.666_666_666_666_666_7,
            ],
        ];
        let expect_1 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [7.0, 3.0, 2.1, 2.0, 1.9, 1.0, -3.0],
            [30.5, 6.5, 3.305, 3.0, 2.705, 0.5, 0.5],
            [
                104.833_333_333_333_33,
                12.166_666_666_666_667,
                4.620_166_666_666_667,
                4.0,
                3.419_833_333_333_333_3,
                -0.166_666_666_666_666_7,
                3.166_666_666_666_666_7,
            ],
        ];
        let expect_h = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [6.5, 2.5, 1.6, 1.5, 1.4, 0.5, -3.5],
            [26.875, 4.875, 2.13, 1.875, 1.63, -0.125, 1.875],
            [
                88.645_833_333_333_33,
                8.479_166_666_666_667,
                2.642_666_666_666_666_7,
                2.1875,
                1.767_333_333_333_333_3,
                -0.604_166_666_666_666_7,
                3.229_166_666_666_666_7,
            ],
        ];
        let expect_nh = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [5.5, 1.5, 0.6, 0.5, 0.4, -0.5, -4.5],
            [20.375, 2.375, 0.53, 0.375, 0.23, -0.625, 5.375],
            [
                61.770_833_333_333_333,
                3.604_166_666_666_666_7,
                0.512_666_666_666_666_7,
                0.3125,
                0.137_333_333_333_333_33,
                -0.479_166_666_666_666_7,
                1.354_166_666_666_666_7,
            ],
        ];
        let expect = [
            (0.0, expect_0),
            (1.0, expect_1),
            (0.5, expect_h),
            (-0.5, expect_nh),
        ];

        for (alpha, expected) in expect {
            for (i, &n) in ns.iter().enumerate() {
                let p_i32 = xs.map(|x| eval_genlaguerre(n, alpha, x));
                np_assert_allclose!(p_i32, expected[i], rtol = 1e-15, atol = f64::EPSILON);

                let p_f64 = xs.map(|x| eval_genlaguerre(f64::from(n), alpha, x));
                np_assert_allclose!(p_f64, expected[i], rtol = 1e-14, atol = f64::EPSILON);
            }
        }
    }

    #[test]
    fn test_eval_laguerre_scipy() {
        // test data obtained from `scipy.special.eval_laguerre`
        let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
        let ns = [0, 1, 2, 3, 4, 5, 6, 7, 8];

        let expect = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [6.0, 2.0, 1.1, 1.0, 0.9, 0.0, -4.0],
            [23.5, 3.5, 1.205, 1.0, 0.805, -0.5, 3.5],
            [
                74.333_333_333_333_33,
                5.666_666_666_666_667,
                1.315_166_666_666_666_7,
                1.0,
                0.714_833_333_333_333_3,
                -0.666_666_666_666_666_7,
                2.666_666_666_666_666_7,
            ],
            [
                205.375,
                8.708_333_333_333_333,
                1.430_670_833_333_333_3,
                1.0,
                0.629_337_5,
                -0.625,
                -1.291_666_666_666_666_7,
            ],
            [
                515.583_333_333_333_3,
                12.883_333_333_333_333,
                1.551_687_583_333_333_3,
                1.0,
                0.548_354_083_333_333_3,
                -0.466_666_666_666_666_7,
                -3.166_666_666_666_666_7,
            ],
            [
                1_203.743_055_555_555_6,
                18.509_722_222_222_222,
                1.678_396_334_722_222_2,
                1.0,
                0.471_728_668_055_555_6,
                -0.256_944_444_444_444_4,
                -2.090_277_777_777_777_8,
            ],
            [
                2_653.410_714_285_714_7,
                25.976_587_301_587_3,
                1.810_980_926_408_730_7,
                1.0,
                0.399_310_759_702_381_06,
                -0.040_476_190_476_190_44,
                0.325_396_825_396_825_6,
            ],
            [
                5_580.251_612_103_175,
                35.757_167_658_730_154,
                1.949_629_705_714_534_2,
                1.0,
                0.330_953_705_397_073_47,
                0.153_993_055_555_555_56,
                2.235_739_087_301_587,
            ],
        ];

        for (i, &n) in ns.iter().enumerate() {
            let p_i32 = xs.map(|x| eval_laguerre(n, x));
            np_assert_allclose!(p_i32, expect[i], rtol = 1e-15, atol = f64::EPSILON);

            let p_f64 = xs.map(|x| eval_laguerre(f64::from(n), x));
            np_assert_allclose!(p_f64, expect[i], rtol = 2e-14, atol = f64::EPSILON);
        }
    }

    fn test_orthogonal_poly<N, F, G>(eval_int: F, eval_f64: G, expected: &[&dyn Fn(f64) -> f64])
    where
        N: num_traits::FromPrimitive,
        F: Fn(N, f64) -> f64,
        G: Fn(f64, f64) -> f64,
    {
        let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
        for (n, &expect_fn) in expected.iter().enumerate() {
            let y_expect = xs.map(expect_fn);
            let y_actual_int = xs.map(|x| eval_int(N::from_usize(n).unwrap(), x));
            let n_f64 = n.to_f64().unwrap();
            let y_actual_f64 = xs.map(|x| eval_f64(n_f64, x));
            np_assert_allclose!(y_actual_int, y_expect, rtol = 1e-15, atol = f64::EPSILON);
            np_assert_allclose!(y_actual_f64, y_expect, rtol = 1e-15, atol = f64::EPSILON);
        }
    }

    #[test]
    fn test_eval_hermite_he() {
        test_orthogonal_poly::<u32, _, _>(
            eval_hermite_he,
            eval_hermite_he,
            &[&|_| 1.0, &|x| x, &|x| x * x - 1.0, &|x| x * x * x - 3.0 * x],
        );
    }

    #[test]
    fn test_eval_hermite_h() {
        test_orthogonal_poly::<u32, _, _>(
            eval_hermite_h,
            eval_hermite_h,
            &[&|_| 1.0, &|x| 2.0 * x, &|x| 4.0 * x * x - 2.0, &|x| {
                8.0 * x * x * x - 12.0 * x
            }],
        );
    }

    #[test]
    fn test_eval_hermite_real_degree_neg_x() {
        // mpmath; these were NaN (the tolerance is limited by the accuracy of `hyp1f1`)
        let h = [
            eval_hermite_h(2.529_868_162_602_888_5, -5.382_735_923_248_713),
            eval_hermite_h(3.3, -0.7),
            eval_hermite_h(-15.3, -8.0),
            eval_hermite_h(40.5, -3.0),
        ];
        let h_expected = [
            -22_618_949_487.496_384,
            5.292_284_105_801_645,
            9.165_934_848_962_248e29,
            -6.302_852_583_554_227e31,
        ];
        np_assert_allclose!(h, h_expected, rtol = 1e-10);

        let he = [
            eval_hermite_he(2.5, -1.5),
            eval_hermite_he(17.25, -6.0),
            eval_hermite_he(-0.5, -9.5),
        ];
        let he_expected = [
            2.234_925_546_810_356,
            110_811_784_506.491_07,
            1.824_056_114_304_247_2e19,
        ];
        np_assert_allclose!(he, he_expected, rtol = 1e-10);
    }

    #[test]
    fn test_eval_chebyshev_t() {
        test_orthogonal_poly::<i32, _, _>(
            eval_chebyshev_t,
            eval_chebyshev_t,
            &[&|_| 1.0, &|x| x, &|x| 2.0 * x * x - 1.0, &|x| {
                4.0 * x * x * x - 3.0 * x
            }],
        );
    }

    #[test]
    fn test_eval_chebyshev_u() {
        test_orthogonal_poly::<i32, _, _>(
            eval_chebyshev_u,
            eval_chebyshev_u,
            &[&|_| 1.0, &|x| 2.0 * x, &|x| 4.0 * x * x - 1.0, &|x| {
                8.0 * x * x * x - 4.0 * x
            }],
        );
    }

    #[test]
    fn test_eval_gegenbauer_negative_alpha() {
        // mpmath: gegenbauer(n, alpha, x)
        let c = crate::eval_gegenbauer(3_i32, -0.3, 0.5);
        np_assert_allclose!([c], [0.1505], rtol = 1e-14);
        let c = crate::eval_gegenbauer(4_i32, -0.25, 2.0);
        np_assert_allclose!([c], [-7.093_75], rtol = 1e-14);
        // C_1^alpha(x) = 2 alpha x, also where binom(2 alpha, 1) = 0
        assert_eq!(crate::eval_gegenbauer(1_i32, -0.5, -0.3), 0.3);
    }

    #[test]
    fn test_eval_gegenbauer_large_alpha() {
        // mpmath: gegenbauer(n, alpha, x)
        let c = [
            eval_gegenbauer(0.5, 100.0, 1.0),
            eval_gegenbauer(2.5, 90.0, -0.7),
            eval_gegenbauer(0.5, 100.0, c64(1.0, 0.0)).re,
            eval_gegenbauer(171.5, 0.01, 1.0),
            eval_gegenbauer(170.624_376_956_302_7, 0.01, 1.0),
        ];
        let expected = [
            15.947_720_785_517_17,
            -3.247_107_273_325_162e25,
            15.947_720_785_517_17,
            1.307_068_781_116_506e-4,
            1.313_641_618_259_32e-4,
        ];
        np_assert_allclose!(c, expected, rtol = 1e-12);
    }

    #[test]
    fn test_eval_jacobi_alpha_m1() {
        // mpmath: jacobi(n, -1 + 1e-35, beta, x)
        np_assert_allclose!([eval_jacobi(3_i32, -1.0, 1.0, 0.3)], [0.1925], rtol = 1e-14);
        np_assert_allclose!([eval_jacobi(3.0, -1.0, 1.0, 0.3)], [0.1925], rtol = 1e-14);
        let p = eval_jacobi(2.5, -1.0, -1.0, 0.3);
        np_assert_allclose!([p], [-0.218_362_414_580_039_7], rtol = 1e-14);
    }

    #[test]
    fn test_shifted_and_scaled() {
        // reference values from scipy.special
        let r = [
            crate::eval_gengegenbauer(3_i32, 2.5, 1.5, 0.3),
            crate::eval_gengegenbauer(2.5, 2.5, 1.5, 0.3),
            crate::eval_sh_legendre(4_i32, 0.3),
            crate::eval_sh_legendre(2.5, 0.7),
            crate::eval_chebyshev_c(5_i32, 1.3),
            crate::eval_chebyshev_s(5_i32, 1.3),
            crate::eval_chebyshev_c(2.5, 1.3),
            crate::eval_chebyshev_s(2.5, 1.3),
            crate::eval_sh_chebyshev_t(4_i32, 0.3),
            crate::eval_sh_chebyshev_u(4_i32, 0.3),
            crate::eval_sh_chebyshev_t(1.5, 0.3),
        ];
        let expected = [
            0.013_587_412_587_412_59,
            0.010_779_963_000_470_787,
            -0.113,
            -0.426_559_547_478_622_85,
            -0.772_07,
            -1.175_07,
            -1.108_120_029_599_681_6,
            0.157_988_300_295_633_32,
            -0.0752,
            -0.5104,
            -0.985_900_603_509_299_1,
        ];
        np_assert_allclose!(r, expected, rtol = 1e-13);

        let c = [
            crate::eval_sh_legendre(3.0, c64(0.2, 0.1)),
            crate::eval_sh_chebyshev_u(1.5, c64(0.3, 0.2)),
        ];
        let expected = [
            c64(0.54, 0.22),
            c64(-1.550_356_818_307_985, 0.051_799_489_871_846_57),
        ];
        np_assert_allclose!(c, expected, rtol = 1e-13);
    }

    #[test]
    fn test_eval_sh_jacobi() {
        // P_n^(alpha, beta)(2x - 1), see scipy/scipy#18988
        let xs = [0.25, 0.5, 0.75];
        let p = xs.map(|x| crate::eval_sh_jacobi(4_i32, 1.0, 2.0, x));
        np_assert_allclose!(p, [-1.210_937_5, 0.625, -0.585_937_5], rtol = 1e-14);
        let p = crate::eval_sh_jacobi(2.5, 1.0, 2.0, 0.3);
        np_assert_allclose!([p], [1.303_920_745_846_762_5], rtol = 1e-14);
    }

    #[test]
    fn test_nan_and_domain() {
        let nan = f64::NAN;
        assert!(eval_legendre(0_i32, nan).is_nan());
        assert!(eval_legendre(0.0, nan).is_nan());
        assert!(crate::eval_chebyshev_t(0_i32, nan).is_nan());
        assert!(crate::eval_chebyshev_u(-1_i32, nan).is_nan());
        assert!(eval_laguerre(0_i32, nan).is_nan());
        assert!(crate::eval_genlaguerre(2_i32, nan, 1.0).is_nan());
        // alpha <= -1
        assert!(crate::eval_genlaguerre(2_i32, -1.5, 0.3).is_nan());
        assert!(crate::eval_genlaguerre(2.0, -1.5, 0.3).is_nan());
    }
}
