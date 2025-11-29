//! Translated from `scipy/special/orthogonal_eval.pxd`
//!
//! <https://github.com/scipy/scipy/blob/c16dc41/scipy/special/orthogonal_eval.pxd>

use crate::ffi::xsf as ffi;
use core::cmp::Ordering;
use core::f64::consts::{PI, SQRT_2};
use core::ops::Range;
use num_traits::Zero;

mod sealed {
    pub trait OrthoPolyArg:
        crate::xsf::Hyp2F1Arg + crate::xsf::Hyp1F1Arg + num_traits::NumOps<f64> + Clone
    {
        const ZERO: Self;

        fn is_finite(&self) -> bool;
        fn is_nan(&self) -> bool;

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
    }
}

/// Helper function for the generalized multiset coefficient
///
/// https://mathworld.wolfram.com/Multichoose.html
#[inline]
fn multiset(n: f64, k: f64) -> f64 {
    unsafe { ffi::binom(n + k - 1.0, k) }
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
    /// Corresponds to `eval_jacobi` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_jacobi(self, n: f64, a: f64, b: f64) -> Self {
        let a1 = a + 1.0;
        self.hyp2f1_sh(-n, n + a1 + b, a1) * multiset(a1, n)
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_legendre(self, n: f64) -> Self {
        self.hyp2f1_sh(-n, n + 1.0, 1.0)
    }
}

impl JacobiArg<i32> for f64 {
    /// Corresponds to `eval_jacobi_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_jacobi(self, n: i32, a: f64, b: f64) -> Self {
        match n.cmp(&0) {
            Ordering::Less => self.eval_jacobi(f64::from(n), a, b),
            Ordering::Equal => 1.0,
            Ordering::Greater => {
                let a1 = a + 1.0;
                // setting u = (x - 1) / 2 simplifies the recurrence
                let u = 0.5 * self - 0.5;
                let d0 = (a + b + 2.0) * u;

                if n == 1 {
                    a1 + d0
                } else {
                    let mut d = d0 / a1;
                    let mut p = d + 1.0;
                    for k in 1..n {
                        let k = f64::from(k);
                        let t = 2.0 * k + a + b;
                        d = (t + 2.0) * ((t + 1.0) * t * u * p + k * (k + b) * d)
                            / ((k + a1) * (k + a1 + b) * t);
                        p += d;
                    }
                    p * multiset(a1, f64::from(n))
                }
            },
        }
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_legendre(self, n: i32) -> Self {
        let n = if n < 0 { -n - 1 } else { n }; // symmetry
        let x = self;

        if n == 0 {
            1.0
        } else if n == 1 {
            x
        } else if x.abs() < 1e-5 {
            // Power series rather than recurrence due to loss of precision
            // http://functions.wolfram.com/Polynomials/LegendreP/02/
            let aa = n / 2;
            let a = f64::from(aa);
            let n_f64 = f64::from(n);

            let mut d = (-1.0_f64).powi(aa) * multiset(a + 1.0, -0.5);
            if n & 1 == 1 {
                // odd n
                d *= n_f64 * x;
            }

            let x2 = x * x;
            let mut p = 0.0;
            for kk in 0..=aa {
                let kk = f64::from(kk);
                p += d;
                // d *= -2 * x**2 * (a - kk) * (2*n + 1 - 2*a + 2*kk) / (
                //     (n + 1 - 2*a + 2*kk) * (n + 2 - 2*a + 2*kk))
                let u = 2.0 * (kk - a);
                let v = 1.0 + u + n_f64;
                d *= u * (v + n_f64) * x2 / (v * (1.0 + v));
            }
            p
        } else {
            let (mut p, mut d) = (x, x - 1.0);
            for k in 1..n {
                let k = f64::from(k);
                // ((2*k+1)/(k+1))*(x-1)*p + (k/(k+1)) * d
                d = ((2.0 * k + 1.0) * (x - 1.0) * p + k * d) / (k + 1.0);
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

///////////////////////////////
// Gegenbauer (Ultraspherical)
//

pub trait GegenbauerArg<N>: sealed::OrthoPolyArg {
    fn eval_gegenbauer(self, n: N, alpha: f64) -> Self;
}

impl<Z: sealed::OrthoPolyArg + From<f64>> GegenbauerArg<f64> for Z {
    /// Corresponds to `eval_gegenbauer` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_gegenbauer(self, n: f64, a: f64) -> Self {
        if a == 0.0 {
            if n.is_nan() || self.clone().is_nan() {
                self * n
            } else {
                Self::ZERO
            }
        } else if !self.clone().is_finite() {
            self * a
        } else {
            let aa = 2.0 * a;
            self.hyp2f1_sh(-n, n + aa, 0.5 + a) * multiset(aa, n)
        }
    }
}

impl GegenbauerArg<i32> for f64 {
    /// Corresponds to `eval_gegenbauer_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_gegenbauer(self, n: i32, alpha: f64) -> Self {
        if alpha.is_nan() || self.is_nan() {
            f64::NAN
        } else if n < 0 || alpha == 0.0 {
            // https://functions.wolfram.com/Polynomials/GegenbauerC3/03/01/02/0001/
            // https://functions.wolfram.com/Polynomials/GegenbauerC3/03/01/03/0012/
            0.0
        } else if n == 0 {
            1.0
        } else if self.abs() < 1e-5 {
            // Power series rather than recurrence due to loss of precision.
            // We use backwards iteration to avoid large powers of z in the initial term.
            // https://functions.wolfram.com/Polynomials/GegenbauerC3/02/0001/

            let m = n / 2;

            let mut term = (-1.0_f64).powi(m) * multiset(alpha, f64::from(m));
            if n == 2 * m + 1 {
                term *= 2.0 * (alpha + f64::from(m)) * self;
            }

            let c = -4.0 * self * self;
            let mut sum = term;
            for k in (1..=m).rev() {
                let k = f64::from(k);
                let nk = f64::from(n) - k;
                let nkk = nk - k;
                term *= c * k * (nk + alpha) / ((nkk + 1.0) * (nkk + 2.0));
                sum += term;
            }
            sum
        } else {
            // Recurrence relation
            // https://functions.wolfram.com/Polynomials/GegenbauerC3/17/01/01/01/0002/

            let aa = alpha * 2.0;
            let d0 = self - 1.0;
            let (mut p, mut d) = (self, d0);
            for k in 1..n {
                let k = f64::from(k);
                d = (k * d + (k + k + aa) * d0 * p) / (k + aa);
                p += d;
            }

            let n = f64::from(n);
            if (alpha / n) < 1e-8 {
                // avoid loss of precision
                aa / n * p
            } else {
                multiset(aa, n) * p
            }
        }
    }
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
    #[inline]
    fn eval_chebyshev_t(self, n: i32) -> Self {
        chebyshev_recurrence(self, 1.0, self, 0..n.abs())
    }

    #[inline]
    fn eval_chebyshev_u(self, n: i32) -> Self {
        if n < -1 {
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

///////////////////////////////
// Generalized Laguerre
// Laguerre
//

pub trait LaguerreArg<N>: sealed::OrthoPolyArg {
    fn eval_genlaguerre(self, n: N, a: f64) -> Self;
}

impl<Z: sealed::OrthoPolyArg> LaguerreArg<f64> for Z {
    /// Corresponds to `eval_genlaguerre` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_genlaguerre(self, n: f64, a: f64) -> Z {
        let a1 = a + 1.0;
        self.hyp1f1(-n, a1) * multiset(a1, n)
    }
}

impl LaguerreArg<i32> for f64 {
    /// Corresponds to `eval_genlaguerre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline]
    fn eval_genlaguerre(self, n: i32, a: f64) -> Self {
        if n < 0 {
            0.0
        } else {
            let a1 = a + 1.0;
            let (mut d, mut p) = (0.0, 1.0);
            for k in 0..n {
                let k = f64::from(k);
                d = (d * k - p * self) / (k + a1);
                p += d;
            }
            p * multiset(a1, f64::from(n))
        }
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
            PI.sqrt() * (scale * n).exp2() / unsafe { ffi::gamma(0.5 - 0.5 * n) }
        } else {
            0.0
        }
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
}
