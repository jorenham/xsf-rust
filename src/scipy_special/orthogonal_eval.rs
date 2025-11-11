//! Translated from `scipy/special/orthogonal_eval.pxd`
//!
//! <https://github.com/scipy/scipy/blob/c16dc41/scipy/special/orthogonal_eval.pxd>

use core::f64::consts::{PI, SQRT_2};
use core::ops::Mul;

use crate::ffi::xsf as ffi;
use num_complex::Complex64;
use num_traits::Zero;

mod sealed {
    pub trait Sealed: crate::xsf::Hyp2F1Arg + crate::xsf::Hyp1F1Arg {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

///////////////////////////////
// Jacobi
// Legendre
//

pub trait JacobiArg<N>: sealed::Sealed {
    // TODO: add the other classical orthogonal polynomials (but skip the confusing shifted ones)
    fn eval_jacobi(&self, n: N, alpha: f64, beta: f64) -> Self;
    fn eval_legendre(&self, n: N) -> Self;
}

impl JacobiArg<f64> for f64 {
    /// Corresponds to `eval_jacobi` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_jacobi(&self, n: f64, a: f64, b: f64) -> Self {
        unsafe {
            ffi::binom(n + a, n) * ffi::hyp2f1(-n, 1.0 + n + a + b, 1.0 + a, 0.5 - 0.5 * self)
        }
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_legendre(&self, n: f64) -> Self {
        unsafe { ffi::hyp2f1(-n, n + 1.0, 1.0, 0.5 - 0.5 * self) }
    }
}

impl JacobiArg<f64> for Complex64 {
    /// Corresponds to `eval_jacobi` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_jacobi(&self, n: f64, a: f64, b: f64) -> Self {
        unsafe {
            ffi::binom(n + a, n) * ffi::hyp2f1_1(-n, 1.0 + n + a + b, 1.0 + a, 0.5 - 0.5 * self)
        }
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_legendre(&self, n: f64) -> Self {
        unsafe { ffi::hyp2f1_1(-n, n + 1.0, 1.0, 0.5 - 0.5 * self) }
    }
}

impl JacobiArg<i32> for f64 {
    /// Corresponds to `eval_jacobi_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_jacobi(&self, n: i32, a: f64, b: f64) -> Self {
        if n < 0 {
            self.eval_jacobi(n as f64, a, b)
        } else if n == 0 {
            1.0
        } else {
            // setting u = (x - 1) / 2 simplifies the recurrence
            let u = 0.5 * self - 0.5;

            if n == 1 {
                (a + 1.0) + (a + b + 2.0) * u
            } else {
                let mut d = (a + b + 2.0) * u / (a + 1.0);
                let mut p = d + 1.0;
                for k in 1..n {
                    let k = k as f64;
                    let t = 2.0 * k + a + b;
                    d = (t + 2.0) * ((t + 1.0) * t * u * p + k * (k + b) * d)
                        / ((k + a + 1.0) * (k + a + b + 1.0) * t);
                    p += d;
                }
                p * unsafe { ffi::binom(n as f64 + a, n as f64) }
            }
        }
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_legendre(&self, n: i32) -> Self {
        // symmetry
        let n = if n < 0 { -n - 1 } else { n };
        let x = self;

        if n == 0 {
            1.0
        } else if n == 1 {
            *self
        } else if x.abs() < 1e-5 {
            // Power series rather than recurrence due to loss of precision
            // http://functions.wolfram.com/Polynomials/LegendreP/02/
            let aa = n / 2;
            let a = aa as f64;
            let n_f64 = n as f64;

            let mut d = if aa % 2 == 0 { 1.0 } else { -1.0 };
            if n == 2 * aa {
                d *= -2.0 / unsafe { ffi::beta(a + 1.0, -0.5) };
            } else {
                d *= 2.0 * x / unsafe { ffi::beta(a + 1.0, 0.5) };
            }

            let x2 = x * x;
            let mut p = 0.0;
            for kk in 0..=aa {
                let kk = kk as f64;
                p += d;
                // d *= -2 * x**2 * (a - kk) * (2*n + 1 - 2*a + 2*kk) / (
                //     (n + 1 - 2*a + 2*kk) * (n + 2 - 2*a + 2*kk))
                let u = 2.0 * (kk - a);
                let v = 1.0 + u + n_f64;
                d *= u * (v + n_f64) * x2 / (v * (1.0 + v));
                if d.abs() < 1e-20 * p.abs() {
                    // converged
                    break;
                }
            }
            p
        } else {
            let mut d = x - 1.0;
            let mut p = *x;
            for k in 1..n {
                let k = k as f64;
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
pub fn eval_jacobi<N, Z>(n: N, alpha: f64, beta: f64, z: Z) -> Z
where
    Z: JacobiArg<N>,
{
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
pub fn eval_legendre<N, Z>(n: N, z: Z) -> Z
where
    Z: JacobiArg<N>,
{
    z.eval_legendre(n)
}

///////////////////////////////
// Generalized Laguerre
// Laguerre
//

pub trait LaguerreArg<N>: sealed::Sealed {
    fn eval_genlaguerre(&self, n: N, a: f64) -> Self;
    fn eval_laguerre(&self, n: N) -> Self;
}

impl<Z> LaguerreArg<f64> for Z
where
    Z: sealed::Sealed + Mul<f64, Output = Z> + Clone,
{
    /// Corresponds to `eval_genlaguerre` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_genlaguerre(&self, n: f64, a: f64) -> Z {
        crate::hyp1f1(-n, 1.0 + a, self.clone()) * unsafe { ffi::binom(n + a, n) }
    }

    #[inline(always)]
    fn eval_laguerre(&self, n: f64) -> Z {
        crate::hyp1f1(-n, 1.0, self.clone())
    }
}

impl LaguerreArg<i32> for f64 {
    /// Corresponds to `eval_genlaguerre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_genlaguerre(&self, n: i32, a: f64) -> Self {
        if n < 0 {
            0.0
        } else {
            let (mut d, mut p) = (0.0, 1.0);
            for k in 0..n {
                let k = k as f64;
                d = (d * k - p * self) / (k + a + 1.0);
                p += d;
            }
            p * unsafe { ffi::binom(n as f64 + a, n as f64) }
        }
    }

    #[inline(always)]
    fn eval_laguerre(&self, n: i32) -> Self {
        if n < 0 {
            0.0
        } else {
            let (mut d, mut p) = (0.0, 1.0);
            for k in 0..n {
                let k = k as f64;
                d = (d * k - p * self) / (k + 1.0);
                p += d;
            }
            p
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
pub fn eval_genlaguerre<N, Z>(n: N, alpha: f64, z: Z) -> Z
where
    Z: LaguerreArg<N>,
{
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
pub fn eval_laguerre<N, Z>(n: N, z: Z) -> Z
where
    Z: LaguerreArg<N>,
{
    z.eval_laguerre(n)
}

///////////////////////////////
// Hermite H  (physicists')
// Hermite He (probabilists')
//

pub trait HermiteArg<N>: sealed::Sealed {
    fn eval_hermite_he(&self, n: N) -> Self;
    fn eval_hermite_h(&self, n: N) -> Self;
}

impl HermiteArg<f64> for f64 {
    #[inline(always)]
    fn eval_hermite_he(&self, n: f64) -> f64 {
        if self.is_nan() || n.is_nan() {
            f64::NAN
        } else if n == 0.0 {
            1.0
        } else if self.is_zero() {
            if n < 0.0 || n % 2.0 == 0.0 {
                (PI * n.exp2()).sqrt() / unsafe { ffi::gamma(0.5 - 0.5 * n) }
            } else {
                0.0
            }
        } else {
            let half_n = 0.5 * n;
            let y = 0.5 * self * self;
            let c = if self.is_sign_positive() {
                half_n.exp2()
            } else {
                (self / y.sqrt()).powf(n)
            };
            c * crate::hypu(-half_n, 0.5, y)
        }
    }

    #[inline(always)]
    fn eval_hermite_h(&self, n: f64) -> f64 {
        if self.is_nan() || n.is_nan() {
            f64::NAN
        } else if n == 0.0 {
            1.0
        } else if self.is_zero() {
            if n < 0.0 || n % 2.0 == 0.0 {
                PI.sqrt() * n.exp2() / unsafe { ffi::gamma(0.5 - 0.5 * n) }
            } else {
                0.0
            }
        } else {
            let y = self * self;
            let c = if self.is_sign_positive() {
                n.exp2()
            } else {
                (2.0 * self / y.sqrt()).powf(n)
            };
            c * crate::hypu(-0.5 * n, 0.5, y)
        }
    }
}

impl HermiteArg<u32> for f64 {
    /// Corresponds to `eval_hermitenorm` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_hermite_he(&self, n: u32) -> f64 {
        if self.is_nan() {
            f64::NAN
        } else {
            let (mut y1, mut y2) = (1.0, 0.0);
            for k in (1..=n).rev() {
                (y1, y2) = (self * y1 - (k as f64) * y2, y1);
            }
            y1
        }
    }

    /// Corresponds to `eval_hermite` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_hermite_h(&self, n: u32) -> f64 {
        (n as f64 / 2.0).exp2() * (SQRT_2 * self).eval_hermite_he(n)
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
        eval_genlaguerre, eval_hermite_h, eval_hermite_he, eval_jacobi, eval_laguerre,
        eval_legendre, np_assert_allclose,
    };
    use num_complex::c64;

    #[test]
    fn test_eval_jacobi_eq_legendre() {
        // manual test cases
        let xs = [-5.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 5.0];
        let ns = [0, 1, 2, 3, 4, 8, 15];
        for n in ns {
            let p0 = xs.map(|x| crate::legendre_p(n, x));

            let pi = xs.map(|x| eval_jacobi(n, 0.0, 0.0, x));
            np_assert_allclose!(pi, p0, rtol = 1e-13);

            let pf = xs.map(|x| eval_jacobi(n as f64, 0.0, 0.0, x));
            np_assert_allclose!(pf, p0, rtol = 1e-12, atol = f64::EPSILON);

            let pc = xs.map(|x| eval_jacobi(n as f64, 0.0, 0.0, c64(x, 0.0)).re);
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
            [-491.0, -1.0, -0.173125, -0.375, -0.539375, 4.0, 584.0],
        ];
        let expect_0_1 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [-8.0, -2.0, -0.65, -0.5, -0.35, 1.0, 7.0],
            [67.0, 3.0, -0.375, -0.5, -0.575, 1.0, 57.0],
            [-584.0, -4.0, 0.539375, 0.375, 0.173125, 1.0, 491.0],
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
                -535.9375, -2.1875, 0.214375, 0.0, -0.214375, 2.1875, 535.9375,
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

                let p_f64 = xs.map(|x| eval_jacobi(n as f64, alpha, beta, x));
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

            let pf = xs.map(|x| eval_legendre(n as f64, x));
            np_assert_allclose!(pf, p0, rtol = 1e-12, atol = f64::EPSILON);
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
                74.33333333333333,
                5.666666666666667,
                1.3151666666666667,
                1.0,
                0.7148333333333333,
                -0.6666666666666667,
                2.6666666666666667,
            ],
        ];
        let expect_1 = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [7.0, 3.0, 2.1, 2.0, 1.9, 1.0, -3.0],
            [30.5, 6.5, 3.305, 3.0, 2.705, 0.5, 0.5],
            [
                104.83333333333333,
                12.166666666666667,
                4.620166666666667,
                4.0,
                3.4198333333333333,
                -0.1666666666666667,
                3.1666666666666667,
            ],
        ];
        let expect_h = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [6.5, 2.5, 1.6, 1.5, 1.4, 0.5, -3.5],
            [26.875, 4.875, 2.13, 1.875, 1.63, -0.125, 1.875],
            [
                88.64583333333333,
                8.479166666666667,
                2.6426666666666667,
                2.1875,
                1.7673333333333333,
                -0.6041666666666667,
                3.2291666666666667,
            ],
        ];
        let expect_nh = [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [5.5, 1.5, 0.6, 0.5, 0.4, -0.5, -4.5],
            [20.375, 2.375, 0.53, 0.375, 0.23, -0.625, 5.375],
            [
                61.770833333333333,
                3.6041666666666667,
                0.5126666666666667,
                0.3125,
                0.13733333333333333,
                -0.4791666666666667,
                1.3541666666666667,
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

                let p_f64 = xs.map(|x| eval_genlaguerre(n as f64, alpha, x));
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
                74.33333333333333,
                5.666666666666667,
                1.3151666666666667,
                1.0,
                0.7148333333333333,
                -0.666666666666667,
                2.6666666666666667,
            ],
            [
                205.375,
                8.708333333333333,
                1.4306708333333333,
                1.0,
                0.6293375,
                -0.625,
                -1.2916666666666667,
            ],
            [
                515.5833333333333,
                12.883333333333333,
                1.5516875833333333,
                1.0,
                0.5483540833333333,
                -0.4666666666666667,
                -3.1666666666666667,
            ],
            [
                1203.7430555555556,
                18.509722222222222,
                1.6783963347222222,
                1.0,
                0.4717286680555556,
                -0.2569444444444444,
                -2.0902777777777778,
            ],
            [
                2653.4107142857147,
                25.9765873015873,
                1.8109809264087307,
                1.0,
                0.39931075970238106,
                -0.04047619047619044,
                0.3253968253968256,
            ],
            [
                5580.251612103175,
                35.757167658730154,
                1.9496297057145342,
                1.0,
                0.33095370539707347,
                0.15399305555555556,
                2.235739087301587,
            ],
        ];

        for (i, &n) in ns.iter().enumerate() {
            let p_i32 = xs.map(|x| eval_laguerre(n, x));
            np_assert_allclose!(p_i32, expect[i], rtol = 1e-15, atol = f64::EPSILON);

            let p_f64 = xs.map(|x| eval_laguerre(n as f64, x));
            np_assert_allclose!(p_f64, expect[i], rtol = 2e-14, atol = f64::EPSILON);
        }
    }

    #[test]
    fn test_eval_hermite_he() {
        let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];

        let y0_expect = xs.map(|_| 1.0);
        let y0_actual_u32 = xs.map(|x| eval_hermite_he(0, x));
        let y0_actual_f64 = xs.map(|x| eval_hermite_he(0.0, x));
        np_assert_allclose!(y0_actual_u32, y0_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y0_actual_f64, y0_expect, rtol = 1e-15, atol = f64::EPSILON);

        let y1_expect = xs.map(|x| x);
        let y1_actual_u32 = xs.map(|x| eval_hermite_he(1, x));
        let y1_actual_f64 = xs.map(|x| eval_hermite_he(1.0, x));
        np_assert_allclose!(y1_actual_u32, y1_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y1_actual_f64, y1_expect, rtol = 1e-15, atol = f64::EPSILON);

        let y2_expect = xs.map(|x| x * x - 1.0);
        let y2_actual_u32 = xs.map(|x| eval_hermite_he(2, x));
        let y2_actual_f64 = xs.map(|x| eval_hermite_he(2.0, x));
        np_assert_allclose!(y2_actual_u32, y2_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y2_actual_f64, y2_expect, rtol = 1e-15, atol = f64::EPSILON);

        let y3_expect = xs.map(|x| x * x * x - 3.0 * x);
        let y3_actual_u32 = xs.map(|x| eval_hermite_he(3, x));
        let y3_actual_f64 = xs.map(|x| eval_hermite_he(3.0, x));
        np_assert_allclose!(y3_actual_u32, y3_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y3_actual_f64, y3_expect, rtol = 1e-15, atol = f64::EPSILON);
    }

    #[test]
    fn test_eval_hermite_h() {
        let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];

        let y0_expect = xs.map(|_| 1.0);
        let y0_actual_u32 = xs.map(|x| eval_hermite_h(0, x));
        let y0_actual_f64 = xs.map(|x| eval_hermite_h(0.0, x));
        np_assert_allclose!(y0_actual_u32, y0_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y0_actual_f64, y0_expect, rtol = 1e-15, atol = f64::EPSILON);

        let y1_expect = xs.map(|x| 2.0 * x);
        let y1_actual_u32 = xs.map(|x| eval_hermite_h(1, x));
        let y1_actual_f64 = xs.map(|x| eval_hermite_h(1.0, x));
        np_assert_allclose!(y1_actual_u32, y1_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y1_actual_f64, y1_expect, rtol = 1e-15, atol = f64::EPSILON);

        let y2_expect = xs.map(|x| 4.0 * x * x - 2.0);
        let y2_actual_u32 = xs.map(|x| eval_hermite_h(2, x));
        let y2_actual_f64 = xs.map(|x| eval_hermite_h(2.0, x));
        np_assert_allclose!(y2_actual_u32, y2_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y2_actual_f64, y2_expect, rtol = 1e-15, atol = f64::EPSILON);

        let y3_expect = xs.map(|x| 8.0 * x * x * x - 12.0 * x);
        let y3_actual_u32 = xs.map(|x| eval_hermite_h(3, x));
        let y3_actual_f64 = xs.map(|x| eval_hermite_h(3.0, x));
        np_assert_allclose!(y3_actual_u32, y3_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y3_actual_f64, y3_expect, rtol = 1e-15, atol = f64::EPSILON);
    }
}
