//! Translated from `scipy/special/orthogonal_eval.pxd`
//!
//! <https://github.com/scipy/scipy/blob/c16dc41/scipy/special/orthogonal_eval.pxd>

use crate::ffi::xsf as ffi;
use num_complex::Complex64;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait OrthoPolyArg<N>: sealed::Sealed {
    // TODO: add the other classical orthogonaal polynomials (but skip the confusing shifted ones)
    fn eval_jacobi(self, n: N, alpha: f64, beta: f64) -> Self;
    fn eval_legendre(self, n: N) -> Self;
}

impl OrthoPolyArg<f64> for f64 {
    /// Corresponds to `eval_jacobi` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_jacobi(self, n: f64, a: f64, b: f64) -> Self {
        unsafe {
            ffi::binom(n + a, n) * ffi::hyp2f1(-n, 1.0 + n + a + b, 1.0 + a, 0.5 - 0.5 * self)
        }
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_legendre(self, n: f64) -> Self {
        unsafe { ffi::hyp2f1(-n, n + 1.0, 1.0, 0.5 * (1.0 - self)) }
    }
}

impl OrthoPolyArg<f64> for Complex64 {
    /// Corresponds to `eval_jacobi` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_jacobi(self, n: f64, a: f64, b: f64) -> Self {
        unsafe {
            ffi::binom(n + a, n) * ffi::hyp2f1_1(-n, 1.0 + n + a + b, 1.0 + a, 0.5 - 0.5 * self)
        }
    }

    /// Corresponds to `eval_legendre_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_legendre(self, n: f64) -> Self {
        unsafe { ffi::hyp2f1_1(-n, n + 1.0, 1.0, 0.5 * (1.0 - self)) }
    }
}

impl OrthoPolyArg<i32> for f64 {
    /// Corresponds to from `eval_jacobi_l` in `scipy/special/orthogonal_eval.pxd`
    #[inline(always)]
    fn eval_jacobi(self, n: i32, a: f64, b: f64) -> Self {
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
    fn eval_legendre(self, n: i32) -> Self {
        // symmetry
        let n = if n < 0 { -n - 1 } else { n };
        let x = self;

        if n == 0 {
            1.0
        } else if n == 1 {
            x
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
            let mut p = x;
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
/// # Notes
///
/// The Jacobi polynomials can be defined via the Gauss hypergeometric function ${}_2F_1$ as
///
/// $$
/// P_n^{(\alpha, \beta)}(z) =
///   \binom{\alpha + n}{n} \\;
///   \hyp{2}{1}{-n,\\, \alpha + \beta + n + 1}{\alpha + 1}{\Bigg\|\\, \frac{1 - z}{2}}
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
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function,
///   $\hyp{2}{1}{a_1,\\ a_2}{b}{\big\|\\,z}$
/// - [`eval_legendre`](crate::eval_legendre): Evaluate Legendre polynomial $P_n$
#[inline]
pub fn eval_jacobi<N, T>(n: N, alpha: f64, beta: f64, z: T) -> T
where
    T: OrthoPolyArg<N>,
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
/// # Notes
///
/// The Legendre polynomials can be defined via the Gauss hypergeometric function ${}_2F_1$ as
///
/// $$
/// P_n(z) = \hyp{2}{1}{-n,\\, n+1}{1}{\Bigg\|\\, \frac{1 - z}{2}}
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
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function,
///   $\hyp{2}{1}{a_1,\\ a_2}{b}{\big\|\\,z}$
/// - [`eval_jacobi`](crate::eval_jacobi): Evaluate Jacobi polynomial $P_n^{(\alpha, \beta)}$,
///   a generalization of Legendre polynomials
#[inline]
pub fn eval_legendre<N, T>(n: N, z: T) -> T
where
    T: OrthoPolyArg<N>,
{
    z.eval_legendre(n)
}

#[cfg(test)]
mod tests {
    use crate::{eval_jacobi, eval_legendre, np_assert_allclose};
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
}
