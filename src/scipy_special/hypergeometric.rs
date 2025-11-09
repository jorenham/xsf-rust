use core::f64::consts::TAU;
use num_complex::Complex64;

const LN_MAX: f64 = 709.782_712_893_384; // f64::MAX.ln()
const LN_MIN: f64 = -708.396_418_532_264_1; // f64::MIN_POSITIVE.ln()

/// Asymptotic expansion for $I_{v-1}(2\sqrt{z}) \Gamma(v)$ for real $z>0$ and $v \to +\infty$.
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_asy`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L59-L100>
#[inline(always)]
fn _hyp0f1_asy(v: f64, z: f64) -> f64 {
    let arg = z.sqrt();
    let v1 = (v - 1.0).abs();
    let x = 2.0 * arg / v1;
    let p1 = (1.0 + x * x).sqrt();
    let eta = p1 + x.ln() - unsafe { crate::ffi::xsf::log1p(p1) };

    let arg_exp_i = unsafe { crate::ffi::xsf::gammaln(v) } - 0.5 * (TAU * p1 * v1).ln();
    let gs = unsafe { crate::ffi::xsf::gammasgn(v) };

    let _v1_eta = v1 * eta;
    let arg_exp_i = arg_exp_i + _v1_eta;
    let arg_exp_k = arg_exp_i - _v1_eta;

    // large-v asymptotic correction, DLMF 10.41.10
    let qv = 1.0 / (24.0 * v1 * p1);
    let q2 = 1.0 / (1.0 + x * x); // p1^-2

    let uv1 = (3.0 - 5.0 * q2) * qv;
    let uv2 = (81.0 - (462.0 + 385.0 * q2) * q2) * qv.powi(2) / 2.0;
    let uv3 = (3037.5 - (36960.3 + (76576.5 - 42542.5 * q2) * q2) * q2) * qv.powi(3) / 3.0;
    let u_corr_i = 1.0 + uv1 + uv2 + uv3;

    let _v1_log_arg = unsafe { crate::ffi::xsf::xlogy(v1, arg) };
    let out = gs * (arg_exp_i - _v1_log_arg).exp() * u_corr_i;

    if v - 1.0 < 0.0 {
        // DLMF 10.27.2: I_{-v} = I_{v} + (2/pi) sin(pi*v) K_v
        let u_corr_k = 1.0 - uv1 + uv2 - uv3;
        out + 2.0
            * gs
            * unsafe { crate::ffi::xsf::sinpi(v1) }
            * (arg_exp_k + _v1_log_arg).exp()
            * u_corr_k
    } else {
        out
    }
}

/// Real-valued kernel
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_real`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L31-L56>
#[inline(always)]
fn _hyp0f1_real(v: f64, z: f64) -> f64 {
    if v <= 0.0 && v.fract() == 0.0 {
        // poles
        f64::NAN
    } else if z == 0.0 && v != 0.0 {
        // zeros
        1.0
    } else if z.abs() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)
        1.0 + z / v + z * z / (2.0 * v * (v + 1.0))
    } else if z > 0.0 {
        let arg = z.sqrt();
        let arg_exp = unsafe { crate::ffi::xsf::xlogy(1.0 - v, arg) + crate::ffi::xsf::gammaln(v) };
        let bess_val = unsafe { crate::ffi::xsf::cyl_bessel_i(v - 1.0, 2.0 * arg) };

        if arg_exp > LN_MAX || bess_val == 0.0 || arg_exp < LN_MIN || bess_val.is_infinite() {
            // overflow or underflow
            _hyp0f1_asy(v, z)
        } else {
            arg_exp.exp() * bess_val * unsafe { crate::ffi::xsf::gammasgn(v) }
        }
    } else {
        let arg = (-z).sqrt();
        let bess_val = unsafe { crate::ffi::xsf::cephes_jv(v - 1.0, 2.0 * arg) };

        arg.powf(1.0 - v) * unsafe { crate::ffi::xsf::gamma(v) } * bess_val
    }
}

/// Complex valued kernel
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_cmplx`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L106-L136>
#[inline(always)]
fn _hyp0f1_cmplx(v: f64, z: Complex64) -> Complex64 {
    if v <= 0.0 && v.fract() == 0.0 {
        // poles
        f64::NAN.into()
    } else if z.re == 0.0 && z.im == 0.0 && v != 0.0 {
        // zeros
        1.0.into()
    } else if z.norm() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)

        // need to do computations in this order, for otherwise $v\approx -z \ll 1$
        // it can lose precision (as was reported for 32-bit linux, see scipy/scipy#6365)
        let t1 = 1.0 + z / v;
        let t2 = z * z / (2.0 * v * (v + 1.0));
        t1 + t2
    } else {
        let (arg, r) = if z.re > 0.0 {
            let arg = z.sqrt();
            let bess_val = unsafe { crate::ffi::xsf::cyl_bessel_i_1(v - 1.0, 2.0 * arg) };
            (arg, bess_val)
        } else {
            let arg = (-z).sqrt();
            let bess_val = unsafe { crate::ffi::xsf::cyl_bessel_j_1(v - 1.0, 2.0 * arg) };
            (arg, bess_val)
        };

        r * unsafe { crate::ffi::xsf::gamma(v) } * arg.powf(1.0 - v)
    }
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait HypergeometricArg: sealed::Sealed {
    fn hyp0f0(self) -> Self;
    fn hyp1f0(self, a: f64) -> Self;
    fn hyp0f1(self, b: f64) -> Self;
}

impl HypergeometricArg for f64 {
    #[inline(always)]
    fn hyp0f0(self) -> Self {
        self.exp()
    }

    #[inline(always)]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self == 1.0 {
            f64::NAN
        } else {
            (1.0 - self).powf(-a)
        }
    }

    #[inline(always)]
    fn hyp0f1(self, b: f64) -> Self {
        _hyp0f1_real(b, self)
    }
}

impl HypergeometricArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn hyp0f0(self) -> Self {
        self.exp()
    }

    #[inline(always)]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self.re == 1.0 && self.im == 0.0 {
            f64::NAN.into()
        } else {
            (1.0 - self).powf(-a)
        }
    }

    #[inline(always)]
    fn hyp0f1(self, b: f64) -> Self {
        _hyp0f1_cmplx(b, self)
    }
}

/// Hypergeometric function $_0F_0\left[ \middle\| z\right]$ for real or complex $z$
///
/// This function is currently not implemented in SciPy, but is straightforward to evaluate.
/// Both [`f64`] and [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _0F_0\left[ \middle\| z\right] = \sum\_{n=0}^\infty {z^n \over n!},
/// $$
///
/// which is just the exponential function, $e^z$.
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function $_1F_0\left[a\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
pub fn hyp0f0<T: HypergeometricArg>(z: T) -> T {
    z.hyp0f0()
}

/// Hypergeometric function $_1F_0\left[a\middle\| z\right]$ for real or complex $z$
///
/// This function is currently not implemented in SciPy, but is straightforward to evaluate.
/// Both [`f64`] and [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _1F_0\left[a\middle\| z\right] = \sum\_{n=0}^\infty \rpow a n {z^n \over n!}
/// $$
///
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
/// It's also equal to $(1 - z)^{-a}$.
///
/// # See also
/// - [`hyp0f0`]: Hypergeometric function $_0F_0\left[ \middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
pub fn hyp1f0<T: HypergeometricArg>(a: f64, z: T) -> T {
    z.hyp1f0(a)
}

/// Confluent hypergeometric limit function $_0F_1\left[b\middle\| z\right]$ for real or complex $z$
///
/// This is a translation of the [`scipy.special.hyp0f1`][hyp0f1] Cython implementation into Rust,
/// using the same Cephes-based [scipy/xsf][xsf] FFI functions as SciPy.
///
/// [hyp0f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp0f1.html
/// [xsf]: https://github.com/scipy/xsf
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _0F_1\left[b\middle\| z\right] = \sum\_{n=0}^\infty {1\over \rpow b n} {z^n \over n!}.
/// $$
///
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
///
/// It's also the limit as $a \to \infty$ of $\hyp 1 1 a b z$, and satisfies the
/// differential equation $f\'\'(z) + b f\'(z) = f(z)$. See [^1] for more information.
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function, $_1F_0\left[a\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
/// [^1]: Weisstein, Eric W. "Confluent Hypergeometric Limit Function." From MathWorld -- A Wolfram
/// Resource. <https://mathworld.wolfram.com/ConfluentHypergeometricLimitFunction.html>
///
pub fn hyp0f1<T: HypergeometricArg>(b: f64, z: T) -> T {
    z.hyp0f1(b)
}

#[cfg(test)]
mod tests {
    use crate::{hyp0f1, np_assert_allclose};
    use core::f64::consts::{E, PI};
    use num_complex::c64;

    #[test]
    fn test_hyp0f0_f64() {
        np_assert_allclose!([crate::hyp0f0(0.0)], [1.0], atol = 1e-15);
        np_assert_allclose!([crate::hyp0f0(1.0)], [E], atol = 1e-15);
        np_assert_allclose!([crate::hyp0f0(-1.0)], [1.0 / E], atol = 1e-15);
    }

    #[test]
    fn test_hyp0f0_c64() {
        np_assert_allclose!(
            [crate::hyp0f0(c64(0.0, 0.0))],
            [c64(1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!([crate::hyp0f0(c64(1.0, 0.0))], [c64(E, 0.0)], atol = 1e-15);
        np_assert_allclose!(
            [crate::hyp0f0(c64(0.0, PI))],
            [c64(-1.0, 0.0)],
            atol = 1e-15
        );
    }

    #[test]
    fn test_hyp1f0_f64() {
        assert_eq!(crate::hyp1f0(-1.0, -1.0), 2.0);
        assert_eq!(crate::hyp1f0(0.0, -1.0), 1.0);
        assert_eq!(crate::hyp1f0(1.0, -1.0), 0.5);
        assert_eq!(crate::hyp1f0(2.0, -1.0), 0.25);

        assert_eq!(crate::hyp1f0(-1.0, 0.0), 1.0);
        assert_eq!(crate::hyp1f0(0.0, 0.0), 1.0);
        assert_eq!(crate::hyp1f0(1.0, 0.0), 1.0);
        assert_eq!(crate::hyp1f0(2.0, 0.0), 1.0);

        assert_eq!(crate::hyp1f0(-1.0, 1.0), 0.0);
        assert_eq!(crate::hyp1f0(0.0, 1.0), 1.0);
        assert!(crate::hyp1f0(1.0, 1.0).is_nan());
        assert!(crate::hyp1f0(2.0, 1.0).is_nan());

        assert_eq!(crate::hyp1f0(-1.0, 2.0), -1.0);
        assert_eq!(crate::hyp1f0(0.0, 2.0), 1.0);
        assert_eq!(crate::hyp1f0(1.0, 2.0), -1.0);
        assert_eq!(crate::hyp1f0(2.0, 2.0), 1.0);
    }

    #[test]
    fn test_hyp1f0_c64() {
        assert_eq!(crate::hyp1f0(-1.0, c64(-1.0, 0.0)), c64(2.0, 0.0));
        assert_eq!(crate::hyp1f0(0.0, c64(-1.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(1.0, c64(-1.0, 0.0)), c64(0.5, 0.0));
        assert_eq!(crate::hyp1f0(2.0, c64(-1.0, 0.0)), c64(0.25, 0.0));

        assert_eq!(crate::hyp1f0(-1.0, c64(0.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(0.0, c64(0.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(1.0, c64(0.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(2.0, c64(0.0, 0.0)), c64(1.0, 0.0));

        assert_eq!(crate::hyp1f0(-1.0, c64(1.0, 0.0)), c64(0.0, 0.0));
        assert_eq!(crate::hyp1f0(0.0, c64(1.0, 0.0)), c64(1.0, 0.0));
        assert!(crate::hyp1f0(1.0, c64(1.0, 0.0)).is_nan());
        assert!(crate::hyp1f0(2.0, c64(1.0, 0.0)).is_nan());

        np_assert_allclose!(
            [crate::hyp1f0(-1.0, c64(2.0, 0.0))],
            [c64(-1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!(
            [crate::hyp1f0(0.0, c64(2.0, 0.0))],
            [c64(1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!(
            [crate::hyp1f0(1.0, c64(2.0, 0.0))],
            [c64(-1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!(
            [crate::hyp1f0(2.0, c64(2.0, 0.0))],
            [c64(1.0, 0.0)],
            atol = 1e-15
        );
    }

    /// Corresponds to `scipy.special.tests.test_basic.TestHyper.test_hyp0f1`
    /// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/tests/test_basic.py#L3144-L3168>
    #[test]
    fn test_hyp0f1() {
        // scalar input
        np_assert_allclose!([hyp0f1(2.5, 0.5)], [1.21482702689997], rtol = 1e-12);
        np_assert_allclose!([hyp0f1(2.5, 0.0)], [1.0], rtol = 1e-15);

        // float input, expected values match mpmath
        let x = [-1.5, -1.0, 0.0, 1.0, 1.5].map(|z| hyp0f1(3.0, z));
        let expected = [
            0.58493659229143,
            0.70566805723127,
            1.0,
            1.37789689539747,
            1.60373685288480,
        ];
        np_assert_allclose!(x, expected, rtol = 1e-12);

        // complex input
        let x = [-1.5, -1.0, 0.0, 1.0, 1.5].map(|re| hyp0f1(3.0, c64(re, 0.0)));
        np_assert_allclose!(x, expected.map(|v| c64(v, 0.0)), rtol = 1e-12);
    }

    /// Corresponds to `scipy.special.tests.test_basic.TestHyper.test_hyp0f1_gh5764`
    /// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/tests/test_basic.py#L3170-L3176>
    /// which looks like a regression test for scipy/scipy#5764
    #[test]
    fn test_hyp0f1_gh5764() {
        // Just checks the point that failed; there's a more systematic test in test_mpmath
        let res = hyp0f1(0.8, c64(0.5, 0.5));
        // The expected value was generated using mpmath
        np_assert_allclose!(
            [res],
            [c64(1.6139719776441115, 0.8089305406179067)],
            atol = 1.5e-7
        );
    }
}
