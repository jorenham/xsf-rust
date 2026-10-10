//! Translated into pure Rust from `xsf/loggamma.h` and `xsf/zlog1.h` (xsf v0.2.2), which were
//! translated from Cython into C++ by SciPy developers in 2024 and 2023, respectively.
//! Original header comment appears below.
//!
//! An implementation of the principal branch of the logarithm of
//! Gamma. Also contains implementations of Gamma and 1/Gamma which are
//! easily computed from log-Gamma.
//!
//! Author: Josh Wilson
//!
//! Distributed under the same license as Scipy.
//!
//! References
//! ----------
//! \[1\] Hare, "Computing the Principal Branch of log-Gamma",
//!     Journal of Algorithms, 1997.
//!
//! \[2\] Julia,
//!     <https://github.com/JuliaLang/julia/blob/v0.6.4/base/special/gamma.jl>
//!     (unlike in xsf, a permalink to the last Julia release that contains it)
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::complex::cdiv;
use crate::xsf::evalpoly::cevalpoly;
use crate::xsf::exp::cexp;
use crate::xsf::log::clog;

use core::f64::consts::PI;
use num_complex::Complex;

const LOGGAMMA_SMALLX: f64 = 7.0;
const LOGGAMMA_SMALLY: f64 = 7.0;
const LOGGAMMA_HLOG2PI: f64 = 0.918_938_533_204_672_8; // log(2*pi)/2
const LOGGAMMA_LOGPI: f64 = 1.144_729_885_849_400_2; // log(pi)
const LOGGAMMA_TAYLOR_RADIUS: f64 = 0.2;

fn loggamma_stirling(z: Complex<f64>) -> Complex<f64> {
    /* Stirling series for log-Gamma
     *
     * The coefficients are B[2*n]/(2*n*(2*n - 1)) where B[2*n] is the
     * (2*n)th Bernoulli number. See (1.1) in [1].
     */
    let coeffs = [
        -0.029_550_653_594_771_242,
        0.006_410_256_410_256_41,
        -0.001_917_526_917_526_917_6,
        0.000_841_750_841_750_841_7,
        -0.000_595_238_095_238_095_3,
        0.000_793_650_793_650_793_7,
        -0.002_777_777_777_777_778,
        0.083_333_333_333_333_33,
    ];
    let rz = cdiv(Complex::new(1.0, 0.0), z);
    let rzz = cdiv(rz, z);

    (z - 0.5) * clog(z) - z + LOGGAMMA_HLOG2PI + rz * cevalpoly(&coeffs, rzz)
}

fn loggamma_recurrence(mut z: Complex<f64>) -> Complex<f64> {
    /* Backward recurrence relation.
     *
     * See Proposition 2.2 in [1] and the Julia implementation [2].
     *
     */
    let mut signflips = 0;
    let mut sb = false;
    let mut shiftprod = z;

    z += 1.0;
    while z.re <= LOGGAMMA_SMALLX {
        shiftprod *= z;
        let nsb = shiftprod.im.is_sign_negative();
        signflips += i32::from(nsb && !sb);
        sb = nsb;
        z += 1.0;
    }
    loggamma_stirling(z) - clog(shiftprod) - f64::from(signflips) * 2.0 * PI * Complex::I
}

fn loggamma_taylor(mut z: Complex<f64>) -> Complex<f64> {
    /* Taylor series for log-Gamma around z = 1.
     *
     * It is
     *
     * loggamma(z + 1) = -gamma*z + zeta(2)*z**2/2 - zeta(3)*z**3/3 ...
     *
     * where gamma is the Euler-Mascheroni constant.
     */
    let coeffs = [
        -0.043_478_266_053_040_26,
        0.045_454_556_293_204_67,
        -0.047_619_070_330_142_226,
        0.050_000_047_698_101_69,
        -0.052_631_679_379_616_66,
        0.055_555_767_627_403_614,
        -0.058_823_978_658_684_585,
        0.062_500_955_141_213_04,
        -0.066_668_705_882_420_46,
        0.071_432_946_295_361_33,
        -0.076_932_516_411_352_2,
        0.083_353_840_546_109,
        -0.090_954_017_145_829_04,
        0.100_099_457_512_781_8,
        -0.111_334_265_869_564_69,
        0.125_509_669_524_743_04,
        -0.144_049_896_768_846_1,
        0.169_557_176_997_408_2,
        -0.207_385_551_028_673_98,
        0.270_580_808_427_784_54,
        -0.400_685_634_386_531_43,
        0.822_467_033_424_113_2,
        -0.577_215_664_901_532_9,
    ];

    z -= 1.0;
    z * cevalpoly(&coeffs, z)
}

/// `xsf::detail::zlog1` from `xsf/zlog1.h` (original author: Josh Wilson, 2016)
fn zlog1(mut z: Complex<f64>) -> Complex<f64> {
    /* Compute log, paying special attention to accuracy around 1. We
     * implement this ourselves because some systems (most notably the
     * Travis CI machines) are weak in this regime. */
    let mut coeff = Complex::new(-1.0, 0.0);
    let mut res = Complex::new(0.0, 0.0);

    if (z - 1.0).norm() > 0.1 {
        return clog(z);
    }

    z -= 1.0;
    for n in 1..17 {
        coeff *= -z;
        res += coeff / f64::from(n);
        if cdiv(res, coeff).norm() < f64::EPSILON {
            break;
        }
    }
    res
}

/// `xsf::loggamma` for complex `z`
#[allow(clippy::float_cmp)]
fn cloggamma(z: Complex<f64>) -> Complex<f64> {
    // Compute the principal branch of log-Gamma

    if z.re.is_nan() || z.im.is_nan() {
        return Complex::new(f64::NAN, f64::NAN);
    }
    if z.re <= 0.0 && z.im == 0.0 && z.re == z.re.floor() {
        return Complex::new(f64::NAN, f64::NAN);
    }
    if z.re > LOGGAMMA_SMALLX || z.im.abs() > LOGGAMMA_SMALLY {
        return loggamma_stirling(z);
    }
    if (z - 1.0).norm() < LOGGAMMA_TAYLOR_RADIUS {
        return loggamma_taylor(z);
    }
    if (z - 2.0).norm() < LOGGAMMA_TAYLOR_RADIUS {
        // Recurrence relation and the Taylor series around 1.
        return zlog1(z - 1.0) + loggamma_taylor(z - 1.0);
    }
    if z.re < 0.1 {
        // Reflection formula; see Proposition 3.1 in [1]
        let tmp = (2.0 * PI).copysign(z.im) * (0.5 * z.re + 0.25).floor();
        // `1.0 - z` as in C++, where the imaginary part is `-z.im` instead of `0.0 - z.im`
        let omz = Complex::new(1.0 - z.re, -z.im);
        return Complex::new(LOGGAMMA_LOGPI, tmp) - clog(crate::sinpi(z)) - cloggamma(omz);
    }
    if !z.im.is_sign_negative() {
        // z.imag() >= 0 but is not -0.0
        return loggamma_recurrence(z);
    }
    loggamma_recurrence(z.conj()).conj()
}

/// `xsf::rgamma` for complex `z`
#[allow(clippy::float_cmp)]
fn crgamma(z: Complex<f64>) -> Complex<f64> {
    // Compute 1/Gamma(z) using loggamma.
    if z.re <= 0.0 && z.im == 0.0 && z.re == z.re.floor() {
        // Zeros at 0, -1, -2, ...
        return Complex::new(0.0, 0.0);
    }
    cexp(-cloggamma(z))
}

pub trait LogGammaArg: crate::sealed::Sealed {
    fn xsf_loggamma(self) -> Self;
    fn xsf_rgamma(self) -> Self;
}

impl LogGammaArg for f64 {
    #[inline]
    fn xsf_loggamma(self) -> f64 {
        if self < 0.0 {
            return f64::NAN;
        }
        crate::xsf::cephes::lgam(self)
    }

    #[inline]
    fn xsf_rgamma(self) -> f64 {
        crate::xsf::cephes::rgamma(self)
    }
}

impl LogGammaArg for Complex<f64> {
    #[inline]
    fn xsf_loggamma(self) -> Self {
        cloggamma(self)
    }

    #[inline]
    fn xsf_rgamma(self) -> Self {
        crgamma(self)
    }
}

/// Principal branch of the logarithm of `gamma(z)`
///
/// Corresponds to [`scipy.special.loggamma`][loggamma] in scipy
///
/// [loggamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.loggamma.html
///
/// # See also
/// - [`gamma`](crate::gamma): Gamma function
#[doc(alias = "lgamma", alias = "ln_gamma", alias = "log_gamma")]
#[inline]
pub fn loggamma<T: LogGammaArg>(z: T) -> T {
    z.xsf_loggamma()
}

/// Reciprocal Gamma function `1 / gamma(z)`
///
/// Corresponds to [`scipy.special.rgamma`][rgamma] in scipy
///
/// [rgamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rgamma.html
///
/// # See also
/// - [`gamma`](crate::gamma): Gamma function
#[inline]
pub fn rgamma<T: LogGammaArg>(z: T) -> T {
    z.xsf_rgamma()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_loggamma_f64() {
        xsref::test("loggamma", "d-d", |x| crate::loggamma(x[0]));
    }

    #[test]
    fn test_loggamma_c64() {
        xsref::test("loggamma", "cd-cd", |x| crate::loggamma(c64(x[0], x[1])));
    }

    #[test]
    fn test_loggamma_c64_negative_axis() {
        // the sign of the zero imaginary part selects the side of the branch cut
        let w = crate::loggamma(c64(-2.5, 0.0));
        assert_eq!(crate::loggamma(c64(-2.5, -0.0)), w.conj());
        // mpmath: loggamma(-2.5 + 1e-30j)
        assert!((w.re / -0.056_243_716_497_674_05 - 1.0).abs() < 1e-13);
        assert_eq!(w.im, -3.0 * core::f64::consts::PI);

        // poles of gamma, zeros of rgamma
        let w = crate::loggamma(c64(-2.0, 0.0));
        assert!(w.re.is_nan() && w.im.is_nan());
        assert_eq!(crate::rgamma(c64(-2.0, 0.0)), c64(0.0, 0.0));
    }

    #[test]
    fn test_rgamma_f64() {
        xsref::test("rgamma", "d-d", |x| crate::rgamma(x[0]));
    }

    #[test]
    fn test_rgamma_c64() {
        xsref::test("rgamma", "cd-cd", |x| crate::rgamma(c64(x[0], x[1])));
    }
}
