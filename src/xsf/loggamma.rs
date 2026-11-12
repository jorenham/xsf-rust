//! Translated into pure Rust from `xsf/loggamma.h` and `xsf/zlog1.h` (scipy/xsf@f2a68ec). Original
//! header comment appears below.
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

use crate::utils::{LN_PI, LN_SQRT_2PI, is_nonpos_int};
use crate::xsf::c64::{C64, ViaC64, c64};

use core::f64::consts::PI;
use num_complex::Complex;

const LOGGAMMA_SMALLX: f64 = 7.0;
const LOGGAMMA_SMALLY: f64 = 7.0;
const LOGGAMMA_TAYLOR_RADIUS: f64 = 0.2;

fn loggamma_stirling(z: C64) -> C64 {
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
    let rz = 1.0 / z;
    let rzz = rz / z;

    (z - 0.5) * z.ln() - z + LN_SQRT_2PI + rz * rzz.evalpoly(&coeffs)
}

fn loggamma_recurrence(mut z: C64) -> C64 {
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
    loggamma_stirling(z) - shiftprod.ln() - f64::from(signflips) * 2.0 * PI * C64::I
}

fn loggamma_taylor(mut z: C64) -> C64 {
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
    z * z.evalpoly(&coeffs)
}

/// [`zlog1_c64`] for [`num_complex::Complex`]
pub(crate) fn zlog1(z: Complex<f64>) -> Complex<f64> {
    z.via(zlog1_c64)
}

/// `xsf::detail::zlog1` from `xsf/zlog1.h` (original author: Josh Wilson, 2016)
fn zlog1_c64(mut z: C64) -> C64 {
    /* Compute log, paying special attention to accuracy around 1. We
     * implement this ourselves because some systems (most notably the
     * Travis CI machines) are weak in this regime. */
    let mut coeff = c64(-1.0, 0.0);
    let mut res = C64::ZERO;

    if (z - 1.0).norm() > 0.1 {
        return z.ln();
    }

    z -= 1.0;
    for n in 1..17 {
        coeff *= -z;
        res += coeff / f64::from(n);
        if (res / coeff).norm() < f64::EPSILON {
            break;
        }
    }
    res
}

/// `xsf::loggamma` for complex `z`
pub(crate) fn cloggamma(z: C64) -> C64 {
    // Compute the principal branch of log-Gamma

    if z.re.is_nan() || z.im.is_nan() {
        return C64::NAN;
    }
    if z.im == 0.0 && is_nonpos_int(z.re) {
        return C64::NAN;
    }
    if z.re > LOGGAMMA_SMALLX || z.im.abs() > LOGGAMMA_SMALLY {
        return loggamma_stirling(z);
    }
    if (z - 1.0).norm() < LOGGAMMA_TAYLOR_RADIUS {
        return loggamma_taylor(z);
    }
    if (z - 2.0).norm() < LOGGAMMA_TAYLOR_RADIUS {
        // Recurrence relation and the Taylor series around 1.
        return zlog1_c64(z - 1.0) + loggamma_taylor(z - 1.0);
    }
    if z.re < 0.1 {
        // Reflection formula; see Proposition 3.1 in [1]
        let tmp = (2.0 * PI).copysign(z.im) * (0.5 * z.re + 0.25).floor();
        return c64(LN_PI, tmp) - z.sinpi().ln() - cloggamma(1.0 - z);
    }
    if !z.im.is_sign_negative() {
        // z.imag() >= 0 but is not -0.0
        return loggamma_recurrence(z);
    }
    loggamma_recurrence(z.conj()).conj()
}

/// `xsf::rgamma` for complex `z`
fn crgamma(z: C64) -> C64 {
    // Compute 1/Gamma(z) using loggamma.
    if z.im == 0.0 && is_nonpos_int(z.re) {
        // Zeros at 0, -1, -2, ...
        return C64::ZERO;
    }
    (-cloggamma(z)).exp()
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
        self.via(cloggamma)
    }

    #[inline]
    fn xsf_rgamma(self) -> Self {
        self.via(crgamma)
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
