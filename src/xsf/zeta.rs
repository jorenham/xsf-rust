//! The complex `riemann_zeta` is translated into pure Rust from `xsf/zeta.h` (xsf v0.2.2).
//!
//! Complex riemann-zeta function implementation based on Python implementation
//! written by Matt Haberland (@mdhaber) in:
//! <https://colab.research.google.com/drive/1zMDSAJlXCLRqMMtJ0e9nDGjQ8iZCnmn5?usp=sharing>
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `zeta_euler_maclaurin` returns NaN if the number of terms doesn't fit in a
//! `size_t` (UB in C++).

use std::f64::consts::{LN_2, PI, TAU};

use num_complex::Complex;

use crate::utils::{LN_PI, LN_TAU};
use crate::xsf::complex::{cdiv, rsub};
use crate::xsf::exp::cexp;
use crate::xsf::log::clog;

const NAN: Complex<f64> = Complex::new(f64::NAN, f64::NAN);

/// `std::pow(x, z)` for real `x > 0` and complex `z`, as in libstdc++
fn rpow(x: f64, z: Complex<f64>) -> Complex<f64> {
    debug_assert!(x > 0.0);
    // std::polar(pow(x, z.real()), z.imag() * log(x))
    let rho = x.powf(z.re);
    let theta = z.im * x.ln();
    Complex::new(rho * theta.cos(), rho * theta.sin())
}

/// Log of absolute value of expansion coefficients for Euler-Maclaurin
/// summation formula. log(|B2k / (2k)!|)
///
/// See <https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms>
///
/// Generated with the script
///
/// ```python
/// import numpy as np
/// from mpmath import mp
///
/// mp.dps = 10000
///
/// results = []
/// for k in range(51):
///     results.append(
///         float(
///             mp.log(abs(mp.bernoulli(2*k)/(mp.factorial(2*k))))
///         )
///     )
/// ```
const ZETA_EM_LOG_ABS_COEFF_LOOKUP: [f64; 51] = [
    0.0,
    -2.484_906_649_788_000_4,
    -6.579_251_212_010_101,
    -1.031_692_083_029_347e1,
    -1.400_580_028_440_740_5e1,
    -1.768_462_940_266_784e1,
    -2.136_113_156_007_322_2e1,
    -2.503_707_050_291_142_3e1,
    -2.871_287_059_984_694_8e1,
    -3.238_863_619_752_229_5e1,
    -3.606_439_319_366_539e1,
    -3.974_014_804_199_518_4e1,
    -4.341_590_235_365_616e1,
    -4.709_165_653_118_148_5e1,
    -5.076_741_067_517_639e1,
    -5.444_316_481_078_909e1,
    -5.811_891_894_430_628e1,
    -6.179_467_307_729_959e1,
    -6.547_042_721_016_194e1,
    -6.914_618_134_299_154e1,
    -7.282_193_547_581_296e1,
    -7.649_768_960_863_234e1,
    -8.017_344_374_145_12e1,
    -8.384_919_787_426_993e1,
    -8.752_495_200_708_863e1,
    -9.120_070_613_990_733e1,
    -9.487_646_027_272_602e1,
    -9.855_221_440_554_472e1,
    -1.022_279_685_383_634e2,
    -1.059_037_226_711_821e2,
    -1.095_794_768_040_007_8e2,
    -1.132_552_309_368_194_7e2,
    -1.169_309_850_696_381_7e2,
    -1.206_067_392_024_568_5e2,
    -1.242_824_933_352_755_5e2,
    -1.279_582_474_680_942_4e2,
    -1.316_340_016_009_129_4e2,
    -1.353_097_557_337_316_3e2,
    -1.389_855_098_665_503e2,
    -1.426_612_639_993_69e2,
    -1.463_370_181_321_877e2,
    -1.500_127_722_650_064e2,
    -1.536_885_263_978_251e2,
    -1.573_642_805_306_437_5e2,
    -1.610_400_346_634_624_5e2,
    -1.647_157_887_962_811_5e2,
    -1.683_915_429_290_998_4e2,
    -1.720_672_970_619_185_4e2,
    -1.757_430_511_947_372_3e2,
    -1.794_188_053_275_559e2,
    -1.830_945_594_603_746e2,
];

/// Complex log of expansion coefficients for Euler-Maclaurin summation formula.
#[allow(clippy::cast_precision_loss)]
fn zeta_em_log_coeff(n: usize) -> Complex<f64> {
    let mut result: Complex<f64> = if n < 50 {
        ZETA_EM_LOG_ABS_COEFF_LOOKUP[n].into()
    } else {
        /* Asymptotic formula
         * Uses https://dlmf.nist.gov/24.11#E1 to approximate B_{2n} and
         * Stirling's approximation for (2n)!.
         */
        (LN_2 - 2.0 * n as f64 * LN_TAU).into()
    };
    if n % 2 == 0 {
        /* B_{2n}/(2n)! is negative for even n. This contributes a term
         * pi*i when taking the log. */
        result += PI * Complex::I;
    }
    result
}

/// Compute riemann_zeta for complex input z using the Euler-Maclaurin formula.
/// Computation of individual terms in expansion are logarithmized to avoid
/// overflow. TODO: only logarithmize when necessary.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]
fn zeta_euler_maclaurin(z: Complex<f64>) -> Complex<f64> {
    if z == Complex::ONE {
        /* Return NaN at pole since value depends on how z approaches 1.0. */
        return NAN;
    }
    let n = (z.im.abs() / 4.0).max(50.0);
    if n >= usize::MAX as f64 {
        // (unlike in xsf) the number of terms doesn't fit in a `size_t`
        return NAN;
    }
    let n = n as usize;
    let m = n;
    let mut result = Complex::ZERO;
    for i in 1..n {
        let term = rpow(i as f64, -z);
        result += term;
        // When z.real() > 1, series converges and we can consider early termination
        if z.re > 1.0 && term.norm() / result.norm() <= f64::EPSILON {
            return result;
        }
    }
    let big_n = n as f64;
    let b = rpow(big_n, -z);
    result += b * (0.5 + cdiv(big_n.into(), z - 1.0));
    /* The terms of the Euler-Maclaurin
     * expansion below are T(k, n) = B2k/(2k)! * n^(1 - z - 2k) * z(z+1)...(z+2k-2).
     * We work with logarithms to avoid overflow in all cases at the expense of
     * some accuracy. At the start of iteration k:
     *     log_poch will equal log(z(z+1)...(z+2k-2))
     *     log_factor will equal log(n^(1 - z - 2k))
     * These are updated one extra time after the loop completes for use in the
     * Euler-Maclaurin error estimate.
     */
    let mut log_poch = clog(z);
    let mut log_factor = -(z + 1.0) * big_n.ln();
    for k in 1..=m {
        let term = cexp(zeta_em_log_coeff(k) + log_factor + log_poch);
        result += term;
        if term.norm() / result.norm() <= f64::EPSILON {
            return result;
        }
        log_poch += clog(z + (2 * k - 1) as f64) + clog(z + (2 * k) as f64);
        log_factor -= 2.0 * big_n.ln();
    }
    /* Euler-maclaurin absolute error estimate.
     * The error is bounded above by |(z + 2m + 1)/(z.real + 2m + 1) * T(m+1, n)|
     * See https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms
     */
    let mut error = cexp(zeta_em_log_coeff(m + 1) + log_factor + log_poch).norm();
    error *= ((z + 2.0 * m as f64 + 1.0) / (z.re + 2.0 * m as f64 + 1.0)).norm();
    // convert to relative error estimate
    error /= result.norm();
    if error > 1e-1 {
        /* If error estimate predicts we don't even get 1 digit of precision, return NaN
         * and signal no result */
        return NAN;
    }
    result
}

/// Lookup table of coefficients for Algorithm 2 from Borwein 1995
/// Borwein, Peter B.. “An efficient algorithm for the Riemann zeta function.” (1995).
///
/// Stores coefficients as dk / dn, where dn is the final coefficient.
///
/// Generated with the Python script:
///
/// ```python
/// import numpy as np
/// import math
/// from mpmath import mp
///
/// mp.dps = 1000
/// n = 50
///
/// coeffs = []
/// S = mp.zero
/// for i in range(n + 1):
///     num = math.factorial(n + i - 1) * 4**i
///     den = math.factorial(n - i) * math.factorial(2*i)
///     S += mp.mpf(num) / mp.mpf(den)
///     coeffs.append(S*n)
///
/// dn = coeffs[-1]
/// coeffs = [float(dk/dn) for dk in coeffs[:-1]]
/// coeffs = np.asarray(coeffs)
/// ```
const ZETA_BORWEIN_COEFF: [f64; 50] = [
    1.055_507_836_138_287_8e-38,
    5.278_594_688_527_578e-35,
    4.401_468_732_204_496_3e-32,
    1.467_453_546_497_519e-29,
    2.617_862_196_688_831e-27,
    2.900_097_799_958_025e-25,
    2.184_440_361_492_933e-23,
    1.189_097_731_291_329_6e-21,
    4.887_139_616_687_227_6e-20,
    1.567_225_369_880_273_4e-18,
    4.022_931_234_264_572e-17,
    8.435_973_533_351_745e-16,
    1.469_296_379_914_116e-14,
    2.154_874_705_457_190_2e-13,
    2.691_953_053_753_512_4e-12,
    2.892_540_916_276_848_4e-11,
    2.695_750_569_369_985_6e-10,
    2.194_772_239_130_839e-9,
    1.570_782_293_700_57e-8,
    9.936_187_220_749_133e-8,
    5.581_721_578_217_702e-7,
    2.796_271_112_037_765e-6,
    1.253_886_254_275_813e-5,
    5.049_261_002_939_051e-5,
    1.831_288_445_970_366_6e-4,
    5.997_690_328_552_426e-4,
    1.778_050_108_246_096_8e-3,
    4.781_802_283_665_968e-3,
    1.169_043_228_713_167_1e-2,
    2.603_430_292_953_596_4e-2,
    5.292_298_247_275_485_6e-2,
    9.842_471_411_648_858e-2,
    1.678_961_079_654_034_4e-1,
    2.635_042_919_476_862_6e-1,
    3.819_442_810_600_665e-1,
    5.137_731_385_068_898e-1,
    6.452_925_385_418_66e-1,
    7.625_449_322_050_362e-1,
    8.556_063_057_019_102e-1,
    9.210_560_628_865_25e-1,
    9.616_100_580_028_147e-1,
    9.835_905_431_606_953e-1,
    9.939_187_229_336_753e-1,
    9.980_782_525_166_648e-1,
    9.994_930_141_459_889e-1,
    9.998_914_788_445_85e-1,
    9.999_819_091_990_203e-1,
    9.999_977_981_288_319e-1,
    9.999_998_260_580_315e-1,
    9.999_999_933_099_243e-1,
];

/// Compute riemann_zeta for complex input z using Algorithm 2 from Borwein 1995.
#[allow(clippy::cast_sign_loss)]
fn zeta_borwein(z: Complex<f64>) -> Complex<f64> {
    let mut result = Complex::ZERO;
    // Sum in reverse order because smaller terms come later.
    for k in (0..50).rev() {
        let sign = (-1.0_f64).powi(k);
        let den = rpow(f64::from(k + 1), z);
        let term = cdiv((sign * (ZETA_BORWEIN_COEFF[k as usize] - 1.0)).into(), den);
        result += term;
    }
    cdiv(-result, rsub(1.0, rpow(2.0, rsub(1.0, z))))
}

/// Compute riemann zeta for complex z and real part >= 0
fn zeta_right_halfplane(z: Complex<f64>) -> Complex<f64> {
    if z == Complex::ONE {
        return NAN;
    }
    /* Cutoff for using Euler-MacLaurin chosen based on cursory empirical search.
     * TODO: Choose cutoffs in a more principled way. */
    if z.re < 50.0 && z.im.abs() > 50.0 {
        if z.re >= 0.0 && z.re < 2.5 && z.im.abs() > 1e9 {
            /* Euler-MacLaurin summation starts to take an unreasonable amount of time in this
             * region, so just give up and return NaN instead. */
            return NAN;
        }
        return zeta_euler_maclaurin(z);
    }
    zeta_borwein(z)
}

/// exp(pi*z) for complex z.
fn exppi(z: Complex<f64>) -> Complex<f64> {
    let x = z.re;
    let y = z.im;
    let factor1 = Complex::new(crate::cospi(y), crate::sinpi(y));
    let factor2 = (PI * x).exp();
    factor1 * factor2
}

/// log(sinpi(z)) using sin(z) = (exp(i*pi*z) - exp(-i*pi*z)) / 2i
///
/// No attempt is made to choose any particular branch of the logarithm.
/// This is an internal function and the intent is that this that the
/// result of log(sinpi(z)) will be added to other terms, and the sum
/// will then be exponentiated, making the choice of a specific branch
/// unnecessary.
fn logsinpi(z: Complex<f64>) -> Complex<f64> {
    let result = clog(crate::sinpi(z));
    // If it doesn't overflow, just do the regular calculation.
    if result.re.is_finite() && !result.im.is_finite() {
        return result;
    }
    /* Otherwise factor before taking log. This is where we may end up
     * taking a branch other than the principal branch. */
    /* Calculating log((exp(i*pi*z) - exp(-i*pi*z)) / 2i). Factor out term
     * with larger magnitude before taking log. */
    if z.im > 0.0 {
        /* if z.imag() > 0 then, exp(-i*pi*z) has greatest magnitude. Factor it
         * out to get:
         * log(exp(-i*pi*z)*((exp(2*i*pi*z) - 1.0)/(2i)) =
         * log(exp(-i*pi*z)) + log((exp(2*i*pi*z) - 1.0)/(2i)) =
         * -i*pi*z + log((exp(2*i*pi*z) - 1.0)/(2i)) */
        return -Complex::I * PI * z
            + clog(cdiv(exppi(2.0 * z * Complex::I) - 1.0, 2.0 * Complex::I));
    }
    /* if z.imag() < 0 then, exp(i*pi*z) has greatest magnitude. Factor similarly
     * to above */
    Complex::I * PI * z
        + clog(cdiv(
            rsub(1.0, exppi(-2.0 * z * Complex::I)),
            2.0 * Complex::I,
        ))
}

/// Leading factor in reflection formula for zeta function.
/// zeta(z) = 2^z * pi^(z-1) * sin(pi*z/2) * gamma(1 - z) * zeta(1 - z)
/// This computes 2^z * pi^(z - 1) * sin(pi*z/2) * gamma(1 - z)
///
/// Computation is logarithimized to prevent overflow.
/// TODO: Complexify the cephes zeta_reflection implementation, which uses
/// the lanczos approximation for the gamma function.
fn zeta_reflection_factor_with_logs(z: Complex<f64>) -> Complex<f64> {
    let t1 = z * LN_2;
    let t2 = (z - 1.0) * LN_PI;
    let t3 = logsinpi(z / 2.0);
    let t4 = crate::loggamma(rsub(1.0, z));
    cexp(t1 + t2 + t3 + t4)
}

fn zeta_reflection(z: Complex<f64>) -> Complex<f64> {
    let mut factor = 2.0 * rpow(TAU, z - 1.0) * crate::sinpi(z / 2.0) * crate::gamma(rsub(1.0, z));
    if !factor.re.is_finite() || !factor.im.is_finite() {
        // Try again with logs if standard calculation had overflow.
        factor = zeta_reflection_factor_with_logs(z);
    }
    let result = zeta_right_halfplane(rsub(1.0, z));
    /* zeta tends to 1.0 as real part tends to +inf. In cases where
     * the real part of zeta tends to -inf, then zeta(1 - z) in the
     * reflection formula will tend to 1.0. Factor overflows then,
     * factor * result below will become NaN. In this case, we just
     * return factor to preserve complex infinity. Only zeta(1 - z) == 1.0
     * is handled because this is the only practical case where we should
     *  expect zeta(1 - z) == x for a real number x when z is not on the
     * real line. */
    if result == Complex::ONE {
        factor
    } else {
        factor * result
    }
}

fn criemann_zeta(z: Complex<f64>) -> Complex<f64> {
    if z.im == 0.0 {
        return crate::xsf::cephes::riemann_zeta(z.re).into();
    }
    if z.re >= 0.5 {
        return zeta_right_halfplane(z);
    }
    zeta_reflection(z)
}

pub trait ZetaArg: crate::sealed::Sealed {
    fn riemann_zeta(self) -> Self;
    fn zeta(self, q: f64) -> Self;
}

impl ZetaArg for f64 {
    #[inline]
    fn riemann_zeta(self) -> Self {
        crate::xsf::cephes::riemann_zeta(self)
    }

    #[inline]
    #[allow(clippy::float_cmp)]
    fn zeta(self, q: f64) -> Self {
        if q == 1.0 {
            return self.riemann_zeta();
        }
        crate::xsf::cephes::zeta(self, q)
    }
}

impl ZetaArg for num_complex::Complex<f64> {
    #[inline]
    fn riemann_zeta(self) -> Self {
        criemann_zeta(self)
    }

    #[inline]
    #[allow(clippy::float_cmp)]
    fn zeta(self, q: f64) -> Self {
        if q == 1.0 {
            return self.riemann_zeta();
        }
        if self.im == 0.0 {
            num_complex::Complex::new(self.re.zeta(q), 0.0)
        } else {
            // Complex input for Hurwitz Zeta is not currently supported.
            num_complex::Complex::new(f64::NAN, f64::NAN)
        }
    }
}

/// Riemann zeta function $\zeta(z)$ for real or complex $z$
///
/// Corresponds to [`scipy.special.zeta(z)`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html
///
/// # Definition
///
/// $$ \zeta(z) = \sum_{k=1}^\infty k^{-z} $$
///
/// # See also
/// - [`zeta`]: Hurwitz zeta function $\zeta(z, q)$
/// - [`zetac`]: Riemann zeta function minus one, $\zeta(x) - 1$
#[must_use]
#[inline]
pub fn riemann_zeta<T: ZetaArg>(z: T) -> T {
    z.riemann_zeta()
}

/// Hurwitz zeta function $\zeta(z, q)$ for real or complex $z$
///
/// Corresponds to [`scipy.special.zeta(z, q)`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html
///
/// # Definition
///
/// $$ \zeta(z, q) = \sum_{k=0}^\infty (k + q)^{-z} $$
///
/// # See also
/// - [`riemann_zeta`]: Riemann zeta function $\zeta(z)$
/// - [`zetac`]: Riemann zeta function minus one, $\zeta(x) - 1$
#[doc(alias = "hurwitz_zeta")]
#[must_use]
#[inline]
pub fn zeta<T: ZetaArg>(z: T, q: f64) -> T {
    z.zeta(q)
}

/// Riemann zeta function minus one, $\zeta(x) - 1$
///
/// Corresponds to [`scipy.special.zetac`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zetac.html
///
/// # Definition
///
/// $$ \zeta(x) - 1 = \sum_{k=2}^\infty k^{-x} $$
///
/// # See also
/// - [`riemann_zeta`]: Riemann zeta function $\zeta(z)$
/// - [`zeta`]: Hurwitz zeta function $\zeta(z, q)$
#[doc(alias = "riemann_zetac")]
#[must_use]
#[inline]
pub fn zetac(x: f64) -> f64 {
    crate::xsf::cephes::zetac(x)
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_riemann_zeta_f64() {
        xsref::test("riemann_zeta", "d-d", |x| crate::riemann_zeta(x[0]));
    }

    #[test]
    fn test_riemann_zeta_c64() {
        xsref::test("riemann_zeta", "cd-cd", |x| {
            crate::riemann_zeta(c64(x[0], x[1]))
        });
    }

    #[test]
    fn test_riemann_zeta_c64_extreme() {
        // (k + 1)^z overflows in the Borwein sum, and 1 / inf = 0 (C99 Annex G)
        let w = crate::riemann_zeta(c64(1000.0, 0.5));
        assert_eq!((w.re, w.im), (1.0, -3.170_083_028_960_87e-302));
        // the reflection factor overflows, and zeta(1 - z) is finite
        assert_eq!(crate::riemann_zeta(c64(-1000.0, 0.5)).re, f64::INFINITY);
        // too many Euler-Maclaurin terms (undefined behavior in xsf)
        for im in [7.4e19, 1e300, f64::INFINITY] {
            let w = crate::riemann_zeta(c64(3.0, im));
            assert!(w.re.is_nan() && w.im.is_nan());
        }
    }

    #[test]
    fn test_zeta_f64() {
        xsref::test("zeta", "d_d-d", |x| crate::zeta(x[0], x[1]));
    }

    #[test]
    fn test_zeta_large_q() {
        // asymptotic expansion for q > 1e8 (reference value from mpmath)
        let (x, q) = (36.0, 1.1e8);
        let expected = 1.016_688_811_413_570_4e-283;
        let rel = (crate::zeta(x, q) - expected).abs() / expected;
        assert!(rel < 1e-15, "rel = {rel:e}");

        assert_eq!(crate::zeta(f64::INFINITY, q), 0.0);
        assert_eq!(crate::zeta(f64::INFINITY, f64::INFINITY), 0.0);
    }

    #[test]
    fn test_zeta_near_1_tiny_q() {
        // the tail sum ~ 1 / (x - 1) isn't negligible (reference values from mpmath)
        for (x, q, expected) in [
            (1.0 + f64::EPSILON, 1e-15, 5.503_599_627_370_504e15),
            (1.000_000_000_001, 1e-17, 1.000_009_999_150_220_6e17),
        ] {
            let rel = (crate::zeta(x, q) - expected).abs() / expected;
            assert!(rel < 1e-15, "zeta({x}, {q}): rel = {rel:e}");
        }
    }

    #[test]
    fn test_zeta_negative_q() {
        // reference values from mpmath, by summing the terms with q + k < 0 directly
        for (x, q, expected) in [
            (4.0, -431.636_697_197_915_6, 63.992_140_081_094_475),
            (2.0, -1e10 - 0.5, 9.869_604_400_989_358),
            (31.0, -7.5, 1.593_135_539_215_446_3e-29),
            (31.0, -0.500_000_000_000_000_1, 3.304_063_430_367_368e-5),
        ] {
            let rel = (crate::zeta(x, q) - expected).abs() / expected;
            assert!(rel < 1e-15, "zeta({x}, {q}): rel = {rel:e}");
        }
        // all terms underflow
        assert_eq!(crate::zeta(1e300, 1.5), 0.0);
        assert_eq!(crate::zeta(f64::INFINITY, 1.5), 0.0);
        assert_eq!(crate::zeta(f64::INFINITY, -0.5), f64::INFINITY);
    }

    #[test]
    fn test_zetac_f64() {
        xsref::test("zetac", "d-d", |x| crate::zetac(x[0]));
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_zetac_large_x() {
        // zetac(x) = 2^-x to double precision, which doesn't underflow before x = 1075
        for x in [127.0, 130.0, 1000.0, 1074.0] {
            assert_eq!(crate::zetac(x), (-x).exp2());
        }
    }
}
