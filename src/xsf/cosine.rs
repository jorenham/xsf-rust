//! Translated into pure Rust from `xsf/cosine.h` (xsf main): translated into C++ by the xsf
//! developers in 2026, original author: Warren Weckesser, 2021.

use crate::xsf::cephes::polevl::polevl;

// M_PI64 is the 64 bit floating point representation of π, e.g., in Python:
//   >>> import math
//   >>> math.pi.hex()
//   '0x1.921fb54442d18p+1'
// It is used in the function cosine_cdf_pade_approx_at_neg_pi,
// which depends on this value being the 64 bit representation.
// Do not replace this with M_PI from math.h or NPY_PI from the
// numpy header files.
const M_PI64: f64 = core::f64::consts::PI; // 0x1.921fb54442d18p+1

// p and q (below) are the coefficients in the numerator and denominator
// polynomials (resp.) of the Pade approximation of
//     f(x) = (pi + x + sin(x))/(2*pi)
// at x=-pi.  The coefficients are ordered from lowest degree to highest.
// These values are used in the function cosine_cdf_pade_approx_at_neg_pi(x).
//
// These coefficients can be derived by using mpmath as follows:
//
//    import mpmath
//
//    def f(x):
//        return (mpmath.pi + x + mpmath.sin(x)) / (2*mpmath.pi)
//
//    # Note: 40 digits might be overkill; a few more digits than the default
//    # might be sufficient.
//    mpmath.mp.dps = 40
//    ts = mpmath.taylor(f, -mpmath.pi, 20)
//    p, q = mpmath.pade(ts, 9, 10)
//
// (A python script with that code is in
//  https://github.com/scipy/scipy/blob/v1.18.0/scipy/special/_precompute/cosine_cdf.py)
//
// The following are the values after converting to 64 bit floating point:
// p = [0.0,
//      0.0,
//      0.0,
//      0.026525823848649224,
//      0.0,
//      -0.0007883197097740538,
//      0.0,
//      1.0235408442872927e-05,
//      0.0,
//      -3.8360369451359084e-08]
// q = [1.0,
//      0.0,
//      0.020281047093125535,
//      0.0,
//      0.00020944197182753272,
//      0.0,
//      1.4162345851873058e-06,
//      0.0,
//      6.498171564823105e-09,
//      0.0,
//      1.6955280904096042e-11]

const COSINE_CDF_PADE_NUMER: [f64; 4] = [
    -3.836_036_945_135_908_4e-8,
    1.023_540_844_287_292_7e-5,
    -0.000_788_319_709_774_053_8,
    0.026_525_823_848_649_224,
];

const COSINE_CDF_PADE_DENOM: [f64; 6] = [
    1.695_528_090_409_604_2e-11,
    6.498_171_564_823_105e-9,
    1.416_234_585_187_305_8e-6,
    0.000_209_441_971_827_532_72,
    0.020_281_047_093_125_535,
    1.0,
];

const COSINE_INVCDF_P2_COEFFS: [f64; 6] = [
    -6.844_846_384_555_272_5e-9,
    3.490_093_422_701_228_4e-6,
    -0.000_305_397_129_071_151_67,
    0.009_350_454_384_541_677,
    -0.116_021_429_402_087_26,
    0.5,
];

const COSINE_INVCDF_Q2_COEFFS: [f64; 6] = [
    -5.579_679_571_562_129e-8,
    1.372_857_015_278_879_3e-5,
    -0.000_891_691_992_732_111_7,
    0.022_927_496_105_281_435,
    -0.252_876_192_137_507_84,
    1.0,
];

const COSINE_INVCDF_POLY_COEFFS: [f64; 7] = [
    1.191_166_794_908_291_5e-8,
    1.683_039_183_039_183e-7,
    43.0 / 17_248_000.0,
    1.0 / 25_200.0,
    1.0 / 1400.0,
    1.0 / 60.0,
    1.0,
];

fn cosine_cdf_pade_approx_at_neg_pi(x: f64) -> f64 {
    // Compute the CDF of the standard cosine distribution for x close to but
    // not less than -π.  A Pade approximant is used to avoid the loss of
    // precision that occurs in the formula 1/2 + (x + sin(x))/(2*pi) when
    // x is near -π.

    // M_PI64 is not exactly π.  In fact, float64(π - M_PI64) is
    // 1.2246467991473532e-16.  h is supposed to be x + π, so to compute
    // h accurately, we write the calculation as:
    let h = (x + M_PI64) + 1.224_646_799_147_353_2e-16;
    let h2 = h * h;
    let h3 = h2 * h;
    let numer = h3 * polevl(h2, &COSINE_CDF_PADE_NUMER);
    let denom = polevl(h2, &COSINE_CDF_PADE_DENOM);
    numer / denom
}

#[inline]
fn cosine_invcdf_p2(t: f64) -> f64 {
    polevl(t, &COSINE_INVCDF_P2_COEFFS)
}

#[inline]
fn cosine_invcdf_q2(t: f64) -> f64 {
    polevl(t, &COSINE_INVCDF_Q2_COEFFS)
}

fn cosine_invcdf_poly_approx(s: f64) -> f64 {
    // Part of the asymptotic expansion of the inverse function at p=0.
    //
    // See, for example, the wikipedia article "Kepler's equation"
    // (https://en.wikipedia.org/wiki/Kepler%27s_equation).  In particular, see the
    // series expansion for the inverse Kepler equation when the eccentricity e is 1.

    // p(s) = s + (1/60) * s**3 + (1/1400) * s**5 + (1/25200) * s**7 +
    //        (43/17248000) * s**9 + (1213/7207200000) * s**11 +
    //        (151439/12713500800000) * s**13 + ...
    //
    // Here we include terms up to s**13.
    let s2 = s * s;
    s * polevl(s2, &COSINE_INVCDF_POLY_COEFFS)
}

/// CDF of the standard (raised) cosine distribution, $\frac{\pi + x + \sin x}{2\pi}$
///
/// The support is $\[-\pi, \pi\]$; this returns 0 for $x < -\pi$ and 1 for $x \ge \pi$.
///
/// Corresponds to `xsf::cosine_cdf` in xsf main, which is used by
/// [`scipy.stats.cosine`][cosine].
///
/// [cosine]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cosine.html
///
/// # See also
/// - [`cosine_invcdf`]
#[must_use]
pub fn cosine_cdf(x: f64) -> f64 {
    // Compute the CDF of the standard cosine distribution.
    if x >= M_PI64 {
        return 1.0;
    }
    if x < -M_PI64 {
        return 0.0;
    }
    if x < -1.6 {
        return cosine_cdf_pade_approx_at_neg_pi(x);
    }
    0.5 + (x + x.sin()) / (2.0 * M_PI64)
}

// The CDF of the cosine distribution is
//     p = (pi + x + sin(x)) / (2*pi),
// We want the inverse of this.
//
// Move the factor 2*pi and the constant pi to express this as
//     pi*(2*p - 1) = x + sin(x)
// Then if f(x) = x + sin(x), and g(x) is the inverse of f, we have
//     x = g(pi*(2*p - 1)).

// The coefficients in the functions _p2 and _q2 are the coefficients in the
// Pade approximation at p=0.5 to the inverse of x + sin(x).
// The following steps were used to derive these:
// 1. Find the coefficients of the Taylor polynomial of x + sin(x) at x = 0.
//    A Taylor polynomial of order 22 was used.
// 2. "Revert" the Taylor coefficients to find the Taylor polynomial of the
//    inverse.  The method for series reversion is described at
//        https://en.wikipedia.org/wiki/Bell_polynomials#Reversion_of_series
// 3. Convert the Taylor coefficients of the inverse to the (11, 10) Pade
//    approximant. The coefficients of the Pade approximant (converted to 64
//    bit floating point) in increasing order of degree are:
//      p = [0.0, 0.5,
//           0.0, -0.11602142940208726,
//           0.0, 0.009350454384541677,
//           0.0, -0.00030539712907115167,
//           0.0, 3.4900934227012284e-06,
//           0.0, -6.8448463845552725e-09]
//      q = [1.0,
//           0.0, -0.25287619213750784,
//           0.0, 0.022927496105281435,
//           0.0, -0.0008916919927321117,
//           0.0, 1.3728570152788793e-05,
//           0.0, -5.579679571562129e-08]
//
// The nonzero values in p and q are used in the functions _p2 and _q2 below.
// The functions assume that the square of the variable is passed in, and
// _p2 does not include the outermost multiplicative factor.
// So to evaluate x = invcdf(p) for a given p, the following is used:
//        y = pi*(2*p - 1)
//        x = y * _p2(y**2) / _q2(y**2)

/// Inverse of [`cosine_cdf`], the quantile function of the standard (raised) cosine distribution
///
/// Returns NaN for `p` outside of $\[0, 1\]$.
///
/// Corresponds to `xsf::cosine_invcdf` in xsf main, which is used by
/// [`scipy.stats.cosine`][cosine].
///
/// [cosine]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cosine.html
///
/// # See also
/// - [`cosine_cdf`]
#[must_use]
#[allow(clippy::float_cmp)]
pub fn cosine_invcdf(mut p: f64) -> f64 {
    // Compute the inverse CDF of the standard cosine distribution.
    let mut sgn = 1.0;

    if p < 0.0 || p > 1.0 {
        return f64::NAN;
    }
    if p <= 1e-48 {
        return -M_PI64;
    }
    if p == 1.0 {
        return M_PI64;
    }

    if p > 0.5 {
        p = 1.0 - p;
        sgn = -1.0;
    }

    let mut x = if p < 0.0925 {
        // `f64::cbrt` instead of the less accurate `cephes::cbrt` used by `xsf::cbrt`
        cosine_invcdf_poly_approx((12.0 * M_PI64 * p).cbrt()) - M_PI64
    } else {
        let y = M_PI64 * (2.0 * p - 1.0);
        let y2 = y * y;
        y * cosine_invcdf_p2(y2) / cosine_invcdf_q2(y2)
    };
    // For p < 0.0018, the asymptotic expansion at p=0 is sufficiently
    // accurate that no more work is needed.  Similarly, for p > 0.42,
    // the Pade approximant is sufficiently accurate.  In between these
    // bounds, we refine the estimate with Halley's method.
    if 0.0018 < p && p < 0.42 {
        // Apply one iteration of Halley's method, with
        //    f(x)   = pi + x + sin(x) - y,
        //    f'(x)  = 1 + cos(x),
        //    f''(x) = -sin(x)
        // where y = 2*pi*p.
        let f0 = M_PI64 + x + x.sin() - 2.0 * M_PI64 * p;
        let f1 = 1.0 + x.cos();
        let f2 = -x.sin();
        x -= 2.0 * f0 * f1 / (2.0 * f1 * f1 - f0 * f2);
    }

    sgn * x
}

#[cfg(test)]
mod tests {
    use super::M_PI64;
    use crate::{cosine_cdf, cosine_invcdf};

    // reference values from scipy/special/tests/test_cosine_distr.py

    #[test]
    #[allow(clippy::approx_constant)]
    fn test_cosine_cdf() {
        assert_eq!(cosine_cdf(-4.0), 0.0);
        assert_eq!(cosine_cdf(0.0), 0.5);
        assert_eq!(cosine_cdf(M_PI64), 1.0);
        assert_eq!(cosine_cdf(4.0), 1.0);
        assert!(cosine_cdf(f64::NAN).is_nan());

        let xs = [
            3.1409, 2.25, -1.599, -1.601, -2.0, -3.0, -3.1415, -3.141_59, -M_PI64,
        ];
        let expected = [
            0.999_999_999_991_185,
            0.981_932_817_328_790_7,
            0.086_419_598_383_825_53,
            0.086_110_582_992_713,
            0.036_970_933_596_161_1,
            7.522_387_241_801_384e-5,
            2.109_869_685_443_648e-14,
            4.956_444_476_505_336e-19,
            4.871_934_450_264_861e-50,
        ];
        let actual = xs.map(cosine_cdf);
        xsref::np_assert_allclose!(&actual, &expected, rtol = 5e-15);
    }

    #[test]
    fn test_cosine_invcdf() {
        assert_eq!(cosine_invcdf(0.0), -M_PI64);
        assert_eq!(cosine_invcdf(0.5), 0.0);
        assert_eq!(cosine_invcdf(1.0), M_PI64);
        assert!(cosine_invcdf(-0.1).is_nan());
        assert!(cosine_invcdf(1.1).is_nan());
        assert!(cosine_invcdf(f64::NAN).is_nan());

        let ps = [
            1e-50,
            1e-14,
            1e-8,
            0.001_800_1,
            0.01,
            0.06,
            0.125,
            0.25,
            0.4,
            0.419,
            0.421,
            0.75,
            0.94,
            0.999_999_999_6,
        ];
        let expected = [
            -M_PI64,
            -3.141_520_413_705_845_4,
            -3.134_368_658_912_452_4,
            -2.732_563_923_138_336,
            -2.412_765_890_086_78,
            -1.788_124_497_533_015_7,
            -1.375_252_366_986_927_4,
            -0.831_711_193_579_736,
            -0.316_795_451_239_528_9,
            -0.255_860_256_269_199_06,
            -0.249_475_707_504_456_63,
            0.831_711_193_579_736,
            1.788_124_497_533_015_3,
            3.139_122_083_991_716_7,
        ];
        let actual = ps.map(cosine_invcdf);
        xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
    }
}
