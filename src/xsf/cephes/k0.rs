//! Translated into pure Rust from `xsf/cephes/k0.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the intervals and the error conditions are corrected in the comments.

use crate::xsf::cephes::chbevl::chbevl;
use crate::xsf::cephes::i0::i0;

/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
 * in the interval [0,2].  The odd order coefficients are all
 * zero; only the even order coefficients are listed.
 *
 * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
 */
const K0_A: [f64; 10] = [
    1.374_465_435_613_523e-16,
    4.259_816_142_796_61e-14,
    1.034_969_525_763_384_2e-11,
    1.904_516_377_220_209e-9,
    2.534_791_079_026_149_4e-7,
    2.286_212_103_119_452e-5,
    1.264_615_411_446_926e-3,
    3.597_993_651_536_15e-2,
    3.442_898_999_246_285e-1,
    -5.353_273_932_339_028e-1,
];

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */
const K0_B: [f64; 25] = [
    5.300_433_772_686_263e-18,
    -1.647_580_430_152_421_2e-17,
    5.210_391_505_039_027_4e-17,
    -1.678_231_096_805_412e-16,
    5.512_055_978_524_319e-16,
    -1.848_593_377_343_779e-15,
    6.340_076_477_405_070_6e-15,
    -2.227_513_326_991_669_8e-14,
    8.032_890_775_363_575e-14,
    -2.980_096_923_172_730_3e-13,
    1.140_340_588_208_475e-12,
    -4.514_597_883_373_944e-12,
    1.855_949_114_954_717_7e-11,
    -7.957_489_244_477_107e-11,
    3.577_397_281_400_301e-10,
    -1.697_534_509_389_06e-9,
    8.574_034_017_414_225e-9,
    -4.660_489_897_687_948e-8,
    2.766_813_639_445_015e-7,
    -1.831_755_522_719_119_5e-6,
    1.394_981_371_887_65e-5,
    -1.284_954_958_162_780_2e-4,
    1.569_883_885_730_053_3e-3,
    -3.144_810_131_196_45e-2,
    2.440_303_082_065_955_5,
];

/// Modified Bessel function, third kind, order zero
///
/// Returns modified Bessel function of the third kind
/// of order zero of the argument.
///
/// The range is partitioned into the two intervals \[0,2\] and
/// (2, infinity).  Chebyshev polynomial expansions are employed
/// in each interval.
///
/// ACCURACY:
///
/// Tested at 20000 random points between 0 and 8.  Peak absolute
/// error (relative when K0 > 1) was 9.1e-16; rms, 9.3e-17.
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       1.3e-15     1.5e-16
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
///  K0 singular        x == 0          INFINITY
///  K0 domain          x < 0           NaN
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn k0(x: f64) -> f64 {
    if x == 0.0 {
        return f64::INFINITY;
    } else if x < 0.0 {
        return f64::NAN;
    }

    if x <= 2.0 {
        let y = x * x - 2.0;
        return chbevl(y, &K0_A) - (0.5 * x).ln() * i0(x);
    }
    let z = 8.0 / x - 2.0;
    (-x).exp() * chbevl(z, &K0_B) / x.sqrt()
}

/// Modified Bessel function, third kind, order zero,
/// exponentially scaled
///
/// Returns exponentially scaled modified Bessel function
/// of the third kind of order zero of the argument.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       1.4e-15     1.4e-16
/// ```
///
/// See k0().
#[allow(clippy::float_cmp)]
pub(crate) fn k0e(x: f64) -> f64 {
    if x == 0.0 {
        return f64::INFINITY;
    } else if x < 0.0 {
        return f64::NAN;
    }

    if x <= 2.0 {
        let y = x * x - 2.0;
        let y = chbevl(y, &K0_A) - (0.5 * x).ln() * i0(x);
        return y * x.exp();
    }

    chbevl(8.0 / x - 2.0, &K0_B) / x.sqrt()
}
