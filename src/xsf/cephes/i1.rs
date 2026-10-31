//! Translated into pure Rust from `xsf/cephes/i1.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1985, 1987, 2000 by Stephen L. Moshier

use crate::utils::MAXLOG;
use crate::xsf::cephes::chbevl::chbevl;

/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */
const I1_A: [f64; 29] = [
    2.777_914_112_761_046_4e-18,
    -2.111_421_214_358_166e-17,
    1.553_631_957_736_200_5e-16,
    -1.105_596_947_735_386_2e-15,
    7.600_684_294_735_408e-15,
    -5.042_185_504_727_912e-14,
    3.223_793_365_945_575e-13,
    -1.983_974_397_764_943_6e-12,
    1.173_618_629_889_090_1e-11,
    -6.663_489_723_502_027e-11,
    3.625_590_281_552_117e-10,
    -1.887_249_751_722_829_4e-9,
    9.381_537_386_495_773e-9,
    -4.445_059_128_796_328e-8,
    2.003_294_753_552_135_3e-7,
    -8.568_720_264_695_455e-7,
    3.470_251_308_137_678_5e-6,
    -1.327_316_365_603_943_6e-5,
    4.781_565_107_550_054e-5,
    -1.617_608_158_258_967_4e-4,
    5.122_859_561_685_758e-4,
    -1.513_572_450_631_253_2e-3,
    4.156_422_944_312_888e-3,
    -1.056_408_489_462_619_7e-2,
    2.472_644_903_062_651_6e-2,
    -5.294_598_120_809_499e-2,
    1.026_436_586_898_471e-1,
    -1.764_165_183_578_340_6e-1,
    2.525_871_864_436_336_5e-1,
];

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */
const I1_B: [f64; 25] = [
    7.517_296_310_842_105e-18,
    4.414_348_323_071_708e-18,
    -4.650_305_368_489_358_6e-17,
    -3.209_525_921_993_424e-17,
    2.962_628_997_645_95e-16,
    3.308_202_310_920_928_5e-16,
    -1.880_354_775_510_782_5e-15,
    -3.814_403_072_437_007_5e-15,
    1.042_027_698_412_880_2e-14,
    4.272_440_016_711_951e-14,
    -2.101_541_842_772_664_3e-14,
    -4.083_551_111_092_197_4e-13,
    -7.198_551_776_245_908e-13,
    2.035_628_544_147_089_6e-12,
    1.412_580_743_661_378_2e-11,
    3.252_603_583_015_488_4e-11,
    -1.897_495_812_350_541_3e-11,
    -5.589_743_462_196_584e-10,
    -3.835_380_385_964_237e-9,
    -2.631_468_846_889_519_6e-8,
    -2.512_236_237_870_209e-7,
    -3.882_564_808_877_691e-6,
    -1.105_889_387_626_237_1e-4,
    -9.761_097_491_361_469e-3,
    7.785_762_350_182_801e-1,
];

/// Modified Bessel function of order one
///
/// Returns modified Bessel function of order one of the
/// argument.
///
/// The function is defined as i1(x) = -i j1( ix ).
///
/// The range is partitioned into the two intervals \[0,8\] and
/// (8, infinity).  Chebyshev polynomial expansions are employed
/// in each interval.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       1.9e-15     2.1e-16
/// ```
pub(crate) fn i1(x: f64) -> f64 {
    let mut z = x.abs();
    if z <= 8.0 {
        let y = (z / 2.0) - 2.0;
        z = chbevl(y, &I1_A) * z * z.exp();
    } else if z > MAXLOG {
        /* exp(z) overflows here even though i1(z) is still finite up to z ~ 713.99,
         * so evaluate the exponential in two halves. */
        let e = (z / 2.0).exp();
        z = e * chbevl(32.0 / z - 2.0, &I1_B) / z.sqrt() * e;
    } else {
        z = z.exp() * chbevl(32.0 / z - 2.0, &I1_B) / z.sqrt();
    }
    if x < 0.0 {
        z = -z;
    }
    z
}

/// Modified Bessel function of order one,
/// exponentially scaled
///
/// Returns exponentially scaled modified Bessel function
/// of order one of the argument.
///
/// The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       1.8e-15     2.0e-16
/// ```
///
/// See i1().
pub(crate) fn i1e(x: f64) -> f64 {
    let mut z = x.abs();
    if z <= 8.0 {
        let y = (z / 2.0) - 2.0;
        z *= chbevl(y, &I1_A);
    } else {
        z = chbevl(32.0 / z - 2.0, &I1_B) / z.sqrt();
    }
    if x < 0.0 {
        z = -z;
    }
    z
}
