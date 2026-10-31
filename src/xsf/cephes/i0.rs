//! Translated into pure Rust from `xsf/cephes/i0.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 2000 by Stephen L. Moshier

use crate::xsf::cephes::chbevl::chbevl;

const MAXLOG: f64 = 709.782_712_893_384; // log(DBL_MAX)

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */
const I0_A: [f64; 30] = [
    -4.415_341_646_479_339_5e-18,
    3.330_794_518_822_238_4e-17,
    -2.431_279_846_547_955e-16,
    1.715_391_285_555_133e-15,
    -1.168_533_287_799_345_1e-14,
    7.676_185_498_604_936e-14,
    -4.856_446_783_111_929e-13,
    2.955_052_663_129_64e-12,
    -1.726_826_291_441_556e-11,
    9.675_809_035_373_237e-11,
    -5.189_795_601_635_263e-10,
    2.659_823_724_682_386_6e-9,
    -1.300_025_009_986_248e-8,
    6.046_995_022_541_919e-8,
    -2.670_793_853_940_612e-7,
    1.117_387_539_120_103_7e-6,
    -4.416_738_358_458_750_5e-6,
    1.644_844_807_072_889_6e-5,
    -5.754_195_010_082_104e-5,
    1.885_028_850_958_416_5e-4,
    -5.763_755_745_385_824e-4,
    1.639_475_616_941_335_7e-3,
    -4.324_309_995_050_576e-3,
    1.054_646_039_459_499_8e-2,
    -2.373_741_480_589_947e-2,
    4.930_528_423_967_071e-2,
    -9.490_109_704_804_764e-2,
    1.716_209_015_222_087_7e-1,
    -3.046_826_723_431_984e-1,
    6.767_952_744_094_761e-1,
];

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */
const I0_B: [f64; 25] = [
    -7.233_180_487_874_754e-18,
    -4.830_504_485_944_182e-18,
    4.465_621_420_296_76e-17,
    3.461_222_867_697_461e-17,
    -2.827_623_980_516_583_6e-16,
    -3.425_485_619_677_219e-16,
    1.772_560_133_056_526_3e-15,
    3.811_680_669_352_622_4e-15,
    -9.554_846_698_828_307e-15,
    -4.150_569_347_287_222e-14,
    1.540_086_217_521_41e-14,
    3.852_778_382_742_142_6e-13,
    7.180_124_451_383_666e-13,
    -1.794_178_531_506_806_2e-12,
    -1.321_581_184_044_771_3e-11,
    -3.149_916_527_963_241_6e-11,
    1.188_914_710_784_643_9e-11,
    4.940_602_388_224_97e-10,
    3.396_232_025_708_386_5e-9,
    2.266_668_990_498_178e-8,
    2.048_918_589_469_063_8e-7,
    2.891_370_520_834_756_7e-6,
    6.889_758_346_916_825e-5,
    3.369_116_478_255_694_3e-3,
    8.044_904_110_141_088e-1,
];

/// Modified Bessel function of order zero
///
/// Returns modified Bessel function of order zero of the
/// argument.
///
/// The function is defined as i0(x) = j0( ix ).
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
///    IEEE      0,30        20000       8.0e-16     1.3e-16
/// ```
pub(crate) fn i0(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }
    if x <= 8.0 {
        let y = (x / 2.0) - 2.0;
        return x.exp() * chbevl(y, &I0_A);
    }

    if x > MAXLOG {
        /* exp(x) overflows here even though i0(x) is still finite up to x ~ 713.99,
         * so evaluate the exponential in two halves. */
        let e = (x / 2.0).exp();
        return e * chbevl(32.0 / x - 2.0, &I0_B) / x.sqrt() * e;
    }

    x.exp() * chbevl(32.0 / x - 2.0, &I0_B) / x.sqrt()
}

/// Modified Bessel function of order zero,
/// exponentially scaled
///
/// Returns exponentially scaled modified Bessel function
/// of order zero of the argument.
///
/// The function is defined as i0e(x) = exp(-|x|) j0( ix ).
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0,30        20000       7.4e-16     1.1e-16
/// ```
///
/// See i0().
pub(crate) fn i0e(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }
    if x <= 8.0 {
        let y = (x / 2.0) - 2.0;
        return chbevl(y, &I0_A);
    }

    chbevl(32.0 / x - 2.0, &I0_B) / x.sqrt()
}
