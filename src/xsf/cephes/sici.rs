//! Translated into pure Rust from `xsf/cephes/sici.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.1:  January, 1989
//! Copyright 1984, 1987, 1989 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Unlike in xsf, `(Si, Ci)` is returned instead of written to output references, the (unused)
//! status code is dropped, and the threshold of the auxiliary functions is corrected in the
//! comments. The `Si` and `Ci` values that xsf computes for finite `x > 1e9` are omitted, because
//! they are always overwritten by the asymptotic expansion.

use crate::utils::SCIPY_EULER;
use crate::xsf::cephes::polevl::{p1evl, polevl};
use core::f64::consts::FRAC_PI_2;

const SICI_SN: [f64; 6] = [
    -8.391_678_279_103_039e-11,
    4.625_917_144_270_128e-8,
    -9.757_593_038_436_328e-6,
    9.769_454_381_704_354e-4,
    -4.134_703_162_294_066e-2,
    1.0,
];

const SICI_SD: [f64; 6] = [
    2.032_692_661_959_519_3e-12,
    1.279_978_911_799_433e-9,
    4.418_278_428_012_189e-7,
    9.964_121_220_438_756e-5,
    1.420_852_393_261_499e-2,
    1.0,
];

const SICI_CN: [f64; 6] = [
    2.025_240_023_891_022_8e-11,
    -1.352_495_049_157_907_6e-8,
    3.593_250_514_199_931e-6,
    -4.740_072_068_734_079e-4,
    2.891_596_526_075_552_3e-2,
    -1.0,
];

const SICI_CD: [f64; 6] = [
    4.077_460_400_618_806e-12,
    3.067_809_975_818_878e-9,
    1.232_103_556_858_834_2e-6,
    3.174_420_247_750_327_5e-4,
    5.100_280_562_364_460_6e-2,
    4.0,
];

const SICI_FN4: [f64; 7] = [
    4.236_128_628_922_166,
    5.459_377_171_618_128_5,
    1.620_832_877_015_383_3,
    1.670_066_118_313_230_4e-1,
    6.810_201_324_725_182e-3,
    1.089_365_806_503_286_7e-4,
    5.489_002_234_213_736e-7,
];

const SICI_FD4: [f64; 7] = [
    /*  1.00000000000000000000E0, */
    8.164_966_342_053_91,
    7.308_288_225_055_645,
    1.867_922_579_501_842,
    1.787_920_529_631_499e-1,
    7.017_106_683_227_897e-3,
    1.100_343_571_539_157_3e-4,
    5.489_002_527_562_557e-7,
];

const SICI_FN8: [f64; 9] = [
    4.558_808_734_704_653e-1,
    7.137_152_741_001_467e-1,
    1.603_001_582_223_194_7e-1,
    1.160_642_294_081_244e-2,
    3.495_564_424_478_590_6e-4,
    4.862_154_308_264_547_5e-6,
    3.200_927_900_910_049e-8,
    9.417_795_761_285_13e-11,
    9.705_071_108_819_52e-14,
];

const SICI_FD8: [f64; 8] = [
    /*  1.00000000000000000000E0, */
    9.174_636_118_736_841e-1,
    1.786_855_453_320_745_4e-1,
    1.222_535_947_719_712_9e-2,
    3.586_964_818_818_515_7e-4,
    4.924_350_643_178_815e-6,
    3.219_569_391_010_46e-8,
    9.437_205_903_502_767e-11,
    9.705_071_108_819_52e-14,
];

const SICI_GN4: [f64; 8] = [
    8.710_016_989_731_142e-2,
    6.113_791_099_522_193e-1,
    3.971_802_963_923_375e-1,
    7.485_277_376_284_691e-2,
    5.388_686_814_621_773e-3,
    1.619_997_945_989_340_3e-4,
    1.979_638_741_409_636_5e-6,
    7.825_790_407_440_903e-9,
];

const SICI_GD4: [f64; 7] = [
    /*  1.00000000000000000000E0, */
    1.644_022_024_133_553_5,
    6.662_967_012_689_88e-1,
    9.887_717_612_776_888e-2,
    6.223_963_454_417_684e-3,
    1.732_210_814_741_771e-4,
    2.026_591_820_863_439_7e-6,
    7.825_792_189_335_346e-9,
];

const SICI_GN8: [f64; 9] = [
    6.973_599_534_432_762e-1,
    3.304_109_793_056_320_7e-1,
    3.848_787_676_499_743e-2,
    1.717_182_390_523_479e-3,
    3.489_411_655_022_794_6e-5,
    3.471_311_670_841_167e-7,
    1.704_044_527_820_445_2e-9,
    3.859_459_254_302_766e-12,
    3.140_400_989_463_633_5e-15,
];

const SICI_GD8: [f64; 9] = [
    /*  1.00000000000000000000E0, */
    1.685_488_988_110_116_5,
    4.878_522_586_953_049_6e-1,
    4.679_131_942_596_258e-2,
    1.902_844_266_743_995_3e-3,
    3.684_755_044_425_611e-5,
    3.570_432_234_437_408_3e-7,
    1.726_937_489_663_161_5e-9,
    3.878_301_660_239_547e-12,
    3.140_400_989_463_633_5e-15,
];

/// Sine and cosine integrals
///
/// Evaluates the integrals
///
/// ```text
///                          x
///                          -
///                         |  cos t - 1
///   Ci(x) = eul + ln x +  |  --------- dt,
///                         |      t
///                        -
///                         0
///             x
///             -
///            |  sin t
///   Si(x) =  |  ----- dt
///            |    t
///           -
///            0
/// ```
///
/// where eul = 0.57721566490153286061 is Euler's constant.
/// The integrals are approximated by rational functions.
/// For x > 4 auxiliary functions f(x) and g(x) are employed
/// such that
///
/// ```text
/// Ci(x) = f(x) sin(x) - g(x) cos(x)
/// Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)
/// ```
///
/// ACCURACY:
///    Test interval = \[0,50\].
/// Absolute error, except relative when > 1:
///
/// ```text
/// arithmetic   function   # trials      peak         rms
///    IEEE        Si        20000       3.8e-16     7.3e-17
///    IEEE        Ci        20000       8.4e-16     4.6e-17
/// ```
#[allow(clippy::float_cmp, clippy::many_single_char_names)]
pub(crate) fn sici(mut x: f64) -> (f64, f64) {
    let sign: i16 = if x < 0.0 {
        x = -x;
        -1
    } else {
        0
    };

    if x == 0.0 {
        return (0.0, f64::NEG_INFINITY);
    }

    if x > 1.0e9 && x.is_infinite() {
        if sign == -1 {
            return (-FRAC_PI_2, f64::NAN);
        }
        return (FRAC_PI_2, 0.0);
    }

    'asympt: {
        if x > 4.0 {
            break 'asympt;
        }

        let z = x * x;
        let mut s = x * polevl(z, &SICI_SN) / polevl(z, &SICI_SD);
        let c = z * polevl(z, &SICI_CN) / polevl(z, &SICI_CD);

        if sign != 0 {
            s = -s;
        }
        let si = s;
        let ci = SCIPY_EULER + x.ln() + c; /* real part if x < 0 */
        return (si, ci);
    }

    /* The auxiliary functions are:
     *
     *
     * *si = *si - M_PI_2;
     * c = cos(x);
     * s = sin(x);
     *
     * t = *ci * s - *si * c;
     * a = *ci * c + *si * s;
     *
     * *si = t;
     * *ci = -a;
     */

    // asympt:

    let s = x.sin();
    let c = x.cos();
    let z = 1.0 / (x * x);
    let (f, g);
    if x < 8.0 {
        f = polevl(z, &SICI_FN4) / (x * p1evl(z, &SICI_FD4));
        g = z * polevl(z, &SICI_GN4) / p1evl(z, &SICI_GD4);
    } else {
        f = polevl(z, &SICI_FN8) / (x * p1evl(z, &SICI_FD8));
        g = z * polevl(z, &SICI_GN8) / p1evl(z, &SICI_GD8);
    }
    let mut si = FRAC_PI_2 - f * c - g * s;
    if sign != 0 {
        si = -(si);
    }
    let ci = f * s - g * c;

    (si, ci)
}
