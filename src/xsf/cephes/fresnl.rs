//! Translated into pure Rust from `xsf/cephes/fresnl.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.1:  January, 1989
//! Copyright 1984, 1987, 1989 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Unlike in xsf, `(S, C)` is returned instead of written to output pointers, the (unused) status
//! code is dropped, and the threshold between the two methods is corrected in the comments.

use crate::xsf::cephes::polevl::{p1evl, polevl};
use crate::xsf::trig::{cephes_cospi as cospi, cephes_sinpi as sinpi};
use core::f64::consts::PI;

/* S(x) for small x */
const FRESNL_SN: [f64; 6] = [
    -2.991_819_194_010_198_3e3,
    7.088_400_452_577_386e5,
    -6.297_414_862_058_625e7,
    2.548_908_805_733_763_7e9,
    -4.429_795_180_596_978e10,
    3.180_162_978_765_678e11,
];

const FRESNL_SD: [f64; 6] = [
    /* 1.00000000000000000000E0, */
    2.813_762_688_899_943e2,
    4.558_478_108_065_326e4,
    5.173_438_887_700_964e6,
    4.193_202_458_981_112e8,
    2.244_117_956_453_409_2e10,
    6.073_663_894_900_846e11,
];

/* C(x) for small x */
const FRESNL_CN: [f64; 6] = [
    -4.988_431_145_735_735_4e-8,
    9.504_280_628_298_596e-6,
    -6.451_914_356_839_651e-4,
    1.888_433_193_967_038_4e-2,
    -2.055_259_009_550_138_8e-1,
    1.0,
];

const FRESNL_CD: [f64; 7] = [
    3.999_829_689_724_96e-12,
    9.154_392_157_746_574e-10,
    1.250_018_624_795_988_2e-7,
    1.222_627_890_241_790_2e-5,
    8.680_295_429_417_843e-4,
    4.121_420_907_221_998e-2,
    1.0,
];

/* Auxiliary function f(x) */
const FRESNL_FN: [f64; 10] = [
    4.215_435_550_436_775e-1,
    1.434_079_197_807_589e-1,
    1.152_209_550_735_857_7e-2,
    3.450_179_397_825_74e-4,
    4.636_137_492_878_673_5e-6,
    3.055_689_837_902_576e-8,
    1.023_045_141_649_072_4e-10,
    1.720_107_432_681_618_3e-13,
    1.342_832_762_330_627_5e-16,
    3.763_297_112_699_879e-20,
];

const FRESNL_FD: [f64; 10] = [
    /*  1.00000000000000000000E0, */
    7.515_863_983_533_789e-1,
    1.168_889_258_591_913_8e-1,
    6.440_515_265_088_586_5e-3,
    1.559_344_091_641_530_1e-4,
    1.846_275_673_489_305_5e-6,
    1.126_992_247_639_990_3e-8,
    3.601_400_295_893_713_6e-11,
    5.887_545_336_215_784e-14,
    4.520_014_340_741_297_3e-17,
    1.254_432_370_900_112_7e-20,
];

/* Auxiliary function g(x) */
const FRESNL_GN: [f64; 11] = [
    5.044_420_736_433_832e-1,
    1.971_028_335_255_234e-1,
    1.876_485_840_925_752_6e-2,
    6.840_793_809_153_931e-4,
    1.151_388_261_118_842_8e-5,
    9.828_524_436_884_223e-8,
    4.453_444_158_617_501_5e-10,
    1.082_680_411_390_208_8e-12,
    1.375_554_606_332_618e-15,
    8.363_544_356_306_774e-19,
    1.869_587_101_627_832_4e-22,
];

const FRESNL_GD: [f64; 11] = [
    /*  1.00000000000000000000E0, */
    1.474_957_599_251_283_3,
    3.377_489_891_200_2e-1,
    2.536_037_414_203_388e-2,
    8.146_791_071_843_061e-4,
    1.275_450_756_677_291_2e-5,
    1.043_145_896_575_719_9e-7,
    4.606_807_281_465_204_3e-10,
    1.102_732_150_662_402_8e-12,
    1.387_965_312_595_788_6e-15,
    8.391_588_162_831_187e-19,
    1.869_587_101_627_832_4e-22,
];

/// Fresnel integral
///
/// Evaluates the Fresnel integrals
///
/// ```text
///           x
///           -
///          | |
/// C(x) =   |   cos(pi/2 t**2) dt,
///        | |
///         -
///          0
///
///           x
///           -
///          | |
/// S(x) =   |   sin(pi/2 t**2) dt.
///        | |
///         -
///          0
/// ```
///
/// The integrals are evaluated by a power series for x < 1.6.
/// For x >= 1.6 auxiliary functions f(x) and g(x) are employed
/// such that
///
/// ```text
/// C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
/// S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
/// ```
///
/// ACCURACY:
///
///  Relative error.
///
/// ```text
/// Arithmetic  function   domain     # trials      peak         rms
///   IEEE       S(x)      0, 10       20000       9.3e-16     2.1e-16
///   IEEE       C(x)      0, 10       20000       8.9e-16     2.1e-16
/// ```
pub(crate) fn fresnl(xxa: f64) -> (f64, f64) {
    let (mut cc, mut ss) = 'done: {
        if xxa.is_infinite() {
            break 'done (0.5, 0.5);
        }

        let x = xxa.abs();
        let x2 = x * x;
        if x2 < 2.5625 {
            let t = x2 * x2;
            let ss = x * x2 * polevl(t, &FRESNL_SN) / p1evl(t, &FRESNL_SD);
            let cc = x * polevl(t, &FRESNL_CN) / polevl(t, &FRESNL_CD);
            break 'done (cc, ss);
        }

        if x > 36974.0 {
            /*
             * http://functions.wolfram.com/GammaBetaErf/FresnelC/06/02/
             * http://functions.wolfram.com/GammaBetaErf/FresnelS/06/02/
             */
            let cc = 0.5 + 1.0 / (PI * x) * sinpi(x * x / 2.0);
            let ss = 0.5 - 1.0 / (PI * x) * cospi(x * x / 2.0);
            break 'done (cc, ss);
        }

        /*             Asymptotic power series auxiliary functions
         *             for large argument
         */
        let x2 = x * x;
        let t = PI * x2;
        let u = 1.0 / (t * t);
        let t = 1.0 / t;
        let f = 1.0 - u * polevl(u, &FRESNL_FN) / p1evl(u, &FRESNL_FD);
        let g = t * polevl(u, &FRESNL_GN) / p1evl(u, &FRESNL_GD);

        let c = cospi(x2 / 2.0);
        let s = sinpi(x2 / 2.0);
        let t = PI * x;
        let cc = 0.5 + (f * s - g * c) / t;
        let ss = 0.5 - (f * c + g * s) / t;
        (cc, ss)
    };

    // done:
    if xxa < 0.0 {
        cc = -cc;
        ss = -ss;
    }

    (ss, cc)
}
