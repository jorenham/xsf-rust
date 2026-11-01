//! Translated into pure Rust from `xsf/cephes/j0.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
//!
//! Note: all coefficients satisfy the relative error criterion
//! except YP, YQ which are designed for absolute error.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the sign of the value of `y0(x) - 2 log(x) j0(x) / pi` at `x = 0`, and the
//! absolute error condition in the `y0` accuracy table, are corrected in the comments.

use crate::utils::{SQRT1OPI, SQRT2OPI};
use crate::xsf::cephes::polevl::{p1evl, polevl};
use core::f64::consts::{FRAC_2_PI, FRAC_PI_4};

const J0_PP: [f64; 7] = [
    7.969_367_292_973_471e-4,
    8.283_523_921_074_408e-2,
    1.239_533_716_464_143,
    5.447_250_030_587_687,
    8.747_165_001_998_17,
    5.303_240_382_353_949,
    1.0,
];

const J0_PQ: [f64; 7] = [
    9.244_088_105_588_637e-4,
    8.562_884_743_544_745e-2,
    1.253_527_439_010_589_5,
    5.470_977_403_304_171,
    8.761_908_832_370_695,
    5.306_052_882_353_947,
    1.0,
];

const J0_QP: [f64; 8] = [
    -1.136_638_388_984_691_6e-2,
    -1.282_527_186_705_093_1,
    -1.955_395_442_577_359_7e1,
    -9.320_601_521_237_683e1,
    -1.776_811_679_804_880_6e2,
    -1.470_775_051_549_511_8e2,
    -5.141_053_267_665_993e1,
    -6.050_143_506_007_285,
];

const J0_QQ: [f64; 7] = [
    /*  1.00000000000000000000E0, */
    6.431_782_561_181_78e1,
    8.564_300_259_769_806e2,
    3.882_401_836_054_016_3e3,
    7.240_467_741_956_525e3,
    5.930_727_011_873_169e3,
    2.062_093_316_603_278_3e3,
    2.420_057_402_402_914e2,
];

const J0_YP: [f64; 8] = [
    1.559_243_678_552_357_4e4,
    -1.466_392_959_039_716e7,
    5.435_264_770_518_765e9,
    -9.821_360_657_179_115e11,
    8.759_063_943_953_67e13,
    -3.466_283_033_847_297e15,
    4.427_332_685_725_698_4e16,
    -1.849_508_004_369_866_8e16,
];

const J0_YQ: [f64; 7] = [
    /*  1.00000000000000000000E0, */
    1.041_283_536_642_598_4e3,
    6.261_073_301_371_35e5,
    2.689_196_333_938_141_5e8,
    8.640_024_871_039_35e10,
    2.029_796_127_501_055_5e13,
    3.171_577_528_429_750_5e15,
    2.505_962_561_726_530_6e17,
];

/*  5.783185962946784521175995758455807035071 */
const J0_DR1: f64 = 5.783_185_962_946_784;

/* 30.47126234366208639907816317502275584842 */
const J0_DR2: f64 = 30.471_262_343_662_087;

const J0_RP: [f64; 4] = [
    -4.794_432_209_782_018e9,
    1.956_174_919_465_565_7e12,
    -2.492_483_443_609_677_2e14,
    9.708_622_510_473_064e15,
];

const J0_RQ: [f64; 8] = [
    /*  1.00000000000000000000E0, */
    4.995_631_471_526_51e2,
    1.737_854_016_763_747e5,
    4.844_096_583_399_621e7,
    1.118_555_370_453_568_3e10,
    2.112_775_201_154_892e12,
    3.105_182_298_574_225_6e14,
    3.181_219_559_432_049_6e16,
    1.710_862_940_810_431_5e18,
];

/// Bessel function of order zero
///
/// Returns Bessel function of order zero of the argument.
///
/// The domain is divided into the intervals \[0, 5\] and
/// (5, infinity). In the first interval the following rational
/// approximation is used:
///
/// ```text
///        2         2
/// (w - r  ) (w - r  ) P (w) / Q (w)
///       1         2    3       8
///
///            2
/// where w = x  and the two r's are zeros of the function.
/// ```
///
/// In the second interval, the Hankel asymptotic expansion
/// is employed with two rational functions of degree 6/6
/// and 7/7.
///
/// ACCURACY:
///
/// ```text
///                      Absolute error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       3.9e-16     5.2e-17
/// ```
pub(crate) fn j0(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }

    if x <= 5.0 {
        let z = x * x;
        if x < 1.0e-5 {
            return 1.0 - z / 4.0;
        }

        let p = (z - J0_DR1) * (z - J0_DR2);
        return p * polevl(z, &J0_RP) / p1evl(z, &J0_RQ);
    }

    let w = 5.0 / x;
    let q = 25.0 / (x * x);
    let p = polevl(q, &J0_PP) / polevl(q, &J0_PQ);
    let q = polevl(q, &J0_QP) / p1evl(q, &J0_QQ);
    if x < 10.0 {
        let xn = x - FRAC_PI_4;
        let p = p * xn.cos() - w * q * xn.sin();
        return p * SQRT2OPI / x.sqrt();
    }
    let p = (p + w * q) * x.cos() + (p - w * q) * x.sin();
    p * SQRT1OPI / x.sqrt()
}

/* Bessel function of second kind, order zero  */

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  y0(x)  -  2 * log(x) * j0(x) / M_PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / M_PI
 * = -0.073804295108687225.
 */

/// Bessel function of the second kind, order zero
///
/// Returns Bessel function of the second kind, of order
/// zero, of the argument.
///
/// The domain is divided into the intervals \[0, 5\] and
/// (5, infinity). In the first interval a rational approximation
/// R(x) is employed to compute
///   y0(x)  = R(x)  +   2 * log(x) * j0(x) / M_PI.
/// Thus a call to j0() is required.
///
/// In the second interval, the Hankel asymptotic expansion
/// is employed with two rational functions of degree 6/6
/// and 7/7.
///
/// ACCURACY:
///
/// ```text
///  Absolute error, when |y0(x)| < 1; else relative error:
///
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       1.3e-15     1.3e-16
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn y0(x: f64) -> f64 {
    if x <= 5.0 {
        if x == 0.0 {
            return f64::NEG_INFINITY;
        } else if x < 0.0 {
            return f64::NAN;
        }
        let z = x * x;
        let mut w = polevl(z, &J0_YP) / p1evl(z, &J0_YQ);
        w += FRAC_2_PI * x.ln() * j0(x);
        return w;
    }

    let w = 5.0 / x;
    let z = 25.0 / (x * x);
    let p = polevl(z, &J0_PP) / polevl(z, &J0_PQ);
    let q = polevl(z, &J0_QP) / p1evl(z, &J0_QQ);
    if x < 10.0 {
        let xn = x - FRAC_PI_4;
        let p = p * xn.sin() + w * q * xn.cos();
        return p * SQRT2OPI / x.sqrt();
    }
    let p = (w * q - p) * x.cos() + (p + w * q) * x.sin();
    p * SQRT1OPI / x.sqrt()
}
