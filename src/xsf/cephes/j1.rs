//! Translated into pure Rust from `xsf/cephes/j1.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the descriptions of the approximations (intervals, kind, and degrees) are
//! corrected in the comments.

use crate::utils::{SQRT_EPS, SQRT1OPI, SQRT2OPI};
use crate::xsf::cephes::polevl::{p1evl, polevl};
use core::f64::consts::{FRAC_1_SQRT_2, FRAC_2_PI};

const THPIO4: f64 = 2.356_194_490_192_345; // 3*pi/4

const J1_RP: [f64; 4] = [
    -8.999_712_257_055_594e8,
    4.522_282_979_981_940_3e11,
    -7.274_942_452_218_183e13,
    3.682_957_328_638_529e15,
];

const J1_RQ: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    6.208_364_781_180_543e2,
    2.569_872_567_577_488_4e5,
    8.351_467_914_319_493e7,
    2.215_115_954_797_925e10,
    4.749_141_220_799_914e12,
    7.843_696_078_762_359e14,
    8.952_223_361_846_274e16,
    5.322_786_203_326_801e18,
];

const J1_PP: [f64; 7] = [
    7.621_256_162_081_731e-4,
    7.313_970_569_409_176e-2,
    1.127_196_081_296_849_3,
    5.112_079_511_468_076,
    8.424_045_901_417_724,
    5.214_515_986_823_615,
    1.0,
];

const J1_PQ: [f64; 7] = [
    5.713_231_280_725_487e-4,
    6.884_559_087_544_954e-2,
    1.105_142_326_340_617,
    5.073_863_861_286_015,
    8.399_855_543_276_042,
    5.209_828_486_823_619,
    1.0,
];

const J1_QP: [f64; 8] = [
    5.108_625_947_501_766e-2,
    4.982_138_729_512_334,
    7.582_382_841_325_453e1,
    3.667_796_093_601_508e2,
    7.108_563_049_989_261e2,
    5.974_896_124_006_136e2,
    2.116_887_571_005_721_3e2,
    2.520_702_058_580_237_2e1,
];

const J1_QQ: [f64; 7] = [
    /* 1.00000000000000000000E0, */
    7.423_732_770_356_752e1,
    1.056_448_860_382_628_3e3,
    4.986_410_583_376_536e3,
    9.562_318_924_047_562e3,
    7.997_041_604_473_507e3,
    2.826_192_785_176_390_8e3,
    3.360_936_078_106_983e2,
];

const J1_YP: [f64; 6] = [
    1.263_204_747_901_780_4e9,
    -6.473_558_763_791_603e11,
    1.145_095_115_418_237_3e14,
    -8.127_702_555_013_251e15,
    2.024_394_757_135_949e17,
    -7.788_771_962_659_501e17,
];

const J1_YQ: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    5.943_015_923_461_282e2,
    2.355_640_929_430_685_6e5,
    7.348_119_444_597_217e7,
    1.876_013_161_087_061_7e10,
    3.882_312_774_962_385_7e12,
    6.205_577_271_469_538e14,
    6.871_410_873_553_005e16,
    3.972_706_081_165_606_4e18,
];

const J1_Z1: f64 = 14.681_970_642_123_893;
const J1_Z2: f64 = 49.218_456_321_694_6;

/// Bessel function of order one
///
/// Returns Bessel function of order one of the argument.
///
/// The domain is divided into the intervals \[0, 5\] and
/// (5, infinity). In the first interval a rational
/// approximation is used. In the second, the asymptotic
/// trigonometric representation is employed using two
/// rational functions of degree 6/6 and 7/7.
///
/// ACCURACY:
///
/// ```text
///                      Absolute error:
/// arithmetic   domain      # trials      peak         rms
///    IEEE      0, 30       20000       3.0e-16     3.6e-17
/// ```
pub(crate) fn j1(x: f64) -> f64 {
    if x < 0.0 {
        return -j1(-x);
    }
    if x < SQRT_EPS {
        return 0.5 * x;
    }
    if x <= 5.0 {
        let z = x * x;
        let w = polevl(z, &J1_RP) / p1evl(z, &J1_RQ);
        return w * x * (z - J1_Z1) * (z - J1_Z2);
    }

    let w = 5.0 / x;
    let z = w * w;
    let p = polevl(z, &J1_PP) / polevl(z, &J1_PQ);
    let q = polevl(z, &J1_QP) / p1evl(z, &J1_QQ);
    if x < 10.0 {
        let xn = x - THPIO4;
        let p = p * xn.cos() - w * q * xn.sin();
        return p * SQRT2OPI / x.sqrt();
    }
    let a = FRAC_1_SQRT_2 * (w * q - p);
    let b = FRAC_1_SQRT_2 * (p + w * q);
    (a * x.cos() + b * x.sin()) * SQRT2OPI / x.sqrt()
}

/// Bessel function of second kind of order one
///
/// Returns Bessel function of the second kind of order one
/// of the argument.
///
/// The domain is divided into the intervals \[0, 5\] and
/// (5, infinity). In the first interval a rational
/// approximation is used, and a call to j1() is required.
/// In the second, the asymptotic trigonometric representation
/// is employed using two rational functions of degree 6/6 and 7/7.
///
/// ACCURACY:
///
/// ```text
///                      Absolute error:
/// arithmetic   domain      # trials      peak         rms
///    IEEE      0, 30       20000       1.0e-15     7.5e-17
///
/// (error criterion relative when |y1| > 1).
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn y1(x: f64) -> f64 {
    if x <= 5.0 {
        if x == 0.0 {
            return f64::NEG_INFINITY;
        } else if x <= 0.0 {
            return f64::NAN;
        }
        let z = x * x;
        let mut w = x * (polevl(z, &J1_YP) / p1evl(z, &J1_YQ));
        w += FRAC_2_PI * (j1(x) * x.ln() - 1.0 / x);
        return w;
    }

    let w = 5.0 / x;
    let z = w * w;
    let p = polevl(z, &J1_PP) / polevl(z, &J1_PQ);
    let q = polevl(z, &J1_QP) / p1evl(z, &J1_QQ);
    if x < 10.0 {
        let xn = x - THPIO4;
        let p = p * xn.sin() + w * q * xn.cos();
        return p * SQRT2OPI / x.sqrt();
    }
    let p = (w * q - p) * x.sin() - (p + w * q) * x.cos();
    p * SQRT1OPI / x.sqrt()
}
