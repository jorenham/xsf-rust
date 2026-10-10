//! Translated into pure Rust from `xsf/cephes/ndtri.h` (xsf v0.2.2).
//!
//! Translated into C++ by SciPy developers in 2024 from the Cephes Math Library Release 2.1
//! (January, 1989), Copyright 1984, 1987, 1989 by Stephen L. Moshier.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::cephes::polevl::{p1evl, polevl};

/// sqrt(2*pi), from `cephes/const.h`
const SQRT2PI: f64 = 2.506_628_274_631_000_7;

/* approximation for 0 <= |y - 0.5| <= 3/8 */
const NDTRI_P0: [f64; 5] = [
    -59.963_350_101_410_79,
    98.001_075_418_599_97,
    -56.676_285_746_907_03,
    13.931_260_938_727_968,
    -1.239_165_838_673_812_5,
];
const NDTRI_Q0: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    1.954_488_583_381_417_6,
    4.676_279_128_988_815,
    86.360_242_139_089_05,
    -225.462_687_854_119_37,
    200.260_212_380_060_66,
    -82.037_225_616_833_34,
    15.905_622_512_621_17,
    -1.183_316_211_213_3,
];

/* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
 * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
 */
const NDTRI_P1: [f64; 9] = [
    4.055_448_923_059_624_5,
    31.525_109_459_989_388,
    57.162_819_224_642_13,
    44.080_507_389_320_08,
    14.684_956_192_885_803,
    2.186_633_068_507_902_5,
    -0.140_256_079_171_354_5,
    -3.504_246_268_278_482e-2,
    -8.574_567_851_546_854e-4,
];
const NDTRI_Q1: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    15.779_988_325_646_675,
    45.390_763_512_887_92,
    41.317_203_825_467_2,
    15.042_538_569_290_75,
    2.504_649_462_083_094,
    -0.142_182_922_854_787_79,
    -3.808_064_076_915_783e-2,
    -9.332_594_808_954_574e-4,
];

/* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
 * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
 */

const NDTRI_P2: [f64; 9] = [
    3.237_748_917_769_460_3,
    6.915_228_890_689_842,
    3.938_810_252_924_744_4,
    1.333_034_608_158_075_5,
    0.201_485_389_549_179_08,
    1.237_166_348_178_200_3e-2,
    3.015_815_535_082_354_3e-4,
    2.658_069_746_867_375_5e-6,
    6.239_745_391_849_833e-9,
];
const NDTRI_Q2: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    6.024_270_393_647_42,
    3.679_835_638_561_608_7,
    1.377_020_994_890_813_2,
    0.216_236_993_594_496_63,
    1.342_040_060_885_431_8e-2,
    3.280_144_646_821_277_4e-4,
    2.892_478_647_453_806_8e-6,
    6.790_194_080_099_813e-9,
];

/// `ndtri(0.5 + y)` for `|y| <= 0.5 - exp(-2)`, i.e. the central branch of [`ndtri`]
///
/// Unlike in xsf, this is a separate function, so that `erfinv` can avoid the rounding error in
/// `0.5 * (y + 1)` for small `|y|`.
pub(super) fn ndtri_central(y: f64) -> f64 {
    let y2 = y * y;
    let x = y + y * (y2 * polevl(y2, &NDTRI_P0) / p1evl(y2, &NDTRI_Q0));
    x * SQRT2PI
}

/// Inverse of Normal distribution function
///
/// Returns the argument, x, for which the area under the
/// Gaussian probability density function (integrated from
/// minus infinity to x) is equal to y.
///
/// For small arguments 0 < y < exp(-2), the program computes
/// z = sqrt( -2.0 * log(y) );  then the approximation is
/// x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
/// There are two rational functions P/Q, one for 0 < y < exp(-32)
/// and the other for y up to exp(-2).  For larger arguments,
/// w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain        # trials      peak         rms
///    IEEE     0.125, 1        20000       6.7e-16     1.2e-16
///    IEEE     3e-308, 0.135   50000       3.1e-16     8.5e-17
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition    value returned
/// ndtri domain       x < 0        NAN
/// ndtri domain       x > 1        NAN
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn ndtri(y0: f64) -> f64 {
    if y0 == 0.0 {
        return f64::NEG_INFINITY;
    }
    if y0 == 1.0 {
        return f64::INFINITY;
    }
    if y0 < 0.0 || y0 > 1.0 {
        return f64::NAN;
    }
    let mut code = 1;
    let mut y = y0;
    if y > (1.0 - 0.135_335_283_236_612_7) {
        /* 0.135... = exp(-2) */
        y = 1.0 - y;
        code = 0;
    }

    if y > 0.135_335_283_236_612_7 {
        return ndtri_central(y - 0.5);
    }

    let mut x = (-2.0 * y.ln()).sqrt();
    let x0 = x - x.ln() / x;

    let z = 1.0 / x;
    let x1 = if x < 8.0 {
        /* y > exp(-32) = 1.2664165549e-14 */
        z * polevl(z, &NDTRI_P1) / p1evl(z, &NDTRI_Q1)
    } else {
        z * polevl(z, &NDTRI_P2) / p1evl(z, &NDTRI_Q2)
    };
    x = x0 - x1;
    if code != 0 {
        x = -x;
    }
    x
}
