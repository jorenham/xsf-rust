//! Translated into pure Rust from `xsf/cephes/ndtr.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.2:  June, 1992
//! Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use core::f64::consts::FRAC_1_SQRT_2;

use crate::utils::MAXLOG;
use crate::xsf::cephes::polevl::{p1evl, polevl};

const NDTR_P: [f64; 9] = [
    2.461_969_814_735_305e-10,
    0.564_189_564_831_068_9,
    7.463_210_564_422_699,
    48.637_197_098_568_14,
    196.520_832_956_077_1,
    526.445_194_995_477_3,
    934.528_527_171_957_6,
    1_027.551_886_895_157_2,
    557.535_335_369_399_4,
];

const NDTR_Q: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    13.228_195_115_474_499,
    86.707_214_088_598_97,
    354.937_778_887_819_9,
    975.708_501_743_205_5,
    1_823.909_166_879_097_3,
    2_246.337_608_187_109_7,
    1_656.663_091_941_613_4,
    557.535_340_817_727_7,
];

const NDTR_R: [f64; 6] = [
    0.564_189_583_547_755_1,
    1.275_366_707_599_781,
    5.019_050_422_511_805,
    6.160_210_979_930_536,
    7.409_742_699_504_489_5,
    2.978_866_653_721_002_2,
];

const NDTR_S: [f64; 6] = [
    /* 1.00000000000000000000E0, */
    2.260_528_632_201_172_6,
    9.396_035_249_380_015,
    12.048_953_980_809_666,
    17.081_445_074_756_59,
    9.608_968_090_632_859,
    3.369_076_451_000_815,
];

const NDTR_T: [f64; 5] = [
    9.604_973_739_870_516,
    90.026_019_720_384_27,
    2_232.005_345_946_843,
    7_003.325_141_128_051,
    55_592.301_301_039_49,
];

const NDTR_U: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    33.561_714_164_750_31,
    521.357_949_780_152_7,
    4_594.323_829_709_801,
    22_629.000_061_389_095,
    49_267.394_260_863_59,
];

/// Complementary error function
///
/// ```text
///  1 - erf(x) =
///
///                           inf.
///                             -
///                  2         | |          2
///   erfc(x)  =  --------     |    exp( - t  ) dt
///               sqrt(pi)   | |
///                           -
///                            x
/// ```
///
/// For small x, erfc(x) = 1 - erf(x); otherwise rational
/// approximations are computed.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0,26.6417   30000       5.7e-14     1.6e-14
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
pub(crate) fn erfc(a: f64) -> f64 {
    if a.is_nan() {
        return f64::NAN;
    }

    let x = if a < 0.0 { -a } else { a };

    if x < 1.0 {
        return 1.0 - erf(a);
    }

    let mut z = -a * a;

    if z >= -MAXLOG {
        z = z.exp();

        let (p, q) = if x < 8.0 {
            (polevl(x, &NDTR_P), p1evl(x, &NDTR_Q))
        } else {
            (polevl(x, &NDTR_R), p1evl(x, &NDTR_S))
        };
        let mut y = (z * p) / q;

        if a < 0.0 {
            y = 2.0 - y;
        }

        if y != 0.0 {
            return y;
        }
    }

    // under:
    if a < 0.0 { 2.0 } else { 0.0 }
}

/// Error function
///
/// The integral is
///
/// ```text
///                           x
///                            -
///                 2         | |          2
///   erf(x)  =  --------     |    exp( - t  ) dt.
///              sqrt(pi)   | |
///                          -
///                           0
/// ```
///
/// For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
/// erf(x) = 1 - erfc(x).
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0,1         30000       3.6e-16     1.0e-16
/// ```
pub(crate) fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }

    if x < 0.0 {
        return -erf(-x);
    }

    if x.abs() > 1.0 {
        return 1.0 - erfc(x);
    }
    let z = x * x;

    x * polevl(z, &NDTR_T) / p1evl(z, &NDTR_U)
}

/// Normal distribution function
///
/// Returns the area under the Gaussian probability density
/// function, integrated from minus infinity to x:
///
/// ```text
///                            x
///                             -
///                   1        | |          2
///    ndtr(x)  = ---------    |    exp( - t /2 ) dt
///               sqrt(2pi)  | |
///                           -
///                          -inf.
///
///             =  ( 1 + erf(z) ) / 2
///             =  erfc(z) / 2
/// ```
///
/// where z = x/sqrt(2). Computation is via the functions
/// erf and erfc.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE     -13,0        30000       3.4e-14     6.8e-15
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition         value returned
/// erfc underflow    x > 37.519379347       0.0
/// ```
pub(crate) fn ndtr(a: f64) -> f64 {
    if a.is_nan() {
        return f64::NAN;
    }

    let x = a * FRAC_1_SQRT_2;
    let z = x.abs();

    if z < 1.0 {
        0.5 + 0.5 * erf(x)
    } else {
        let y = 0.5 * erfc(z);
        if x > 0.0 { 1.0 - y } else { y }
    }
}
