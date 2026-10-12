//! Translated into pure Rust from `xsf/cephes/ellpk.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library, Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the domain error condition and return value are corrected in the comments.

use crate::xsf::cephes::polevl::polevl;

const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53

const ELLPK_P: [f64; 11] = [
    1.379_828_646_062_732_5e-4,
    2.280_257_240_058_756e-3,
    7.974_040_132_204_152e-3,
    9.858_213_790_212_26e-3,
    6.874_896_874_499_499e-3,
    6.189_010_336_376_876e-3,
    8.790_782_739_527_438e-3,
    1.493_804_489_168_052_6e-2,
    3.088_514_652_467_12e-2,
    9.657_359_028_116_902e-2,
    1.386_294_361_119_890_6,
];

const ELLPK_Q: [f64; 11] = [
    2.940_789_550_485_985e-5,
    9.141_847_238_659_173e-4,
    5.940_583_037_531_678e-3,
    1.548_505_166_497_624e-2,
    2.390_896_027_159_248_8e-2,
    3.012_047_152_276_040_4e-2,
    3.737_743_141_738_232_6e-2,
    4.882_803_475_709_983e-2,
    7.031_249_969_639_575e-2,
    1.249_999_999_998_708_3e-1,
    0.5,
];

const ELLPK_C1: f64 = 1.386_294_361_119_890_6; /* log(4) */

/// Complete elliptic integral of the first kind
///
/// Approximates the integral
///
/// ```text
///            pi/2
///             -
///            | |
///            |           dt
/// K(m)  =    |    ------------------
///            |                   2
///          | |    sqrt( 1 - m sin t )
///           -
///            0
/// ```
///
/// where m = 1 - m1, using the approximation
///
/// ```text
///     P(x)  -  log x Q(x).
/// ```
///
/// The argument m1 is used internally rather than m so that the logarithmic
/// singularity at m = 1 will be shifted to the origin; this
/// preserves maximum accuracy.
///
/// K(0) = pi/2.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE       0,1        20000       2.6e-16     6.8e-17
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// ellpk domain         x<0               NaN
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn ellpk(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    if x > 1.0 {
        if x.is_infinite() {
            return 0.0;
        }
        return ellpk(1.0 / x) / x.sqrt();
    }

    if x > MACHEP {
        polevl(x, &ELLPK_P) - x.ln() * polevl(x, &ELLPK_Q)
    } else if x == 0.0 {
        f64::INFINITY
    } else {
        ELLPK_C1 - 0.5 * x.ln()
    }
}
