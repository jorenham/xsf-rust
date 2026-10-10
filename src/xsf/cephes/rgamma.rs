//! Translated into pure Rust from `xsf/cephes/rgamma.h` (xsf v0.2.2).
//!
//! Cephes Math Library Release 2.0:  April, 1987,
//! Copyright 1985, 1987 by Stephen L. Moshier.

use crate::xsf::cephes::chbevl::chbevl;
use crate::xsf::cephes::gamma;

/* Chebyshev coefficients for reciprocal Gamma function
 * in interval 0 to 1.  Function is 1/(x Gamma(x)) - 1
 */
const RGAMMA_R: [f64; 16] = [
    3.131_734_582_312_3e-17,
    -6.707_186_064_779_08e-16,
    2.200_390_781_722_595_4e-15,
    2.476_916_303_482_541_4e-13,
    -6.600_741_004_112_952e-12,
    5.138_501_863_242_27e-11,
    1.089_653_864_544_186_7e-9,
    -3.339_646_306_868_369_4e-8,
    2.689_759_964_405_954_6e-7,
    2.960_011_775_188_017e-6,
    -8.048_141_249_784_711e-5,
    4.166_091_387_096_889e-4,
    5.065_798_640_286_087e-3,
    -6.419_254_361_091_582e-2,
    -4.985_587_286_840_036e-3,
    0.127_546_015_610_523_95,
];

/// Reciprocal Gamma function
///
/// Returns one divided by the Gamma function of the argument.
///
/// The function is approximated by a Chebyshev expansion in
/// the interval \[0,1\].  Range reduction is by recurrence
/// for arguments between -4 and +4.  Outside this range,
/// 1 / Gamma(x) is returned.  (The original Cephes version
/// used the recurrence between -34.034 and +34.84425627277176174,
/// and the cosecant reflection formula below -34.034.)
///
/// The reciprocal Gamma function has no singularities,
/// but overflow and underflow may occur for large arguments.
/// These conditions return either INFINITY or 0 with
/// appropriate sign.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE     -30,+30      30000       9.7e-16     2.1e-16
/// ```
///
/// For arguments less than -34.034 the peak error is 6.9e-16
/// (IEEE, -170,-34.034, 10000 trials), excepting overflow or underflow.
#[allow(clippy::float_cmp)]
pub(crate) fn rgamma(x: f64) -> f64 {
    if x == 0.0 {
        // This case is separate from below to get correct sign for zero.
        return x;
    }

    if x < 0.0 && x == x.floor() {
        // Gamma poles.
        return 0.0;
    }

    if x.abs() > 4.0 {
        return 1.0 / gamma(x);
    }

    let mut z = 1.0;
    let mut w = x;

    while w > 1.0 {
        /* Downward recurrence */
        w -= 1.0;
        z *= w;
    }
    while w < 0.0 {
        /* Upward recurrence */
        z /= w;
        w += 1.0;
    }
    if w == 0.0 {
        /* Nonpositive integer */
        return 0.0;
    }
    if w == 1.0 {
        /* Other integer */
        return 1.0 / z;
    }

    w * (1.0 + chbevl(4.0 * w - 2.0, &RGAMMA_R)) / z
}
