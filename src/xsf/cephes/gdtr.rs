//! Translated into pure Rust from `xsf/cephes/gdtr.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.3:  March,1995
//! Copyright 1984, 1987, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the values returned on domain errors are documented as NaN instead of 0.0.

use crate::xsf::cephes::igam::{igam, igamc};
use crate::xsf::cephes::igami::igamci;

/// Gamma distribution function
///
/// Returns the integral from zero to x of the Gamma probability
/// density function:
///
/// ```text
///                x
///        b       -
///       a       | |   b-1  -at
/// y =  -----    |    t    e    dt
///       -     | |
///      | (b)   -
///               0
/// ```
///
///  The incomplete Gamma integral is used, according to the
/// relation
///
/// ```text
/// y = igam( b, ax ).
/// ```
///
/// ACCURACY:
///
/// See igam().
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// gdtr domain         x < 0            NaN
/// ```
pub(crate) fn gdtr(a: f64, b: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    igam(b, a * x)
}

/// Complemented Gamma distribution function
///
/// Returns the integral from x to infinity of the Gamma
/// probability density function:
///
/// ```text
///               inf.
///        b       -
///       a       | |   b-1  -at
/// y =  -----    |    t    e    dt
///       -     | |
///      | (b)   -
///               x
/// ```
///
///  The incomplete Gamma integral is used, according to the
/// relation
///
/// ```text
/// y = igamc( b, ax ).
/// ```
///
/// ACCURACY:
///
/// See igamc().
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// gdtrc domain         x < 0            NaN
/// ```
pub(crate) fn gdtrc(a: f64, b: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    igamc(b, a * x)
}

#[allow(dead_code)] // unused, also in xsf
pub(crate) fn gdtri(a: f64, b: f64, y: f64) -> f64 {
    if (y < 0.0) || (y > 1.0) || (a <= 0.0) || (b < 0.0) {
        return f64::NAN;
    }

    igamci(b, 1.0 - y) / a
}
