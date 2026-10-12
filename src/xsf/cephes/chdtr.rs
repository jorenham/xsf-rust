//! Translated into pure Rust from `xsf/cephes/chdtr.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, some formulas and domain conditions in the comments are corrected.

use crate::xsf::cephes::igam::{igam, igamc};
use crate::xsf::cephes::igami::igamci;

/// Complemented Chi-square distribution
///
/// Returns the area under the right hand tail (from x to
/// infinity) of the Chi square probability density function
/// with v degrees of freedom:
///
/// ```text
///                                  inf.
///                                    -
///                        1          | |  v/2-1  -t/2
///  P( x | v )   =   -----------     |   t      e     dt
///                    v/2  -       | |
///                   2    | (v/2)   -
///                                   x
/// ```
///
/// where x is the Chi-square variable.
///
/// The incomplete Gamma integral is used, according to the
/// formula
///
/// ```text
///     y = chdtrc( v, x ) = igamc( v/2.0, x/2.0 ).
/// ```
///
/// The arguments must both be positive.
///
/// ACCURACY:
///
/// See igamc().
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// chdtrc domain       x < 0               NaN
/// ```
pub(crate) fn chdtrc(df: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    igamc(df / 2.0, x / 2.0)
}

/// Chi-square distribution
///
/// Returns the area under the left hand tail (from 0 to x)
/// of the Chi square probability density function with
/// v degrees of freedom.
///
/// ```text
///                                   x
///                                    -
///                        1          | |  v/2-1  -t/2
///  P( x | v )   =   -----------     |   t      e     dt
///                    v/2  -       | |
///                   2    | (v/2)   -
///                                   0
/// ```
///
/// where x is the Chi-square variable.
///
/// The incomplete Gamma integral is used, according to the
/// formula
///
/// ```text
///     y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
/// ```
///
/// The arguments must both be positive.
///
/// ACCURACY:
///
/// See igam().
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// chdtr domain        x < 0               NaN
/// ```
pub(crate) fn chdtr(df: f64, x: f64) -> f64 {
    if x < 0.0 {
        // || (df < 1.0)
        return f64::NAN;
    }
    igam(df / 2.0, x / 2.0)
}

/// Inverse of complemented Chi-square distribution
///
/// Finds the Chi-square argument x such that the integral
/// from x to infinity of the Chi-square density is equal
/// to the given cumulative probability y.
///
/// This is accomplished using the inverse Gamma integral
/// function and the relation
///
/// ```text
///    x/2 = igamci( df/2, y );
/// ```
///
/// ACCURACY:
///
/// See igami.c.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// chdtri domain   y < 0 or y > 1          NaN
/// ```
pub(crate) fn chdtri(df: f64, y: f64) -> f64 {
    if (y < 0.0) || (y > 1.0) {
        // || (df < 1.0)
        return f64::NAN;
    }

    let x = igamci(0.5 * df, y);
    2.0 * x
}
