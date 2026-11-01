//! Translated into pure Rust from `xsf/cephes/fdtr.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.3:  March, 1995
//! Copyright 1984, 1987, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, some formulas and domain conditions in the comments are corrected.

use crate::xsf::cephes::incbet::{betainc_xy, incbet};
use crate::xsf::cephes::incbi::incbi;

/// Complemented F distribution
///
/// Returns the area from x to infinity under the F density
/// function (also known as Snedcor's density or the
/// variance ratio density).
///
/// The incomplete beta integral is used, according to the
/// formula
///
/// ```text
///  1-P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).
/// ```
///
/// ACCURACY:
///
/// Tested at random points (a,b,x) in the indicated intervals.
///
/// ```text
///                x     a,b                     Relative error:
/// arithmetic  domain  domain     # trials      peak         rms
///    IEEE      0,1    1,100        20000      4.5e-14     3.5e-15
///    IEEE      1,5    1,100        20000      4.0e-14     1.3e-14
///    IEEE      0,1    1,10000      20000      2.3e-11     7.1e-13
///    IEEE      1,5    1,10000      20000      4.7e-11     6.1e-12
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
///
/// See also incbet.c.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// fdtrc domain    a<=0, b<=0, x<0       NaN
/// ```
///
/// (Unlike in xsf, where the value returned on a domain error is documented as 0.0.)
pub(crate) fn fdtrc(a: f64, b: f64, x: f64) -> f64 {
    if (a <= 0.0) || (b <= 0.0) || (x < 0.0) {
        return f64::NAN;
    }
    let w = b / (b + a * x);
    // Unlike in xsf (scipy/xsf#358), `1 - w` is passed as well (see `betainc_xy`), e.g.
    // fdtrc(0.00118, 1.46e11, 0.00286) was 1 instead of 0.007.
    betainc_xy(0.5 * b, 0.5 * a, w, a * x / (b + a * x))
}

/// F distribution
///
/// Returns the area from zero to x under the F density
/// function (also known as Snedcor's density or the
/// variance ratio density).  This is the density
/// of x = (u1/df1)/(u2/df2), where u1 and u2 are random
/// variables having Chi square distributions with df1
/// and df2 degrees of freedom, respectively.
///
/// The incomplete beta integral is used, according to the
/// formula
///
/// ```text
///     P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).
/// ```
///
/// The arguments a and b are greater than zero, and x is
/// nonnegative.
///
/// ACCURACY:
///
/// Tested at random points (a,b,x).
///
/// ```text
///                x     a,b                     Relative error:
/// arithmetic  domain  domain     # trials      peak         rms
///    IEEE      0,1    0,100        20000      8.9e-14     1.3e-14
///    IEEE      1,5    0,100        20000      9.6e-14     2.7e-15
///    IEEE      0,1    1,10000      20000      4.0e-11     6.7e-12
///    IEEE      1,5    1,10000      20000      1.1e-11     2.9e-13
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
///
/// See also incbet.c.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// fdtr domain     a<=0, b<=0, x<0       NaN
/// ```
///
/// (Unlike in xsf, where the value returned on a domain error is documented as 0.0.)
pub(crate) fn fdtr(a: f64, b: f64, x: f64) -> f64 {
    if (a <= 0.0) || (b <= 0.0) || (x < 0.0) {
        return f64::NAN;
    }
    let mut w = a * x;
    w = w / (b + w);
    // Unlike in xsf (scipy/xsf#358), `1 - w` is passed as well (see `betainc_xy`), e.g.
    // fdtr(1.38, 0.112, 5.77e62) was 1 instead of 0.99974.
    betainc_xy(0.5 * a, 0.5 * b, w, b / (b + a * x))
}

/// Inverse of F distribution
///
/// Finds the F density argument x such that the integral
/// from -infinity to x of the F density is equal to the
/// given probability p.
///
/// This is accomplished using the inverse beta integral
/// function and the relations
///
/// ```text
///      z = incbi( df2/2, df1/2, 1-p )
///      x = df2 (1-z) / (df1 z).
/// ```
///
/// Note: the following relations hold for the inverse of
/// the uncomplemented F distribution:
///
/// ```text
///      z = incbi( df1/2, df2/2, p )
///      x = df2 z / (df1 (1-z)).
/// ```
///
/// ACCURACY:
///
/// Tested at random points (a,b,p).
///
/// ```text
///              a,b                     Relative error:
/// arithmetic  domain     # trials      peak         rms
///  For p between .001 and 1:
///    IEEE     1,100        20000      2.9e-14     3.3e-15
///    IEEE     1,10000      20000      9.1e-12     2.3e-13
///  For p between 10^-6 and 10^-3:
///    IEEE     1,100        20000      2.6e-11     2.6e-13
///    IEEE     1,10000      20000      2.8e-12     8.4e-14
/// ```
///
/// See also fdtrc.c.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// fdtri domain   p <= 0 or p > 1       NaN
///                a<=0, b<=0
/// ```
pub(crate) fn fdtri(a: f64, b: f64, mut y: f64) -> f64 {
    if (a <= 0.0) || (b <= 0.0) || (y <= 0.0) || (y > 1.0) {
        return f64::NAN;
    }
    y = 1.0 - y;
    /* Compute probability for x = 0.5.  */
    let mut w = incbet(0.5 * b, 0.5 * a, 0.5);
    /* If that is greater than y, then the solution w < .5.
     * Otherwise, solve at 1-y to remove cancellation in (b - b*w).  */
    if w > y || y < 0.001 {
        w = incbi(0.5 * b, 0.5 * a, y);
        (b - b * w) / (a * w)
    } else {
        w = incbi(0.5 * a, 0.5 * b, 1.0 - y);
        b * w / (a * (1.0 - w))
    }
}
