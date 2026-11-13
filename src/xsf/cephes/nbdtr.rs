//! Translated into pure Rust from `xsf/cephes/nbdtr.h` (scipy/xsf@09c93c7).
//!
//! Cephes Math Library Release 2.3:  March, 1995
//! Copyright 1984, 1987, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `k + 1` is computed in `f64`, so that it doesn't overflow for `k = i32::MAX`.

use crate::xsf::cephes::incbet::{incbet, incbet_xy};
use crate::xsf::cephes::incbi::incbi;

/// Complemented negative binomial distribution
///
/// Returns the sum of the terms k+1 to infinity of the negative
/// binomial distribution:
///
/// ```text
///   inf
///   --  ( n+j-1 )   n      j
///   >   (       )  p  (1-p)
///   --  (   j   )
///  j=k+1
/// ```
///
/// The terms are not computed individually; instead the incomplete
/// beta integral is employed, according to the formula
///
/// ```text
/// y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p ).
/// ```
///
/// The arguments must be positive, with p ranging from 0 to 1.
///
/// ACCURACY:
///
/// Tested at random points (a,b,p), with p between 0 and 1.
///
/// ```text
///               a,b                     Relative error:
/// arithmetic  domain     # trials      peak         rms
///    IEEE     0,100         2000      3.1e-13     9.3e-15
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
///
/// See also incbet.c.
pub(crate) fn nbdtrc(k: i32, n: i32, p: f64) -> f64 {
    if (p < 0.0) || (p > 1.0) || k < 0 {
        // domerr
        return f64::NAN;
    }

    let dk = f64::from(k) + 1.0;
    let dn = f64::from(n);
    // Unlike in xsf (scipy/xsf#358), `p` is passed as the exact `1 - (1 - p)`.
    incbet_xy(dk, dn, 1.0 - p, p)
}

/// Negative binomial distribution
///
/// Returns the sum of the terms 0 through k of the negative
/// binomial distribution:
///
/// ```text
///   k
///   --  ( n+j-1 )   n      j
///   >   (       )  p  (1-p)
///   --  (   j   )
///  j=0
/// ```
///
/// In a sequence of Bernoulli trials, this is the probability
/// that k or fewer failures precede the nth success.
///
/// The terms are not computed individually; instead the incomplete
/// beta integral is employed, according to the formula
///
/// ```text
/// y = nbdtr( k, n, p ) = incbet( n, k+1, p ).
/// ```
///
/// The arguments must be positive, with p ranging from 0 to 1.
///
/// ACCURACY:
///
/// Tested at random points (a,b,p), with p between 0 and 1.
///
/// ```text
///               a,b                     Relative error:
/// arithmetic  domain     # trials      peak         rms
///    IEEE     0,100         2000      1.6e-13     4.1e-15
/// ```
///
/// See also incbet.c.
pub(crate) fn nbdtr(k: i32, n: i32, p: f64) -> f64 {
    if (p < 0.0) || (p > 1.0) || k < 0 {
        // domerr
        return f64::NAN;
    }
    let dk = f64::from(k) + 1.0;
    let dn = f64::from(n);
    incbet(dn, dk, p)
}

/// Functional inverse of negative binomial distribution
///
/// Finds the argument p such that nbdtr(k,n,p) is equal to y.
///
/// ACCURACY:
///
/// Tested at random points (a,b,y), with y between 0 and 1.
///
/// ```text
///               a,b                     Relative error:
/// arithmetic  domain     # trials      peak         rms
///    IEEE     0,100         2000      3.1e-15     3.6e-16
/// ```
///
/// See also incbi.c.
pub(crate) fn nbdtri(k: i32, n: i32, p: f64) -> f64 {
    if (p < 0.0) || (p > 1.0) || k < 0 {
        // domerr
        return f64::NAN;
    }
    let dk = f64::from(k) + 1.0;
    let dn = f64::from(n);
    incbi(dn, dk, p).0
}
