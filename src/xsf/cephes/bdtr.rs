//! Translated into pure Rust from `xsf/cephes/bdtr.h` (scipy/xsf@bd9efcc).
//!
//! Cephes Math Library Release 2.3:  March, 1995
//! Copyright 1984, 1987, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! The Cephes `expm1` from `cephes/unity.h` is replaced by `f64::exp_m1` from the Rust standard
//! library, like in `xsf/exp.rs`. In `bdtri`, `log1p(y - 1)` is replaced by `ln(y)`, which is
//! more accurate, since `y - 1` is exact for `y > 0.8`.

use crate::xsf::cephes::incbet::{incbet, incbet_xy};
use crate::xsf::cephes::incbi::incbi;

/// Complemented binomial distribution
///
/// Returns the sum of the terms k+1 through n of the Binomial
/// probability density:
///
/// ```text
///   n
///   --  ( n )   j      n-j
///   >   (   )  p  (1-p)
///   --  ( j )
///  j=k+1
/// ```
///
/// The terms are not summed directly; instead the incomplete
/// beta integral is employed, according to the formula
///
/// ```text
/// y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
/// ```
///
/// The arguments must be positive, with p ranging from 0 to 1.
///
/// ACCURACY:
///
/// Tested at random points (a,b,p).
///
/// ```text
///               a,b                     Relative error:
/// arithmetic  domain     # trials      peak         rms
///  For p between 0.001 and 1:
///    IEEE     0,100         1000      2.0e-13     6.8e-15
///  For p between 0 and .001:
///    IEEE     0,100       100000      4.0e-13     2.7e-14
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// bdtrc domain      x<0, x>1, n<k       NaN
/// ```
pub(crate) fn bdtrc(k: f64, n: i32, p: f64) -> f64 {
    let fk = k.floor();

    if p.is_nan() || k.is_nan() {
        return f64::NAN;
    }

    if p < 0.0 || p > 1.0 || f64::from(n) < fk {
        return f64::NAN;
    }

    if fk < 0.0 {
        return 1.0;
    }

    if fk == f64::from(n) {
        return 0.0;
    }

    let dn = f64::from(n) - fk;
    if k == 0.0 {
        if p < 0.01 {
            -(dn * (-p).ln_1p()).exp_m1()
        } else {
            1.0 - (1.0 - p).powf(dn)
        }
    } else {
        incbet(fk + 1.0, dn, p)
    }
}

/// Binomial distribution
///
/// Returns the sum of the terms 0 through k of the Binomial
/// probability density:
///
/// ```text
///   k
///   --  ( n )   j      n-j
///   >   (   )  p  (1-p)
///   --  ( j )
///  j=0
/// ```
///
/// The terms are not summed directly; instead the incomplete
/// beta integral is employed, according to the formula
///
/// ```text
/// y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
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
///  For p between 0.001 and 1:
///    IEEE     0,100         2000      4.9e-14     2.2e-15
/// ```
///
/// See also incbet.c.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// bdtr domain         k < 0            NaN
///                     n < k
///                     x < 0, x > 1
/// ```
pub(crate) fn bdtr(k: f64, n: i32, p: f64) -> f64 {
    let fk = k.floor();

    if p.is_nan() || k.is_nan() {
        return f64::NAN;
    }

    if p < 0.0 || p > 1.0 || fk < 0.0 || f64::from(n) < fk {
        return f64::NAN;
    }

    if fk == f64::from(n) {
        return 1.0;
    }

    let dn = f64::from(n) - fk;
    // Unlike in xsf (scipy/xsf#358), the rounding error of `1 - p` (exact if `p >= 0.5`) is
    // accounted for.
    //
    // For `k = 0`, `(1 - p)^n = q^n (1 + e/q)^n ~ q^n exp(n e/q)`, with the rounded `q = 1 - p`
    // and its exact rounding error `e = (1 - q) - p`.
    if fk == 0.0 {
        let q = 1.0 - p;
        let qn = q.powf(dn);
        if p < 0.5 {
            qn * (dn * ((1.0 - q) - p) / q).exp()
        } else {
            qn
        }
    } else {
        incbet_xy(dn, fk + 1.0, 1.0 - p, p)
    }
}

/// Inverse binomial distribution
///
/// Finds the event probability p such that the sum of the
/// terms 0 through k of the Binomial probability density
/// is equal to the given cumulative probability y.
///
/// This is accomplished using the inverse beta integral
/// function and the relation
///
/// ```text
/// 1 - p = incbi( n-k, k+1, y ).
/// ```
///
/// ACCURACY:
///
/// Tested at random points (a,b,p).
///
/// ```text
///               a,b                     Relative error:
/// arithmetic  domain     # trials      peak         rms
///  For p between 0.001 and 1:
///    IEEE     0,100          500      1.2e-15     1.6e-16
///    IEEE     0,10000        500      8.4e-16     9.1e-17
///  For p between 10^-6 and 0.001:
///    IEEE     0,100          500      2.9e-16     7.7e-17
///    IEEE     0,10000        500      5.0e-16     7.4e-17
/// ```
///
/// Here y = bdtr(k, n, p) for random p in the given range, and the
/// error is relative to the exact inverse of that (rounded) y, which
/// can lie outside that range, e.g. if y rounds to almost 1. Cases
/// where y underflows or rounds to 1 were excluded (e.g. 171 of the
/// 500 trials for a,b up to 10000 and p between 0.001 and 1).
///
/// See also incbi.c.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// bdtri domain     k < 0, n <= k         NaN
///                  x < 0, x > 1
/// ```
pub(crate) fn bdtri(k: f64, n: i32, y: f64) -> f64 {
    let fk = k.floor();

    if k.is_nan() {
        return f64::NAN;
    }

    if y < 0.0 || y > 1.0 || fk < 0.0 || f64::from(n) <= fk {
        return f64::NAN;
    }

    let dn = f64::from(n) - fk;

    if fk == f64::from(n) {
        return 1.0;
    }

    if fk == 0.0 {
        if y > 0.8 {
            -(y.ln() / dn).exp_m1()
        } else {
            1.0 - y.powf(1.0 / dn)
        }
    } else {
        let dk = fk + 1.0;
        // 1 - y is exact for y >= 0.5
        if y >= 0.5 {
            incbi(dk, dn, 1.0 - y).0
        } else {
            1.0 - incbi(dn, dk, y).0
        }
    }
}
