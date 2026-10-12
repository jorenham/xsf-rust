//! Translated into pure Rust from `xsf/cephes/pdtr.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.3:  March, 1995
//! Copyright 1984, 1987, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `k + 1` in `pdtri` is computed in `f64`, so that it doesn't overflow for
//! `k = i32::MAX`, and the value returned on a domain error is documented as NaN instead of 0.0.

use crate::xsf::cephes::igam::{igam, igamc};
use crate::xsf::cephes::igami::igamci;

/// Complemented poisson distribution
///
/// Returns the sum of the terms k+1 to infinity of the Poisson
/// distribution:
///
/// ```text
///  inf.       j
///   --   -m  m
///   >   e    --
///   --       j!
///  j=k+1
/// ```
///
/// The terms are not summed directly; instead the incomplete
/// Gamma integral is employed, according to the formula
///
/// ```text
/// y = pdtrc( k, m ) = igam( k+1, m ).
/// ```
///
/// The arguments must both be nonnegative.
///
/// ACCURACY:
///
/// See igam.c.
#[allow(clippy::float_cmp)]
pub(crate) fn pdtrc(k: f64, m: f64) -> f64 {
    if k < 0.0 || m < 0.0 {
        return f64::NAN;
    }
    if m == 0.0 {
        return 0.0;
    }
    let v = k.floor() + 1.0;
    igam(v, m)
}

/// Poisson distribution
///
/// Returns the sum of the first k terms of the Poisson
/// distribution:
///
/// ```text
///   k         j
///   --   -m  m
///   >   e    --
///   --       j!
///  j=0
/// ```
///
/// The terms are not summed directly; instead the incomplete
/// Gamma integral is employed, according to the relation
///
/// ```text
/// y = pdtr( k, m ) = igamc( k+1, m ).
/// ```
///
/// The arguments must both be nonnegative.
///
/// ACCURACY:
///
/// See igamc().
#[allow(clippy::float_cmp)]
pub(crate) fn pdtr(k: f64, m: f64) -> f64 {
    if k < 0.0 || m < 0.0 {
        return f64::NAN;
    }
    if m == 0.0 {
        return 1.0;
    }
    let v = k.floor() + 1.0;
    igamc(v, m)
}

/// Inverse Poisson distribution
///
/// Finds the Poisson variable x such that the integral
/// from 0 to x of the Poisson density is equal to the
/// given probability y.
///
/// This is accomplished using the inverse Gamma integral
/// function and the relation
///
/// ```text
///    m = igamci( k+1, y ).
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
/// pdtri domain    y < 0 or y >= 1       NaN
///                     k < 0
/// ```
pub(crate) fn pdtri(k: i32, y: f64) -> f64 {
    if (k < 0) || (y < 0.0) || (y >= 1.0) {
        return f64::NAN;
    }
    let v = f64::from(k) + 1.0;
    igamci(v, y)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_pdtri_k_max() {
        // `k + 1` overflows `int` in xsf
        let m = super::pdtri(i32::MAX, 0.5);
        assert!(m.is_finite() && m > f64::from(i32::MAX), "{m}");
    }
}
