//! Translated into pure Rust from `xsf/cephes/unity.h` (scipy/xsf@e58388a).
//!
//! Relative error approximations for function arguments near
//! unity.
//!
//!    lgam1p(x) = lgam(1+x)
//!
//! Scipy changes:
//! - 06-10-2016: added lgam1p
//!
//! The other functions in `unity.h` are translated in `xsf/log.rs` (`log1pmx`) and `xsf/trig.rs`
//! (`cosm1`), or replaced by `f64::ln_1p` and `f64::exp_m1` from the Rust standard library.

use crate::utils::EGAMMA;
use crate::xsf::cephes::consts::MACHEP;
use crate::xsf::cephes::{lgam, zeta};

// cephes/const.h
const MAXITER: u32 = 500;

/// Compute lgam(x + 1) around x = 0 using its Taylor series.
fn lgam1p_taylor(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let mut res = -EGAMMA * x;
    let mut xfac = -x;
    for n in 2..MAXITER {
        xfac *= -x;
        let coeff = zeta(f64::from(n), 1.0) * xfac / f64::from(n);
        res += coeff;
        if coeff.abs() <= MACHEP * res.abs() {
            break;
        }
    }

    res
}

/// Compute lgam(x + 1).
///
/// see `scipy.special._ufuncs._lgam1p`
#[doc(hidden)]
#[must_use]
#[inline]
pub fn lgam1p(x: f64) -> f64 {
    if x.abs() <= 0.5 {
        lgam1p_taylor(x)
    } else if (x - 1.0).abs() < 0.5 {
        x.ln() + lgam1p_taylor(x - 1.0)
    } else {
        lgam(x + 1.0)
    }
}
