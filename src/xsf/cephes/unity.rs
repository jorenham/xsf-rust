//! Translated into pure Rust from `xsf/cephes/unity.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
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
//!
//! Unlike in xsf (scipy/xsf#348), the Taylor series of `lgam1p` is summed until it converges (at
//! most `MAXITER` terms, as in `log1pmx`), instead of stopping after 41 terms.

use crate::utils::{MACHEP, SCIPY_EULER};
use crate::xsf::cephes::{lgam, zeta};

// cephes/const.h
const MAXITER: u32 = 500;

/// Compute lgam(x + 1) around x = 0 using its Taylor series.
#[allow(clippy::float_cmp)]
fn lgam1p_taylor(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let mut res = -SCIPY_EULER * x;
    let mut xfac = -x;
    // (unlike in xsf, scipy/xsf#348) MAXITER instead of 42, since 41 terms don't converge for |x|
    // near 1/2
    for n in 2..MAXITER {
        xfac *= -x;
        let coeff = zeta(f64::from(n), 1.0) * xfac / f64::from(n);
        res += coeff;
        if coeff.abs() < MACHEP * res.abs() {
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
