//! Translated into pure Rust from `xsf/ellint_carlson/_rc.hh` (xsf main).
//!
//! ```text
//! References:
//! [1] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
//!     Functions," NIST, US Dept. of Commerce.
//!     https://dlmf.nist.gov/19.16.E6
//! [2] B. C. Carlson, "Numerical computation of real or complex elliptic
//!     integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
//!     https://arxiv.org/abs/math/9409227
//!     https://doi.org/10.1007/BF02198293
//! ```

use super::arithmetic::{comp_horner, ocrt};
use super::common::{ExitStatus, config, constants, is_horrible};
use super::typing::Num;
use super::{argcheck, common};

pub(super) fn rc<T: Num>(x: T, y: T, rerr: f64) -> (ExitStatus, T) {
    let mut status = ExitStatus::Success;
    if argcheck::invalid_rerr(rerr, 2.0e-4) {
        return (ExitStatus::BadRerr, T::NAN);
    }

    /* Cauchy principal value with real negative y */
    if argcheck::too_small(y.imag()) && (y.real() < 0.0) {
        /* Ref[2], Eq. 2.14 or Eq. (21) in the arXiv preprint */
        let (status, tmpres) = rc(x - y, -y, rerr);
        let res = if is_horrible(status) {
            T::NAN
        } else {
            tmpres.prod(x.quot(x - y).sqrt())
        };
        return (status, res);
    } else if y.too_small() || !x.ph_good() {
        return (ExitStatus::BadArgs, T::NAN);
    } else if x.isinf() || y.isinf() {
        return (ExitStatus::Success, T::from(0.0));
    }

    let mut am = (x + y * 2.0) / 3.0;
    let mut fterm = (am - x).abs() / ocrt(3.0 * rerr);

    let mut xm = x;
    let mut ym = y;
    let mut sm = y - am;

    let mut m = 0; /* loop variable */
    while common::max(&[(xm - ym).abs(), fterm]) >= am.abs() {
        if m > config::MAX_ITER {
            status = ExitStatus::NIter;
            break;
        }

        let lam = (xm.sqrt() * 2.0).prod(ym.sqrt()) + ym;
        am = (am + lam) * 0.25;
        xm = (xm + lam) * 0.25;
        ym = (ym + lam) * 0.25;
        sm = sm * 0.25;
        fterm *= 0.25;

        m += 1;
    }
    am = (xm + ym + ym) / 3.0;
    sm = sm.quot(am);
    /* Eq. (20) of Ref[2] */
    let res = comp_horner(sm, &constants::RC_C).quot(am.sqrt() * constants::RC_C[0]);
    (status, res)
}
