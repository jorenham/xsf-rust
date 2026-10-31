//! Translated into pure Rust from `xsf/ellint_carlson/_rf.hh` (xsf main).
//!
//! ```text
//! References
//! [1] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
//!     Functions," NIST, US Dept. of Commerce.
//!     https://dlmf.nist.gov/19.16.E1
//! [2] B. C. Carlson, "Numerical computation of real or complex elliptic
//!     integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
//!     https://arxiv.org/abs/math/9409227
//!     https://doi.org/10.1007/BF02198293
//! ```

use super::arithmetic::{comp_horner, dot2, ocrt, sum2};
use super::common::{ExitStatus, abscmp, config, constants, sort_by};
use super::typing::Num;
use super::{argcheck, common};

#[inline]
pub(super) fn agm_update<T: Num>(x: &mut T, y: &mut T) {
    let xnext = (*x + *y) * 0.5;
    let ynext = x.prod(*y).sqrt();

    *x = xnext;
    *y = ynext;
}

fn rf0<T: Num>(x: T, y: T, rerr: f64) -> (ExitStatus, T) {
    let mut status = ExitStatus::Success;
    let rsq = 2.0 * rerr.sqrt();

    let mut xm = x.sqrt();
    let mut ym = y.sqrt();
    let mut m = 0;
    while (xm - ym).abs() >= rsq * xm.abs().min(ym.abs()) {
        if m > config::MAX_ITER {
            status = ExitStatus::NIter;
            break;
        }

        agm_update(&mut xm, &mut ym);

        m += 1;
    }

    let res = T::from(constants::PI).quot(xm + ym);
    (status, res)
}

pub(super) fn rf<T: Num>(x: T, y: T, z: T, rerr: f64) -> (ExitStatus, T) {
    let mut status = ExitStatus::Success;
    if argcheck::invalid_rerr(rerr, 3.0e-4) {
        return (ExitStatus::BadRerr, T::NAN);
    }

    if x.ph_good() && y.ph_good() && z.ph_good() {
        if x.isinf() || y.isinf() || z.isinf() {
            return (ExitStatus::Success, T::from(0.0));
        }
    } else {
        return (ExitStatus::BadArgs, T::NAN);
    }

    let mut cct1 = [x, y, z];
    sort_by(&mut cct1, abscmp);
    let mut xm = cct1[0];
    let mut ym = cct1[1];
    let mut zm = cct1[2];
    if xm.too_small() {
        if ym.too_small() {
            return (ExitStatus::Singular, T::HUGE);
        }
        let (status, tmpres) = rf0(ym, zm, rerr * 0.5);
        /* Correction for non-zero x, see Eq. 19.27.3 in
         * https://dlmf.nist.gov/19.27.E3 */
        let res = tmpres - xm.quot(ym.prod(zm)).sqrt();
        return (status, res);
    }

    let mut am = sum2(&cct1) / 3.0;
    let mut xxm = am - xm;
    let mut yym = am - ym;
    let mut fterm = common::max(&[xxm.abs(), yym.abs(), (am - zm).abs()]) / ocrt(3.0 * rerr);
    let mut m = 0;
    loop {
        let a_am = am.abs();
        if !(a_am <= fterm || a_am <= common::max(&[xxm.abs(), yym.abs(), (am - zm).abs()])) {
            break;
        }
        if m > config::MAX_ITER {
            status = ExitStatus::NIter;
            break;
        }
        cct1 = [xm.sqrt(), ym.sqrt(), zm.sqrt()];
        let cct2 = [cct1[1], cct1[2], cct1[0]];
        let lam = dot2(&cct1, &cct2);
        am = (am + lam) * 0.25;
        xm = (xm + lam) * 0.25;
        ym = (ym + lam) * 0.25;
        zm = (zm + lam) * 0.25;
        xxm = xxm * 0.25;
        yym = yym * 0.25;
        fterm *= 0.25;

        m += 1;
    }
    /* Burn some extra cycles re-balancing Am as the "true" centroid. */
    cct1 = [xm, ym, zm];
    am = sum2(&cct1) / 3.0;
    xxm = xxm.quot(am);
    yym = yym.quot(am);
    /* Prepare the E_2 and E_3 terms used in the expansion */
    let zzm = -(xxm + yym);
    let e2 = xxm.prod(yym) - zzm.prod(zzm);
    let e3 = xxm.prod(yym.prod(zzm));
    /* Evaluate the 7th-degree expansion using the E_2 and E3 terms, following
     * Eq. 19.36.1 of [1], https://dlmf.nist.gov/19.36.E1
     * The order of expansion is higher than that in Eq. (14) of Ref. [2]. */
    let mut s = comp_horner(e2, &constants::RF_C1);
    s = s + e3.prod(comp_horner(e2, &constants::RF_C2) + e3 * constants::RF_C33);
    s = s / constants::RF_DENOM;
    s = s + 1.0;

    let res = s.quot(am.sqrt());
    (status, res)
}
