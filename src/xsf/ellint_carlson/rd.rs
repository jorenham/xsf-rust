//! Translated into pure Rust from `xsf/ellint_carlson/_rd.hh` (xsf main).
//!
//! ```text
//! References
//! [1] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
//!     Functions," NIST, US Dept. of Commerce.
//!     https://dlmf.nist.gov/19.16.E5
//! [2] B. C. Carlson, "Numerical computation of real or complex elliptic
//!     integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
//!     https://arxiv.org/abs/math/9409227
//!     https://doi.org/10.1007/BF02198293
//! ```

use super::arithmetic::{comp_horner, dot2, ndot2, nsum2, ocrt, sum2_acc};
use super::common::{ExitStatus, config, constants};
use super::typing::Num;
use super::{argcheck, common};

#[allow(clippy::too_many_lines)]
pub(super) fn rd<T: Num>(x: T, y: T, z: T, rerr: f64) -> (ExitStatus, T) {
    let mut cct1 = [T::from(0.0); 6];
    let mut cct2 = [T::from(0.0); 6];

    let mut status = ExitStatus::Success;
    if argcheck::invalid_rerr(rerr, 1.0e-4) {
        return (ExitStatus::BadRerr, T::NAN);
    }

    // Unlike in xsf (scipy/xsf#376), huge or tiny arguments are rescaled
    if let Some(f) = common::scale_factor(&[x, y, z], f64::MIN_POSITIVE) {
        let (status, res) = rd(x * f * f, y * f * f, z * f * f, rerr);
        return (status, res * f * f * f);
    }

    if !(x.ph_good() && y.ph_good() && z.ph_good()) {
        return (ExitStatus::BadArgs, T::NAN);
    }

    if z.too_small() {
        return (ExitStatus::Singular, T::HUGE);
    }

    if x.isinf() || y.isinf() || z.isinf() {
        return (status, T::from(0.0));
    }

    if x.too_small() && y.too_small() {
        return (ExitStatus::Singular, T::HUGE);
    }

    cct1[0] = x;
    cct1[1] = y;
    cct1[2] = z;
    cct1[3] = z;
    cct1[4] = z;
    let mut am = nsum2(&cct1, 5) / 5.0;
    let mut xm = x;
    let mut ym = y;
    let mut zm = z;
    let mut xxm = am - xm;
    let mut yym = am - ym;
    let mut fterm = common::max(&[xxm.abs(), yym.abs(), (am - z).abs()]) / ocrt(rerr / 5.0);
    let mut d4m = 1.0;
    let mut adt = T::from(0.0);
    let mut ade = T::from(0.0);

    let mut m = 0;
    let mut tmp;
    loop {
        let a_am = am.abs();
        if !(a_am <= fterm || a_am <= common::max(&[xxm.abs(), yym.abs(), (am - zm).abs()])) {
            break;
        }
        if m > config::MAX_ITER {
            status = ExitStatus::NIter;
            break;
        }
        cct1[0] = xm.sqrt();
        cct2[2] = cct1[0];
        cct1[1] = ym.sqrt();
        cct2[0] = cct1[1];
        cct1[2] = zm.sqrt();
        cct2[1] = cct1[2];
        let lam = ndot2(&cct1, &cct2, 3);

        /* The cumulative sum term in Ref. [2], Eq. (41), implemented by
         * accumulating the summation term into adt with compensation
         * term ade */
        tmp = T::from(d4m).quot(cct1[2].prod(zm + lam));
        sum2_acc(tmp, &mut adt, &mut ade);

        am = (am + lam) * 0.25;
        xm = (xm + lam) * 0.25;
        ym = (ym + lam) * 0.25;
        zm = (zm + lam) * 0.25;
        xxm = xxm * 0.25;
        yym = yym * 0.25;
        fterm *= 0.25;
        d4m *= 0.25;

        m += 1;
    }
    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    cct1[3] = zm;
    cct1[4] = zm;
    am = nsum2(&cct1, 5) / 5.0;
    xxm = xxm.quot(am);
    yym = yym.quot(am);
    let zzm = (xxm + yym) / -3.0;
    /* Prepare the "elementary" terms E_n as in Eqs. (39-40) in Ref. [2] */
    let xy = xxm.prod(yym);
    let zz2 = zzm.prod(zzm);
    let e2 = xy - zz2 * 6.0;
    let e3 = (xy * 3.0 - zz2 * 8.0).prod(zzm);
    let e4 = (xy - zz2).prod(zz2) * 3.0;
    let e5 = xy.prod(zz2).prod(zzm);
    let mut t = am.sqrt();
    tmp = T::from(d4m).quot(t.prod(t).prod(t));
    /* Evaluate the 7th-degree expansion using the E_n terms, following
     * Eq. 19.36.2 of [1], https://dlmf.nist.gov/19.36#E2
     * The order of expansion is higher than that in Eq. (41) of Ref. [2]. */
    cct1[0] = comp_horner(e2, &constants::RDJ_C1);
    cct1[1] = comp_horner(e3, &constants::RDJ_C2);
    cct1[2] = comp_horner(e2, &constants::RDJ_C3);
    cct1[3] = comp_horner(e2, &constants::RDJ_C4);
    cct1[4] = comp_horner(e2, &constants::RDJ_C5);
    cct1[5] = e3 * constants::RDJ_C5[1];

    cct2[0] = T::from(1.0);
    cct2[1] = T::from(1.0);
    cct2[2] = e3;
    cct2[3] = e4;
    cct2[4] = e5;
    cct2[5] = e4;
    t = dot2(&cct1, &cct2) / constants::RDJ_DENOM + T::from(1.0);

    /* Combine in to final result using compensated dot. See Eq. (41) of Ref.
     * [2]. */
    cct1[0] = tmp;
    cct1[1] = T::from(3.0);
    cct1[2] = T::from(3.0);
    cct2[0] = t;
    cct2[1] = adt;
    cct2[2] = ade;
    let res = ndot2(&cct1, &cct2, 3);
    (status, res)
}
