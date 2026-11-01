//! Translated into pure Rust from `xsf/ellint_carlson/_rg.hh` (xsf main).
//!
//! ```text
//! References
//! [1] B. C. Carlson, "Numerical computation of real or complex elliptic
//!     integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
//!     https://arxiv.org/abs/math/9409227
//!     https://doi.org/10.1007/BF02198293
//! [2] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
//!     Functions," NIST, US Dept. of Commerce.
//!     https://dlmf.nist.gov/19.16.E1
//!     https://dlmf.nist.gov/19.20.ii
//! ```

use num_complex::Complex64;

use super::arithmetic::{dot2, eft_prod, eft_sum};
use super::common::{self, ExitStatus, abscmp, config, constants, is_horrible, sort_by};
use super::rd::rd;
use super::rf::agm_update;
use super::typing::Num;

/// The complex `arithmetic::aux::rg_dot2_acc`; the real one is `arithmetic::fdot2_acc`
#[inline]
pub(super) fn rg_dot2_acc_cplx(
    fac: f64,
    term: Complex64,
    acc: &mut Complex64,
    cor: &mut Complex64,
) {
    let (ar, cr) = eft_prod(fac, term.re);
    let (ai, ci) = eft_prod(fac, term.im);
    let (p, q) = eft_sum(*acc, Complex64::new(ar, ai));
    *acc = p;
    *cor += q + Complex64::new(cr, ci);
}

fn rg0<T: Num>(x: T, y: T, rerr: f64) -> (ExitStatus, T) {
    let mut status = ExitStatus::Success;
    let rsq = 2.0 * rerr.sqrt();
    let mut fac = 0.25;

    let mut xm = x.sqrt();
    let mut ym = y.sqrt();
    let mut dm = (xm + ym) * 0.5;
    let mut sum = (-dm).prod(dm);
    let mut cor = T::from(0.0);
    dm = xm - ym;
    let mut m = 0;
    while dm.abs() >= rsq * xm.abs().min(ym.abs()) {
        if m > config::MAX_ITER {
            status = ExitStatus::NIter;
            break;
        }

        agm_update(&mut xm, &mut ym);
        /* Ref[1], Eq. 2.39 or Eq. (46) in the arXiv preprint */
        fac *= 2.0;
        dm = xm - ym;
        T::rg_dot2_acc(fac, dm.prod(dm), &mut sum, &mut cor);

        m += 1;
    }

    let mut res = T::from(constants::PI).quot(xm + ym);
    res = res.prod((sum + cor) * -0.5);
    (status, res)
}

pub(super) fn rg<T: Num>(x: T, y: T, z: T, rerr: f64) -> (ExitStatus, T) {
    let mut status = ExitStatus::Success;
    if super::argcheck::invalid_rerr(rerr, 1.0e-4) {
        return (ExitStatus::BadRerr, T::NAN);
    }

    // Unlike in xsf (scipy/xsf#376), huge or tiny arguments are rescaled
    if let Some(f) = common::scale_factor(&[x, y, z], f64::MIN_POSITIVE) {
        let (status, res) = rg(x * f * f, y * f * f, z * f * f, rerr);
        return (status, res / f);
    }

    let mut cct = [x, y, z];
    sort_by(&mut cct, abscmp);

    if (cct[0].isinf() || cct[1].isinf() || cct[2].isinf())
        && cct[0].ph_good()
        && cct[1].ph_good()
        && cct[2].ph_good()
    {
        return (ExitStatus::Singular, T::HUGE);
    }

    if cct[0].too_small() {
        if cct[1].too_small() {
            /* Special case -- also covers the case of z ~ zero. */
            return (status, cct[2].sqrt() * 0.5);
        }
        /* Special case -- use the AGM algorithm. */
        return rg0(cct[1], cct[2], rerr);
    }

    /* Ref[2], Eq. 19.21.11 (second identity) <https://dlmf.nist.gov/19.21.E11>
     * i.e. 6R_G() = sum of [ x * (y + z) * R_D() ] over cyclic permutations of
     * (x, y, z).
     * This cyclic form is manifestly symmetric and is preferred over
     * Eq. 19.21.10 ibid.
     * Here we put the three R_D terms in the buffer cct2, and the
     * x * (y + z) = dot({x, x}, {y, z}) terms in cct1. */
    let mut cct1 = [T::from(0.0); 3];
    let mut cct2 = [T::from(0.0); 3];

    macro_rules! check_status_or_fail {
        ($status_tmp:expr) => {
            if $status_tmp != ExitStatus::Success {
                status = $status_tmp;
            }
            if is_horrible(status) {
                return (status, T::NAN);
            }
        };
    }

    let mut status_tmp;

    (status_tmp, cct2[0]) = rd(y, z, x, rerr);
    check_status_or_fail!(status_tmp);
    (status_tmp, cct2[1]) = rd(z, x, y, rerr);
    check_status_or_fail!(status_tmp);
    (status_tmp, cct2[2]) = rd(x, y, z, rerr);
    check_status_or_fail!(status_tmp);

    /* Fill the cct1 buffer via the intermediate buffers tm1, tm2
     * that are dotted together (using the compensation algorithm). */
    let mut tm1 = [x, x];
    let mut tm2 = [y, z];
    cct1[0] = dot2(&tm1, &tm2);
    tm1 = [y, y];
    tm2[0] = x;
    cct1[1] = dot2(&tm1, &tm2);
    tm1 = [z, z];
    tm2[1] = y;
    cct1[2] = dot2(&tm1, &tm2);

    let res = dot2(&cct1, &cct2) / 6.0;
    (status, res)
}
