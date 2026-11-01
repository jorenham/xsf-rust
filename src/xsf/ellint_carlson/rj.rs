//! Translated into pure Rust from `xsf/ellint_carlson/_rj.hh` (xsf main).
//!
//! ```text
//! References
//! [1] B. C. Carlson, "Numerical computation of real or complex elliptic
//!     integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
//!     https://arxiv.org/abs/math/9409227
//!     https://doi.org/10.1007/BF02198293
//! [2] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
//!     Functions," NIST, US Dept. of Commerce.
//!     https://dlmf.nist.gov/19.20.iii
//! [3] B. C. Carlson, J. FitzSimmons, "Reduction Theorems for Elliptic
//!     Integrands with the Square Root of Two Quadratic Factors," J.
//!     Comput. Appl. Math., vol. 118, nos. 1-2, pp. 71-85, 2000.
//!     https://doi.org/10.1016/S0377-0427(00)00282-X
//! [4] F. Johansson, "Numerical Evaluation of Elliptic Functions, Elliptic
//!     Integrals and Modular Forms," in J. Blumlein, C. Schneider, P.
//!     Paule, eds., "Elliptic Integrals, Elliptic Functions and Modular
//!     Forms in Quantum Field Theory," pp. 269-293, 2019 (Cham,
//!     Switzerland: Springer Nature Switzerland)
//!     https://arxiv.org/abs/1806.06725
//!     https://doi.org/10.1007/978-3-030-04480-0
//! [5] B. C. Carlson, J. L. Gustafson, "Asymptotic Approximations for
//!     Symmetric Elliptic Integrals," SIAM J. Math. Anls., vol. 25, no. 2,
//!     pp. 288-303, 1994.
//!     https://arxiv.org/abs/math/9310223
//!     https://doi.org/10.1137/S0036141092228477
//! [6] J. L. Gustafson, "Asymptotic Formulas for Elliptic Integrals", Doctoral
//!     dissertataion, Dept. of Mathematics, Iowa State University, 1982
//!     http://johngustafson.net/pubs/pub03/JGDissertation.pdf
//! ```
//!
//! The complex `std::atan` is `catan` from `xsf/complex.rs`, which is a translation of FreeBSD's
//! `catan`, so it isn't bit-identical to the glibc one that is used by g++.

use num_complex::Complex64;

use super::arithmetic::{comp_horner, dot2, ndot2, nsum2, ocrt};
use super::common::{ExitStatus, config, constants, is_horrible, is_troublesome, sort_by};
use super::rc::rc;
use super::rf::rf;
use super::rg::rg;
use super::typing::Num;
use super::{argcheck, common};
use crate::xsf::complex::{atanh, catan, csqrt};

mod rjimpl {
    use super::{ExitStatus, Num, argcheck, config, is_horrible, is_troublesome, ndot2, nsum2};
    use super::{rc, rf, rj};

    #[inline]
    pub(super) fn asymp_zero(r: f64) -> bool {
        (0.0 < r) && (r <= config::ASYM_ZERO_UL)
    }

    #[inline]
    pub(super) fn abs_close_zero(r: f64) -> bool {
        (0.0 < r) && (r <= config::ASYM_CLOSE_UL)
    }

    #[derive(Clone, Copy, PartialEq, Eq)]
    pub(super) enum AsymFlag {
        Nothing,
        Hugep, /* x, y, z << p */
        Tinyp, /* p << geom. mean of x, y, z */
        Hugey, /* max(x, p) << min(y, z) == y */
        #[allow(dead_code)] // see the commented-out code in `rj_asym_conf`
        Tinyx, /* x << min(y, z, p) == min(y, p) */
        Tinyy, /* max(x, y) == y << min(z, p) */
        Hugez, /* max(x, y, p) == max(y, p) << z */
    }

    #[allow(clippy::struct_excessive_bools)]
    #[derive(Default)]
    pub(super) struct ArgCases {
        pub(super) maybe_asymp: bool,   /* might be good for an asympt. case */
        pub(super) retry_caupv: bool,   /* should use Cauchy principal value */
        pub(super) good_infinity: bool, /* the "good" kind of directed infinity */
        pub(super) hit_pole: bool,      /* singular */
    }

    /* Comparison function by the real part. */
    #[inline]
    pub(super) fn rcmp<T: Num>(a: &T, b: &T) -> bool {
        a.real() < b.real()
    }

    /* Check whether the input arguments for RJ are out of domain, while set
     * corresponding flags in the class(ification) variable.
     *
     * NOTE: x, y, z must be in-order by real parts.
     */
    pub(super) fn good_args<T: Num>(x: T, y: T, z: T, p: T, classify: &mut ArgCases) -> bool {
        classify.hit_pole = x.too_small() && y.too_small() && z.ph_good() && !p.too_small();
        if classify.hit_pole {
            return false;
        }
        classify.good_infinity = (x.isinf() || y.isinf() || z.isinf() || p.isinf())
            && (x.ph_good() && y.ph_good() && z.ph_good());
        if classify.good_infinity {
            return false;
        }
        let xr = x.real();
        let xi = x.imag();

        let yr = y.real();
        let yi = y.imag();

        let zi = z.imag();

        let pr = p.real();
        let pi = p.imag();

        /* In the text above Eq. 2.26 of Ref[1]:
         * "If x, y, z are real and nonnegative, at most one of them is 0, and the
         * fourth variable of RJ is negative, the Cauchy principal value ..." */
        let xyzreal_nonneg_atmost1z = argcheck::too_small(xi)
            && argcheck::too_small(yi)
            && argcheck::too_small(zi)
            && (xr >= 0.0)
            && (yr > 0.0);
        if argcheck::too_small(pi) && xyzreal_nonneg_atmost1z {
            classify.retry_caupv = pr < 0.0;
            if classify.retry_caupv {
                return false;
            }
            /* Ref[2], Sec. 19.27 <https://dlmf.nist.gov/19.27.iv>:
             * "Assume x, y, and z are real and nonnegative, at most one of them is
             * 0, and p > 0" */
            classify.maybe_asymp = pr > 0.0;
            if classify.maybe_asymp {
                return true;
            }
        }
        /* Ref[1], in text above Eq. 2.17 ("Algorithm for RJ")
         * "Let x, y, z have nonnegative real part and at most one of them* be
         * 0, while Re p > 0."
         *     * By "them", Carlson seems to have meant the numbers x, y, z
         *       themselves, rather than their "real parts". */
        let x0 = x.too_small();
        let y0 = y.too_small();
        let z0 = z.too_small();
        if (pr > 0.0) && (xr >= 0.0) && !((y0 && (x0 || z0)) || (x0 && z0)) {
            return true;
        }
        /* "Alternatively, if p != 0 and |ph p| < pi ..." */
        if !p.too_small() && p.ph_good() {
            /* "... either let x, y, z be real and nonnegative and at most one of
             * them be 0, ..." */
            let mut flag = false;
            if xyzreal_nonneg_atmost1z {
                flag ^= true;
            }
            /* "... or else let two of the variables x, y, z be nonzero and
             * conjugate complex with phase less in magnitude than pi and the
             * third variable be real and nonnegative." */
            if argcheck::r1conj2(x, y, z)
                || argcheck::r1conj2(y, z, x)
                || argcheck::r1conj2(z, x, y)
            {
                flag ^= true;
            }
            return flag;
        }
        false
    }

    /* Cauchy principal value dispatcher
     * Ref[1], Eq. 2.26 */
    pub(super) fn rj_cpv_dispatch(x: f64, y: f64, z: f64, p: f64, rerr: f64) -> (ExitStatus, f64) {
        /* Retry with principal value evaluation, valid for reals. */
        let mut xct1 = [x, y, -p, z];
        let mut xct2 = [0.0; 4];

        let xy = x * y;
        xct2[3] = 1.0 - p / z;
        let pn = (nsum2(&xct1, 3) - xy / xct1[3]) / xct2[3];

        let mut status = ExitStatus::Success;
        let mut status_tmp;
        (status_tmp, xct2[0]) = rj(x, y, z, pn, rerr, false);
        if is_horrible(status_tmp) {
            return (status_tmp, 0.0);
        } else if is_troublesome(status_tmp) {
            status = status_tmp;
        }

        (status_tmp, xct2[1]) = rf(x, y, z, rerr);
        if is_horrible(status_tmp) {
            return (status_tmp, 0.0);
        } else if is_troublesome(status_tmp) {
            status = status_tmp;
        }

        let pq = pn * xct1[2];
        let xypq = xy + pq;
        (status_tmp, xct2[2]) = rc(xypq, pq, rerr);
        if is_horrible(status_tmp) {
            return (status_tmp, 0.0);
        } else if is_troublesome(status_tmp) {
            status = status_tmp;
        }
        xct1[0] = pn - z;
        xct1[1] = -3.0;
        xct1[2] = 3.0 * (xy * z / xypq).sqrt();

        let mut tmpres = ndot2(&xct1, &xct2, 3);
        tmpres /= z - p;
        let res = tmpres;

        (status, res)
    }

    /* The notation follows the one used in the DLMF:
     * https://dlmf.nist.gov/19.27.i
     */
    #[derive(Default)]
    pub(super) struct AsymConfig {
        pub(super) a: f64,
        pub(super) b: f64,
        pub(super) c: f64,
        pub(super) f: f64,
        pub(super) g: f64,
        pub(super) h: f64,
    }

    #[allow(clippy::manual_midpoint)]
    pub(super) fn rj_asym_conf(x: f64, y: f64, z: f64, p: f64, conf: &mut AsymConfig) -> AsymFlag {
        let mut t;

        t = z / p;
        if asymp_zero(t) {
            conf.c = (x + y + z) / 3.0;
            return AsymFlag::Hugep;
        }

        /* this bound is sharp. RJ in this case behaves with logarithmic
         * singularity as p -> +0 */
        if abs_close_zero(p) || ((x != 0.0) && asymp_zero(p / x)) {
            conf.f = (x * y * z).sqrt();
            // Unlike in xsf, the square roots are multiplied if the computed `f` isn't normal,
            // e.g. R_J(2^-540, 2^-540, 2^250, 2^-550) was 6.3e126 (not 8.8e125).
            if !conf.f.is_normal() {
                conf.f = x.sqrt() * y.sqrt() * z.sqrt();
            }
            return AsymFlag::Tinyp;
        }

        /* XXX until RJ(0, y, z, p) is implemented, this is useless */
        /*
        t = (*x) / fmin((*y), (*p));
        if ( ( 0.0 < (*x) && (*x) <= 1e-26 ) || ASYMP_ZERO(t) )
        {
        *h = sqrt((*y) * (*z));
        if ( (*h) / (*p) + 0.5 * ((*y) + (*z)) / (*h) <= sqrt((*h) / (*x)) )
        {
            return asymp_tinyx;
        }
        }
        */

        t = y / z.min(p);
        if ((0.0 < y) && (y <= 1e-26)) || asymp_zero(t) {
            /* bound fairly sharp even if p is large */
            conf.a = 0.5 * (x + y);
            conf.g = (x * y).sqrt();
            if (conf.a / z + conf.a / p) * (p / conf.a).ln().abs() <= 1.0 {
                return AsymFlag::Tinyy;
            }
        }

        if (x != 0.0) && asymp_zero(z.max(p) / x) {
            /* bound might not be sharp if x + 2p much larger than (yz) ** 2,
             * but this is unlikely to be true anyway. */
            return AsymFlag::Hugey;
        }

        if (z != 0.0) && asymp_zero(y.max(p) / z) {
            conf.b = 0.5 * (x + y);
            conf.h = (x * y).sqrt();
            /* when bounds are sharp */
            if (z / (conf.b + conf.h)).ln().abs() <= z.sqrt() {
                return AsymFlag::Hugez;
            }
        }

        AsymFlag::Nothing
    }
}

/* Prevent division by zero due to underflow in atan(sqrt(z)) / sqrt(z) and
 * square root of negative number in the real context. */
#[inline]
pub(super) fn safe_atan_sqrt_div_cplx(z: Complex64) -> Complex64 {
    if z.too_small() {
        return Complex64::from(1.0);
    }
    let s = csqrt(z);
    catan(s).quot(s)
}

#[inline]
pub(super) fn safe_atan_sqrt_div(x: f64) -> f64 {
    if argcheck::too_small(x) {
        return 1.0;
    }
    let s;
    if x < 0.0 {
        s = (-x).sqrt();
        return atanh(s) / s;
    }
    s = x.sqrt();
    s.atan() / s
}

#[allow(clippy::too_many_lines)]
pub(super) fn rj<T: Num>(x: T, y: T, z: T, p: T, rerr: f64, noasymp: bool) -> (ExitStatus, T) {
    use rjimpl::AsymFlag;

    let mut classify = rjimpl::ArgCases::default();
    let mut cct1 = [T::from(0.0); 6];
    let mut cct2 = [T::from(0.0); 6];

    let mut status = ExitStatus::Success;
    if argcheck::invalid_rerr(rerr, 1.0e-4) {
        return (ExitStatus::BadRerr, T::NAN);
    }

    // Unlike in xsf, huge or tiny arguments are rescaled (homogeneity), see `common::scale_factor`,
    // but not downscaled below 2^-254, so that products such as `beta` below don't underflow.
    if let Some(f) = common::scale_factor(&[x, y, z, p], 2.0_f64.powi(-254)) {
        let (status, res) = rj(x * f * f, y * f * f, z * f * f, p * f * f, rerr, noasymp);
        return (status, res * f * f * f);
    }

    /* Put the symmetric arguments in the order of real parts. */
    cct1[0] = x;
    cct1[1] = y;
    cct1[2] = z;
    sort_by(&mut cct1[..3], rjimpl::rcmp);
    let mut xm = cct1[0];
    let mut ym = cct1[1];
    let mut zm = cct1[2];
    if !rjimpl::good_args(xm, ym, zm, p, &mut classify) {
        if classify.good_infinity {
            return (ExitStatus::Success, T::from(0.0));
        }

        if classify.retry_caupv {
            /* Retry with principal value evaluation, valid for reals. */
            let xr = xm.real();
            let yr = ym.real();
            let zr = zm.real();
            let pr = p.real();
            let (status, tmpres) = rjimpl::rj_cpv_dispatch(xr, yr, zr, pr, rerr);
            let res = if is_horrible(status) {
                T::NAN
            } else {
                T::from(tmpres)
            };
            return (status, res);
        } else if classify.hit_pole {
            return (ExitStatus::Singular, T::HUGE);
        }
        return (ExitStatus::BadArgs, T::NAN);
    }

    if classify.maybe_asymp && !noasymp {
        /* might be dealt with by asymptotic expansion of real-arg RJ */
        let mut tmpres = 0.0;
        let mut config = rjimpl::AsymConfig::default();
        let xr = xm.real();
        let yr = ym.real();
        let zr = zm.real();
        let pr = p.real();
        let cres = rjimpl::rj_asym_conf(xr, yr, zr, pr, &mut config);
        match cres {
            AsymFlag::Nothing => {},
            AsymFlag::Hugep => {
                /* Ref[2], Eq. 19.27.11 <https://dlmf.nist.gov/19.27.E11> */
                (status, tmpres) = rf(xr, yr, zr, rerr * 0.5);
                tmpres = 3.0 * (tmpres - 0.5 * constants::PI / pr.sqrt()) / pr;
            },
            AsymFlag::Tinyp => {
                /* Ref[6], Table 2 and Theorem 11 (Eq. 4.12) */
                let mut xct1 = [0.0; 3];
                let mut xct2 = [0.0; 3];
                let r = rerr * 0.5;
                xct1[1] = xr.sqrt();
                xct2[0] = xct1[1];
                xct1[2] = yr.sqrt();
                xct2[1] = xct1[2];
                xct1[0] = zr.sqrt();
                xct2[2] = xct1[0];
                let lamt = dot2(&xct1, &xct2);
                let xct3 = [pr; 3];
                let mut alpha = dot2(&xct1, &xct3) + config.f;
                alpha *= alpha;
                let mut beta = pr + lamt;
                beta = beta * beta * pr;
                // Unlike in xsf, NaN is returned if alpha or beta overflow, instead of dropping
                // the R_C term, so that e.g. R_J(2^-995, 2^692, 2^913, 2^-579) isn't 0.
                if !(alpha.is_finite() && beta.is_finite()) {
                    return (ExitStatus::Overflow, T::NAN);
                }
                let status_tmp;
                (status_tmp, xct2[0]) = rc(alpha, beta, r);
                (status, xct2[1]) = rj(xr + lamt, yr + lamt, zr + lamt, pr + lamt, r, false);
                if status_tmp != ExitStatus::Success {
                    status = status_tmp;
                }
                xct1[0] = 3.0;
                xct1[1] = 2.0;
                tmpres = ndot2(&xct1, &xct2, 2);
            },
            AsymFlag::Hugey => {
                /* Ref[2], Eq. 19.27.14 <https://dlmf.nist.gov/19.27.E14> */
                let r = rerr / 3.0;
                tmpres = 1.0 / (yr * zr).sqrt();
                let (status_tmp, t1) = rc(xr, pr, r);
                let t2;
                (status, t2) = rg(0.0, yr, zr, r);
                if status_tmp != ExitStatus::Success {
                    status = status_tmp;
                }
                tmpres *= 3.0 * t1 - 2.0 * t2 * tmpres;
            },
            AsymFlag::Tinyx => {
                /* Ref[2], Eq. 18.27.15 <https://dlmf.nist.gov/19.27.E15> */
                (status, tmpres) = rj(0.0, yr, zr, pr, rerr, false);
                tmpres -= 3.0 * xr.sqrt() / (config.h * pr);
            },
            AsymFlag::Tinyy => {
                /* Ref[5], Proposition J3
                 * The bounds given in the ref is better than the one in
                 * Ref[6]. */
                let mut tx;
                (status, tx) = rc(1.0, pr / zr, rerr);
                tmpres = (8.0 * zr / (config.a + config.g)).ln() - 2.0 * tx;
                tx = (2.0 * pr / (config.a + config.g)).ln() / pr;
                let rt_h = tx * config.a * (1.0 + 0.5 * pr / zr) / (1.0 - config.a / pr);
                tx = 1.5 / (zr.sqrt() * pr);
                tmpres *= tx;
                tmpres += rt_h * tx;
            },
            AsymFlag::Hugez => {
                /* Try to switch between Ref[5], Propositino J6 and direct
                 * computation without asymptotics */
                let mut tt = config.h + pr;
                tt *= tt;
                let tm = config.b + config.h;
                (status, tmpres) = rc(tt, 2.0 * tm * pr, rerr * 0.5);
                tt = (2.0 * zr / tm).ln();
                let rt_h = 0.75 * (2.0 * pr * tmpres / zr - tt / (zr - config.h));
                tt = 3.0 / zr.sqrt();
                tmpres *= tt;
                /* Also try if direct computation gives reasonable approx. */
                let (status_direct, tmp_direct) = rj(xr, yr, zr, pr, rerr * 0.5, true);
                if (tmp_direct > tmpres) || (tmp_direct < tmpres + rt_h * tt) {
                    tmpres += tt * 0.3 * rt_h;
                } else {
                    tmpres = tmp_direct;
                    status = status_direct;
                }
            },
        } /* end of switch */
        if cres != AsymFlag::Nothing {
            return (status, T::from(tmpres));
        }
    } /* end of condition on whether one should try asymptotic approx. */

    cct1[3] = p;
    cct1[4] = p;
    let mut am = nsum2(&cct1, 5) / 5.0;
    let delta = (p - xm).prod(p - ym).prod(p - zm);
    let mut xxm = am - xm;
    let mut yym = am - ym;
    let mut zzm = am - zm;
    let mut fterm =
        common::max(&[xxm.abs(), yym.abs(), zzm.abs(), (am - p).abs()]) / ocrt(rerr / 5.0);
    /* The algorithm is a direct implementation of Ref[3], Theorem A.1., where
     * a purely arithmetic recurrence relation is used instead of evaluating
     * special functions in the inner loop. */

    /* m = 0; */
    let mut d4m = 1.0;
    let mut pm = p;
    let mut prm;
    let mut lam;
    let mut dm;

    macro_rules! ellint_rj_calc {
        () => {
            prm = pm.sqrt();
            cct1[0] = xm.sqrt();
            cct2[1] = cct1[0];
            cct1[1] = ym.sqrt();
            cct2[2] = cct1[1];
            cct1[2] = zm.sqrt();
            cct2[0] = cct1[2];
            lam = ndot2(&cct1, &cct2, 3);
            dm = (prm + cct1[0]).prod(prm + cct1[1]).prod(prm + cct1[2]);
        };
    }

    macro_rules! ellint_rj_updt {
        () => {
            am = (am + lam) * 0.25;
            xm = (xm + lam) * 0.25;
            ym = (ym + lam) * 0.25;
            zm = (zm + lam) * 0.25;
            pm = (pm + lam) * 0.25;
            xxm = xxm * 0.25;
            yym = yym * 0.25;
            zzm = zzm * 0.25;
            d4m *= 0.25;
            fterm *= 0.25;
        };
    }

    ellint_rj_calc!();
    let mut sm = dm * 0.5;
    /* next */
    ellint_rj_updt!();

    let mut m = 1;
    loop {
        let a_am = am.abs();
        if !(a_am <= fterm
            || a_am <= common::max(&[xxm.abs(), yym.abs(), zzm.abs(), (am - pm).abs()]))
        {
            break;
        }
        if m > config::MAX_ITER {
            status = ExitStatus::NIter;
            break;
        }
        let rm = sm.prod(((delta * d4m).quot(sm.prod(sm)) + 1.0).sqrt() + 1.0);
        ellint_rj_calc!();
        sm = ((rm.prod(dm) - delta * (d4m * d4m)) * 0.5).quot(dm + rm * d4m);

        /* next */
        ellint_rj_updt!();
        m += 1;
    }
    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    cct1[3] = pm;
    cct1[4] = pm;
    am = nsum2(&cct1, 5) / 5.0;
    xxm = xxm.quot(am);
    yym = yym.quot(am);
    zzm = zzm.quot(am);
    /* pp = -0.5 * (xxm + yym + zzm) */
    cct1[0] = xxm;
    cct2[2] = xxm;
    cct1[1] = yym;
    cct2[0] = yym;
    cct1[2] = zzm;
    cct2[1] = zzm;
    let pp = nsum2(&cct1, 3) * -0.5;
    cct1[3] = pp;
    cct2[3] = pp;
    cct1[3] = cct1[3] * -3.0;
    let pp2 = pp.prod(pp);
    let xyz = yym.prod(zzm).prod(xxm);
    /* Prepare the "elementary" terms E_n as in Eqs. (31-32) in Ref. [2] */
    /* e2 = xxm * yym + zzm * xxm + yym * zzm - pp2 * 3.0 */
    let e2 = ndot2(&cct1, &cct2, 4);
    /* e3 = xyz + 2.0 * pp * (e2 + 2.0 * pp2) */
    let e3 = xyz + (pp * 2.0).prod(e2 + pp2 * 2.0);
    /* e4 = (2.0 * xyz + (e2 + 3.0 * pp2) * pp) * pp */
    let e4 = (xyz * 2.0 + (e2 + pp2 * 3.0).prod(pp)).prod(pp);
    let e5 = xyz.prod(pp2);
    /* tmp = d4m * pow(sqrt(Am), -3) */
    let mut t = am.sqrt();
    let mut tmp = T::from(d4m).quot(t.prod(t).prod(t));
    /* Evaluate the 7th-degree expansion using the E_n terms, following
     * Eq. 19.36.2 of [1], https://dlmf.nist.gov/19.36#E2
     * The order of expansion is higher than that in Eq. (32) of Ref. [2]. */
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
    t = dot2(&cct1, &cct2) / constants::RDJ_DENOM + 1.0;
    tmp = tmp.prod(t);
    t = (delta * d4m).quot(sm.prod(sm));
    /* Use the library atan function, following Ref[4], Sec. 1.6.1. */
    tmp = tmp + (t.safe_atan_sqrt_div() * 3.0).quot(sm);

    (status, tmp)
}
