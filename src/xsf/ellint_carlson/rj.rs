//! Translated into pure Rust from `xsf/ellint_carlson/_rj.hh` (scipy/xsf@39019a4).
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

use core::cell::Cell;
use core::f64::consts::FRAC_PI_2;

use super::arithmetic::{comp_horner, dot2, ndot2, nsum2, ocrt};
use super::common::{ExitStatus, config, constants, is_horrible, is_troublesome, sort_by};
use super::rc::rc;
use super::rf::rf;
use super::rg::rg;
use super::typing::Num;
use super::{ELLINT_CARLSON_RERR, argcheck, common};
use crate::utils::{SQRT_EPSILON, pow2};
use crate::xsf::c64::{C64, c64};
use crate::xsf::complex::atanh;

mod rjimpl {
    use super::{Cell, pow2, rc, rf, rj_nested};
    use super::{ExitStatus, Num, argcheck, config, is_horrible, is_troublesome, ndot2, nsum2};

    // Unlike in xsf, also if `r` underflows to 0 (it is a ratio of positive numbers)
    #[inline]
    pub(super) fn asymp_zero(r: f64) -> bool {
        (0.0..=config::ASYM_ZERO_UL).contains(&r)
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
        // Unlike in xsf, subnormal arguments are not treated as 0 here and below
        let zero = T::from(0.0);
        classify.hit_pole = x == zero && y == zero && z.ph_good() && p != zero;
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
        let x0 = x == zero;
        let y0 = y == zero;
        let z0 = z == zero;
        if (pr > 0.0) && (xr >= 0.0) && !((y0 && (x0 || z0)) || (x0 && z0)) {
            return true;
        }
        /* "Alternatively, if p != 0 and |ph p| < pi ..." */
        if p != zero && p.ph_good() {
            /* "... either let x, y, z be real and nonnegative and at most one of
             * them be 0, ..." */
            /* "... or else let two of the variables x, y, z be nonzero and
             * conjugate complex with phase less in magnitude than pi and the
             * third variable be real and nonnegative." */
            return xyzreal_nonneg_atmost1z
                || argcheck::r1conj2(x, y, z)
                || argcheck::r1conj2(y, z, x)
                || argcheck::r1conj2(z, x, y);
        }
        false
    }

    /* Cauchy principal value dispatcher
     * Ref[1], Eq. 2.26 */
    pub(super) fn rj_cpv_dispatch(
        x: f64,
        y: f64,
        z: f64,
        p: f64,
        rerr: f64,
        left: &Cell<Option<u32>>,
    ) -> (ExitStatus, f64) {
        /* Retry with principal value evaluation, valid for reals. */
        let mut xct1 = [x, y, -p, z];
        let mut xct2 = [0.0; 4];

        let xy = x * y;
        xct2[3] = 1.0 - p / z;
        let pn = (nsum2(&xct1, 3) - xy / xct1[3]) / xct2[3];

        let mut status = ExitStatus::Success;
        let mut status_tmp;
        (status_tmp, xct2[0]) = rj_nested(x, y, z, pn, rerr, false, left);
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
        // Unlike in xsf, in closed form if `x y z` or `pn q` under- or overflows
        let closed = !((xy * z).is_normal() && pq.is_normal());
        (status_tmp, xct2[2]) = if closed {
            (ExitStatus::Success, rc_term(x, y, z, pn, -p))
        } else {
            rc(xypq, pq, rerr)
        };
        if is_horrible(status_tmp) {
            return (status_tmp, 0.0);
        } else if is_troublesome(status_tmp) {
            status = status_tmp;
        }
        xct1[0] = pn - z;
        xct1[1] = -3.0;
        xct1[2] = 3.0 * if closed { 1.0 } else { (xy * z / xypq).sqrt() };

        let mut tmpres = ndot2(&xct1, &xct2, 3);
        // Unlike in xsf, NaN if the terms cancel to less than half of the digits, e.g. for tiny
        // `x`, `y` and `p` (instead of a wrong value)
        let size: f64 = (0..3).map(|i| (xct1[i] * xct2[i]).abs()).sum();
        if size > const { pow2(26) } * tmpres.abs() {
            return (ExitStatus::PrecLoss, f64::NAN);
        }
        tmpres /= z - p;
        let res = tmpres;

        (status, res)
    }

    /// `sqrt(x y z / (x y + pn q)) R_C(x y + pn q, pn q)` of `rj_cpv_dispatch` for `x <= y`, i.e.
    /// `sqrt(z) asinh(t) / sqrt(x y + pn q)` with `t = sqrt(x y / (pn q))`
    /// ([DLMF 19.2.19](https://dlmf.nist.gov/19.2.E19)), from the square roots so that it can't
    /// spuriously overflow
    pub(super) fn rc_term(x: f64, y: f64, z: f64, pn: f64, q: f64) -> f64 {
        if x == 0.0 {
            return 0.0;
        }
        let [sx, sy, sp, sq] = [x, y, pn, q].map(f64::sqrt);
        // `sx / sp <= 1`, as `pn - x = (z - x) (y + q) / (z + q) >= 0`
        let t = (sx / sp) * (sy / sq);
        // `n / (a b)`, dividing by the larger one first
        let div = |n: f64, a: f64, b: f64| n / a.max(b) / a.min(b);
        if t <= 1.0 {
            div(z.sqrt() * (t.asinh() / t.hypot(1.0)), sp, sq)
        } else {
            // asinh(t) = ln(2 t) to double precision if `t` overflows
            let asinh = if t.is_finite() {
                t.asinh()
            } else {
                core::f64::consts::LN_2 + (sx / sp).ln() + (sy.ln() - sq.ln())
            };
            div(z.sqrt() * (asinh / t.recip().hypot(1.0)), sx, sy)
        }
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

    /// `sqrt(x * y * ...)`
    ///
    /// Unlike in xsf (scipy/xsf#376), as `sqrt(x) sqrt(y) ...` if `x * y`, `x * y * z`, ...
    /// underflows or overflows.
    fn sqrt_prod(args: &[f64]) -> f64 {
        args[1..]
            .iter()
            .try_fold(args[0], |prod, &a| Some(prod * a).filter(|p| p.is_normal()))
            .map_or_else(|| args.iter().copied().map(f64::sqrt).product(), f64::sqrt)
    }

    /// `ln(c x / y)`
    ///
    /// Unlike in xsf (scipy/xsf#376), as `ln(c) + ln(x) - ln(y)` if `c x / y` overflows.
    pub(super) fn ln_ratio(c: f64, x: f64, y: f64) -> f64 {
        let r = c * x / y;
        if r.is_finite() {
            r.ln()
        } else {
            c.ln() + x.ln() - y.ln()
        }
    }

    /// `R_C(a^2, p b^2)` for when `a^2` or `p b^2` under- or overflows, where
    /// `a = p S + sqrt(x y z)`, `b = p + sqrt(x y) + sqrt(y z) + sqrt(z x)`, and `S` is the sum
    /// of the square roots `[sz, sx, sy]`:
    ///
    /// - R_C(a^2, p b^2) = (s / b) R_C((s a / b)^2, s^2 p), by homogeneity
    /// - s = 2^64 keeps s^2 p normal for subnormal p
    /// - a / b = p (S / b) + sqrt(x) sqrt(y) (sqrt(z) / b) can't overflow, as p <= b and
    ///   sqrt(y z) <= b
    pub(super) fn rc_rescaled(
        [sz, sx, sy]: [f64; 3],
        p: f64,
        b: f64,
        rerr: f64,
    ) -> (ExitStatus, f64) {
        let s = const { pow2(64) };
        let t = (p * ((sx + sy + sz) / b) + sx * (sy * (sz / b))) * s;
        let (status, res) = rc(t * t, s * s * p, rerr);
        (status, res * s / b)
    }

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
            conf.f = sqrt_prod(&[x, y, z]);
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
            conf.g = sqrt_prod(&[x, y]);
            // Unlike in xsf, only if the relative error bound is below the rounding error (not 1),
            // e.g. not for `y = 6e-29, z = 1e-21, p = 3e-9`, where the error is about `a / z`. And
            // also if `a` underflows to 0, where the bound tends to 0.
            if conf.a == 0.0
                || (conf.a / z + conf.a / p) * ln_ratio(1.0, p, conf.a).abs() <= f64::EPSILON
            {
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
            conf.h = sqrt_prod(&[x, y]);
            /* when bounds are sharp */
            if ln_ratio(1.0, z, conf.b + conf.h).abs() <= z.sqrt() {
                return AsymFlag::Hugez;
            }
        }

        AsymFlag::Nothing
    }
}

/* Prevent division by zero due to underflow in atan(sqrt(z)) / sqrt(z) and
 * square root of negative number in the real context. */
#[inline]
pub(super) fn safe_atan_sqrt_div_cplx(z: C64) -> C64 {
    if z.too_small() {
        return C64::ONE;
    }
    let s = z.sqrt();
    s.atan() / s
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

// The Carlson integrals give RJ = 3 RF/p with relative remainder
// at most pi/(2 cos(arg(p)/2)) * sqrt(max(x,y,z)/|p|).
#[cold]
fn rj_far_p<T: Num>(xyz: [T; 3], p: T) -> Option<T> {
    let max_p = p.max_norm();
    let real_p = p.real() > 0.0 && p.imag() == 0.0;
    let balanced_p = p.real().abs().min(p.imag().abs()) >= 0.5 * max_p;
    // For near-real p, each projection of 1/(p+t) differs relatively from
    // that of 1/p by at most 2t/(Re(p)+t).
    let near_real_p = p.real() > 0.0 && !real_p && p.imag().abs() < 0.5 * p.real();
    // For balanced p, cos(arg(p)/2) > 1/5 and each projection of 1/p
    // is at least 1/sqrt(5) of its magnitude.
    let spatial_range = 0.0..=pow2(-120) * max_p;
    if !p.is_finite()
        || !(real_p || balanced_p || near_real_p)
        || xyz
            .iter()
            .any(|a| a.imag() != 0.0 || !spatial_range.contains(&a.real()))
    {
        return None;
    }
    let args = [xyz[0], xyz[1], xyz[2], p];
    if !common::needs_duplication(&args, realness(args).1) {
        return None;
    }
    let [x, y, z] = xyz.map(T::real);
    let (status, f) = rf(x, y, z, ELLINT_CARLSON_RERR * 0.5);
    let (man, exp) = common::frexp(p);
    let value = if near_real_p {
        // Unlike in xsf, keep the smaller projection's exponent through division.
        let (im, im_exp) = common::frexp(p - p.real());
        let den = man.real().mul_add(man.real(), man.imag() * man.imag());
        let scaled = 3.0 * f / den;
        T::from(scaled * man.real()).ldexp(-exp) - (im * scaled).ldexp(im_exp - 2 * exp)
    } else {
        (T::from(3.0 * f) / man).ldexp(-exp)
    };
    let clear_normal = |x: f64| x.is_finite() && x.abs() >= 2.0 * f64::MIN_POSITIVE;
    (status == ExitStatus::Success
        && clear_normal(value.real())
        && (real_p || clear_normal(value.imag())))
    .then_some(value)
}

/// `R_J(x, y, z, p) = R_J(x', y', z', p') / 4 + 3 R_C(alpha^2, beta^2)` of the duplicated
/// arguments ([DLMF 19.26.22](https://dlmf.nist.gov/19.26.E22)), with
/// `R_C(alpha^2, beta^2) = 2 R_C(1, 2 beta / d) / d` by duplication
/// ([DLMF 19.26.18](https://dlmf.nist.gov/19.26.E18)) and
/// `d = alpha + beta = (sqrt(p) + sqrt(x)) (sqrt(p) + sqrt(y)) (sqrt(p) + sqrt(z))`
///
/// The arguments are scaled by `h^2` first (`R_J` is homogeneous of degree -3/2), i.e. quartered
/// for `h = 1 / 2`.
fn rj_duplicated<T: Num>(
    x: T,
    y: T,
    z: T,
    p: T,
    h: f64,
    rerr: f64,
    noasymp: bool,
    left: &Cell<Option<u32>>,
) -> (ExitStatus, T) {
    // `sqrt(a) h`, exact for subnormal `a` if `h = 1 / 2`
    let sq = [x, y, z].map(|a| a.sqrt() * h);
    let sp = p.sqrt() * h;
    let xyz4 = [x, y, z].map(|a| a * (h * h));
    let p4 = p * (h * h);
    let lam = common::lambda([x, y, z], sq);
    let xyz = common::add_lambda(xyz4, sq, lam);
    let mut pd = p4 + lam;
    // Like in the loop, from the nearest shifted argument if it cancels, but only the parts that
    // cancel themselves: in the other one, the shift can lose a tiny part of `p`, and with it the
    // side of the branch cut that `p` is on
    if pd.abs() < p4.abs() {
        let pl = common::add_lambda_p(p4, xyz4, xyz);
        let cancels = |part: fn(T) -> f64| part(pd).abs() < part(p4).abs();
        let re = if cancels(T::real) {
            pl.real()
        } else {
            pd.real()
        };
        let imag_source = if cancels(T::imag) { pl } else { pd };
        pd = imag_source + -imag_source.real() + re;
    }
    // Unlike in xsf, not quartered (`h = 1`) if that rounds a subnormal `p - a` and the square
    // root sum of `p` and `a` cancels, which is formed from it, or `p + lam` cancels, which is
    // formed from that of the nearest `a` (any `a` counts here). Likewise if one of these cancels
    // and it leaves the duplicated arguments with only subnormal imaginary parts, which count as 0.
    let rounded = |i: usize| {
        let d = p4 - xyz4[i];
        d.too_small() && d * 4.0 != p - [x, y, z][i]
    };
    let cancels = |i: usize| pd.abs() < p4.abs() || (sp + sq[i]).abs() < sp.abs();
    let dup = [xyz[0], xyz[1], xyz[2], pd];
    let lost_imag = realness(dup).0 && dup.iter().any(|a| a.imag() != 0.0);
    if h < 1.0 && (0..3).any(|i| (lost_imag || rounded(i)) && cancels(i)) {
        return rj_duplicated(x, y, z, p, 1.0, rerr, noasymp, left);
    }
    // NaN if the unquartered duplicated arguments overflow
    if h >= 1.0 && !(xyz.iter().all(|a| a.is_finite()) && pd.is_finite()) {
        return (ExitStatus::PrecLoss, T::NAN);
    }
    // `d h^3 = d1 d23 2^de` and `2 beta / d` with all factors as `frexp`, so that no partial
    // product or quotient under- or overflows (near the branch cut, a factor of `d` and `p + lam`
    // can be tiny next to a huge `sqrt(p)`)
    let factor = |i: usize| common::sqrt_sum_frexp(p4, xyz4[i], sp, sq[i]);
    let [(d1, e1), (d2, e2), (d3, e3)] = core::array::from_fn(factor);
    let (d23, de) = (d2 * d3, e1 + e2 + e3);
    let [(spm, spe), (pdm, pde)] = [sp, pd].map(common::frexp);
    let u = (spm / d1 * (pdm / d23)).ldexp(spe + pde - de + 1);
    // Unlike in xsf, NaN if `p + lam` cancels and `u` isn't normal, i.e. under- or overflows
    if pd.abs() < p4.abs() && !u.max_norm().is_normal() {
        return (ExitStatus::PrecLoss, T::NAN);
    }
    let (status_rc, res_rc) = rc(T::from(1.0), u, rerr);
    let (mut status, res) = rj_nested(xyz[0], xyz[1], xyz[2], pd, rerr, noasymp, left);
    // NaN if rounding made two of the duplicated arguments 0, i.e. a pole
    if status == ExitStatus::Singular {
        return (ExitStatus::PrecLoss, T::NAN);
    }
    if status_rc != ExitStatus::Success {
        status = status_rc;
    }
    let c = 2.0 * h * h * h;
    (status, res * c + (res_rc * (3.0 * c) / d1 / d23).ldexp(-de))
}

/// Whether `a` is below `sqrt(f64::MIN_POSITIVE)`, the bound for `delta` and `s^2` of the loop:
/// step `m` needs `delta 4^-2m`, which underflows long before `delta` does
fn underflows_in_loop<T: Num>(a: T) -> bool {
    a.max_norm() < f64::MIN_POSITIVE.sqrt()
}

/// The most evaluations that a call of [`rj`] nests, i.e. its rescalings, duplication steps and
/// real fallbacks: a hard bound on the recursion. A huge `p` next to tiny `x, y, z` takes 384
/// steps, which quarter it from at most 2^1024 to 2^256.
const MAX_NESTED: u32 = 512;

pub(super) fn rj<T: Num>(x: T, y: T, z: T, p: T) -> T {
    // Unlike in xsf, use the large-p expansion before repeated duplication.
    // Apply the component range check at the original argument scale.
    if p.max_norm() > const { pow2(256) } {
        if let Some(value) = rj_far_p([x, y, z], p) {
            return value;
        }
    }
    rj_nested(x, y, z, p, ELLINT_CARLSON_RERR, false, &Cell::new(Some(MAX_NESTED))).1
}

pub(super) fn rj_cplx(x: C64, y: C64, z: C64, p: C64) -> C64 {
    // Retain the stock result when both components are already within
    // the accuracy tolerance of the cut tests.
    const AGREEMENT_TOL: f64 = 1e-13;

    let old = rj(x, y, z, p);
    if !(p.re > 0.0 && (x.re < 0.0 || y.re < 0.0 || z.re < 0.0)) {
        return old;
    }
    let mut xyz = [x, y, z];
    sort_by(&mut xyz, rjimpl::rcmp);
    // Unlike in xsf, recover the same-lip coalescence limit or an exact repeat.
    rj_cut_confluent(xyz, p, old)
        .filter(|value| {
            let agrees = |old: f64, new: f64| (old - new).abs() <= AGREEMENT_TOL * new.abs();
            !(agrees(old.re, value.re) && agrees(old.im, value.im))
        })
        .unwrap_or(old)
}

fn rj_cut_confluent([x, y, z]: [C64; 3], p: C64, old: C64) -> Option<C64> {
    fn rc_value<T: Num>(x: T, y: T) -> Option<T> {
        let (status, value) = rc(x, y, ELLINT_CARLSON_RERR);
        (status == ExitStatus::Success).then_some(value)
    }

    let a = -x.re;
    let denominator = p + a;
    let cut_distance = z.re + a;
    // The caller sorts the spatial arguments and has already found a negative real part.
    let same_lip = x.re == y.re
        && x.im.abs().min(y.im.abs()) > 0.0
        && x.im.is_sign_negative() == y.im.is_sign_negative();
    let finite_geometry = z.re > 0.0
        && z.is_real()
        && [x, y, z, p, denominator]
            .iter()
            .all(|value| value.is_finite())
        && cut_distance.is_finite();
    if !(same_lip && finite_geometry) {
        return None;
    }

    let offset = x.im.abs() + y.im.abs();
    // The contour bound |RJ| <= 3*pi/(a*sqrt(Re p)) excludes larger offsets.
    let near_cut = 0.5 * offset <= f64::EPSILON * a;
    let gap = (x.im - y.im).abs();
    // The quadratic error bound cannot pass for wider gaps: |RJ| <= 3*pi/(a*sqrt(Re p)).
    let try_repeated_formula = old.re.is_nan()
        && old.im.is_nan()
        && x.im.is_normal()
        && y.im.is_normal()
        && (x == y || (!near_cut && gap <= a * SQRT_EPSILON));
    if !(near_cut || try_repeated_formula) {
        return None;
    }

    let right = rc_value(z, p)?;
    let project = |left: C64, denominator: C64| {
        let value = (left - right) / denominator * 3.0;
        let scale = 3.0 * ((left.norm() + right.norm()) / denominator.max_norm());
        // Allow at most four bits of cancellation in either component.
        (value.is_finite()
            && value.max_norm() <= 0.5 * f64::MAX
            && scale / 16.0 <= value.re.abs().min(value.im.abs()))
        .then_some(value)
    };

    let limit = || {
        if !near_cut {
            return None;
        }
        // At the shared lip RJ(-a,-a,z,p) = 3 (RC(z,-a) - RC(z,p))/(p+a).
        // The negative real RC returns a principal value, so restore the lip residue.
        let re = rc_value(z.re, -a)?;
        let im = -(FRAC_PI_2 / cut_distance.sqrt()).copysign(x.im);
        let value = project(c64(re, im), denominator)?;
        // The rotated contour bounds the omitted offsets by 8*offset/(a² sqrt(Re p)).
        // The floor keeps accepted components normal in the relative check.
        let perturbation = [a, a, p.re.sqrt()]
            .into_iter()
            .fold(8.0 * offset, |error, divisor| (error / divisor).max(f64::MIN_POSITIVE));
        (perturbation <= f64::EPSILON * value.re.abs().min(value.im.abs())).then_some(value)
    };
    limit().or_else(|| {
        if !try_repeated_formula {
            return None;
        }
        // Unlike in xsf, use DLMF 19.20.8 on exact repeats: RJ(x,x,z,p)
        // = 3 (RC(z,x) - RC(z,p))/(p-x).
        let repeated = |w: C64| {
            let denominator = p - w;
            if !denominator.is_finite() {
                return None;
            }
            project(rc_value(z, w)?, denominator).filter(|v| v.re.is_normal() && v.im.is_normal())
        };
        let first = repeated(x)?;
        if x == y {
            return Some(first);
        }
        // Unlike in xsf, average the two repeated-argument integrals.
        // On the rotated contour the error is at most 16*(Im x-Im y)^2/(a^3 sqrt(Re p)).
        let separation = gap / a;
        let perturbation = [a, p.re.sqrt()]
            .into_iter()
            .fold((16.0 * separation * separation).max(f64::MIN_POSITIVE), |error, divisor| {
                (error / divisor).max(f64::MIN_POSITIVE)
            });
        let second = repeated(y)?;
        let value = 0.5 * (first + second);
        (first.re.is_sign_negative() == second.re.is_sign_negative()
            && first.im.is_sign_negative() == second.im.is_sign_negative()
            && perturbation <= f64::EPSILON * value.re.abs().min(value.im.abs()))
        .then_some(value)
    })
}

/// [`rj`] with the number of evaluations that are `left` of [`MAX_NESTED`], `None` if it ran out
fn rj_nested<T: Num>(
    x: T,
    y: T,
    z: T,
    p: T,
    rerr: f64,
    noasymp: bool,
    left: &Cell<Option<u32>>,
) -> (ExitStatus, T) {
    if argcheck::invalid_rerr(rerr, 1.0e-4) {
        return (ExitStatus::BadRerr, T::NAN);
    }
    // Unlike in xsf, NaN if none is left
    left.set(left.get().and_then(|n| n.checked_sub(1)));
    if left.get().is_none() {
        return (ExitStatus::PrecLoss, T::NAN);
    }

    // Unlike in xsf (scipy/xsf#376), huge or tiny arguments are rescaled, and if negligible
    // imaginary parts prevent that and the result is not finite, rescaled with those flushed, and
    // then also with negligible real parts flushed if none is negative and the arguments are in
    // the domain (flushed, two of `x, y, z` can become conjugate)
    let args = [x, y, z, p];
    let (real, min) = realness(args);
    // A real part that is flushed by downscaling becomes `MIN_POSITIVE`, not 0: on the imaginary
    // axis, `p` can leave the domain, and `p + lam` cancel to 0, e.g. for
    // `x, y, z, p = i, i, i, 1 - 15i` scaled (the other real parts are at least that already)
    let rescaled = |f: f64| {
        let [x, y, z, p] = args.map(|a| {
            let s = a * f * f;
            if f < 1.0 && a.real() > 0.0 && s.real() < f64::MIN_POSITIVE {
                s + -s.real() + f64::MIN_POSITIVE
            } else {
                s
            }
        });
        let (status, res) = rj_nested(x, y, z, p, rerr, noasymp, left);
        (status, res * f * f * f)
    };
    if let Some(f) = common::scale_factor(&args, min) {
        return rescaled(f);
    }
    // Unlike in xsf, complex arguments are also upscaled if `delta` underflows in the loop: its
    // `R_C` terms need the precision of `delta`. Not if that takes them out of the domain, where a
    // subnormal imaginary part counts as 0.
    if !real && underflows_in_loop((p - x) * (p - y) * (p - z)) {
        let out = common::upscale_factor(&args).map(rescaled);
        if let Some(out) = out.filter(|out| out.0 != ExitStatus::BadArgs) {
            return out;
        }
    }
    let in_domain = || {
        let mut xyz = [x, y, z];
        sort_by(&mut xyz, rjimpl::rcmp);
        rjimpl::good_args(xyz[0], xyz[1], xyz[2], p, &mut rjimpl::ArgCases::default())
    };
    let mut out = rj_unscaled(x, y, z, p, rerr, noasymp, left);
    for flush_re in [false, true] {
        if out.1.is_finite() || (flush_re && !in_domain()) {
            break;
        }
        out = common::relaxed_scale_factor(&args, flush_re, min).map_or(out, rescaled);
    }
    // NaN as well if a fallback made up for an evaluation that ran out, like the asymptotic one
    // for a huge `z`
    left.get().map_or((ExitStatus::PrecLoss, T::NAN), |_| out)
}

/// Whether all arguments are real (imaginary parts 0 or subnormal), and the smallest argument size
/// of the rescaling: 2^-511 so
/// that products of two arguments don't underflow, or 2^-254 if complex (wider spans are
/// inaccurate there)
fn realness<T: Num>(args: [T; 4]) -> (bool, f64) {
    let real = args.iter().all(|a| argcheck::too_small(a.imag()));
    let min = if real {
        f64::MIN_POSITIVE.sqrt()
    } else {
        const { pow2(-254) }
    };
    (real, min)
}

/// The scale of the `|Im p|` that `rj_unscaled` raises a tinier one to for real `lam` and
/// `Re p < 0`, relative to the distance from `p` to `0, x, y, z`: chosen to keep `4^-m Im p` from
/// underflowing in the steps, while the result changes only by about as much
const MIN_IM_P: f64 = pow2(-60);

fn rj_unscaled<T: Num>(
    x: T,
    y: T,
    z: T,
    p: T,
    rerr: f64,
    noasymp: bool,
    left: &Cell<Option<u32>>,
) -> (ExitStatus, T) {
    use rjimpl::AsymFlag;

    let (real, min) = realness([x, y, z, p]);
    let mut classify = rjimpl::ArgCases::default();
    let mut cct1 = [T::from(0.0); 6];
    let mut cct2 = [T::from(0.0); 6];

    let mut status = ExitStatus::Success;

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
            let (status, tmpres) = rjimpl::rj_cpv_dispatch(xr, yr, zr, pr, rerr, left);
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

    // Unlike in xsf, keep `lam` real for real `x, y, z` or an exact conjugate pair and a real
    // third, which the steps keep: its rounding error can outweigh a tiny `Im p`
    let pair = |a: T, b: T, c: T| a.real() == b.real() && a.imag() == -b.imag() && c.imag() == 0.0;
    let real_lam = [xm, ym, zm].iter().all(|a| a.imag() == 0.0)
        || pair(xm, ym, zm)
        || pair(ym, zm, xm)
        || pair(zm, xm, ym);
    // Unlike in xsf, then a tiny `Im p` left of the imaginary axis is raised to 1 to 2 times
    // `MIN_IM_P` times the distance from `p` to `0, x, y, z` (in the max norm, so that the one to
    // `0` is finite): it only selects the side of the cut, which `p + lam` then keeps in every
    // loop step
    let tiny_im = |dist: f64| p.imag().abs() < MIN_IM_P * dist;
    let p = if real_lam && p.real() < 0.0 && p.imag() != 0.0 && tiny_im(p.max_norm()) {
        let dist = [xm, ym, zm]
            .map(|a| (p - a).max_norm())
            .into_iter()
            .fold(p.max_norm(), f64::min);
        if tiny_im(dist) {
            common::frexp(p - p.real()).0 * (2.0 * MIN_IM_P * dist) + p.real()
        } else {
            p
        }
    } else {
        p
    };

    // Unlike in xsf (scipy/xsf#376), if a tiny argument prevents the downscaling of a huge one, a
    // duplication step first raises the tiny one
    if common::needs_duplication(&[xm, ym, zm, p], min) {
        return rj_duplicated(xm, ym, zm, p, 0.5, rerr, noasymp, left);
    }

    // Unlike in xsf, also for complex input with `p` tiny relative to `x, y, z`, where the loop's
    // `R_C` terms lose accuracy (`1 + delta / s^2` cancels to about `sqrt(p / x)`), as the step
    // gets `2 beta / d` from `sqrt(p)` directly. Real input takes the "tiny p" form below, which
    // is the same identity.
    //
    // (not if rounding moves the shifted arguments out of the domain, e.g. a real part to below 0)
    let tiny_p = |a: &T| p.max_norm() <= config::ASYM_ZERO_UL * a.max_norm();
    if !real && [xm, ym, zm].iter().all(tiny_p) {
        let (status, res) = rj_duplicated(xm, ym, zm, p, 0.5, rerr, noasymp, left);
        if status != ExitStatus::BadArgs {
            return (status, res);
        }
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
                // Unlike in xsf, with the next terms `6 R_G / p^2 - 9 pi c / (4 p^(5/2))` of the
                // expansion of `1 / (t + p)` in the integral of `R_J`, with `c` the mean of
                // `x, y, z` (relative error `O((c / p)^2)` instead of `O(c / p)`). The `R_G` term
                // is at most `2 c / p` relative to `R_F`, as `R_G <= sqrt(c)` and
                // `R_F >= 1 / sqrt(c)`, so it is skipped if that is below 2^-55.
                if config.c / pr > const { pow2(-56) } {
                    let (status_g, res_g) = rg(xr, yr, zr, rerr * 0.5);
                    if status_g != ExitStatus::Success {
                        status = status_g;
                    }
                    tmpres += 2.0 * res_g / pr;
                }
                let pi_term = 0.5 * constants::PI / pr.sqrt() * (1.0 + 1.5 * config.c / pr);
                tmpres = 3.0 * (tmpres - pi_term) / pr;
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
                let status_tmp;
                (status_tmp, xct2[0]) = if alpha.is_normal() && beta.is_normal() {
                    rc(alpha, beta, r)
                } else {
                    // Unlike in xsf (scipy/xsf#376), rescaled if `alpha` or `beta` under- or
                    // overflows.
                    rjimpl::rc_rescaled(xct1, pr, pr + lamt, r)
                };
                (status, xct2[1]) =
                    rj_nested(xr + lamt, yr + lamt, zr + lamt, pr + lamt, r, false, left);
                if status_tmp != ExitStatus::Success {
                    status = status_tmp;
                }
                xct1[0] = 3.0;
                xct1[1] = 2.0;
                tmpres = ndot2(&xct1, &xct2, 2);
                // Unlike in xsf, inf (not NaN) if a term overflows
                if tmpres.is_nan() {
                    tmpres = 3.0 * xct2[0] + 2.0 * xct2[1];
                }
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
                (status, tmpres) = rj_nested(0.0, yr, zr, pr, rerr, false, left);
                tmpres -= 3.0 * xr.sqrt() / (config.h * pr);
            },
            AsymFlag::Tinyy => {
                /* Ref[5], Proposition J3
                 * The bounds given in the ref is better than the one in
                 * Ref[6]. */
                // Unlike in xsf, ln(c v / (a + g)) as 2 ln(sqrt(2 c) sqrt(v) / (sqrt(x) + sqrt(y)))
                // if `a + g = (sqrt(x) + sqrt(y))^2 / 2` is subnormal
                let ag = config.a + config.g;
                let ln_ag = |c: f64, v: f64| {
                    if ag.is_normal() {
                        rjimpl::ln_ratio(c, v, ag)
                    } else {
                        2.0 * rjimpl::ln_ratio((2.0 * c).sqrt(), v.sqrt(), xr.sqrt() + yr.sqrt())
                    }
                };
                let mut tx;
                (status, tx) = rc(1.0, pr / zr, rerr);
                tmpres = ln_ag(8.0, zr) - 2.0 * tx;
                tx = ln_ag(2.0, pr) / pr;
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
                tt = rjimpl::ln_ratio(2.0, zr, tm);
                let rt_h = 0.75 * (2.0 * pr * tmpres / zr - tt / (zr - config.h));
                tt = 3.0 / zr.sqrt();
                tmpres *= tt;
                /* Also try if direct computation gives reasonable approx. */
                let (status_direct, tmp_direct) = rj_nested(xr, yr, zr, pr, rerr * 0.5, true, left);
                let asym = tmpres + tt * 0.3 * rt_h;
                // Unlike in xsf (scipy/xsf#376), `asym` is also used if `tmp_direct` is NaN, unless
                // `tmpres` is 0 or subnormal (e.g. if `rc` overflows), or `asym` is not finite
                if (tmp_direct > tmpres)
                    || (tmp_direct < tmpres + rt_h * tt)
                    || (tmp_direct.is_nan() && tmpres.is_normal() && asym.is_finite())
                {
                    tmpres = asym;
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
    let delta = (p - xm) * (p - ym) * (p - zm);
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
    // Whether `diff`, a sum that cancels `a`, lost more than 4 bits
    let cancels = |a: T, diff: T| a.max_norm() > diff.max_norm() * const { pow2(4) };
    // Whether the arguments are in the right half-plane, where `Re lam >= 0`
    let right = xm.real() >= 0.0 && p.real() > 0.0;
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
            if real_lam {
                lam = T::from(lam.real());
            }
            dm = common::sqrt_sum(pm, xm, prm, cct1[0])
                * common::sqrt_sum(pm, ym, prm, cct1[1])
                * common::sqrt_sum(pm, zm, prm, cct1[2]);
        };
    }

    macro_rules! ellint_rj_updt {
        () => {
            am = (am + lam) * 0.25;
            let xyz = [xm, ym, zm];
            let xyzl = common::add_lambda(xyz, [cct1[0], cct1[1], cct1[2]], lam);
            let mut pl = pm + lam;
            // Unlike in xsf (scipy/xsf#409), `p + lam` as `(a + lam) + (p - a)` if it cancels, with
            // `a` the closest of `x, y, z`.
            if pl.abs() < pm.abs() {
                pl = common::add_lambda_p(pm, xyz, xyzl);
                // Unlike in xsf, restore `Im(p + lam) = 4^-m Im p` for real `lam`: a complex `a`
                // can round away a tiny `Im p`
                if real_lam {
                    pl = (p - p.real()) * d4m + pl.real();
                }
                // Unlike in xsf, in the right half-plane the real part of `p + lam` is at least its
                // lower bound `4^-m Re p`: rounded to the left of the imaginary axis, `p + lam` and
                // the `R_C` terms can be on different branches.
                let re = p.real() * d4m;
                if right && pl.real() < re {
                    pl = pl + -pl.real() + re;
                }
            }
            [xm, ym, zm] = xyzl.map(|a| a * 0.25);
            pm = pl * 0.25;
            xxm = xxm * 0.25;
            yym = yym * 0.25;
            zzm = zzm * 0.25;
            d4m *= 0.25;
            fterm *= 0.25;
        };
    }

    ellint_rj_calc!();
    let mut sm = dm * 0.5;
    // Unlike in xsf, the step also comes first for complex arguments if both `delta` and `s^2`
    // still underflow in the loop after upscaling: its `R_C` terms need both, the step's neither
    //
    // (not if the shifted arguments leave the domain: a subnormal imaginary part counts as 0, but
    // can make theirs normal)
    if !real && underflows_in_loop(delta) && underflows_in_loop(sm * sm) {
        let (status, res) = rj_duplicated(xm, ym, zm, p, 0.5, rerr, noasymp, left);
        if status != ExitStatus::BadArgs {
            return (status, res);
        }
    }
    /* next */
    ellint_rj_updt!();
    // `s^2 + delta / 4 = s beta`, see below
    let mut w2 = sm * prm * pm * 4.0;
    // `w r` and `d + r 4^-m` of the previous step, for `u` below (`u = 2 beta / d` at first)
    let mut wr = T::from(0.5);
    let mut den = T::from(1.0);

    // `R_C(1, u) = atan(sqrt(t)) / sqrt(t)` with `t = u - 1`
    //
    // Unlike in xsf, evaluate `R_C` from `u` if `1 + t` cancels, e.g. after a step where `p + lam`
    // cancels, as `R_C(1, u) ~ ln(4 / u) / 2` there, but only where `u` keeps the side of the cut:
    // in the right half-plane, and for real `lam`, where `Im(p + lam) = 4^-m Im p` (see
    // `MIN_IM_P`). NaN if `u` underflows.
    //
    // (elsewhere `p + lam` can lose a tiny `Im p`, and `u` its side of the cut)
    let rc1 = |t: T, u: T| {
        if !((right || real_lam) && cancels(T::from(1.0), u)) {
            t.safe_atan_sqrt_div()
        } else if u.max_norm().is_normal() {
            rc(T::from(1.0), u, rerr).1
        } else {
            T::NAN
        }
    };
    // The sum of the `R_C` terms of the steps before `m`, `3 atan(sqrt(t)) / (sqrt(t) s)` with
    // `t = delta 4^-m / s^2`, as in xsf's last step
    let s_term = |s: T, d4m: f64| (delta * d4m / (s * s)).safe_atan_sqrt_div() * 3.0 / s;
    // Unlike in xsf, the `R_C` terms are summed directly once the update of `s` cancels by more
    // than 4 bits, e.g. for `x, y, z << p`, where `s` is lost while `r d ~ delta 4^-2m`. The term
    // of a step is the `R_C` term of `rj_duplicated`, `6 4^-m R_C(1, 1 + e) / d` with
    // `e = delta 4^-3m / d^2` and `1 + e = 2 beta / d`, with the `p + lam` of `beta` as the last
    // factor (it can be subnormal).
    let step_rc_term = |d: T, d4m: f64, prm: T, pl: T| {
        let e = delta * d4m / d * (T::from(d4m) / d) * d4m;
        rc1(e, prm / d * (pl * 2.0)) * (d4m * 6.0) / d
    };
    let mut rc_sum = None;

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
        if let Some(sum) = rc_sum {
            ellint_rj_calc!();
            ellint_rj_updt!();
            rc_sum = Some(sum + step_rc_term(dm, d4m * 4.0, prm, pm * 4.0));
            m += 1;
            continue;
        }
        let ss = sm * sm;
        let t = delta * d4m / ss;
        let u = t + 1.0;
        // Unlike in xsf, `w = s sqrt(u)` as `sqrt(w2)` if `|u| < 2`, since `u` cancels as `u -> 0`,
        // e.g. for `x, y << p << z`. The branch is that of `sqrt(u)`, i.e. `Re(w / s) >= 0` (`s` is
        // inaccurate if `|u| >> 1`). But not if `w2` differs from `s^2 u` by more than the
        // cancellation can explain, e.g. after a step where `d + r 4^-m` cancels, or if `s^2`
        // underflows (in the max norm, which is cheaper).
        let by_w2 = u.real() * u.real() + u.imag() * u.imag() < 4.0
            && ss.max_norm().is_normal()
            && (w2 - ss * u).max_norm() <= ss.max_norm() * 0.5;
        let (wm, rm) = if by_w2 {
            let w = w2.sqrt();
            let w = if w.real() * sm.real() + w.imag() * sm.imag() < 0.0 {
                -w
            } else {
                w
            };
            (w, sm + w)
        } else {
            let su = u.sqrt();
            (sm * su, sm * (su + 1.0))
        };
        // `u = w2 / s^2` if that is the accurate one, with the `p + lam` of `beta` as the last
        // factor (`w2` can be subnormal for a normal `u`)
        let pl = pm * 4.0;
        let uw =
            (by_w2 && cancels(T::from(1.0), u)).then(|| wr / ss * (dm / den) * (prm / den) * pl);
        ellint_rj_calc!();
        den = dm + rm * d4m;
        let rd = rm * dm;
        let num = rd - delta * (d4m * d4m);
        // The `R_C` terms before this step
        let before =
            cancels(rd, num).then(|| uw.map_or_else(|| s_term(sm, d4m), |u| rc1(t, u) * 3.0 / sm));
        sm = num * 0.5 / den;

        /* next */
        ellint_rj_updt!();
        // `w2 = s^2 + delta 4^-m` (of the next step) as `w r d beta / (d + r 4^-m)^2`, without
        // cancellation, from `d^2 + delta 4^-3m = 2 d beta` with `beta = sqrt(p) (p + lam)`
        // (https://dlmf.nist.gov/19.26.E23), where `p + lam = 4 p'`
        wr = wm * rm;
        w2 = wr * (dm / den) * (prm * pm * 4.0 / den);
        rc_sum = before.map(|sum| sum + step_rc_term(dm, d4m * 4.0, prm, pm * 4.0));
        m += 1;
    }
    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    cct1[3] = pm;
    cct1[4] = pm;
    am = nsum2(&cct1, 5) / 5.0;
    xxm = xxm / am;
    yym = yym / am;
    zzm = zzm / am;
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
    let pp2 = pp * pp;
    let xyz = yym * zzm * xxm;
    /* Prepare the "elementary" terms E_n as in Eqs. (31-32) in Ref. [2] */
    /* e2 = xxm * yym + zzm * xxm + yym * zzm - pp2 * 3.0 */
    let e2 = ndot2(&cct1, &cct2, 4);
    /* e3 = xyz + 2.0 * pp * (e2 + 2.0 * pp2) */
    let e3 = xyz + pp * 2.0 * (e2 + pp2 * 2.0);
    /* e4 = (2.0 * xyz + (e2 + 3.0 * pp2) * pp) * pp */
    let e4 = (xyz * 2.0 + (e2 + pp2 * 3.0) * pp) * pp;
    let e5 = xyz * pp2;
    /* tmp = d4m * pow(sqrt(Am), -3) */
    let mut t = am.sqrt();
    let mut tmp = T::from(d4m) / (t * t * t);
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
    tmp = tmp * t;
    /* Use the library atan function, following Ref[4], Sec. 1.6.1. */
    tmp = tmp + rc_sum.unwrap_or_else(|| s_term(sm, d4m));

    (status, tmp)
}

#[cfg(test)]
mod tests {
    use super::rjimpl::rc_term;
    use super::{C64, Cell, ExitStatus, Num, pow2, rj_far_p, rj_nested};
    use crate::utils::DENORM_MIN;
    use crate::xsf::c64::c64;

    /// Whether `rj` nests `n` evaluations: it is finite with that many, and NaN with fewer
    fn nests<T: Num>(x: T, y: T, z: T, p: T, n: u32) -> bool {
        let rj = |n| rj_nested(x, y, z, p, 5e-16, false, &Cell::new(Some(n)));
        let ran_out = |n| matches!(rj(n), (ExitStatus::PrecLoss, res) if res.abs().is_nan());
        rj(n).1.is_finite() && (0..n).all(ran_out)
    }

    #[test]
    fn test_rj_nested() {
        // 4 duplication steps quarter a huge `p` next to tiny `x, y, z` to 2^256
        let t = f64::MIN_POSITIVE;
        assert!(nests(t, t, t, pow2(264), 5));
        // the asymptotic fallback for a huge `z` doesn't make up for its direct evaluation
        assert!(nests(1.0, 1.0, 1e50, 1.0, 2));
        // the duplication step and the upscaling alternate
        let [x, y, z] = [0.0, DENORM_MIN, pow2(253)].map(C64::from);
        assert!(nests(x, y, z, c64(pow2(-1020), pow2(-1020)), 5));
    }

    #[test]
    fn test_rj_far_p_declines_subnormal() {
        // Arb rounds both components to -minsub at the original scale.
        let a = c64(1.273_270_726_927_572_3e106, 0.0);
        let p = c64(-pow2(900), pow2(900));
        assert!(rj_far_p([C64::ZERO, a, a], p).is_none());
    }

    #[test]
    fn test_rc_term() {
        let actual = [
            // `t` overflows
            rc_term(1.0, 1e300, 1e308, 1e300, 5e-324),
            // `sqrt(z / (pn q))` overflows
            rc_term(1.0, 1.0, 1e300, 2.0, 5e-324),
            rc_term(1e-200, 1e-200, 2.0, 3e-200, 1e-210),
            rc_term(0.5, 4.0, 9.0, 5.0, 1e-300),
        ];
        // Arb
        let expected = [
            3_729_131.831_412_505_7,     // 3729131.8314125056872
            3.725_666_095_509_706e152,   // 3.7256660955097061359e152
            1.648_515_727_704_121_5e201, // 1.6485157277041215091e201
            733.176_604_211_949,         // 733.17660421194900136
        ];
        for (a, e) in actual.into_iter().zip(expected) {
            assert!((a - e).abs() <= 1e-15 * e, "{a} != {e}");
        }
        assert_eq!(rc_term(0.0, 1e308, 1.0, 1e308, 5e-324), 0.0);
    }
}
