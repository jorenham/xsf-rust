//! Inverses of distribution functions with respect to a parameter, which `scipy.special` computes
//! with cdflib or Boost, and noncentral distribution functions from cdflib's series. The roots are
//! found with the bracket search and Chandrupatla solver of [`gdtrib`](crate::gdtrib) to machine
//! precision, whereas cdflib's `dinvr` stops at a relative tolerance of `1e-10` (`1e-8` for the
//! noncentral t). If the probability can't be reached, 0 or inf is returned, whichever limit is
//! closest, instead of cdflib's search bound. If it is 0 or 1, the limit where the CDF is smallest
//! or largest is returned, instead of a root where the CDF underflows or rounds to it.

use core::f64::consts::{LN_2, SQRT_2};

use crate::ext::quadpack::quad_unimodal;
use crate::utils::{SQRT_FRAC_2_PI, SQRT_PI, pow2};
use crate::xsf::cephes::dd_real::{ldexp_f64, two_sum};
use crate::xsf::cephes::{
    betainc_xy, fdtr, fdtr_w, igam, igam_exp_scaled, igam_fac, igamc, incbet_cf_factor,
    incbet_ln_pow_asymp, lbeta, ln_a_beta, ln_gamma_term,
};
use crate::xsf::{
    bracket_root_for_cdf_inversion, erfcx, find_positive_root, find_root_chandrupatla, log1pmx,
};

/// Inverse of [`stdtr`](crate::stdtr) with respect to `nu`
///
/// Corresponds to [`scipy.special.stdtridf`][stdtridf].
///
/// [stdtridf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtridf.html
///
/// # See also
/// - [`stdtr`](crate::stdtr): Student's t distribution CDF
/// - [`stdtri`](crate::stdtri): Inverse of [`stdtr`](crate::stdtr) with respect to `x`
#[must_use]
pub fn stdtridf(p: f64, t: f64) -> f64 {
    // stdtr(df, t) is 1/2 if t = 0, and 1 or 0 if t = +-inf, independent of df
    if !((0.0..=1.0).contains(&p) && t != 0.0 && t.is_finite()) {
        return f64::NAN;
    }
    // stdtr(df, -t) = 1 - stdtr(df, t) is used for the upper tail
    let func = cdf_minus_p(p, |df| crate::stdtr(df, t), |df| crate::stdtr(df, -t));
    // stdtr(df, t) increases from 1/2 at df = 0 to ndtr(t) at df = inf if t > 0
    limit_if_p_is_0_or_1(p, t > 0.0).unwrap_or_else(|| find_df_root(func, t > 0.0))
}

/// Inverse of [`fdtr`](crate::fdtr) with respect to `dfd`
///
/// Corresponds to [`scipy.special.fdtridfd`][fdtridfd].
///
/// For `f < 1`, `fdtr` can peak in `dfd`, so that there can be two roots, of which one is returned.
/// The same goes for [`ncfdtridfd`] and [`ncfdtridfn`].
///
/// [fdtridfd]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtridfd.html
///
/// # See also
/// - [`fdtr`](crate::fdtr): F distribution CDF
/// - [`fdtri`](crate::fdtri): Inverse of [`fdtr`](crate::fdtr) with respect to `x`
#[must_use]
pub fn fdtridfd(dfn: f64, p: f64, f: f64) -> f64 {
    // fdtr(dfn, dfd, f) is 0 or 1 if f = 0 or inf, independent of dfd
    if !(dfn > 0.0 && (0.0..=1.0).contains(&p) && is_pos_finite(f)) {
        return f64::NAN;
    }
    let func = cdf_minus_p(
        p,
        |dfd| crate::fdtr(dfn, dfd, f),
        |dfd| crate::fdtrc(dfn, dfd, f),
    );
    // fdtr(dfn, dfd, f) goes from 0 at dfd = 0 to chdtr(dfn, dfn f) at dfd = inf
    limit_if_p_is_0_or_1(p, true).unwrap_or_else(|| find_df_root(func, true))
}

/// Inverse of [`bdtr`](crate::bdtr) with respect to `k`
///
/// Corresponds to [`scipy.special.bdtrik`][bdtrik].
///
/// The CDF $\I_{1-p}(n - k, k + 1)$ is continuous in `k`, so that `k` is real. If `y` lies below
/// the CDF at $k = 0$, i.e. $(1 - p)^n$, 0 is returned.
///
/// [bdtrik]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrik.html
///
/// # See also
/// - [`bdtr`](crate::bdtr): Binomial distribution function
/// - [`bdtri`](crate::bdtri): Inverse of [`bdtr`](crate::bdtr) with respect to `p`
/// - [`bdtrin`]: Inverse of [`bdtr`](crate::bdtr) with respect to `n`
#[must_use]
pub fn bdtrik(y: f64, n: f64, p: f64) -> f64 {
    if !((0.0..=1.0).contains(&y) && (0.0..f64::INFINITY).contains(&n) && (0.0..=1.0).contains(&p))
    {
        return f64::NAN;
    }
    // the CDF is 0 below k = n if p = 1
    if y > 0.0 && (y == 1.0 || p == 1.0) {
        return n;
    }
    // The CDF `betaincc(k + 1, n - k, p)` increases from (1 - p)^n at k = 0 to 1 at k = n, beyond
    // which it stays 1.
    invert_betainc(y, true, n, 1.0 - p, p, |k| ((n - k).max(0.0), k + 1.0))
}

/// Inverse of [`bdtr`](crate::bdtr) with respect to `n`
///
/// Corresponds to [`scipy.special.bdtrin`][bdtrin].
///
/// The CDF $\I_{1-p}(n - k, k + 1)$ is continuous in `n`, so that `n` is real.
///
/// [bdtrin]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrin.html
///
/// # See also
/// - [`bdtr`](crate::bdtr): Binomial distribution function
/// - [`bdtri`](crate::bdtri): Inverse of [`bdtr`](crate::bdtr) with respect to `p`
/// - [`bdtrik`]: Inverse of [`bdtr`](crate::bdtr) with respect to `k`
#[must_use]
pub fn bdtrin(k: f64, y: f64, p: f64) -> f64 {
    // As in SciPy, NaN if p is 0 or 1, where the CDF is 1 or 0 for all n > k, unless y is 0 or 1,
    // which give the limits n = inf and n = k.
    let p_is_valid = p > 0.0 && p < 1.0 || (p == 0.0 || p == 1.0) && (y == 0.0 || y == 1.0);
    if !((0.0..f64::INFINITY).contains(&k) && (0.0..=1.0).contains(&y) && p_is_valid) {
        return f64::NAN;
    }
    // The CDF `betaincc(k + 1, n - k, p)` decreases from 1 at n = k to 0 at n = inf. It is solved
    // for `n - k`, which keeps its precision if `n` is close to `k`.
    k + invert_betainc(y, false, f64::MAX, 1.0 - p, p, |nk| (nk, k + 1.0))
}

/// Inverse of [`betainc`](crate::betainc) with respect to `a`
///
/// Corresponds to [`scipy.special.btdtria`][btdtria].
///
/// Unlike in SciPy, NaN is returned for `x > 1` even if `p` is 0 or 1, instead of the limit.
///
/// [btdtria]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.btdtria.html
///
/// # See also
/// - [`betainc`](crate::betainc): Regularized incomplete beta function
/// - [`betaincinv`](crate::betaincinv): Inverse of [`betainc`](crate::betainc) with respect to `x`
/// - [`btdtrib`]: Inverse of [`betainc`](crate::betainc) with respect to `b`
#[must_use]
pub fn btdtria(p: f64, b: f64, x: f64) -> f64 {
    if !btdtr_is_valid(p, b, x) {
        return f64::NAN;
    }
    // the CDF decreases from 1 at a = 0 to 0 at a = inf
    invert_betainc(p, false, f64::MAX, x, 1.0 - x, |a| (a, b))
}

/// Inverse of [`betainc`](crate::betainc) with respect to `b`
///
/// Corresponds to [`scipy.special.btdtrib`][btdtrib].
///
/// Unlike in SciPy, NaN is returned for `x > 1` even if `p` is 0 or 1, instead of the limit.
///
/// [btdtrib]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.btdtrib.html
///
/// # See also
/// - [`betainc`](crate::betainc): Regularized incomplete beta function
/// - [`betaincinv`](crate::betaincinv): Inverse of [`betainc`](crate::betainc) with respect to `x`
/// - [`btdtria`]: Inverse of [`betainc`](crate::betainc) with respect to `a`
#[must_use]
pub fn btdtrib(a: f64, p: f64, x: f64) -> f64 {
    if !btdtr_is_valid(p, a, x) {
        return f64::NAN;
    }
    // the CDF increases from 0 at b = 0 to 1 at b = inf
    invert_betainc(p, true, f64::MAX, x, 1.0 - x, |b| (a, b))
}

/// Whether `p`, the other shape and `x` are valid arguments of [`btdtria`] and [`btdtrib`]. As in
/// SciPy, the limits are also returned for `p` = 0 or 1 if the other shape is inf or `x = 1`, where
/// the CDF doesn't depend on the shape that is solved for.
fn btdtr_is_valid(p: f64, shape: f64, x: f64) -> bool {
    let is_limit = p == 0.0 || p == 1.0;
    (0.0..=1.0).contains(&p)
        && shape > 0.0
        && x > 0.0
        && x <= 1.0
        && (is_limit || shape.is_finite() && x < 1.0)
}

/// Inverse of [`chdtr`](crate::chdtr) with respect to `v`
///
/// Corresponds to [`scipy.special.chdtriv`][chdtriv].
///
/// [chdtriv]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtriv.html
///
/// # See also
/// - [`chdtr`](crate::chdtr): Chi-squared distribution function
/// - [`chdtri`](crate::chdtri): Inverse of [`chdtrc`](crate::chdtrc) with respect to `x`
/// - [`gdtrib`](crate::gdtrib): Inverse of [`gdtr`](crate::gdtr) with respect to `b`
#[must_use]
pub fn chdtriv(p: f64, x: f64) -> f64 {
    if !(x > 0.0) {
        return f64::NAN;
    }
    // chdtr(v, x) = gdtr(1/2, v/2, x) decreases from 1 at v = 0 to 0 at v = inf, also if x = inf
    limit_if_p_is_0_or_1(p, false).unwrap_or_else(|| {
        // as in SciPy, NaN if x/2 underflows to 0, where the CDF would be 0 for all v
        if 0.5 * x == 0.0 {
            return f64::NAN;
        }
        2.0 * crate::gdtrib(0.5, p, x)
    })
}

/// Inverse of [`nbdtr`](crate::nbdtr) with respect to `k`
///
/// Corresponds to [`scipy.special.nbdtrik`][nbdtrik].
///
/// The CDF $\I_p(n, k + 1)$ is continuous in `k`, so that `k` is real. If `y` lies below the CDF at
/// $k = 0$, i.e. $p^n$, 0 is returned.
///
/// [nbdtrik]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtrik.html
///
/// # See also
/// - [`nbdtr`](crate::nbdtr): Negative binomial distribution function
/// - [`nbdtri`](crate::nbdtri): Inverse of [`nbdtr`](crate::nbdtr) with respect to `p`
/// - [`nbdtrin`]: Inverse of [`nbdtr`](crate::nbdtr) with respect to `n`
#[must_use]
pub fn nbdtrik(y: f64, n: f64, p: f64) -> f64 {
    if !((0.0..=1.0).contains(&y) && is_pos_finite(n) && (0.0..=1.0).contains(&p)) {
        return f64::NAN;
    }
    // the CDF increases from p^n at k = 0 to 1 at k = inf
    invert_betainc(y, true, f64::MAX, p, 1.0 - p, |k| (n, k + 1.0))
}

/// Inverse of [`nbdtr`](crate::nbdtr) with respect to `n`
///
/// Corresponds to [`scipy.special.nbdtrin`][nbdtrin].
///
/// The CDF $\I_p(n, k + 1)$ is continuous in `n`, so that `n` is real. Unlike in SciPy, the limits
/// inf and 0 are returned if `y` is 0 or 1, instead of NaN. As in SciPy, NaN is returned if `p` is
/// 0 or 1.
///
/// [nbdtrin]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtrin.html
///
/// # See also
/// - [`nbdtr`](crate::nbdtr): Negative binomial distribution function
/// - [`nbdtri`](crate::nbdtri): Inverse of [`nbdtr`](crate::nbdtr) with respect to `p`
/// - [`nbdtrik`]: Inverse of [`nbdtr`](crate::nbdtr) with respect to `k`
#[must_use]
pub fn nbdtrin(k: f64, y: f64, p: f64) -> f64 {
    // NaN if p is 0 or 1, where the CDF is 0 or 1 for all n > 0
    if !((0.0..f64::INFINITY).contains(&k) && (0.0..=1.0).contains(&y) && p > 0.0 && p < 1.0) {
        return f64::NAN;
    }
    // Unlike in SciPy, the root is found for n, instead of for n + k, which loses n if k >> n.
    // The CDF decreases from 1 at n = 0 to 0 at n = inf.
    invert_betainc(y, false, f64::MAX, p, 1.0 - p, |n| (n, k + 1.0))
}

/// Inverse of [`pdtr`](crate::pdtr) with respect to `k`
///
/// Corresponds to [`scipy.special.pdtrik`][pdtrik].
///
/// The CDF $Q(k + 1, m)$ is continuous in `k`, so that `k` is real. If `p` lies below the CDF at
/// $k = 0$, i.e. $e^{-m}$, 0 is returned. As in SciPy, NaN is returned if $p = 1$.
///
/// [pdtrik]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtrik.html
///
/// # See also
/// - [`pdtr`](crate::pdtr): Poisson distribution function
/// - [`pdtri`](crate::pdtri): Inverse of [`pdtr`](crate::pdtr) with respect to `m`
#[must_use]
pub fn pdtrik(p: f64, m: f64) -> f64 {
    // as in SciPy, 0 if p = 0, also if m = inf
    if !((0.0..1.0).contains(&p) && m >= 0.0 && (m.is_finite() || p == 0.0)) {
        return f64::NAN;
    }
    // the CDF increases from exp(-m) at k = 0 to 1 at k = inf
    invert_cdf(
        p,
        true,
        f64::MAX,
        |k| igamc(k + 1.0, m),
        |k| igam(k + 1.0, m),
    )
}

/// `cdf(x) - p`, which is computed as `(1 - p) - sf(x)` with the survival function `sf` if
/// `p > 1/2`, where `1 - p` is exact, so that the root keeps the precision of the smaller tail.
fn cdf_minus_p(p: f64, cdf: impl Fn(f64) -> f64, sf: impl Fn(f64) -> f64) -> impl Fn(f64) -> f64 {
    let q = 1.0 - p;
    move |x| if p <= q { cdf(x) - p } else { q - sf(x) }
}

/// 0 or inf, wherever a CDF that is `increasing` in its parameter is smallest if `p = 0`, or
/// largest if `p = 1`.
fn limit_if_p_is_0_or_1(p: f64, increasing: bool) -> Option<f64> {
    (p == 0.0 || p == 1.0).then_some(if (p == 0.0) == increasing {
        0.0
    } else {
        f64::INFINITY
    })
}

/// The root in `[MIN_POSITIVE, upper]` of `cdf(x) = p`, with the survival function `sf`, for a CDF
/// that is `increasing` in `x`, or the limit if `p` is 0 or 1
fn invert_cdf(
    p: f64,
    increasing: bool,
    upper: f64,
    cdf: impl Fn(f64) -> f64,
    sf: impl Fn(f64) -> f64,
) -> f64 {
    limit_if_p_is_0_or_1(p, increasing).unwrap_or_else(|| {
        find_positive_root(
            cdf_minus_p(p, cdf, sf),
            f64::MIN_POSITIVE,
            upper,
            increasing,
        )
    })
}

/// [`invert_cdf`] of `betainc(a, b, x)`, with `(a, b) = ab(v)`, given `y = 1 - x`, so that the CDF
/// can also be a complement `betaincc(b, a, y)`, by swapping `x` and `y`
fn invert_betainc(
    p: f64,
    increasing: bool,
    upper: f64,
    x: f64,
    y: f64,
    ab: impl Fn(f64) -> (f64, f64),
) -> f64 {
    let cdf = |v| {
        let (a, b) = ab(v);
        betainc_xy(a, b, x, y)
    };
    let sf = |v| {
        let (a, b) = ab(v);
        betainc_xy(b, a, y, x)
    };
    invert_cdf(p, increasing, upper, cdf, sf)
}

/// The degrees of freedom are searched for in `[8^-DF_K, 8^DF_K] ~ [5e-100, 2e99]`, like in cdflib
/// from `1e-100`, instead of in all of `f64`, which saves ~460 evaluations of the CDF if `p` can't
/// be reached.
const DF_K: i32 = 110;

/// [`find_positive_root`], also searching the other way from the start (1) if `func` doesn't
/// change sign in the given direction, and then around the extremum of `func` with
/// [`find_extremum_root`].
///
/// The CDFs need not be monotonic in the degrees of freedom: e.g. `fdtr(dfn, dfd, f)` can peak in
/// `dfd` if `f < 1`, so that if `p` lies above its limit at `dfd = inf`, there can be two roots.
fn find_df_root(func: impl Fn(f64) -> f64, increasing: bool) -> f64 {
    let (lower, upper) = (8f64.powi(-DF_K), 8f64.powi(DF_K));
    let x = find_positive_root(&func, lower, upper, increasing);
    if x == 0.0 || x == f64::INFINITY {
        let x_other = find_positive_root(&func, lower, upper, !increasing);
        if x_other > 0.0 && x_other < f64::INFINITY {
            return x_other;
        }
        return scan_same_sign(&func)
            .and_then(|points| find_extremum_root(&func, &points, increasing))
            .unwrap_or(x);
    }
    x
}

/// The number of search points below 1 over which [`scan_same_sign`] needs `func` to stay constant
/// to assume that it reached its limit, since the CDFs can also round to a constant before a peak
/// at small `df`, e.g. the noncentral t CDF with `t > 0` between `df = 1` and its peak at
/// `df ~ 1 / nc^2`, which `8^14 > 1e12` spans for `|nc| <= 1e6`.
const PLATEAU_K: i32 = 14;

/// `(k, func(8^k))` at the search points `8^k` of [`find_df_root`] from 1 outwards, until `func`
/// reaches its limits, i.e. stops changing, if `func` has the same sign at all of them.
fn scan_same_sign(func: impl Fn(f64) -> f64) -> Option<Vec<(i32, f64)>> {
    let mut points = vec![(0, func(1.0))];
    for step in [-1, 1] {
        let mut y_prev = points[0].1;
        let mut plateau = 0;
        let mut k = step;
        while (-DF_K..=DF_K).contains(&k) {
            let yk = func(8f64.powi(k));
            points.push((k, yk));
            plateau = if yk == y_prev { plateau + 1 } else { 0 };
            if plateau == if step < 0 { PLATEAU_K } else { 1 } {
                break;
            }
            y_prev = yk;
            k += step;
        }
    }
    points.sort_by_key(|&(k, _)| k);
    let sign = points[0].1.signum();
    points
        .iter()
        .all(|&(_, y)| y != 0.0 && y.signum() == sign)
        .then_some(points)
}

/// Root of `func` next to its extremum, if it has the same sign at all `points` of
/// [`scan_same_sign`], but the opposite sign at an extremum between them. Of the two roots, the
/// one where `func` is `increasing` is returned.
fn find_extremum_root(
    func: impl Fn(f64) -> f64,
    points: &[(i32, f64)],
    increasing: bool,
) -> Option<f64> {
    // `sign * func`, which is positive at the points, and negative at the root's extremum
    let sign = points[0].1.signum();
    let g = |x: f64| sign * func(x);
    let gs: Vec<f64> = points.iter().map(|&(_, y)| sign * y).collect();
    // the smallest local minimum, which is bracketed by its neighbors, unless it's rounding noise
    let below = |i: usize, j: usize| gs[i] < (1.0 - 8.0 * f64::EPSILON) * gs[j];
    let i = (1..gs.len() - 1)
        .filter(|&i| below(i, i - 1) && below(i, i + 1))
        .min_by(|&i, &j| gs[i].total_cmp(&gs[j]))?;
    let (xl, xr) = (8f64.powi(points[i - 1].0), 8f64.powi(points[i + 1].0));

    // Golden-section search for the minimum of `g` in `ln x`, until `g` is negative. 100 steps
    // shrink the interval by `0.618^100 ~ 1e-21`, below the spacing of `ln x`.
    let inv_phi = (5f64.sqrt() - 1.0) / 2.0;
    let (mut a, mut b) = (xl.ln(), xr.ln());
    let (mut c, mut d) = (b - inv_phi * (b - a), a + inv_phi * (b - a));
    let (mut gc, mut gd) = (g(c.exp()), g(d.exp()));
    for _ in 0..100 {
        if !(gc >= 0.0 && gd >= 0.0) {
            break;
        }
        if gc < gd {
            (b, d, gd) = (d, c, gc);
            c = b - inv_phi * (b - a);
            gc = g(c.exp());
        } else {
            (a, c, gc) = (c, d, gd);
            d = a + inv_phi * (b - a);
            gd = g(d.exp());
        }
    }
    let xm = if gc < 0.0 {
        c.exp()
    } else if gd < 0.0 {
        d.exp()
    } else {
        return None;
    };

    // `func` increases from the search points to the extremum if it's negative there
    let (x1, x2) = if increasing == (sign < 0.0) {
        (xl, xm)
    } else {
        (xm, xr)
    };
    // only a relative tolerance, as in `find_positive_root`
    let (x, status) =
        find_root_chandrupatla(&func, x1, x2, func(x1), func(x2), f64::EPSILON, 0.0, 100);
    (status == 0).then_some(x)
}

/// The largest noncentrality for [`cumchn`] and [`cumfnc`]: as in cdflib, where `nc / 2` must fit
/// an `int`, since they sum `~17 sqrt(nc)` terms.
const NC_MAX: f64 = 2.0 * i32::MAX as f64;

/// Whether `x` is positive and finite
fn is_pos_finite(x: f64) -> bool {
    x > 0.0 && x < f64::INFINITY
}

/// The exponent of the largest power of 2 by which [`rescale`] divides at once.
const RESCALE: i32 = 500;

/// `ln(x^a y^b / (a B(a, b)))`, with `y = 1 - x`, i.e. `ln(I_x(a, b) - I_x(a + 1, b))`
/// (A&S 26.5.16), and whether it avoids the cancellation of `lbeta(a, b)`
fn ln_beta_term(a: f64, b: f64, x: f64, y: f64, lnx: f64, lny: f64) -> (f64, bool) {
    // Unlike in cdflib, the limit if `a = dfn / 2` rounds to 0, as `a B(a, b) -> 1`
    if a == 0.0 {
        return (b * lny, false);
    }
    let t = -lbeta(a, b) - a.ln() + a * lnx + b * lny;
    // Unlike in cdflib, without the cancellation of `lbeta(a, b)` against `a ln(x)` or `b ln(y)`
    // for a huge `b` or `a`, like the `I_x(a, b)` that the terms start from, if `x` or `y` is
    // normal (a subnormal one is rounded here)
    let ln_pow = |p, q, v: f64| incbet_ln_pow_asymp(p, q, v, t).filter(|_| v >= f64::MIN_POSITIVE);
    if let Some(p) = ln_pow(a, b, x) {
        (p - a.ln() + b * lny, true)
    } else if let Some(p) = ln_pow(b, a, y) {
        (p - a.ln() + a * lnx, true)
    } else {
        (t, false)
    }
}

/// `(m, k)` with `I_x(a, b) = m 2^k`, where `m ~ 1` and `k < 0` if `I_x(a, b)` underflows, so that
/// the terms of a series that starts there can be carried as multiples of `2^k`
fn betainc_scaled(a: f64, b: f64, x: f64, y: f64, lnx: f64, lny: f64) -> (f64, i32) {
    let ix = betainc_xy(a, b, x, y);
    if ix >= f64::MIN_POSITIVE {
        return (ix, 0);
    }
    let (ln_term, _) = ln_beta_term(a, b, x, y, lnx, lny);
    if !ln_term.is_finite() {
        return (ix, 0);
    }
    // I_x(a, b) = x^a y^b / (a B(a, b)) w, with x below the mean if I_x underflows
    let (w, _) = incbet_cf_factor(a, b, x, y);
    if !(w > 0.0 && w.is_finite()) {
        return (ix, 0);
    }
    scaled_exp(ln_term, w)
}

/// `(m, k)` with `exp(ln_term) w = m 2^k` and `m ~ 1`
fn scaled_exp(ln_term: f64, w: f64) -> (f64, i32) {
    let k = ((ln_term + w.ln()) / LN_2).round() as i32;
    ((ln_term - f64::from(k) * LN_2).exp() * w, k)
}

/// `(m, k)` with `P(a, x) = m 2^k`, as [`betainc_scaled`] for `I_x(a, b)`
fn igam_scaled(a: f64, x: f64) -> (f64, i32) {
    let p = igam(a, x);
    if p >= f64::MIN_POSITIVE {
        return (p, 0);
    }
    igam_exp_scaled(a, x).map_or((p, 0), |(m, ln_f)| scaled_exp(ln_f, m))
}

/// Divides the `values` by powers of 2, adding their exponents to `scale`, while
/// `peak * num / den > 2^RESCALE` and `scale < 0`, so that the scaled terms stay finite and below
/// their unscaled values. Returns the sum of the exponents.
///
/// `num / den` is separate, as it can overflow.
fn rescale(mut peak: f64, num: f64, den: f64, scale: &mut i32, values: &mut [&mut f64]) -> i32 {
    let mut shift = 0;
    while *scale < 0 && peak * num > pow2(RESCALE) * den {
        let k = RESCALE.min(-*scale);
        peak *= pow2(-k);
        for v in values.iter_mut() {
            **v *= pow2(-k);
        }
        *scale += k;
        shift += k;
    }
    shift
}

/// `lambda^k exp(-lambda) / k!` for `k >= 1`
fn poisson_weight(k: f64, lambda: f64) -> f64 {
    // `igam_fac` loses `~|ln(lambda)| eps` if lambda < 1
    if k == 1.0 {
        lambda * (-lambda).exp()
    } else {
        igam_fac(k + 1.0, lambda) / lambda
    }
}

/// CDF of the noncentral chi-square distribution
///
/// Corresponds to [`scipy.special.chndtr`][chndtr].
///
/// [chndtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtr.html
///
/// # See also
/// - [`chdtr`](crate::chdtr): CDF of the (central) chi-square distribution
/// - [`chndtri`]: Inverse of [`chndtr`] with respect to `x`
/// - [`chndtridf`]: Inverse of [`chndtr`] with respect to `df`
/// - [`chndtrinc`]: Inverse of [`chndtr`] with respect to `nc`
#[must_use]
pub fn chndtr(x: f64, df: f64, nc: f64) -> f64 {
    // the edge cases of SciPy
    if !(x >= 0.0 && df >= 0.0 && nc >= 0.0) {
        f64::NAN
    } else if nc == 0.0 && (df == 0.0 || x.is_infinite()) {
        1.0
    } else if nc == 0.0 && df.is_infinite() {
        0.0
    } else if df == 0.0 || df.is_infinite() && nc.is_finite() {
        f64::NAN
    } else if nc.is_infinite() {
        if x.is_infinite() { f64::NAN } else { 0.0 }
    } else if x.is_infinite() {
        1.0
    } else if nc > NC_MAX {
        f64::NAN
    } else {
        // Unlike in SciPy, cdflib's series (A&S 26.4.25), as Boost's returns 0 where its central
        // term underflows, e.g. `chndtr(20, 63.5, 774) = 1.9e-147`.
        let p = cumchn(x, df, nc);
        // the rounding of the sum can exceed 1, but keep NaN
        if p > 1.0 { 1.0 } else { p }
    }
}

/// Inverse of [`chndtr`] with respect to `x`
///
/// Corresponds to [`scipy.special.chndtrix`][chndtrix].
///
/// [chndtrix]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtrix.html
///
/// # See also
/// - [`chndtridf`]: Inverse of [`chndtr`] with respect to `df`
/// - [`chndtrinc`]: Inverse of [`chndtr`] with respect to `nc`
#[must_use]
pub fn chndtri(p: f64, df: f64, nc: f64) -> f64 {
    if !((0.0..=1.0).contains(&p) && is_pos_finite(df) && (0.0..=NC_MAX).contains(&nc)) {
        return f64::NAN;
    }
    limit_if_p_is_0_or_1(p, true).unwrap_or_else(|| {
        // as low as the smallest subnormal, where `cumchn` is accurate too
        find_positive_root(|x| cumchn(x, df, nc) - p, f64::from_bits(1), f64::MAX, true)
    })
}

/// Inverse of [`chndtr`] with respect to `df`
///
/// Corresponds to [`scipy.special.chndtridf`][chndtridf].
///
/// [chndtridf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtridf.html
///
/// # See also
/// - [`chndtri`]: Inverse of [`chndtr`] with respect to `x`
/// - [`chndtrinc`]: Inverse of [`chndtr`] with respect to `nc`
#[must_use]
pub fn chndtridf(x: f64, p: f64, nc: f64) -> f64 {
    if !(is_pos_finite(x) && (0.0..=1.0).contains(&p) && (0.0..=NC_MAX).contains(&nc)) {
        return f64::NAN;
    }
    // the CDF decreases in df
    limit_if_p_is_0_or_1(p, false)
        .unwrap_or_else(|| find_df_root(|df| cumchn(x, df, nc) - p, false))
}

/// Inverse of [`chndtr`] with respect to `nc`
///
/// Corresponds to [`scipy.special.chndtrinc`][chndtrinc].
///
/// [chndtrinc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtrinc.html
///
/// # See also
/// - [`chndtri`]: Inverse of [`chndtr`] with respect to `x`
/// - [`chndtridf`]: Inverse of [`chndtr`] with respect to `df`
#[must_use]
pub fn chndtrinc(x: f64, df: f64, p: f64) -> f64 {
    if !(is_pos_finite(x) && is_pos_finite(df) && (0.0..=1.0).contains(&p)) {
        return f64::NAN;
    }
    // the CDF decreases in nc, from chdtr(df, x) at nc = 0
    limit_if_p_is_0_or_1(p, false).unwrap_or_else(|| {
        find_positive_root(|nc| cumchn(x, df, nc) - p, f64::MIN_POSITIVE, NC_MAX, false)
    })
}

/// CDF of the noncentral chi-square distribution, from cdflib's `cumchn`.
///
/// It sums the series A&S 26.4.25 backward and forward from the term with the largest Poisson
/// weight, with the incomplete gamma functions from the recurrence A&S 6.5.21.
///
/// Unlike in cdflib, as in [`cumfnc`], with `ln_igam_fac` for the first terms of the recurrence,
/// since the `lgam` differences lose `~a ln(a) eps`, and `igam` instead of `cumchi`.
fn cumchn(x: f64, df: f64, pnonc: f64) -> f64 {
    let eps = f64::EPSILON;

    if !(x > 0.0) {
        return 0.0;
    }
    if !(pnonc >= 2.0 * f64::MIN_POSITIVE) {
        return crate::chdtr(df, x);
    }

    let xnonc = pnonc / 2.0;
    // Unlike in cdflib, the first term `exp(-nc / 2) chdtr(df, x)` for a subnormal `x`, as `x / 2`
    // rounds, to `~(1 + nc) x` relative
    if x < 2.0 * f64::MIN_POSITIVE {
        return crate::chdtr(df, x) * (-xnonc).exp();
    }
    let icent = xnonc.floor().max(1.0);
    let centwt = poisson_weight(icent, xnonc);
    let chid2 = x / 2.0;
    let mut adn = df / 2.0 + icent;
    // in units of `2^scale` where `P(adn, chid2)` underflows, as the terms below `icent` need not
    let (mut pdn, mut scale) = igam_scaled(adn, chid2);
    let scale_up = scale;
    let mut pup = pdn;
    let mut ssum = centwt * pdn;

    let mut wt = centwt;
    let mut i = icent;

    // chid2^adn exp(-chid2) / Gamma(adn + 1)
    let mut dnterm = (ln_gamma_term(adn, chid2) - f64::from(scale) * LN_2).exp();

    // as in `cumfnc`, and in cdflib
    loop {
        wt *= i / xnonc;
        i -= 1.0;
        // `dnterm` is scaled via the denominator, as scaling it first can make it subnormal
        let shift = rescale(pdn, adn, chid2, &mut scale, &mut [&mut pdn, &mut ssum]);
        dnterm = dnterm * adn / ldexp_f64(chid2, shift);
        adn -= 1.0;
        pdn += dnterm;
        ssum += wt * pdn;
        if !(wt * pdn > eps * ssum && i > 0.0 && wt >= f64::MIN_POSITIVE) {
            break;
        }
    }

    i = icent;
    wt = centwt;

    let mut aup = df / 2.0 + icent;
    let mut upterm = (ln_gamma_term(aup, chid2) - f64::from(scale) * LN_2).exp();
    pup = ldexp_f64(pup, scale_up - scale);

    loop {
        i += 1.0;
        wt *= xnonc / i;
        pup -= upterm;
        aup += 1.0;
        upterm *= chid2 / aup;
        ssum += wt * pup;
        if !(wt * pup > eps * ssum) {
            break;
        }
    }

    ldexp_f64(ssum, scale)
}

/// CDF of the noncentral F distribution
///
/// Corresponds to [`scipy.special.ncfdtr`][ncfdtr].
///
/// [ncfdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtr.html
///
/// # See also
/// - [`ncfdtri`]: Inverse of [`ncfdtr`] with respect to `f`
/// - [`fdtr`](crate::fdtr): CDF of the (central) F distribution
#[must_use]
pub fn ncfdtr(dfn: f64, dfd: f64, nc: f64, f: f64) -> f64 {
    // as in SciPy, also for infinite parameters
    if f == f64::INFINITY && dfn > 0.0 && dfd > 0.0 && nc >= 0.0 {
        return 1.0;
    }
    if !(is_pos_finite(dfn) && is_pos_finite(dfd) && (0.0..=NC_MAX).contains(&nc) && f >= 0.0) {
        return f64::NAN;
    }
    // Unlike in SciPy, cdflib's series (A&S 26.6.20), as Boost's is inaccurate where its central
    // term underflows, e.g. `ncfdtr(5, 50, 1500, 5) = 1.7e-188`, not `8.2e-161`.
    let p = cumfnc(f, dfn, dfd, nc);
    // the rounding of the sum can exceed 1, but keep NaN
    if p > 1.0 { 1.0 } else { p }
}

/// Inverse of [`ncfdtr`] with respect to `f`
///
/// Corresponds to [`scipy.special.ncfdtri`][ncfdtri].
///
/// [ncfdtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtri.html
///
/// # See also
/// - [`ncfdtridfd`]: Inverse of [`ncfdtr`] with respect to `dfd`
/// - [`ncfdtridfn`]: Inverse of [`ncfdtr`] with respect to `dfn`
/// - [`ncfdtrinc`]: Inverse of [`ncfdtr`] with respect to `nc`
#[must_use]
pub fn ncfdtri(dfn: f64, dfd: f64, nc: f64, p: f64) -> f64 {
    if !(is_pos_finite(dfn)
        && is_pos_finite(dfd)
        && (0.0..=NC_MAX).contains(&nc)
        && (0.0..=1.0).contains(&p))
    {
        return f64::NAN;
    }
    limit_if_p_is_0_or_1(p, true).unwrap_or_else(|| {
        // as low as the smallest subnormal, like `chndtri`
        find_positive_root(
            |f| cumfnc(f, dfn, dfd, nc) - p,
            f64::from_bits(1),
            f64::MAX,
            true,
        )
    })
}

/// CDF of the noncentral F distribution, from cdflib's `cumfnc`.
///
/// It sums the series A&S 26.6.20 forward and backward from the term with the largest Poisson
/// weight, with the incomplete beta functions from the recurrences A&S 26.5.16.
///
/// Unlike in SciPy,
/// - the series is summed to machine precision, instead of to a relative tolerance of `1e-4` or
///   until the sum is below `1e-300`,
/// - `nc` is treated as 0 below `2 MIN_POSITIVE` instead of `1e-10`,
/// - the central Poisson weight is computed with `igam_fac`, since `exp` loses `~nc ln(nc) eps`,
/// - the first terms of the recurrences are computed with `lbeta`, also if `a < 2`, where the
///   `lgam` differences lose `~b ln(b) eps`,
/// - `ln(1 - x)` is computed with `ln_1p`,
/// - subnormal `x` is handled separately, since it rounds,
/// - the beta terms are rescaled by powers of 2 if the central one underflows,
/// - `cumf` and `bratio` are replaced by [`fdtr`](crate::fdtr) and `betainc_xy`.
fn cumfnc(f: f64, dfn: f64, dfd: f64, pnonc: f64) -> f64 {
    let eps = f64::EPSILON;

    if !(f > 0.0) {
        return 0.0;
    }
    // the weights below are subnormal if `nc / 2 < MIN_POSITIVE`, where `nc` changes the CDF by
    // less than `nc / 2` relative
    if !(pnonc >= 2.0 * f64::MIN_POSITIVE) {
        return fdtr(dfn, dfd, f);
    }
    if !(pnonc <= NC_MAX) {
        return f64::NAN;
    }

    let xnonc = pnonc / 2.0;
    // Unlike in cdflib, the first term `exp(-nc / 2) I_x(dfn / 2, dfd / 2)` separately if
    // `dfn / 2 + dfd / 2 < 1 / 2`, as it can exceed the next by `~1 / (dfn / 2 + dfd / 2)`
    // relative, past where the loop below stops, which then starts above it
    let i_min = if dfn + dfd < 1.0 { 1.0 } else { 0.0 };
    let icent = xnonc.floor().max(i_min + 1.0);
    let centwt = poisson_weight(icent, xnonc);
    let prod = dfn * f;
    let dsum = dfd + prod;
    // Unlike in cdflib, a subnormal `x = dfn f / dsum`, which rounds: if `dfn f < eps`, the first
    // term `exp(-nc / 2) fdtr(dfn, dfd, f)`, to `~(1 + nc) dfn f` relative, else (`dfd > 1e291`)
    // the noncentral chi-square CDF of `dfn f`
    if prod / dsum < 2.0 * f64::MIN_POSITIVE {
        if prod >= eps {
            return cumchn(prod, dfn, pnonc);
        }
        return fdtr(dfn, dfd, f) * (-xnonc).exp();
    }
    // `dfd / dsum`, also if `prod` overflows
    let mut yy = if prod.is_finite() {
        dfd / dsum
    } else {
        dfd / dfn / f
    };
    let xx;
    // `ln` of the rounded complement would lose ~eps / min(xx, yy), which `b` or `adn` amplify
    let (lnxx, lnyy);
    if yy > 0.5 {
        // Unlike in cdflib, without the rounding of a subnormal `prod`
        xx = fdtr_w(dfn, dfd, f).0;
        yy = 1.0 - xx;
        (lnxx, lnyy) = (xx.ln(), (-xx).ln_1p());
    } else {
        xx = 1.0 - yy;
        // without the rounding of a subnormal `yy`
        let ln_y = if !prod.is_finite() {
            dfd.ln() - dfn.ln() - f.ln()
        } else if yy < f64::MIN_POSITIVE {
            dfd.ln() - dsum.ln()
        } else {
            yy.ln()
        };
        (lnxx, lnyy) = ((-yy).ln_1p(), ln_y);
    }
    let mut adn = dfn / 2.0 + icent;
    let mut aup = adn;
    // Unlike in cdflib, if `dfd / 2` is subnormal, where it rounds (e.g. to 0), the series for
    // `b = MIN_POSITIVE` scaled by `dfd / (2 b)`, as its terms are `O(b)` (but the first if
    // `dfn / 2 ~ b`, which is separate then)
    let b = (dfd / 2.0).max(f64::MIN_POSITIVE);
    // If `I_x(adn, b)` underflows, the beta terms are carried as multiples of `2^scale`, since the
    // terms below `icent` can still be significant, e.g. if x is tiny.
    // Unlike in cdflib, if `yy` is subnormal or 0, as it rounds: `1 - I_y(b, adn)`, with (DLMF
    // 8.17.7) `I_y(b, adn) = y^b / (b B(adn, b)) F(b, 1 - adn; b + 1; y)`. If `adn > 1 / eps`, its
    // limit `Q(b, z)` as `adn -> inf`, with `z = adn y`, up to `O(1 / adn)`, which avoids
    // cancelling `ln(y)` and `ln(adn)`. Else, or if `z` is subnormal, `z < eps` and the first term
    // suffices: the next is `b (1 - adn) y / (b + 1)` relative, so `< z / 35` relative to `1 - I`.
    // `z` avoids the rounding of `yy`, and the overflow of `adn / dsum` for a tiny `dsum` (where
    // `adn dfd = z dsum` is normal if `z >= eps`).
    let z = if prod.is_finite() {
        let r = adn / dsum;
        if r.is_finite() {
            r * dfd
        } else {
            adn * dfd / dsum
        }
    } else {
        adn / dfn / f * dfd
    };
    let (mut betdn, mut scale) = if yy >= f64::MIN_POSITIVE {
        betainc_scaled(adn, b, xx, yy, lnxx, lnyy)
    } else if adn * eps > 1.0 && z >= f64::MIN_POSITIVE {
        (igamc(b, z), 0)
    } else {
        (-(b * lnyy - ln_a_beta(b, adn)).exp_m1(), 0)
    };
    let scale_up = scale;
    let mut betup = betdn;
    let mut ssum = centwt * betdn;

    let mut xmult = centwt;
    let mut i = icent;

    // x^adn y^b / (adn B(adn, b))
    let (ln_cent, asymp) = ln_beta_term(adn, b, xx, yy, lnxx, lnyy);
    let cent_term = (ln_cent - f64::from(scale) * LN_2).exp();
    let mut dnterm = cent_term;

    // Unlike in cdflib, at least once, as the central term underflows for a tiny `nc` while the
    // next need not, with `>` instead of `>=`, which never ends if `eps * ssum` underflows, and not
    // past the underflow of the Poisson weight, which gets stuck at a subnormal while the scaled
    // `betdn` grows
    loop {
        xmult *= i / xnonc;
        i -= 1.0;
        // Unlike in cdflib, `dfn / 2` last, not `adn - 1`, which loses up to `~icent eps` of it
        // (all if `dfn < eps`), while `I_x(dfn / 2, b)` depends on it relatively
        adn = if i > 0.0 { adn - 1.0 } else { dfn / 2.0 };
        let num = adn + 1.0;
        let shift = rescale(
            betdn,
            num,
            (adn + b) * xx,
            &mut scale,
            &mut [&mut betdn, &mut ssum],
        );
        // Unlike in cdflib, `xx` scaled first, as `(adn + b) xx` can round to a subnormal
        dnterm = dnterm * num / ((adn + b) * ldexp_f64(xx, shift));
        betdn += dnterm;
        ssum += xmult * betdn;
        if !(xmult * betdn > eps * ssum && i > i_min && xmult >= f64::MIN_POSITIVE) {
            break;
        }
    }

    i = icent + 1.0;
    xmult = centwt;

    // x^(aup - 1) y^b / ((aup - 1) B(aup - 1, b)), with `aup - 1` not rounded to 0 if dfn < eps
    let a = dfn / 2.0 + (icent - 1.0);
    // Unlike in cdflib, from the first term of the down sum if that avoids the cancellation: the
    // result depends on the difference of both sums, so their first terms should share their
    // rounding error, which is independent for `a` and `adn` otherwise
    let mut upterm = if asymp {
        ldexp_f64(cent_term, scale_up - scale) * aup / ((a + b) * xx)
    } else {
        (ln_beta_term(a, b, xx, yy, lnxx, lnyy).0 - f64::from(scale) * LN_2).exp()
    };
    betup = ldexp_f64(betup, scale_up - scale);

    loop {
        xmult *= xnonc / i;
        i += 1.0;
        aup += 1.0;
        upterm *= (aup + b - 2.0) * xx / (aup - 1.0);
        betup -= upterm;
        ssum += xmult * betup;
        if !(xmult * betup > eps * ssum) {
            break;
        }
    }

    let first = if i_min > 0.0 {
        fdtr(dfn, dfd, f) * (-xnonc).exp()
    } else {
        0.0
    };
    ldexp_f64(ssum * (dfd / (2.0 * b)), scale) + first
}

/// Inverse of the noncentral F CDF with respect to `dfd`
///
/// Corresponds to [`scipy.special.ncfdtridfd`][ncfdtridfd].
///
/// [ncfdtridfd]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtridfd.html
///
/// # See also
/// - [`ncfdtridfn`]: Inverse of the noncentral F CDF with respect to `dfn`
/// - [`ncfdtrinc`]: Inverse of the noncentral F CDF with respect to `nc`
/// - [`fdtridfd`]: Inverse of [`fdtr`](crate::fdtr) with respect to `b`
#[must_use]
pub fn ncfdtridfd(dfn: f64, p: f64, nc: f64, f: f64) -> f64 {
    if !(dfn > 0.0 && (0.0..=1.0).contains(&p) && (0.0..=NC_MAX).contains(&nc) && is_pos_finite(f))
    {
        return f64::NAN;
    }
    limit_if_p_is_0_or_1(p, true)
        .unwrap_or_else(|| find_df_root(|dfd| cumfnc(f, dfn, dfd, nc) - p, true))
}

/// Inverse of the noncentral F CDF with respect to `dfn`
///
/// Corresponds to [`scipy.special.ncfdtridfn`][ncfdtridfn].
///
/// [ncfdtridfn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtridfn.html
///
/// # See also
/// - [`ncfdtridfd`]: Inverse of the noncentral F CDF with respect to `dfd`
/// - [`ncfdtrinc`]: Inverse of the noncentral F CDF with respect to `nc`
#[must_use]
pub fn ncfdtridfn(p: f64, dfd: f64, nc: f64, f: f64) -> f64 {
    if !(dfd > 0.0 && (0.0..=1.0).contains(&p) && (0.0..=NC_MAX).contains(&nc) && is_pos_finite(f))
    {
        return f64::NAN;
    }
    // ncfdtr goes from exp(-nc/2) at dfn = 0 to chdtrc(dfd, dfd / f) at dfn = inf
    let increasing = crate::chdtrc(dfd, dfd / f) > (-nc / 2.0).exp();
    limit_if_p_is_0_or_1(p, increasing)
        .unwrap_or_else(|| find_df_root(|dfn| cumfnc(f, dfn, dfd, nc) - p, increasing))
}

/// Inverse of the noncentral F CDF with respect to `nc`
///
/// Corresponds to [`scipy.special.ncfdtrinc`][ncfdtrinc].
///
/// [ncfdtrinc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtrinc.html
///
/// # See also
/// - [`ncfdtridfd`]: Inverse of the noncentral F CDF with respect to `dfd`
/// - [`ncfdtridfn`]: Inverse of the noncentral F CDF with respect to `dfn`
#[must_use]
pub fn ncfdtrinc(dfn: f64, dfd: f64, p: f64, f: f64) -> f64 {
    if !(dfn > 0.0 && dfd > 0.0 && (0.0..=1.0).contains(&p) && is_pos_finite(f)) {
        return f64::NAN;
    }
    // the CDF decreases in nc
    limit_if_p_is_0_or_1(p, false).unwrap_or_else(|| {
        find_positive_root(
            |nc| cumfnc(f, dfn, dfd, nc) - p,
            f64::MIN_POSITIVE,
            NC_MAX,
            false,
        )
    })
}

/// CDF of the noncentral t distribution
///
/// Corresponds to [`scipy.special.nctdtr`][nctdtr]. It is slower, but also accurate where SciPy's
/// loses all digits or returns NaN, e.g. for `t < 0 < nc` or a large `nc`.
///
/// [nctdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtr.html
///
/// # See also
/// - [`stdtr`](crate::stdtr): CDF of the (central) Student's t distribution
/// - [`nctdtridf`]: Inverse of [`nctdtr`] with respect to `df`
/// - [`nctdtrinc`]: Inverse of [`nctdtr`] with respect to `nc`
#[must_use]
pub fn nctdtr(df: f64, nc: f64, t: f64) -> f64 {
    if !(df > 0.0 && !nc.is_nan() && !t.is_nan()) {
        f64::NAN
    } else if t.is_infinite() {
        // the limits of SciPy, where infinite `t` goes before infinite `nc`
        if t > 0.0 { 1.0 } else { 0.0 }
    } else if nc.is_infinite() {
        f64::NAN
    } else if df.is_infinite() {
        crate::ndtr(t - nc)
    } else if t == 0.0 || 0.5 * df == 0.0 {
        // the `df -> 0` limit, also if `df / 2` underflows
        crate::ndtr(-nc)
    } else {
        // Unlike in SciPy, quadrature, as Boost's series cancels for `t < 0 < nc`, e.g.
        // `nctdtr(10, 10, -1) = 2.6e-27`
        nct_quad(t, df, nc)
    }
}

/// Inverse of the noncentral t CDF with respect to `t`
///
/// Corresponds to [`scipy.special.nctdtrit`][nctdtrit]. It is slower, but also accurate where
/// SciPy's loses all digits or returns NaN, e.g. for a small `df` or a large `nc`.
///
/// [nctdtrit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtrit.html
///
/// # See also
/// - [`nctdtr`]: CDF of the noncentral t distribution
/// - [`stdtri`](crate::stdtri): Inverse of the (central) Student's t CDF
#[must_use]
pub fn nctdtrit(df: f64, nc: f64, p: f64) -> f64 {
    if !(df > 0.0 && nc.is_finite() && (0.0..=1.0).contains(&p)) {
        f64::NAN
    } else if df.is_infinite() || p == 0.0 || p == 1.0 {
        // Unlike in SciPy, -inf if `p = 0`, instead of inf
        nc + crate::ndtri(p)
    } else if p > 0.5 {
        // `F(t; nc) = 1 - F(-t; -nc)`, since `1 - p` is exact, but not `nctdtr` near 1
        -nctdtrit(df, -nc, 1.0 - p)
    } else {
        // Unlike in SciPy, bracketed, with Chandrupatla in `asinh(t)`, as the tails can be as
        // heavy as `|t|^-df`
        let func = |x: f64| nctdtr(df, nc, x.sinh()) - p;
        let x_max = f64::MAX.asinh();
        let x0 = (nc + crate::ndtri(p)).asinh();
        let (xl, xr, f_xl, f_xr, bracket_status) =
            bracket_root_for_cdf_inversion(func, x0, -x_max, x_max, -1.0, 1.0, 2.0, 2.0, true, 64);
        // relative and absolute, as the rounding of `nctdtr` can jump at `t = 0`
        let tol = f64::EPSILON;
        match bracket_status {
            0 => match find_root_chandrupatla(func, xl, xr, f_xl, f_xr, tol, tol, 100) {
                (x, 0) => x.sinh(),
                _ => f64::NAN,
            },
            1 => f64::NEG_INFINITY,
            2 => f64::INFINITY,
            _ => f64::NAN,
        }
    }
}

/// CDF of the noncentral t distribution, from cdflib's `cumtnc`.
///
/// It sums the series for the upper tail on p. 532 of Johnson, Kotz and Balakrishnan, Continuous
/// Univariate Distributions, Vol. 2, 2nd ed. (1995), forward and backward from the term with the
/// largest Poisson weight.
///
/// Unlike in SciPy, and as in [`cumfnc`],
/// - the series is summed to machine precision, instead of to a relative tolerance of `1e-7`,
/// - `t` and `nc` are treated as 0 below `eps` and `1e-102` instead of `1e-10`, and the shortcuts
///   for `B(2 cent)` and `B(2 cent + 1)` near 0 or 1 are dropped, since they cut off the sum at
///   `1e-10`, or lose the terms below `cent`,
/// - the Poisson weights are computed with `igam_fac` (directly if `cent = 1`), and `ln(1 - x)`
///   with `ln_1p`,
/// - `s(cent)` is computed with `lbeta`, since `lgam` loses `~df ln(df) eps`,
/// - the beta terms are rescaled by powers of 2 if the central ones underflow,
/// - `cumt`, `cumnor` and `bratio` are replaced by [`stdtr`](crate::stdtr),
///   [`ndtr`](crate::ndtr) and `betainc_xy`.
///
/// As in SciPy, the lower tail for `t > 0` is computed as 1 minus the upper tail, so that it is
/// lost below `~eps`, and the series cancels if `nc` and `t` have opposite signs.
fn cumtnc(t: f64, df: f64, pnonc: f64) -> f64 {
    let tiny = f64::EPSILON;

    // the weights below, `~lambda^1.5` if `lambda = nc^2 / 2 < 1`, are subnormal if `|nc| < 1e-102`
    if pnonc.abs() < 1e-102 {
        return crate::stdtr(df, t);
    }
    let qrevs = t < 0.0;
    let tt = t.abs();
    let dpnonc = if qrevs { -pnonc } else { pnonc };
    let pnonc2 = dpnonc * dpnonc;
    let t2 = tt * tt;
    // t is effectively 0, also if 1 - x underflows
    if tt <= tiny || t2 / df == 0.0 {
        return crate::ndtr(-pnonc);
    }

    let lmbda = 0.5 * pnonc2;
    let x = df / (df + t2);
    // not t2 / (df + t2), which is NaN if t2 overflows
    let omx = 1.0 / (1.0 + df / t2);
    let (lnx, lnomx) = if x > 0.5 {
        ((-omx).ln_1p(), omx.ln())
    } else {
        (x.ln(), (-x).ln_1p())
    };
    // Case : i = lambda
    let cent = lmbda.floor().max(1.0);
    // Compute d=T(2i) and e=t(2i+1), offset by exp(-lambda)
    // lambda^cent exp(-lambda) / cent! and lambda^(cent + 1/2) exp(-lambda) / Gamma(cent + 3/2),
    // where `igam_fac` loses `~|ln(lambda)| eps` if lambda < 1
    let (dcent, mut ecent) = if cent == 1.0 {
        let dcent = lmbda * (-lmbda).exp();
        (dcent, dcent * lmbda.sqrt() / (0.75 * SQRT_PI))
    } else {
        (
            igam_fac(cent + 1.0, lmbda) / lmbda,
            igam_fac(cent + 1.5, lmbda) / lmbda,
        )
    };
    if dpnonc < 0.0 {
        ecent = -ecent;
    }
    // Compute bcent=B(2*cent) and bbcent=B(2*cent+1). If they underflow, the terms are carried as
    // multiples of `2^scale`, since the later ones can still be significant.
    let (bcent, scale_b) = betainc_scaled(0.5 * df, cent + 0.5, x, omx, lnx, lnomx);
    let (bbcent, scale_bb) = betainc_scaled(0.5 * df, cent + 1.0, x, omx, lnx, lnomx);
    let scale_cent = scale_b.max(scale_bb);
    let bcent = ldexp_f64(bcent, scale_b - scale_cent);
    let bbcent = ldexp_f64(bbcent, scale_bb - scale_cent);
    let mut scale = scale_cent;

    // First term in ccum is D*B + E*BB
    let mut ccum = dcent * bcent + ecent * bbcent;
    // Compute s(cent) = B(2*(cent+1)) - B(2*cent))
    let scent =
        (-lbeta(0.5 * df, cent + 0.5) - (cent + 0.5).ln() + 0.5 * df * lnx + (cent + 0.5) * lnomx
            - f64::from(scale) * LN_2)
            .exp();
    // Compute ss(cent) = B(2*cent+3) - B(2*cent+1)
    let sscent =
        (-lbeta(0.5 * df, cent + 1.0) - (cent + 1.0).ln() + 0.5 * df * lnx + (cent + 1.0) * lnomx
            - f64::from(scale) * LN_2)
            .exp();
    // Sum forward. If nc < 0, the sum can cancel to 0 or below, and terms that are NaN end it too.
    let mut xi = cent + 1.0;
    let mut twoi = 2.0 * xi;
    let mut d = dcent;
    let mut e = ecent;
    let mut b = bcent;
    let mut s = scent;
    let mut bb = bbcent;
    let mut ss = sscent;
    loop {
        rescale(
            (b + s).max(bb + ss),
            1.0,
            1.0,
            &mut scale,
            &mut [&mut b, &mut bb, &mut s, &mut ss, &mut ccum],
        );
        b += s;
        bb += ss;
        d *= lmbda / xi;
        e *= lmbda / (xi + 0.5);
        let term = d * b + e * bb;
        ccum += term;
        s *= omx * (df + twoi - 1.0) / (twoi + 1.0);
        ss *= omx * (df + twoi) / (twoi + 2.0);
        xi += 1.0;
        twoi = 2.0 * xi;
        // if scaled, not past the underflow of the Poisson weights, which would get stuck at a
        // subnormal, while `b` and `bb` grow
        if !(term.abs() > f64::EPSILON * ccum.abs())
            || scale < 0 && d < f64::MIN_POSITIVE && e.abs() < f64::MIN_POSITIVE
        {
            break;
        }
    }
    // Sum Backward
    xi = cent;
    twoi = 2.0 * xi;
    d = dcent;
    e = ecent;
    let shift = scale_cent - scale;
    b = ldexp_f64(bcent, shift);
    bb = ldexp_f64(bbcent, shift);
    s = ldexp_f64(scent, shift) * (1.0 + twoi) / ((df + twoi - 1.0) * omx);
    ss = ldexp_f64(sscent, shift) * (2.0 + twoi) / ((df + twoi) * omx);
    loop {
        b -= s;
        bb -= ss;
        d *= xi / lmbda;
        e *= (xi + 0.5) / lmbda;
        let term = d * b + e * bb;
        ccum += term;
        xi -= 1.0;
        if xi < 0.5 {
            break;
        }
        twoi = 2.0 * xi;
        s *= (1.0 + twoi) / ((df + twoi - 1.0) * omx);
        ss *= (2.0 + twoi) / ((df + twoi) * omx);
        if !(term.abs() > f64::EPSILON * ccum.abs()) {
            break;
        }
    }
    // Due to roundoff error the answer may not lie between zero and one
    // Force it to do so
    let ccum = ldexp_f64(0.5 * ccum, scale).clamp(0.0, 1.0);
    if qrevs { ccum } else { 1.0 - ccum }
}

/// CDF of the noncentral t distribution, as `E[ndtr(t W - nc)]` over the chi distribution of
/// `W = sqrt(V / df)` with `V ~ chi^2(df)`, by quadrature in `s = ln(W)`.
///
/// The integrand is unimodal in `s`. Below `w0` it is split into `ndtr(-nc) P(W < w0)` and the
/// difference, whose tail is far shorter if `df` is small.
fn nct_quad(t: f64, df: f64, nc: f64) -> f64 {
    // above the `50 eps` floor of the error estimate of `quad`
    const EPSREL: f64 = 1e-13;
    let a = 0.5 * df;
    // the density `2 x^a e^-x / Gamma(a)` of `s`, with `x = a e^(2 s)` left unrounded
    let peak = 2.0 * igam_fac(a, a);
    let density = |s: f64| {
        let u = 2.0 * s;
        // `ln(x / a) - (x / a - 1)`
        let ln_ratio = if u.abs() < 1.0 {
            log1pmx(u.exp_m1())
        } else {
            u - u.exp_m1()
        };
        peak * (a * ln_ratio).exp()
    };

    let phi0 = crate::ndtr(-nc);
    let w0 = (0.5 / (t.abs() * (nc.abs() + 2.0))).min(0.5);
    let s0 = w0.ln();
    // `P(W < w0)`, by its leading term `y^a / Gamma(a + 1)` where `y = a w0^2` underflows
    let y = a * w0 * w0;
    let head = phi0
        * if y < f64::MIN_POSITIVE {
            (a * (a.ln() + 2.0 * s0)).exp() * crate::rgamma(a + 1.0)
        } else {
            igam(a, y)
        };
    // `t e^s - nc`, from a finite `t - nc` if `e^s > 1/2`, as `e^s` rounds away a small `s`, with
    // the rounding error of `t - nc`, which would otherwise bias all nodes alike
    let (t_nc, t_nc_err) = two_sum(t, -nc);
    let z = |s: f64| {
        if s > -LN_2 && t_nc.is_finite() {
            t_nc + (t * s.exp_m1() + t_nc_err)
        } else {
            t * s.exp() - nc
        }
    };
    // the difference has the sign of `t`, and its tail decays as `w^(df + 1)`
    let diff = |s: f64| density(s) * (crate::ndtr(z(s)) - phi0).abs();
    let below = quad_unimodal(diff, f64::NEG_INFINITY, s0, s0, 1.0 / (df + 1.0), EPSREL);

    // the derivatives of `ln(density(s) ndtr(z))`, with `z = t w - nc` and `w = e^s`
    let slopes = |s: f64| {
        let w = s.exp();
        let (x, z) = (t * w, z(s));
        // ndtr'(z) / ndtr(z)
        let h = SQRT_FRAC_2_PI / erfcx(-z / SQRT_2);
        let xh = x * h;
        let d1 = xh - df * (2.0 * s).exp_m1();
        // clamped, as `x z` can overflow where `x h` is 0
        let d2 = -2.0 * df * w * w + xh * (1.0 - x * (z + h)).max(-f64::MAX);
        (d1, d2)
    };
    let d1 = |s: f64| slopes(s).0;
    // the mode, to well within the width `~1 / sqrt(2 df)` of the density
    let mode = if d1(s0) > 0.0 {
        // bracketed in doubling steps from 0, since the mode is `~x h / (2 df)` if `df` is large
        let width = 1.0 / (df + 1.0).sqrt();
        let up = d1(0.0) > 0.0;
        let (mut near, mut step) = (0.0, width);
        let far = loop {
            let far = if up { near + step } else { s0.max(near - step) };
            if (d1(far) > 0.0) != up {
                break far;
            }
            (near, step) = (far, 2.0 * step);
        };
        let (lo, hi) = if up { (near, far) } else { (far, near) };
        find_root_chandrupatla(d1, lo, hi, d1(lo), d1(hi), 1e-8, 1e-3 * width, 100).0
    } else {
        s0
    };
    // at most 1, the scale of the density if `df` is small, and at most `1 / |x|`, the scale of the
    // step of `ndtr(x - nc)` that the curvature misses beyond it
    let x = t * mode.exp();
    let width = 1.0 / (-slopes(mode).1).max(1.0).max(x * x).min(f64::MAX).sqrt();
    let integrand = |s: f64| density(s) * crate::ndtr(z(s));
    let above = quad_unimodal(integrand, s0, f64::INFINITY, mode, width, EPSREL);

    // the rounding of the sum can exceed 1, but keep NaN
    let p = head + t.signum() * below + above;
    if p > 1.0 { 1.0 } else { p }
}

/// The largest `|nc|` for which [`cumtnc`] is computed, as in cdflib, since it sums `~12 |nc|`
/// terms.
const NCT_MAX: f64 = 1e6;

/// Inverse of the noncentral t CDF with respect to `df`
///
/// Corresponds to [`scipy.special.nctdtridf`][nctdtridf].
///
/// [nctdtridf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtridf.html
///
/// # See also
/// - [`nctdtrinc`]: Inverse of the noncentral t CDF with respect to `nc`
/// - [`stdtridf`]: Inverse of [`stdtr`](crate::stdtr) with respect to `nu`
#[must_use]
pub fn nctdtridf(p: f64, nc: f64, t: f64) -> f64 {
    // the CDF is ndtr(-nc) if t = 0, and 1 or 0 if t = +-inf, independent of df
    if !((0.0..=1.0).contains(&p) && nc.abs() <= NCT_MAX && t != 0.0 && t.is_finite()) {
        return f64::NAN;
    }
    let increasing = t > 0.0;
    if let Some(df) = limit_if_p_is_0_or_1(p, increasing) {
        return df;
    }
    let func = |df: f64| cumtnc(t, df, nc) - p;
    // The CDF goes from ndtr(-nc) at df = 0 to ndtr(t - nc) at df = inf. If p lies outside, and
    // also beyond the CDF at the search points, only an extremum in between can reach it, which
    // saves the searches of `find_df_root` up to its bounds, of ~110 CDFs of ~12 |nc| terms each.
    let (at_0, at_inf) = (crate::ndtr(-nc) - p, crate::ndtr(t - nc) - p);
    if (at_0 > 0.0 && at_inf > 0.0) || (at_0 < 0.0 && at_inf < 0.0) {
        let closest = if at_0.abs() < at_inf.abs() {
            0.0
        } else {
            f64::INFINITY
        };
        // With `W = sqrt(V / df)` and `V ~ chi^2(df)`, the CDF is `E[ndtr(t W - nc)]`, so that its
        // distance from 0 if `t > 0`, or from 1 if `t < 0`, is `E[ndtr(|t| W - s)]` with
        // `s = sign(t) nc`. This is at least `ndtr(-s)`, and by Markov's inequality
        // `P(W >= w) <= 1 / w^2` with `w = s / (2 |t|)`, at most `4 t^2 / s^2 + ndtr(-s / 2)` if
        // `s > 0`, and the same distance `q` of `p` must lie in between.
        let (q, s) = if t > 0.0 { (p, nc) } else { (1.0 - p, -nc) };
        // `q < ndtr(-s)`, without the rounding of `1 - p`
        let below_min = (t > 0.0) == (at_0 > 0.0);
        if below_min || s > 0.0 && q > 4.0 * t * t / (s * s) + crate::ndtr(-s / 2.0) {
            return closest;
        }
        if let Some(points) = scan_same_sign(func) {
            return find_extremum_root(func, &points, increasing).unwrap_or(closest);
        }
    }
    find_df_root(func, increasing)
}

/// Inverse of the noncentral t CDF with respect to `nc`
///
/// Corresponds to [`scipy.special.nctdtrinc`][nctdtrinc].
///
/// [nctdtrinc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtrinc.html
///
/// # See also
/// - [`nctdtridf`]: Inverse of the noncentral t CDF with respect to `df`
/// - [`stdtridf`]: Inverse of [`stdtr`](crate::stdtr) with respect to `nu`
#[must_use]
pub fn nctdtrinc(df: f64, p: f64, t: f64) -> f64 {
    if !(df > 0.0 && (0.0..=1.0).contains(&p) && t.is_finite()) {
        return f64::NAN;
    }
    // The CDF decreases in nc, from 1 at -inf to 0 at inf.
    if p == 0.0 || p == 1.0 {
        return if p == 0.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
    }
    let func = |nc: f64| cumtnc(t, df, nc) - p;
    // The bracket is expanded in steps of 1, 2, 4, ... from 0, and 64 steps are more than enough
    // to reach NCT_MAX.
    let (xl, xr, f_xl, f_xr, bracket_status) = bracket_root_for_cdf_inversion(
        func, 0.0, -NCT_MAX, NCT_MAX, -1.0, 1.0, 2.0, 2.0, false, 64,
    );
    match bracket_status {
        0 => match find_root_chandrupatla(func, xl, xr, f_xl, f_xr, f64::EPSILON, 1e-100, 100) {
            (nc, 0) => nc,
            _ => f64::NAN,
        },
        1 => f64::NEG_INFINITY,
        2 => f64::INFINITY,
        _ => f64::NAN,
    }
}
