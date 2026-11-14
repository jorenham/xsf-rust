//! Inverses of distribution functions with respect to a parameter, which `scipy.special` computes
//! with cdflib or Boost. The roots are found with the bracket search and Chandrupatla solver of
//! [`gdtrib`](crate::gdtrib) to machine precision, whereas cdflib's `dinvr` stops at a relative
//! tolerance of `1e-10` (`1e-8` for the noncentral t). If the probability can't be reached, 0 or
//! inf is returned, whichever limit is closest, instead of cdflib's search bound. If it is 0 or 1,
//! the limit where the CDF is smallest or largest is returned, instead of a root where the CDF
//! underflows or rounds to it.

use core::f64::consts::{LN_2, PI};

use crate::utils::pow2;
use crate::xsf::cephes::dd_real::ldexp_f64;
use crate::xsf::cephes::{betainc_xy, igam_fac, incbet_cf_factor, lbeta};
use crate::xsf::{bracket_root_for_cdf_inversion, find_positive_root, find_root_chandrupatla};
use crate::{betainc, betaincc};

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
    if !(0.0..=1.0).contains(&p) || t == 0.0 || !t.is_finite() {
        return f64::NAN;
    }
    // stdtr(df, -t) = 1 - stdtr(df, t) is used for the upper tail
    let func = cdf_minus_p(p, |df| crate::stdtr(df, t), |df| crate::stdtr(df, -t));
    // stdtr(df, t) increases from 1/2 at df = 0 to ndtr(t) at df = inf if t > 0
    limit_if_p_is_0_or_1(p, t > 0.0).unwrap_or_else(|| find_df_root(func, t > 0.0))
}

/// Inverse of [`fdtr`](crate::fdtr) with respect to `b`
///
/// Corresponds to [`scipy.special.fdtridfd`][fdtridfd].
///
/// For `f < 1`, `fdtr` can peak in `b`, so that there can be two roots, of which one is returned.
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
    if !(dfn > 0.0) || !(0.0..=1.0).contains(&p) || !(f > 0.0) || f.is_infinite() {
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
/// The CDF `I_{1-p}(n - k, k + 1)` is continuous in `k`, so that `k` is real. If `y` lies below the
/// CDF at `k = 0`, i.e. `(1 - p)^n`, 0 is returned.
///
/// [bdtrik]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrik.html
///
/// # See also
/// - [`bdtr`](crate::bdtr): Binomial distribution function
/// - [`bdtri`](crate::bdtri): Inverse of [`bdtr`](crate::bdtr) with respect to `p`
/// - [`bdtrin`]: Inverse of [`bdtr`](crate::bdtr) with respect to `n`
#[must_use]
pub fn bdtrik(y: f64, n: f64, p: f64) -> f64 {
    if !(0.0..=1.0).contains(&y) || !(n >= 0.0) || n.is_infinite() || !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    if y == 0.0 {
        return 0.0;
    }
    // the CDF is 0 below k = n if p = 1
    if y == 1.0 || p == 1.0 {
        return n;
    }
    // the CDF increases from (1 - p)^n at k = 0 to 1 at k = n, beyond which it stays 1
    let func = cdf_minus_p(
        y,
        |k| {
            if k < n {
                betaincc(k + 1.0, n - k, p)
            } else {
                1.0
            }
        },
        |k| {
            if k < n {
                betainc(k + 1.0, n - k, p)
            } else {
                0.0
            }
        },
    );
    find_positive_root(func, f64::MIN_POSITIVE, n, true)
}

/// Inverse of [`bdtr`](crate::bdtr) with respect to `n`
///
/// Corresponds to [`scipy.special.bdtrin`][bdtrin].
///
/// The CDF `I_{1-p}(n - k, k + 1)` is continuous in `n`, so that `n` is real.
///
/// [bdtrin]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrin.html
///
/// # See also
/// - [`bdtr`](crate::bdtr): Binomial distribution function
/// - [`bdtri`](crate::bdtri): Inverse of [`bdtr`](crate::bdtr) with respect to `p`
/// - [`bdtrik`]: Inverse of [`bdtr`](crate::bdtr) with respect to `k`
#[must_use]
pub fn bdtrin(k: f64, y: f64, p: f64) -> f64 {
    if !(k.is_finite() && k >= 0.0 && (0.0..=1.0).contains(&y) && (0.0..=1.0).contains(&p)) {
        return f64::NAN;
    }
    // the limits n = inf and n = k, as in SciPy also if p is 0 or 1
    if y == 0.0 {
        return f64::INFINITY;
    }
    if y == 1.0 {
        return k;
    }
    // as in SciPy, NaN if p is 0 or 1, where the CDF is 1 or 0 for all n > k
    if p == 0.0 || p == 1.0 {
        return f64::NAN;
    }
    // The CDF decreases from 1 at n = k to 0 at n = inf. It is solved for `n - k`, which keeps
    // its precision if `n` is close to `k`.
    let func = cdf_minus_p(
        y,
        |nk| betaincc(k + 1.0, nk, p),
        |nk| betainc(k + 1.0, nk, p),
    );
    k + find_positive_root(func, f64::MIN_POSITIVE, f64::MAX, false)
}

/// Inverse of [`betainc`] with respect to `a`
///
/// Corresponds to [`scipy.special.btdtria`][btdtria].
///
/// Unlike in SciPy, NaN is returned if `x > 1`, also if `p` is 0 or 1.
///
/// [btdtria]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.btdtria.html
///
/// # See also
/// - [`betainc`]: Regularized incomplete beta function
/// - [`betaincinv`](crate::betaincinv): Inverse of [`betainc`] with respect to `x`
/// - [`btdtrib`]: Inverse of [`betainc`] with respect to `b`
#[must_use]
pub fn btdtria(p: f64, b: f64, x: f64) -> f64 {
    if !((0.0..=1.0).contains(&p) && b > 0.0 && x > 0.0 && x <= 1.0) {
        return f64::NAN;
    }
    // The CDF decreases from 1 at a = 0 to 0 at a = inf. As in SciPy, the limits are also returned
    // if b = inf or x = 1, where the CDF is 1 for all a.
    limit_if_p_is_0_or_1(p, false).unwrap_or_else(|| {
        if b.is_infinite() || x == 1.0 {
            return f64::NAN;
        }
        let func = cdf_minus_p(p, |a| betainc(a, b, x), |a| betaincc(a, b, x));
        find_positive_root(func, f64::MIN_POSITIVE, f64::MAX, false)
    })
}

/// Inverse of [`betainc`] with respect to `b`
///
/// Corresponds to [`scipy.special.btdtrib`][btdtrib].
///
/// Unlike in SciPy, NaN is returned if `x > 1`, also if `p` is 0 or 1.
///
/// [btdtrib]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.btdtrib.html
///
/// # See also
/// - [`betainc`]: Regularized incomplete beta function
/// - [`betaincinv`](crate::betaincinv): Inverse of [`betainc`] with respect to `x`
/// - [`btdtria`]: Inverse of [`betainc`] with respect to `a`
#[must_use]
pub fn btdtrib(a: f64, p: f64, x: f64) -> f64 {
    if !(a > 0.0 && (0.0..=1.0).contains(&p) && x > 0.0 && x <= 1.0) {
        return f64::NAN;
    }
    // The CDF increases from 0 at b = 0 to 1 at b = inf. As in SciPy, the limits are also returned
    // if a = inf or x = 1, where the CDF is 0 or 1 for all b.
    limit_if_p_is_0_or_1(p, true).unwrap_or_else(|| {
        if a.is_infinite() || x == 1.0 {
            return f64::NAN;
        }
        let func = cdf_minus_p(p, |b| betainc(a, b, x), |b| betaincc(a, b, x));
        find_positive_root(func, f64::MIN_POSITIVE, f64::MAX, true)
    })
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

/// The largest noncentrality for which [`cumfnc`] is computed: like in cdflib, where `nc / 2` must
/// fit an `int`, since it sums `~17 sqrt(nc)` terms.
const NC_MAX: f64 = 2.0 * i32::MAX as f64;

/// The exponent of the largest power of 2 by which [`rescale`] divides at once.
const RESCALE: i32 = 500;

/// `ln(x^a y^b / (a B(a, b)))`, with `y = 1 - x`, i.e. `ln(I_x(a, b) - I_x(a + 1, b))`
/// (A&S 26.5.16)
fn ln_beta_term(a: f64, b: f64, lnx: f64, lny: f64) -> f64 {
    -lbeta(a, b) - a.ln() + a * lnx + b * lny
}

/// `(m, k)` with `I_x(a, b) = m 2^k`, where `m ~ 1` and `k < 0` if `I_x(a, b)` underflows, so that
/// the terms of a series that starts there can be carried as multiples of `2^k`
fn betainc_scaled(a: f64, b: f64, x: f64, y: f64, lnx: f64, lny: f64) -> (f64, i32) {
    let ix = betainc_xy(a, b, x, y);
    if ix >= f64::MIN_POSITIVE {
        return (ix, 0);
    }
    let ln_term = ln_beta_term(a, b, lnx, lny);
    if !ln_term.is_finite() {
        return (ix, 0);
    }
    // I_x(a, b) = x^a y^b / (a B(a, b)) w, with x below the mean if I_x underflows
    let (w, _) = incbet_cf_factor(a, b, x, y);
    if !(w > 0.0 && w.is_finite()) {
        return (ix, 0);
    }
    let k = ((ln_term + w.ln()) / LN_2).round() as i32;
    ((ln_term - f64::from(k) * LN_2).exp() * w, k)
}

/// Divides the `values` by powers of 2, adding their exponents to `scale`, while `peak * factor`
/// lies above `2^RESCALE`, so that the scaled terms of [`betainc_scaled`] stay finite, but not past
/// `scale = 0`, so that they stay below their unscaled values.
///
/// `factor` is separate, since `peak * factor` can overflow before the division.
fn rescale(mut peak: f64, factor: f64, scale: &mut i32, values: &mut [&mut f64]) {
    while *scale < 0 && peak * factor > pow2(RESCALE) {
        let k = RESCALE.min(-*scale);
        peak *= pow2(-k);
        for v in values.iter_mut() {
            **v *= pow2(-k);
        }
        *scale += k;
    }
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
        return crate::fdtr(dfn, dfd, f);
    }
    if !(pnonc <= NC_MAX) {
        return f64::NAN;
    }

    let xnonc = pnonc / 2.0;
    let icent = xnonc.floor().max(1.0);
    // xnonc^icent exp(-xnonc) / icent!, where `igam_fac` loses `~|ln(xnonc)| eps` if xnonc < 1
    let centwt = if icent == 1.0 {
        xnonc * (-xnonc).exp()
    } else {
        igam_fac(icent + 1.0, xnonc) / xnonc
    };
    let prod = dfn * f;
    let dsum = dfd + prod;
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
        xx = prod / dsum;
        yy = 1.0 - xx;
        (lnxx, lnyy) = (xx.ln(), (-xx).ln_1p());
    } else {
        xx = 1.0 - yy;
        // without the rounding of a subnormal `yy`
        let ln_y = if prod.is_finite() {
            yy.ln()
        } else {
            dfd.ln() - dfn.ln() - f.ln()
        };
        (lnxx, lnyy) = ((-yy).ln_1p(), ln_y);
    }
    let mut adn = dfn / 2.0 + icent;
    let mut aup = adn;
    let b = dfd / 2.0;
    // If `I_x(adn, b)` underflows, the beta terms are carried as multiples of `2^scale`, since the
    // terms below `icent` can still be significant, e.g. if x is tiny.
    let (mut betdn, mut scale) = betainc_scaled(adn, b, xx, yy, lnxx, lnyy);
    let scale_up = scale;
    let mut betup = betdn;
    let mut ssum = centwt * betdn;

    let mut xmult = centwt;
    let mut i = icent;

    // x^adn y^b / (adn B(adn, b))
    let mut dnterm = (ln_beta_term(adn, b, lnxx, lnyy) - f64::from(scale) * LN_2).exp();

    // `>` instead of `>=`, which never ends if `eps * ssum` underflows, and not past the underflow
    // of the Poisson weight, which would get stuck at a subnormal, while the scaled `betdn` grows
    while xmult * betdn > eps * ssum && i > 0.0 && xmult >= f64::MIN_POSITIVE {
        xmult *= i / xnonc;
        i -= 1.0;
        adn -= 1.0;
        let ratio = (adn + 1.0) / ((adn + b) * xx);
        rescale(
            betdn,
            ratio,
            &mut scale,
            &mut [&mut betdn, &mut dnterm, &mut ssum],
        );
        dnterm *= ratio;
        betdn += dnterm;
        ssum += xmult * betdn;
    }

    i = icent + 1.0;
    xmult = centwt;

    // x^(aup - 1) y^b / ((aup - 1) B(aup - 1, b)), with `aup - 1` not rounded to 0 if dfn < eps
    let a = dfn / 2.0 + (icent - 1.0);
    let mut upterm = (ln_beta_term(a, b, lnxx, lnyy) - f64::from(scale) * LN_2).exp();
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

    ldexp_f64(ssum, scale)
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
    if !(dfn > 0.0)
        || !(0.0..=1.0).contains(&p)
        || !(0.0..=NC_MAX).contains(&nc)
        || !(f > 0.0)
        || f.is_infinite()
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
    if !(dfd > 0.0)
        || !(0.0..=1.0).contains(&p)
        || !(0.0..=NC_MAX).contains(&nc)
        || !(f > 0.0)
        || f.is_infinite()
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
    if !(dfn > 0.0) || !(dfd > 0.0) || !(0.0..=1.0).contains(&p) || !(f > 0.0) || f.is_infinite() {
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
        (dcent, dcent * lmbda.sqrt() / (0.75 * PI.sqrt()))
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
    if !(0.0..=1.0).contains(&p) || !(nc.abs() <= NCT_MAX) || t == 0.0 || !t.is_finite() {
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
    if !(df > 0.0) || !(0.0..=1.0).contains(&p) || !t.is_finite() {
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
