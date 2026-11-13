//! Inverses of distribution functions with respect to a parameter, which `scipy.special` computes
//! with cdflib. Unlike in SciPy, the roots are found with the bracket search and Chandrupatla
//! solver of [`gdtrib`](crate::gdtrib) to machine precision, instead of cdflib's `dinvr` with a
//! relative tolerance of `1e-10`. If `p` can't be reached, 0 or inf is returned, whichever limit
//! is closest, instead of cdflib's search bound. If `p` is 0 or 1, the limit where the CDF is
//! smallest or largest is returned, instead of a root where the CDF underflows or rounds to `p`.

use crate::xsf::{find_positive_root, find_root_chandrupatla};

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
    let q = 1.0 - p;
    // stdtr(df, -t) = 1 - stdtr(df, t) is used for the upper tail
    let func = |df: f64| {
        if p <= q {
            crate::stdtr(df, t) - p
        } else {
            q - crate::stdtr(df, -t)
        }
    };
    // stdtr(df, t) increases from 1/2 at df = 0 to ndtr(t) at df = inf if t > 0
    limit_if_p_is_0_or_1(p, t > 0.0).unwrap_or_else(|| find_df_root(func, t > 0.0))
}

/// Inverse of [`fdtr`](crate::fdtr) with respect to `b`
///
/// Corresponds to [`scipy.special.fdtridfd`][fdtridfd].
///
/// For `f < 1`, `fdtr` can peak in `b`, so that there can be two roots, of which one is returned.
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
    let q = 1.0 - p;
    let func = |dfd: f64| {
        if p <= q {
            crate::fdtr(dfn, dfd, f) - p
        } else {
            q - crate::fdtrc(dfn, dfd, f)
        }
    };
    // fdtr(dfn, dfd, f) goes from 0 at dfd = 0 to chdtr(dfn, dfn f) at dfd = inf
    limit_if_p_is_0_or_1(p, true).unwrap_or_else(|| find_df_root(func, true))
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
    let (x, status) =
        find_root_chandrupatla(&func, x1, x2, func(x1), func(x2), f64::EPSILON, 1e-100, 100);
    (status == 0).then_some(x)
}
