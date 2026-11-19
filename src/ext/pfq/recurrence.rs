//! 3F2 for `|z| < 1` by a recurrence in a parameter
//!
//! The series can cancel at a lower parameter `b` and not at `b + m` for a large `m`, e.g. for
//! large `a_i` and `z < 0`. With `u = b - 1`, `d_i = a_i - u` and the other lower parameter `b'`,
//!
//! ```text
//! c_0 F(b) + c_1 F(b - 1) + c_2 F(b - 2) + c_3 F(b - 3) = 0,
//! c_0 = z d_1 d_2 d_3,
//! c_1 = u (z (d_1 d_2 + d_1 d_3 + d_2 d_3 + d_1 + d_2 + d_3 + 1) - (u - 1) (u - b')),
//! c_2 = u (u - 1) (z (d_1 + d_2 + d_3 + 3) + 2 (u - 1) - b'),
//! c_3 = u (u - 1) (u - 2) (z - 1),
//! ```
//!
//! which is the differential equation (DLMF 16.8.3) with `(theta + b - 1) F(b) = (b - 1) F(b - 1)`
//! for `theta = z d/dz` (DLMF 16.3.4 with `n = 1`).
//!
//! It can also cancel at an upper parameter `a < 0` and not at `a + m`, e.g. where its terms
//! alternate up to `n = -a`, as for `z > 0` and positive other parameters. With `w = a + 1`,
//! `e_j = b_j - w`, and `h_i = w - a_i` for the other two upper parameters,
//!
//! ```text
//! c_0 F(a + 3) + c_1 F(a + 2) + c_2 F(a + 1) + c_3 F(a) = 0,
//! c_0 = w (w + 1) (1 - z),
//! c_1 = w (e_1 + e_2 - w - 2 + z (h_1 + h_2 + 1)),
//! c_2 = w (1 - e_1 - e_2) + e_1 e_2 - z h_1 h_2,
//! c_3 = -e_1 e_2,
//! ```
//!
//! by `(theta + a) F(a) = a F(a + 1)` (DLMF 16.3.3 with `n = 1`).
//!
//! `F` follows in `m` steps from the series at that parameter plus `m`, `m + 1` and `m + 2`, in
//! double-double, in which each of those is exact.
//!
//! The steps can make the errors of the series and their own rounding errors larger by many
//! orders of magnitude, so the error of the result is estimated. An error `e` of the F at `k`
//! steps from the result adds `g_k e` to it, where the `g_k` follow from the adjoint recurrence,
//! upwards from `g_0 = 1`. The estimate is not a bound, as the `g_k` are computed in `f64` from
//! the rounded coefficients.

use core::iter::zip;

use super::{
    Dd, Estimate, PFQ_DD_MAXITER, PFQ_MAXITER, Params, SeriesReal, descending, first_zero,
    pfq_extreme, pfq_ratio, pfq_sum_terms, pfq_terms_of,
};

/// The `m` for which the recurrence starts from the series at `b + m`, in the order in which
/// they are tried; for an upper parameter `a`, it starts from where `a + m` is about minus
/// these, in reverse
///
/// Of 7136 results on 91004 points, 426 are from the last one and 537 from the one before it,
/// and 8192 would give 273 more (measured).
const RECURRENCE_STARTS: [u32; 10] = [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];

/// Bound of the rounding error of a step, relative to the sizes of the terms of its sum:
/// `Dd::EPSILON` for each of its at most 17 inexact operations in a row (`u - 1`, `u - 2`,
/// `z - 1` and the products with `minus` are exact)
const RECURRENCE_EPSILON: f64 = 17.0 * Dd::EPSILON;

/// The coefficients `[c_0, c_1, c_2, c_3]`, for `minus = -1`
///
/// For `minus = 1` and the sizes of the arguments, no terms cancel: each is then the sum of the
/// sizes of its terms, which is the scale of its rounding error.
fn lower_coefficients<T: SeriesReal>(
    u: T,
    [d1, d2, d3]: [T; 3],
    z: T,
    other: T,
    minus: f64,
) -> [T; 4] {
    let (m, sum) = (T::from_f64(minus), d1 + d2 + d3);
    // u - 1 and u - 2
    let (u1, u2) = (u + minus, u + 2.0 * minus);
    [
        d1 * d2 * d3 * z,
        u * ((d1 * d2 + d1 * d3 + d2 * d3 + sum + 1.0) * z + m * u1 * (u + m * other)),
        u * u1 * ((sum + 3.0) * z + u1 + u1 + m * other),
        u * u1 * u2 * (z + minus),
    ]
}

/// As `lower_coefficients`, for the relation in an upper parameter
fn upper_coefficients<T: SeriesReal>(w: T, e: [T; 2], h: [T; 2], z: T, minus: f64) -> [T; 4] {
    let (m, sum) = (T::from_f64(minus), e[0] + e[1]);
    [
        w * (w + 1.0) * (m * z + 1.0),
        w * (sum + m * w + 2.0 * minus + z * (h[0] + h[1] + 1.0)),
        w * (m * sum + 1.0) + e[0] * e[1] + m * z * h[0] * h[1],
        m * e[0] * e[1],
    ]
}

/// The series in double-double for the ratios `ratio(n)` of its terms, which can be those of
/// more exact parameters than `a` and `b`, and an estimate of its absolute error, or `None` if
/// it is not accurate
fn recurrence_start(
    a: [f64; 3],
    b: [f64; 3],
    z: f64,
    ratio: impl Fn(f64) -> Dd,
) -> Option<(Dd, f64)> {
    let terms = pfq_terms_of(first_zero(&a), z, |n, t: Dd| t * (ratio(n) * Dd::new(z)));
    let (a, b) = (descending(a), descending(b));
    // (summed to the accuracy of double-double, as the steps can make its error larger)
    let (y, err, _) = pfq_sum_terms(&a, &b, z, PFQ_DD_MAXITER, terms, Dd::EPSILON);
    let accurate = Estimate::absolute((y.hi, err)).is_accurate();
    accurate.then_some((y, err))
}

/// F in `m` steps from `start(m + k)` for `k = 0, 1, 2`, in double-double, the part of the
/// estimate of its absolute error from the rounding errors of the steps, and that estimate, or
/// `None` if a start is `None`
///
/// `step(k)` gives the `c_i` of the relation `c_0 F_{k+3} + c_1 F_{k+2} + c_2 F_{k+1} + c_3 F_k
/// = 0` between the F at `k` steps from the result, and the sizes of the terms of each `c_i`.
fn recurrence_from(
    m: u32,
    start: impl Fn(f64) -> Option<(Dd, f64)>,
    step: impl Fn(f64) -> ([Dd; 4], [f64; 4]),
) -> Option<(Dd, f64, f64)> {
    let start = |k: u32| start(f64::from(m + k));
    let starts = [start(0)?, start(1)?, start(2)?];
    // `[F_{k+1}, F_{k+2}, F_{k+3}]` for the next step `k`
    let mut f = starts.map(|(y, _)| y);
    // of each step: the factors of those three in its F, and the bound of its rounding error
    let mut steps = Vec::with_capacity(m as usize);
    for k in (0..m).rev().map(f64::from) {
        // (the sizes of the terms of c_0, c_1 and c_2 give those of the sum)
        let ([c0, c1, c2, c3], [e0, e1, e2, _]) = step(k);
        let sum = c0 * f[2] + c1 * f[1] + c2 * f[0];
        let next = -sum / c3;
        let terms = e0 * f[2].hi.abs() + e1 * f[1].hi.abs() + e2 * f[0].hi.abs();
        // an underflow adds at most `f64::MIN_POSITIVE` to the sum, to the products c_0
        // and c_3, and to the quotient; in c_1 and c_2, less than the rounding of their
        // terms without z
        let underflows = 1.0 + f[2].hi.abs() + next.hi.abs() + c3.hi.abs();
        let bound = (RECURRENCE_EPSILON * terms + f64::MIN_POSITIVE * underflows) / c3.hi.abs();
        steps.push(([c2, c1, c0].map(|c| -c.hi / c3.hi), bound));
        f = [next, f[0], f[1]];
    }
    // the `g_k` of the three F above each step, upwards from the result, and the bound of the
    // rounding errors
    let (mut g, mut rounding) = ([1.0_f64, 0.0, 0.0], 0.0);
    for (factors, bound) in steps.iter().rev() {
        rounding += g[0].abs() * bound;
        g = [g[0] * factors[0] + g[1], g[0] * factors[1] + g[2], g[0] * factors[2]];
    }
    let series: f64 = zip(g, starts).map(|(g, (_, err))| g.abs() * err).sum();
    Some((f[0], rounding, rounding + series))
}

/// `recurrence_from` in the lower parameter `bj`, where `other` is the other one
fn lower_recurrence(a: [f64; 3], bj: f64, other: f64, z: f64, m: u32) -> Option<(Dd, f64, f64)> {
    let checked = pfq_extreme(&a, &[bj, other, 1.0]);
    let start = |k: f64| {
        // the terms for `u` itself, which `pfq_terms` would round
        let u = Dd::new(bj) + k;
        let ratio = |n| pfq_ratio::<Dd>(&a, &[other, 1.0], n, checked) / (u + n);
        recurrence_start(a, [u.hi, other, 1.0], z, ratio)
    };
    let step = |k: f64| {
        // F(b_j + k) from the relation at b = b_j + k + 3
        let u = Dd::new(bj) + (k + 2.0);
        let size = u.hi.abs();
        (
            lower_coefficients(u, a.map(|ai| ai - u), Dd::new(z), Dd::new(other), -1.0),
            lower_coefficients(size, a.map(|ai| ai.abs() + size), z.abs(), other.abs(), 1.0),
        )
    };
    recurrence_from(m, start, step)
}

/// `upper_coefficients` of the relation in the upper parameter `a[i]` at `a = w - 1`, and the
/// sizes of their terms
fn upper_step(a: [f64; 3], i: usize, b: [f64; 2], z: f64, w: Dd) -> ([Dd; 4], [f64; 4]) {
    let others = [a[(i + 1) % 3], a[(i + 2) % 3]];
    let size = |v: f64| v.abs() + w.hi.abs();
    (
        upper_coefficients(w, b.map(|bj| bj - w), others.map(|ai| w - ai), Dd::new(z), -1.0),
        upper_coefficients(w.hi.abs(), b.map(size), others.map(size), z.abs(), 1.0),
    )
}

/// `recurrence_from` down in the upper parameter `a[i]`, where `b` are the lower ones
fn upper_recurrence(a: [f64; 3], i: usize, b: [f64; 2], z: f64, m: u32) -> Option<(Dd, f64, f64)> {
    let (others, b1) = ([a[(i + 1) % 3], a[(i + 2) % 3]], [b[0], b[1], 1.0]);
    let checked = pfq_extreme(&a, &b1);
    let start = |k: f64| {
        // the terms for `u` itself, as for a lower parameter
        let u = Dd::new(a[i]) + k;
        let ratio = |n| pfq_ratio::<Dd>(&others, &b1, n, checked) * (u + n);
        recurrence_start([u.hi, others[0], others[1]], b1, z, ratio)
    };
    // F(a_i + k) from the relation at a = a_i + k
    recurrence_from(m, start, |k| upper_step(a, i, b, z, Dd::new(a[i]) + (k + 1.0)))
}

/// `y` as `f64` and the estimate `err` of its absolute error plus its rounding, if that is
/// accurate
fn accurate(y: Dd, err: f64) -> Option<(f64, f64)> {
    let y = (y.hi, err + 0.5 * f64::EPSILON * y.hi.abs());
    Estimate::absolute(y).is_accurate().then_some(y)
}

/// The first accurate result of `run(m)` for the `m` of `starts`, which ascend, and the estimate
/// of its absolute error
fn recurrence_first(
    starts: impl IntoIterator<Item = u32>,
    run: impl Fn(u32) -> Option<(Dd, f64, f64)>,
) -> Option<(f64, f64)> {
    for m in starts {
        // (a series that is not accurate may be for a larger `m`)
        let Some((y, rounding, err)) = run(m) else {
            continue;
        };
        if let Some(result) = accurate(y, err) {
            return Some(result);
        }
        // a larger `m` ends with the same steps, after more of them
        if accurate(y, rounding).is_none() {
            break;
        }
    }
    None
}

/// 3F2 for `|z| < 1` by the recurrence in a parameter, and an estimate of its absolute error,
/// or `None` if it is not 3F2 or not accurate
pub(super) fn pfq_recurrence<P: Params>(a: &P, b: &P, z: f64) -> Option<(f64, f64)> {
    let a: [f64; 3] = a.as_ref().try_into().ok()?;
    // `b` is in descending order, and has the 1 of `n!`
    let &([1.0, large, small] | [large, 1.0, small] | [large, small, 1.0]) = b.as_ref() else {
        return None;
    };
    let lower =
        |(bj, other)| recurrence_first(RECURRENCE_STARTS, |m| lower_recurrence(a, bj, other, z, m));
    let upper = |i: usize| {
        // from the series where a_i + m is about -4096, ..., -8, the fewest steps first (after
        // `steps` of them, it would be in (0, 1))
        let steps = (-a[i]).ceil() as u32;
        let offsets = RECURRENCE_STARTS.iter().rev().filter(|&&p| p < steps);
        let run = |m| upper_recurrence(a, i, [large, small], z, m);
        recurrence_first(offsets.map(|p| steps - p), run)
    };
    // first in the smaller lower parameter, and last in the negative upper ones, in at most as
    // many steps as the series has terms
    let mut negative = (0..3).filter(|&i| (-f64::from(PFQ_MAXITER)..0.0).contains(&a[i]));
    let lowers = [(small, large), (large, small)].into_iter().find_map(lower);
    lowers.or_else(|| negative.find_map(upper))
}

#[cfg(test)]
mod tests {
    use super::super::pfq_series;
    use super::*;

    /// Asserts that the result is accurate, and that it differs from the series, which is
    /// accurate too, by at most their estimated errors
    fn assert_estimate(result: Option<(f64, f64)>, a: [f64; 3], b: [f64; 3], z: f64) {
        let ((y, err), (expected, series_err)) = (result.unwrap(), pfq_series(&a, &b, z));
        assert!(err <= 1e-13 * y.abs(), "({a:?}, {b:?}, {z}): {err}");
        assert!((y - expected).abs() <= err + series_err, "{y} != {expected} ({err})");
    }

    #[test]
    fn test_pfq_recurrence_estimate() {
        for (a, b, z) in [
            ([1.5, 2.5, 3.5], [4.25, 1.0, 0.75], -0.5),
            ([10.5, -7.25, 0.5], [3.5, 1.0, -2.75], 0.3),
            ([30.5, 20.25, 0.125], [2.5, 1.0, -0.999_999], -0.9),
            // z > 1/2, where the errors grow in each step
            ([0.5, 1.25, -0.75], [12.5, 2.5, 1.0], 0.75),
        ] {
            assert_estimate(pfq_recurrence(&a, &b, z), a, b, z);
        }
    }

    #[test]
    fn test_upper_recurrence_estimate() {
        // in `m` steps down in the last upper parameter
        for (a, [b1, b2], z, m) in [
            ([10.5, 0.5, -20.25], [3.5, -2.75], 0.3, 12),
            ([1.25, 0.75, -40.5], [2.5, 1.5], 0.75, 33),
            ([2.5, 1.5, -30.5], [4.25, 0.75], -0.5, 23),
            ([-3.5, -12.25, -24.75], [-1.5, -30.5], 0.9, 17),
        ] {
            let result = recurrence_first([m], |m| upper_recurrence(a, 2, [b1, b2], z, m));
            assert_estimate(result, a, descending([b1, b2, 1.0]), z);
        }
    }
}
