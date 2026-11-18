//! 3F2 for `|z| < 1` by a recurrence in a lower parameter `b`
//!
//! Where the series cancels, e.g. for large `a_i` and `z < 0`, it may not for `b + m` with a
//! large `m`. With `(theta + b - 1) F(b) = (b - 1) F(b - 1)` for `theta = z d/dz` (DLMF 16.3.4
//! with `n = 1`), the differential equation (DLMF 16.8.3) becomes
//!
//! ```text
//! c_0 F(b) + c_1 F(b - 1) + c_2 F(b - 2) + c_3 F(b - 3) = 0,
//! c_0 = z d_1 d_2 d_3,
//! c_1 = u (z (d_1 d_2 + d_1 d_3 + d_2 d_3 + d_1 + d_2 + d_3 + 1) - (u - 1) (u - b')),
//! c_2 = u (u - 1) (z (d_1 + d_2 + d_3 + 3) + 2 u - 2 - b'),
//! c_3 = u (u - 1) (u - 2) (z - 1),
//! ```
//!
//! with `u = b - 1`, `d_i = a_i - u` and the other lower parameter `b'`. For `b -> inf`, `F(b)`
//! tends to 1, and the solutions grow per step down by about 1, 1 and `z / (z - 1)`, so the
//! recurrence is unstable downwards for `z > 1/2`. For `b` below the `a_i`, other solutions can
//! also grow by many orders of magnitude.
//!
//! It is run in double-double from the series at `b + m + k` for `k = 0, 1, 2`. All `b + k` are
//! exact in double-double, so that the only errors are those of the three series and the
//! rounding errors of the steps. The result is linear in them: an error `e` of the `F(b + k)`
//! of a step adds `g_k e` to it, where the `g_k` follow from the adjoint recurrence, upwards
//! from `g_0 = 1`. The estimate of the absolute error of the result is the sum of the `|g_k|`
//! times the bounds of the rounding errors and the estimates of the errors of the series. It
//! is not a bound: the `g_k` are computed in `f64` from the rounded coefficients, and can be
//! too small where they cancel.

use core::iter::zip;

use super::{
    Dd, PFQ_DD_MAXITER, Params, SeriesReal, accurate, descending, first_zero, pfq_extreme,
    pfq_ratio, pfq_sum_terms, pfq_terms_of,
};

/// Smallest `m` for which the recurrence starts from the series at `b + m`. The next ones are
/// `2 m`, `4 m`, ...
const RECURRENCE_MIN_START: u32 = 8;

/// Largest such `m`
///
/// Of 7136 results on 91004 points, 426 are from it and 537 from half of it, and twice it gave
/// 8 more results on 1351 points (measured).
const RECURRENCE_MAX_START: u32 = 4096;

/// Bound of the rounding error of a step, relative to the sizes of the terms of its sum:
/// `Dd::EPSILON` for each of its at most 17 inexact operations in a row (`u - 1`, `u - 2`,
/// `z - 1` and the products with `minus` are exact)
const RECURRENCE_EPSILON: f64 = 17.0 * Dd::EPSILON;

/// The coefficients `[c_0, c_1, c_2, c_3]`, for `minus = -1`
///
/// For `minus = 1` and the sizes of the arguments, no terms cancel: each is then the sum of the
/// sizes of its terms, which bounds its rounding error.
fn coefficients<T: SeriesReal>(u: T, [d1, d2, d3]: [T; 3], z: T, other: T, minus: f64) -> [T; 4] {
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

/// F at `b_j + k` by the series in double-double and an estimate of its absolute error, or
/// `None` if it is not accurate
fn recurrence_start(a: [f64; 3], bj: f64, other: f64, z: f64, k: f64) -> Option<(Dd, f64)> {
    let u = Dd::new(bj) + k;
    // the terms for `u` itself, which `pfq_terms` would round
    let checked = pfq_extreme(&a, &[bj, other, 1.0]);
    let ratio = |n| pfq_ratio::<Dd>(&a, &[other, 1.0], n, checked) / (u + n);
    let terms = pfq_terms_of(first_zero(&a), z, ratio);
    let b = descending([u.hi, other, 1.0]);
    // (summed to the accuracy of double-double, as the steps can make its error larger)
    let (y, err, _) = pfq_sum_terms(&a, &b, z, PFQ_DD_MAXITER, terms, Dd::EPSILON);
    accurate((y.hi, err))?;
    Some((y, err))
}

/// F at `b_j` in `m` steps from the series at `b_j + m + k` for `k = 0, 1, 2`, the part of
/// the estimate of its absolute error from the rounding errors, and that estimate, or `None`
/// if one of the series is not accurate
fn recurrence_from(a: [f64; 3], bj: f64, other: f64, z: f64, m: u32) -> Option<(f64, f64, f64)> {
    let start = |k: u32| recurrence_start(a, bj, other, z, f64::from(m + k));
    let starts = [start(0)?, start(1)?, start(2)?];
    // F at the three `b` above that of the next step, in ascending order
    let mut f = starts.map(|(y, _)| y);
    // of each step: the factors of those three in its F, and the bound of its rounding error
    let mut steps = Vec::with_capacity(m as usize);
    for k in (0..m).rev().map(f64::from) {
        // F(b_j + k) from the relation at b = b_j + k + 3
        let u = Dd::new(bj) + (k + 2.0);
        let d = a.map(|ai| ai - u);
        let [c0, c1, c2, c3] = coefficients(u, d, Dd::new(z), Dd::new(other), -1.0);
        let sum = c0 * f[2] + c1 * f[1] + c2 * f[0];
        let next = -sum / c3;
        // the sizes of the terms of c_0, c_1 and c_2, and with them those of the sum
        let size = u.hi.abs();
        let sizes = a.map(|ai| ai.abs() + size);
        let [e0, e1, e2, _] = coefficients(size, sizes, z.abs(), other.abs(), 1.0);
        let terms = e0 * f[2].hi.abs() + e1 * f[1].hi.abs() + e2 * f[0].hi.abs();
        // an underflow adds at most `f64::MIN_POSITIVE` to the sum, to the products c_0
        // and c_3, and to the quotient; in c_1 and c_2, less than the rounding of their
        // terms without z
        let factors = 1.0 + f[2].hi.abs() + next.hi.abs() + c3.hi.abs();
        let rounding = (RECURRENCE_EPSILON * terms + f64::MIN_POSITIVE * factors) / c3.hi.abs();
        steps.push(([c2, c1, c0].map(|c| -c.hi / c3.hi), rounding));
        f = [next, f[0], f[1]];
    }
    // the `g_k` of the three F above each step, upwards from the result, and the bound of the
    // rounding errors, which starts with that of the result
    let (mut g, mut rounding) = ([1.0_f64, 0.0, 0.0], 0.5 * f64::EPSILON * f[0].hi.abs());
    for (factors, step) in steps.iter().rev() {
        rounding += g[0].abs() * step;
        g = [g[0] * factors[0] + g[1], g[0] * factors[1] + g[2], g[0] * factors[2]];
    }
    let series: f64 = zip(g, starts).map(|(g, (_, err))| g.abs() * err).sum();
    Some((f[0].hi, rounding, rounding + series))
}

/// 3F2 for `|z| < 1` by the recurrence in a lower parameter, and an estimate of its absolute
/// error, or `None` if it is not 3F2 or not accurate
pub(super) fn pfq_recurrence<P: Params>(a: &P, b: &P, z: f64) -> Option<(f64, f64)> {
    let a: [f64; 3] = a.as_ref().try_into().ok()?;
    // `b` is in descending order, and has the 1 of `n!`
    let &([1.0, large, small] | [large, 1.0, small] | [large, small, 1.0]) = b.as_ref() else {
        return None;
    };
    // first in the smaller one
    for (bj, other) in [(small, large), (large, small)] {
        let mut m = RECURRENCE_MIN_START;
        while m <= RECURRENCE_MAX_START {
            // (`None` if a series is not accurate, which it may be for a larger `m`)
            if let Some((y, rounding, err)) = recurrence_from(a, bj, other, z, m) {
                if accurate((y, err)).is_some() {
                    return Some((y, err));
                }
                // a larger `m` ends with the same steps, after more of them
                if accurate((y, rounding)).is_none() {
                    break;
                }
            }
            m *= 2;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::super::pfq_series;
    use super::*;

    #[test]
    fn test_pfq_recurrence_estimate() {
        // where the series is accurate too, they differ by at most their estimated errors
        for (a, b, z) in [
            ([1.5, 2.5, 3.5], [4.25, 1.0, 0.75], -0.5),
            ([10.5, -7.25, 0.5], [3.5, 1.0, -2.75], 0.3),
            ([30.5, 20.25, 0.125], [2.5, 1.0, -0.999_999], -0.9),
            // z > 1/2, where the errors grow in each step
            ([0.5, 1.25, -0.75], [12.5, 2.5, 1.0], 0.75),
        ] {
            let (y, err) = pfq_recurrence(&a, &b, z).unwrap();
            let (expected, series_err) = pfq_series(&a, &b, z);
            assert!(err <= 1e-13 * y.abs(), "pfq_recurrence({a:?}, {b:?}, {z}): {err}");
            assert!((y - expected).abs() <= err + series_err, "{y} != {expected} ({err})");
        }
    }
}
