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

use core::iter::{successors, zip};

use super::{
    Dd, PFQ_DD_MAXITER, Params, SeriesReal, accurate, descending, first_zero, pfq_extreme,
    pfq_ratio, pfq_sum_terms, pfq_terms_of,
};

/// Largest `m` for which the recurrence starts from the series at `b + m`
///
/// Of 7136 results on 91004 points, 426 are from it and 537 from half of it, and twice it gave
/// 8 more results on 1351 points (measured).
const RECURRENCE_MAX_START: f64 = 4096.0;

/// Bound of the rounding error of a step, relative to the sizes of the terms of its sum:
/// `Dd::EPSILON` for each of its at most 17 inexact operations in a row (`u - 1`, `u - 2`,
/// `z - 1` and the products with a sign are exact)
const RECURRENCE_EPSILON: f64 = 17.0 * Dd::EPSILON;

/// The coefficients `c_0, ..., c_3` for `sign = -1`. For `sign = 1` and the sizes of the
/// arguments, each is the sum of the sizes of its terms instead, which do not cancel.
fn coefficients<T: SeriesReal>(u: T, [d1, d2, d3]: [T; 3], z: T, other: T, sign: f64) -> [T; 4] {
    let (s, sum) = (T::from_f64(sign), d1 + d2 + d3);
    // u - 1 and u - 2
    let (u1, u2) = (u + sign, u + 2.0 * sign);
    [
        d1 * d2 * d3 * z,
        u * ((d1 * d2 + d1 * d3 + d2 * d3 + sum + 1.0) * z + s * u1 * (u + s * other)),
        u * u1 * ((sum + 3.0) * z + u1 + u1 + s * other),
        u * u1 * u2 * (z + sign),
    ]
}

/// 3F2 for `|z| < 1` by the recurrence in a lower parameter `b_j`, and an estimate of its
/// absolute error, or `None` if it is not 3F2 or the series that start it are not accurate
///
/// It starts at `b_j + m` for `m = 8, 16, ...`, until the result is accurate or its rounding
/// errors alone are too large, first for the smaller `b_j`.
pub(super) fn pfq_recurrence<P: Params>(a: &P, b: &P, z: f64) -> Option<(f64, f64)> {
    let (&[a1, a2, a3], &[_, _, _]) = (a.as_ref(), b.as_ref()) else {
        return None;
    };
    let one = b.as_ref().iter().position(|&v| v == 1.0)?;
    // b_j and the other lower parameter
    let lower = |j: usize| (b.as_ref()[j], b.as_ref()[3 - one - j]);
    let (degree, checked) = (first_zero(a.as_ref()), pfq_extreme(a.as_ref(), b.as_ref()));
    // F at `b_j + k` by the series in double-double and its absolute error, if it is accurate
    let series = |j: usize, k: f64| {
        let (bj, other) = lower(j);
        let u = Dd::new(bj) + k;
        // the terms for `u` itself, which `pfq_terms` would round
        let ratio = |n| pfq_ratio::<Dd>(a.as_ref(), &[other, 1.0], n, checked) / (u + n);
        let terms = pfq_terms_of(degree, z, ratio);
        let mut rounded = b.clone();
        rounded.as_mut()[j] = u.hi;
        let rounded = descending(rounded);
        // (summed to the accuracy of double-double, as the steps can make its error larger)
        let (y, err, _) = pfq_sum_terms(a, &rounded, z, PFQ_DD_MAXITER, terms, Dd::EPSILON);
        accurate((y.hi, err)).map(|_| (y, err))
    };
    // F at `b_j` from the series at `b_j + m + k`, the part of the estimate of its error from
    // the rounding errors, and that estimate
    let run = |j: usize, m: f64| {
        let (bj, other) = lower(j);
        let seeds = [series(j, m)?, series(j, m + 1.0)?, series(j, m + 2.0)?];
        let mut f = seeds.map(|(y, _)| y);
        // of each step, F(b_j + k) in terms of the three above it, and its rounding error
        let mut steps = Vec::with_capacity(m as usize);
        for k in (0..m as u32).rev().map(f64::from) {
            // the relation at b = b_j + k + 3
            let u = Dd::new(bj) + (k + 2.0);
            let d = [a1, a2, a3].map(|ai| ai - u);
            let [c0, c1, c2, c3] = coefficients(u, d, Dd::new(z), Dd::new(other), -1.0);
            // the sizes of the terms of c_0, c_1 and c_2, which bound their rounding errors
            let size = u.hi.abs();
            let sizes = [a1, a2, a3].map(|ai| ai.abs() + size);
            let [e0, e1, e2, _] = coefficients(size, sizes, z.abs(), other.abs(), 1.0);
            let terms = e0 * f[2].hi.abs() + e1 * f[1].hi.abs() + e2 * f[0].hi.abs();
            let next = -(c0 * f[2] + c1 * f[1] + c2 * f[0]) / c3;
            // an underflow adds at most `f64::MIN_POSITIVE` to the sum, to the products c_0
            // and c_3, and to the quotient; in c_1 and c_2, less than the rounding of their
            // terms without z
            let factors = 1.0 + f[2].hi.abs() + next.hi.abs() + c3.hi.abs();
            let rounding = (RECURRENCE_EPSILON * terms + f64::MIN_POSITIVE * factors) / c3.hi.abs();
            steps.push(([c2, c1, c0].map(|c| -c.hi / c3.hi), rounding));
            f = [next, f[0], f[1]];
        }
        // the factors of the errors of the last three F(b_j + k) in the result, upwards, and
        // the bound of the rounding errors, with that of the result
        let (mut g, mut rounding) = ([1.0_f64, 0.0, 0.0], 0.5 * f64::EPSILON * f[0].hi.abs());
        for (row, step) in steps.iter().rev() {
            rounding += g[0].abs() * step;
            g = [g[0] * row[0] + g[1], g[0] * row[1] + g[2], g[0] * row[2]];
        }
        let series_err: f64 = zip(g, seeds).map(|(g, (_, err))| g.abs() * err).sum();
        Some((f[0].hi, rounding, rounding + series_err))
    };
    let starts = || {
        let doubled = successors(Some(8.0), |m| Some(2.0 * m));
        doubled.take_while(|&m| m <= RECURRENCE_MAX_START)
    };
    // the smaller lower parameter first, as `b` is in descending order
    let mut results = (0..3).rev().filter(|&j| j != one).flat_map(|j| {
        let runs = starts().filter_map(move |m| run(j, m));
        // a farther start ends with the same steps, after more of them
        runs.take_while(|&(y, rounding, _)| accurate((y, rounding)).is_some())
    });
    results.find_map(|(y, _, err)| accurate((y, err)).map(|y| (y, err)))
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
