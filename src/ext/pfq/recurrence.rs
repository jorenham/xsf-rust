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
//! For `z < 0` it also cancels at large upper parameters `a > 0`, and less at `a - m`. The last
//! relation then gives `F(a)` in `m` steps up from the series at `a - m`, `a - m - 1` and
//! `a - m - 2`. The limiting characteristic roots of a step up, for large `w` and fixed other
//! parameters, are 1, 1 and `1 / (1 - z)`: none is above 1, but the repeated root at 1 leaves
//! room for the errors to grow. For more than one large `a_i`, the starts of the steps in one
//! come from the steps in the next, not from the series.
//!
//! For a large negative `b` and `z < 0`, the steps down from `b + m > 0` can make the errors
//! larger by many orders of magnitude where `b + k` is positive, and much less where it is
//! negative. They then start where `b + m` is about -8, -16, ..., from the steps up in the large
//! upper parameters there. Where the series cancel there too, as for a negative other lower
//! parameter, those starts come from the steps down in that one from where it is positive.
//!
//! The steps can make the errors of the series and their own rounding errors larger by many
//! orders of magnitude, so the error of the result is estimated. An error `e` of the F at `k`
//! steps from the result adds `g_k e` to it, where the `g_k` follow from the adjoint recurrence,
//! upwards from `g_0 = 1`. The estimate is not a bound, as the `g_k` are computed in `f64` from
//! the rounded coefficients.

use core::array::from_fn;
use core::cell::Cell;
use core::iter::{once, zip};

use super::{
    Dd, Estimate, PFQ_DD_MAXITER, PFQ_LONG_MAXITER, PFQ_MAXITER, PFQ_TOL, Params, SeriesReal,
    descending, first_zero, pfq_extreme, pfq_ratio, pfq_sum_terms, pfq_terms_of,
};

/// The `m` for which the recurrence steps down from the series at `b + m`, in the order in which
/// they are tried; in a negative parameter `v`, it steps down from where `v + m` is about minus
/// these, in reverse
const RECURRENCE_STARTS: [u32; 10] = [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];

/// Smallest upper parameter `a` that the recurrence goes up in: in a step from `a - 1`, `a - 2`
/// and `a - 3`, which are not negative, as a step up divides by the `c_0` of the lowest one,
/// which is small near -1 and -2
const RAISED_MIN: f64 = 3.0;

/// Most steps of the recurrences up in the upper parameters plus terms of the series that they
/// start from, of all numbers of steps that are tried together: up in three of them, each try
/// is 13 recurrences from 27 series, which can take much longer than the other methods
///
/// Of 1927 results on 126504 points, 1880 are within it, and twice as much gives 37 more
/// (measured).
const RAISED_WORK: u32 = 2 * PFQ_MAXITER;

/// Most steps and series terms of all tries of `recurrence_balanced` together: the terms of
/// three starts of `PFQ_LONG_MAXITER` each
///
/// Of its 82 results without it on 1059 points where the other methods are not accurate, 80 are
/// within it (measured).
const BALANCED_WORK: u32 = 3 * PFQ_LONG_MAXITER;

/// Most steps and series terms of all tries of `recurrence_nested` together
///
/// Of its 127 results without it on 979 points where the other methods are not accurate, 108 are
/// within it, and half as much again gives 8 more (measured).
const NESTED_WORK: u32 = 4 * PFQ_MAXITER;

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
/// it is not accurate within `maxiter` terms or `work` runs out (each term takes one)
fn recurrence_start(
    a: [f64; 3],
    b: [f64; 2],
    z: f64,
    maxiter: u32,
    work: &Cell<u32>,
    ratio: impl Fn(f64) -> Dd,
) -> Option<(Dd, f64)> {
    let limit = maxiter.min(work.get());
    let terms = pfq_terms_of(first_zero(&a), z, |n, t: Dd| {
        work.set(work.get() - 1);
        t * (ratio(n) * Dd::new(z))
    });
    let (a, b) = (descending(a), descending(b));
    // (summed to the accuracy of double-double, as the steps can make its error larger)
    let (y, err, _) = pfq_sum_terms(&a, &b, z, limit, terms, Dd::EPSILON);
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

/// `recurrence_from` in the lower parameter `bj`, where `other` is the other one, or `None` if a
/// series is not accurate or `work` runs out, of which each step and term takes one
fn lower_recurrence(
    a: [f64; 3],
    bj: f64,
    other: f64,
    z: f64,
    m: u32,
    work: &Cell<u32>,
) -> Option<(Dd, f64, f64)> {
    let checked = pfq_extreme(&a, &[bj, other]);
    work.set(work.get().checked_sub(m)?);
    let start = |k: f64| {
        // the terms for `u` itself, which `pfq_terms` would round
        let u = Dd::new(bj) + k;
        let ratio = |n| pfq_ratio::<Dd>(&a, &[other], n, checked) / (u + n);
        recurrence_start(a, [u.hi, other], z, PFQ_DD_MAXITER, work, ratio)
    };
    lower_steps(a, bj, other, z, m, start)
}

/// `recurrence_from` in the lower parameter `bj` from `start(k)`, the F at `bj + k`
fn lower_steps(
    a: [f64; 3],
    bj: f64,
    other: f64,
    z: f64,
    m: u32,
    start: impl Fn(f64) -> Option<(Dd, f64)>,
) -> Option<(Dd, f64, f64)> {
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

/// `recurrence_from` down in the upper parameter `a[i]`, where `b` are the lower ones, from series
/// of at most `maxiter` terms, or `None` if a series is not accurate or `work` runs out, of which
/// each step and term takes one
fn upper_recurrence(
    a: [f64; 3],
    i: usize,
    b: [f64; 2],
    z: f64,
    m: u32,
    maxiter: u32,
    work: &Cell<u32>,
) -> Option<(Dd, f64, f64)> {
    let others = [a[(i + 1) % 3], a[(i + 2) % 3]];
    let checked = pfq_extreme(&a, &b);
    work.set(work.get().checked_sub(m)?);
    let start = |k: f64| {
        // the terms for `u` itself, as for a lower parameter
        let u = Dd::new(a[i]) + k;
        let ratio = |n| pfq_ratio::<Dd>(&others, &b, n, checked) * (u + n);
        recurrence_start([u.hi, others[0], others[1]], b, z, maxiter, work, ratio)
    };
    // F(a_i + k) from the relation at a = a_i + k
    recurrence_from(m, start, |k| upper_step(a, i, b, z, Dd::new(a[i]) + (k + 1.0)))
}

/// F and an estimate of its absolute error by `recurrence_from` up in each upper parameter from
/// `a[i]` on that is at least `RAISED_MIN`, in at most `m` steps, or by the series if there is
/// none, where `b` are the lower parameters
///
/// It is `None` if a series is not accurate, or if `work` runs out: what is left of
/// `RAISED_WORK`.
fn raised_recurrence(
    a: [f64; 3],
    i: usize,
    b: [f64; 2],
    z: f64,
    m: u32,
    work: &Cell<u32>,
) -> Option<(Dd, f64)> {
    let Some(i) = (i..3).find(|&i| a[i] >= RAISED_MIN) else {
        let checked = pfq_extreme(&a, &b);
        // near z = -1, the series can take more than `PFQ_DD_MAXITER` terms, up to all work left
        let ratio = |n| pfq_ratio::<Dd>(&a, &b, n, checked);
        return recurrence_start(a, b, z, u32::MAX, work, ratio);
    };
    // (so that the lowest start is not negative)
    let m = m.min(a[i] as u32 - 2);
    work.set(work.get().checked_sub(m)?);
    let start = |k: f64| {
        // exact: a_i < 2^53 and the integer k <= a_i are multiples of the spacing of the
        // `f64` at a_i, and so is a_i - k, which is in [0, a_i]
        let lowered = from_fn(|j| if j == i { a[i] - k } else { a[j] });
        raised_recurrence(lowered, i + 1, b, z, m, work)
    };
    let step = |k: f64| {
        // F(a_i - k) from the relation at a = a_i - k - 3, whose F are in the reverse order
        let (mut c, mut sizes) = upper_step(a, i, b, z, Dd::new(a[i]) - (k + 2.0));
        c.reverse();
        sizes.reverse();
        (c, sizes)
    };
    let (y, _, err) = recurrence_from(m, start, step)?;
    Some((y, err))
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
) -> Option<(Dd, f64)> {
    for m in starts {
        // (a series that is not accurate may be for a larger `m`)
        let Some((y, rounding, err)) = run(m) else {
            continue;
        };
        if accurate(y, err).is_some() {
            return Some((y, err));
        }
        // a larger `m` ends with the same steps, after more of them
        if accurate(y, rounding).is_none() {
            break;
        }
    }
    None
}

/// The `ceil(-v)` steps after which `v` would be in (0, 1), or 0 if `v >= 0` or if they are more
/// than `PFQ_MAXITER`
fn steps_to_zero(v: f64) -> u32 {
    // (the cast saturates at 0 for v >= 0)
    let steps = (-v).ceil() as u32;
    if steps <= PFQ_MAXITER { steps } else { 0 }
}

/// The `m` for which `v + m` is about -4096, ..., -8, the fewest first, if `v` is in reach of
/// `PFQ_MAXITER` steps
fn toward_zero(v: f64) -> impl DoubleEndedIterator<Item = u32> {
    let steps = steps_to_zero(v);
    let offsets = RECURRENCE_STARTS.iter().rev().filter(move |&&p| p < steps);
    offsets.map(move |p| steps - p)
}

/// F by `raised_recurrence` for `z < 0`, first in all steps of the largest upper parameter,
/// `a[0]`, then in fewer, which make the errors less large and the series cancel more, while
/// there is `work` left, and the estimate of its absolute error, if that is accurate
fn raised_first(a: [f64; 3], b: [f64; 2], z: f64, work: &Cell<u32>) -> Option<(Dd, f64)> {
    if z >= 0.0 || !(RAISED_MIN..f64::from(PFQ_MAXITER)).contains(&a[0]) {
        return None;
    }
    let longest = a[0] as u32 - 2;
    let fewer = RECURRENCE_STARTS.into_iter().rev().filter(|&m| m < longest);
    let tries = once(longest).chain(fewer).take_while(|_| work.get() > 0);
    let run = |m| raised_recurrence(a, 0, b, z, m, work);
    let mut results = tries.filter_map(run);
    results.find(|&(y, err)| accurate(y, err).is_some())
}

/// The upper and lower parameters of a 3F2 in descending order, or `None` if it is not 3F2
fn params_3f2<A: Params, B: Params>(a: &A, b: &B) -> Option<([f64; 3], [f64; 2])> {
    Some((a.as_ref().try_into().ok()?, b.as_ref().try_into().ok()?))
}

/// 3F2 for `|z| < 1` by the recurrence in a parameter, and an estimate of its absolute error,
/// or `None` if it is not 3F2 or not accurate
pub(super) fn pfq_recurrence<A: Params, B: Params>(a: &A, b: &B, z: f64) -> Option<(f64, f64)> {
    let (a, [large, small]) = params_3f2(a, b)?;
    let lower = |(bj, other)| {
        let run = |m| lower_recurrence(a, bj, other, z, m, &Cell::new(u32::MAX));
        recurrence_first(RECURRENCE_STARTS, run)
    };
    let raised = || raised_first(a, [large, small], z, &Cell::new(RAISED_WORK));
    // first in the smaller lower parameter, then the larger, and last up in the large upper ones
    // for z < 0
    let (y, err) = [(small, large), (large, small)]
        .into_iter()
        .find_map(lower)
        .or_else(raised)?;
    accurate(y, err)
}

/// 3F2 for `|z| < 1` by the slower recurrences (lowered, balanced upper, and nested), and an
/// estimate of its absolute error, or `None` if it is not 3F2 or not accurate
pub(super) fn pfq_recurrence_extended<A: Params, B: Params>(
    a: &A,
    b: &B,
    z: f64,
) -> Option<(f64, f64)> {
    let (a, b) = params_3f2(a, b)?;
    recurrence_lowered(a, b, z)
        .or_else(|| recurrence_balanced(a, b, z))
        .or_else(|| recurrence_nested(a, b, z))
}

/// Steps down in a negative smaller lower parameter from where `raised_first` gives F
fn recurrence_lowered(a: [f64; 3], [large, small]: [f64; 2], z: f64) -> Option<(f64, f64)> {
    if z >= 0.0 || a[0] < RAISED_MIN {
        return None;
    }
    // one `RAISED_WORK` for each of the three starts, shared by all `m`
    let work = Cell::new(3 * RAISED_WORK);
    // (exact, as `small + k` is in `(small, 0)`)
    let start = |k: f64| raised_first(a, [small + k, large], z, &work);
    let (y, err) =
        recurrence_first(toward_zero(small), |m| lower_steps(a, small, large, z, m, start))?;
    accurate(y, err)
}

/// Steps down in a negative smaller lower parameter from where it is about -8, -16, ..., whose
/// starts are by the steps down in the other one from where that is positive
///
/// It tries the starts nearest to 0 first, where the series that their F come from cancel least.
/// All tries share `NESTED_WORK`.
fn recurrence_nested(a: [f64; 3], [large, small]: [f64; 2], z: f64) -> Option<(f64, f64)> {
    if z >= 0.0 {
        return None;
    }
    let work = Cell::new(NESTED_WORK);
    let start = |k: f64| {
        let positive = RECURRENCE_STARTS
            .into_iter()
            .filter(|&m| f64::from(m) > -large);
        // (exact, as `small + k` is in `(small, 0)`)
        recurrence_first(positive, |m| lower_recurrence(a, large, small + k, z, m, &work))
    };
    let run = |m| {
        work.set(work.get().checked_sub(m)?);
        lower_steps(a, small, large, z, m, start)
    };
    let mut results = toward_zero(small).rev().filter_map(run);
    results.find_map(|(y, _, err)| accurate(y, err))
}

/// Steps down in a negative upper parameter from about the most steps whose rounding errors are
/// within half of the tolerance
///
/// With fewer steps, the series at their start cancels more, and with more, their rounding errors
/// grow more. The series have up to `PFQ_LONG_MAXITER` terms, as near `z = 1` they can need more
/// than `PFQ_DD_MAXITER` where `a_i + m` is near 0. All of them and the steps share
/// `BALANCED_WORK`.
///
/// The search is a heuristic: it bisects as if the estimate of the rounding errors grew with the
/// steps and only too few steps made a start inaccurate, which the long series near `z = 1` and
/// the rounding of the starts can break, and it returns the first accurate result it finds.
fn recurrence_balanced(a: [f64; 3], b: [f64; 2], z: f64) -> Option<(f64, f64)> {
    let work = Cell::new(BALANCED_WORK);
    let balanced = |i| {
        let run = |m, maxiter| upper_recurrence(a, i, b, z, m, maxiter, &work);
        // try near 0 (`high - 8`) first with a short `PFQ_DD_MAXITER` budget so large negative
        // `a[i]` finish in one attempt without starving bisection, then bisect `[1, steps_to_zero)`
        let (mut low, mut high) = (0, steps_to_zero(a[i]));
        if high > 2 * RECURRENCE_STARTS[0] {
            let m = high - RECURRENCE_STARTS[0];
            if let Some((y, rounding, err)) = run(m, PFQ_DD_MAXITER) {
                if let Some(y) = accurate(y, err) {
                    return Some(y);
                }
                if 2.0 * rounding <= PFQ_TOL * y.hi.abs() {
                    low = m;
                } else {
                    high = m;
                }
            }
        }
        while low + 1 < high {
            let m = low.midpoint(high);
            let y = run(m, PFQ_LONG_MAXITER);
            if let Some(y) = y.and_then(|(y, _, err)| accurate(y, err)) {
                return Some(y);
            }
            // within half of the tolerance, which leaves the other half to the series
            let within = y.is_none_or(|(y, rounding, _)| 2.0 * rounding <= PFQ_TOL * y.hi.abs());
            if within { low = m } else { high = m }
        }
        None
    };
    (0..3).find_map(balanced)
}

#[cfg(test)]
mod tests {
    use super::super::pfq_series;
    use super::*;

    /// Asserts that the result is accurate, and that it differs from the series, which is
    /// accurate too, by at most their estimated errors
    fn assert_estimate(result: Option<(f64, f64)>, a: [f64; 3], b: [f64; 2], z: f64) {
        let ((y, err), (expected, series_err)) = (result.unwrap(), pfq_series(&a, &b, z));
        assert!(err <= 1e-13 * y.abs(), "({a:?}, {b:?}, {z}): {err}");
        assert!((y - expected).abs() <= err + series_err, "{y} != {expected} ({err})");
    }

    #[test]
    fn test_pfq_recurrence_estimate() {
        for (a, b, z) in [
            ([1.5, 2.5, 3.5], [4.25, 0.75], -0.5),
            ([10.5, -7.25, 0.5], [3.5, -2.75], 0.3),
            ([30.5, 20.25, 0.125], [2.5, -0.999_999], -0.9),
            // z > 1/2, where the errors grow in each step
            ([0.5, 1.25, -0.75], [12.5, 2.5], 0.75),
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
            let result = recurrence_first([m], |m| {
                upper_recurrence(a, 2, [b1, b2], z, m, PFQ_DD_MAXITER, &Cell::new(u32::MAX))
            });
            let result = result.and_then(|(y, err)| accurate(y, err));
            assert_estimate(result, a, descending([b1, b2]), z);
        }
    }

    #[test]
    fn test_raised_recurrence_estimate() {
        // up in one, two and three upper parameters, in at most `m` steps each; the last one
        // from F = 1 at a_3 = 0
        for (a, [b1, b2], z, m) in [
            ([12.5, 0.75, 0.25], [2.5, 1.5], -0.5, 8),
            ([14.25, 9.5, 0.5], [0.75, -2.25], -0.3, 16),
            ([8.5, 6.25, 4.0], [3.5, 1.25], -0.9, 4),
        ] {
            let result = raised_recurrence(a, 0, [b1, b2], z, m, &Cell::new(RAISED_WORK));
            let result = result.and_then(|(y, err)| accurate(y, err));
            assert_estimate(result, a, descending([b1, b2]), z);
        }
    }

    #[test]
    fn test_raised_recurrence_work() {
        // the work that it takes is enough, and one less is not
        let run = |work| {
            let work = Cell::new(work);
            let result = raised_recurrence([12.5, 0.75, 0.25], 0, [2.5, 1.5], -0.5, 8, &work);
            (result.map(|(y, err)| (y.hi.to_bits(), err.to_bits())), work.get())
        };
        let (expected, left) = run(RAISED_WORK);
        assert!(expected.is_some());
        assert_eq!(run(RAISED_WORK - left), (expected, 0));
        assert_eq!(run(RAISED_WORK - left - 1).0, None);
    }

    #[test]
    fn test_recurrence_balanced_probe_budget() {
        // probe at `high - 8` fails with `PFQ_DD_MAXITER` and must leave enough of `BALANCED_WORK`
        // for the bisection midpoint `m = 18` to succeed
        let a = [16.983_936_871_431_336, -0.114_022_201_911_358_77, -36.194_302_717_639_12];
        let b = [0.871_800_762_335_526_7, 0.001_410_048_592_030_205_8];
        let (y, err) = recurrence_balanced(a, b, 0.996).unwrap();
        let expected = 194.547_589_493_083_65;
        assert!(
            err <= PFQ_TOL * y.abs() && (y - expected).abs() <= 1e-13 * expected.abs(),
            "{y} vs {expected} (err {err})"
        );
    }
}
