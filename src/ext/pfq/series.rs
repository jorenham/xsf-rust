//! Power series evaluation, ratio bounds, compensated summation, Weniger acceleration, and
//! overflow detection for `pFq`.

use core::iter::{once, zip};

use super::{
    Dd, Estimate, PFQ_DD_MAXITER, PFQ_LONG_MAXITER, PFQ_MAX_UNCHECKED, PFQ_MAXITER, PFQ_MODERATE,
    PFQ_TOL, Scratch, SeriesReal, dd_hi_down, dd_hi_up, first_zero, is_normal, ldexp_i64,
    next_up_nonneg, pfq_settled, sort_descending,
};
use crate::utils::pow2;
use crate::xsf::cephes::dd_real::{frexp_f64, ldexp_f64};

/// Most parameters, with the 1 of `n!`, for which the rounding errors of the methods in `f64`
/// are measured
const PFQ_MEASURED_PARAMS: f64 = 6.0;

/// Largest size of terms that add nothing to a sum, even if there are `PFQ_MAXITER` of them
const PFQ_NEGLIGIBLE: f64 = 1e-150;

/// Below the logarithm of the smallest positive normal number
const PFQ_LN_UNDERFLOW: f64 = -709.0;

/// The factor of the measured rounding errors for more than `PFQ_MEASURED_PARAMS` parameters,
/// as each of them adds its own to a term
fn pfq_rounding_scale(a: &[f64], b: &[f64]) -> f64 {
    ((a.len() + b.len() + 1) as f64 / PFQ_MEASURED_PARAMS).max(1.0)
}

/// Maximum order of `pfq_weniger`
const PFQ_WENIGER_MAXITER: usize = 40;

/// Number of the last steps of `pfq_weniger` in its estimate where its noise stopped it
/// (measured)
const PFQ_WENIGER_STEPS: usize = 4;
/// `b` with the `1` of `n!` inserted in descending order (for a descending `b`)
fn with_one(b: &[f64]) -> impl Iterator<Item = f64> + Clone + '_ {
    let k = b.partition_point(|&v| v > 1.0);
    b[..k]
        .iter()
        .copied()
        .chain(once(1.0))
        .chain(b[k..].iter().copied())
}

/// Bound on `|t_{m+1} / t_m|` for all `m >= n`, if all `a_i + n` and `b_j + n` are positive,
/// except that for a longer `a` the unpaired `a_i + m` are taken at `m = n` (and they grow)
///
/// `b` is the complete denominator (including the `1` of `n!`), of length `b_len`. Each `a_i` is
/// paired with `b_i`: `(a_i + m) / (b_i + m)` decreases to 1, or stays below 1. The other
/// `1 / (b_j + m)` decrease. This is tightest if `a` and `b` are sorted alike.
///
/// The factors are multiplied as mantissas and binary exponents so intermediate products cannot
/// overflow to `inf` and divide the bound down to `0`, and any inexact product or quotient is
/// rounded up (including subnormal conversion) so the result never falls below the exact bound
/// for the given `a`, `b`, `z`, and `n`.
fn pfq_ratio_bound_iter(
    a: &[f64],
    b: impl Iterator<Item = f64> + Clone,
    b_len: usize,
    z: f64,
    n: f64,
) -> f64 {
    let az = z.abs();
    let scale = |(m, e, r): (f64, i64, usize), f: f64, divide: bool| {
        let ok = f.is_finite() && if divide { f > 0.0 } else { f >= 0.0 };
        ok.then(|| {
            let (fm, fe) = frexp_f64(f);
            let (m, me) = frexp_f64(if divide { m / fm } else { m * fm });
            let e = if divide {
                e - i64::from(fe)
            } else {
                e + i64::from(fe)
            };
            (m, e + i64::from(me), r + 2)
        })
    };
    let acc = (az > 0.0 && az.is_finite()).then(|| {
        let (zm, ze) = frexp_f64(az);
        let init = (zm, i64::from(ze), 0_usize);
        let acc = zip(a, b.clone()).try_fold(init, |acc, (&ai, bj)| {
            if ai > bj {
                scale(scale(acc, ai + n, false)?, bj + n, true)
            } else {
                (ai + n >= 0.0 && bj + n > 0.0 && (bj + n).is_finite()).then_some(acc)
            }
        })?;
        let acc = a
            .iter()
            .skip(b_len)
            .try_fold(acc, |acc, &ai| scale(acc, ai + n, false))?;
        b.skip(a.len())
            .try_fold(acc, |acc, bj| scale(acc, bj + n, true))
    });
    if az == 0.0 {
        0.0
    } else if let Some((m, e, roundings)) = acc.flatten() {
        if roundings == 0 {
            az
        } else if m == 0.0 {
            0.0
        } else {
            let slack = 1.0 + (roundings as f64 + 2.0) * f64::EPSILON;
            let bound = ldexp_i64(m * slack, e);
            if !bound.is_finite() {
                f64::INFINITY
            } else if bound < f64::MIN_POSITIVE {
                next_up_nonneg(bound)
            } else {
                bound
            }
        }
    } else {
        f64::INFINITY
    }
}

/// `pfq_ratio_bound` for a complete descending denominator `b` that already includes the `n + 1`
/// factorial factor (e.g. from `Residues` or already-shifted `pfq_alg_series_dd`)
fn pfq_ratio_bound_full(a: &[f64], b: &[f64], z: f64, n: f64) -> f64 {
    pfq_ratio_bound_iter(a, b.iter().copied(), b.len(), z, n)
}

/// Bound on `|t_{m+1} / t_m|` for all `m >= n` (if all `a_i + n` and `b_j + n` are positive,
/// except that for `a.len() > b.len() + 1` the unpaired `a_i + m` are taken at `m = n` and grow),
/// for descending `a` and `b` without the `1` of `n!`, which is inserted in descending order
fn pfq_ratio_bound(a: &[f64], b: &[f64], z: f64, n: f64) -> f64 {
    pfq_ratio_bound_iter(a, with_one(b), b.len() + 1, z, n)
}

/// Certifies that for `up_s.len() > low_s.len()`, term ratios stay at least 1 once they reach it
///
/// The ratios are at most 1 up to `m0 = n + delta`, where `base` bounds `|z|` times the paired
/// ratios and `delta = base^(-1/k) - x_s` for the largest unpaired `x_s = up_s[low_s.len()]`. From
/// `m0` on they grow if the derivative of their log is positive: the `k` unpaired factors add at
/// least `k / (x_s + m - n)`, and each pair with `u > l` subtracts
/// `(u - l) / ((u + m - n) (l + m - n))`. Times `x_s + m - n`, that is at most
/// `(u - l) / (l + delta)` times `max(1, (x_s + delta) / (u + delta))` (for descending `up_s`,
/// `(x_s + m - n) / (u + m - n)` is at most 1 and increases toward 1).
fn regrowth_certified(
    up_s: &[f64],
    low_s: &[f64],
    unpaired: &[f64],
    z: f64,
    bound: f64,
) -> Option<()> {
    bound.is_finite().then_some(())?;
    let (k, x_s) = (unpaired.len() as f64, unpaired[0]);
    // conservative upper bound on `base`, giving a conservative lower bound on `delta`
    let base = pfq_ratio_bound_full(&up_s[..low_s.len()], low_s, z, 0.0);
    base.is_normal().then_some(())?;
    // certify `base * root^k <= 1` via an upward-rounded scaled product
    let fits = |r: f64| {
        (r > 0.0 && r.is_finite()).then(|| {
            let ((bm, be), (rm, re)) = (frexp_f64(base), frexp_f64(r));
            let (m, e) = (0..unpaired.len()).fold((bm, i64::from(be)), |(m, e), _| {
                let (nm, ne) = frexp_f64(m * rm);
                (nm, e + i64::from(re) + i64::from(ne))
            });
            let slack = 1.0 + (2 * unpaired.len() + 2) as f64 * f64::EPSILON;
            ldexp_i64(m * slack, e) <= 1.0
        }) == Some(true)
    };
    let r0 = base.powf(-k.recip()) * (1.0 - 128.0 * f64::EPSILON);
    let root = [r0, r0 * (1.0 - 1024.0 * f64::EPSILON)]
        .into_iter()
        .find(|&r| fits(r))?;
    let delta = root - x_s;
    let delta = if delta > 0.0 {
        f64::from_bits(delta.to_bits() - 1)
    } else {
        delta
    };
    let pulls = zip(up_s, low_s).filter(|&(u, l)| u > l);
    let pull = pulls.map(|(&u, &l)| (u - l) / (l + delta) * ((x_s + delta) / (u + delta)).max(1.0));
    let pull_slack = 1.0 + (4.0 * low_s.len() as f64 + 4.0) * f64::EPSILON;
    (delta >= 0.0 && pull.sum::<f64>() * pull_slack <= k).then_some(())
}

/// Reusable scratch buffers for shifted `Dd` ratio/regrowth bounds
pub(super) struct BoundScratch {
    up: Scratch<f64>,
    low: Scratch<f64>,
}

impl BoundScratch {
    pub(super) fn new(up_len: usize, low_len: usize) -> Self {
        Self {
            up: Scratch::new(up_len, 0.0),
            low: Scratch::new(low_len, 0.0),
        }
    }

    fn shift<'a>(&'a mut self, up: &[Dd], low: &[Dd], n: f64) -> (&'a [f64], &'a [f64]) {
        let up_s = self.up.fill(up.len(), up.iter().map(|&x| dd_hi_up(x + n)));
        let low_s = self
            .low
            .fill(low.len(), low.iter().map(|&x| dd_hi_down(x + n)));
        (up_s, low_s)
    }

    /// `pfq_ratio_bound_full` on `up + n` (rounded up) and `low + n` (rounded down) in `Dd`
    pub(super) fn ratio_bound(&mut self, up: &[Dd], low: &[Dd], z: f64, n: f64) -> f64 {
        let (up_s, low_s) = self.shift(up, low, n);
        pfq_ratio_bound_full(up_s, low_s, z, 0.0)
    }

    /// `ratio_bound`, or infinite unless ratios for `up.len() > low.len()` stay at least 1 once
    /// they reach it
    pub(super) fn regrowth_bound(&mut self, up: &[Dd], low: &[Dd], z: f64, n: f64) -> f64 {
        let (up_s, low_s) = self.shift(up, low, n);
        let bound = pfq_ratio_bound_full(up_s, low_s, z, 0.0);
        match up_s.get(low_s.len()..).filter(|v| !v.is_empty()) {
            None => bound,
            Some(unpaired) => regrowth_certified(up_s, low_s, unpaired, z, bound)
                .map_or(f64::INFINITY, |()| bound),
        }
    }
}

/// The products of the `a_i + n` and of `(1 + n)` times the `b_j + n`, whose quotient is
/// `t_{n+1} / (z t_n)`
///
/// They and their quotient are normal for up to `PFQ_MAX_UNCHECKED` factors with sizes in
/// `PFQ_MODERATE`, as each nonzero `|a_i + n|` is then in `[2^-100, 2^101]`.
fn pfq_products<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> [T; 2] {
    let shifted = |p: T, &v: &f64| p * (T::from_f64(v) + n);
    let one = T::from_f64(1.0);
    [a.iter().fold(one, shifted), b.iter().fold(one + n, shifted)]
}

/// `t_{n+1} / (z t_n)` for any parameters and `n`, as a mantissa and an `i64` exponent, so that
/// no partial product over- or underflows
pub(super) fn pfq_frexp_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> (T, i64) {
    let product = |params: &[f64], init: (T, i64)| {
        let times = |(m, e): (T, i64), &v: &f64| {
            let (vm, ve) = (T::from_f64(v) + n).frexp();
            let (m, me) = (m * vm).frexp();
            (m, e + i64::from(ve) + i64::from(me))
        };
        params.iter().fold(init, times)
    };
    let (fm, fe) = (T::from_f64(1.0) + n).frexp();
    let ((num, ne), (den, de)) =
        (product(a, (T::from_f64(1.0), 0)), product(b, (fm, i64::from(fe))));
    (num / den, ne - de)
}

/// `t_{n+1} / (z t_n)` for any parameters and `n`, which is zero if it underflows, and infinite
/// if it overflows or would be inexact, as it is not normal
#[cold]
fn pfq_checked_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64) -> T {
    let (ratio, e) = pfq_frexp_ratio::<T>(a, b, n);
    let ratio = ldexp_i64(ratio, e);
    let ln_abs = |params: &[f64]| params.iter().map(|v| (v + n).abs().ln()).sum::<f64>();
    if is_normal(ratio) {
        ratio
    } else if ln_abs(a) - ln_abs(b) - (1.0 + n).ln() < PFQ_LN_UNDERFLOW {
        T::from_f64(0.0)
    } else {
        T::from_f64(f64::INFINITY)
    }
}

/// Bound on the size of the term after `term` if it underflows, or if its ratio does
#[cold]
#[inline(never)]
fn pfq_underflow_bound<T: SeriesReal>(term: T, z: f64) -> T {
    T::from_f64(T::MIN_POSITIVE * term.to_f64().abs().max(1.0) * z.abs().max(1.0))
}

/// Whether `pfq_products` is not exact enough for these parameters
pub(super) fn pfq_extreme(a: &[f64], b: &[f64]) -> bool {
    let many = a.len().max(b.len() + 1) > PFQ_MAX_UNCHECKED;
    many || a.iter().chain(b).any(|v| !PFQ_MODERATE.contains(&v.abs()))
}

/// `t_{n+1} / (z t_n)`, which is `checked` if the parameters are extreme
#[inline]
pub(super) fn pfq_ratio<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, checked: bool) -> T {
    if checked {
        pfq_checked_ratio(a, b, n)
    } else {
        let [num, den]: [T; 2] = pfq_products(a, b, n);
        num / den
    }
}

/// `t_{n+1} / t_n` as `Ok(ratio)` if the products, their quotient, and the product with `z` are
/// normal, else `Err(ratio)` with the already-computed `q * z`
///
/// With `checked = false`, the parameters must not be `pfq_extreme` so that partial products of
/// `pfq_products` stay normal.
#[inline]
fn pfq_normal_ratio<T: SeriesReal>(
    a: &[f64],
    b: &[f64],
    n: f64,
    z: f64,
    checked: bool,
) -> Result<T, T> {
    let (q, products_normal) = if checked {
        (pfq_checked_ratio(a, b, n), true)
    } else {
        let [num, den]: [T; 2] = pfq_products(a, b, n);
        (num / den, is_normal(num) && is_normal(den))
    };
    let ratio = q * T::from_f64(z);
    if products_normal && is_normal(q) && is_normal(ratio) {
        Ok(ratio)
    } else {
        Err(ratio)
    }
}

/// The term `t_{n+1}` after `term = t_n`, with a ratio that is `checked` if the parameters are
/// extreme
#[inline]
fn pfq_next<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, z: f64, term: T, checked: bool) -> T {
    match pfq_normal_ratio::<T>(a, b, n, z, checked) {
        Ok(ratio) => term * ratio,
        Err(ratio) => pfq_rescaled_next(a, b, n, z, term).unwrap_or_else(|| term * ratio),
    }
}

/// `pfq_next` without its fallback to the rounded ratio, which can have more rounding errors, so
/// `None` where `t_{n+1}` is not normal
fn pfq_normal_next(a: &[f64], b: &[f64], n: f64, z: f64, term: f64, checked: bool) -> Option<f64> {
    match pfq_normal_ratio::<f64>(a, b, n, z, checked) {
        Ok(ratio) => Some(term * ratio).filter(|next| next.is_normal()),
        Err(_) => pfq_rescaled_next(a, b, n, z, term),
    }
}

/// `pfq_next` from the mantissas and exponents of its factors, as an intermediate product or
/// quotient can over- or underflow where `t_{n+1}` is normal; `None` if `t_{n+1}` is not normal
/// either
#[cold]
#[inline(never)]
fn pfq_rescaled_next<T: SeriesReal>(a: &[f64], b: &[f64], n: f64, z: f64, term: T) -> Option<T> {
    let (ratio, e) = pfq_frexp_ratio::<T>(a, b, n);
    let ((zm, ze), (tm, te)) = (frexp_f64(z), term.frexp());
    let rescaled = ldexp_i64(tm * (ratio * T::from_f64(zm)), e + i64::from(ze + te));
    is_normal(rescaled).then_some(rescaled)
}

/// The terms `t_1, t_2, ...` of the series after `t_0 = 1`, in the arithmetic of `T`
///
/// A term that underflows is a bound on it. The terms end before one that overflows or is NaN,
/// or that is such a bound and not negligible.
pub(super) fn pfq_terms<T: SeriesReal>(a: &[f64], b: &[f64], z: f64) -> impl Iterator<Item = T> {
    let checked = pfq_extreme(a, b);
    pfq_terms_of(first_zero(a), z, move |n, term| pfq_next(a, b, n, z, term, checked))
}

/// `pfq_terms` with `next(n, t_n) = t_{n+1}`, for the degree of a polynomial, else infinity
pub(super) fn pfq_terms_of<T: SeriesReal>(
    degree: f64,
    z: f64,
    next: impl Fn(f64, T) -> T,
) -> impl Iterator<Item = T> {
    // the largest term, which is small once the terms are only upper bounds
    let mut limit = f64::MAX;
    let indices = (0_u32..).map(f64::from);
    indices.scan(T::from_f64(1.0), move |term, n| {
        let next = next(n, *term);
        // (the terms after the last one of a polynomial are zero; tested last, as that is faster)
        *term = if next.to_f64().abs() < T::MIN_POSITIVE && n < degree {
            // the term or its ratio underflows, after which the terms can become large again:
            // follow a bound on them
            limit = PFQ_NEGLIGIBLE;
            pfq_underflow_bound(*term, z)
        } else {
            next
        };
        (term.to_f64().abs() <= limit).then_some(*term)
    })
}

/// The series `sum_n t_n` (DLMF 16.2.1) and an estimate of its relative error
///
/// The estimate is infinite if the series did not converge within `PFQ_MAXITER` terms, or if
/// its terms overflow, or underflow before they become large again. If it converged, but its
/// terms cancel, then it is summed again in double-double arithmetic.
pub(super) fn pfq_series(a: &[f64], b: &[f64], z: f64) -> Estimate {
    let (y, err, _) = pfq_sum::<f64>(a, b, z, PFQ_MAXITER);
    // (not in double-double, whose `EPSILON` leaves room for many more parameters)
    let err = err * pfq_rounding_scale(a, b);
    let y = if err <= PFQ_TOL * y.abs() || err == f64::INFINITY {
        (y, err)
    } else {
        let (y, err, _) = pfq_sum::<Dd>(a, b, z, PFQ_DD_MAXITER);
        (f64::from(y), err)
    };
    Estimate::absolute(y)
}

/// The long series in double-double, if it is accurate
pub(super) fn pfq_series_long(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    // A 3F2 close to 1 can need more than 20,000 terms to bound the remainder.
    let maxiter = if a.len() == 3 && z >= 0.99 {
        PFQ_LONG_MAXITER + PFQ_MAXITER
    } else {
        PFQ_LONG_MAXITER
    };
    let (y, err, _) = pfq_sum::<Dd>(a, b, z, maxiter);
    let y = Estimate::absolute((f64::from(y), err));
    y.is_accurate().then_some(y)
}

/// Largest partial sum, in units of `unit`, whose square is added to the sums of squares
const PFQ_SQUARE_MAX: f64 = 1e150;

/// The factor 2^-400 by which the units of those sums shrink once it is exceeded
const PFQ_UNIT: f64 = pow2(-400);

/// Maximum number of roundings of a ratio `t_{k+1} / t_k`: each `a_i + k`, `b_j + k`, and
/// `1 + k`, the product, the quotient, and the factors `z` and `t_k`
fn pfq_ratio_roundings(a: &[f64], b: &[f64]) -> f64 {
    2.0 * (a.len() + b.len() + 1) as f64 + 3.0
}

/// `pfq_series` in the arithmetic of `T`, with at most `maxiter` terms, and with the sum
/// and the root of its squared tails in that arithmetic
///
/// The rounding error of each ratio `t_{n+1} / t_n` affects the terms after it, which sum to
/// `s - s_n` for the partial sum `s_n` and the sum `s`. These errors are added in quadrature,
/// like those of the partial sums. A part of them adds up instead, as `a + n` rounds alike
/// for all `n` between two powers of 2. The constants are measured.
pub(super) fn pfq_sum<T: SeriesReal>(a: &[f64], b: &[f64], z: f64, maxiter: u32) -> (T, f64, f64) {
    let terms = pfq_terms::<T>(a, b, z);
    pfq_sum_terms(a, b, z, maxiter, terms, f64::EPSILON)
}

/// `pfq_sum` of the given terms `t_1, t_2, ...`, up to a bound on the rest of them of `tol` times
/// the sum. The terms can be those of a more exact `b`.
#[inline]
pub(super) fn pfq_sum_terms<T: SeriesReal>(
    a: &[f64],
    b: &[f64],
    z: f64,
    maxiter: u32,
    mut terms: impl Iterator<Item = T>,
    tol: f64,
) -> (T, f64, f64) {
    let settled = pfq_settled(a, b);
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    let (mut sum, mut comp) = (T::from_f64(1.0), T::from_f64(0.0));
    // the compensated partial sum `s_n` and its offset `d_n = s_n - 1` as `f64`
    let (mut s, mut d) = T::sum_and_delta(sum, comp);
    // the sums of `d_n`, `d_n^2`, and `s_n^2`, in units of `unit` so that squares don't overflow;
    // shifting by `t_0 = 1` keeps `s - s_n = d - d_n` accurate near a flat sum
    let (mut dsum1, mut dsum2, mut sum2, mut unit) = (0.0, 0.0_f64, 0.0_f64, 1.0);
    // the sum of |t_k|, in units of `unit`
    let mut abs1 = 1.0;
    let ratio_roundings = pfq_ratio_roundings(a, b);
    // the rounding error, and the root of the sum of (s - s_n)^2 = (d - d_n)^2
    let estimates = |n: f64, d: f64, dsum1: f64, dsum2: f64, sum2: f64, unit: f64, abs1: f64| {
        // scaled by max(|d * unit|, sqrt(dsum2)) so that neither tiny `d` nor huge sums overflow
        let root2 = sum2.sqrt();
        let (du, droot2) = (d * unit, dsum2.sqrt());
        let scale = du.abs().max(droot2);
        let tails = if scale == 0.0 {
            0.0
        } else {
            let (r, r2) = (du / scale, droot2 / scale);
            ((n + 1.0) * (r * r) - 2.0 * r * (dsum1 / scale) + r2 * r2)
                .max(0.0)
                .sqrt()
                * scale
        };
        let mut rounding = T::EPSILON * ((4.0 + 0.25 * (n + 1.0).sqrt()) * tails + root2);
        if sum2 > f64::MAX * unit * unit {
            // The measured estimate was only accepted where the squares don't overflow. Beyond,
            // t_k has at most k such roundings that all go the same way, so the sum has at most
            // sum_k k |t_k| <= (n + 1) sum_k |t_k| of them.
            rounding += T::EPSILON * (n + 1.0) * ratio_roundings * abs1;
        }
        (rounding / unit, tails / unit)
    };
    for n in (0..maxiter).map(f64::from) {
        while s.abs().max(d.abs()) * unit > PFQ_SQUARE_MAX {
            // the squares would overflow, so the sums are rescaled (exactly, by a power of 2)
            let u2 = PFQ_UNIT * PFQ_UNIT;
            (dsum1, dsum2, sum2, unit) = (dsum1 * PFQ_UNIT, dsum2 * u2, sum2 * u2, unit * PFQ_UNIT);
            abs1 *= PFQ_UNIT;
        }
        let (su, du) = (s * unit, d * unit);
        dsum1 += du;
        dsum2 += du * du;
        sum2 += su * su;
        if n == degree {
            // the series terminated
            let (err, tails) = estimates(n, d, dsum1, dsum2, sum2, unit, abs1);
            return (sum + comp, err, tails);
        }
        let Some(term) = terms.next() else { break };
        let t = term.to_f64();
        abs1 += t.abs() * unit;
        (sum, comp) = T::add_compensated(sum, comp, term);
        (s, d) = T::sum_and_delta(sum, comp);
        if !s.is_finite() {
            break;
        }
        if n >= settled && t.abs() <= tol * s.abs() {
            let rho = pfq_ratio_bound(a, b, z, n + 1.0);
            // geometric bound on the remaining terms, if the series converges (p <= q + 1)
            let tail = t.abs() * rho / (1.0 - rho);
            if a.len() <= b.len() + 1 && rho < 1.0 && tail <= tol * s.abs() {
                let (rounding, tails) = estimates(n, d, dsum1, dsum2, sum2, unit, abs1);
                return (sum + comp, rounding + tail, tails);
            }
        }
    }
    (sum + comp, f64::INFINITY, f64::INFINITY)
}

/// Weniger's delta transformation of the series, and an estimate of its absolute error
///
/// This is a factorial Levin-type sequence transformation \[1\], which is used for
/// pFq in HypergeometricFunctions.jl \[2\]. For `z < 0` it converges much faster than the
/// series, and also for `z < -1`, where the series diverges. With the term `t_{n+1}` as the
/// remainder estimate of the partial sum `s_n`, its order `k` is
///
/// ```text
/// delta_k = (sum_{j=0}^k c_j s_j / t_{j+1}) / (sum_{j=0}^k c_j / t_{j+1}),
/// c_j = (-1)^j binom(k, j) (j + 1)_{k-1},
/// ```
///
/// by the three-term recursion of its numerator and denominator. It is applied from the `n`
/// where the terms start to alternate, so that the weights of the `s_j` are positive. The same
/// recursion for `|s_j|` then gives their weighted mean, which scales the rounding errors.
///
/// It is computed in the arithmetic of `T`. The estimate is infinite if the orders up to
/// `PFQ_WENIGER_MAXITER` did not converge, or if a term overflows, or underflows before the terms
/// become large again. It is `None` if the terms only alternate after more than `PFQ_MAXITER`
/// terms, or if a term before that overflows or underflows.
///
/// 1. E. J. Weniger, "Nonlinear sequence transformations for the acceleration of convergence
///    and the summation of divergent series", Comput. Phys. Rep. 10 (1989), 189-371
/// 2. R. M. Slevinsky, "Fast and stable rational approximation of generalized hypergeometric
///    functions", Numer. Algor. 98 (2025), 587-624, <https://doi.org/10.1007/s11075-024-01808-w>
pub(super) fn pfq_weniger<T: SeriesReal>(a: &[f64], b: &[f64], z: f64) -> Option<(f64, f64)> {
    let settled = pfq_settled(a, b).ceil();
    if settled > f64::from(PFQ_MAXITER) {
        return None;
    }
    // the degree of a polynomial, else infinity
    let degree = first_zero(a);
    // the transformation starts at this n
    let start = settled.min(degree);
    let checked = pfq_extreme(a, b);
    // the rounding errors of the partial sums and of the terms (the factor is measured)
    let eps = 4.0 * T::EPSILON * pfq_rounding_scale(a, b);
    let one = T::from_f64(1.0);
    let (mut term, mut sum) = (one, one);
    // the largest |s_n| before `start`
    let mut scale = 1.0_f64;
    for n in (0..start as u32).map(f64::from) {
        term = pfq_next(a, b, n, z, term, checked);
        if !is_normal(term) {
            // overflow or underflow
            return None;
        }
        sum = sum + term;
        scale = scale.max(sum.to_f64().abs());
    }
    if start == degree {
        // the series terminated
        return Some((sum.to_f64(), eps * scale));
    }
    let mut rows = [[T::from_f64(0.0); 3]; PFQ_WENIGER_MAXITER];
    let mut delta = sum;
    // the last steps, the newest first
    let mut steps = [f64::INFINITY; PFQ_WENIGER_STEPS];
    for k in 0..PFQ_WENIGER_MAXITER {
        let n = start + k as f64;
        term = pfq_next(a, b, n, z, term, checked);
        if !is_normal(term) {
            // overflow, or underflow, after which the terms only matter if they can increase
            let negligible = term.to_f64().abs() < 1.0
                && a.len() <= b.len() + 1
                && pfq_ratio_bound(a, b, z, n) <= 1.0;
            let noise = eps * scale;
            return Some((sum.to_f64(), if negligible { noise } else { f64::INFINITY }));
        }
        // the numerator, the denominator, and the numerator for |s_j|
        rows[k] = [sum / term, one / term, sum.abs() / term];
        for j in (0..k).rev() {
            // the order of `rows[j]`, which starts at s_j
            let (m, i) = (j as f64, (k - j - 1) as f64);
            let c = if i == 0.0 {
                one
            } else {
                let (s, t) = (m + i, m + 2.0 * i);
                T::from_f64((s + 1.0) * s) / T::from_f64((t + 1.0) * t)
            };
            rows[j] = core::array::from_fn(|l| rows[j + 1][l] - c * rows[j][l]);
        }
        let [num, den, abs] = rows[0];
        let next = num / den;
        let noise = eps * ((abs / den).to_f64().abs() + next.to_f64().abs() + scale);
        // the steps need not be below the rounding of the `f64` result, which only double-double
        // gets under (a NaN noise stays NaN)
        let rounding = f64::EPSILON * next.to_f64().abs();
        let tol = if noise < rounding { rounding } else { noise };
        steps.copy_within(..PFQ_WENIGER_STEPS - 1, 1);
        steps[0] = (next - delta).to_f64().abs();
        delta = next;
        // delta is a weighted mean of the last one and of the one of that order from s_1
        // (NaN for k = 0). For a large |z| the weight of the latter is small, and so are the
        // steps, also while the two still differ, so it has to agree as well.
        let other = (rows[1][0] / rows[1][1] - delta).to_f64().abs();
        if steps[0] <= tol && steps[1] <= tol && other <= tol {
            // The steps also count, as they can be smaller than the error once the noise grows
            // past them (e.g. by 4x per order for z = -16). Then they can also dip far below the
            // error where it changes sign, so the largest of the last few (that it took) counts.
            let counted = if noise < rounding { 2 } else { steps.len() }.min(k + 1);
            let err = noise + steps[..counted].iter().fold(other, |m, &s| m.max(s));
            return Some((delta.to_f64(), err));
        }
        sum = sum + term;
    }
    Some((delta.to_f64(), f64::INFINITY))
}

/// Double-double Weniger estimate including the rounding errors accumulated by the terms
/// up to `settled`, as `a + n` can round alike.
pub(super) fn pfq_weniger_dd(a: &[f64], b: &[f64], z: f64, settled: f64) -> Option<Estimate> {
    let (y, err) = pfq_weniger::<Dd>(a, b, z)?;
    // (each scaled, as their sum can overflow)
    let drift = once(1.0)
        .chain(pfq_terms(a, b, z))
        .zip(0..=settled as u32)
        .map(|(t, n)| 4.0 * Dd::EPSILON * f64::from(n) * t.abs());
    Some(Estimate::absolute((y, err + drift.sum::<f64>())))
}

/// Most terms that `pfq_overflow` follows
///
/// Of its 57 results without it on 1657 points where the other methods are not accurate, 55 are
/// within it, and 54 within half of it (measured).
const PFQ_OVERFLOW_MAXITER: u32 = 4 * PFQ_MAXITER;

/// `+-inf` for `0 < z < 1` where the series shows that `_{q+1}F_q` overflows, else `None`
///
/// From the first `n >= pfq_settled` on, the terms have the sign of `t_n`. Once the sum of their
/// sizes is at least `2^1025`, and that of the terms before at most `2^1023`, the result has that
/// sign and a size of at least `2^1024`. Each term has at most `pfq_ratio_roundings` roundings per
/// step before it, and each sum one per term.
pub(super) fn pfq_overflow(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    if z <= 0.0 {
        return None;
    }
    let (settled, checked) = (pfq_settled(a, b), pfq_extreme(a, b));
    let roundings = pfq_ratio_roundings(a, b) + 1.0;
    let next = |n, term| pfq_normal_next(a, b, n, z, term, checked);
    let indices = (0..PFQ_OVERFLOW_MAXITER).map(f64::from);

    // the sum of the sizes of the terms before `settled`, and the first term after them
    let (before, mut t) = indices
        .clone()
        .take_while(|&n| n < settled)
        .try_fold((0.0, 1.0), |(sum, t): (f64, f64), n| Some((sum + t.abs(), next(n, t)?)))?;

    // the sum of the sizes of the terms from `settled` on, and `t`, both in units of `2^e`
    let mut after = 0.0;
    let mut e = 0;
    for n in indices.skip_while(|&n| n < settled) {
        // the sum in [1/2, 1), exactly, and t with it, which is exact unless it underflows
        let (m, k) = frexp_f64(after + t.abs());
        after = m;
        e += k;
        t = ldexp_f64(t, -k);
        if !t.is_normal() {
            return None;
        }
        // (a bound of the relative rounding errors, with a factor 2 to spare)
        let err = (n + 1.0) * roundings * f64::EPSILON;
        // (the lower bound `>= 2^1025` by its exponent, as it can be beyond `f64`)
        if before * (1.0 + err) <= pow2(1023) && frexp_f64(after * (1.0 - err)).1 + e > 1025 {
            return Some(Estimate::exact(f64::INFINITY.copysign(t)));
        }
        let next = next(n, t)?;
        if next.abs() < t.abs() {
            // the rest is at most a geometric series
            let rho = pfq_ratio_bound(a, b, z, n + 1.0);
            if rho < 1.0 && before + ldexp_f64(after + next.abs() / (1.0 - rho), e) < pow2(1023) {
                return None;
            }
        }
        t = next;
    }
    None
}

/// `pfq_series` after sorting the parameters, which tightens `pfq_ratio_bound`
pub(super) fn pfq_series_sorted(a: &[f64], b: &[f64], z: f64) -> Estimate {
    let (mut a_s, mut b_s) = (Scratch::new(a.len(), 0.0), Scratch::new(b.len(), 0.0));
    let (a, b) = (a_s.fill(a.len(), a.iter().copied()), b_s.fill(b.len(), b.iter().copied()));
    sort_descending(a);
    sort_descending(b);
    pfq_series(a, b, z)
}

#[cfg(test)]
mod tests {
    use super::super::recip_up;
    use super::*;

    #[test]
    fn test_pfq_sum_rescales_repeatedly() {
        // the partial sum jumps from 1 to 1e300, whose square needs two rescalings
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-1.0], &[], 1.0, 2, [1e300].into_iter(), 1e-13);
        assert_eq!(y, 1e300);
        assert!(err.is_finite() && tails.is_finite());
    }

    #[test]
    fn test_pfq_ratio_bound_safe() {
        // huge denominator factors with huge `z`: exact ratio is 2^-7, not 0 from `inf` division
        let rho = pfq_ratio_bound(&[], &[pow2(520), pow2(510)], pow2(1023), 0.0);
        assert!(rho >= pow2(-7) && rho < 1.0, "{rho}");
        let rho = pfq_ratio_bound_full(&[], &[pow2(520), pow2(510)], pow2(1023), 0.0);
        assert!(rho >= pow2(-7) && rho < 1.0, "{rho}");
        // full bound at n > 0 does not insert an extra (1 + n) factor (exact 1/3, not 1/6)
        let rho = pfq_ratio_bound_full(&[], &[2.0], 1.0, 1.0);
        assert!((1.0 / 3.0..0.34).contains(&rho), "{rho}");
        // exact ratio exceeds 1 by ~4.4e-17 while the unpadded f64 product rounds below 1
        let b = 5.0 * pow2(78);
        let a = b + 3.0 * pow2(28);
        let z = -(1.0 + 7.0 * f64::EPSILON);
        let rho = pfq_ratio_bound(&[a, a, a, 1.0], &[b, b, b], 1.0 / z, 9999.0);
        assert!(rho >= 1.0, "{rho}");
        // exact bound of 1 when every a_i <= b_i (including implicit 1) at |z| = 1
        assert_eq!(pfq_ratio_bound(&[1.0, 1.0, 1.0], &[1e300, 2.0], -1.0, 1.0), 1.0);
        // implicit 1 is inserted in descending order: [0.75] -> [1.0, 0.75] pairs with [1.5, 0.5]
        let rho = pfq_ratio_bound(&[1.5, 0.5], &[0.75], -1.0, 1.0);
        assert!((1.25..1.3).contains(&rho), "{rho}");
        // subnormal scaling rounds upward, never down or to zero
        let u = f64::from_bits(1);
        assert!(pfq_ratio_bound(&[], &[4.0], 5.0 * u, 0.0) >= 2.0 * u);
        assert!(pfq_ratio_bound(&[], &[4.0], u, 0.0) >= u);
        // reciprocal bound rounds upward even on the exact-z fast path
        assert!(pfq_ratio_bound(&[1.0], &[], recip_up(-3.0), 0.0) > 1.0 / 3.0);
        assert!(pfq_ratio_bound_full(&[1.0], &[2.0], recip_up(-3.0), 0.0) > 1.0 / 3.0);
        // before settled (a_i + n < 0), the bound is infinite
        assert_eq!(pfq_ratio_bound(&[-2.0], &[], 0.5, 0.0), f64::INFINITY);
    }

    #[test]
    fn test_pfq_sum_flat_compensated() {
        // 1000 terms of 1e-17 each: preserved by Neumaier compensation, with nonzero tail norm
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-1000.0], &[], 1.0, 1001, core::iter::repeat(1e-17), 1e-13);
        assert!((y - (1.0 + 1e-14)).abs() <= 2.0 * f64::EPSILON, "{y}");
        assert!(tails > 1e-13 && err.is_finite(), "{tails}, {err}");
        // cancelling partial sums leaving a tiny final offset: tail norm stays finite near 1
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-3.0], &[], 1.0, 4, [1.0, -1.0, 1e-200].into_iter(), 1e-13);
        assert_eq!(y, 1.0);
        assert!(err.is_finite() && (tails - 1.0).abs() <= 4.0 * f64::EPSILON, "{err}, {tails}");
        let (y, err, tails) = pfq_sum_terms::<Dd>(
            &[-3.0],
            &[],
            1.0,
            4,
            [1.0, -1.0, 1e-200].into_iter().map(Dd::new),
            1e-13,
        );
        assert_eq!(y.hi, 1.0);
        assert!(err.is_finite() && (tails - 1.0).abs() <= 4.0 * f64::EPSILON, "{err}, {tails}");
        // non-finite partial sum bails out cleanly with infinite error
        let (y, err, tails) =
            pfq_sum_terms::<f64>(&[-2.0], &[], 1.0, 3, [f64::MAX, f64::MAX].into_iter(), 1e-13);
        assert!(y.is_infinite() && err == f64::INFINITY && tails == f64::INFINITY);
    }

    #[test]
    fn test_pfq_next_rescaled_subnormal_quotient() {
        // extreme parameters whose ratio underflows while `ratio * z` is normal (checked mode)
        let next = pfq_next::<f64>(&[1e-160, 1e-160], &[], 0.0, 1e200, 1.0, true);
        assert!((next / 1e-120 - 1.0).abs() <= 4.0 * f64::EPSILON, "{next}");
        // smallest subnormal term is not double-rounded to zero by the rescaled path
        let a = [f64::from_bits(0x0010_0000_0000_0001)];
        let z = f64::from_bits(0x3c9f_ffff_ffff_ffff);
        assert_eq!(pfq_next::<f64>(&a, &[], 0.0, z, 1.0, true), f64::from_bits(1));
    }

    #[test]
    fn test_pfq_series_sorted_boundary() {
        // unsorted `(p, q)` across the 8-element stack/heap boundary of `Scratch<f64>`
        for (p, q) in [(8, 8), (8, 9), (9, 9)] {
            let a: Vec<f64> = (0..p).map(|i| 0.5 + 0.75 * f64::from(i)).collect();
            let b: Vec<f64> = (0..q).map(|j| 1.25 + 0.5 * f64::from(j)).collect();
            let mut a_sorted = a.clone();
            let mut b_sorted = b.clone();
            sort_descending(&mut a_sorted);
            sort_descending(&mut b_sorted);
            let actual = pfq_series_sorted(&a, &b, -2.5);
            let expected = pfq_series(&a_sorted, &b_sorted, -2.5);
            assert!(actual.is_accurate(), "({p}, {q}): {actual:?}");
            assert_eq!(actual.y.to_bits(), expected.y.to_bits(), "({p}, {q})");
            assert_eq!(actual.err.to_bits(), expected.err.to_bits(), "({p}, {q})");
        }
    }
}
