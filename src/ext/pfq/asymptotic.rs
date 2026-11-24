//! Asymptotic, algebraic, exponential, and analytic-continuation expansions for `pFq`.

use core::f64::consts::{PI, TAU};
use core::iter::{from_fn, once, successors, zip};
use core::ops::ControlFlow;

use super::gamma::{LogGammaCache, gamma_alg_term, gamma_series};
use super::residue::{self, Residues};
use super::series::{BoundScratch, pfq_terms_of};
use super::{
    Dd, Estimate, OrBetter, PFQ_MAX_UNCHECKED, PFQ_MAXITER, PFQ_MODERATE, Scaled, Scratch,
    SeriesReal, add_exp, add_exp_dd, dd_hi_up, first_zero, first_zero_dd, is_normal, is_pole,
    ldexp_i64, ln_dd, next_up_nonneg, positive_from, recip_up, root_dd, skip, sort_descending_dd,
    sum_dd,
};
use crate::utils::{EPSILON_SQ, is_nonpos_int, neg_one_pow};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::{FRAC_2_PI_DD, LN_SQRT_2PI_DD, frexp_f64};
use crate::xsf::{cospi, sinpi};

/// The analytic continuation to `z < -1` (DLMF 16.8.8) and an estimate of its relative error
/// (which may exceed `PFQ_TOL`), or `None` if evaluation fails or the result is not normal
///
/// This is for `p = q + 1`. It evaluates the algebraic sum via `pfq_alg` (`Alg::Value`), falling
/// back to `pfq_alg_dd` on failure or inaccuracy when at least one `a_i` is not in an
/// integer-difference class. It can still return `Some` when neither method meets `PFQ_TOL`, so
/// callers must check the returned error bound.
pub(super) fn pfq_continuation(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let to_est = |s: Scaled| {
        let (y, err) = s.to_f64();
        y.is_normal().then(|| Estimate::new(y, err))
    };
    let best = pfq_alg(a, b, z, Alg::Value)
        .and_then(add_exp)
        .and_then(to_est);
    let has_singleton =
        zip(a, 0..).any(|(&aj, j)| skip(a, j).all(|ak| !(Dd::new(ak) - aj).is_integer()));
    if has_singleton {
        best.or_better(|| pfq_alg_dd(a, b, z).and_then(to_est))
    } else {
        best
    }
}

/// A term of a series with a bound on its size
pub(super) trait Term: Copy {
    /// The arithmetic in which the series is summed
    type Value: SeriesReal;
    fn value(self) -> Self::Value;
    fn bound(self) -> f64;
    /// Additional absolute rounding error carried by the generated term itself
    fn error(self) -> f64 {
        0.0
    }
}

impl<V: SeriesReal> Term for V {
    type Value = Self;

    fn value(self) -> V {
        self
    }

    fn bound(self) -> f64 {
        self.abs_up()
    }
}

/// A value and a bound on its size
impl Term for (f64, f64) {
    type Value = f64;

    fn value(self) -> f64 {
        self.0
    }

    fn bound(self) -> f64 {
        self.1
    }
}

/// A double-double term and a bound on its generation error
impl Term for (Dd, f64) {
    type Value = Dd;

    fn value(self) -> Dd {
        self.0
    }

    fn bound(self) -> f64 {
        let mag = dd_hi_up(self.0.abs());
        let b = if mag > 0.0 && self.1 > 0.0 {
            next_up_nonneg(mag + self.1)
        } else {
            mag + self.1
        };
        if b.is_finite() { b } else { f64::INFINITY }
    }

    fn error(self) -> f64 {
        self.1
    }
}

/// Incremental state for `asymptotic_sum` and each residue class of `pfq_exp_cos`
#[derive(Clone, Copy)]
struct AsymptoticSum<V> {
    settled: f64,
    z_abs: f64,
    remainder: f64,
    tol: f64,
    sum: V,
    abs_sum: f64,
    term_err: f64,
    prev: f64,
    min: f64,
    best: (V, f64),
    done: bool,
}

impl<V: SeriesReal> AsymptoticSum<V> {
    fn new(settled: f64, z: f64, tol: f64) -> Self {
        let zero = V::from_f64(0.0);
        let z_abs = z.abs();
        Self {
            settled,
            z_abs,
            remainder: (TAU * z_abs).sqrt().max(1.0),
            tol,
            sum: zero,
            abs_sum: 0.0,
            term_err: 0.0,
            prev: f64::INFINITY,
            min: f64::INFINITY,
            best: (zero, f64::INFINITY),
            done: false,
        }
    }

    /// Advances by term `n` (with `next_bound` from the lookahead term); breaks when done
    fn step(
        &mut self,
        n: f64,
        t: V,
        bound: f64,
        term_err: f64,
        next_bound: f64,
        mut ratio_bound: impl FnMut(f64) -> f64,
    ) -> ControlFlow<()> {
        let size = bound.max(next_bound);
        let stop = !t.to_f64().is_finite()
            || !bound.is_finite()
            || n >= self.settled + self.z_abs && bound >= self.prev;
        let converged = !stop && n >= self.settled && {
            // a later sum can disqualify the current best even when its term is larger
            let err = self.remainder * size + n * V::EPSILON * self.abs_sum + self.term_err;
            self.best.1 = self
                .best
                .1
                .max((self.sum - self.best.0).to_f64().abs() - err);
            if err < self.best.1 {
                self.best = (self.sum, err);
            }
            if size < self.min {
                self.min = size;
                self.remainder * self.min <= self.tol * self.sum.to_f64().abs()
                    && ratio_bound(n) < 1.0
            } else {
                false
            }
        };
        if stop || converged {
            self.done = true;
            ControlFlow::Break(())
        } else {
            self.sum = self.sum + t;
            self.abs_sum += bound;
            self.term_err += term_err;
            self.prev = bound;
            // (at the cap, terms that grew after the smallest may still shrink to smaller ones)
            if n + 1.0 == f64::from(PFQ_MAXITER) && size > self.min {
                self.best.1 = f64::INFINITY;
            }
            ControlFlow::Continue(())
        }
    }
}

/// `asymptotic_sum` stopping once the remainder bound is at most `tol * |sum|`
fn asymptotic_sum_tol<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    z: f64,
    tol: f64,
    mut ratio_bound: impl FnMut(f64) -> f64,
) -> (T::Value, f64) {
    let mut state = AsymptoticSum::<T::Value>::new(settled, z, tol);
    let mut terms = terms.peekable();
    let _ = (0..PFQ_MAXITER)
        .map(f64::from)
        .try_for_each(|n| match terms.next() {
            Some(term) => {
                let next = terms.peek().map_or(0.0, |next| next.bound());
                state.step(n, term.value(), term.bound(), term.error(), next, &mut ratio_bound)
            },
            None => ControlFlow::Break(()),
        });
    state.best
}

/// The sum of an asymptotic series in `1/z` up to about its smallest term from term `settled`
/// on, and an estimate of its absolute error
///
/// From term `settled` on, all parameters plus `n` must be positive. Before it, the terms may
/// vanish or grow. From term `settled + |z|` on, the term ratios are about `n / |z|`, so that
/// growing terms keep growing. The size of a term `t_n` is the larger of the bounds on `t_n` and
/// on `t_{n+1}`, as a `t_n` of a recurrence of order 2 can vanish while `t_{n+1}` does not.
///
/// The error is the smallest size times `sqrt(2 pi |z|)`, about the sum of the terms near
/// it, plus rounding. For terms of one sign, the remainder can be half of that, as for
/// `sum_n n! / z^(n+1)` with `z > 0`, whose term ratios are also about `n / z`. The sum is the
/// one whose error is smallest, as the rounding grows with `n`. It stops where the remainder is
/// negligible if `ratio_bound(n) < 1`, so that no later term is larger before they grow for good.
///
/// Before term `settled + |z|`, the terms can grow again after a small one, and then shrink to
/// smaller ones, e.g. after a tiny `a_m` or a small positive `1 + a_m - a_i + n` of a term of
/// the algebraic part of pFq. So the error of a sum is taken to be at least its distance to any
/// later sum, less the error of that one (a heuristic, as that error is an estimate too), and
/// infinite if `PFQ_MAXITER` ends the terms while the last is larger than the smallest.
pub(super) fn asymptotic_sum<T: Term>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    z: f64,
    ratio_bound: impl FnMut(f64) -> f64,
) -> (T::Value, f64) {
    asymptotic_sum_tol(terms, settled, z, T::Value::EPSILON, ratio_bound)
}

/// The sum of a convergent series, in the arithmetic of its terms, and an estimate of its
/// absolute error
///
/// Stops once the geometric tail bound is at most `tol * |sum|` (or returns the bound at
/// `PFQ_MAXITER` if `ratio_bound < 1`, else infinite error). From term `settled` on, the bound
/// on each term is at most `ratio_bound(n)` times the one before it, and `ratio_bound` does not
/// increase. `step_roundings` bounds the rounding errors accumulated per recurrence step when
/// generating each term (beyond the summation rounding itself).
fn convergent_sum<T: Term>(
    mut terms: impl Iterator<Item = T>,
    settled: f64,
    tol: f64,
    step_roundings: f64,
    mut ratio_bound: impl FnMut(f64) -> f64,
) -> (T::Value, f64) {
    let eps = T::Value::EPSILON;
    let init = (T::Value::from_f64(0.0), 0.0, 0.0, f64::INFINITY);
    let step = |(sum, abs_sum, gen_err, _): (T::Value, f64, f64, f64), n: f64| match terms.next() {
        Some(term) if term.value().to_f64().is_finite() => {
            let next = sum + term.value();
            if !next.to_f64().is_finite() {
                ControlFlow::Break((next, f64::INFINITY))
            } else {
                let bound = term.bound();
                let term_gen = n * step_roundings * eps * bound;
                let abs_sum = abs_sum + bound;
                let gen_err = gen_err + term_gen;
                let rho = if n >= settled { ratio_bound(n) } else { 1.0 };
                if rho < 1.0 {
                    // geometric bound on the remaining terms, inflating the last term by its
                    // generation uncertainty and rounding upward so `rho` near 1 or tiny `Dd`
                    // terms cannot understate the tail
                    let tail = if bound == 0.0 || rho <= 0.0 {
                        0.0
                    } else {
                        let tail_base = if term_gen > 0.0 {
                            next_up_nonneg(bound + term_gen)
                        } else {
                            bound
                        };
                        next_up_nonneg(next_up_nonneg(tail_base * rho) / (1.0 - rho))
                    };
                    let err = tail + n * eps * abs_sum + gen_err;
                    if tail <= tol * next.to_f64().abs() {
                        ControlFlow::Break((next, err))
                    } else {
                        ControlFlow::Continue((next, abs_sum, gen_err, err))
                    }
                } else {
                    ControlFlow::Continue((next, abs_sum, gen_err, f64::INFINITY))
                }
            }
        },
        _ => ControlFlow::Break((sum, f64::INFINITY)),
    };
    match (0..PFQ_MAXITER).map(f64::from).try_fold(init, step) {
        ControlFlow::Break(res) => res,
        ControlFlow::Continue((sum, _, _, last)) => (sum, last),
    }
}

/// Sum of a terminating series of known `degree`, with its rounding error estimate
///
/// `step_roundings` bounds the rounding errors accumulated per recurrence step when generating
/// each term (beyond the summation rounding itself).
fn polynomial_sum<T: SeriesReal>(
    mut terms: impl Iterator<Item = T>,
    degree: f64,
    step_roundings: f64,
) -> (T, f64) {
    let zero = T::from_f64(0.0);
    let step = |(sum, err_sum): (T, f64), k: u32| match terms.next() {
        Some(t) if t.to_f64().is_finite() => {
            let next = sum + t;
            if next.to_f64().is_finite() {
                let weight = ((degree + 1.0) + f64::from(k) * step_roundings) * T::EPSILON;
                ControlFlow::Continue((next, err_sum + weight * t.to_f64().abs()))
            } else {
                ControlFlow::Break(next)
            }
        },
        _ => ControlFlow::Break(sum),
    };
    if degree >= f64::from(PFQ_MAXITER) {
        (zero, f64::INFINITY)
    } else {
        match (0..=degree as u32).try_fold((zero, 0.0), step) {
            ControlFlow::Continue((sum, err_sum)) => (sum, err_sum),
            ControlFlow::Break(sum) => (sum, f64::INFINITY),
        }
    }
}

/// The sum of the series in `1 / z` of a term of the algebraic part, where the term ratios are
/// about `n^kappa / |z|`: by `asymptotic_sum`, or by `convergent_sum` for `kappa = 0`, which
/// also stands for `p > q + 1`
fn pfq_alg_sum<T: Term<Value = f64>>(
    terms: impl Iterator<Item = T>,
    settled: f64,
    kappa: usize,
    z: f64,
    step_roundings: f64,
    ratio_bound: impl FnMut(f64) -> f64,
) -> (f64, f64) {
    if kappa == 0 {
        convergent_sum(terms, settled, f64::EPSILON, step_roundings, ratio_bound)
    } else {
        let az = z.abs();
        let scale = match kappa {
            1 => az,
            2 => az.sqrt(),
            _ => az.powf(1.0 / kappa as f64),
        };
        asymptotic_sum(terms, settled, scale, ratio_bound)
    }
}

/// Whether `up` and `low` (complete `Dd` denominator) need `frexp` rescaling in `pfq_next_dd`
fn pfq_extreme_dd(up: &[Dd], low: &[Dd]) -> bool {
    let many = up.len().max(low.len()) > PFQ_MAX_UNCHECKED;
    let extreme = |x: &Dd| {
        !PFQ_MODERATE.contains(&x.hi.abs()) || (x.lo != 0.0 && !PFQ_MODERATE.contains(&x.lo.abs()))
    };
    many || up.iter().chain(low).any(extreme)
}

/// `t_{n+1}` after `term = t_n` for complete `Dd` parameters `up` and `low`, preserving `Dd` low
/// parts through `(x + n).frexp()` when intermediate products are extreme or non-normal
#[cold]
#[inline(never)]
fn pfq_rescaled_next_dd(up: &[Dd], low: &[Dd], n: f64, w: f64, term: f64) -> f64 {
    let product = |params: &[Dd]| {
        params.iter().fold((1.0_f64, 0_i64), |(m, e), &x| {
            let (vm, ve) = (x + n).frexp();
            let (m, me) = frexp_f64(m * vm.hi);
            (m, e + i64::from(ve) + i64::from(me))
        })
    };
    let ((num, ne), (den, de)) = (product(up), product(low));
    let ((wm, we), (tm, te)) = (frexp_f64(w), frexp_f64(term));
    ldexp_i64(tm * ((num / den) * wm), (ne - de) + i64::from(we + te))
}

/// `t_{n+1}` after `term = t_n` for complete `Dd` parameters `up` and `low`
#[inline]
fn pfq_next_dd(up: &[Dd], low: &[Dd], n: f64, w: f64, term: f64, checked: bool) -> f64 {
    let normal = |v: f64| is_normal(v).then_some(v);
    let fast = (!checked)
        .then(|| {
            let num = normal(up.iter().fold(1.0_f64, |p, &x| p * (x + n).hi))?;
            let den = normal(low.iter().fold(1.0_f64, |p, &x| p * (x + n).hi))?;
            let ratio = normal(normal(num / den)? * w)?;
            Some(term * ratio)
        })
        .flatten();
    fast.unwrap_or_else(|| pfq_rescaled_next_dd(up, low, n, w, term))
}

/// The series of the term of `a_m` of the algebraic part (DLMF 16.11.2), as its sum and
/// absolute error
///
/// The shifted parameters `a_m`, `1 + a_m - b_j`, and `1 + a_m - a_i` (including `1` for `i = m`)
/// are kept in double-double so cancellations in `x + n`, pole checks, settlement, and bounds do
/// not lose their low parts. With `sizes`, it is the sum of the sizes of its terms.
fn pfq_alg_series(a: &[f64], b: &[f64], m: usize, z: f64, sizes: bool) -> (f64, f64) {
    let am = Dd::new(a[m]);
    let (mut up_s, mut low_s) =
        (Scratch::new(b.len() + 1, Dd::ZERO), Scratch::new(a.len(), Dd::ZERO));
    let up = up_s.fill(b.len() + 1, once(am).chain(b.iter().map(|&bj| (am - bj) + 1.0)));
    let low = low_s.fill(a.len(), a.iter().map(|&ai| (am - ai) + 1.0));
    sort_descending_dd(up);
    sort_descending_dd(low);
    let degree = first_zero_dd(up);
    let (p, q) = (a.len(), b.len());
    let w = neg_one_pow(p + q + 1) / z;
    let checked = pfq_extreme_dd(up, low);
    let series = once(1.0)
        .chain(pfq_terms_of(degree, w, |n, term| pfq_next_dd(up, low, n, w, term, checked)))
        .map(|t| if sizes { t.abs() } else { t });
    // each step does `up.len() + low.len()` roundings and multiplications plus division and
    // scaling (`<= 2 * (up + low) + 3` `f64` roundings of `0.5 * f64::EPSILON` each)
    let step_roundings = (up.len() + low.len() + 2) as f64;
    if degree.is_finite() {
        polynomial_sum(series, degree, step_roundings)
    } else {
        let settled = up
            .iter()
            .chain(&*low)
            .map(positive_from)
            .fold(0.0, f64::max);
        let kappa = (q + 1).saturating_sub(p);
        let mut scratch = BoundScratch::new(up.len(), low.len());
        let bound = |n| scratch.regrowth_bound(up, low, recip_up(z), n);
        pfq_alg_sum(series, settled, kappa, z, step_roundings, bound)
    }
}

/// `pfq_alg_series` for `p >= q` in double-double, with the exact parameters, and an infinite
/// error if no term bounds the rest
///
/// A product of the factors of a term ratio that isn't normal is less exact, and ends the terms.
fn pfq_alg_series_dd(a: &[f64], b: &[f64], m: usize, z: f64) -> (Dd, f64) {
    let [mut up, mut low] = residue::series_params(a, b, m);
    // (so that the factors of a term ratio pair by size, with the smallest of `up` unpaired)
    for v in [&mut up, &mut low] {
        sort_descending_dd(v);
    }
    let degree = first_zero_dd(&up);
    let ratio = |n: f64| {
        let normal = |p: Dd| if is_normal(p) { p } else { Dd::NAN };
        let prod = |v: &[Dd]| v.iter().fold(Dd::ONE, |p, &x| normal(p * (x + n)));
        // in (-1)^(p+q+1) / z
        if n == degree {
            Dd::ZERO
        } else {
            prod(&up) / normal(prod(&low) * (neg_one_pow(a.len() + b.len() + 1) * z))
        }
    };
    let terms = once(Dd::ONE).chain(pfq_terms_of(degree, 1.0, |n, t| t * ratio(n)));
    // each double-double step does `up.len() + low.len()` additions (`<= 3 * 2^-106`) and
    // multiplications (`<= 7 * 2^-106`) plus 1 division (`<= 15 * 2^-106`), bounded in units of
    // `Dd::EPSILON = 64 * 2^-106` by `0.25 * (up + low + 2)` (`(16 * (up + low) + 32) * 2^-106`)
    let step_roundings = 0.25 * (up.len() + low.len() + 2) as f64;
    if degree.is_finite() {
        // (all of its terms, with per-step factor-product roundings in double-double)
        polynomial_sum(terms, degree, step_roundings)
    } else if up.len() <= low.len() {
        // (it converges for p > q + 1, and for p = q + 1 with |z| > 1)
        let mut scratch = BoundScratch::new(up.len(), low.len());
        let bound = |n: f64| scratch.ratio_bound(&up, &low, recip_up(z), n);
        let settled = up.iter().chain(&low).map(positive_from).fold(0.0, f64::max);
        convergent_sum(terms, settled, Dd::EPSILON, step_roundings, bound)
    } else {
        // Where its factors are positive, each paired one of a term ratio is monotone with the limit
        // 1, and the unpaired one grows. So for `m >= n` the ratio is at most a `w(m)` that is linear
        // in `m`, with `w(n) = above`, and at least `w(m) / spread`.
        //
        // The product of `w(n)` to `w(m)` is at most 1 while `w(n) + w(m) <= 2`: no term up to there
        // is larger than term `n`. And from `w(m) >= spread` on, the terms grow for good. The sum
        // only ends at an `n` with `above + spread <= 2` (above their rounding), and so not at a
        // term that a later one exceeds.
        let ends = |n: f64| {
            let shift = |x: &Dd| (*x + n).hi;
            let unpaired = shift(&up[low.len()]) / z.abs();
            let quotients = zip(&up, &low).map(|(u, l)| shift(u) / shift(l));
            let [above, spread] = quotients.fold([unpaired, 1.0], |[above, spread], f| {
                [above * f.max(1.0), spread * f.max(f.recip())]
            });
            let positive = up.iter().chain(&low).all(|x| shift(x) > 0.0);
            positive
                && unpaired.is_normal()
                && (above + spread) * (1.0 + f64::EPSILON.sqrt()) <= 2.0
        };
        // (the unpaired factor alone is above 1 after that)
        let mut ns = (0..PFQ_MAXITER)
            .map(f64::from)
            .take_while(|&n| (up[low.len()] + n).hi <= z.abs());
        match ns.find(|&n| ends(n)) {
            Some(start) => {
                let terms = zip(terms, (0_u32..).map(f64::from))
                    .take_while(|&(_, n)| n < start || ends(n))
                    .map(|(term, _)| term);
                asymptotic_sum(terms, start, z, |_| 0.0)
            },
            None => (Dd::ONE, f64::INFINITY),
        }
    }
}

/// Appends the term of `pfq_alg` of the class of `a_m` (the `a_i` that differ from `a_m` by
/// integers) to `terms`
///
/// The terms of these `a_i` are infinite, unlike the limit of their sum, which is the sum of
/// the `Residues` at the poles that they share. It pushes one term for the first of the smallest
/// of them, with the residue `(-1)^k / k!` of the Gamma function in place of each `Gamma(-k)`
/// of its factor, and leaves `terms` unchanged for the other `a_i` or if the class has no
/// residues.
///
/// The residues are complex for `z > 0`, where it is the sum of the bounds on their sizes, as
/// with `sizes`.
///
/// For a polynomial, returns `None` for the class of its non-positive integer `a_i`, which the
/// `Residues` exclude, and leaves `terms` unchanged (returning `Some(())`) for the other classes,
/// whose terms vanish with the `1 / Gamma(a_i)` of their factor.
fn pfq_class_term(
    a: &[f64],
    b: &[f64],
    m: usize,
    z: f64,
    sizes: bool,
    terms: &mut Vec<Scaled>,
    cache: LogGammaCache<'_>,
) -> Option<()> {
    let am = a[m];
    let others = skip(a, m).map(Dd::new);
    let diffs = others.clone().map(|ai| ai - am);
    let before = |(d, i): (Dd, usize)| d.is_integer() && (d.hi < 0.0 || d.hi == 0.0 && i < m);
    if diffs.clone().zip(0..).any(before) {
        return Some(());
    }
    if first_zero(a).is_finite() {
        return (!is_nonpos_int(am)).then_some(());
    }
    let residues = Residues::of_class(a, b, m, z);
    let [up, low] = residues.params();
    // (their order is final after the last zero factor, and their bounds are valid once all
    // parameters plus n exceed 1)
    let settled = up
        .iter()
        .chain(&low)
        .map(|&x| positive_from(&(x - 1.0)))
        .fold(0.0, f64::max);
    let mut scratch = BoundScratch::new(up.len(), low.len());
    let bound =
        |n| scratch.ratio_bound(&up, &low, recip_up(z), n) * Residues::growth_bound(&up, &low, n);
    let sizes = sizes || z > 0.0;
    let series = residues.map(|(t, size)| (if sizes { size } else { t }, size));
    let kappa = (b.len() + 1).saturating_sub(a.len());
    let (sum, err) = pfq_alg_sum(series, settled, kappa, z, 0.0, bound);
    if (sum, err) != (0.0, 0.0) {
        let bs = b.iter().copied().map(Dd::new);
        let args = bs.clone().map(|bj| bj - am);
        let poles = args.clone().filter(is_pole);
        let sign = neg_one_pow(poles.clone().filter(|g| g.hi % 2.0 != 0.0).count());
        // (Gamma(a_i - a_m) has the residue 1 for a_i = a_m)
        let num = bs
            .chain(diffs.filter(|d| d.hi != 0.0))
            .chain(poles.map(|g| -g + 1.0));
        let den = others.chain(args.filter(|g| !is_pole(g)));
        let l = -ln_dd(z.abs()) * am;
        let mode = if sizes { Alg::Sizes } else { Alg::Value };
        gamma_alg_term(num, den, (Dd::new(sign * sum), err), l, mode, terms, cache)?;
    }
    Some(())
}

/// What `pfq_alg` gives for each term
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Alg {
    /// Its value
    Value,
    /// A bound on its size
    Sizes,
    /// Its value for `p >= q`, with its series and its factor in double-double
    Dd,
}

/// The non-vanishing terms of the algebraic part of the asymptotic expansion of pFq for `p <= q`
/// (DLMF 16.11.2), with `|z|` instead of `-z`, or `None` if a factor of a term overflows or
/// underflows, or if a non-positive integer `a_i` differs from another `a_i` by an integer
///
/// For `p > q + 1` and `z < 0`, their series converge, and they sum to `pfq_borel`.
///
/// The term of `a_m` is `|z|^-a_m` times a `_{q+1}F_{p-1}` series in `(-1)^(p+q+1) / z`, with the
/// parameters `a_m` and `1 + a_m - b_j`, and `1 + a_m - a_i` (`i != m`). It is omitted if it
/// vanishes, where `1 / Gamma(b_j - a_m)` does. The `a_i` that differ by integers have one term,
/// by `pfq_class_term`, which for `z > 0` is a bound on its size.
///
/// With `Alg::Sizes`, each term is the sum of the sizes of the terms of its series, times the
/// size of its factor. Their sum bounds the algebraic part at any phase of `z` and of its terms.
///
/// With `Alg::Dd`, the series and the factor of each term are in double-double, except for the
/// `a_i` that differ by integers.
pub(super) fn pfq_alg(a: &[f64], b: &[f64], z: f64, mode: Alg) -> Option<Vec<Scaled>> {
    pfq_alg_with_cache(a, b, z, mode, &mut [])
}

fn pfq_alg_with_cache(
    a: &[f64],
    b: &[f64],
    z: f64,
    mode: Alg,
    cache: LogGammaCache<'_>,
) -> Option<Vec<Scaled>> {
    let sizes = mode == Alg::Sizes;
    let ln_z = ln_dd(z.abs());
    let bs = b.iter().copied().map(Dd::new);
    let mut terms = Vec::with_capacity(a.len());
    for (m, &am) in a.iter().enumerate() {
        let others = skip(a, m).map(Dd::new);
        let diffs = others.clone().map(|ai| ai - am);
        if diffs.clone().any(Dd::is_integer) {
            pfq_class_term(a, b, m, z, sizes, &mut terms, cache)?;
        } else {
            let num = bs.clone().chain(diffs);
            let den = others.chain(bs.clone().map(|bj| bj - am));
            if !den.clone().any(|g| is_pole(&g)) {
                let l = -ln_z * am;
                let (sum, err) = if mode == Alg::Dd {
                    pfq_alg_series_dd(a, b, m, z)
                } else {
                    let (sum, err) = pfq_alg_series(a, b, m, z, sizes);
                    (Dd::new(sum), err)
                };
                gamma_alg_term(num, den, (sum, err), l, mode, &mut terms, cache)?;
            }
        }
    }
    Some(terms)
}

/// Relative error bound of `x = kappa * root_dd(|z|, kappa)`
pub(super) fn pfq_exp_x_rel(kappa: usize) -> f64 {
    let kf = kappa as f64;
    if kappa == 1 {
        0.0
    } else if kappa == 2 {
        4.0 * EPSILON_SQ
    } else if kappa <= 1024 {
        // `root_dd` scales `u` into `[2^-513, 2^512]`, keeping `u` normal and `|ln u| / kappa <= ln 2`
        (4.0 * kf + 4.0) * EPSILON_SQ
    } else {
        // unscaled `root_dd` can lose double-double precision on extreme or subnormal inputs
        f64::INFINITY
    }
}

/// The terms `t_k = c_k x^-k` of the exponential part of the asymptotic expansion of
/// pFq for `p <= q` (DLMF 16.11.1)
///
/// With `D = x d/dx` and `kappa = q - p + 1`, the differential equation (DLMF 16.8.3) is
/// `(D prod_j (D + kappa (b_j - 1)) - x^kappa prod_i (D + kappa a_i)) w = 0`, and each factor
/// `D + g` maps `e^x x^s` to `e^x (x^(s+1) + (s + g) x^s)`. So `c_0 = 1` and
/// `kappa k c_k = sum_{m=1}^q d_{m+1}(nu - k + m) c_{k-m}`, where `d_j(s)` is the
/// coefficient of `e^x x^(s+q+1-j)` in the image of `e^x x^s`.
fn pfq_exp_terms(
    a: &[f64],
    b: &[f64],
    nu: Dd,
    nu_err: f64,
    x: Dd,
    x_rel: f64,
) -> impl Iterator<Item = (Dd, f64)> {
    let (q, kappa) = (b.len(), (b.len() + 1 - a.len()) as f64);
    let (xm, xe) = x.frexp();
    let left: Vec<(Dd, f64)> = b
        .iter()
        .map(|&bj| ((Dd::new(bj) - 1.0) * kappa, 2.0 * Dd::EPSILON * (bj.abs() + 1.0) * kappa))
        .chain(once((Dd::ZERO, 0.0)))
        .collect();
    let right: Vec<(Dd, f64)> = a
        .iter()
        .map(|&ai| (Dd::new(ai) * kappa, Dd::EPSILON * ai.abs() * kappa))
        .collect();
    // the coefficients of e^x x^(s+n-j) in prod_(i<n) (D + g_i) e^x x^s, for j <= q + 1,
    // propagating absolute rounding uncertainty from factor construction through cancellation
    let fill_coefs = |g: &[(Dd, f64)], s: Dd, s_err: f64, c: &mut [Dd], err_c: &mut [f64]| {
        c[0] = Dd::ONE;
        c[1..].fill(Dd::ZERO);
        err_c[0] = 0.0;
        err_c[1..].fill(0.0);
        for (i, &(gi, gi_err)) in g.iter().enumerate() {
            for j in (1..=i + 1).rev() {
                let shift = (i + 1 - j) as f64;
                let f = s + shift + gi;
                let f_err = s_err + gi_err + 3.0 * Dd::EPSILON * (s.hi.abs() + shift + gi.hi.abs());
                let prod = c[j - 1] * f;
                let prod_err = err_c[j - 1] * (f.hi.abs() + f_err)
                    + c[j - 1].hi.abs() * f_err
                    + Dd::EPSILON * prod.hi.abs();
                let sum_mag = c[j].hi.abs() + prod.hi.abs();
                c[j] = c[j] + prod;
                err_c[j] += prod_err + Dd::EPSILON * sum_mag;
            }
        }
    };
    // t_k, its error bound, and the sums/error bounds so far of t_(k+1), ..., t_(k+q)
    let (mut k, mut t, mut t_err) = (0.0, Dd::ONE, 0.0);
    let (mut sums, mut sums_err) = (vec![Dd::ZERO; q], vec![0.0; q]);
    let (mut l, mut le) = (vec![Dd::ZERO; q + 2], vec![0.0; q + 2]);
    let (mut r, mut re) = (vec![Dd::ZERO; q + 2], vec![0.0; q + 2]);
    from_fn(move || {
        let term = (t, t_err);
        let s = nu - k;
        let s_err = nu_err + Dd::EPSILON * (nu.hi.abs() + k);
        fill_coefs(&left, s, s_err, &mut l, &mut le);
        fill_coefs(&right, s, s_err, &mut r, &mut re);
        let (mut scale_m, mut scale_e) = (Dd::ONE, 0_i64);
        let (tm, te) = t.frexp();
        let (tem, tee) = frexp_f64(t_err);
        for (m, (sum, sum_err)) in zip(&mut sums, &mut sums_err).enumerate() {
            // keep `x^-(m+1)` split into mantissa and exponent so large `|x|` does not underflow
            // before multiplying by large `d * t`
            let (sm, se) = (scale_m / xm).frexp();
            scale_m = sm;
            scale_e += i64::from(se) - i64::from(xe);
            let d = l[m + 2] - r[m + 2];
            let d_err =
                le[m + 2] + re[m + 2] + Dd::EPSILON * (l[m + 2].hi.abs() + r[m + 2].hi.abs());
            let (dm, de) = d.frexp();
            let prod_m = dm * scale_m * tm;
            let prod_e = i64::from(de) + scale_e + i64::from(te);
            let contrib = ldexp_i64(prod_m, prod_e);
            let (dem, dee) = frexp_f64(d_err);
            let sm_abs = scale_m.hi.abs();
            let err_dt =
                ldexp_i64(dem * tm.hi.abs() * sm_abs, i64::from(dee) + i64::from(te) + scale_e);
            let err_ee = ldexp_i64(dem * tem * sm_abs, i64::from(dee) + i64::from(tee) + scale_e);
            let err_td =
                ldexp_i64(dm.hi.abs() * tem * sm_abs, i64::from(de) + i64::from(tee) + scale_e);
            let scale_rel = 2.0 * (m + 1) as f64 * (x_rel + Dd::EPSILON) + 3.0 * Dd::EPSILON;
            let subnormal_err = if prod_m.hi != 0.0 && contrib.hi.abs() < Dd::MIN_POSITIVE {
                f64::EPSILON * contrib.hi.abs() + f64::MIN_POSITIVE
            } else {
                0.0
            };
            let contrib_err =
                err_dt + err_ee + err_td + scale_rel * contrib.hi.abs() + subnormal_err;
            let sum_mag = sum.hi.abs() + contrib.hi.abs();
            *sum = *sum + contrib;
            *sum_err += contrib_err + Dd::EPSILON * sum_mag;
        }
        k += 1.0;
        let denom = kappa * k;
        t = sums[0] / denom;
        t_err = (sums_err[0] + Dd::EPSILON * sums[0].hi.abs()) / denom;
        sums.rotate_left(1);
        sums[q - 1] = Dd::ZERO;
        sums_err.rotate_left(1);
        sums_err[q - 1] = 0.0;
        Some(term)
    })
}

/// `cos(pi / k)` and `sin(pi / k)` in double-double precision, for `k >= 2`
fn cos_sin_pi_div(k: usize) -> (Dd, Dd) {
    if k == 2 {
        (Dd::ZERO, Dd::ONE)
    } else {
        let phi = 2.0 / (FRAC_2_PI_DD * k as f64);
        // n, and the terms n - 1 of their Taylor series, of the degrees 2 n - 2 and 2 n - 1, which
        // are below 2^-106 from n = 21 on
        let terms = successors(Some((1.0, Dd::ONE, phi)), |&(n, c, _)| {
            let next = -(c * phi * phi) / ((2.0 * n - 1.0) * (2.0 * n));
            Some((n + 1.0, next, next * phi / (2.0 * n + 1.0)))
        });
        terms
            .take(20)
            .fold((Dd::ZERO, Dd::ZERO), |(cos, sin), (_, c, s)| (cos + c, sin + s))
    }
}

/// The sum of the two exponential parts of `pfq_exp` for `z < 0` and `p < q`, which have
/// `x e^(+-pi i / kappa)` as `x` (DLMF 16.11.8, 16.11.9), relative to
/// `e^(x cos(pi / kappa))`:
///
/// ```text
/// 2 sum_k t_k cos(theta - pi k / kappa),
/// theta = x sin(pi / kappa) + pi nu / kappa
/// ```
///
/// with its absolute error, and `cos(pi / kappa)`, where `terms` are the `t_k = c_k x^-k`
/// and `root` is `|z|^(1 / kappa)`.
fn pfq_exp_cos(
    terms: impl Iterator<Item = (Dd, f64)>,
    settled: f64,
    kappa: usize,
    root: Dd,
    nu: Dd,
    nu_err: f64,
    x_rel: f64,
) -> ((f64, f64), Dd) {
    let mut terms = terms.fuse();
    let kf = kappa as f64;
    let x = kf * root.hi;
    let (cos, sin) = cos_sin_pi_div(kappa);
    // theta / pi in double-double, reduced to [-1, 1], and a bound of the error of theta
    let turns = root * kf * sin * (FRAC_2_PI_DD * 0.5) + nu / kf;
    let t = turns - 2.0 * (0.5 * turns.hi).round();
    let dtheta =
        4.0 * EPSILON_SQ * (x + PI * nu.hi.abs()) + (x * sin.hi).abs() * x_rel + PI * nu_err / kf;
    // single pass over `terms`: each residue class `k = r mod kappa` tracks its own
    // `AsymptoticSum` state in double-double, with a `kappa`-term lookahead window
    let mut classes = vec![AsymptoticSum::<Dd>::new(settled / kf, x / kf, Dd::EPSILON); kappa];
    let mut win: Vec<Option<(Dd, f64)>> = (0..kappa).map(|_| terms.next()).collect();
    let mut active = kappa;
    let _ = (0..PFQ_MAXITER).try_for_each(|j| {
        let (nf, sign) = (f64::from(j), if j % 2 == 0 { 1.0 } else { -1.0 });
        for (cls, slot) in zip(&mut classes, &mut win) {
            let raw = core::mem::replace(slot, terms.next());
            if !cls.done {
                let next = slot.map_or(0.0, Term::bound);
                match raw {
                    None => {
                        cls.done = true;
                        active -= 1;
                    },
                    Some(term @ (v, v_err))
                        if cls
                            .step(nf, v * sign, term.bound(), v_err, next, |_| 0.0)
                            .is_break() =>
                    {
                        active -= 1;
                    },
                    Some(_) => {},
                }
            }
        }
        if active == 0 {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    });
    let (cos0, sin0) = if kappa == 2 {
        (cospi(t.hi), sinpi(t.hi))
    } else {
        (0.0, 0.0)
    };
    let (sum, err, size) =
        classes
            .iter()
            .enumerate()
            .fold((Dd::ZERO, 0.0, 0.0), |(sum, err, size), (r, cls)| {
                let (s, s_err) = cls.best;
                let (c, u_err) = if kappa == 2 {
                    let c = if r == 0 {
                        Dd::new(cos0) - Dd::new(PI) * t.lo * sin0
                    } else {
                        Dd::new(sin0) + Dd::new(PI) * t.lo * cos0
                    };
                    (c, 0.0)
                } else {
                    let u = t.hi - r as f64 / kf;
                    (Dd::new(cospi(u)) - Dd::new(PI) * t.lo * sinpi(u), 2.0 * f64::EPSILON)
                };
                let c_err = dtheta + u_err;
                let rounding = c_err * s.hi.abs() + 3.0 * f64::EPSILON * (s.hi * c.hi).abs();
                (
                    sum + s * c,
                    err + s_err * (c.hi.abs() + c_err) + rounding,
                    size + s.hi.abs() + s_err,
                )
            });
    // the other exponentials of DLMF 16.11.9, with x e^(+-3 pi i / kappa), relative to these
    let next = 2.0 * (x * (cospi(3.0 / kf) - cos.hi)).exp();
    let others = if kappa == 2 { 0.0 } else { next };
    ((2.0 * sum.hi, 2.0 * (err + others * size)), cos)
}

/// The exponential part of the asymptotic expansion of pFq for `p <= q` (DLMF 16.11.1),
/// with `|z|^nu` instead of `z^nu`:
///
/// ```text
/// prod_j Gamma(b_j) / prod_i Gamma(a_i)
///     (2 pi)^((1 - kappa) / 2) kappa^(-nu - 1/2) e^x |x|^nu sum_k c_k x^-k,
/// nu = sum_i a_i - sum_j b_j + (kappa - 1) / 2
/// ```
///
/// with `kappa = q - p + 1` and `x = kappa z^(1 / kappa)`, which is real for `z > 0` or
/// `p = q`. Otherwise, it is the sum of the two parts with
/// `x = kappa |z|^(1 / kappa) e^(+-pi i / kappa)`.
fn pfq_exp(a: &[f64], b: &[f64], z: f64, cache: LogGammaCache<'_>) -> Option<Scaled> {
    let kappa = b.len() + 1 - a.len();
    let x_rel = pfq_exp_x_rel(kappa);
    if !x_rel.is_finite() {
        return None;
    }
    let kf = kappa as f64;
    let nu = sum_dd(a) - sum_dd(b) + 0.5 * (kf - 1.0);
    let nu_abs = a.iter().chain(b).map(|v| v.abs()).sum::<f64>() + 0.5 * (kf - 1.0);
    let nu_err = (a.len() + b.len() + 2) as f64 * Dd::EPSILON * nu_abs;
    let real = z > 0.0 || kappa == 1;
    let root = root_dd(z.abs(), kappa) * if real { z.signum() } else { 1.0 };
    let x_dd = root * kf;
    let x = x_dd.hi;
    // (the factors of the last coefficient are negative from here on)
    let beta = b.iter().fold(0.0_f64, |m, bj| m.max(kf * (bj - 1.0)));
    let settled = (nu.hi + b.len() as f64 + beta).max(0.0);
    let terms = pfq_exp_terms(a, b, nu, nu_err, x_dd, x_rel);
    let (series, growth) = if real {
        let (sum, err) = asymptotic_sum_tol(terms, settled, x, f64::EPSILON, |_| 0.0);
        // the other exponentials of DLMF 16.11.8 and 16.11.9, with x e^(+-2 pi i / kappa),
        // relative to this one
        let next = 2.0 * (x * (cospi(2.0 / kf) - 1.0)).exp();
        let others = if kappa == 1 { 0.0 } else { next };
        ((sum.hi, err + (others + f64::EPSILON) * sum.hi.abs()), Dd::ONE)
    } else {
        pfq_exp_cos(terms, settled, kappa, root, nu, nu_err, x_rel)
    };
    // (2 pi)^((1 - kappa) / 2) kappa^(-1/2) |z|^(nu / kappa) e^x
    let ln_power = nu / kf * ln_dd(z.abs()) - ln_dd(kf) * 0.5;
    let l = ln_power - LN_SQRT_2PI_DD * (kf - 1.0) + root * kf * growth;
    let num = b.iter().copied().map(Dd::new);
    let den = a.iter().copied().map(Dd::new);
    gamma_series(num, den, series, l, cache)
}

/// pFq for `p <= q` by its asymptotic expansion (DLMF 16.11.7 to 16.11.9)
///
/// Only the exponentially large part is used for `z > 0`, and only the algebraic part for
/// `z < 0` and `p = q`. The size of the other part is added to the error, for `z > 0` as the
/// sizes of the terms of the algebraic part. For `z < 0` and `p < q`, it is the sum of both
/// parts.
pub(super) fn pfq_asymptotic(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let kappa = b.len() + 1 - a.len();
    // A 2F2 evaluation reuses the same Gamma arguments across the algebraic and exponential terms.
    let mut log_gamma;
    let cache: LogGammaCache<'_> = if a.len() == 2 && b.len() == 2 {
        log_gamma = [None; 8];
        &mut log_gamma
    } else {
        &mut []
    };
    let mode = if z > 0.0 { Alg::Sizes } else { Alg::Value };
    let alg = pfq_alg_with_cache(a, b, z, mode, cache)?;
    let exp = pfq_exp(a, b, z, cache)?;
    if z > 0.0 {
        // distinguish an empty algebraic part (zero error) from a failed `add_exp`
        let other = if alg.is_empty() {
            None
        } else {
            Some(add_exp(alg)?)
        };
        Some(Estimate::scaled(exp, other))
    } else if kappa == 1 {
        Some(Estimate::scaled(add_exp(alg)?, Some(exp)))
    } else {
        // (DLMF 16.11.9 for p <= q - 2 has the exponential part only; with the algebraic part,
        // it is theorem 3 of Paris, https://arxiv.org/abs/1708.04824)
        Some(Estimate::scaled(add_exp(alg.into_iter().chain([exp]))?, None))
    }
}

/// The sum of the terms of `pfq_alg` in double-double, for `p >= q`
pub(super) fn pfq_alg_dd(a: &[f64], b: &[f64], z: f64) -> Option<Scaled> {
    add_exp_dd(pfq_alg(a, b, z, Alg::Dd)?)
}

/// pFq for `p = q` and `z < 0` by the algebraic part of `pfq_asymptotic` in double-double, where
/// its terms cancel, e.g. for nearly equal `a_i`, or the terms of their series do
pub(super) fn pfq_asymptotic_dd(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    if z > 0.0 || a.len() != b.len() {
        return None;
    }
    let alg = pfq_alg_dd(a, b, z)?;
    Some(Estimate::scaled(alg, Some(pfq_exp(a, b, z, &mut [])?)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::pow2;

    #[test]
    fn test_asymptotic_sum_cap_after_hump() {
        // small terms, a hump before term `settled + |z|`, then terms that shrink too slowly to
        // end before the cap
        let term = |n: u32| {
            let t = match n {
                0..5 => 10.0_f64.powi(-i32::try_from(n).unwrap()),
                5..10 => 10.0_f64.powi(i32::try_from(n).unwrap() - 7),
                _ => 1e3 * 0.9999_f64.powi(i32::try_from(n).unwrap() - 10),
            };
            (t, t)
        };
        let (_, err) = asymptotic_sum((0..).map(term), 0.0, 20.0, |_| 0.0);
        assert_eq!(err, f64::INFINITY);
    }

    #[test]
    fn test_pfq_alg_series_dd_growing() {
        // the terms of the series of a_1 grow from the fourth one on, long before `n = |z|`
        let (a, b) = ([pow2(-70), -19.75, -29.625], [1.0, -219.75, -229.625]);
        assert!(pfq_alg_series_dd(&a, &b, 0, -100.0).1.is_infinite());
        // the low parts of the parameters make the term after the one where they settle larger
        let (s, u) = (pow2(-49), pow2(-52));
        let (a, b) = ([s - 16.0, 1.0 + 3.0 * u], [1.0 + u, 2.0 * s - 16.0]);
        assert!(pfq_alg_series_dd(&a, &b, 0, -1.1 * s).1.is_infinite());
        // the bound of the unpaired factor underflows, and the terms grow to ~1e72 from the fifth
        let c = (1..=10).map(|i| f64::from(i) / 16.0);
        let a: Vec<_> = once(f64::from_bits(1)).chain(c.clone()).collect();
        let b: Vec<_> = once(3.0).chain(c.map(|c| c - pow2(49))).collect();
        assert!(pfq_alg_series_dd(&a, &b, 0, -2.0).1.is_infinite());
    }

    #[test]
    fn test_pfq_alg_series_dd_polynomial() {
        // of degree 0, which is below `|z|`
        let (sum, err) = pfq_alg_series_dd(&[0.5, 0.75], &[1.5, 0.25], 0, -0.5);
        assert!(sum.hi == 1.0 && err <= Dd::EPSILON);
    }

    #[test]
    fn test_pfq_asymptotic_integer_difference() {
        // the bound on the algebraic part is too large for equal a_i, and near them the
        // expansion is accurate or `None` (Arb reference)
        let c = 1.0 + pow2(-10);
        let (b, z) = ([c, c, 11.25], 21.25);
        let accurate = |a: &[f64; 3]| pfq_asymptotic(a, &b, z).and_then(Estimate::accurate);
        assert_eq!(accurate(&[1.0, 1.0, 10.5]), None);
        let near = accurate(&[1.0, 1.0 + pow2(-20), 10.5]);
        let expected = 739_413_030.005_544_4;
        assert!(near.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected));
    }

    #[test]
    fn test_cos_sin_pi_div_two() {
        let (cos, sin) = cos_sin_pi_div(2);
        assert_eq!((cos.hi, cos.lo, sin.hi, sin.lo), (0.0, 0.0, 1.0, 0.0));
    }

    #[test]
    fn test_pfq_asymptotic_nearly_cancelled() {
        // equal a_i within sqrt(eps) of b_2 and b_3, where the algebraic part is about 1 and
        // changes much with a_i: accurate or `None` (Arb reference)
        let t = -pow2(-40);
        let (root, root_2) = (f64::EPSILON.sqrt(), (2.0 * f64::EPSILON).sqrt());
        let y = pfq_asymptotic(&[t; 3], &[t - 1.0, t + root, t + root_2], 30.0)
            .and_then(Estimate::accurate);
        let expected = -816_436.169_307_593_7_f64;
        assert!(y.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected.abs()));
    }

    #[test]
    fn test_pfq_alg_dd_params() {
        // exact terminating parameter (-2) beside one with identical hi word and nonzero lo
        // (-2 - 2^-52): degree is 2, not infinite or 0
        let [up, _] = residue::series_params(&[-2.0, 0.5], &[1.0 + pow2(-52)], 0);
        assert_eq!(up[0].hi, -2.0);
        assert_eq!(up[1].hi, -2.0);
        assert_ne!(up[1].lo, 0.0);
        assert_eq!(first_zero_dd(&up), 2.0);
        let (sum, err) = pfq_alg_series(&[-2.0, 0.5], &[1.0 + pow2(-52)], 0, -10.0, false);
        // 1 + (-2)(-2 - 2^-52)/(1 * -1.5) * (-0.1) + ...
        let expected = 1.0 + 4.0 / 15.0 + 4.0 / 150.0;
        assert!((sum - expected).abs() <= 1e-13 && err <= 1e-13, "{sum}, {err}");
        // parameter with hi = -3 and nonzero lo = 2^-60 does not falsely terminate
        let [up_nonpole, _] = residue::series_params(&[pow2(-60), 0.5], &[4.0], 0);
        assert_eq!(up_nonpole[1].hi, -3.0);
        assert_eq!(up_nonpole[1].lo, pow2(-60));
        assert!(first_zero_dd(&up_nonpole).is_infinite());
        // |a_m| >> 1 with v ~ a_m preserves +1 in (a_m - v) + 1
        let big = pow2(60);
        let [up_big, _] = residue::series_params(&[big, 0.5], &[big - 256.0], 0);
        assert_eq!(up_big[1], Dd::new(257.0));
        // genuine b_j = 1.0 with a_j = 2^-60 keeps exact 2^-60 in up (not rounded to 0)
        let tiny = pow2(-60);
        let [up_one, _] = residue::series_params(&[tiny, 0.5], &[1.0], 0);
        assert_eq!(up_one[0], Dd::new(tiny));
        assert_eq!(up_one[1], Dd::new(tiny));
        let (sum, err) = pfq_alg_series(&[tiny, 0.5], &[1.0], 0, -10.0, false);
        assert!(
            (sum - 1.0).abs() <= 4.0 * f64::EPSILON && err <= 4.0 * f64::EPSILON,
            "{sum}, {err}"
        );
        // exact termination at huge |z| does not overflow `TAU * |z|`
        let (sum, err) = pfq_alg_series(&[-1.0], &[1.0, 8.0], 0, -f64::MAX, false);
        assert_eq!(sum, 1.0);
        assert!(err.is_finite() && err <= 4.0 * f64::EPSILON, "{err}");
        // regrowth guard rounds conservatively on threshold pull = 1 + eps
        let mut scratch = BoundScratch::new(2, 1);
        let rho = scratch.regrowth_bound(
            &[Dd::new(3.0), Dd::new(3.0)],
            &[Dd::new(1.0)],
            next_up_nonneg(1.0 / 12.0),
            0.0,
        );
        assert_eq!(rho, f64::INFINITY);
        // six-factor threshold where `1/5` rounding in `powf` would exceed a small fixed margin
        let mut scratch = BoundScratch::new(6, 1);
        let rho = scratch.regrowth_bound(
            &[Dd::new(8.590_444_340_720_38e49); 6],
            &[Dd::new(1.0)],
            next_up_nonneg(1e-300),
            0.0,
        );
        assert_eq!(rho, f64::INFINITY);
        // repeated-pole shifted-factor case: early f64 rounding drops the combined bound below 1,
        // while preserving Dd low parts keeps it strictly above 1
        let am = 1.0 + pow2(-41);
        let a = [am, am, 8191.5 + pow2(-40)];
        let b = [8191.5, 1.0];
        let z = -5.002_442_002_406_743;
        let residues = Residues::of_class(&a, &b, 0, z);
        let [up_r, low_r] = residues.params();
        let mut scratch = BoundScratch::new(up_r.len(), low_r.len());
        let bound = scratch.ratio_bound(&up_r, &low_r, recip_up(z), 8191.0)
            * Residues::growth_bound(&up_r, &low_r, 8191.0);
        assert!(bound > 1.0, "{bound}");
        // terminating algebraic series with many parameters accounts for term-generation rounding
        let mut b_many = vec![2.5];
        b_many.resize(51, pow2(-53));
        let (sum, err) = pfq_alg_series(&[0.5], &b_many, 0, pow2(20), false);
        let expected = 305.041_624_171_279_54;
        assert!(
            err.is_finite() && (sum - expected).abs() <= err,
            "{sum} vs {expected} (err {err})"
        );
        // huge finite degree-1 term does not overflow error accumulation before scaling by eps
        let (sum_big, err_big) = pfq_alg_series(&[0.5, 1.0], &[2.5, pow2(1023)], 0, -2.0, false);
        assert!(sum_big.is_finite() && err_big.is_finite(), "{sum_big}, {err_big}");
        // subnormal scaled denominator in double-double rejects the lossy term
        let (_, err_sub) =
            pfq_alg_series_dd(&[0.5, 0.75], &[2.5, 1.5 + pow2(-52)], 0, -3.0 * f64::from_bits(1));
        assert!(err_sub.is_infinite(), "{err_sub}");
    }

    #[test]
    fn test_pfq_exp_cos_single_pass() {
        // kappa = 2: class 0 has a zero followed by a nonzero `k + kappa` lookahead term,
        // while class 1 decays and stops independently; nu = 0 gives phase theta = 200 rad
        // and cos(pi / 2) = 0; verify the underlying iterator is never polled past its first None
        let mut inner = [1.0, 0.5, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0].into_iter();
        let mut exhausted = false;
        let seq = from_fn(move || {
            assert!(!exhausted, "polled past first None");
            let next = inner.next();
            if next.is_none() {
                exhausted = true;
            }
            next
        });
        let ((sum, err), cos) =
            pfq_exp_cos(seq.map(|v| (Dd::new(v), 0.0)), 0.0, 2, Dd::new(100.0), Dd::ZERO, 0.0, 0.0);
        assert_eq!(cos, Dd::ZERO);
        // class 0 exact alternating sum is 1.0 - 0.0 + 0.25 = 1.25; must not stop at the 0.0
        let expected = 2.0 * 1.25 * 200.0_f64.cos() + 2.0 * 0.5 * 200.0_f64.sin();
        assert!((sum - expected).abs() <= err && err <= 1e-10, "{sum} vs {expected} (err {err})");
        // phase with integer and half-integer hi word and nonzero lo word preserves both
        // `-PI * t.lo * sin0` and `+PI * t.lo * cos0` corrections (omitting them would give 0)
        let root = Dd::new(100.0);
        let base = root * FRAC_2_PI_DD;
        let lo = base.lo;
        assert!(lo.abs() > 1e-18);
        let nu_int = Dd::new(2.0 * (64.0 - base.hi));
        let ((sum_int, err_int), _) = pfq_exp_cos(
            [0.0, 1.0, 0.0, 0.0].into_iter().map(|v| (Dd::new(v), 0.0)),
            0.0,
            2,
            root,
            nu_int,
            0.0,
            0.0,
        );
        let exp_int = 2.0 * PI * lo;
        assert!(
            (sum_int - exp_int).abs() <= err_int && (sum_int / exp_int - 1.0).abs() <= 1e-12,
            "{sum_int} vs {exp_int} (err {err_int})"
        );
        let nu_half = Dd::new(2.0 * (64.5 - base.hi));
        let ((sum_half, err_half), _) = pfq_exp_cos(
            [1.0, 0.0, 0.0, 0.0].into_iter().map(|v| (Dd::new(v), 0.0)),
            0.0,
            2,
            root,
            nu_half,
            0.0,
            0.0,
        );
        let exp_half = -2.0 * PI * lo;
        assert!(
            (sum_half - exp_half).abs() <= err_half && (sum_half / exp_half - 1.0).abs() <= 1e-12,
            "{sum_half} vs {exp_half} (err {err_half})"
        );
        // large cancelling `l[m+2] - r[m+2]` in `pfq_exp_terms` stays accurate in double-double
        let (b1, b2) = (-4000.6_f64, -9000.46_f64);
        let est = super::super::hyp2f2_estimate(b1 + 1.0, b2 + 1.0, b1, b2, 700.0).unwrap();
        let expected = 7.717_091_476_174_813e303;
        assert!(
            est.is_accurate() && (est.y - expected).abs() <= 1e-14 * expected,
            "{est:?} vs {expected}"
        );
        // Term::bound includes the low word and rounds outward
        assert!(
            (
                Dd {
                    hi: 1.0,
                    lo: pow2(-54)
                },
                0.0
            )
                .bound()
                > 1.0
        );
        // uncertain-zero lookahead preceded by an exact zero in the same residue class
        let (_, err_unc) = asymptotic_sum(
            [(Dd::ONE, 0.0), (Dd::ZERO, 0.0), (Dd::ZERO, 1.0), (Dd::ZERO, 0.0)].into_iter(),
            0.0,
            100.0,
            |_| 0.0,
        );
        assert!(err_unc >= 1.0, "{err_unc}");
        let ((_, err_cos_unc), _) = pfq_exp_cos(
            [
                (Dd::ONE, 0.0),
                (Dd::ZERO, 0.0),
                (Dd::ZERO, 0.0),
                (Dd::ZERO, 0.0),
                (Dd::ZERO, 1.0),
                (Dd::ZERO, 0.0),
                (Dd::ZERO, 0.0),
                (Dd::ZERO, 0.0),
            ]
            .into_iter(),
            0.0,
            2,
            root,
            nu_int,
            0.0,
            0.0,
        );
        assert!(err_cos_unc >= 1.0, "{err_cos_unc}");
        // `x_rel` propagates into `pfq_exp_cos` phase error budget for kappa = 3
        let root_276 = root_dd(pow2(276), 3);
        let ((_, err_phase), _) = pfq_exp_cos(
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                .into_iter()
                .map(|v| (Dd::new(v), 0.0)),
            0.0,
            3,
            root_276,
            Dd::new(-2.0),
            0.0,
            1e-29,
        );
        // 2 * |x * sin(pi / 3)| * 1e-29 ~ 0.257
        assert!(err_phase >= 0.25, "{err_phase}");
        // cross-term `s_err * c_err` with nonzero `nu_err` and an uncertain coefficient at
        // quadrature (`t = 5/6` puts class 1 at `cos(pi / 2) = 0`)
        let (_, sin3) = cos_sin_pi_div(3);
        let nu_quad = (Dd::new(5.0 / 6.0) - root * 3.0 * sin3 * (FRAC_2_PI_DD * 0.5)) * 3.0;
        let mut quad_terms = [(Dd::ZERO, 0.0); 9];
        quad_terms[1] = (Dd::ZERO, 1.0);
        let ((_, err_quad), _) =
            pfq_exp_cos(quad_terms.into_iter(), 0.0, 3, root, nu_quad, 1e-6, 0.0);
        // 2 * s_err (>= 1) * (PI * 1e-6 / 3) ~ 2.09e-6
        assert!(err_quad >= 2e-6, "{err_quad:e}");
        // pre-cancellation rounding in `nu - k + g_i` (DLMF 13.7.1: 1F1(2^-60; 2^60; 2^60))
        // is bounded by `t_err`: exact t_2 = 1 - 2^-61 + O(2^-120)
        let (a_tiny, b_huge) = (pow2(-60), pow2(60));
        let nu_k = Dd::new(a_tiny) - b_huge;
        let nu_k_err = 4.0 * Dd::EPSILON * (a_tiny + b_huge);
        let (t2, t2_err) =
            pfq_exp_terms(&[a_tiny], &[b_huge], nu_k, nu_k_err, Dd::new(b_huge), 0.0)
                .nth(2)
                .unwrap();
        let true_err = ((t2 - Dd::ONE) + pow2(-61)).hi.abs();
        assert!(t2_err >= true_err, "{t2_err} < {true_err} (t2={t2:?})");
        // reciprocal-power underflow: 2F2(1, 2; 2^64, 2^64; 2^540) has x^-2 < f64::MIN_POSITIVE
        // while d_3 * x^-2 ~ 4.85e-268 is normal
        let (b64, z540) = (pow2(64), pow2(540));
        let nu_2f2 = Dd::new(3.0) - Dd::new(b64) - Dd::new(b64);
        let nu_2f2_err = 6.0 * Dd::EPSILON * (3.0 + 2.0 * b64);
        let (t2_u, t2_u_err) =
            pfq_exp_terms(&[1.0, 2.0], &[b64, b64], nu_2f2, nu_2f2_err, Dd::new(z540), 0.0)
                .nth(2)
                .unwrap();
        let exact_t2_u = Dd {
            hi: 4.469_444_793_151_709_3e-249,
            lo: -1.453_734_526_361_732_8e-267,
        };
        let err_2f2 = (t2_u - exact_t2_u).hi.abs();
        assert!(
            err_2f2 <= t2_u_err && err_2f2 <= 1e-275,
            "{err_2f2:e} vs bound {t2_u_err:e} (t2={t2_u:?})"
        );
        // root_dd error propagation: 0F2(; 1, 2; 2^900) has exact t_1 = -2 / (9 * 2^300)
        let z900 = pow2(900);
        let x900 = root_dd(z900, 3) * 3.0;
        let nu_0f2 = Dd::new(-2.0);
        let nu_0f2_err = 4.0 * Dd::EPSILON * 4.0;
        let (t1_r, t1_r_err) =
            pfq_exp_terms(&[], &[1.0, 2.0], nu_0f2, nu_0f2_err, x900, pfq_exp_x_rel(3))
                .nth(1)
                .unwrap();
        let exact_t1_r = Dd {
            hi: -1.090_909_658_955_050_3e-91,
            lo: -6.055_765_105_789_37e-108,
        };
        let err_0f2 = (t1_r - exact_t1_r).hi.abs();
        assert!(err_0f2 <= t1_r_err, "{err_0f2:e} > bound {t1_r_err:e} (t1={t1_r:?})");
        // real exponential series stops at `f64::EPSILON` without extra double-double tail terms
        let mut drawn = 0usize;
        let counted = pfq_exp_terms(&[-0.5], &[0.5], Dd::new(-1.0), 0.0, Dd::new(100.0), 0.0)
            .inspect(|_| drawn += 1);
        let (sum_f64, err_f64) = asymptotic_sum_tol(counted, 0.0, 100.0, f64::EPSILON, |_| 0.0);
        let (sum_dd, _) = asymptotic_sum(
            pfq_exp_terms(&[-0.5], &[0.5], Dd::new(-1.0), 0.0, Dd::new(100.0), 0.0),
            0.0,
            100.0,
            |_| 0.0,
        );
        assert!(drawn <= 20, "took {drawn} terms");
        assert!(
            (sum_f64 - sum_dd).hi.abs() <= err_f64 && err_f64 <= 4.0 * f64::EPSILON * sum_dd.hi,
            "{sum_f64:?} vs {sum_dd:?} (err {err_f64:e})"
        );
    }

    #[test]
    fn test_convergent_sum_edge_cases() {
        // reaches PFQ_MAXITER with rho < 1: returns finite tail bound rather than infinity
        let r = 0.999_5_f64;
        let geom = successors(Some(1.0_f64), |&t| Some(t * r));
        let (sum, err) = convergent_sum(geom, 0.0, f64::EPSILON, 0.0, |_| r);
        let exact = 1.0 / (1.0 - r);
        assert!(err.is_finite() && (sum - exact).abs() <= err, "{sum} vs {exact} ({err})");
        // early None (e.g. from overflow/NaN inside pfq_terms_of) returns infinite error
        let (_, err) = convergent_sum([1.0, 0.5].into_iter(), 0.0, f64::EPSILON, 0.0, |_| 0.5);
        assert_eq!(err, f64::INFINITY);
        // underflow followed by regrowth: yields the first negligible bound, then stops -> inf err
        let make_regrow =
            || pfq_terms_of(f64::INFINITY, 2e155, |n, t| if n == 0.0 { 1e-320 } else { t * 2e155 });
        assert_eq!(make_regrow().count(), 1);
        let (_, err) =
            convergent_sum(once(1.0_f64).chain(make_regrow()), 10.0, f64::EPSILON, 0.0, |_| 0.5);
        assert_eq!(err, f64::INFINITY);
        // sum overflow on finite terms returns infinite error
        let (_, err) =
            convergent_sum([f64::MAX, f64::MAX].into_iter(), 0.0, f64::EPSILON, 0.0, |_| 0.5);
        assert_eq!(err, f64::INFINITY);
        // last-term generation error inflates the geometric tail when stopping at PFQ_MAXITER with
        // rho very close to 1
        let rho_near1 = 1.0 - 0.5 * f64::EPSILON;
        let (_, err_tail) =
            convergent_sum(core::iter::repeat(1.0_f64), 0.0, f64::EPSILON, 1.0, |_| rho_near1);
        let n_last = f64::from(PFQ_MAXITER - 1);
        let uninflated = rho_near1 / (1.0 - rho_near1)
            + n_last * f64::EPSILON * f64::from(PFQ_MAXITER)
            + 0.5 * n_last * f64::from(PFQ_MAXITER) * f64::EPSILON;
        assert!(
            err_tail >= uninflated + 0.5 * n_last * f64::EPSILON * rho_near1 / (1.0 - rho_near1),
            "err_tail={err_tail:e} vs uninflated={uninflated:e}"
        );
        // double-double tail rounds both the last term magnitude (1 + 2^-54) and tail upward
        let dd_term = Dd {
            hi: 1.0,
            lo: pow2(-54),
        };
        let (_, err_tail_dd) =
            convergent_sum(core::iter::repeat(dd_term), 0.0, Dd::EPSILON, 1.0, |_| rho_near1);
        assert!(
            err_tail_dd >= 9_007_199_254_740_992.0,
            "err_tail_dd={err_tail_dd:.1} understated exact tail ~9007199254740991.5"
        );
        // many-parameter 31F30 reciprocal series accounts for per-step factor-product roundings
        let a31: Vec<f64> = (0..31).map(|i| 0.8 + 0.005 * f64::from(i)).collect();
        let b30: Vec<f64> = (0..30).map(|i| 0.7 + 0.005 * f64::from(i)).collect();
        let (s_f64, err_f64) = pfq_alg_series(&a31, &b30, 0, -25.0, false);
        let (s_dd, err_dd) = pfq_alg_series_dd(&a31, &b30, 0, -25.0);
        let actual = (s_f64 - s_dd.hi).abs();
        assert!(
            err_dd <= 1e-25 && actual <= err_f64 && err_f64 >= 20.0 * f64::EPSILON,
            "actual={actual:e} reported={err_f64:e} dd_err={err_dd:e}"
        );
    }

    #[test]
    fn test_pfq_continuation_edge_cases() {
        // zero inner polynomial (a_0 = 0.5 gives 1 - 1 = 0) must not abort the full sum
        let a = [0.5, 1.5 - pow2(-40), 1.75];
        let b = [2.5, 1.0];
        let z = -pow2(40);
        let e = pfq_continuation(&a, &b, z).unwrap();
        let expected = 7.959_621_827_055_401e-7;
        assert!(
            e.is_accurate() && (e.y - expected).abs() <= 1e-13 * expected.abs(),
            "{e:?} vs {expected}"
        );
        // mixed exact-integer class {2.75, 3.75} + cancelling singletons {0.5, 0.5 + 2^-32}
        let a = [0.5, 0.5 + pow2(-32), 2.75, 3.75];
        let b = [2.25, 3.5, 4.125];
        let e = pfq_continuation(&a, &b, -1000.0).unwrap();
        let expected = 0.191_031_737_123_466_5;
        assert!(
            e.is_accurate() && (e.y - expected).abs() <= 1e-13 * expected.abs(),
            "{e:?} vs {expected}"
        );
        // direct continuation just past z = -1 (bypasses Weniger routing)
        let a = [0.5, 1.25, 2.75];
        let b = [2.25, 3.5];
        let e = pfq_continuation(&a, &b, -1.01).unwrap();
        let expected = 0.842_648_061_384_124_5;
        assert!(
            e.is_accurate() && (e.y - expected).abs() <= 1e-13 * expected.abs(),
            "{e:?} vs {expected}"
        );
    }
}
