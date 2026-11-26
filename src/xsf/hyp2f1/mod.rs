//! The complex `hyp2f1` is translated into pure Rust from `xsf/hyp2f1.h` (scipy/xsf@f2a68ec).
//!
//! ```text
//! Implementation of Gauss's hypergeometric function for complex values.
//!
//! This implementation is based on the Fortran implementation by Shanjie Zhang and
//! Jianming Jin included in specfun.f [1]_.  Computation of Gauss's hypergeometric
//! function involves handling a patchwork of special cases. By default the Zhang and
//! Jin implementation has been followed as closely as possible except for situations where
//! an improvement was obvious. We've attempted to document the reasons behind decisions
//! made by Zhang and Jin and to document the reasons for deviating from their implementation
//! when this has been done. References to the NIST Digital Library of Mathematical
//! Functions [2]_ have been added where they are appropriate. The review paper by
//! Pearson et al [3]_ is an excellent resource for best practices for numerical
//! computation of hypergeometric functions. We have followed this review paper
//! when making improvements to and correcting defects in Zhang and Jin's
//! implementation. When Pearson et al propose several competing alternatives for a
//! given case, we've used our best judgment to decide on the method to use.
//!
//! Author: Albert Steppi
//!
//! Distributed under the same license as Scipy.
//!
//! References
//! ----------
//! .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996
//! .. [2] NIST Digital Library of Mathematical Functions. http://dlmf.nist.gov/,
//!        Release 1.1.1 of 2021-03-15. F. W. J. Olver, A. B. Olde Daalhuis,
//!        D. W. Lozier, B. I. Schneider, R. F. Boisvert, C. W. Clark, B. R. Miller,
//!        B. V. Saunders, H. S. Cohl, and M. A. McClain, eds.
//! .. [3] Pearson, J.W., Olver, S. & Porter, M.A.
//!        "Numerical methods for the computation of the confluent and Gauss
//!        hypergeometric functions."
//!        Numer Algor 74, 821-866 (2017). https://doi.org/10.1007/s11075-016-0173-0
//! .. [4] Raimundas Vidunas, "Degenerate Gauss Hypergeometric Functions",
//!        Kyushu Journal of Mathematics, 2007, Volume 61, Issue 1, Pages 109-135,
//! .. [5] López, J.L., Temme, N.M. New series expansions of the Gauss hypergeometric
//!        function. Adv Comput Math 39, 349-365 (2013).
//!        https://doi.org/10.1007/s10444-012-9283-y
//! ```
//!
//! The C++ generator classes are functions that return `FnMut` closures here. `set_error` is a
//! no-op in our build of xsf, so the calls to it are omitted.

mod beta_endpoints;
mod polynomial;
mod power;

use core::f64::consts::{FRAC_1_PI, PI, TAU};

use crate::ext::{HYP2F1_RECUR_TOL, Scaled, fold_exp, hyp2f1_poly_dd, ln_dd, mul_gamma_ratio};
use crate::utils::{is_nonpos_int, neg_one_pow, next_down_nonneg, next_up_nonneg, pow2};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, two_sum};
use crate::xsf::cephes::hyp2f1::{hyp2f1_at_1, ln_max_growth, max_term_index};
use crate::xsf::cephes::lanczos::LANCZOS_G;
use crate::xsf::cephes::{gamma, gammasgn, lanczos_sum_expg_scaled, lgam, poch, rgamma};
use crate::xsf::complex::cpow_integral;
use crate::xsf::tools::{self, series_eval_fixed_length};
use crate::{cospi, digamma, sinpi};

use power::{ScaledC64, ldexp_sat};

const HYP2F1_EPS: f64 = 1e-15;
/* The original implementation in SciPy from Zhang and Jin used 1500 for the
 * maximum number of series iterations in some cases and 500 in others.
 * Through the empirical results on the test cases in
 * scipy/special/_precompute/hyp2f1_data.py, it was determined that these values
 * can lead to early termination of series which would have eventually converged
 * at a reasonable level of accuracy. We've bumped the iteration limit to 3000,
 * and may adjust it again based on further analysis. */
const HYP2F1_MAXITER: u64 = 3000;
// Unlike in xsf, bound the work of the optional Lopez-Temme retry in `hyp2f1_inverse` by the
// smaller limit of Zhang and Jin.
const LOPEZ_TEMME_RETRY_MAXITER: u64 = 500;

/// `tools::series_eval` with the tolerance and term limit of hyp2f1
fn series_eval(g: impl FnMut() -> C64) -> C64 {
    tools::series_eval(g, C64::ZERO, HYP2F1_EPS, HYP2F1_MAXITER)
}

fn four_gammas_lanczos(u: f64, v: f64, w: f64, x: f64) -> f64 {
    /* Compute ratio of gamma functions using lanczos approximation.
     *
     * Computes gamma(u)*gamma(v)/(gamma(w)*gamma(x))
     *
     * It is assumed that x = u + v - w, but it is left to the user to
     * ensure this.
     *
     * The lanczos approximation takes the form
     *
     * gamma(x) = factor(x) * lanczos_sum_expg_scaled(x)
     *
     * where factor(x) = ((x + lanczos_g - 0.5)/e)**(x - 0.5).
     *
     * The formula above is only valid for x >= 0.5, but can be extended to
     * x < 0.5 with the reflection principle.
     *
     * Using the lanczos approximation when computing this ratio of gamma functions
     * allows factors to be combined analytically to avoid underflow and overflow
     * and produce a more accurate result. The condition x = u + v - w makes it
     * possible to cancel the factors in the expression
     *
     * factor(u) * factor(v) / (factor(w) * factor(x))
     *
     * by taking one factor and absorbing it into the others. Currently, this
     * implementation takes the factor corresponding to the argument with largest
     * absolute value and absorbs it into the others.
     *
     * Since this is only called internally by four_gammas. It is assumed that
     * |u| >= |v| and |w| >= |x|.
     */

    /* The below implementation may incorrectly return finite results
     * at poles of the gamma function. Handle these cases explicitly. */
    if is_nonpos_int(u) || is_nonpos_int(v) {
        /* Return nan if numerator has pole. Diverges to +- infinity
         * depending on direction so value is undefined. */
        return f64::NAN;
    }
    if is_nonpos_int(w) || is_nonpos_int(x) {
        // Return 0 if denominator has pole but not numerator.
        return 0.0;
    }

    let (u_num, u_den, ugh, u_prime) = lanczos_parts(u);
    let (v_num, v_den, vgh, v_prime) = lanczos_parts(v);
    let (w_num, w_den, wgh, w_prime) = lanczos_parts(w);
    let (x_num, x_den, xgh, x_prime) = lanczos_parts(x);
    let mut result = u_num / u_den * v_num / v_den * w_den / w_num * x_den / x_num;

    if u.abs() >= w.abs() {
        // u has greatest absolute value. Absorb ugh into the others.
        result *= pow_ratio(vgh, ugh, v_prime - u_prime, v, v > 100.0);
        result *= pow_ratio(ugh, wgh, u_prime - w_prime, w, u > 100.0);
        result *= pow_ratio(ugh, xgh, u_prime - x_prime, x, u > 100.0);
    } else {
        // w has greatest absolute value. Absorb wgh into the others.
        result *= pow_ratio(ugh, wgh, u_prime - w_prime, u, u > 100.0);
        result *= pow_ratio(vgh, wgh, v_prime - w_prime, v, v > 100.0);
        result *= pow_ratio(wgh, xgh, w_prime - x_prime, x, x > 100.0);
    }
    // This exhausts all cases because we assume |u| >= |v| and |w| >= |x|.

    result
}

/// `(num, den, gh, t')`, with `num / den` the Lanczos sum factor of `gamma(t)`, `t'` either `t` or
/// `1 - t` if `t < 0.5` (reflection), and `gh = t' + g - 1/2`
///
/// `num` or `den` is 1, so that `r * num / den` rounds as a single product or quotient.
fn lanczos_parts(t: f64) -> (f64, f64, f64, f64) {
    if t >= 0.5 {
        (lanczos_sum_expg_scaled(t), 1.0, t + LANCZOS_G - 0.5, t)
    } else {
        let den = lanczos_sum_expg_scaled(1.0 - t) * (PI * t).sin() * FRAC_1_PI;
        (1.0, den, 0.5 - t + LANCZOS_G, 1.0 - t)
    }
}

/// `(num_gh / den_gh)^(e - 1/2)`, where `d` is the difference of the primed arguments of `num_gh`
/// and `den_gh` (before the Lanczos shift), via `ln_1p(d / den_gh)` if `big` and the base is close
/// to 1
fn pow_ratio(num_gh: f64, den_gh: f64, d: f64, e: f64, big: bool) -> f64 {
    if (d * (e - 0.5)).abs() < 100.0 * den_gh && big {
        /* Special case where base is close to 1. Condition taken from
         * Boost's beta function implementation. */
        ((e - 0.5) * (d / den_gh).ln_1p()).exp()
    } else {
        (num_gh / den_gh).powf(e - 0.5)
    }
}

/// Whether `four_gammas` tries the direct ratio of gamma functions at these arguments
fn in_direct_gamma_range(args: impl IntoIterator<Item = f64>) -> bool {
    /* Direct ratio tends to be more accurate for arguments in this range. Range
     * chosen empirically based on the relevant benchmarks in
     * scipy/special/_precompute/hyp2f1_data.py */
    args.into_iter().all(|x| x.abs() <= 100.0)
}

fn four_gammas(mut u: f64, mut v: f64, mut w: f64, mut x: f64) -> f64 {
    // Without loss of generality, ensure |u| >= |v| and |w| >= |x|.
    if v.abs() > u.abs() {
        core::mem::swap(&mut u, &mut v);
    }
    if x.abs() > w.abs() {
        core::mem::swap(&mut x, &mut w);
    }
    if in_direct_gamma_range([u, v, w, x]) {
        let result = gamma(u) * gamma(v) * (rgamma(w) * rgamma(x));
        if result.is_finite() && result != 0.0 {
            return result;
        }
    }
    let mut result = four_gammas_lanczos(u, v, w, x);
    if result.is_finite() && result != 0.0 {
        return result;
    }
    // If overflow or underflow, try again with logs.
    result = (lgam(v) - lgam(x) + lgam(u) - lgam(w)).exp();
    result *= gammasgn(u) * gammasgn(w) * gammasgn(v) * gammasgn(x);
    result
}

fn hypergeometric_series_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> C64 {
    /* Maclaurin series for hyp2f1.
     *
     * Series is convergent for |z| < 1 but is only practical for numerical
     * computation when |z| < 0.9.
     */
    let mut term = C64::ONE;
    let mut k: u64 = 0;
    move || {
        let output = term;
        let k_ = k as f64;
        term = term * (a + k_) * (b + k_) / ((k + 1) as f64 * (c + k_)) * z;
        k += 1;
        output
    }
}

/// `series_eval` of the Maclaurin series for |z| < 1, where `rounding` is the relative error of
/// the later terms from rounding `a` and `b`
//
// Unlike in xsf, a small term only ends the sum if the later ones can add up to the tolerance,
// e.g. after a tiny first term from a tiny `a`. They are then summed as a series of their own, as
// they can be far below the sum, and kept only if their estimated error (each of the n terms errs
// by about n eps) implies an accurate result, and their sum exceeds that error plus the rounding
// of the earlier terms, so that the result is closer than without them.
fn hypergeometric_series(a: f64, b: f64, c: f64, z: C64, rounding: f64) -> C64 {
    let mut series_generator = hypergeometric_series_generator(a, b, c, z);
    let (mut n, mut last) = (0.0, C64::ZERO);
    let head = series_eval_cancellation(|| {
        (n, last) = (n + 1.0, series_generator());
        last
    });
    // Unlike in xsf, continue the series analytically if it cancels.
    if head.is_continued() {
        if let Some(value) = hyp2f1_continuation(a, b, c, z, head.value.into()) {
            return value;
        }
    }
    let x = z.norm();
    let growth = ln_max_growth(a, b, c, x, n - 1.0);
    // A subnormal term would pass its rounding on to the larger later ones.
    if !(growth > 0.0 && last.norm() >= f64::MIN_POSITIVE) {
        return head.value;
    }
    let peak = (last.norm().ln() + growth).exp();
    // As a heuristic, the later terms have the magnitudes of the real ones at |z|, which grow up to
    // the largest one, at m, and then fall like a Gaussian, from the slope of the logarithm of the
    // term ratio there, but no faster than a geometric series in the limit |z| of that ratio (also
    // if the slope is negative, as `max` drops the NaN).
    let m = max_term_index(a, b, c, x);
    let slope = 1.0 / (c + m) + 1.0 / (m + 1.0) - 1.0 / (a + m) - 1.0 / (b + m);
    let width = (TAU / slope).sqrt().max(1.0 / (1.0 - x));
    let mass = peak * (m - n + 1.0 + width);
    let negligible = head.value.norm() * HYP2F1_EPS;
    if !(mass >= negligible) {
        return head.value;
    }
    // Their sum, as they turn by the angle of z (the Fourier transform of the Gaussian), unless
    // negligible, and the number of them until one is below `HYP2F1_EPS` times that sum
    let angle = z.im.atan2(z.re);
    let sum = (peak * width * (-angle * angle / (2.0 * slope)).exp()).max(negligible);
    let length = m - n + (2.0 * (peak / (HYP2F1_EPS * sum)).ln() / slope).sqrt();
    // Unlike in xsf, NaN if z is real, as then (with z.re >= 0) the later terms have the same
    // sign, so that `head` is wrong.
    let fallback = if z.is_real() && peak > head.value.norm() {
        C64::NAN
    } else {
        head.value
    };
    // Give up if the sum would take too long, or be below the rounding of the terms.
    if !(length <= HYP2F1_MAXITER as f64 && mass * f64::EPSILON <= sum) {
        return fallback;
    }
    let n_head = n as u64;
    let tail = series_eval_cancellation(|| {
        n += 1.0;
        series_generator()
    });
    if tail.value.is_nan() {
        return fallback;
    }
    // Add the first term 1 last, so that the rounding of the rest of `head` isn't lost.
    let mut head_terms = hypergeometric_series_generator(a, b, c, z);
    let value = head_terms() + series_eval_fixed_length(head_terms, tail.value, n_head - 1);
    let error = (n * f64::EPSILON + rounding) * tail.cancellation * tail.value.norm();
    // The first term 1 is exact.
    let head_rounding = f64::EPSILON * (head.cancellation * head.value.norm() - 1.0);
    if error <= HYP2F1_RECUR_TOL * value.norm() && error + head_rounding <= tail.value.norm() {
        value
    } else {
        head.value
    }
}

fn hyp2f1_transform1_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> C64 {
    /* 1 -z transformation of standard series.*/
    let factor1 = c64(four_gammas(c, c - a - b, c - a, c - b), 0.0);
    let factor2 = (1.0 - z).powf(c - a - b) * four_gammas(c, a + b - c, a, b);
    let mut generator1 = hypergeometric_series_generator(a, b, a + b - c + 1.0, 1.0 - z);
    let mut generator2 = hypergeometric_series_generator(c - a, c - b, c - a - b + 1.0, 1.0 - z);
    move || factor1 * generator1() + factor2 * generator2()
}

fn hyp2f1_transform1_limit_series_generator(a: f64, b: f64, m: f64, z: C64) -> impl FnMut() -> C64 {
    /* 1 - z transform in limit as c - a - b approaches an integer m. */
    let mut d1 = digamma(a);
    let mut d2 = digamma(b);
    let mut d3 = digamma(1.0 + m);
    let mut d4 = digamma(1.0);
    let log_1_z = (1.0 - z).ln();
    let mut factor = c64(rgamma(m + 1.0), 0.0);
    let mut k: i32 = 0;
    move || {
        let term = (d1 + d2 - d3 - d4 + log_1_z) * factor;
        let k_ = f64::from(k);
        // Use digamma(x + 1) = digamma(x) + 1/x
        d1 += 1.0 / (a + k_); // d1 = digamma(a + k)
        d2 += 1.0 / (b + k_); // d2 = digamma(b + k)
        d3 += 1.0 / (1.0 + m + k_); // d3 = digamma(1 + m + k)
        d4 += 1.0 / (1.0 + k_); // d4 = digamma(1 + k)
        factor *= (1.0 - z) * ((a + k_) * (b + k_) / ((k_ + 1.0) * (m + k_ + 1.0)));
        k += 1;
        term
    }
}

fn hyp2f1_transform2_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> (C64, C64) {
    /* 1/z transformation of standard series.*/
    let factor1 = (-z).powf(-a) * four_gammas(c, b - a, b, c - a);
    let factor2 = (-z).powf(-b) * four_gammas(c, a - b, a, c - b);
    let mut generator1 = hypergeometric_series_generator(a, a - c + 1.0, a - b + 1.0, 1.0 / z);
    let mut generator2 = hypergeometric_series_generator(b, b - c + 1.0, b - a + 1.0, 1.0 / z);
    // Unlike in xsf, the terms of both series are returned apart, for `series_eval_pairs`.
    move || (factor1 * generator1(), factor2 * generator2())
}

/// `hyp2f1_transform2_generator(c - a, c - b, c, z)`, but without `c - (c - a)` losing `a` if tiny
fn hyp2f1_transform2_euler_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> (C64, C64) {
    let factor1 = (-z).powf(a - c) * four_gammas(c, a - b, c - b, a);
    let factor2 = (-z).powf(b - c) * four_gammas(c, b - a, c - a, b);
    let mut generator1 = hypergeometric_series_generator(c - a, 1.0 - a, b - a + 1.0, 1.0 / z);
    let mut generator2 = hypergeometric_series_generator(c - b, 1.0 - b, a - b + 1.0, 1.0 / z);
    // (the terms of both series apart, for `series_eval_pairs`)
    move || (factor1 * generator1(), factor2 * generator2())
}

/// The estimated cancellation `abs_sum / |sum|`, or infinity if NaN (e.g. if `sum` is NaN)
fn cancellation(abs_sum: f64, sum: C64) -> f64 {
    let ratio = abs_sum / sum.norm();
    if ratio.is_nan() { f64::INFINITY } else { ratio }
}

/// A value with its estimated absolute error
#[derive(Clone, Copy)]
struct Estimate {
    value: C64,
    error: f64,
}

// A value whose only error is its final rounding to f64
impl From<C64> for Estimate {
    fn from(value: C64) -> Self {
        Self { value, error: 0.0 }
    }
}

impl Estimate {
    // Whether the estimate is finite with an estimated error below its magnitude.
    fn resolved(&self) -> bool {
        self.value.is_finite() && self.error < self.value.linf_norm()
    }

    // Whether the estimate leaves at least half the significand.
    fn accurate(&self) -> bool {
        self.value.is_finite() && self.error <= pow2(-26) * self.value.linf_norm()
    }
}

/// A value with its estimated cancellation
#[derive(Clone, Copy)]
struct Cancelled {
    value: C64,
    cancellation: f64,
}

impl Cancelled {
    /// Whether the estimated relative roundoff exceeds `HYP2F1_EPS`
    fn is_lossy(self) -> bool {
        self.cancellation * f64::EPSILON > HYP2F1_EPS
    }

    /// Whether the estimated relative roundoff exceeds `CONTINUED_ROUNDOFF`
    fn is_continued(self) -> bool {
        self.cancellation * f64::EPSILON > CONTINUED_ROUNDOFF
    }

    /// Whether the estimated relative roundoff is at most `USABLE_ERROR`
    fn is_usable(self) -> bool {
        self.cancellation * f64::EPSILON <= USABLE_ERROR
    }
}

/// The estimated relative roundoff of a series above which to try `hyp2f1_continuation`; below
/// it the continuation gains little accuracy at much more cost
const CONTINUED_ROUNDOFF: f64 = 100.0 * HYP2F1_EPS;

/// The largest estimated relative error of a usable sum
const USABLE_ERROR: f64 = 0.125;

/// `series_eval` of the generator, and its estimated cancellation sum(|term|) / |sum|
fn series_eval_cancellation(g: impl FnMut() -> C64) -> Cancelled {
    series_eval_cancellation_to(g, HYP2F1_EPS)
}

/// `series_eval_cancellation` to the relative tolerance `tol`
fn series_eval_cancellation_to(mut g: impl FnMut() -> C64, tol: f64) -> Cancelled {
    let mut abs_sum = 0.0;
    let g = || {
        let term = g();
        abs_sum += term.norm();
        term
    };
    let result = tools::series_eval(g, C64::ZERO, tol, HYP2F1_MAXITER);
    Cancelled {
        value: result,
        cancellation: cancellation(abs_sum, result),
    }
}

/// `series_eval_cancellation` of the sums of the pairs of terms, but with both their norms in the
/// cancellation if `pairs`
fn series_eval_pairs(mut g: impl FnMut() -> (C64, C64), pairs: bool) -> Cancelled {
    let mut abs_sum = 0.0;
    let sum = series_eval(|| {
        let (first, second) = g();
        let term = first + second;
        abs_sum += if pairs {
            first.norm() + second.norm()
        } else {
            term.norm()
        };
        term
    });
    Cancelled {
        value: sum,
        cancellation: cancellation(abs_sum, sum),
    }
}

/// `series_eval_cancellation`, and the estimated relative size of the omitted tail, a geometric
/// series in the ratio of the last two terms (infinite unless that ratio is below one)
fn series_eval_tail(mut g: impl FnMut() -> C64) -> (Cancelled, f64) {
    let (mut previous, mut last) = (0.0, 0.0);
    let series = series_eval_cancellation(|| {
        let term = g();
        (previous, last) = (last, term.norm());
        term
    });
    let ratio = last / previous;
    let tail = if ratio < 1.0 {
        last * ratio / (1.0 - ratio)
    } else {
        f64::INFINITY
    };
    (series, tail / series.value.norm())
}

/// `series_eval_cancellation`, but converged only after three successive small terms
//
// One term of the limit or paired series can vanish with its logarithmic coefficient.
fn series_eval_settled(mut g: impl FnMut() -> C64) -> Cancelled {
    let (mut sum, mut abs_sum, mut small) = (C64::ZERO, 0.0, 0);
    for _ in 0..HYP2F1_MAXITER {
        let term = g();
        sum += term;
        abs_sum += term.norm();
        // (a NaN sum never settles)
        if sum.norm().is_nan() {
            break;
        }
        small = if term.norm() < sum.norm() * HYP2F1_EPS {
            small + 1
        } else {
            0
        };
        if small == 3 {
            return Cancelled {
                value: sum,
                cancellation: cancellation(abs_sum, sum),
            };
        }
    }
    Cancelled {
        value: C64::NAN,
        cancellation: f64::INFINITY,
    }
}

/// `series_eval_cancellation` of a series whose ratio of terms `n + 1` and `n` divides by `c + n`
///
/// A `retry` gives the term magnitudes that make it NaN.
//
// Unlike in xsf, stop only after two successive terms are small relative to the sum and the
// nearest denominator ahead, since the terms of one parity can be small and a small denominator can
// grow them again. A retry can replace a NaN, so it waits until all denominators ahead are positive
// and the last pair plus its geometric tail is small.
fn series_eval_lopez_temme(mut g: impl FnMut() -> C64, c: f64, retry: Option<f64>) -> Cancelled {
    let (mut sum, mut abs_sum) = (C64::ZERO, 0.0);
    // The previous term magnitude, and the magnitudes of the pairs of terms ending 2 and 1 before
    let (mut previous, mut older_pair, mut last_pair) = (0.0, 0.0, 0.0);
    for n in 0..HYP2F1_MAXITER {
        let term = g();
        sum += term;
        let magnitude = term.norm();
        abs_sum += magnitude;
        if retry.is_some_and(|limit| !(abs_sum < limit)) {
            break;
        }
        let ahead = c + n as f64;
        let nearest = match retry {
            _ if ahead > 0.0 => ahead,
            None => (c - c.round()).abs(),
            Some(_) => 0.0,
        };
        let tol = sum.norm() * HYP2F1_EPS * nearest.min(1.0);
        // Skip pairs with the first term, since the terms after it can be tiny and still grow.
        let pair = previous + magnitude;
        let settled = match retry {
            None => previous < tol && magnitude < tol,
            Some(_) => n > 3 && 0.0 < pair && pair < tol * (1.0 - pair / older_pair),
        };
        if settled {
            return Cancelled {
                value: sum,
                cancellation: cancellation(abs_sum, sum),
            };
        }
        (previous, older_pair, last_pair) = (magnitude, last_pair, pair);
    }
    Cancelled {
        value: C64::NAN,
        cancellation: f64::INFINITY,
    }
}

const HYP2F1_MAX_FIXED_TERMS: u64 = 1 << 20;

/// `series_eval_fixed_length`, or NaN if that would sum more than `HYP2F1_MAX_FIXED_TERMS` terms
//
// Unlike in xsf (scipy/xsf#388), cap num_terms.
fn series_eval_finite(g: impl FnMut() -> C64, num_terms: u64) -> C64 {
    if num_terms > HYP2F1_MAX_FIXED_TERMS {
        return C64::NAN;
    }
    series_eval_fixed_length(g, C64::ZERO, num_terms)
}

fn hyp2f1_transform2_limit_series_generator(
    a: f64,
    b: f64,
    c: f64,
    m: f64,
    factor: f64,
    z: C64,
) -> impl FnMut() -> C64 {
    /* 1/z transform in limit as a - b approaches a non-negative integer m. (Can swap a and b to
     * handle the m a negative integer case. */
    let mut d1 = digamma(1.0);
    let mut d2 = digamma(1.0 + m);
    let mut d3 = digamma(a);
    let mut d4 = digamma(c - a);
    let log_neg_z = (-z).ln();
    // Unlike in xsf, start from the caller's factor, which includes the gammas of the prefactor.
    let mut factor = c64(factor, 0.0);
    let mut k: u64 = 0;
    move || {
        let term = (d1 + d2 - d3 - d4 + log_neg_z) * factor;
        let k_ = k as f64;
        // Use digamma(x + 1) = digamma(x) + 1/x
        d1 += 1.0 / (1.0 + k_); // d1 = digamma(1 + k)
        d2 += 1.0 / (1.0 + m + k_); // d2 = digamma(1 + m + k)
        let ak = a + k_;
        // Unlike in xsf, restart d3 = digamma(a + k) when a + k is within 1/2 of its pole at 0, as
        // the recurrence keeps the rounding error of 1/(a + k) in every later term.
        d3 = if (-0.5..0.5).contains(&ak) {
            digamma(ak + 1.0)
        } else {
            d3 + 1.0 / ak
        };
        let s = c - a - k_ - 1.0;
        d4 -= 1.0 / s; // d4 = digamma(c - a - k)
        // Unlike in xsf, take 1 - c + b + m as 1 - c + a, so that the zero of the factor cancels
        // the pole of d4 exactly.
        let f = -(b + m + k_) * s / ((k + 1) as f64 * (m + k_ + 1.0));
        factor *= f / z;
        k += 1;
        term
    }
}

fn hyp2f1_transform2_limit_series_cminus_a_int_generator(
    a: f64,
    b: f64,
    c: f64,
    m: f64,
    n: f64,
    scaled_factor: ScaledC64,
    z: C64,
) -> impl FnMut() -> C64 {
    /* 1/z transform in limit as a - b approaches a non-negative integer m, and c - a approaches
     * a positive integer n. */
    let mut d1 = digamma(1.0);
    let mut d2 = digamma(1.0 + m);
    let mut d3 = digamma(a);
    let mut d4 = digamma(n);
    let log_neg_z = (-z).ln();
    // Unlike in xsf, start from the caller's factor, which includes the gammas of the prefactor,
    // relative to its scale.
    let mut factor = scaled_factor.mantissa;
    let mut k: u64 = 0;
    move || {
        let k_ = k as f64;
        if k_ < n {
            let term = (d1 + d2 - d3 - d4 + log_neg_z) * factor;
            // Use digamma(x + 1) = digamma(x) + 1/x
            d1 += 1.0 / (1.0 + k_); // d1 = digamma(1 + k)
            d2 += 1.0 / (1.0 + m + k_); // d2 = digamma(1 + m + k)
            d3 += 1.0 / (a + k_); // d3 = digamma(a + k)
            d4 -= 1.0 / (n - k_ - 1.0); // d4 = digamma(c - a - k)
            let f = (b + m + k_) * (1.0 - c + b + m + k_) / ((k + 1) as f64 * (m + k_ + 1.0));
            factor *= f / z;
            k += 1;
            return term;
        }
        if k_ == n {
            /* When c - a approaches a positive integer and k_ >= c - a = n then
             * poch(1 - c + b + m + k) = poch(1 - c + a + k) = approaches zero and
             * digamma(c - a - k) approaches a pole. However we can use the limit
             * digamma(-n + epsilon) / gamma(-n + epsilon) -> (-1)**(n + 1) * (n+1)! as epsilon -> 0
             * to continue the series.
             *
             * poch(1 - c + b, m + k) = gamma(1 - c + b + m + k)/gamma(1 - c + b)
             *
             * If a - b is an integer and c - a is an integer, then a and b must both be integers, so assume
             * a and b are integers and take the limit as c approaches an integer.
             *
             * gamma(1 - c + epsilon + a + k)/gamma(1 - c - epsilon + b) =
             * (gamma(c + epsilon - b) / gamma(c + epsilon - a - k)) *
             * (sin(pi * (c + epsilon - b)) / sin(pi * (c + epsilon - a - k))) (reflection principle)
             *
             * In the limit as epsilon goes to zero, the ratio of sines will approach
             * (-1)**(a - b + k) = (-1)**(m + k)
             *
             * We may then replace
             *
             * poch(1 - c - epsilon + b, m + k)*digamma(c + epsilon - a - k)
             *
             * with
             *
             * (-1)**(a - b + k)*gamma(c + epsilon - b) * digamma(c + epsilon - a - k) / gamma(c + epsilon - a - k)
             *
             * and taking the limit epsilon -> 0 gives
             *
             * (-1)**(a - b + k) * gamma(c - b) * (-1)**(k + a - c + 1)(k + a - c)!
             * = (-1)**(c - b - 1)*Gamma(k + a - c + 1)
             */
            // Unlike in xsf, times the gammas of the prefactor, relative to the scale.
            let gammas = gamma_ratio([c.into()], [a.into(), (c - b).into()]);
            let f = (-1.0_f64).powf(m + n)
                * crate::binom(c - 1.0, b - 1.0)
                * poch(c - a + 1.0, m - 1.0)
                * ldexp_sat(gammas.mantissa.re, gammas.exponent - scaled_factor.exponent);
            // Unlike in xsf, use an integer power for the rounded index to retain small phases.
            factor = if n != c - a {
                (1.0 / z).powi(k as i32) * f
            } else {
                f / z.powf(k_)
            };
        }
        let term = factor;
        let f = (b + m + k_) * (k_ + a - c + 1.0) / ((k + 1) as f64 * (m + k_ + 1.0));
        factor *= f / z;
        k += 1;
        term
    }
}

fn hyp2f1_transform2_limit_finite_part_generator(
    b: f64,
    c: f64,
    m: f64,
    z: C64,
) -> impl FnMut() -> C64 {
    /* Initial finite sum in limit as a - b approaches a non-negative integer m. The limiting series
     * for the 1 - z transform also has an initial finite sum, but it is a standard hypergeometric
     * series. */
    // Unlike in xsf, start from 1, as the caller includes Gamma(m) / Gamma(c - b).
    let mut term = C64::ONE;
    let mut k: u64 = 0;
    move || {
        let output = term;
        let k_ = k as f64;
        term = term * (b + k_) * (c - b - k_ - 1.0) / ((k + 1) as f64 * (m - k_ - 1.0)) / z;
        k += 1;
        output
    }
}

#[allow(non_snake_case)]
fn lopez_temme_series_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> C64 {
    /* Lopez-Temme Series for Gaussian hypergeometric function [4].
     *
     * Converges for all z with real(z) < 1, including in the regions surrounding
     * the points exp(+- i*pi/3) that are not covered by any of the standard
     * transformations.
     */
    let mut n: u64 = 0;
    let mut phi_previous = 1.0;
    let mut phi = 1.0 - 2.0 * b / c;
    let mut Z = z * a / (z - 2.0);
    move || {
        if n == 0 {
            n += 1;
            return C64::ONE;
        }
        if n > 1 {
            // Update phi and Z for n>=2
            let n1 = (n - 1) as f64;
            let new_phi = (n1 * phi_previous - (2.0 * b - c) * phi) / (c + n1);
            phi_previous = phi;
            phi = new_phi;
            Z = Z * z / (z - 2.0) * ((a + n1) / n as f64);
        }
        n += 1;
        Z * phi
    }
}

/// The series, and its estimated cancellation as in `series_eval_cancellation`
///
/// Unless summed in pairs, a `retry` gives the term magnitudes that make it NaN.
fn lopez_temme_series(a: f64, b: f64, c: f64, z: C64, retry: Option<f64>) -> Cancelled {
    // Unlike in xsf, sum the even series (DLMF 15.8.13) when c = 2b. The odd
    // terms vanish, so the single-term rule would stop at the first one.
    let even = (c.is_finite() && c == 2.0 * b).then(|| {
        let w = z / (2.0 - z);
        let half_a = 0.5 * a;
        let terms = hypergeometric_series_generator(half_a, half_a + 0.5, b + 0.5, w * w);
        series_eval_lopez_temme(terms, b + 0.5, retry)
    });
    // Large negative powers after the Euler transform can still cancel, so
    // prefer the paired sum for a symmetric positive beta weight.
    if let Some(even) = even.filter(|even| even.cancellation <= 16.0) {
        return even;
    }
    // Unlike in xsf, keep the even sum if the paired sum fails, or else restart term by term.
    lopez_temme_paired(a, b, c, z, HYP2F1_MAXITER)
        .or(even)
        .unwrap_or_else(|| {
            series_eval_lopez_temme(lopez_temme_series_generator(a, b, c, z), c, retry)
        })
}

/// `prefactor` times the Lopez-Temme series of `p, q, c, z`, re-summed more accurately when the
/// sum `series` loses precision, or NaN
//
// Unlike in xsf, if the series loses more than HYP2F1_EPS to cancellation, also sum it with the
// upper parameters swapped, and prefer that if it cancels 4 times less. Take the preferred finite
// order with a small estimated error, smaller if swapped since that can replace a NaN, or a swapped
// sum that checks xsf's, or else NaN.
fn lopez_temme_reordered(p: f64, q: f64, c: f64, z: C64, prefactor: C64, series: Cancelled) -> C64 {
    // A usable sum estimates the swapped one, so stop that at twice the term magnitudes it can
    // have if preferred.
    let limit = if series.is_usable() {
        0.5 * series.cancellation * series.value.norm() * (1.0 - z * 0.5).norm().powf(q - p)
    } else {
        f64::INFINITY
    };
    let swapped_series = lopez_temme_series(q, p, c, z, Some(limit));
    let swapped = (1.0 - z * 0.5).powf(-q) * swapped_series.value; // Lopez-Temme prefactor.
    let swapped_error = swapped_series.cancellation * f64::EPSILON;
    let primary =
        Some(prefactor * series.value).filter(|value| series.is_usable() && value.is_finite());
    let certified = Some(swapped)
        .filter(|value| swapped_error <= crate::ext::HYP2F1_RECUR_TOL && value.is_finite());
    // xsf's sum, if the usable swapped one confirms it within its estimated error, or else the
    // swapped one, if within USABLE_ERROR of xsf's. Both are finite if the distance is that small.
    let checked = || {
        let xsf = prefactor * series_eval(lopez_temme_series_generator(p, q, c, z));
        let distance = (xsf - swapped).norm() / swapped.norm();
        let confirmed = if distance <= swapped_error {
            xsf
        } else {
            swapped
        };
        (swapped_series.is_usable() && distance <= USABLE_ERROR).then_some(confirmed)
    };
    let preferred = if 4.0 * swapped_series.cancellation < series.cancellation {
        certified.or_else(checked).or(primary)
    } else {
        primary.or(certified).or_else(checked)
    };
    preferred.unwrap_or(C64::NAN)
}

/// The Lopez-Temme series and its estimated cancellation, once a tail bound holds for a pair of
/// terms within `max_terms` terms, or `None` (always unless a, b > 0 and b < c)
fn lopez_temme_paired(a: f64, b: f64, c: f64, z: C64, max_terms: u64) -> Option<Cancelled> {
    if !(a > 0.0 && a.is_finite() && b > 0.0 && b < c && c.is_finite()) {
        return None;
    }
    let mut terms = lopez_temme_series_generator(a, b, c, z);
    // Unlike in xsf, bound both tails using decreasing absolute beta moments.
    // For odd moments, reflection gives a beta-weight difference of fixed sign.
    let w2 = (z / (z - 2.0)).norm_sqr();
    let upper = a.max(1.0);
    let ratio = |n: f64| w2 * ((upper + n) / (n + 1.0)) * ((upper + n + 1.0) / (n + 2.0));
    // The ratio does not increase with n, so no pair can stop if the last one cannot.
    if ratio((max_terms - 2) as f64) >= 1.0 {
        return None;
    }
    let (mut sum, mut abs_sum) = (C64::ZERO, 0.0);
    for n in (0..max_terms - 1).step_by(2) {
        let even = terms();
        let odd = terms();
        sum += even;
        sum += odd;
        // Unlike in xsf, skip Euclidean norms when the pair is too large to stop.
        // Two bounds sqrt(2) with rounding room; half-MAX keeps the norms finite.
        let pair_abs = even.l1_norm() + odd.l1_norm();
        abs_sum += pair_abs;
        let sum_abs = sum.l1_norm();
        if pair_abs.max(sum_abs) <= 0.5 * f64::MAX && pair_abs > 2.0 * HYP2F1_EPS * sum_abs {
            continue;
        }
        let last = even.norm() + odd.norm();
        let tol = sum.norm() * HYP2F1_EPS;
        if !last.is_finite() || !tol.is_finite() {
            return None;
        }
        if last < tol {
            let ratio = ratio(n as f64);
            if last * ratio < tol * (1.0 - ratio) {
                return Some(Cancelled {
                    value: sum,
                    cancellation: abs_sum / sum_abs,
                });
            }
        }
    }
    None
}

/// The parameters (p, q) of the Lopez-Temme series F(p, q; c; z) for F(a, b; c; z), its
/// prefactor, and the Euler transform prefactor if that applies
fn lopez_temme_form(a: f64, b: f64, c: f64, z: C64) -> (f64, f64, C64, Option<C64>) {
    /* This condition for applying Euler Transformation (DLMF 15.8.1)
     * was determined empirically to work better for this case than that
     * used in Zhang and Jin's implementation for |z| < 0.9,
     *  real(z) >= 0. */
    if (c - a <= a && c - b < b) || (c - a < a && c - b <= b) {
        let prefactor = (1.0 - z * 0.5).powf(a - c); // Lopez-Temme prefactor
        let euler = (1.0 - z).powf(c - a - b); // Euler transform prefactor.
        return (c - a, c - b, prefactor, Some(euler));
    }
    (a, b, (1.0 - z * 0.5).powf(-a), None) // Lopez-Temme prefactor.
}

/// The result, and its estimated cancellation as in `series_eval_cancellation`
//
// Unlike in xsf, `settled` or a cancelled gamma ratio sums the series with `series_eval_settled`.
fn hyp2f1_transform1_limiting_case(a: f64, b: f64, c: f64, z: C64, settled: bool) -> Cancelled {
    /* 1 - z transform in limiting case where c - a - b approaches an integer m. */
    let m = (c - a - b).round();
    // gamma(c) / (gamma(a) gamma(b)) for c near a + b + m and m >= 0.
    // Unlike in xsf, cancel it as four gammas if the ratio over- or underflows.
    // (in the order of a and b whose rising factorial is representable, if any)
    let mut recovered = false;
    let mut cancelled = |ratio: f64, a: f64, b: f64, m: f64| {
        if ratio.is_finite() && ratio != 0.0 {
            return ratio;
        }
        let product = [(a, b), (b, a)]
            .into_iter()
            .map(|(a, b)| four_gammas(c, 1.0, a, b + m + 1.0) * poch(b, m + 1.0))
            .find(|product| product.is_finite() && *product != 0.0);
        recovered = product.is_some();
        product.unwrap_or(ratio)
    };
    let mut result = C64::ZERO;
    let (prefactor, series_generator) = if m >= 0.0 {
        if m != 0.0 {
            let series_generator = hypergeometric_series_generator(a, b, 1.0 - m, 1.0 - z);
            result +=
                series_eval_finite(series_generator, m as u64) * four_gammas(m, c, a + m, b + m);
        }
        let gammas = cancelled(gamma(c) / (gamma(a) * gamma(b)), a, b, m);
        let prefactor = (1.0 - z).powf(m) * (neg_one_pow(m + 1.0) * gammas);
        let series_generator = hyp2f1_transform1_limit_series_generator(a + m, b + m, m, z);
        (prefactor, series_generator)
    } else {
        result = (1.0 - z).powf(m) * four_gammas(-m, c, a, b);
        let series_generator1 = hypergeometric_series_generator(a + m, b + m, 1.0 + m, 1.0 - z);
        result *= series_eval_finite(series_generator1, (-m) as u64);
        let gammas = cancelled(gamma(c) * (rgamma(a + m) * rgamma(b + m)), a + m, b + m, -m);
        let prefactor = neg_one_pow(m + 1.0) * gammas;
        let series_generator2 = hyp2f1_transform1_limit_series_generator(a, b, -m, z);
        (c64(prefactor, 0.0), series_generator2)
    };
    let series = if settled || recovered {
        series_eval_settled(series_generator)
    } else {
        series_eval_cancellation(series_generator)
    };
    let abs_sum = result.norm() + prefactor.norm() * series.value.norm() * series.cancellation;
    result += prefactor * series.value;
    let lost = cancellation(abs_sum, result);
    // Unlike in xsf, the sum with a cancelled ratio is NaN if it loses two thirds of the digits.
    if recovered && f64::EPSILON * lost >= f64::EPSILON.cbrt() {
        result = C64::NAN;
    }
    Cancelled {
        value: result,
        cancellation: lost,
    }
}

/// The product of the gamma functions at `num` divided by that at `den`, or 0 at a pole
fn gamma_ratio<const N: usize, const M: usize>(num: [Dd; N], den: [Dd; M]) -> ScaledC64 {
    // (directly where that tends to be more accurate, as in `four_gammas`)
    let (num_f64, den_f64) = (num.map(f64::from), den.map(f64::from));
    if in_direct_gamma_range(num_f64.into_iter().chain(den_f64)) {
        let direct = num_f64.map(gamma).iter().product::<f64>()
            * den_f64.map(rgamma).iter().product::<f64>();
        if direct.is_normal() {
            return ScaledC64::from(c64(direct, 0.0)).normalized();
        }
    }
    mul_gamma_ratio(num, den, || Scaled::new(1.0, 0, 0.0)).map_or(ScaledC64::from(C64::ZERO), |y| {
        let Scaled { m, e, .. } = fold_exp(y);
        ScaledC64 {
            mantissa: c64(m, 0.0),
            exponent: e.into(),
        }
    })
}

/// The generator, adding the norm of its k-th term (from 1) times k to `abs_sum`
fn weighted<'s>(mut g: impl FnMut() -> C64 + 's, abs_sum: &'s mut f64) -> impl FnMut() -> C64 + 's {
    let mut k = 0.0;
    move || {
        let term = g();
        k += 1.0;
        *abs_sum += k * term.norm();
        term
    }
}

/// The result, and its estimated cancellation as in `series_eval_cancellation`
fn hyp2f1_transform2_limiting_case(a: f64, b: f64, c: f64, m: f64, z: C64) -> (ScaledC64, f64) {
    /* 1 / z transform in limiting case where a - b approaches a non-negative integer m. Negative integer case
     * can be handled by swapping a and b. */
    // Unlike in xsf, take each ratio of gammas, including those of the first term of the series,
    // as one, as its factors can overflow or underflow alone. The second uses
    // (1 - x)_m / Gamma(x) = (-1)^m / Gamma(x - m) for x = c - b, so it is 0 at a pole of
    // Gamma(x - m), not 0 times infinity.
    //
    // (x is c - b rounded as in xsf, which in the Euler transform is exactly the original a.)
    let (a_dd, b_dd, c_dd, cb) = (Dd::new(a), Dd::new(b), Dd::new(c), Dd::new(c - b));
    let series_generator1 = hyp2f1_transform2_limit_finite_part_generator(b, c, m, z);
    // Unlike in xsf, also return the estimated cancellation, with the norm of the k-th term of
    // each sum weighted by k, as each step of its recurrence rounds.
    let (mut finite_abs, mut series_abs) = (0.0, 0.0);
    // Unlike in xsf, keep a separate scale for each power and ratio, and for each part, as they
    // can underflow or overflow alone.
    let finite_factor = power::neg_z_pow(z, -b) * gamma_ratio([c_dd, m.into()], [a_dd, cb]);
    let finite_sum = series_eval_finite(weighted(series_generator1, &mut finite_abs), m as u64);
    let finite = finite_factor.mantissa * finite_sum;
    let n = c - a;
    let n_int = n.round();
    let near_int = (n - n_int).abs() < HYP2F1_EPS;
    // Unlike in xsf, take x - m as c - a for c - a below 1/2 and not near an integer, so that a
    // zero of 1 / Gamma(c - a) cancels the pole of digamma(c - a) exactly.
    let x_m = if near_int || n_int > 0.0 {
        cb - m
    } else {
        Dd::new(n)
    };
    let mut factor = gamma_ratio([c_dd, b_dd + m], [a_dd, b_dd, x_m, (m + 1.0).into()]);
    factor.mantissa *= neg_one_pow(m);
    let series = if near_int {
        // Unlike in xsf, pass an integer at the positive-parameter coefficient switch.
        let n = if n_int > 0.0 && b > 0.0 { n_int } else { n };
        let g = hyp2f1_transform2_limit_series_cminus_a_int_generator(a, b, c, m, n, factor, z);
        series_eval(weighted(g, &mut series_abs))
    } else {
        let g = hyp2f1_transform2_limit_series_generator(a, b, c, m, factor.mantissa.re, z);
        series_eval(weighted(g, &mut series_abs))
    };
    let prefactor = power::neg_z_pow(z, -a);
    sum_scaled_parts([
        (finite, finite_factor.mantissa.norm() * finite_abs, finite_factor.exponent),
        (
            prefactor.mantissa * series,
            prefactor.mantissa.norm() * series_abs,
            prefactor.exponent + factor.exponent,
        ),
    ])
}

/// The sum of the (value, sum of absolute terms, exponent) parts at the largest exponent of a
/// nonzero part, and its estimated cancellation as in `series_eval_cancellation`
fn sum_scaled_parts(parts: [(C64, f64, f64); 2]) -> (ScaledC64, f64) {
    let exponent = parts
        .iter()
        .filter_map(|&(part, _, exponent)| (part != C64::ZERO).then_some(exponent))
        .reduce(f64::max)
        .unwrap_or(0.0);
    let [(first, first_abs), (second, second_abs)] =
        parts.map(|(part, part_abs, part_exponent)| {
            let shift = part_exponent - exponent;
            let part = ScaledC64 {
                mantissa: part,
                exponent: shift,
            };
            (part.to_c64(), ldexp_sat(part_abs, shift))
        });
    let mantissa = first + second;
    (ScaledC64 { mantissa, exponent }, cancellation(first_abs + second_abs, mantissa))
}

#[cold]
fn hyp2f1_euler_cut_lip(a: f64, b: f64, c: f64, z: C64) -> Option<C64> {
    // Bound the cost of the extra compensated polynomial evaluation.
    const MAX_DEGREE: f64 = 32.0;
    let terminating = |a: f64, b: f64| {
        let (p, low) = two_sum(c, -b);
        let n = -p;
        ((1.0..=MAX_DEGREE).contains(&n) && n.fract() == 0.0 && low == 0.0 && 0.0 < a && a < c)
            .then_some((n, a))
    };
    let (n, a) = terminating(a, b).or_else(|| terminating(b, a))?;
    let d = Dd::new(-n) - a;
    let nearest = d.hi.round();
    let delta = (d - nearest).hi;
    if !(delta != 0.0 && delta.abs() < 0.5 && delta.abs() <= f64::EPSILON * d.hi.abs()) {
        return None;
    }
    let phase = c64(cospi(delta), -z.im.signum() * sinpi(delta));
    let polynomial = hyp2f1_poly_dd(-n, Dd::new(c) - a, c, z.re)?;
    let r = Dd::new(z.re) - 1.0;
    let certified_lip = z.im == 0.0 || {
        // For these ordinates, |1 - z*t| <= x on the vertical segment for 0 <= t <= 1.
        // The beta integral (DLMF 15.6.1) then bounds |P'| by n*x^(n-1).
        // The normalized power contributes |d P|/(x-1).
        let (p, _) = polynomial.to_f64();
        let shift = z.im.abs() * (n * z.re.powf(n - 1.0) + d.hi.abs() * p.abs() / r.hi);
        let limit = f64::EPSILON / 4.0 * p.abs() * phase.re.abs().min(phase.im.abs());
        p.is_normal() && limit.is_normal() && shift <= limit
    };
    let polynomial = certified_lip.then_some(polynomial)?;
    let (v, err) = polynomial.with_l(d * ln_dd(r)).to_f64();
    let value = phase * (neg_one_pow(nearest) * v);
    ([value.re.abs(), value.im.abs()]
        .into_iter()
        .all(|a| (2.0 * f64::MIN_POSITIVE..=0.5 * f64::MAX).contains(&a))
        && err <= 8.0 * f64::EPSILON)
        .then_some(value)
}

const CERT_U: f64 = f64::EPSILON / 2.0;
const OUTWARD: f64 = 1.0 + 32.0 * CERT_U;
// Floors keep the products of nonnegative enclosure factors normal.
const CERT_FLOOR: f64 = pow2(-400);

fn seed_series(a: f64, b: f64, c: f64, z: C64) -> Option<C64> {
    let (mut sum, mut term) = (C64::ONE, C64::ONE);
    let radius = z.norm() * OUTWARD;
    for j in 0..HYP2F1_MAXITER {
        let j = j as f64;
        let p = (a + j) / (1.0 + j);
        let q = (b + j) / (c + j);
        let rho = radius * p.max(1.0) * q.max(1.0);
        if rho < 1.0 && term.l1_norm() * rho <= CERT_U * sum.l1_norm() * (1.0 - rho) {
            return Some(sum);
        }
        term = term * (p * q) * z;
        if !term.l1_norm().is_normal() {
            return None;
        }
        sum += term;
    }
    None
}

// Only used for bounded inputs. The normal floor also covers a subnormal squared modulus.
fn mod_upper(z: C64) -> f64 {
    (z.re.mul_add(z.re, z.im * z.im).max(f64::MIN_POSITIVE)).sqrt() * OUTWARD
}

fn steps_for_c_down(a: f64, b: f64, c: f64) -> Option<usize> {
    let shift = (a.max(b) - c).floor() + 1.0;
    if shift >= HYP2F1_MAXITER as f64 {
        return None;
    }
    let n = shift as usize;
    let (_, low) = two_sum(c, (n + 1) as f64);
    let t = c + n as f64;
    // Lower-index arguments have the same dyadic fractional part and no coarser
    // grid than c or c+n+1. Check the outer endpoint exactly, without a lattice cutoff.
    (low == 0.0 && t >= a.max(b)).then_some(n)
}

// DLMF 15.5.18, with t=c+k and F_t=F(a,b;t;z):
// F_(t-1) = ((1-l*z)*F_t + beta*z*F_(t+1))/(1-z).
// Carry an absolute bound over the vertical segment of width `vertical` from z.
#[cold]
fn c_down(a: f64, b: f64, c: f64, z: C64, vertical: f64) -> Option<(C64, f64)> {
    let n = steps_for_c_down(a, b, c)?;
    let r = (mod_upper(z) + vertical) * OUTWARD;
    let k = (1..=4).find(|&k| r <= 1.0 - pow2(-k))?;
    // At a seed with t >= max(a,b), the absolute series is at most
    // (1-r)^(-min(a,b)), by comparison with the binomial series.
    let e = f64::from(k) * a.min(b).ceil();
    if e > 1023.0 {
        return None;
    }
    let (mut u, mut v) = (pow2(e as i32), pow2(e as i32));
    let t0 = c + n as f64;
    let (mut f, mut fp) = (seed_series(a, b, t0, z)?, seed_series(a, b, t0 + 1.0, z)?);
    let w = 1.0 - z;
    let inv = C64::ONE / w;
    // On a vertical segment from the real axis, |1-z| is smallest at its base.
    let inv_bound = OUTWARD / w.re.mul_add(w.re, w.im * w.im).sqrt();
    for k in (1..=n).rev() {
        let t = c + k as f64;
        let d = t - 1.0;
        let l = (2.0 * t - 1.0 - a - b) / d;
        let beta = (t - a) * (t - b) / (t * d);
        let p = 1.0 - l * z;
        let next = (p * f + (beta * z) * fp) * inv;
        // t and d are exact, and nonzero beta factors stay normal under the
        // shift cap. The coefficient 1-l*z needs an additive allowance for cancellation.
        let l_error = 16.0 * CERT_U * ((2.0 * t).abs() + 1.0 + a + b) / d.abs();
        let p_error = l_error * r + 8.0 * CERT_U * (1.0 + l.abs() * r);
        let pa =
            ((mod_upper(p) + p_error + l.abs() * vertical) * inv_bound * OUTWARD).max(CERT_FLOOR);
        let qa = (beta.abs() * OUTWARD * r * inv_bound * OUTWARD).max(CERT_FLOOR);
        let new_u = ((pa * u + qa * v) * OUTWARD).max(CERT_FLOOR);
        (f, fp, u, v) = (next, f, new_u, u);
        if !f.l1_norm().is_normal() || !u.is_normal() {
            return None;
        }
    }
    Some((f, u))
}

// The true component is at most `bound` in magnitude. If the
// candidates are separated by more than twice the resulting error radius,
// the triangle inequality proves that the replacement is closer.
fn closer_component(candidate: f64, legacy: f64, bound: f64) -> Option<f64> {
    if candidate == legacy {
        Some(legacy)
    } else {
        (candidate.is_finite()
            && next_down_nonneg((candidate - legacy).abs())
                > 2.0 * next_up_nonneg(candidate.abs() + bound))
        .then_some(candidate)
    }
}

#[cold]
fn certified_c_down(a: f64, b: f64, c: f64, z: C64, current: C64) -> Option<C64> {
    if !current.is_finite() {
        return None;
    }
    let tiny = z.im.abs();
    let scaled = tiny != 0.0 && tiny < f64::MIN_POSITIVE.sqrt() * f64::EPSILON.sqrt();
    let (value, re_bound, im_bound) = if scaled {
        let real_z = c64(z.re, 0.0);
        // Enlarge the segment to keep the coefficient bounds in normal range.
        let vertical = tiny.max(CERT_FLOOR);
        // DLMF 15.5.1: F'=(ab/c)F(a+1,b+1;c+1). Use the same
        // enclosure on the shifted function over the segment.
        let [(ap, da), (bp, db), (cp, dc)] = [a, b, c].map(|v| two_sum(v, 1.0));
        if da != 0.0 || db != 0.0 || dc != 0.0 {
            return None;
        }
        let (f, u) = c_down(a, b, c, real_z, vertical)?;
        closer_component(f.re, current.re, u)?;
        let (shifted, shifted_u) = c_down(ap, bp, cp, real_z, vertical)?;
        let scale = (a * b).max(CERT_FLOOR) / c.abs() * OUTWARD;
        let im_bound = next_up_nonneg((scale * shifted_u * OUTWARD) * tiny);
        (c64(f.re, ((a * b) / c) * shifted.re * z.im), u, im_bound)
    } else {
        let (f, u) = c_down(a, b, c, z, 0.0)?;
        (f, u, if tiny == 0.0 { 0.0 } else { u })
    };
    Some(c64(
        closer_component(value.re, current.re, re_bound)?,
        closer_component(value.im, current.im, im_bound)?,
    ))
}

// DLMF 15.10.17 solved for w1 = F(a,b;c;z), with w2 = z^(1-c) (1-z)^(c-a-b) F(1-a,1-b;2-c;z)
// (DLMF 15.10.12, 15.8.1), w3 = F(a,b;a+b+1-c;1-z), and Gamma(c-1) by reflection (DLMF 5.5.3),
// and an estimated bound of its absolute error. For a, b > 0 > c, the gamma ratios are positive
// except for the sine.
fn hyp2f1_negative_c(a: f64, b: f64, c: f64, z: C64) -> Option<(C64, f64)> {
    let w = 1.0 - z;
    // For |1 - z| > 1, this is NaN unless the terms first decrease below the tolerance.
    let (w3, w3_tail) = series_eval_tail(hypergeometric_series_generator(a, b, a + b + 1.0 - c, w));
    if !(w3.value.is_finite() && w3_tail.is_finite()) {
        return None;
    }
    let (f2, f2_tail) =
        series_eval_tail(hypergeometric_series_generator(1.0 - a, 1.0 - b, 2.0 - c, z));
    // The gamma ratios can over- or underflow where the terms do not, so keep their logs.
    let one_c = -Dd::new(c) + 1.0;
    let (a1, b1) = (one_c + a, one_c + b);
    let one = || Scaled::new(1.0, 0, 0.0);
    let g1 = mul_gamma_ratio([a1, b1], [a1 + b, one_c], one)?;
    // (with Gamma(a) = Gamma(a + 1) / a, as `ln_gamma` fails at a subnormal a)
    let (a_1, b_1) = (Dd::new(a) + 1.0, Dd::new(b) + 1.0);
    let g2 = mul_gamma_ratio([a1, b1], [a_1, b_1, one_c, one_c + 1.0], one)?;
    let sin = sinpi(c);
    let ln_g2 = g2.l.hi + a.ln() + b.ln() + (PI / sin.abs()).ln();
    let g = g1.l.exp().hi;
    let t1 = if g.is_normal() {
        g * w3.value
    } else {
        (w3.value.ln() + g1.l.hi).exp()
    };
    let (ln_z, ln_w) = ((1.0 - c) * z.ln(), (c - a - b) * w.ln());
    let t2 = (ln_z + ln_w + ln_g2).exp() * sin.signum() * f2.value;
    // (with the tail of each series, and the rounding of the exponent of t2)
    let exponent_cancellation = ln_z.norm() + ln_w.norm() + ln_g2.abs();
    let t1_error = f64::EPSILON * w3.cancellation + w3_tail + g1.err;
    let t2_error = f64::EPSILON * (f2.cancellation + exponent_cancellation) + f2_tail + g2.err;
    Some((t1 + t2, t1.norm() * t1_error + t2.norm() * t2_error))
}

fn hyp2f1_pfaff(a: f64, b: f64, c: f64, z: C64) -> C64 {
    let (a, b) = if 0.0 < b && b < a && a < c {
        (b, a)
    } else {
        (a, b)
    };
    let z_pfaff = z / (z - 1.0);
    let series_generator = hypergeometric_series_generator(a, c - b, c, z_pfaff);
    let positive_parameters = a > 0.0 && a.is_finite() && b > 0.0 && b < c && c.is_finite();
    if z.is_real() || !positive_parameters {
        return (1.0 - z).powf(-a) * series_eval(series_generator);
    }
    // Unlike in xsf, evaluate the Pfaff prefactor without rounding |1-z| first.
    let power = pfaff_power(a, z);
    // Unlike in xsf, sum this positive-parameter Pfaff series to machine epsilon, and continue it
    // analytically if it cancels.
    let series = series_eval_cancellation_to(series_generator, f64::EPSILON);
    // (only for a normal power, as the continuation's error is then relative to its value)
    let sum = if series.is_continued() && power.norm().is_normal() {
        hyp2f1_continuation(a, c - b, c, z_pfaff, series.value.into()).unwrap_or(series.value)
    } else {
        series.value
    };
    power * sum
}

/// (1 - z)^(-a), from ln |1 - z| without rounding |1 - z| first
fn pfaff_power(a: f64, z: C64) -> C64 {
    let log_modulus = 0.5 * z.re.mul_add(z.re - 2.0, z.im * z.im).ln_1p();
    let (sin, cos) = (a * z.im.atan2(1.0 - z.re)).sin_cos();
    c64(cos, sin) * (-a * log_modulus).exp()
}

// Unlike in xsf, use the complementary incomplete beta integral (DLMF 8.17.4, 8.17.7)
// when c = a + 1 and d = 1 - b is small. The constant terms of the two series
// then combine as [d B(d,a) - 1] - [w^d - 1], for w = 1 - z.
fn hyp2f1_incomplete_beta(a: f64, b: f64, c: f64, z: C64) -> Option<C64> {
    // c - 1 is exact here, so the numerator comparison tests the identity exactly.
    if !(1.0 < c && c < 2.0) {
        return None;
    }
    let shape = c - 1.0;
    let d = if shape == a {
        1.0 - b
    } else if shape == b {
        1.0 - a
    } else {
        return None;
    };
    if !(d > 0.0 && d <= 1.0 / 64.0) {
        return None;
    }
    // Below the cut, a subnormal d*Im(z) can erase the phase of expm1(d*log(w)).
    // Use z F'(z) = a [(1-z)^(-b) - F(z)] where (Im(z)/(1-Re(z)))^2 underflows.
    let y = z.im;
    let tiny_phase = z.re < 1.0 && y != 0.0 && y.abs() < f64::MIN_POSITIVE.sqrt() * (1.0 - z.re);
    let z = if tiny_phase { c64(z.re, 0.0) } else { z };
    let w = 1.0 - z;
    let log_power = w.ln() * d;
    // log(d B(d,a)) = log Gamma(1+d) - integral of psi from a to a+d.
    let log_beta = crate::lgam1p(d) - digamma_integral(shape, d);
    let constant = log_beta.exp_m1() - crate::xsf::cexpm1(log_power);

    // Factor the first nonconstant term out of F(1-a,d;1+d;w) - 1. The normalized
    // series starts at one even for a subnormal w, and its term ratios are at most |w|.
    let mut term = C64::ONE;
    let mut n = 1.0;
    let tail = series_eval(|| {
        let value = term;
        term *= w * ((1.0 - shape + n) * (d + n) / ((n + 1.0) * (d + n + 1.0)));
        n += 1.0;
        value
    });
    let value = z.powf(-shape)
        * shape
        * (constant / d - (log_power.exp() * w) * ((1.0 - shape) / (1.0 + d)) * tail);
    if tiny_phase {
        let derivative = (shape / z.re) * ((1.0 - z.re).powf(d - 1.0) - value.re);
        return Some(c64(value.re, y * derivative));
    }
    Some(value)
}

/// `ln Gamma(x + h) - ln Gamma(x)`, the integral of digamma from `x` to `x + h`, or NaN if
/// `x < -HYP2F1_MAXITER`
//
// Gamma(x + 1) = x Gamma(x) (DLMF 5.5.1) moves x to at least max(1, 128 |h|), where
// three-point Gauss quadrature errs by less than 4.2e-4 (h/x)^6 |h| < 1e-16 |h|.
fn digamma_integral(mut x: f64, h: f64) -> f64 {
    if x < -(HYP2F1_MAXITER as f64) {
        return f64::NAN;
    }
    let mut integral = 0.0;
    while x < (128.0 * h.abs()).max(1.0) {
        integral -= (h / x).ln_1p();
        x += 1.0;
    }
    let mid = x + 0.5 * h;
    let step = 0.5 * (3.0_f64 / 5.0).sqrt() * h;
    integral + h / 18.0 * (5.0 * (digamma(mid - step) + digamma(mid + step)) + 8.0 * digamma(mid))
}

/// The 1 - z series (DLMF 15.8.4) for a balance c - a - b near zero, uniformly in the balance,
/// and its estimated cancellation as in `series_eval_cancellation`, where it is safe and more
/// accurate than the two separate series
//
// The k-th terms of the two series, whose poles in the balance d cancel, are paired: with
// A_k = Gamma(a+k) Gamma(b+k) / (Gamma(a+d) Gamma(b+d) k! Gamma(1-d+k)) and
// B_k = (a+d)_k (b+d)_k / (k! Gamma(1+d+k)), it is -Gamma(c) Gamma(-d) / (Gamma(a) Gamma(b))
// w^d / B_0 times the sum of B_k w^k expm1(log(A_k / B_k) - d log(w)), w = 1 - z.
fn hyp2f1_transform1_near_balance(a: f64, b: f64, c: f64, z: C64) -> Option<Cancelled> {
    let d = f64::from(Dd::new(c) - a - b);
    // Only where the unpaired series lose more than HYP2F1_EPS to their cancelling poles.
    let pole_loss = f64::EPSILON / d.abs();
    if pole_loss <= HYP2F1_EPS {
        return None;
    }
    // Beyond its direct ratio, `four_gammas` reflects x < 1/2 with a relative error near
    // eps |x| / dist(x, Z), worth it only where both the unpaired series, which lose eps / |d|
    // to their cancelling poles, and the limit for an integer balance, which is off by about
    // |d|, lose more.
    let tolerance = d.abs().max(pole_loss);
    let reflected = |x: f64| x < 0.5 && (x - x.round()).abs() <= (tolerance * x).abs();
    if !in_direct_gamma_range([c, -d, a, b]) && [a, b, c].into_iter().any(reflected) {
        return None;
    }
    let w = 1.0 - z;
    let log_power = w.ln() * d;
    // log(A_0 / B_0)
    let log_ratio =
        || crate::lgam1p(d) - crate::lgam1p(-d) - digamma_integral(a, d) - digamma_integral(b, d);
    let prefactor = log_power.exp() * -four_gammas(c, -d, a, b);
    hyp2f1_paired((a, b), 0.0, d, w, log_ratio, log_power, prefactor)
}

/// `prefactor` times the sum of w^k B_k expm1(log_ratio() + log(A_k / B_k) - log_power) over
/// k, with A_k = (a)_k (b)_k / ((1-d)_k (1+m)_k) and B_k = (a+d)_k (b+d)_k / (k! (1+m+d)_k),
/// and its estimated cancellation as in `series_eval_cancellation`, where that is accurate
fn hyp2f1_paired(
    (a, b): (f64, f64),
    m: f64,
    d: f64,
    w: C64,
    log_ratio: impl FnOnce() -> f64,
    log_power: C64,
    prefactor: C64,
) -> Option<Cancelled> {
    // The phase of exp(log_power) must not underflow before the prefactor restores its scale,
    // and an over- or underflowed prefactor loses the scale.
    let scaled = prefactor.is_finite() && prefactor.norm() >= f64::MIN_POSITIVE;
    if w.im != 0.0 && log_power.im.abs() < f64::MIN_POSITIVE || !scaled {
        return None;
    }
    // (only then, as its digamma integrals can be slow)
    let mut log_ratio = Some(log_ratio()).filter(|x| x.is_finite())?;
    let mut upper = C64::ONE; // w^k B_k
    let mut k = 0.0;
    let series = series_eval_settled(|| {
        let term = upper * crate::xsf::cexpm1(log_ratio - log_power);
        let (p, q, n) = (a + k, b + k, k + 1.0);
        let (shifted, nd) = ((p + d) * (q + d), n - d);
        // A_k+1 / B_k+1 = (A_k / B_k) (1 + d step)
        let step = ((1.0 + n / (n + m)) * p * q - (p + q + d) * nd) / (nd * shifted);
        log_ratio += (d * step).ln_1p();
        upper *= w * (shifted / (n * (n + m + d)));
        k = n;
        term
    });
    let value = series.value * prefactor;
    // A product whose imaginary part cancels by more than 1 / |d| loses more than the unpaired
    // series.
    let rotated = (series.value.re * prefactor.im).abs() + (series.value.im * prefactor.re).abs();
    let accurate = value.im.abs() >= d.abs() * rotated;
    (accurate && value.is_finite()).then_some(Cancelled { value, ..series })
}

/// The limit of the 1 - z series for an integer balance, if it is closer than `result`.
#[cold]
fn hyp2f1_transform1_near_limit(a: f64, b: f64, c: f64, z: C64, offset: f64, result: C64) -> C64 {
    // At the integer balance, which its gamma ratios assume.
    let limit = hyp2f1_transform1_limiting_case(a, b, a + b + (c - a - b).round(), z, true);
    // A heuristic for the relative derivative of the terms of the limit in the balance, from
    // their logarithm and digammas (DLMF 15.8.10). The sum amplifies it by its cancellation.
    let sensitivity =
        (1.0 - z).ln().norm() + [a, b, c].map(|x| digamma(x).abs()).iter().sum::<f64>();
    let error = (offset * sensitivity + f64::EPSILON) * limit.cancellation;
    // The margin covers the heuristic estimate.
    let agree = (result - limit.value).norm() <= 32.0 * error * limit.value.norm();
    if error < f64::EPSILON.cbrt() && !agree {
        limit.value
    } else {
        result
    }
}

// A step is at most this fraction of the distance to the nearest singular point, 0 or 1, for
// Taylor terms that decay asymptotically as 4^-k.
const TAYLOR_STEP_RATIO: f64 = 0.25;
// A step h also has |kappa h| at most this, for the local exponents kappa of the equation, as
// the Taylor terms of a solution near exp(kappa h) cancel by up to exp(|kappa h| - Re(kappa h)).
// (measured for large a, b and c near a + b: 4 sums up to 45% more terms, 16 cancels too much)
const TAYLOR_STEP_PHASE: f64 = 8.0;
// The Taylor terms of a continuation, as many as of 20 series: about that many were measured for
// a, b <= 400, c near a + b and |1 - z| near 1, as the steps shrink with the distance to 0
const TAYLOR_MAX_TERMS: u64 = 20 * HYP2F1_MAXITER;

/// F(p + h) and F'(p + h) from F(p) and F'(p), for each solution in `x`, by the Taylor series
/// at p of the hypergeometric equation (DLMF 15.10.1), and the sums of the absolute values of
/// the terms of the first, if they converge within `terms`, which it reduces by the number of
/// terms
//
// With e_k = F^(k)(p) h^k / k!, the equation gives p (1 - p) (k + 2) (k + 1) e_(k+2) =
// (k + a) (k + b) h^2 e_k - ((1 - 2p) k + c - (a + b + 1) p) (k + 1) h e_(k+1).
fn hyp2f1_taylor_step<const N: usize>(
    a: f64,
    b: f64,
    c: f64,
    p: C64,
    h: C64,
    x: [[C64; 2]; N],
    terms: &mut u64,
) -> Option<([[C64; 2]; N], [f64; 2])> {
    let g = h / (p * (1.0 - p));
    // c - (a + b + 1) p, without rounding a + b near p = 1 or cancelling near p = 0
    let t = if p.norm() < (1.0 - p).norm() {
        c - (a + b + 1.0) * p
    } else {
        let (a_plus_b, low) = two_sum(a, b);
        (c - a_plus_b - low - 1.0) + (1.0 - p) * (a + b + 1.0)
    };
    let (u, v, s) = (g * h, g * (1.0 - 2.0 * p), g * t);
    // For each solution, [e_k, e_(k+1)], the sums of the terms of F and h F', their absolute
    // values, and whether they have yet to converge
    let mut solutions = x.map(|[f, fp]| {
        let (e0, e1) = (f, fp * h);
        ([e0, e1], [e0 + e1, e1], [e0.l1_norm() + e1.l1_norm(), e1.l1_norm()], true)
    });
    let mut k = 0.0;
    while solutions.iter().any(|solution| solution.3) {
        let (ab, vs, w) = ((k + a) * (k + b), v * k + s, (k + 2.0) * (k + 1.0));
        for ([e0, e1], [sum, derivative], [abs_sum, abs_derivative], pending) in
            solutions.iter_mut().filter(|solution| solution.3)
        {
            *terms = terms.checked_sub(1)?;
            let e2 = (u * *e0 * ab - vs * *e1 * (k + 1.0)) / w;
            *sum += e2;
            *derivative += e2 * (k + 2.0);
            *abs_sum += e2.l1_norm();
            *abs_derivative += e2.l1_norm() * (k + 2.0);
            // (also for the derivative, whose terms have an extra factor k + 2)
            let tail = (k + 2.0) * (e1.l1_norm() + e2.l1_norm());
            // (stop at a non-finite sum, which no later term repairs)
            if !(sum.is_finite() && derivative.is_finite()) {
                return None;
            }
            *pending = !(tail <= f64::EPSILON * abs_sum.min(*abs_derivative));
            (*e0, *e1) = (*e1, e2);
        }
        k += 1.0;
    }
    let [abs_sum, abs_derivative] = solutions[0].2;
    let next = solutions.map(|(_, [sum, derivative], _, _)| [sum, derivative / h]);
    Some((next, [abs_sum, abs_derivative / h.norm()]))
}

/// F(a, b; c; z) for a, b, c > 0 and |1 - z| < 1 by analytic continuation of the Maclaurin
/// series at x0 = min(1/2, 1 - |1 - z|), along the real axis to 1 - |1 - z| and then along the
/// circle about 1 through z, if its error estimate is below half its distance from `result`, which
/// is then less accurate, or below the estimated error of `result`. For |1 - z| > 1 and Re z < 1,
/// the same for its Pfaff transform, if c - b > 0. For Re z > 1, where a small component may stem
/// from the cut, a `result` within sqrt(eps) relative to it is also kept if the estimate exceeds a
/// component, or only that component if the error of `result` is larger.
//
// On that path, the nearest singular point is 1, or 0 near x0 if |1 - z| > 1/2.
fn hyp2f1_continuation(
    a: f64,
    b: f64,
    c: f64,
    z: C64,
    Estimate {
        value: result,
        error: result_error,
    }: Estimate,
) -> Option<C64> {
    // Unlike in xsf, also continue the Pfaff transform: F(a, b; c; z) = (1 - z)^(-a)
    // F(a, c - b; c; z / (z - 1)) (DLMF 15.8.1), whose argument is within the unit circle about 1
    // if z is not, and off the cut for Re z < 1.
    let (power, b, z) = if (1.0 - z).norm() > 1.0 && z.re < 1.0 {
        // (only if normal, as the error estimate omits the rounding of a subnormal power)
        let power = Some(pfaff_power(a, z)).filter(|p| p.norm().is_normal())?;
        (Some(power), c - b, z / (z - 1.0))
    } else {
        (None, b, z)
    };
    let w = 1.0 - z;
    let r = w.norm();
    let (end, angle) = (1.0 - r, w.im.atan2(w.re));
    // (and only if the real axis reaches the circle, and the first Maclaurin term keeps a b)
    if !(a > 0.0 && b > 0.0 && c > 0.0 && r < 1.0 && end < 1.0 && (a * b).is_normal()) {
        return None;
    }
    let x0 = end.min(0.5);
    let start = c64(x0, 0.0);
    // The Maclaurin terms at x0 are positive, and grow up to `max_term_index`: sum past it before
    // stopping at a small term. Each of its n terms errs by about eps times the sum, by rounding
    // it, the roundings of the recurrence for it, or its part of the omitted tail.
    let maclaurin = |a: f64, b: f64, c: f64| {
        let head = (max_term_index(a, b, c, x0) + 1.0).max(2.0);
        if !(head <= HYP2F1_MAXITER as f64) {
            return None;
        }
        let (mut series_generator, mut n) = (hypergeometric_series_generator(a, b, c, start), 0.0);
        let mut terms = || {
            n += 1.0;
            series_generator()
        };
        let head = series_eval_fixed_length(&mut terms, C64::ZERO, head as u64);
        let sum = tools::series_eval(terms, head, f64::EPSILON, HYP2F1_MAXITER);
        Some((sum, n * sum.l1_norm()))
    };
    // F' = (a b / c) F(a + 1, b + 1; c + 1) (DLMF 15.5.1)
    let [f, g] = [(a, b, c), (a + 1.0, b + 1.0, c + 1.0)].map(|(a, b, c)| maclaurin(a, b, c));
    let ((f, f_error), (g, g_error)) = (f?, g?);
    let ab_c = a * b / c;
    let x = [f, g * ab_c];
    if !(x[0].is_finite() && x[1].is_finite()) {
        return None;
    }
    // The local exponents are the roots of p (1 - p) kappa^2 + (c - (a + b + 1) p) kappa = a b.
    let step = |p: C64, distance: f64| {
        let (s, t) = ((p * (1.0 - p)).norm(), (c - (a + b + 1.0) * p).norm());
        let kappa = (t + t.hypot(2.0 * (s * a * b).sqrt())) / (2.0 * s);
        (TAYLOR_STEP_RATIO * distance).min(TAYLOR_STEP_PHASE / kappa)
    };
    // The error of F at z is estimated by the rounding errors of F and F' of the Maclaurin sums
    // and of each step, times the first row of the product of the transfer matrices of the later
    // steps, whose columns start from F = 1, F' = 0 and F = 0, F' = 1. The estimate does not count
    // the roundings of the partial sums, so it may fall short by about 10%, which
    // the factor 2 in the comparison below absorbs.
    // (the Maclaurin sums first, before any step)
    let mut steps = vec![([[C64::ZERO; 2]; 2], [f_error, g_error * ab_c])];
    let (mut x, mut p, mut turned, mut terms) = (x, start, 0.0, TAYLOR_MAX_TERMS);
    let columns = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
    while p != z {
        let q = if p.re < end {
            c64((p.re + step(p, 1.0 - p.re)).min(end), 0.0)
        } else {
            turned = (turned + step(p, r.min(p.norm())) / r).min(angle.abs());
            if turned < angle.abs() {
                1.0 - r * C64::cis(turned.copysign(angle))
            } else {
                z
            }
        };
        // (not if it rounds to no step)
        if q == p {
            return None;
        }
        let h = q - p;
        let ([next, t0, t1], abs) =
            hyp2f1_taylor_step(a, b, c, p, h, [x, columns[0], columns[1]], &mut terms)?;
        steps.push(([t0, t1], abs));
        (x, p) = (next, q);
    }
    let (mut row, mut error) = ([C64::ONE, C64::ZERO], 0.0);
    for (transfer, abs) in steps.into_iter().rev() {
        error += (row[0].l1_norm() * abs[0] + row[1].l1_norm() * abs[1]) * f64::EPSILON;
        row = transfer.map(|t| row[0] * t[0] + row[1] * t[1]);
    }
    // (for the Pfaff transform, also the rounding of its argument z / (z - 1), times |F'| there)
    let (f, error) = power.map_or((x[0], error), |power| {
        let argument = x[1].norm() * z.norm() * f64::EPSILON;
        (power * x[0], power.norm() * (error + argument))
    });
    // (a first-order estimate, so only if well below one; NaN `result` is less accurate)
    let (tolerance, distance) = (f64::EPSILON.sqrt() * f.norm(), (result - f).norm());
    let less_accurate = result_error > error;
    let accurate = error < tolerance && (less_accurate || !(distance <= 2.0 * error));
    if !accurate {
        return None;
    }
    // (keep a close `result` if the estimate does not resolve a component)
    let unresolved = z.re > 1.0 && error >= f.re.abs().min(f.im.abs()) && distance < tolerance;
    if !unresolved {
        return Some(f);
    }
    // (only that smaller component if `result` is less accurate)
    less_accurate.then(|| {
        if f.re.abs() < f.im.abs() {
            c64(result.re, f.im)
        } else {
            c64(f.re, result.im)
        }
    })
}

fn hyp2f1_near_one(a: f64, b: f64, c: f64, z: C64) -> C64 {
    if let Some(value) = hyp2f1_incomplete_beta(a, b, c, z) {
        return value;
    }
    let balance = c - a - b;
    let offset = (balance - balance.round()).abs();
    let mut estimate = if offset < HYP2F1_EPS {
        // Removable singularity when c - a - b is an integer. Need to use limiting formula.
        // Unlike in xsf (scipy/xsf#388), use the Maclaurin series if c is so large that its
        // terms vanish fast.
        if c > HYP2F1_MAXITER as f64 * (1.0 + a.abs()) * (1.0 + b.abs()) {
            let series_generator = hypergeometric_series_generator(a, b, c, z);
            return series_eval(series_generator);
        }
        hyp2f1_transform1_limiting_case(a, b, c, z, false)
    } else {
        // Unlike in xsf, pair the two series for a balance near zero.
        hyp2f1_transform1_near_balance(a, b, c, z)
            .unwrap_or_else(|| series_eval_cancellation(hyp2f1_transform1_generator(a, b, c, z)))
    };
    // Unlike in xsf, if the 1 - z series lose more than HYP2F1_EPS to cancellation, also sum
    // the Maclaurin series, and keep the one that cancels less: for large c - a - b the 1 - z
    // terms are huge, but the first Maclaurin terms decay fast for large c > 0 (DLMF 15.12.2),
    // e.g. for hyp2f1(4.5, 2.5, 57, 0.64+0.65i).
    if c > 0.0 && estimate.is_lossy() {
        let series_generator = hypergeometric_series_generator(a, b, c, z);
        let maclaurin = series_eval_cancellation(series_generator);
        if maclaurin.cancellation < estimate.cancellation {
            estimate = maclaurin;
        }
    }
    let mut result = estimate.value;
    // Unlike in xsf, check the result against the limit for a nearly integer balance.
    // Rounding in c - a - b moves the poles that cancel between the two 1 - z series by a
    // relative eps / offset, for an error near eps / offset^2, where the limit is off by its
    // sensitivity to the balance: the two meet at an offset near cbrt(eps).
    // (not for a balance beyond the term limit, which the limit sums as many terms for)
    let near = offset >= HYP2F1_EPS && offset.powi(3) < f64::EPSILON;
    if near && balance.abs() <= HYP2F1_MAXITER as f64 {
        result = hyp2f1_transform1_near_limit(a, b, c, z, offset, result);
    }
    // Unlike in xsf, continue the Maclaurin series analytically if both still cancel.
    if estimate.is_lossy() {
        return hyp2f1_continuation(a, b, c, z, result.into()).unwrap_or(result);
    }
    result
}

/// The offset of a - b from an integer below which `hyp2f1_inverse` pairs the two 1/z series,
/// whose poles cancel between them: an empirical cutoff, beyond which the pairing was measured
/// to be less accurate
const TRANSFORM2_POLE_OFFSET: f64 = 1.0 / 1024.0;

/// The 1/z series (DLMF 15.8.2) for a - b near an integer m, uniformly in the offset, and its
/// estimated cancellation as in `series_eval_cancellation`, without and with the rounding of the
/// powers of -z (which the unpaired estimate omits too), where that is accurate
//
// With a = b + m + d, the k-th term of the first series is paired with the (k + m)-th term of
// the second, whose poles in d cancel. With p = b + m, their sum is -Gamma(c) Gamma(-m-d) /
// (Gamma(b) Gamma(c-a)) (-z)^(-a) times the sum of `hyp2f1_paired` for p, p - c + 1, w = 1/z,
// log_power = -d log(-z) and log_ratio = log(Gamma(1+m+d) Gamma(p) Gamma(c-a) / (m! Gamma(1-d)
// Gamma(a) Gamma(c-a+d))).
//
// Each parameter is passed with its complement to c, and a - b exactly, so that for the Euler
// transform c - (c - a) does not lose a tiny a, as in `hyp2f1_transform2_euler_generator`.
fn hyp2f1_transform2_near_integer(
    (a, ca): (f64, f64),
    (b, cb): (f64, f64),
    a_minus_b: Dd,
    c: f64,
    z: C64,
) -> Option<(Cancelled, f64)> {
    let ((a, ca), (b, cb), a_minus_b) = if a_minus_b < 0.0 {
        ((b, cb), (a, ca), -a_minus_b)
    } else {
        ((a, ca), (b, cb), a_minus_b)
    };
    let m = f64::from(a_minus_b).round();
    let d = f64::from(a_minus_b - m);
    // Not if the first series terminates, which the pairing divides by, or if c - b is an integer
    // of at most m, where the second series stops before its poles or vanishes, and those of the
    // first cancel in its gamma ratio.
    let terminates = ca > 0.0 && ca.fract() == 0.0;
    let stops = cb.fract() == 0.0 && cb <= m;
    if d.abs() >= TRANSFORM2_POLE_OFFSET || m > HYP2F1_MAXITER as f64 || terminates || stops {
        return None;
    }
    let w = 1.0 / z;
    let p = b + m;
    let log_ratio = || {
        let log_ratio = digamma_integral(m + 1.0, d) - crate::lgam1p(-d);
        log_ratio - digamma_integral(p, d) - digamma_integral(ca, d)
    };
    // Gamma(-m-d) / Gamma(b) = (-1)^m (b)_m / (1+d)_m Gamma(-d) / Gamma(p), a balanced ratio.
    let rising = neg_one_pow(m) * (poch(b, m) / poch(1.0 + d, m));
    let prefactor = (-z).powf(-a) * (rising * -four_gammas(c, -d, p, ca));
    let log_z = (-z).ln();
    let log_power = log_z * -d;
    // (1 - (c - b) keeps the exact zero of the second series if c - b is an integer)
    let tail = hyp2f1_paired((p, 1.0 - cb + m), m, d, w, log_ratio, log_power, prefactor)?;
    // The pairing rebuilds a, 1 - (c - a) and c - b - m from d, with its rounding error.
    let rebuilt = a.abs().min((1.0 - ca).abs()).min((cb - m).abs());
    let tail = Cancelled {
        // (from zero like the unpaired series, so that a zero imaginary part is +0)
        value: C64::ZERO + tail.value,
        cancellation: tail.cancellation.max(d.abs() / rebuilt),
    };
    if m == 0.0 {
        return Some((tail, tail.cancellation));
    }
    // The first m terms of the second series.
    let mut terms = hypergeometric_series_generator(b, 1.0 - cb, 1.0 - f64::from(a_minus_b), w);
    let (head, head_abs) = (0..m as u64)
        .map(|_| terms())
        .fold((C64::ZERO, 0.0), |(sum, abs), term| (sum + term, abs + term.norm()));
    let factor = (-z).powf(-b) * four_gammas(c, f64::from(a_minus_b), a, cb);
    let value = factor * head + tail.value;
    let abs_sum = factor.norm() * head_abs + tail.value.norm() * tail.cancellation;
    // (-z)^(-b) and (-z)^(-a) lose about |b log(-z)| and |a log(-z)| to rounding.
    let powers =
        (factor * head).norm() * (b * log_z).norm() + tail.value.norm() * (a * log_z).norm();
    let series = Cancelled {
        value,
        cancellation: cancellation(abs_sum, value),
    };
    Some((series, cancellation(abs_sum + powers, value)))
}

/// The factor by which the paired 1/z series must cancel less than the unpaired ones to be used,
/// as their estimates are rough: an empirical choice
const TRANSFORM2_PAIRED_MARGIN: f64 = 1.25;

/// The `unpaired` sum and its cancellation, or, if that loses more than HYP2F1_EPS, the `paired`
/// sum and its full cancellation if its comparable one is clearly less and keeps some digits
fn hyp2f1_transform2_pick(
    unpaired: Cancelled,
    paired: impl FnOnce() -> Option<(Cancelled, f64)>,
) -> Cancelled {
    if !unpaired.is_lossy() {
        return unpaired;
    }
    paired()
        .filter(|(paired, _)| {
            let lost = paired.cancellation;
            TRANSFORM2_PAIRED_MARGIN * lost < unpaired.cancellation && lost * f64::EPSILON < 1.0
        })
        .map_or(unpaired, |(paired, full)| Cancelled {
            cancellation: full,
            ..paired
        })
}

/// `hyp2f1_transform2_limiting_case` of the Euler transform (DLMF 15.8.1), times its power.
/// None unless c - a and c - b are exact and both the result before the power and the product
/// are nonzero and finite
fn hyp2f1_inverse_euler_limit(a: f64, b: f64, c: f64, m: f64, z: C64) -> Option<Cancelled> {
    // (only if c - (c - a) and c - (c - b) are exact, and b is an integer or not near one)
    let b_near_int = b != b.round() && (b - b.round()).abs() < HYP2F1_EPS;
    if c - (c - a) != a || c - (c - b) != b || b_near_int {
        return None;
    }
    let (result, lost) = hyp2f1_transform2_limiting_case(c - b, c - a, c, m, z);
    // Unlike in xsf, multiply by the power with a separate scale, as it can underflow or
    // overflow alone. (The result before it must be in range, as elsewhere the Euler form
    // tends to be less accurate.) Its degree is rounded, so correct it to first order.
    let degree = Dd::new(c) - a - b;
    let rounding = 1.0 + degree.lo * (1.0 - z).ln();
    let value = (power::one_minus_z_pow(z, degree.hi) * result * rounding).to_c64();
    let in_range = |x: C64| x.norm().is_finite() && x != C64::ZERO;
    let cancelled = Cancelled {
        value,
        cancellation: lost,
    };
    (in_range(result.to_c64()) && in_range(value)).then_some(cancelled)
}

/// `hyp2f1_inverse` for a - b near an integer, from the limit of the 1/z transform, also on the
/// Euler transform if `euler`
fn hyp2f1_inverse_limit(a: f64, b: f64, c: f64, z: C64, euler: bool) -> C64 {
    let (a, b) = if b > a { (b, a) } else { (a, b) };
    let m = (a - b).round();
    let euler_limit = if euler {
        hyp2f1_inverse_euler_limit(a, b, c, m, z)
    } else {
        None
    };
    // Unlike in xsf, if the Euler result is not used or loses more than HYP2F1_EPS to
    // cancellation, use the one that cancels less, or NaN if that overflows or is not usable
    // (unless it is zero, as if it underflows).
    let best = euler_limit.filter(|l| !l.is_lossy()).unwrap_or_else(|| {
        let (result, lost) = hyp2f1_transform2_limiting_case(a, b, c, m, z);
        let limit = Cancelled {
            value: result.to_c64(),
            cancellation: lost,
        };
        euler_limit
            .filter(|l| l.cancellation < lost)
            .unwrap_or(limit)
    });
    let usable = best.value == C64::ZERO || best.is_usable();
    let result = if best.value.norm().is_finite() && usable {
        best.value
    } else {
        C64::NAN
    };
    // Unlike in xsf, continue the Maclaurin series analytically if the limit is NaN or lossy and
    // the continuation is more accurate.
    let estimate = Estimate {
        value: result,
        error: f64::EPSILON * best.cancellation * result.norm(),
    };
    let continued = (result.is_nan() || best.is_lossy())
        .then(|| hyp2f1_continuation(a, b, c, z, estimate))
        .flatten();
    continued.unwrap_or_else(|| {
        // Unlike in xsf, if that fails for a NaN limit, use the Euler form even if it does not
        // reduce a and b, if usable.
        if euler || !result.is_nan() {
            return result;
        }
        hyp2f1_inverse_euler_limit(a, b, c, m, z)
            .filter(|l| l.is_usable())
            .map_or(result, |l| l.value)
    })
}

fn hyp2f1_inverse(a: f64, b: f64, c: f64, z: C64) -> C64 {
    // Unlike in xsf (scipy/xsf#387), also try it on the Euler transform (DLMF 15.8.1) if that
    // reduces a and b, as for |z| < 0.9.
    let euler = c - a < a && c - b < b;
    if (a - b - (a - b).round()).abs() < HYP2F1_EPS {
        return hyp2f1_inverse_limit(a, b, c, z, euler);
    }
    // Unlike in xsf, pair the two series for a - b near an integer if that clearly cancels less,
    // counting for a near b both parts of each unpaired term, whose poles cancel within it.
    let near = (a - b).abs() < TRANSFORM2_POLE_OFFSET;
    let mut series = hyp2f1_transform2_pick(
        series_eval_pairs(hyp2f1_transform2_generator(a, b, c, z), near),
        || hyp2f1_transform2_near_integer((a, c - a), (b, c - b), Dd::new(a) - b, c, z),
    );
    if euler {
        // use the one that cancels less (also paired)
        let euler_series = hyp2f1_transform2_pick(
            series_eval_pairs(hyp2f1_transform2_euler_generator(a, b, c, z), near),
            || hyp2f1_transform2_near_integer((c - a, a), (c - b, b), Dd::new(b) - a, c, z),
        );
        if euler_series.cancellation < series.cancellation {
            // (times (1 - z)^(c - a - b), if both are normal)
            let prefactor = (1.0 - z).powf(c - a - b);
            if prefactor.norm().is_normal() && euler_series.value.norm().is_normal() {
                series = Cancelled {
                    value: prefactor * euler_series.value,
                    ..euler_series
                };
            }
        }
    }
    // Unlike in xsf, if the 1/z series lose more than HYP2F1_EPS to cancellation, also sum the
    // Lopez-Temme series for Re z < 1 where its tail is bounded.
    if z.re < 1.0 && series.is_lossy() {
        if let Some(value) = hyp2f1_lopez_temme_retry(a, b, c, z, series.cancellation) {
            return value;
        }
    }
    // Unlike in xsf, continue the Maclaurin series analytically if the 1/z series cancel.
    if series.is_lossy() {
        return hyp2f1_continuation(a, b, c, z, series.value.into()).unwrap_or(series.value);
    }
    series.value
}

/// F(a, b; c; z) from the Lopez-Temme series, if it cancels less than half as much as the 1/z
/// series lost to `cancellation` (the estimates use different norms) and its estimated summation
/// roundoff is below sqrt(EPSILON)
fn hyp2f1_lopez_temme_retry(a: f64, b: f64, c: f64, z: C64, cancellation: f64) -> Option<C64> {
    let (p, q, prefactor, euler) = lopez_temme_form(a, b, c, z);
    // (not if a prefactor component is subnormal or zero, as the result would lose it)
    let normal = |f: C64| f.re.is_normal() && (f.im.is_normal() || z.is_real());
    if !normal(prefactor) || !euler.is_none_or(normal) {
        return None;
    }
    let paired = lopez_temme_paired(p, q, c, z, LOPEZ_TEMME_RETRY_MAXITER)?;
    let result = prefactor * paired.value;
    let value = euler.map_or(result, |euler| euler * result);
    let roundoff = f64::EPSILON * paired.cancellation;
    let accurate = 2.0 * paired.cancellation < cancellation && roundoff < f64::EPSILON.sqrt();
    (accurate && value.norm().is_normal()).then_some(value)
}

// Unlike in xsf, raise the Ferrers order for F(1/2,1/2;1/2-n;z) (DLMF 14.3.1,
// 14.10.1). This avoids cancellation in the Euler-transformed polynomial.
fn hyp2f1_ferrers(n: u64, z: C64) -> C64 {
    let y = z.im;
    let tiny = f64::MIN_POSITIVE.sqrt() * (1.0 - z.re).abs().min(1.0);
    let tiny_phase = y != 0.0 && y.abs() < tiny;
    let z = if tiny_phase {
        c64(z.re, 0.0_f64.copysign(y))
    } else {
        z
    };
    let w = 1.0 - z;
    let ratio = z / w;
    let mut previous = 1.0 / w.sqrt();
    let n_f = n as f64;
    // Sum the differences when successive values can be nearly equal. On the real
    // axis below the equator this also preserves the difference needed for F'_n.
    let (mut value, diff) =
        if (1.0 - ratio).norm_sqr() >= 0.25 || (tiny_phase && (1.0 - 2.0 * z.re) * n_f > 1.0) {
            let mut diff = -ratio * previous;
            let mut value = previous;
            let mut comp = C64::ZERO;
            tools::kahan_step(&mut value, &mut comp, diff);
            for j in 1..n {
                let k = (2.0 * j as f64).powi(2);
                diff = ratio * (previous / (k - 1.0) - diff);
                previous = value;
                tools::kahan_step(&mut value, &mut comp, diff);
            }
            (value, diff)
        } else {
            // At the equator F_n vanishes for odd n. Form 1-r from 1-2z to retain
            // the distance to the root.
            let linear = (1.0 - 2.0 * z) / w;
            let mut value = linear * previous;
            for j in 1..n {
                let k = (2.0 * j as f64).powi(2);
                let next = linear * value + ratio * (previous + previous / (k - 1.0));
                (previous, value) = (value, next);
            }
            (value, value - previous)
        };
    if tiny_phase {
        // z F'_n = -(n+1/2)(F_(n+1)-F_n). Use the next difference divided by z
        // before scaling by Im(z), so even a subnormal ordinate retains its phase.
        let h = diff - previous / ((2.0 * n_f).powi(2) - 1.0);
        let factor = (n_f + 0.5) / w.re;
        let component = if z.re > 1.0 { -h.im } else { h.re };
        let derivative = factor * component;
        let offset = if derivative.is_infinite() {
            (y * component) * factor
        } else {
            y * derivative
        };
        if z.re > 1.0 {
            value.re = offset;
        } else {
            value.im = offset;
            // At the odd equatorial root the leading real part is quadratic.
            // The hypergeometric ODE gives F''(1/2) = (4n+2) F'(1/2).
            if z.re == 0.5 && n % 2 == 1 {
                value.re = -((2.0 * n_f + 1.0) * y) * offset;
            }
        }
    }
    value
}

fn hyp2f1_c64(a: f64, b: f64, c: f64, z: C64) -> C64 {
    /* Special Cases
     * nan input */
    if a.is_nan() || b.is_nan() || c.is_nan() || z.is_nan() {
        return C64::NAN;
    }

    /* Takes constant value 1 when a = 0 or b = 0, even if c is a non-positive
     * integer. This follows mpmath. */
    if a == 0.0 || b == 0.0 {
        return C64::ONE;
    }
    let z_abs = z.norm();
    // Equals 1 when z i 0, unless c is 0.
    if z_abs == 0.0 {
        if c != 0.0 {
            return C64::ONE;
        }
        // Returning real part NAN and imaginary part 0 follows mpmath.
        return C64::NAN_RE;
    }
    let a_neg_int = a == a.trunc() && a < 0.0;
    let b_neg_int = b == b.trunc() && b < 0.0;
    let c_non_pos_int = is_nonpos_int(c);
    /* Diverges when c is a non-positive integer unless a is an integer with
     * c <= a <= 0 or b is an integer with c <= b <= 0, (or z equals 0 with
     * c != 0) Cases z = 0, a = 0, or b = 0 have already been handled. We follow
     * mpmath in handling the degenerate cases where any of a, b, c are
     * non-positive integers. See [3] for a treatment of degenerate cases. */
    if c_non_pos_int && !((a_neg_int && c <= a && a < 0.0) || (b_neg_int && c <= b && b < 0.0)) {
        return C64::INFINITY;
    }
    // Unlike in xsf, cancel a finite common upper and lower parameter:
    // F(a,b;a;z) = (1-z)^(-b) (DLMF 15.4.6).
    // At a pole of c, the other upper parameter has to terminate the series first.
    for (common, other, terminates) in [(a, b, b_neg_int), (b, a, a_neg_int)] {
        let cancels = common == c && c.is_finite() && other.is_finite();
        if cancels && (!c_non_pos_int || (terminates && c <= other)) {
            return power::pow_1mz(z, -other);
        }
    }
    /* Reduces to a polynomial when a or b is a negative integer.
     * If a and b are both negative integers, we take care to terminate
     * the series at a or b of smaller magnitude. This is to ensure proper
     * handling of situations like a < c < b <= 0, a, b, c all non-positive
     * integers, where terminating at a would lead to a term of the form 0 / 0. */
    if a_neg_int || b_neg_int {
        return polynomial::terminating(a, b, c, z).value;
    }
    // Unlike in xsf, try a beta contour for a large negative nonintegral upper parameter.
    if let Some(value) = beta_endpoints::evaluate(a, b, c, z, z_abs) {
        return value;
    }
    // Kummer's Theorem for z = -1; c = 1 + a - b (DLMF 15.4.26)
    if (z + 1.0).norm() < HYP2F1_EPS && (1.0 + a - b - c).abs() < HYP2F1_EPS && !c_non_pos_int {
        return c64(four_gammas(a - b + 1.0, 0.5 * a + 1.0, a + 1.0, 0.5 * a - b + 1.0), 0.0);
    }
    // Unlike in xsf, use Gauss summation before checking for a rounded Euler polynomial.
    if z.re == 1.0 && z.is_real() {
        // The Euler factor can give a signed real zero on the lower lip.
        let exact_zero = |y| if z.im.is_sign_negative() { y } else { 0.0 };
        if let Some(y) = hyp2f1_at_1(a, b, c, || f64::INFINITY, exact_zero) {
            return c64(y, 0.0);
        }
    }
    let c_minus_a_neg_int = c - a == (c - a).trunc() && c - a < 0.0;
    let c_minus_b_neg_int = c - b == (c - b).trunc() && c - b < 0.0;
    /* If one of c - a or c - b is a negative integer, reduces to evaluating
     * a polynomial through an Euler hypergeometric transformation.
     * (DLMF 15.8.1) */
    if c_minus_a_neg_int || c_minus_b_neg_int {
        let (terminating, other) = if c_minus_b_neg_int { (b, a) } else { (a, b) };
        let max_degree = (c - terminating).abs();
        if a == 0.5
            && b == 0.5
            && c.fract() == -0.5
            && (2.0..=(HYP2F1_MAX_FIXED_TERMS - 2) as f64).contains(&max_degree)
            && z.is_finite()
        {
            let value = hyp2f1_ferrers(max_degree as u64, z);
            if value.is_finite() {
                return value;
            }
        }
        if max_degree <= u64::MAX as f64 {
            // Unlike in xsf, retain the exact cut phase for certified tiny ordinates.
            // The power term alone limits |Im z| to about eps² (Re z - 1).
            if z.re > 1.0
                && z.re.is_finite()
                && z.im.abs() <= f64::EPSILON.powi(2) * (z.re - 1.0).min(1.0)
            {
                if let Some(value) = hyp2f1_euler_cut_lip(a, b, c, z) {
                    return value;
                }
            }
            let w = 1.0 - z;
            let regular = w.is_finite() && w != C64::ZERO;
            // Unlike in xsf, retain the phase of exactly integral Euler powers across the cut.
            let integral = other.fract() == 0.0 && regular && two_sum(c, -terminating).1 == 0.0;
            let prefactor = integral
                .then(|| cpow_integral(w, Dd::new(-max_degree) - other))
                .filter(|power| power.is_finite())
                .unwrap_or_else(|| w.powf(c - a - b));
            // Unlike in xsf, sum the polynomial like a directly terminating one.
            let terminating = polynomial::terminating(c - a, c - b, c, z);
            // Unlike in xsf, keep a small other parameter that c - other loses.
            let estimate = polynomial::complement(other, c, max_degree as u64, z, terminating);
            let polynomial = estimate.value;
            let value = prefactor * polynomial;
            // Unlike in xsf, carry the scale of a prefactor or product that under- or overflows.
            let peak = prefactor.linf_norm();
            let in_range = peak.is_normal();
            // (a prefactor below the range needs a resolved polynomial, as the error of a
            // cancelled sum would surface instead of underflowing; otherwise a finite nonzero
            // one is rescaled when the prefactor overflows or the product is not finite)
            let scale = if peak < f64::MIN_POSITIVE {
                estimate.resolved()
            } else {
                !(value.is_finite() && in_range)
                    && polynomial.is_finite()
                    && polynomial != C64::ZERO
            };
            if regular && scale {
                // (with the exponent as the sum of two words, which keeps an integral one exact)
                let Dd { hi, lo } = Dd::new(c) - a - b;
                let power = power::one_minus_z_pow_sqrt(z, hi) * power::one_minus_z_pow_sqrt(z, lo);
                let scaled = (power * polynomial).to_c64();
                // Keep a finite component of a product with the prefactor in range, and an
                // infinite one where the scaled power loses a component far below the other.
                let pick = |scaled: f64, value: f64| {
                    if (in_range && value.is_finite()) || (scaled == 0.0 && value.is_infinite()) {
                        value
                    } else {
                        scaled
                    }
                };
                return c64(pick(scaled.re, value.re), pick(scaled.im, value.im));
            }
            return value;
        }
        return C64::NAN;
    }
    /* Diverges as real(z) -> 1 when c <= a + b.
     * Todo: Actually check for overflow instead of using a fixed tolerance for
     * all parameter combinations like in the Fortran original. */
    if (1.0 - z.re).abs() < HYP2F1_EPS && z.is_real() && c - a - b <= 0.0 && !c_non_pos_int {
        return C64::INFINITY;
    }
    // Gauss's Summation Theorem for z = 1; c - a - b > 0 (DLMF 15.4.20).
    if z.re == 1.0 && z.is_real() && c - a - b > 0.0 && !c_non_pos_int {
        return c64(four_gammas(c, c - a - b, c - a, c - b), 0.0);
    }
    /* |z| < 0.9, z.real() >= 0. Use the Maclaurin Series.
     * -----------------------------------------------------------------------
     * Apply Euler Hypergeometric Transformation (DLMF 15.8.1) to reduce
     * size of a and b if possible. We follow Zhang and Jin's
     * implementation [1] although there is very likely a better heuristic
     * to determine when this transformation should be applied. As it
     * stands, this hurts precision in some cases. */
    if z_abs < 0.9 && z.re >= 0.0 {
        if c - a < a && c - b < b {
            let prefactor = (1.0 - z).powf(c - a - b);
            // Unlike in xsf, track cancellation and try a certified recurrence down in c, or else
            // `hyp2f1_negative_c`.
            if a > 0.0 && b > 0.0 && c < 0.0 {
                let mut series_generator = hypergeometric_series_generator(c - a, c - b, c, z);
                let (mut mass, mut terms) = (0.0, 0.0);
                let sum = series_eval(|| {
                    let term = series_generator();
                    mass += term.l1_norm();
                    terms += 1.0;
                    term
                });
                let value = prefactor * sum;
                if !(8.0 * CERT_U * terms * mass <= crate::ext::HYP2F1_RECUR_TOL * sum.l1_norm()) {
                    let current = certified_c_down(a, b, c, z, value).unwrap_or(value);
                    // Take each component of the connection, if it keeps some digits, that differs
                    // from `current` by more than twice the estimated error, so that `current` is
                    // estimated to be less accurate there.
                    return match hyp2f1_negative_c(a, b, c, z) {
                        Some((alt, error)) if error < alt.norm() => {
                            let pick = |alt: f64, current: f64| {
                                if (alt - current).abs() <= 2.0 * error {
                                    current
                                } else {
                                    alt
                                }
                            };
                            c64(pick(alt.re, current.re), pick(alt.im, current.im))
                        },
                        _ => current,
                    };
                }
                return value;
            }
            // Unlike in xsf, account for the rounding of `c - a` or `c - b` near a nonpositive
            // integer, as the later terms are proportional to their distance to it.
            let relative_rounding = |p: f64| {
                if p.round() > 0.0 {
                    0.0
                } else {
                    f64::EPSILON * p.abs() / (p - p.round()).abs()
                }
            };
            let rounding = relative_rounding(c - a).max(relative_rounding(c - b));
            return prefactor * hypergeometric_series(c - a, c - b, c, z, rounding);
        }
        return hypergeometric_series(a, b, c, z, 0.0);
    }
    /* Points near exp(iπ/3), exp(-iπ/3) not handled by any of the standard
     * transformations. Use series of López and Temme [5]. These regions
     * were not correctly handled by Zhang and Jin's implementation.
     * -------------------------------------------------------------------------*/
    if (0.9..1.1).contains(&z_abs) && (1.0 - z).norm() >= 0.9 && z.re >= 0.0 {
        let (p, q, prefactor, euler) = lopez_temme_form(a, b, c, z);
        let series = lopez_temme_series(p, q, c, z, None);
        let result = if series.is_lossy() {
            lopez_temme_reordered(p, q, c, z, prefactor, series)
        } else {
            prefactor * series.value
        };
        let value = euler.map_or(result, |euler| euler * result);
        // Unlike in xsf, continue the Maclaurin series analytically if this series cancels.
        if series.is_continued() {
            return hyp2f1_continuation(a, b, c, z, value.into()).unwrap_or(value);
        }
        return value;
    }
    /* z/(z - 1) transformation (DLMF 15.8.1). Avoids cancellation issues that
     * occur with Maclaurin series for real(z) < 0.
     * -------------------------------------------------------------------------*/
    if z_abs < 1.1 && z.re < 0.0 {
        return hyp2f1_pfaff(a, b, c, z);
    }
    /* 1 - z transformation (DLMF 15.8.4). */
    if (0.9..1.1).contains(&z_abs) {
        return hyp2f1_near_one(a, b, c, z);
    }
    /* 1/z transformation (DLMF 15.8.2). */
    hyp2f1_inverse(a, b, c, z)
}

pub trait Hyp2F1Arg: crate::sealed::Sealed {
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self;
}

// Unlike in xsf, use the symmetric quadratic transform (DLMF 15.8.21)
// for moderate positive real arguments.
fn hyp2f1_symmetric_real(h: f64, b: f64, c: f64, z: f64) -> Option<f64> {
    if h != b || !(0.5..=32.0).contains(&h) || c != 2.0 * h || !(z > 0.0 && z <= 0.75) {
        return None;
    }
    let u = 1.0 + (1.0 - z).sqrt();
    let rho = z / (u * u);
    let w = rho * rho;
    // Successive positive terms have ratio at most w <= 1/9.
    let tail = w / (1.0 - w);
    let (mut sum, mut term) = (1.0, 1.0);
    let mut n = 0.0;
    while term * tail > f64::EPSILON * 0.25 * sum {
        term *= w * (h + n) * (n + 0.5) / ((h + n + 0.5) * (n + 1.0));
        sum += term;
        n += 1.0;
    }
    Some((c * rho.ln_1p()).exp() * sum)
}

impl Hyp2F1Arg for f64 {
    #[inline]
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self {
        // Unlike in xsf, try the symmetric real fast path before Cephes.
        hyp2f1_symmetric_real(a, b, c, self)
            .unwrap_or_else(|| crate::xsf::cephes::hyp2f1(a, b, c, self))
    }
}

impl Hyp2F1Arg for num_complex::Complex<f64> {
    #[inline]
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self {
        // Unlike in xsf, exactly real arguments can use the symmetric real fast path.
        if self.im == 0.0 {
            if let Some(value) = hyp2f1_symmetric_real(a, b, c, self.re) {
                return Self::new(value, 0.0);
            }
        }
        self.via(|z| hyp2f1_c64(a, b, c, z))
    }
}

/// Gauss' hypergeometric function $_2F_1$
///
/// Corresponds to [`scipy.special.hyp2f1`][hyp2f1] in SciPy, and accepts both `f64` and
/// `num_complex::Complex<f64>` inputs for `z`.
///
/// [hyp2f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp2f1.html
///
/// # Definition
///
/// For $\abs z < 1$,
///
/// $$
/// \hyp 2 1 {a_1\enspace a_2} b z
/// = \sum\_{n=0}^\infty
///     {\rpow{a_1}{n} \ \rpow{a_2}{n} \over \rpow b n}
///     {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)), and the
/// analytic continuation elsewhere [^1]. A non-positive integer $a_1$ or $a_2$ equal to $-n$ makes
/// it a polynomial of degree at most $n$. At $z = 1$, away from poles of $b$, the series also
/// converges when $b - a_1 - a_2 > 0$ [^2].
///
/// For `f64` $z$, it is:
///
/// - NaN where the internal error estimates say that none of its digits are correct;
/// - usually NaN for $z > 1$, where it is generally complex, unless $a_1$, $a_2$, $b - a_1$ or
///   $b - a_2$ is a non-positive integer, so that it is a polynomial, or a polynomial times a power
///   of $1 - z$, if that is real;
/// - $\pm\infty$ at $z = 1$ where it diverges, and $+\infty$ e.g. at a pole of $b$, or where
///   $b - a_1$ or $b - a_2$ rounds to $b$.
///
/// Here, a parameter within about $10^{-13}$ of an integer can count as one.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/15.2>
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/15.4.E20>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function,
///   $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
#[must_use]
#[inline]
pub fn hyp2f1<T: Hyp2F1Arg>(a: f64, b: f64, c: f64, z: T) -> T {
    z.hyp2f1(a, b, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hyp2f1_continuation_maclaurin_growth() {
        // The Maclaurin terms at 1/2 first decay below eps and then grow.
        let half = c64(0.5, 0.0);
        let value = hyp2f1_continuation(1e-20, 700.0, 240.125, half, C64::NAN.into()).unwrap();
        xsref::np_assert_allclose!(
            [value.re, value.im],
            [1.000_007_275_456_677, 0.0],
            rtol = 1e-14
        );
        // Its terms stay below eps, so the sum rounds to 1, which the error estimate covers.
        let exact = c64(1.000_000_000_000_007_3, 0.0);
        assert_eq!(hyp2f1_continuation(1e-29, 700.0, 240.125, half, exact.into()), None);
    }
}
