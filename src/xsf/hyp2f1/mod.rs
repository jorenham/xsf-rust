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

use core::f64::consts::{FRAC_1_PI, PI};

use crate::ext::{Scaled, hyp2f1_poly_dd, ln_dd, mul_gamma_ratio};
use crate::utils::{is_nonpos_int, neg_one_pow, next_down_nonneg, next_up_nonneg, pow2};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, two_sum};
use crate::xsf::cephes::hyp2f1::hyp2f1_at_1;
use crate::xsf::cephes::lanczos::LANCZOS_G;
use crate::xsf::cephes::{gamma, gammasgn, lanczos_sum_expg_scaled, lgam, poch, rgamma};
use crate::xsf::complex::cpow_integral;
use crate::xsf::tools::{self, series_eval_fixed_length};
use crate::{cospi, digamma, sinpi};

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
fn in_direct_gamma_range(args: [f64; 4]) -> bool {
    /* Direct ratio tends to be more accurate for arguments in this range. Range
     * chosen empirically based on the relevant benchmarks in
     * scipy/special/_precompute/hyp2f1_data.py */
    args.iter().all(|x| x.abs() <= 100.0)
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

fn hyp2f1_transform2_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> C64 {
    /* 1/z transformation of standard series.*/
    let factor1 = (-z).powf(-a) * four_gammas(c, b - a, b, c - a);
    let factor2 = (-z).powf(-b) * four_gammas(c, a - b, a, c - b);
    let mut generator1 = hypergeometric_series_generator(a, a - c + 1.0, a - b + 1.0, 1.0 / z);
    let mut generator2 = hypergeometric_series_generator(b, b - c + 1.0, b - a + 1.0, 1.0 / z);
    move || factor1 * generator1() + factor2 * generator2()
}

/// `hyp2f1_transform2_generator(c - a, c - b, c, z)`, but without `c - (c - a)` losing `a` if tiny
fn hyp2f1_transform2_euler_generator(a: f64, b: f64, c: f64, z: C64) -> impl FnMut() -> C64 {
    let factor1 = (-z).powf(a - c) * four_gammas(c, a - b, c - b, a);
    let factor2 = (-z).powf(b - c) * four_gammas(c, b - a, c - a, b);
    let mut generator1 = hypergeometric_series_generator(c - a, 1.0 - a, b - a + 1.0, 1.0 / z);
    let mut generator2 = hypergeometric_series_generator(c - b, 1.0 - b, a - b + 1.0, 1.0 / z);
    move || factor1 * generator1() + factor2 * generator2()
}

/// The estimated cancellation `abs_sum / |sum|`, or infinity if NaN (e.g. if `sum` is NaN)
fn cancellation(abs_sum: f64, sum: C64) -> f64 {
    let ratio = abs_sum / sum.norm();
    if ratio.is_nan() { f64::INFINITY } else { ratio }
}

/// `series_eval` of the generator, and its estimated cancellation sum(|term|) / |sum|
fn series_eval_cancellation(mut g: impl FnMut() -> C64) -> (C64, f64) {
    let mut abs_sum = 0.0;
    let g = || {
        let term = g();
        abs_sum += term.norm();
        term
    };
    let result = series_eval(g);
    (result, cancellation(abs_sum, result))
}

/// `series_eval_cancellation`, and the estimated relative size of the omitted tail, a geometric
/// series in the ratio of the last two terms (infinite unless that ratio is below one)
fn series_eval_tail(mut g: impl FnMut() -> C64) -> (C64, f64, f64) {
    let (mut previous, mut last) = (0.0, 0.0);
    let (sum, cancellation) = series_eval_cancellation(|| {
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
    (sum, cancellation, tail / sum.norm())
}

/// `series_eval_cancellation`, but converged only after three successive small terms
//
// One term of the limit or paired series can vanish with its logarithmic coefficient.
fn series_eval_settled(mut g: impl FnMut() -> C64) -> (C64, f64) {
    let (mut sum, mut abs_sum, mut small) = (C64::ZERO, 0.0, 0);
    for _ in 0..HYP2F1_MAXITER {
        let term = g();
        sum += term;
        abs_sum += term.norm();
        small = if term.norm() < sum.norm() * HYP2F1_EPS {
            small + 1
        } else {
            0
        };
        if small == 3 {
            return (sum, cancellation(abs_sum, sum));
        }
    }
    (C64::NAN, f64::INFINITY)
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
    z: C64,
) -> impl FnMut() -> C64 {
    /* 1/z transform in limit as a - b approaches a non-negative integer m. (Can swap a and b to
     * handle the m a negative integer case. */
    let mut d1 = digamma(1.0);
    let mut d2 = digamma(1.0 + m);
    let mut d3 = digamma(a);
    let mut d4 = digamma(c - a);
    let log_neg_z = (-z).ln();
    let mut factor = c64(poch(b, m) * poch(1.0 - c + b, m) * rgamma(m + 1.0), 0.0);
    let mut k: u64 = 0;
    move || {
        let term = (d1 + d2 - d3 - d4 + log_neg_z) * factor;
        let k_ = k as f64;
        // Use digamma(x + 1) = digamma(x) + 1/x
        d1 += 1.0 / (1.0 + k_); // d1 = digamma(1 + k)
        d2 += 1.0 / (1.0 + m + k_); // d2 = digamma(1 + m + k)
        d3 += 1.0 / (a + k_); // d3 = digamma(a + k)
        d4 -= 1.0 / (c - a - k_ - 1.0); // d4 = digamma(c - a - k)
        let f = (b + m + k_) * (1.0 - c + b + m + k_) / ((k + 1) as f64 * (m + k_ + 1.0));
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
    z: C64,
) -> impl FnMut() -> C64 {
    /* 1/z transform in limit as a - b approaches a non-negative integer m, and c - a approaches
     * a positive integer n. */
    let mut d1 = digamma(1.0);
    let mut d2 = digamma(1.0 + m);
    let mut d3 = digamma(a);
    let mut d4 = digamma(n);
    let log_neg_z = (-z).ln();
    let mut factor = c64(poch(b, m) * poch(1.0 - c + b, m) * rgamma(m + 1.0), 0.0);
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
            let f = (-1.0_f64).powf(m + n)
                * crate::binom(c - 1.0, b - 1.0)
                * poch(c - a + 1.0, m - 1.0);
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
    let mut term = c64(gamma(m) * rgamma(c - b), 0.0);
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

fn lopez_temme_series(a: f64, b: f64, c: f64, z: C64) -> C64 {
    // Unlike in xsf, sum the even series (DLMF 15.8.13) for a symmetric
    // positive beta weight. The odd terms vanish when c = 2b.
    if c.is_finite() && c > 0.0 && c == 2.0 * b {
        let w = z / (2.0 - z);
        let half_a = 0.5 * a;
        let terms = hypergeometric_series_generator(half_a, half_a + 0.5, b + 0.5, w * w);
        let (sum, cancellation) = series_eval_cancellation(terms);
        // Large negative powers after the Euler transform can still cancel.
        if cancellation <= 16.0 {
            return sum;
        }
    }
    // Unlike in xsf, restart with the single-term rule if the paired sum fails.
    lopez_temme_paired(a, b, c, z, HYP2F1_MAXITER)
        .map_or_else(|| series_eval(lopez_temme_series_generator(a, b, c, z)), |(sum, _)| sum)
}

/// The Lopez-Temme series and its estimated cancellation, once a tail bound holds for a pair of
/// terms within `max_terms` terms, or `None` (always unless a, b > 0 and b < c)
fn lopez_temme_paired(a: f64, b: f64, c: f64, z: C64, max_terms: u64) -> Option<(C64, f64)> {
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
                return Some((sum, abs_sum / sum_abs));
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
fn hyp2f1_transform1_limiting_case(a: f64, b: f64, c: f64, z: C64, settled: bool) -> (C64, f64) {
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
    let (series, series_cancellation) = if settled || recovered {
        series_eval_settled(series_generator)
    } else {
        series_eval_cancellation(series_generator)
    };
    let abs_sum = result.norm() + prefactor.norm() * series.norm() * series_cancellation;
    result += prefactor * series;
    let lost = cancellation(abs_sum, result);
    // Unlike in xsf, the sum with a cancelled ratio is NaN if it loses two thirds of the digits.
    if recovered && f64::EPSILON * lost >= f64::EPSILON.cbrt() {
        result = C64::NAN;
    }
    (result, lost)
}

fn hyp2f1_transform2_limiting_case(a: f64, b: f64, c: f64, m: f64, z: C64) -> C64 {
    /* 1 / z transform in limiting case where a - b approaches a non-negative integer m. Negative integer case
     * can be handled by swapping a and b. */
    let series_generator1 = hyp2f1_transform2_limit_finite_part_generator(b, c, m, z);
    let mut result = (-z).powf(-b) * (gamma(c) * rgamma(a));
    result *= series_eval_finite(series_generator1, m as u64);
    let prefactor = (-z).powf(-a) * (rgamma(a) * rgamma(c - b)) * gamma(c);
    let n = c - a;
    let n_int = n.round();
    if (n - n_int).abs() < HYP2F1_EPS {
        // Unlike in xsf, pass an integer at the positive-parameter coefficient switch.
        let n = if n_int > 0.0 && b > 0.0 { n_int } else { n };
        let series_generator2 =
            hyp2f1_transform2_limit_series_cminus_a_int_generator(a, b, c, m, n, z);
        result += prefactor * series_eval(series_generator2);
        return result;
    }
    let series_generator2 = hyp2f1_transform2_limit_series_generator(a, b, c, m, z);
    result += prefactor * series_eval(series_generator2);
    result
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
    let (w3, w3_cancellation, w3_tail) =
        series_eval_tail(hypergeometric_series_generator(a, b, a + b + 1.0 - c, w));
    if !(w3.is_finite() && w3_tail.is_finite()) {
        return None;
    }
    let (f2, f2_cancellation, f2_tail) =
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
        g * w3
    } else {
        (w3.ln() + g1.l.hi).exp()
    };
    let (ln_z, ln_w) = ((1.0 - c) * z.ln(), (c - a - b) * w.ln());
    let t2 = (ln_z + ln_w + ln_g2).exp() * sin.signum() * f2;
    // (with the tail of each series, and the rounding of the exponent of t2)
    let exponent_cancellation = ln_z.norm() + ln_w.norm() + ln_g2.abs();
    let t1_error = f64::EPSILON * w3_cancellation + w3_tail + g1.err;
    let t2_error = f64::EPSILON * (f2_cancellation + exponent_cancellation) + f2_tail + g2.err;
    Some((t1 + t2, t1.norm() * t1_error + t2.norm() * t2_error))
}

fn hyp2f1_pfaff(a: f64, b: f64, c: f64, z: C64) -> C64 {
    let (a, b) = if 0.0 < b && b < a && a < c {
        (b, a)
    } else {
        (a, b)
    };
    let series_generator = hypergeometric_series_generator(a, c - b, c, z / (z - 1.0));
    let positive_parameters = a > 0.0 && a.is_finite() && b > 0.0 && b < c && c.is_finite();
    if z.is_real() || !positive_parameters {
        return (1.0 - z).powf(-a) * series_eval(series_generator);
    }
    // Unlike in xsf, evaluate the Pfaff prefactor without rounding |1-z| first.
    let log_modulus = 0.5 * z.re.mul_add(z.re - 2.0, z.im * z.im).ln_1p();
    let (sin, cos) = (a * z.im.atan2(1.0 - z.re)).sin_cos();
    let power = c64(cos, sin) * (-a * log_modulus).exp();
    // Unlike in xsf, sum this positive-parameter Pfaff series to machine epsilon.
    power * tools::series_eval(series_generator, C64::ZERO, f64::EPSILON, HYP2F1_MAXITER)
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
fn hyp2f1_transform1_near_balance(a: f64, b: f64, c: f64, z: C64) -> Option<(C64, f64)> {
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
    // The phase of w^d must not underflow before the prefactor restores its scale.
    if z.im != 0.0 && log_power.im.abs() < f64::MIN_POSITIVE {
        return None;
    }
    let log_ratio = crate::lgam1p(d) - crate::lgam1p(-d) - digamma_integral(a, d);
    let mut log_ratio = log_ratio - digamma_integral(b, d); // log(A_k / B_k)
    if !log_ratio.is_finite() {
        return None;
    }
    let mut upper = C64::ONE; // B_k w^k / B_0
    let mut k = 0.0;
    let (series, cancellation) = series_eval_settled(|| {
        let term = upper * crate::xsf::cexpm1(log_ratio - log_power);
        let (p, q, n) = (a + k, b + k, k + 1.0);
        let (shifted, nd) = ((p + d) * (q + d), n - d);
        // A_k+1 / B_k+1 = (A_k / B_k) (1 + d step)
        let step = (2.0 * p * q - (p + q + d) * nd) / (nd * shifted);
        log_ratio += (d * step).ln_1p();
        upper *= w * (shifted / (n * (n + d)));
        k = n;
        term
    });
    let prefactor = log_power.exp() * -four_gammas(c, -d, a, b);
    let value = series * prefactor;
    // An underflowed prefactor loses the scale, and a product whose imaginary part cancels by
    // more than 1 / |d| loses more than the unpaired series.
    let rotated = (series.re * prefactor.im).abs() + (series.im * prefactor.re).abs();
    let accurate = prefactor.norm() >= f64::MIN_POSITIVE && value.im.abs() >= d.abs() * rotated;
    (accurate && value.is_finite()).then_some((value, cancellation))
}

/// The limit of the 1 - z series for an integer balance, if it is closer than `result`.
#[cold]
fn hyp2f1_transform1_near_limit(a: f64, b: f64, c: f64, z: C64, offset: f64, result: C64) -> C64 {
    // At the integer balance, which its gamma ratios assume.
    let (limit, cancellation) =
        hyp2f1_transform1_limiting_case(a, b, a + b + (c - a - b).round(), z, true);
    // A heuristic for the relative derivative of the terms of the limit in the balance, from
    // their logarithm and digammas (DLMF 15.8.10). The sum amplifies it by its cancellation.
    let sensitivity =
        (1.0 - z).ln().norm() + [a, b, c].map(|x| digamma(x).abs()).iter().sum::<f64>();
    let error = (offset * sensitivity + f64::EPSILON) * cancellation;
    // The margin covers the heuristic estimate.
    let agree = (result - limit).norm() <= 32.0 * error * limit.norm();
    if error < f64::EPSILON.cbrt() && !agree {
        limit
    } else {
        result
    }
}

fn hyp2f1_near_one(a: f64, b: f64, c: f64, z: C64) -> C64 {
    if let Some(value) = hyp2f1_incomplete_beta(a, b, c, z) {
        return value;
    }
    let balance = c - a - b;
    let offset = (balance - balance.round()).abs();
    let (mut result, cancellation) = if offset < HYP2F1_EPS {
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
    if c > 0.0 && cancellation * f64::EPSILON > HYP2F1_EPS {
        let series_generator = hypergeometric_series_generator(a, b, c, z);
        let (maclaurin, maclaurin_cancellation) = series_eval_cancellation(series_generator);
        if maclaurin_cancellation < cancellation {
            result = maclaurin;
        }
    }
    // Unlike in xsf, check the result against the limit for a nearly integer balance.
    // Rounding in c - a - b moves the poles that cancel between the two 1 - z series by a
    // relative eps / offset, for an error near eps / offset^2, where the limit is off by its
    // sensitivity to the balance: the two meet at an offset near cbrt(eps).
    // (not for a balance beyond the term limit, which the limit sums as many terms for)
    let near = offset >= HYP2F1_EPS && offset.powi(3) < f64::EPSILON;
    if near && balance.abs() <= HYP2F1_MAXITER as f64 {
        return hyp2f1_transform1_near_limit(a, b, c, z, offset, result);
    }
    result
}

fn hyp2f1_inverse(a: f64, b: f64, c: f64, z: C64) -> C64 {
    // Unlike in xsf (scipy/xsf#387), also try it on the Euler transform (DLMF 15.8.1) if that
    // reduces a and b, as for |z| < 0.9.
    let euler = c - a < a && c - b < b;
    // times (1 - z)^(c - a - b), if both are normal
    let euler_scaled = move |result: C64| {
        let prefactor = (1.0 - z).powf(c - a - b);
        (prefactor.norm().is_normal() && result.norm().is_normal()).then(|| prefactor * result)
    };
    if (a - b - (a - b).round()).abs() < HYP2F1_EPS {
        let (a, b) = if b > a { (b, a) } else { (a, b) };
        let m = (a - b).round();
        // (only if c - (c - a) and c - (c - b) are exact, and b is an integer or not near one)
        let b_near_int = b != b.round() && (b - b.round()).abs() < HYP2F1_EPS;
        if euler && c - (c - a) == a && c - (c - b) == b && !b_near_int {
            let result = hyp2f1_transform2_limiting_case(c - b, c - a, c, m, z);
            if let Some(result) = euler_scaled(result) {
                return result;
            }
        }
        return hyp2f1_transform2_limiting_case(a, b, c, m, z);
    }
    let series_generator = hyp2f1_transform2_generator(a, b, c, z);
    let (mut result, mut cancellation) = series_eval_cancellation(series_generator);
    if euler {
        // use the one that cancels less
        let series_generator = hyp2f1_transform2_euler_generator(a, b, c, z);
        let (euler_result, euler_cancellation) = series_eval_cancellation(series_generator);
        if euler_cancellation < cancellation {
            if let Some(euler_result) = euler_scaled(euler_result) {
                (result, cancellation) = (euler_result, euler_cancellation);
            }
        }
    }
    // Unlike in xsf, if the 1/z series lose more than HYP2F1_EPS to cancellation, also sum the
    // Lopez-Temme series for Re z < 1 where its tail is bounded.
    if z.re < 1.0 && cancellation * f64::EPSILON > HYP2F1_EPS {
        if let Some(value) = hyp2f1_lopez_temme_retry(a, b, c, z, cancellation) {
            return value;
        }
    }
    result
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
    let (sum, lt_cancellation) = lopez_temme_paired(p, q, c, z, LOPEZ_TEMME_RETRY_MAXITER)?;
    let result = prefactor * sum;
    let value = euler.map_or(result, |euler| euler * result);
    let roundoff = f64::EPSILON * lt_cancellation;
    let accurate = 2.0 * lt_cancellation < cancellation && roundoff < f64::EPSILON.sqrt();
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
            let estimate = polynomial::terminating(c - a, c - b, c, z);
            let polynomial = estimate.value;
            let value = prefactor * polynomial;
            // Unlike in xsf, carry the scale of a prefactor or product that under- or overflows.
            let peak = prefactor.re.abs().max(prefactor.im.abs());
            let in_range = (f64::MIN_POSITIVE..=f64::MAX).contains(&peak);
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
            let mut series_generator = hypergeometric_series_generator(c - a, c - b, c, z);
            // Unlike in xsf, track cancellation and try a certified recurrence down in c, or else
            // `hyp2f1_negative_c`.
            if a > 0.0 && b > 0.0 && c < 0.0 {
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
            return prefactor * series_eval(series_generator);
        }
        let series_generator = hypergeometric_series_generator(a, b, c, z);
        return series_eval(series_generator);
    }
    /* Points near exp(iπ/3), exp(-iπ/3) not handled by any of the standard
     * transformations. Use series of López and Temme [5]. These regions
     * were not correctly handled by Zhang and Jin's implementation.
     * -------------------------------------------------------------------------*/
    if (0.9..1.1).contains(&z_abs) && (1.0 - z).norm() >= 0.9 && z.re >= 0.0 {
        let (p, q, prefactor, euler) = lopez_temme_form(a, b, c, z);
        let result = prefactor * lopez_temme_series(p, q, c, z);
        return euler.map_or(result, |euler| euler * result);
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
