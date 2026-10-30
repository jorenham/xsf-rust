//! The complex `hyp2f1` is translated into pure Rust from `xsf/hyp2f1.h` (xsf v0.2.2).
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

use core::f64::consts::{FRAC_1_PI, PI};

use num_complex::{Complex64, c64};

use crate::xsf::cephes::{LANCZOS_G, gamma, gammasgn, lanczos_sum_expg_scaled, lgam, poch, rgamma};
use crate::xsf::complex::{cdiv, cmul, rsub};
use crate::xsf::log::clog;
use crate::xsf::tools::{series_eval, series_eval_fixed_length};

const HYP2F1_EPS: f64 = 1e-15;
/* The original implementation in SciPy from Zhang and Jin used 1500 for the
 * maximum number of series iterations in some cases and 500 in others.
 * Through the empirical results on the test cases in
 * scipy/special/_precompute/hyp2f1_data.py, it was determined that these values
 * can lead to early termination of series which would have eventually converged
 * at a reasonable level of accuracy. We've bumped the iteration limit to 3000,
 * and may adjust it again based on further analysis. */
const HYP2F1_MAXITER: u64 = 3000;

/// `std::pow(z, y)` for complex `z` and real `y`, as in libstdc++
fn cpow(z: Complex64, y: f64) -> Complex64 {
    if z.im == 0.0 && z.re > 0.0 {
        return c64(z.re.powf(y), 0.0);
    }
    let t = clog(z);
    // std::polar(exp(y * t.real()), y * t.imag())
    let rho = (y * t.re).exp();
    let theta = y * t.im;
    c64(rho * theta.cos(), rho * theta.sin())
}

/// `x / z` for real `x` and complex `z`, as in C++
fn rdiv(x: f64, z: Complex64) -> Complex64 {
    cdiv(c64(x, 0.0), z)
}

fn digamma(x: f64) -> f64 {
    crate::digamma(x)
}

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names
)]
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
    if (u == u.trunc() && u <= 0.0) || (v == v.trunc() && v <= 0.0) {
        /* Return nan if numerator has pole. Diverges to +- infinity
         * depending on direction so value is undefined. */
        return f64::NAN;
    }
    if (w == w.trunc() && w <= 0.0) || (x == x.trunc() && x <= 0.0) {
        // Return 0 if denominator has pole but not numerator.
        return 0.0;
    }

    let mut result = 1.0;
    let (ugh, vgh, wgh, xgh, u_prime, v_prime, w_prime, x_prime);

    if u >= 0.5 {
        result *= lanczos_sum_expg_scaled(u);
        ugh = u + LANCZOS_G - 0.5;
        u_prime = u;
    } else {
        result /= lanczos_sum_expg_scaled(1.0 - u) * (PI * u).sin() * FRAC_1_PI;
        ugh = 0.5 - u + LANCZOS_G;
        u_prime = 1.0 - u;
    }

    if v >= 0.5 {
        result *= lanczos_sum_expg_scaled(v);
        vgh = v + LANCZOS_G - 0.5;
        v_prime = v;
    } else {
        result /= lanczos_sum_expg_scaled(1.0 - v) * (PI * v).sin() * FRAC_1_PI;
        vgh = 0.5 - v + LANCZOS_G;
        v_prime = 1.0 - v;
    }

    if w >= 0.5 {
        result /= lanczos_sum_expg_scaled(w);
        wgh = w + LANCZOS_G - 0.5;
        w_prime = w;
    } else {
        result *= lanczos_sum_expg_scaled(1.0 - w) * (PI * w).sin() * FRAC_1_PI;
        wgh = 0.5 - w + LANCZOS_G;
        w_prime = 1.0 - w;
    }

    if x >= 0.5 {
        result /= lanczos_sum_expg_scaled(x);
        xgh = x + LANCZOS_G - 0.5;
        x_prime = x;
    } else {
        result *= lanczos_sum_expg_scaled(1.0 - x) * (PI * x).sin() * FRAC_1_PI;
        xgh = 0.5 - x + LANCZOS_G;
        x_prime = 1.0 - x;
    }

    if u.abs() >= w.abs() {
        // u has greatest absolute value. Absorb ugh into the others.
        if ((v_prime - u_prime) * (v - 0.5)).abs() < 100.0 * ugh && v > 100.0 {
            /* Special case where base is close to 1. Condition taken from
             * Boost's beta function implementation. */
            result *= ((v - 0.5) * ((v_prime - u_prime) / ugh).ln_1p()).exp();
        } else {
            result *= (vgh / ugh).powf(v - 0.5);
        }

        if ((u_prime - w_prime) * (w - 0.5)).abs() < 100.0 * wgh && u > 100.0 {
            result *= ((w - 0.5) * ((u_prime - w_prime) / wgh).ln_1p()).exp();
        } else {
            result *= (ugh / wgh).powf(w - 0.5);
        }

        if ((u_prime - x_prime) * (x - 0.5)).abs() < 100.0 * xgh && u > 100.0 {
            result *= ((x - 0.5) * ((u_prime - x_prime) / xgh).ln_1p()).exp();
        } else {
            result *= (ugh / xgh).powf(x - 0.5);
        }
    } else {
        // w has greatest absolute value. Absorb wgh into the others.
        if ((u_prime - w_prime) * (u - 0.5)).abs() < 100.0 * wgh && u > 100.0 {
            result *= ((u - 0.5) * ((u_prime - w_prime) / wgh).ln_1p()).exp();
        } else {
            result *= (ugh / wgh).powf(u - 0.5);
        }
        if ((v_prime - w_prime) * (v - 0.5)).abs() < 100.0 * wgh && v > 100.0 {
            result *= ((v - 0.5) * ((v_prime - w_prime) / wgh).ln_1p()).exp();
        } else {
            result *= (vgh / wgh).powf(v - 0.5);
        }
        if ((w_prime - x_prime) * (x - 0.5)).abs() < 100.0 * xgh && x > 100.0 {
            result *= ((x - 0.5) * ((w_prime - x_prime) / xgh).ln_1p()).exp();
        } else {
            result *= (wgh / xgh).powf(x - 0.5);
        }
    }
    // This exhausts all cases because we assume |u| >= |v| and |w| >= |x|.

    result
}

#[allow(clippy::float_cmp, clippy::many_single_char_names)]
fn four_gammas(mut u: f64, mut v: f64, mut w: f64, mut x: f64) -> f64 {
    // Without loss of generality, ensure |u| >= |v| and |w| >= |x|.
    if v.abs() > u.abs() {
        core::mem::swap(&mut u, &mut v);
    }
    if x.abs() > w.abs() {
        core::mem::swap(&mut x, &mut w);
    }
    /* Direct ratio tends to be more accurate for arguments in this range. Range
     * chosen empirically based on the relevant benchmarks in
     * scipy/special/_precompute/hyp2f1_data.py */
    if u.abs() <= 100.0 && v.abs() <= 100.0 && w.abs() <= 100.0 && x.abs() <= 100.0 {
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

#[allow(clippy::cast_precision_loss)]
fn hypergeometric_series_generator(
    a: f64,
    b: f64,
    c: f64,
    z: Complex64,
) -> impl FnMut() -> Complex64 {
    /* Maclaurin series for hyp2f1.
     *
     * Series is convergent for |z| < 1 but is only practical for numerical
     * computation when |z| < 0.9.
     */
    let mut term = c64(1.0, 0.0);
    let mut k: u64 = 0;
    move || {
        let output = term;
        let k_ = k as f64;
        term = cmul(term * (a + k_) * (b + k_) / ((k + 1) as f64 * (c + k_)), z);
        k += 1;
        output
    }
}

fn hyp2f1_transform1_generator(a: f64, b: f64, c: f64, z: Complex64) -> impl FnMut() -> Complex64 {
    /* 1 -z transformation of standard series.*/
    let factor1 = c64(four_gammas(c, c - a - b, c - a, c - b), 0.0);
    let factor2 = cpow(rsub(1.0, z), c - a - b) * four_gammas(c, a + b - c, a, b);
    let mut generator1 = hypergeometric_series_generator(a, b, a + b - c + 1.0, rsub(1.0, z));
    let mut generator2 =
        hypergeometric_series_generator(c - a, c - b, c - a - b + 1.0, rsub(1.0, z));
    move || cmul(factor1, generator1()) + cmul(factor2, generator2())
}

fn hyp2f1_transform1_limit_series_generator(
    a: f64,
    b: f64,
    m: f64,
    z: Complex64,
) -> impl FnMut() -> Complex64 {
    /* 1 - z transform in limit as c - a - b approaches an integer m. */
    let mut d1 = digamma(a);
    let mut d2 = digamma(b);
    let mut d3 = digamma(1.0 + m);
    let mut d4 = digamma(1.0);
    let log_1_z = clog(rsub(1.0, z));
    let mut factor = c64(rgamma(m + 1.0), 0.0);
    let mut k: i32 = 0;
    move || {
        let term = cmul(d1 + d2 - d3 - d4 + log_1_z, factor);
        let k_ = f64::from(k);
        // Use digamma(x + 1) = digamma(x) + 1/x
        d1 += 1.0 / (a + k_); // d1 = digamma(a + k)
        d2 += 1.0 / (b + k_); // d2 = digamma(b + k)
        d3 += 1.0 / (1.0 + m + k_); // d3 = digamma(1 + m + k)
        d4 += 1.0 / (1.0 + k_); // d4 = digamma(1 + k)
        factor = cmul(
            factor,
            rsub(1.0, z) * ((a + k_) * (b + k_) / ((k_ + 1.0) * (m + k_ + 1.0))),
        );
        k += 1;
        term
    }
}

fn hyp2f1_transform2_generator(a: f64, b: f64, c: f64, z: Complex64) -> impl FnMut() -> Complex64 {
    /* 1/z transformation of standard series.*/
    let factor1 = cpow(-z, -a) * four_gammas(c, b - a, b, c - a);
    let factor2 = cpow(-z, -b) * four_gammas(c, a - b, a, c - b);
    let mut generator1 = hypergeometric_series_generator(a, a - c + 1.0, a - b + 1.0, rdiv(1.0, z));
    let mut generator2 = hypergeometric_series_generator(b, b - c + 1.0, b - a + 1.0, rdiv(1.0, z));
    move || cmul(factor1, generator1()) + cmul(factor2, generator2())
}

#[allow(clippy::cast_precision_loss)]
fn hyp2f1_transform2_limit_series_generator(
    a: f64,
    b: f64,
    c: f64,
    m: f64,
    z: Complex64,
) -> impl FnMut() -> Complex64 {
    /* 1/z transform in limit as a - b approaches a non-negative integer m. (Can swap a and b to
     * handle the m a negative integer case. */
    let mut d1 = digamma(1.0);
    let mut d2 = digamma(1.0 + m);
    let mut d3 = digamma(a);
    let mut d4 = digamma(c - a);
    let log_neg_z = clog(-z);
    let mut factor = c64(poch(b, m) * poch(1.0 - c + b, m) * rgamma(m + 1.0), 0.0);
    let mut k: u64 = 0;
    move || {
        let term = cmul(d1 + d2 - d3 - d4 + log_neg_z, factor);
        let k_ = k as f64;
        // Use digamma(x + 1) = digamma(x) + 1/x
        d1 += 1.0 / (1.0 + k_); // d1 = digamma(1 + k)
        d2 += 1.0 / (1.0 + m + k_); // d2 = digamma(1 + m + k)
        d3 += 1.0 / (a + k_); // d3 = digamma(a + k)
        d4 -= 1.0 / (c - a - k_ - 1.0); // d4 = digamma(c - a - k)
        let f = (b + m + k_) * (1.0 - c + b + m + k_) / ((k + 1) as f64 * (m + k_ + 1.0));
        factor = cmul(factor, rdiv(f, z));
        k += 1;
        term
    }
}

#[allow(clippy::cast_precision_loss, clippy::float_cmp)]
fn hyp2f1_transform2_limit_series_cminus_a_int_generator(
    a: f64,
    b: f64,
    c: f64,
    m: f64,
    n: f64,
    z: Complex64,
) -> impl FnMut() -> Complex64 {
    /* 1/z transform in limit as a - b approaches a non-negative integer m, and c - a approaches
     * a positive integer n. */
    let mut d1 = digamma(1.0);
    let mut d2 = digamma(1.0 + m);
    let mut d3 = digamma(a);
    let mut d4 = digamma(n);
    let log_neg_z = clog(-z);
    let mut factor = c64(poch(b, m) * poch(1.0 - c + b, m) * rgamma(m + 1.0), 0.0);
    let mut k: u64 = 0;
    move || {
        let k_ = k as f64;
        if k_ < n {
            let term = cmul(d1 + d2 - d3 - d4 + log_neg_z, factor);
            // Use digamma(x + 1) = digamma(x) + 1/x
            d1 += 1.0 / (1.0 + k_); // d1 = digamma(1 + k)
            d2 += 1.0 / (1.0 + m + k_); // d2 = digamma(1 + m + k)
            d3 += 1.0 / (a + k_); // d3 = digamma(a + k)
            d4 -= 1.0 / (n - k_ - 1.0); // d4 = digamma(c - a - k)
            let f = (b + m + k_) * (1.0 - c + b + m + k_) / ((k + 1) as f64 * (m + k_ + 1.0));
            factor = cmul(factor, rdiv(f, z));
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
            factor = rdiv(f, cpow(z, k_));
        }
        let term = factor;
        let f = (b + m + k_) * (k_ + a - c + 1.0) / ((k + 1) as f64 * (m + k_ + 1.0));
        factor = cmul(factor, rdiv(f, z));
        k += 1;
        term
    }
}

#[allow(clippy::cast_precision_loss)]
fn hyp2f1_transform2_limit_finite_part_generator(
    b: f64,
    c: f64,
    m: f64,
    z: Complex64,
) -> impl FnMut() -> Complex64 {
    /* Initial finite sum in limit as a - b approaches a non-negative integer m. The limiting series
     * for the 1 - z transform also has an initial finite sum, but it is a standard hypergeometric
     * series. */
    let mut term = c64(gamma(m) * rgamma(c - b), 0.0);
    let mut k: u64 = 0;
    move || {
        let output = term;
        let k_ = k as f64;
        term = cdiv(
            term * (b + k_) * (c - b - k_ - 1.0) / ((k + 1) as f64 * (m - k_ - 1.0)),
            z,
        );
        k += 1;
        output
    }
}

#[allow(clippy::cast_precision_loss, non_snake_case)]
fn lopez_temme_series_generator(a: f64, b: f64, c: f64, z: Complex64) -> impl FnMut() -> Complex64 {
    /* Lopez-Temme Series for Gaussian hypergeometric function [4].
     *
     * Converges for all z with real(z) < 1, including in the regions surrounding
     * the points exp(+- i*pi/3) that are not covered by any of the standard
     * transformations.
     */
    let mut n: u64 = 0;
    let mut phi_previous = 1.0;
    let mut phi = 1.0 - 2.0 * b / c;
    let mut Z = cdiv(z * a, z - 2.0);
    move || {
        if n == 0 {
            n += 1;
            return c64(1.0, 0.0);
        }
        if n > 1 {
            // Update phi and Z for n>=2
            let n1 = (n - 1) as f64;
            let new_phi = (n1 * phi_previous - (2.0 * b - c) * phi) / (c + n1);
            phi_previous = phi;
            phi = new_phi;
            Z = cdiv(cmul(Z, z), z - 2.0) * ((a + n1) / n as f64);
        }
        n += 1;
        Z * phi
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn hyp2f1_transform1_limiting_case(a: f64, b: f64, c: f64, m: f64, z: Complex64) -> Complex64 {
    /* 1 - z transform in limiting case where c - a - b approaches an integer m. */
    let mut result = Complex64::ZERO;
    if m >= 0.0 {
        if m != 0.0 {
            let series_generator = hypergeometric_series_generator(a, b, 1.0 - m, rsub(1.0, z));
            result += series_eval_fixed_length(series_generator, Complex64::ZERO, m as u64)
                * four_gammas(m, c, a + m, b + m);
        }
        let prefactor =
            cpow(rsub(1.0, z), m) * ((-1.0_f64).powf(m + 1.0) * gamma(c) / (gamma(a) * gamma(b)));
        let series_generator = hyp2f1_transform1_limit_series_generator(a + m, b + m, m, z);
        result += cmul(
            prefactor,
            series_eval(
                series_generator,
                Complex64::ZERO,
                HYP2F1_EPS,
                HYP2F1_MAXITER,
            ),
        );
        result
    } else {
        result = cpow(rsub(1.0, z), m) * four_gammas(-m, c, a, b);
        let series_generator1 =
            hypergeometric_series_generator(a + m, b + m, 1.0 + m, rsub(1.0, z));
        result = cmul(
            result,
            series_eval_fixed_length(series_generator1, Complex64::ZERO, (-m) as u64),
        );
        let prefactor = (-1.0_f64).powf(m + 1.0) * gamma(c) * (rgamma(a + m) * rgamma(b + m));
        let series_generator2 = hyp2f1_transform1_limit_series_generator(a, b, -m, z);
        result += series_eval(
            series_generator2,
            Complex64::ZERO,
            HYP2F1_EPS,
            HYP2F1_MAXITER,
        ) * prefactor;
        result
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn hyp2f1_transform2_limiting_case(a: f64, b: f64, c: f64, m: f64, z: Complex64) -> Complex64 {
    /* 1 / z transform in limiting case where a - b approaches a non-negative integer m. Negative integer case
     * can be handled by swapping a and b. */
    let series_generator1 = hyp2f1_transform2_limit_finite_part_generator(b, c, m, z);
    let mut result = cpow(-z, -b) * (gamma(c) * rgamma(a));
    result = cmul(
        result,
        series_eval_fixed_length(series_generator1, Complex64::ZERO, m as u64),
    );
    let prefactor = cpow(-z, -a) * (rgamma(a) * rgamma(c - b)) * gamma(c);
    let n = c - a;
    if (n - n.round()).abs() < HYP2F1_EPS {
        let series_generator2 =
            hyp2f1_transform2_limit_series_cminus_a_int_generator(a, b, c, m, n, z);
        result += cmul(
            prefactor,
            series_eval(
                series_generator2,
                Complex64::ZERO,
                HYP2F1_EPS,
                HYP2F1_MAXITER,
            ),
        );
        return result;
    }
    let series_generator2 = hyp2f1_transform2_limit_series_generator(a, b, c, m, z);
    result += cmul(
        prefactor,
        series_eval(
            series_generator2,
            Complex64::ZERO,
            HYP2F1_EPS,
            HYP2F1_MAXITER,
        ),
    );
    result
}

#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]
fn hyp2f1_c64(mut a: f64, mut b: f64, c: f64, z: Complex64) -> Complex64 {
    /* Special Cases
     * nan input */
    if a.is_nan() || b.is_nan() || c.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return c64(f64::NAN, f64::NAN);
    }

    /* Takes constant value 1 when a = 0 or b = 0, even if c is a non-positive
     * integer. This follows mpmath. */
    if a == 0.0 || b == 0.0 {
        return c64(1.0, 0.0);
    }
    let z_abs = z.norm();
    // Equals 1 when z i 0, unless c is 0.
    if z_abs == 0.0 {
        if c != 0.0 {
            return c64(1.0, 0.0);
        }
        // Returning real part NAN and imaginary part 0 follows mpmath.
        return c64(f64::NAN, 0.0);
    }
    let a_neg_int = a == a.trunc() && a < 0.0;
    let b_neg_int = b == b.trunc() && b < 0.0;
    let c_non_pos_int = c == c.trunc() && c <= 0.0;
    /* Diverges when c is a non-positive integer unless a is an integer with
     * c <= a <= 0 or b is an integer with c <= b <= 0, (or z equals 0 with
     * c != 0) Cases z = 0, a = 0, or b = 0 have already been handled. We follow
     * mpmath in handling the degenerate cases where any of a, b, c are
     * non-positive integers. See [3] for a treatment of degenerate cases. */
    if c_non_pos_int && !((a_neg_int && c <= a && a < 0.0) || (b_neg_int && c <= b && b < 0.0)) {
        return c64(f64::INFINITY, 0.0);
    }
    /* Reduces to a polynomial when a or b is a negative integer.
     * If a and b are both negative integers, we take care to terminate
     * the series at a or b of smaller magnitude. This is to ensure proper
     * handling of situations like a < c < b <= 0, a, b, c all non-positive
     * integers, where terminating at a would lead to a term of the form 0 / 0. */
    let max_degree;
    if a_neg_int || b_neg_int {
        if a_neg_int && b_neg_int {
            max_degree = if a > b { a.abs() } else { b.abs() };
        } else if a_neg_int {
            max_degree = a.abs();
        } else {
            max_degree = b.abs();
        }
        if max_degree <= u64::MAX as f64 {
            let series_generator = hypergeometric_series_generator(a, b, c, z);
            return series_eval_fixed_length(
                series_generator,
                Complex64::ZERO,
                (max_degree + 1.0) as u64,
            );
        }
        return c64(f64::NAN, f64::NAN);
    }
    // Kummer's Theorem for z = -1; c = 1 + a - b (DLMF 15.4.26)
    if (z + 1.0).norm() < HYP2F1_EPS && (1.0 + a - b - c).abs() < HYP2F1_EPS && !c_non_pos_int {
        return c64(
            four_gammas(a - b + 1.0, 0.5 * a + 1.0, a + 1.0, 0.5 * a - b + 1.0),
            0.0,
        );
    }
    let mut result;
    let c_minus_a_neg_int = c - a == (c - a).trunc() && c - a < 0.0;
    let c_minus_b_neg_int = c - b == (c - b).trunc() && c - b < 0.0;
    /* If one of c - a or c - b is a negative integer, reduces to evaluating
     * a polynomial through an Euler hypergeometric transformation.
     * (DLMF 15.8.1) */
    if c_minus_a_neg_int || c_minus_b_neg_int {
        let max_degree = if c_minus_b_neg_int {
            (c - b).abs()
        } else {
            (c - a).abs()
        };
        if max_degree <= u64::MAX as f64 {
            result = cpow(rsub(1.0, z), c - a - b);
            let series_generator = hypergeometric_series_generator(c - a, c - b, c, z);
            result = cmul(
                result,
                series_eval_fixed_length(
                    series_generator,
                    Complex64::ZERO,
                    (max_degree + 2.0) as u64,
                ),
            );
            return result;
        }
        return c64(f64::NAN, f64::NAN);
    }
    /* Diverges as real(z) -> 1 when c <= a + b.
     * Todo: Actually check for overflow instead of using a fixed tolerance for
     * all parameter combinations like in the Fortran original. */
    if (1.0 - z.re).abs() < HYP2F1_EPS && z.im == 0.0 && c - a - b <= 0.0 && !c_non_pos_int {
        return c64(f64::INFINITY, 0.0);
    }
    // Gauss's Summation Theorem for z = 1; c - a - b > 0 (DLMF 15.4.20).
    if z.re == 1.0 && z.im == 0.0 && c - a - b > 0.0 && !c_non_pos_int {
        return c64(four_gammas(c, c - a - b, c - a, c - b), 0.0);
    }
    /* |z| < 0, z.real() >= 0. Use the Maclaurin Series.
     * -----------------------------------------------------------------------
     * Apply Euler Hypergeometric Transformation (DLMF 15.8.1) to reduce
     * size of a and b if possible. We follow Zhang and Jin's
     * implementation [1] although there is very likely a better heuristic
     * to determine when this transformation should be applied. As it
     * stands, this hurts precision in some cases. */
    if z_abs < 0.9 && z.re >= 0.0 {
        if c - a < a && c - b < b {
            result = cpow(rsub(1.0, z), c - a - b);
            let series_generator = hypergeometric_series_generator(c - a, c - b, c, z);
            result = cmul(
                result,
                series_eval(
                    series_generator,
                    Complex64::ZERO,
                    HYP2F1_EPS,
                    HYP2F1_MAXITER,
                ),
            );
            return result;
        }
        let series_generator = hypergeometric_series_generator(a, b, c, z);
        return series_eval(
            series_generator,
            Complex64::ZERO,
            HYP2F1_EPS,
            HYP2F1_MAXITER,
        );
    }
    /* Points near exp(iπ/3), exp(-iπ/3) not handled by any of the standard
     * transformations. Use series of López and Temme [5]. These regions
     * were not correctly handled by Zhang and Jin's implementation.
     * -------------------------------------------------------------------------*/
    if (0.9..1.1).contains(&z_abs) && rsub(1.0, z).norm() >= 0.9 && z.re >= 0.0 {
        /* This condition for applying Euler Transformation (DLMF 15.8.1)
         * was determined empirically to work better for this case than that
         * used in Zhang and Jin's implementation for |z| < 0.9,
         *  real(z) >= 0. */
        if (c - a <= a && c - b < b) || (c - a < a && c - b <= b) {
            let series_generator = lopez_temme_series_generator(c - a, c - b, c, z);
            result = cpow(rsub(1.0, z * 0.5), a - c); // Lopez-Temme prefactor
            result = cmul(
                result,
                series_eval(
                    series_generator,
                    Complex64::ZERO,
                    HYP2F1_EPS,
                    HYP2F1_MAXITER,
                ),
            );
            return cmul(cpow(rsub(1.0, z), c - a - b), result); // Euler transform prefactor.
        }
        let series_generator = lopez_temme_series_generator(a, b, c, z);
        result = series_eval(
            series_generator,
            Complex64::ZERO,
            HYP2F1_EPS,
            HYP2F1_MAXITER,
        );
        return cmul(cpow(rsub(1.0, z * 0.5), -a), result); // Lopez-Temme prefactor.
    }
    /* z/(z - 1) transformation (DLMF 15.8.1). Avoids cancellation issues that
     * occur with Maclaurin series for real(z) < 0.
     * -------------------------------------------------------------------------*/
    if z_abs < 1.1 && z.re < 0.0 {
        if 0.0 < b && b < a && a < c {
            core::mem::swap(&mut a, &mut b);
        }
        let series_generator = hypergeometric_series_generator(a, c - b, c, cdiv(z, z - 1.0));
        return cmul(
            cpow(rsub(1.0, z), -a),
            series_eval(
                series_generator,
                Complex64::ZERO,
                HYP2F1_EPS,
                HYP2F1_MAXITER,
            ),
        );
    }
    /* 1 - z transformation (DLMF 15.8.4). */
    if (0.9..1.1).contains(&z_abs) {
        if (c - a - b - (c - a - b).round()).abs() < HYP2F1_EPS {
            // Removable singularity when c - a - b is an integer. Need to use limiting formula.
            let m = (c - a - b).round();
            return hyp2f1_transform1_limiting_case(a, b, c, m, z);
        }
        let series_generator = hyp2f1_transform1_generator(a, b, c, z);
        return series_eval(
            series_generator,
            Complex64::ZERO,
            HYP2F1_EPS,
            HYP2F1_MAXITER,
        );
    }
    /* 1/z transformation (DLMF 15.8.2). */
    if (a - b - (a - b).round()).abs() < HYP2F1_EPS {
        if b > a {
            core::mem::swap(&mut a, &mut b);
        }
        let m = (a - b).round();
        return hyp2f1_transform2_limiting_case(a, b, c, m, z);
    }
    let series_generator = hyp2f1_transform2_generator(a, b, c, z);
    series_eval(
        series_generator,
        Complex64::ZERO,
        HYP2F1_EPS,
        HYP2F1_MAXITER,
    )
}

pub trait Hyp2F1Arg: crate::sealed::Sealed {
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self;
}

impl Hyp2F1Arg for f64 {
    #[inline]
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self {
        crate::xsf::cephes::hyp2f1(a, b, c, self)
    }
}

impl Hyp2F1Arg for num_complex::Complex<f64> {
    #[inline]
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self {
        hyp2f1_c64(a, b, c, self)
    }
}

/// Gauss' hypergeometric function $_2F_1$
///
/// Corresponds to [`scipy.special.hyp2f1`][hyp2f1] in SciPy, and accepts both `f64` and
/// `num_complex::Complex<f64>` inputs for `z`.
///
/// [hyp2f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp2f1.html
///
/// # Notes
///
/// This function is defined for $\abs z < 1$ as
///
/// $$
/// \hyp 2 1 {a_1\enspace a_2} b z
/// = \sum\_{n=0}^\infty
///     {\rpow{a_1}{n} \ \rpow{a_2}{n} \over \rpow b n}
///     {z^n \over n!}
/// ,
/// $$
///
/// and defined on the rest of the complex $z$-plane by analytic continuation [^1].
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
/// When $n$ is a non-negative integer the result is a polynomial of degree $n$.
///
/// The implementation for complex values of $z$ is described in [^2], except for $z$ in the region
/// defined by
///
/// $$
/// 0.9 \le \abs z < 1.1 , \\
/// \abs{1 - z} \ge 0.9 , \\
/// \mathrm{Re}(z) \ge 0
/// $$
///
/// in which the implementation follows [^4].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/15.2>
/// [^2]: S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996
/// [^4]: J.L. Lopez and N.M. Temme, "New series expansions of the Gauss hypergeometric function",
/// Adv Comput Math 39, 349-365 (2013). <https://doi.org/10.1007/s10444-012-9283-y>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function,
///   $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
///
#[must_use]
#[inline]
pub fn hyp2f1<T: Hyp2F1Arg>(a: f64, b: f64, c: f64, z: T) -> T {
    z.hyp2f1(a, b, c)
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_hyp2f1_f64() {
        xsref::test("hyp2f1", "d_d_d_d-d", |x| {
            crate::hyp2f1(x[0], x[1], x[2], x[3])
        });
    }

    #[test]
    fn test_hyp2f1_c64() {
        xsref::test("hyp2f1", "d_d_d_cd-cd", |x| {
            crate::hyp2f1(x[0], x[1], x[2], c64(x[3], x[4]))
        });
    }
}
