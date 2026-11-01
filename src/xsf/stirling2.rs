//! Translated into pure Rust from `xsf/stirling2.h` (xsf main).
//!
//! Stirling numbers of the second kind
//!
//! SYNOPSIS: Stirling numbers of the second kind count the number of ways to
//!  make a partition of n distinct elements into k non-empty subsets.
//!
//! DESCRIPTION: n is the number of distinct elements of the set to be
//!  partitioned and k is the number of non-empty subsets. Values for n < 0 or
//!  k < 0 are interpreted as 0.
//!
//! NOTE: Original implementation by Lucas Roberts. Ported to xsf by SciPy
//!  developers.
//!
//! References:
//!  - N. M. Temme, "Asymptotic estimates of Stirling numbers", Studies in
//!    Applied Mathematics, 89(3):233-243, 1993.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! This is only called with non-negative integer n and k (through
//! [`stirling2::<f64>`](crate::stirling2())), which rules out the out-of-bounds reads in
//! `stirling2_dp` that xsf has for `NaN` k, or for `0 < k < 1`.

use crate::xsf::binom::binom;
use crate::xsf::lambertw::lambertw;
use num_complex::Complex;

/// Dynamic programming method for small n.
/// Only valid for n <= 219; the stack array is sized accordingly (max 110 elements).
#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::cast_lossless
)]
fn stirling2_dp(n: f64, k: f64) -> f64 {
    // Maximum arraySize for n <= 219 is 110. Guard against misuse.
    const MAX_ARRAY_SIZE: usize = 111;

    if (n == 0.0 && k == 0.0) || (n == 1.0 && k == 1.0) {
        return 1.0;
    }
    if k <= 0.0 || k > n || n < 0.0 {
        return 0.0;
    }

    let ni = n as i32;
    let ki = k as i32;
    let array_size = if ki <= ni - ki + 1 { ki } else { ni - ki + 1 } as usize;

    if array_size > MAX_ARRAY_SIZE {
        return f64::NAN;
    }

    let mut curr = [1.0; MAX_ARRAY_SIZE];

    if ki <= ni - ki + 1 {
        for _ in 1..=ni - ki {
            for j in 1..ki as usize {
                curr[j] = (j + 1) as f64 * curr[j] + curr[j - 1];
                if curr[j].is_infinite() {
                    return f64::INFINITY;
                }
            }
        }
    } else {
        for i in 1..ki {
            for j in 1..(ni - ki + 1) as usize {
                curr[j] += (i + 1) as f64 * curr[j - 1];
                if curr[j].is_infinite() {
                    return f64::INFINITY;
                }
            }
        }
    }
    curr[array_size - 1]
}

/// Second-order Temme asymptotic approximation for large n.
/// Reference: Temme (1993).
#[allow(clippy::float_cmp, clippy::many_single_char_names)]
fn stirling2_temme(n: f64, k: f64) -> f64 {
    if (n == k && n >= 0.0) || (n > 0.0 && k == 1.0) {
        return 1.0;
    }
    if k <= 0.0 || k > n || n < 0.0 {
        return 0.0;
    }

    let mu = k / n;
    let d = (-1.0 / mu).exp() / mu;
    let delta = Complex::new(-d, 0.0);

    // lambertw returns complex; we only need the real part (branch k=0)
    let lwv = lambertw(delta, 0, 1e-8);
    let x0 = lwv.re + 1.0 / mu;
    let t0 = 1.0 / mu - 1.0;
    #[allow(non_snake_case)]
    let F = (t0 / ((1.0 + t0) * (x0 - t0))).sqrt();
    #[allow(non_snake_case)]
    let A = -n * x0.ln() + k * (x0.exp() - 1.0).ln() - k * t0 + (n - k) * t0.ln();

    // F1 correction term (Horner scheme applied to numerator)
    let xt = x0 * t0;
    let t0power3 = t0 * t0 * t0;
    let mut num = -2.0 * x0 * x0 * x0;
    num += ((t0 + 2.0) * t0 + 2.0) * (2.0 * t0power3);
    num +=
        (-6.0 * t0power3 + (8.0 * t0 - 6.0 * x0 - 5.0) * xt + ((2.0 * x0 + 1.0) * x0 + 3.0) * x0)
            * xt;
    let denom = 24.0 * F * (1.0 + t0) * (1.0 + t0) * (x0 - t0) * (x0 - t0) * (x0 - t0) * (x0 - t0);
    #[allow(non_snake_case)]
    let F1 = num / denom;

    A.exp() * k.powf(n - k) * binom(n, k) * (F - F1 / k)
}

/// Stirling numbers of the second kind, using dynamic programming for n <= 219, and the
/// second-order Temme asymptotic approximation for n > 219.
pub(crate) fn stirling2(n: f64, k: f64) -> f64 {
    // Unlike in xsf (scipy/xsf#378), up to n = 219 (not 50), the largest n where no S(n, k)
    // overflows
    if n <= 219.0 {
        return stirling2_dp(n, k);
    }
    stirling2_temme(n, k)
}
