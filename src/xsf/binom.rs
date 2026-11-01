//! Translated into pure Rust from `xsf/binom.h` (xsf v0.2.2), which was translated from Cython
//! into C++ by SciPy developers in 2024.
//!
//! Original authors: Pauli Virtanen, Eric Moore

use crate::xsf::cephes::{beta, gamma, lbeta, poch};
use crate::xsf::trig::cephes_sinpi;

use core::f64::consts::PI;

/// Binomial coefficient considered as a function of two real variables
///
/// For real arguments, the binomial coefficient is defined as
///
/// $$
/// \begin{align*}
/// \binom{n}{k}
/// &= {\Gamma(n+1) \over \Gamma(k+1) \\ \Gamma(n-k+1) } \\\\
/// &= {1 \over (n+1) \\ \Beta(k+1,\\, n-k+1)}
/// \end{align*}
/// $$
///
/// Where $\Gamma$ is the Gamma function ([`gamma`](crate::gamma)) and $\Beta$ the Beta function
/// ([`beta`](crate::beta)).
///
/// Corresponds to [`scipy.special.binom`][binom].
///
/// [binom]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.binom.html
///
/// # See also
/// - [`comb`](crate::comb) integer version "n choose k"
#[must_use]
#[inline]
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub fn binom(n: f64, k: f64) -> f64 {
    let mut kx;
    let mut nx;
    let mut num;
    let mut den;

    if n < 0.0 {
        nx = n.floor();
        if n == nx {
            // Undefined
            return f64::NAN;
        }
    }

    kx = k.floor();
    if k == kx && (n.abs() > 1e-8 || n == 0.0) {
        /* Integer case: use multiplication formula for less rounding
         * error for cases where the result is an integer.
         *
         * This cannot be used for small nonzero n due to loss of
         * precision. */
        nx = n.floor();
        if nx == n && kx > nx / 2.0 && nx > 0.0 {
            // Reduce kx by symmetry
            kx = nx - kx;
        }

        if (0.0..20.0).contains(&kx) {
            num = 1.0;
            den = 1.0;
            for i in 1..=kx as i32 {
                num *= f64::from(i) + n - kx;
                den *= f64::from(i);
                if num.abs() > 1e50 {
                    num /= den;
                    den = 1.0;
                }
            }
            return num / den;
        }
    }

    // general case
    if n >= 1e10 * k && k > 0.0 {
        // avoid under/overflows intermediate results
        return (-lbeta(1.0 + n - k, 1.0 + k) - (n + 1.0).ln()).exp();
    }
    if k > 1e8 * n.abs() {
        // avoid loss of precision
        // Unlike in xsf (scipy/xsf#352), the second term includes the factor n + 1.
        num = gamma(1.0 + n) / k.abs() + gamma(1.0 + n) * n / (2.0 * k * k) * (n + 1.0); // + ...
        let den = PI * k.abs().powf(n);
        if den.is_normal() {
            num /= den;
        } else {
            // Unlike in xsf (scipy/xsf#352), |k|^n is split in two if `den` over- or underflows
            let kn2 = k.abs().powf(0.5 * n);
            num = num / (PI * kn2) / kn2;
        }
        if k > 0.0 {
            kx = k.floor();
            // Unlike in xsf (scipy/xsf#352), the parity of kx is also determined for kx >= 2^31
            let dk = k - kx;
            let sgn = if kx % 2.0 == 0.0 { 1.0 } else { -1.0 };
            // Unlike in xsf (scipy/xsf#352), sinpi is used instead of sin
            return num * cephes_sinpi(dk - n) * sgn;
        }
        kx = k.floor();
        if kx == f64::from(kx as i32) {
            return 0.0;
        }
        return num * (k * PI).sin();
    }
    if k.is_finite() && k.abs() > 2e4 && n.abs() < 170.0 {
        // Unlike in xsf (scipy/xsf#352), the reflection formula
        // binom(n, k) = Gamma(1 + n) sin(pi (n - k + 1)) / (pi poch(k - n, n + 1))
        // is used for large |k|. The cutoff ensures that the `poch` base exceeds 1e4.
        let (a, s) = if k > 0.0 {
            kx = k.round();
            nx = n.round();
            // sin(pi (n - k + 1)) = (-1)^(kx - nx) sin(pi ((k - kx) - (n - nx)))
            let sgn = if (kx - nx) % 2.0 == 0.0 { 1.0 } else { -1.0 };
            (k - n, cephes_sinpi((k - kx) - (n - nx)) * sgn)
        } else {
            // the reflection of Gamma(1 + k) instead of Gamma(1 + n - k)
            (-k, cephes_sinpi(-k))
        };
        if s == 0.0 {
            return 0.0;
        }
        // poch(a, n + 1) = poch(a, n) (a + n), because the rounding error of n + 1 would be
        // amplified by ~ln(a)
        let p = poch(a, n);
        let q = p * (a + n);
        if p.is_normal() && q.is_normal() {
            return gamma(1.0 + n) / q * (s / PI);
        }
        // poch(a, n) is split in two if it over- or underflows
        let m = 0.5 * n;
        return gamma(1.0 + n) / poch(a, m) * (s / PI) / (poch(a + m, m) * (a + n));
    }
    1.0 / (n + 1.0) / beta(1.0 + n - k, 1.0 + k)
}
