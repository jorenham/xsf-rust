//! Translated into pure Rust from `xsf/binom.h` (xsf v0.2.2), which was translated from Cython
//! into C++ by SciPy developers in 2024.
//!
//! Original authors: Pauli Virtanen, Eric Moore

use crate::xsf::cephes::{beta, gamma, lbeta};
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
        // Unlike in xsf, the second term includes the factor n + 1 (from the expansion
        // Gamma(k - n) / Gamma(k + 1) ~ k^(-n-1) (1 + n (n + 1) / (2 k) + ...)), without which the
        // relative error is ~n^2 / (2 k).
        num = gamma(1.0 + n) / k.abs() + gamma(1.0 + n) * n / (2.0 * k * k) * (n + 1.0); // + ...
        num /= PI * k.abs().powf(n);
        if k > 0.0 {
            kx = k.floor();
            // Unlike in xsf, the parity of kx is also determined for kx >= 2^31, where the C++
            // `static_cast<int>(kx)` is UB (with GCC on x86-64, sin((k - n) * pi) is then used).
            let dk = k - kx;
            let sgn = if kx % 2.0 == 0.0 { 1.0 } else { -1.0 };
            // Unlike in xsf, sinpi avoids the rounding error of (dk - n) * pi near its zeros
            return num * cephes_sinpi(dk - n) * sgn;
        }
        kx = k.floor();
        if kx == f64::from(kx as i32) {
            return 0.0;
        }
        return num * (k * PI).sin();
    }
    1.0 / (n + 1.0) / beta(1.0 + n - k, 1.0 + k)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_binom_f64() {
        xsref::test("binom", "d_d-d", |x| crate::binom(x[0], x[1]));
    }

    #[test]
    fn test_binom_large_k() {
        // unlike in xsf: the n (n + 1) / (2 k) term, the parity of k >= 2^31, and sinpi
        for (n, k, expected) in [
            (2.5, 300_000_000.7, 1.329_598_785_965_806_2e-30),
            (0.5, 10_000_000_001.25, 1.994_711_401_707_956_7e-16),
        ] {
            let actual = crate::binom(n, k);
            assert!(
                (actual / expected - 1.0).abs() < 1e-14,
                "{actual:e} != {expected:e}"
            );
        }
        assert_eq!(crate::binom(1.0, 2_147_483_647.0), 0.0);
        assert_eq!(crate::binom(169.0, 20_000_000_000.0), 0.0);
    }
}
