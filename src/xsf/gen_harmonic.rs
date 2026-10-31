//! Translated into pure Rust from `xsf/gen_harmonic.h` (xsf main).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! In xsf, the integer arguments are templated on `T`, which is either `long long` or `double`,
//! where the *values* are expected to be positive integers less than 2**53. Here they are `u64`.

#![allow(clippy::cast_precision_loss)]

use crate::xsf::cephes::zeta;

//
// If zeta(a, n + 1) / zeta(a, 1) exceeds zeta_ratio_threshold, the two
// values are close enough that the loss of precision in the subtraction
// zeta(a, 1) - zeta(a, n + 1) should be avoided by using the direct
// sum of powers instead.  The most conservative value would 0.5, but
// experimentation shows that 0.9 maintains a relative error of less
// than 6e-15.
//
// (Measured against mpmath on ~4400 sampled inputs with at most 2e5 terms in the direct
// sums, indices below 2**53, and results above 1e-300: max relative error 5.5e-15 for
// gen_harmonic, and 6.1e-14 for normalized_gen_harmonic when its numerator is a direct sum
// of 6e4 nearly equal terms. Results near the underflow threshold are less accurate, e.g.
// 2.6e-12 for normalized_gen_harmonic(74, 133, 439, 163.53) = 2.4e-306.)
//
const ZETA_RATIO_THRESHOLD: f64 = 0.9;

/// `(n + 1) as f64` without overflowing at `u64::MAX` (which rounds to 2**64 either way)
#[inline]
fn succ_f64(n: u64) -> f64 {
    // Unlike in xsf, avoid integer overflow in n + 1 and k + 1 (UB for LLONG_MAX in C++)
    n.saturating_add(1) as f64
}

//
// Compute sum_{i=m}^{n} i**-a.
//
// This function assumes 1 <= m <= n and `a` is finite.
//
fn sum_powers(m: u64, n: u64, a: f64) -> f64 {
    let mut sum = 0.0;
    if a >= 0.0 {
        for i in (m..=n).rev() {
            sum += (i as f64).powf(-a);
        }
    } else {
        for i in m..=n {
            sum += (i as f64).powf(-a);
        }
    }
    sum
}

//
// Computes
//
//      sum_{i=j}^{k} i**-a
//      -------------------
//      sum_{i=1}^{n} i**-a
//
// Requires 1 <= j <= k <= n and finite a.
//
fn normalized_sum_powers(j: u64, k: u64, n: u64, a: f64) -> f64 {
    let mut numer = 0.0;
    let mut denom = 0.0;
    if a >= 0.0 {
        for i in (1..=n).rev() {
            let term = (i as f64).powf(-a);
            denom += term;
            if i >= j && i <= k {
                numer += term;
            }
        }
    } else {
        for i in 1..=n {
            let term = (i as f64).powf(-a);
            denom += term;
            if i >= j && i <= k {
                numer += term;
            }
        }
    }
    numer / denom
}

/// Generalized harmonic number $H_{n,a} = \sum_{i=1}^n i^{-a}$
///
/// If $a > 1$ and not too close to $1$, the sum is computed using the Hurwitz zeta function as
/// $\zeta(a, 1) - \zeta(a, n + 1)$. Otherwise it is computed as a simple sum of the powers.
///
/// Returns NaN if $n = 0$.
///
/// Corresponds to `xsf::gen_harmonic` in xsf main, which is used by
/// [`scipy.stats.zipfian`][zipfian].
///
/// [zipfian]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.zipfian.html
///
/// # See also
/// - [`normalized_gen_harmonic`]: Ratio of partial sums of $i^{-a}$
/// - [`zeta`](crate::zeta): Hurwitz zeta function
#[must_use]
#[allow(clippy::float_cmp)]
pub fn gen_harmonic(n: u64, a: f64) -> f64 {
    if n < 1 {
        return f64::NAN;
    }
    if n == 1 {
        // IEEE: pow(1.0, _) is 1.0.
        return 1.0;
    }
    if a.is_nan() {
        return f64::NAN;
    }
    if a.is_infinite() {
        return if a > 0.0 {
            // a = +inf
            1.0
        } else {
            // a = -inf
            f64::INFINITY
        };
    }
    if a == 0.0 {
        return n as f64;
    }
    if a <= 1.0 {
        sum_powers(1, n, a)
    } else {
        // If here, we know a is finite and a > 1, and n > 1.
        //
        // For a > 1, we can use the formula zeta(a, 1) - zeta(a, n + 1),
        // where zeta(a, k) is the Hurwitz zeta function.
        // But if zeta(a, 1) and zeta(a, n + 1) are close, precision is lost
        // in the subtraction, so we use the explicit sum instead. We consider
        // the values "close" if zeta(a, n + 1)/zeta(a, 1) > zeta_ratio_threshold.
        //
        let z1 = zeta(a, 1.0);
        let znp1 = zeta(a, succ_f64(n));
        if znp1 / z1 <= ZETA_RATIO_THRESHOLD {
            z1 - znp1
        } else {
            sum_powers(1, n, a)
        }
    }
}

/// Normalized generalized harmonic number $\sum_{i=j}^k i^{-a} \big/ \sum_{i=1}^n i^{-a}$
///
/// Requires $1 \le j \le k \le n$; returns NaN if that is not true.
///
/// Corresponds to `xsf::normalized_gen_harmonic` in xsf main, which is used by
/// [`scipy.stats.zipfian`][zipfian].
///
/// [zipfian]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.zipfian.html
///
/// # See also
/// - [`gen_harmonic`]: Generalized harmonic number $H_{n,a}$
#[must_use]
#[allow(clippy::float_cmp)]
pub fn normalized_gen_harmonic(j: u64, k: u64, n: u64, a: f64) -> f64 {
    if j < 1 || k < j || n < k {
        return f64::NAN;
    }
    //
    // Now we know 1 <= j <= k <= n
    //
    if n == 1 {
        // IEEE: pow(1.0, _) is 1.0.
        // n == 1 implies j == k == 1.
        return 1.0;
    }
    if a.is_nan() {
        return f64::NAN;
    }
    if a.is_infinite() {
        return if a > 0.0 {
            // a = +inf
            if j == 1 {
                // Numerator and denominator are both 1.
                1.0
            } else {
                // Numerator is 0, denominator is 1.
                0.0
            }
        } else {
            // a = -inf
            if k == 1 {
                // Numerator is 1, denominator is +inf.
                0.0
            } else {
                // Numerator and denominator are both +inf.
                f64::NAN
            }
        };
    }
    if a == 0.0 {
        return (k - j + 1) as f64 / n as f64;
    }
    if a <= 1.0 {
        normalized_sum_powers(j, k, n, a)
    } else {
        // If here, we know a is finite and a > 1, n > 1, and 1 <= j <= k <= n.
        // See the comments in _gen_harmonic() for an explanation of why
        // we check the ratios of the zeta functions before using them
        // to compute the result.
        let zj = zeta(a, j as f64);
        let zkp1 = zeta(a, succ_f64(k));
        let z1 = zeta(a, 1.0);
        let znp1 = zeta(a, succ_f64(n));
        let zeta_numer_ok = zkp1 / zj < ZETA_RATIO_THRESHOLD;
        let zeta_denom_ok = znp1 / z1 < ZETA_RATIO_THRESHOLD;
        if zeta_numer_ok {
            if zeta_denom_ok {
                // OK to use the zeta formula for the numerator
                // and denominator.
                (zj - zkp1) / (z1 - znp1)
            } else {
                (zj - zkp1) / sum_powers(1, n, a)
            }
        } else if zeta_denom_ok {
            sum_powers(j, k, a) / (z1 - znp1)
        } else {
            normalized_sum_powers(j, k, n, a)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{gen_harmonic, normalized_gen_harmonic};

    #[test]
    fn test_gen_harmonic() {
        // mpmath: sum(mpf(i)**-a for i in range(1, n + 1))
        let cases: [(u64, f64, f64); 8] = [
            (2, 2.0, 1.25),
            (10, 0.5, 5.020_997_899_292_666_5),
            (100, -1.5, 40_501.224_515_318_936),
            (1000, 1.0, 7.485_470_860_550_345),
            (1000, 2.0, 1.643_934_566_681_56),
            (10, 1.01, 2.902_261_421_419_068_7),
            (1_000_000_000_000_000, 3.0, 1.202_056_903_159_594_2),
            (7, -0.25, 9.606_504_613_320_602),
        ];
        xsref::np_assert_allclose!(
            cases.map(|(n, a, _)| gen_harmonic(n, a)),
            cases.map(|(.., h)| h),
            rtol = 1e-15
        );

        assert!(gen_harmonic(0, 2.0).is_nan());
        assert_eq!(gen_harmonic(1, f64::NAN), 1.0);
        assert!(gen_harmonic(2, f64::NAN).is_nan());
        assert_eq!(gen_harmonic(5, f64::INFINITY), 1.0);
        assert_eq!(gen_harmonic(5, f64::NEG_INFINITY), f64::INFINITY);
        assert_eq!(gen_harmonic(5, 0.0), 5.0);
    }

    #[test]
    fn test_normalized_gen_harmonic() {
        // mpmath: H(j..k, a) / H(1..n, a)
        let cases: [(u64, u64, u64, f64, f64); 7] = [
            (1, 1, 2, 2.0, 0.8),
            (2, 5, 10, 0.5, 0.444_467_552_195_255_47),
            (3, 50, 100, -1.5, 0.178_880_381_185_630_27),
            (10, 100, 1000, 1.0, 0.315_064_918_106_964_85),
            (1, 10, 1000, 2.0, 0.942_718_623_098_786_7),
            (5, 8, 10, 1.01, 0.214_675_893_759_013_25),
            (10, 20, 1_000_000_000_000, 3.0, 0.003_607_030_809_312_731),
        ];
        xsref::np_assert_allclose!(
            cases.map(|(j, k, n, a, _)| normalized_gen_harmonic(j, k, n, a)),
            cases.map(|(.., r)| r),
            rtol = 1e-15
        );

        assert!(normalized_gen_harmonic(0, 1, 1, 2.0).is_nan());
        assert!(normalized_gen_harmonic(2, 1, 3, 2.0).is_nan());
        assert!(normalized_gen_harmonic(1, 3, 2, 2.0).is_nan());
        assert_eq!(normalized_gen_harmonic(1, 1, 1, f64::NAN), 1.0);
        assert_eq!(normalized_gen_harmonic(1, 2, 3, f64::INFINITY), 1.0);
        assert_eq!(normalized_gen_harmonic(2, 2, 3, f64::INFINITY), 0.0);
        assert_eq!(normalized_gen_harmonic(1, 1, 3, f64::NEG_INFINITY), 0.0);
        assert!(normalized_gen_harmonic(1, 2, 3, f64::NEG_INFINITY).is_nan());
        assert_eq!(normalized_gen_harmonic(2, 3, 4, 0.0), 0.5);
    }
}
