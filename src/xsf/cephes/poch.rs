//! Translated into pure Rust from `xsf/cephes/poch.h` (xsf v0.2.2).
//!
//! Pochhammer symbol (a)_m = gamma(a + m) / gamma(a)

use crate::xsf::cephes::{gammasgn, lgam, lgam_ratio};
use crate::xsf::trig::cephes_sinpi;

#[allow(clippy::float_cmp)]
fn is_nonpos_int(x: f64) -> bool {
    x <= 0.0 && x == x.ceil() && x.abs() < 1e13
}

#[allow(clippy::float_cmp)]
pub(crate) fn poch(a: f64, mut m: f64) -> f64 {
    let mut r = 1.0;

    /*
     * 1. Reduce magnitude of `m` to |m| < 1 by using recurrence relations.
     *
     * This may end up in over/underflow, but then the function itself either
     * diverges or goes to zero. In case the remainder goes to the opposite
     * direction, we end up returning 0*INF = NAN, which is OK.
     */

    /* Recurse down */
    while m >= 1.0 {
        if a + m == 1.0 {
            break;
        }
        m -= 1.0;
        r *= a + m;
        if !r.is_finite() || r == 0.0 {
            break;
        }
    }

    /* Recurse up */
    while m <= -1.0 {
        if a + m == 0.0 {
            break;
        }
        r /= a + m;
        m += 1.0;
        if !r.is_finite() || r == 0.0 {
            break;
        }
    }

    /*
     * 2. Evaluate function with reduced `m`
     *
     * Now either `m` is not big, or the `r` product has over/underflown.
     * If so, the function itself does similarly.
     */

    if m == 0.0 {
        /* Easy case */
        return r;
    } else if a > 1e4 && m.abs() <= 1.0 {
        /* Avoid loss of precision */
        return r
            * a.powf(m)
            * (1.0
                + m * (m - 1.0) / (2.0 * a)
                + m * (m - 1.0) * (m - 2.0) * (3.0 * m - 1.0) / (24.0 * a * a)
                + m * m * (m - 1.0) * (m - 1.0) * (m - 2.0) * (m - 3.0) / (48.0 * a * a * a));
    } else if a > 10.0 && m.abs() <= 1.0 {
        // Unlike in xsf (scipy/xsf#357), `lgam_ratio` avoids the cancellation in
        // lgam(a + m) - lgam(a) (below), e.g. a relative error of 8e-12 for
        // poch(5878.8618471964755, -0.029080395984806877)
        return r * lgam_ratio(a, m).exp();
    }

    /* Check for infinity */
    if is_nonpos_int(a + m) && !is_nonpos_int(a) && a + m != m {
        return f64::INFINITY;
    }

    /* Check for zero */
    if !is_nonpos_int(a + m) && is_nonpos_int(a) {
        return 0.0;
    }

    if a < -10.0 && m.abs() <= 1.0 {
        // Unlike in xsf (scipy/xsf#357), the reflection formula avoids the cancellation in
        // lgam(a + m) - lgam(a), e.g. a relative error of 3.6e-4 for poch(-96483701341.52873,
        // -0.0010688260996181021). The sines use the exact fractional part of `a` instead of the
        // rounded `a + m`.
        let f = a - a.floor();
        return r * (cephes_sinpi(f) / cephes_sinpi(f + m)) * poch(1.0 - a - m, m);
    }

    r * (lgam(a + m) - lgam(a)).exp() * gammasgn(a + m) * gammasgn(a)
}

/// Rising factorial $\rpow x m$
///
/// $$\rpow x m = {\Gamma(x+m) \over \Gamma(x)}$$
///
/// Corresponds to [`scipy.special.poch`][poch] in SciPy.
///
/// [poch]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.poch.html
///
/// # See also
/// - [`pow_falling`]: falling factorial $\fpow x m$
/// - [`gamma`](crate::gamma): gamma function $\Gamma(x)$
#[doc(alias = "poch")]
#[must_use]
#[inline]
pub fn pow_rising(x: f64, m: f64) -> f64 {
    poch(x, m)
}

/// Falling factorial $\fpow x m$
///
/// $$\fpow x m = {\Gamma(x+1) \over \Gamma(x-m+1)}$$
///
/// Note that there is no `scipy.special` analogue for this function, but it can be expressed in
/// terms of the rising factorial as `pow_rising(x - m + 1, m)`.
///
/// # See also
/// - [`pow_rising`]: rising factorial $\rpow x m$
/// - [`gamma`](crate::gamma): gamma function $\Gamma(x)$
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
pub fn pow_falling(x: f64, m: f64) -> f64 {
    let a = x - m + 1.0;
    if a < -10.0 && x < 0.0 && m != m.floor() {
        // `poch(a, m)` would reflect, which amplifies the rounding error of `a` near the poles of
        // Gamma(x + 1), e.g. 3e-5 for pow_falling(-59.0000000001, -7.654166228426843)
        let q = poch(x + 1.0, -m);
        if q.is_normal() {
            return 1.0 / q;
        }
    }
    poch(a, m)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_pow_rising() {
        xsref::test("poch", "d_d-d", |x| crate::pow_rising(x[0], x[1]));
    }

    #[test]
    fn test_pow_rising_lgam_ratio() {
        // unlike in xsf (scipy/xsf#357): `lgam_ratio` for 10 < a <= 1e4
        for (a, m, expected) in [
            (
                5_878.861_847_196_475_5,
                -0.029_080_395_984_806_877,
                0.776_941_938_775_515_5,
            ),
            (150.5, 0.75, 42.941_972_832_542_55),
            (50.25, -0.5, 0.142_132_904_591_272_65),
        ] {
            let actual = crate::pow_rising(a, m);
            assert!(
                (actual / expected - 1.0).abs() < 5e-15,
                "{actual:e} != {expected:e}"
            );
        }
    }

    #[test]
    fn test_pow_rising_reflection() {
        // unlike in xsf (scipy/xsf#357): the reflection formula for a < -10
        for (a, m, expected) in [
            (
                -96_483_701_341.528_73,
                -0.001_068_826_099_618_102_1,
                0.973_630_092_687_070_2,
            ),
            (-10_000.25, 0.5, -100.002_499_984_374_61),
            (-50.75, 2.5, 17_673.802_508_037_614),
        ] {
            let actual = crate::pow_rising(a, m);
            assert!(
                (actual / expected - 1.0).abs() < 5e-15,
                "{actual:e} != {expected:e}"
            );
        }
    }

    #[test]
    fn test_pow_falling_near_pole() {
        for (x, m, expected) in [
            (
                -59.000_000_000_1,
                -7.654_166_228_426_843,
                0.000_141_662_338_550_097_54,
            ),
            (-0.999_999_999_9, 20.5, -3.527_763_913_031_247e28),
            (-1_000.25, -104.5, -9.343_131_111_265e-312), // subnormal
        ] {
            let actual = crate::pow_falling(x, m);
            assert!(
                (actual / expected - 1.0).abs() < 1e-12,
                "{actual:e} != {expected:e}"
            );
        }
    }
}
