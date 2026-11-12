//! Translated into pure Rust from `xsf/cephes/poch.h` (scipy/xsf@09c93c7).
//!
//! Pochhammer symbol (a)_m = gamma(a + m) / gamma(a)

use core::f64::consts::FRAC_1_PI;

use crate::xsf::cephes::dd_real::two_sum;
use crate::xsf::cephes::{gammasgn, lgam, lgam_ratio};
use crate::xsf::trig::cephes_sinpi;

/// `xsf::cephes::detail::is_nonpos_int`, which is false for `|x| >= 1e13`
fn is_small_nonpos_int(x: f64) -> bool {
    x <= 0.0 && x == x.ceil() && x.abs() < 1e13
}

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
        // lgam(a + m) - lgam(a) (below)
        return r * lgam_ratio(a, m).exp();
    }

    /* Check for infinity */
    if is_small_nonpos_int(a + m) && !is_small_nonpos_int(a) && a + m != m {
        return f64::INFINITY;
    }

    /* Check for zero */
    if !is_small_nonpos_int(a + m) && is_small_nonpos_int(a) {
        return 0.0;
    }

    if a < -10.0 && m.abs() <= 1.0 {
        // Unlike in xsf (scipy/xsf#357), the reflection formula avoids the cancellation in
        // lgam(a + m) - lgam(a), with the sines of the exact fractional part of `a`.
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
pub fn pow_falling(x: f64, m: f64) -> f64 {
    // `poch(a, m)` amplifies the rounding error of `a` near the poles of Gamma(a) and Gamma(x + 1),
    // so it is only the last resort
    let a = x - m + 1.0;
    // Reflection (DLMF 5.5.3), where the sine of x - m isn't affected by the rounding of `a`:
    // -x poch(1 - x, m - 1) sin(pi (x - m)) / sin(pi x)
    let reflected = (a < 0.5 && x != x.floor())
        .then(|| poch(1.0 - x, m - 1.0))
        .filter(|p| p.is_normal());
    if let Some(p) = reflected {
        // -x / sin(pi x), which is -1 / pi if sin(pi x) is subnormal
        let s = cephes_sinpi(x);
        let r = if s.is_normal() { -x / s } else { -FRAC_1_PI };
        p * (sinpi_sub(x, m) * r)
    } else {
        // Gamma(x + 1) / Gamma(x + 1 - m), where x + 1 is exact near the poles of Gamma(x + 1)
        let q = poch(x + 1.0, -m);
        if q.is_normal() { 1.0 / q } else { poch(a, m) }
    }
}

/// sin(pi (x - m)), without the rounding error of x - m near an integer
fn sinpi_sub(x: f64, m: f64) -> f64 {
    let (s, e) = two_sum(x, -m);
    let n = s.round();
    let t = cephes_sinpi((s - n) + e);
    if n % 2.0 == 0.0 { t } else { -t }
}
