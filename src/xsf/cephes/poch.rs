//! Translated into pure Rust from `xsf/cephes/poch.h` (xsf v0.2.2).
//!
//! Pochhammer symbol (a)_m = gamma(a + m) / gamma(a)

use crate::xsf::cephes::{gammasgn, lgam};

#[allow(clippy::float_cmp)]
fn is_nonpos_int(x: f64) -> bool {
    x <= 0.0 && x == x.ceil() && x.abs() < 1e13
}

#[allow(clippy::float_cmp)]
fn poch(a: f64, mut m: f64) -> f64 {
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
    }

    /* Check for infinity */
    if is_nonpos_int(a + m) && !is_nonpos_int(a) && a + m != m {
        return f64::INFINITY;
    }

    /* Check for zero */
    if !is_nonpos_int(a + m) && is_nonpos_int(a) {
        return 0.0;
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
    poch(x - m + 1.0, m)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_pow_rising() {
        xsref::test("poch", "d_d-d", |x| crate::pow_rising(x[0], x[1]));
    }
}
