//! Translated from scipy/special/_agm.pxd

use core::f64::consts::FRAC_PI_4;

const SQRT_HALF_MAX: f64 = 9.480_751_908_109_176e+153; // sqrt(np.finfo(float).max / 2.0)
const INV_SQRT_HALF_MAX: f64 = 1.054_768_661_486_3e-154; // 1.0 / SQRT_HALF_MAX;

/// Arithmetic-geometric mean, iterative implementation
/// a and b must be positive (not zero, not nan).
#[inline]
fn agm_iter(a: f64, b: f64) -> f64 {
    let (mut a, mut b) = (a, b);
    let mut count = 20;
    let mut amean = 0.5 * a + 0.5 * b;
    #[allow(clippy::float_cmp)]
    while count > 0 && amean != a && amean != b {
        (a, b) = (amean, a.sqrt() * b.sqrt());
        amean = 0.5 * a + 0.5 * b;
        count -= 1;
    }

    // return amean
    amean
}

/// Arithmetic-geometric mean
///
/// Corresponds to [`scipy.special.agm`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.agm.html
#[must_use]
#[inline]
pub fn agm(a: f64, b: f64) -> f64 {
    #[allow(clippy::float_cmp)]
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) {
        // a and b have opposite sign
        f64::NAN
    } else if (a.is_infinite() || b.is_infinite()) && (a == 0.0 || b == 0.0) {
        // one value is inf and the other is 0
        f64::NAN
    } else if a == 0.0 || b == 0.0 {
        // at least one of the arguments is 0
        0.0
    } else if a == b {
        a
    } else {
        let (sgn, a, b) = if a < 0.0 { (-1.0, -a, -b) } else { (1.0, a, b) };

        // At this point, a and b are both positive and not nan.

        if INV_SQRT_HALF_MAX < a && a < SQRT_HALF_MAX && INV_SQRT_HALF_MAX < b && b < SQRT_HALF_MAX
        {
            let e = 4.0 * a * b / (a + b).powi(2);
            sgn * FRAC_PI_4 * (a + b) / crate::ellipkm1(e)
        } else {
            // At least one value is "extreme" (very big or very small).
            // Use the iteration to avoid overflow or underflow.
            sgn * agm_iter(a, b)
        }
    }
}
