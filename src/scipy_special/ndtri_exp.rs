//! Ported from `scipy/special/_ndtri_exp.pxd`

use core::f64::consts::SQRT_2;

use crate::xsf::cephes::polevl::{p1evl, polevl};
use crate::xsf::cephes::{NDTRI_P1, NDTRI_P2, NDTRI_Q1, NDTRI_Q2};

/// Return inverse of log CDF of normal distribution for very small y
///
/// For p sufficiently small, the inverse of the CDF of the normal
/// distribution can be approximated to high precision as a rational function
/// in sqrt(-2.0 * log(p)).
#[inline]
fn ndtri_exp_small_y(y: f64) -> f64 {
    // sqrt(-2 * y) is faster and has more precision but overflows when y < -DBL_MAX * 0.5
    let x = if y >= -f64::MAX * 0.5 {
        (-2.0 * y).sqrt()
    } else {
        SQRT_2 * (-y).sqrt()
    };
    let x0 = x - x.ln() / x;
    let z = x.recip();
    let x1 = z * if x < 8.0 {
        polevl(z, &NDTRI_P1) / p1evl(z, &NDTRI_Q1)
    } else {
        polevl(z, &NDTRI_P2) / p1evl(z, &NDTRI_Q2)
    };
    x1 - x0
}

/// Inverse of [`log_ndtr`](crate::log_ndtr) vs x.
///
/// Allows for greater precision than [`ndtri`](crate::ndtri) composed with [`f64::exp`] for very
/// small values of `y` and for `y` close to 0.
///
/// Corresponds to [`scipy.special.ndtri_exp`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri_exp.html
///
/// # See also
/// - [`log_ndtr`](crate::log_ndtr): log of the standard normal cumulative distribution function
/// - [`ndtr`](crate::ndtr): standard normal cumulative distribution function
/// - [`ndtri`](crate::ndtri): standard normal percentile function
#[must_use]
#[inline]
pub fn ndtri_exp(y: f64) -> f64 {
    if y < -f64::MAX {
        f64::NEG_INFINITY
    } else if y < -2.0 {
        ndtri_exp_small_y(y)
    } else if y > -0.145_413_457_868_859_06 {
        // y > log1p(-exp(-2))
        -crate::ndtri(-y.exp_m1())
    } else {
        crate::ndtri(y.exp())
    }
}
