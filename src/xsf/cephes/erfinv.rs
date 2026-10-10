//! Translated into pure Rust from `xsf/cephes/erfinv.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::cephes::ndtri;
use crate::xsf::cephes::ndtri::ndtri_central;

use core::f64::consts::{FRAC_1_SQRT_2, FRAC_2_SQRT_PI, LN_2};

/// Inverse of the error function [*erf(x)*](crate::erf)
///
/// In the complex domain, there is no unique complex number w satisfying erf(w)=z.
/// This indicates a true inverse function would be multivalued. When the domain restricts to the
/// real, -1 < x < 1, there is a unique real number satisfying erf(erfinv(x))=x.
///
/// Note that unlike [`scipy.special.erfinv`][scipy], which uses the Boost implementation,
/// this function uses the Cephes implementation, which can be less accurate in certain regions.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfinv.html
///
/// # See also
/// - [`erfcinv`]: Inverse of the complementary error function
/// - [`erf`](crate::erf): Error function
/// - [`erfc`](crate::erfc): Error function
#[doc(alias = "erf_inv")]
#[must_use]
#[inline]
pub fn erfinv(y: f64) -> f64 {
    /*
     * Inverse of the error function.
     *
     * Computes the inverse of the error function on the restricted domain
     * -1 < y < 1. This restriction ensures the existence of a unique result
     * such that erf(erfinv(y)) = y.
     */
    const DOMAIN_LB: f64 = -1.0;
    const DOMAIN_UB: f64 = 1.0;

    // Unlike in xsf (1e-7), so that the neglected cubic term stays below 1/4 ulp
    const THRESH: f64 = 1e-8;

    /*
     * For small arguments, use the Taylor expansion
     * erf(y) = 2/\sqrt{\pi} (y - y^3 / 3 + O(y^5)),    y\to 0
     * where we only retain the linear term.
     * Otherwise, y + 1 loses precision for |y| << 1.
     */
    if (-THRESH < y) && (y < THRESH) {
        return y / FRAC_2_SQRT_PI;
    }
    /*
     * Unlike in xsf, avoid the rounding error of y + 1, which loses precision for |y| << 1 and
     * for y -> 1 (e.g. erfinv(1 - 2^-53) would be inf):
     *   ndtri(0.5 * (y + 1)) = -ndtri(0.5 * (1 - y)), and = ndtri_central(0.5 * y) for |y| < 0.5,
     * where 0.5 * y is exact, and y + 1 and 1 - y are exact for y <= -0.5 and y >= 0.5.
     */
    if (-0.5 < y) && (y < 0.5) {
        return ndtri_central(0.5 * y) * FRAC_1_SQRT_2;
    }
    if (DOMAIN_LB < y) && (y < 0.0) {
        ndtri(f64::midpoint(y, 1.0)) * FRAC_1_SQRT_2
    } else if (0.0 < y) && (y < DOMAIN_UB) {
        -ndtri(0.5 * (1.0 - y)) * FRAC_1_SQRT_2
    } else if y == DOMAIN_LB {
        f64::NEG_INFINITY
    } else if y == DOMAIN_UB {
        f64::INFINITY
    } else if y.is_nan() {
        y
    } else {
        f64::NAN
    }
}

/// Inverse of the complementary error function [*erfc(x)*](crate::erfc)
///
/// In the complex domain, there is no unique complex number w satisfying erfc(w)=z.
/// This indicates a true inverse function would be multivalued. When the domain restricts to the
/// real, 0 < x < 2, there is a unique real number satisfying erfc(erfcinv(x))=erfcinv(erfc(x)).
///
/// It is related to inverse of the error function by erfcinv(1 - y) = [erfinv(y)](fn.erfinv.html).
///
/// Note that [`scipy.special.erfcinv`][scipy] also uses the Cephes implementation.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcinv.html
///
/// # See also
/// - [`erfinv`]: Inverse of the error function
/// - [`erf`](crate::erf): Error function
/// - [`erfc`](crate::erfc): Error function
#[doc(alias = "erfc_inv")]
#[must_use]
#[inline]
pub fn erfcinv(y: f64) -> f64 {
    /*
     * Inverse of the complementary error function.
     *
     * Computes the inverse of the complimentary error function on the restricted
     * domain 0 < y < 2. This restriction ensures the existence of a unique result
     * such that erfc(erfcinv(y)) = y.
     */
    const DOMAIN_LB: f64 = 0.0;
    const DOMAIN_UB: f64 = 2.0;

    if (DOMAIN_LB < y) && (y < 2.0 * f64::MIN_POSITIVE) {
        // Unlike in xsf, where 0.5 * y can lose precision (or underflow to 0) for these y:
        // ndtri(0.5 * y) = ndtri_exp(ln(y) - ln(2)).
        return -crate::ndtri_exp(y.ln() - LN_2) * FRAC_1_SQRT_2;
    }
    if (DOMAIN_LB < y) && (y < DOMAIN_UB) {
        -ndtri(0.5 * y) * FRAC_1_SQRT_2
    } else if y == DOMAIN_LB {
        f64::INFINITY
    } else if y == DOMAIN_UB {
        f64::NEG_INFINITY
    } else if y.is_nan() {
        y
    } else {
        f64::NAN
    }
}

#[cfg(test)]
mod tests {
    /// Based on `scipy.special.tests.test_TestInverseErrorFunction.test_literal_values`
    #[test]
    fn test_erfinv_literal_values() {
        let y = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
        let actual = y.map(crate::erfinv);
        let expected = [
            0.0,
            0.088_855_990_494_257_69,
            0.179_143_454_621_291_7,
            0.272_462_714_726_754_3,
            0.370_807_158_593_557_95,
            0.476_936_276_204_469_9,
            0.595_116_081_449_994_8,
            0.732_869_077_959_216_7,
            0.906_193_802_436_823_3,
            1.163_087_153_676_674_3,
        ];
        crate::np_assert_allclose!(&actual, &expected, atol = 1e-15);
    }

    #[test]
    fn test_erfcinv() {
        xsref::test("erfcinv", "d-d", |x| crate::erfcinv(x[0]));
    }
}
