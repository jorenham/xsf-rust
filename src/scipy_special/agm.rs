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

/// Arithemtic-geometric mean
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
            sgn * FRAC_PI_4 * (a + b) / unsafe { crate::ffi::xsf::ellipkm1(e) }
        } else {
            // At least one value is "extreme" (very big or very small).
            // Use the iteration to avoid overflow or underflow.
            sgn * agm_iter(a, b)
        }
    }
}

#[cfg(test)]
mod test {
    //! Translated from `scipy.special.tests.test_basic.test_agm_simple`

    use crate::np_assert_allclose;
    use core::f64::consts::SQRT_2;

    #[test]
    #[allow(clippy::too_many_lines)]
    fn test_agm_simple() {
        let rtol = 1e-13;

        // Gauss's constant
        np_assert_allclose!(
            [1.0 / crate::agm(1.0, SQRT_2)],
            [0.834_626_841_674_073_2],
            rtol = rtol
        );

        // These values were computed using Wolfram Alpha, with the
        // function ArithmeticGeometricMean[a, b].
        let agm13 = 1.863_616_783_244_897;
        let agm15 = 2.604_008_190_530_94;
        let agm35 = 3.936_235_503_649_555;
        np_assert_allclose!(
            [
                (1.0, 1.0),
                (1.0, 3.0),
                (1.0, 5.0),
                (3.0, 1.0),
                (3.0, 3.0),
                (3.0, 5.0),
            ]
            .map(|ab| crate::agm(ab.0, ab.1)),
            [1.0, agm13, agm15, agm13, 3.0, agm35],
            rtol = rtol
        );

        // Computed by the iteration formula using mpmath, with mpmath.mp.prec = 1000:
        let agm12 = 1.456_791_031_046_906_8;
        np_assert_allclose!([crate::agm(1.0, 2.0)], [agm12], rtol = rtol);
        np_assert_allclose!([crate::agm(2.0, 1.0)], [agm12], rtol = rtol);
        np_assert_allclose!([crate::agm(-1.0, -2.0)], [-agm12], rtol = rtol);
        np_assert_allclose!(
            [crate::agm(24.0, 6.0)],
            [13.458_171_481_725_614],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(13.0, 123_456_789.5)],
            [11_111_458.498_599_306],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(1e30, 1.0)],
            [2.229_223_055_945_383e+28],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(1e-22, 1.0)],
            [0.030_182_566_420_169_886],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(1e150, 1e180)],
            [2.229_223_055_945_383e+178],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(1e180, 1e-150)],
            [2.063_472_251_016_267_7e+177],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(1e-150, 1e-170)],
            [3.311_261_967_046_375_6e-152],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::MIN_POSITIVE, f64::MAX)],
            [1.989_207_205_001_547_3e+305],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(0.75 * f64::MAX, f64::MAX)],
            [1.564_904_312_298_045e+308],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::MIN_POSITIVE, 3.0 * f64::MIN_POSITIVE)],
            [4.146_684_986_673_500_5e-308],
            rtol = rtol
        );

        // zero, nan and inf cases.
        np_assert_allclose!([crate::agm(0.0, 0.0)], [0.0], rtol = rtol);
        np_assert_allclose!([crate::agm(99.0, 0.0)], [0.0], rtol = rtol);

        np_assert_allclose!([crate::agm(99.0, 0.0)], [0.0], rtol = rtol);
        np_assert_allclose!([crate::agm(0.0, f64::INFINITY)], [f64::NAN], rtol = rtol);
        np_assert_allclose!([crate::agm(f64::INFINITY, 0.0)], [f64::NAN], rtol = rtol);
        np_assert_allclose!(
            [crate::agm(0.0, f64::NEG_INFINITY)],
            [f64::NAN],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::NEG_INFINITY, 0.0)],
            [f64::NAN],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::INFINITY, f64::NEG_INFINITY)],
            [f64::NAN],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::NEG_INFINITY, f64::INFINITY)],
            [f64::NAN],
            rtol = rtol
        );
        np_assert_allclose!([crate::agm(1.0, f64::NAN)], [f64::NAN], rtol = rtol);
        np_assert_allclose!([crate::agm(f64::NAN, -1.0)], [f64::NAN], rtol = rtol);

        np_assert_allclose!(
            [crate::agm(1.0, f64::INFINITY)],
            [f64::INFINITY],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::INFINITY, 1.0)],
            [f64::INFINITY],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(-1.0, f64::NEG_INFINITY)],
            [f64::NEG_INFINITY],
            rtol = rtol
        );
        np_assert_allclose!(
            [crate::agm(f64::NEG_INFINITY, -1.0)],
            [f64::NEG_INFINITY],
            rtol = rtol
        );
    }
}
