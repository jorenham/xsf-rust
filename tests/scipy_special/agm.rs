//! Translated from `scipy.special.tests.test_basic.test_agm_simple`

use core::f64::consts::SQRT_2;
use xsref::np_assert_allclose;

#[test]
#[allow(clippy::too_many_lines)]
fn test_agm_simple() {
    let rtol = 1e-13;

    // Gauss's constant
    np_assert_allclose!(
        [1.0 / xsf::agm(1.0, SQRT_2)],
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
        .map(|ab| xsf::agm(ab.0, ab.1)),
        [1.0, agm13, agm15, agm13, 3.0, agm35],
        rtol = rtol
    );

    // Computed by the iteration formula using mpmath, with mpmath.mp.prec = 1000:
    let agm12 = 1.456_791_031_046_906_8;
    np_assert_allclose!([xsf::agm(1.0, 2.0)], [agm12], rtol = rtol);
    np_assert_allclose!([xsf::agm(2.0, 1.0)], [agm12], rtol = rtol);
    np_assert_allclose!([xsf::agm(-1.0, -2.0)], [-agm12], rtol = rtol);
    np_assert_allclose!([xsf::agm(24.0, 6.0)], [13.458_171_481_725_614], rtol = rtol);
    np_assert_allclose!(
        [xsf::agm(13.0, 123_456_789.5)],
        [11_111_458.498_599_306],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(1e30, 1.0)],
        [2.229_223_055_945_383e+28],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(1e-22, 1.0)],
        [0.030_182_566_420_169_886],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(1e150, 1e180)],
        [2.229_223_055_945_383e+178],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(1e180, 1e-150)],
        [2.063_472_251_016_267_7e+177],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(1e-150, 1e-170)],
        [3.311_261_967_046_375_6e-152],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(f64::MIN_POSITIVE, f64::MAX)],
        [1.989_207_205_001_547_3e+305],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(0.75 * f64::MAX, f64::MAX)],
        [1.564_904_312_298_045e+308],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(f64::MIN_POSITIVE, 3.0 * f64::MIN_POSITIVE)],
        [4.146_684_986_673_500_5e-308],
        rtol = rtol
    );

    // zero, nan and inf cases.
    np_assert_allclose!([xsf::agm(0.0, 0.0)], [0.0], rtol = rtol);
    np_assert_allclose!([xsf::agm(99.0, 0.0)], [0.0], rtol = rtol);

    np_assert_allclose!([xsf::agm(99.0, 0.0)], [0.0], rtol = rtol);
    np_assert_allclose!([xsf::agm(0.0, f64::INFINITY)], [f64::NAN], rtol = rtol);
    np_assert_allclose!([xsf::agm(f64::INFINITY, 0.0)], [f64::NAN], rtol = rtol);
    np_assert_allclose!([xsf::agm(0.0, f64::NEG_INFINITY)], [f64::NAN], rtol = rtol);
    np_assert_allclose!([xsf::agm(f64::NEG_INFINITY, 0.0)], [f64::NAN], rtol = rtol);
    np_assert_allclose!(
        [xsf::agm(f64::INFINITY, f64::NEG_INFINITY)],
        [f64::NAN],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(f64::NEG_INFINITY, f64::INFINITY)],
        [f64::NAN],
        rtol = rtol
    );
    np_assert_allclose!([xsf::agm(1.0, f64::NAN)], [f64::NAN], rtol = rtol);
    np_assert_allclose!([xsf::agm(f64::NAN, -1.0)], [f64::NAN], rtol = rtol);

    np_assert_allclose!([xsf::agm(1.0, f64::INFINITY)], [f64::INFINITY], rtol = rtol);
    np_assert_allclose!([xsf::agm(f64::INFINITY, 1.0)], [f64::INFINITY], rtol = rtol);
    np_assert_allclose!(
        [xsf::agm(-1.0, f64::NEG_INFINITY)],
        [f64::NEG_INFINITY],
        rtol = rtol
    );
    np_assert_allclose!(
        [xsf::agm(f64::NEG_INFINITY, -1.0)],
        [f64::NEG_INFINITY],
        rtol = rtol
    );
}
