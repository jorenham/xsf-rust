//! Translated from `scipy.special.tests.test_basic.test_agm_simple`

use core::f64::consts::SQRT_2;
use xsref::{np_assert_allclose, np_assert_equal};

#[test]
fn test_agm_simple() {
    let rtol = 1e-13;

    // Gauss's constant
    np_assert_allclose!(
        [1.0 / xsf::agm(1.0, SQRT_2)],
        [0.834_626_841_674_073_2],
        rtol = rtol
    );

    // Wolfram Alpha, ArithmeticGeometricMean[a, b]
    let agm13 = 1.863_616_783_244_897;
    let agm15 = 2.604_008_190_530_94;
    let agm35 = 3.936_235_503_649_555;
    // mpmath iteration with mp.prec = 1000
    let agm12 = 1.456_791_031_046_906_8;
    let tiny = f64::MIN_POSITIVE;
    let huge = f64::MAX;
    for (a, b, expected) in [
        (1.0, 1.0, 1.0),
        (1.0, 3.0, agm13),
        (1.0, 5.0, agm15),
        (3.0, 1.0, agm13),
        (3.0, 3.0, 3.0),
        (3.0, 5.0, agm35),
        (1.0, 2.0, agm12),
        (2.0, 1.0, agm12),
        (-1.0, -2.0, -agm12),
        (24.0, 6.0, 13.458_171_481_725_614),
        (13.0, 123_456_789.5, 11_111_458.498_599_306),
        (1e30, 1.0, 2.229_223_055_945_383e+28),
        (1e-22, 1.0, 0.030_182_566_420_169_886),
        (1e150, 1e180, 2.229_223_055_945_383e+178),
        (1e180, 1e-150, 2.063_472_251_016_267_7e+177),
        (1e-150, 1e-170, 3.311_261_967_046_375_6e-152),
        (tiny, huge, 1.989_207_205_001_547_3e+305),
        (0.75 * huge, huge, 1.564_904_312_298_045e+308),
        (tiny, 3.0 * tiny, 4.146_684_986_673_500_5e-308),
    ] {
        np_assert_allclose!([xsf::agm(a, b)], [expected], rtol = rtol);
    }

    // zero, nan and inf cases
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    for (a, b, expected) in [
        (0.0, 0.0, 0.0),
        (99.0, 0.0, 0.0),
        (-1.0, 10.0, nan),
        (0.0, inf, nan),
        (inf, 0.0, nan),
        (0.0, -inf, nan),
        (-inf, 0.0, nan),
        (inf, -inf, nan),
        (-inf, inf, nan),
        (1.0, nan, nan),
        (nan, -1.0, nan),
        (1.0, inf, inf),
        (inf, 1.0, inf),
        (-1.0, -inf, -inf),
        (-inf, -1.0, -inf),
    ] {
        np_assert_equal!([xsf::agm(a, b)], [expected]);
    }
}
