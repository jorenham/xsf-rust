/// Based on `scipy.special.tests.test_TestInverseErrorFunction.test_literal_values`
#[test]
fn test_erfinv_literal_values() {
    let y = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
    let actual = y.map(xsf::erfinv);
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
    xsref::np_assert_allclose!(&actual, &expected, atol = 1e-15);
}

#[test]
fn test_erfcinv() {
    xsref::test("erfcinv", "d-d", |x| xsf::erfcinv(x[0]));
}

/// Based on the "cephes::erfcinv subnormal arguments" test in xsf main
#[test]
fn test_erfcinv_subnormal() {
    let y = [
        5e-324,
        1e-323,
        1.5e-323,
        f64::MIN_POSITIVE,
        2.0 * f64::MIN_POSITIVE,
    ];
    let actual = y.map(xsf::erfcinv);
    let expected = [
        27.213_293_210_812_95,
        27.200_563_366_536_258,
        27.193_114_126_203_97,
        26.543_258_454_250_98,
        26.530_207_562_276_377,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 2e-16);
}
