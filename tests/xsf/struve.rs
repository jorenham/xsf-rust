#[test]
fn test_itstruve0_f64() {
    xsref::test("itstruve0", "d-d", |x| xsf::itstruve0(x[0]));
}

#[test]
fn test_it2struve0_f64() {
    xsref::test("it2struve0", "d-d", |x| xsf::it2struve0(x[0]));
}

#[test]
fn test_itstruve0_neumann() {
    // Neumann series for 0.5 <= x < 40 (reference values from mpmath)
    for (x, expected) in [
        (5.5, 1.937_558_487_381_615_6),
        (29.746_969_855_043_805, 3.083_870_170_792_807),
        (37.5, 3.046_057_878_828_238_6),
    ] {
        let rel = (xsf::itstruve0(x) - expected).abs() / expected;
        assert!(rel < 1e-15, "itstruve0({x}): rel = {rel:e}");
    }
}

#[test]
fn test_it2struve0_neumann_asymptotic() {
    // Neumann series for 0.5 <= x < 40, asymptotic expansion beyond (reference values from mpmath)
    for (x, expected, tol) in [
        (3.5, 0.103_660_422_008_807_57, 2e-15),
        (24.450_718_139_337_205, 0.026_268_893_388_791_609, 1e-13),
        (43.234_469_249_914_99, 0.014_722_474_596_212_911, 1e-15),
        (1000.0, 0.000_661_413_224_079_251_6, 1e-15),
    ] {
        let rel = (xsf::it2struve0(x) - expected).abs() / expected;
        assert!(rel < tol, "it2struve0({x}): rel = {rel:e}");
    }
}

#[test]
fn test_itmodstruve0_f64() {
    xsref::test("itmodstruve0", "d-d", |x| xsf::itmodstruve0(x[0]));
}

#[test]
fn test_itmodstruve0_series_asymptotic() {
    // series for x <= 40, asymptotic expansion beyond (reference values from mpmath)
    for (x, expected) in [
        (12.5, 32_045.561_159_494_297),
        (21.281_112_146_816_426, 155_893_845.553_037_7),
        (37.5, 1.280_693_999_522_079e15),
        (45.0, 2.107_522_588_189_009_7e18),
        (110.0, 2.265_187_556_512_247_2e46),
    ] {
        let rel = (xsf::itmodstruve0(x) - expected).abs() / expected;
        assert!(rel < 4e-15, "itmodstruve0({x}): rel = {rel:e}");
    }
}

#[test]
fn test_struve_h_f64() {
    xsref::test("struve_h", "d_d-d", |x| xsf::struve_h(x[0], x[1]));
}

#[test]
fn test_struve_l_f64() {
    xsref::test("struve_l", "d_d-d", |x| xsf::struve_l(x[0], x[1]));
}

#[test]
fn test_struve_l_near_half_integer() {
    // asymptotic expansion, which is exact at half-integers for L_v - I_{-v} (reference
    // values from mpmath)
    for (v, z, expected) in [
        (0.5, 13.0, 48_951.351_990_355_41),
        (-7.499_999_999_7, 15.0, 50_900.655_128_588_08),
        (-15.500_000_1, 15.0, 162.742_950_300_856_3),
    ] {
        let rel = (xsf::struve_l(v, z) - expected).abs() / expected;
        assert!(rel < 2e-12, "struve_l({v}, {z}): rel = {rel:e}");
    }
}
