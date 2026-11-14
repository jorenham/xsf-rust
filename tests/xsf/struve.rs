xsref::tests! {
    test_itstruve0_f64: "itstruve0", "d-d", |x| xsf::itstruve0(x[0]);
    test_it2struve0_f64: "it2struve0", "d-d", |x| xsf::it2struve0(x[0]);
}

#[test]
fn test_itstruve0_it2struve0_large_x() {
    // the Neumann series for 0.5 <= x < 40, and the asymptotic expansion beyond (mpmath)
    let actual = [
        xsf::itstruve0(29.746_969_855_043_805),
        xsf::it2struve0(24.450_718_139_337_205),
    ];
    let expected = [3.083_870_170_792_807, 0.026_268_893_388_791_609];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);

    let actual = [
        xsf::it2struve0(43.234_469_249_914_99),
        xsf::it2struve0(1000.0),
    ];
    let expected = [0.014_722_474_596_212_911, 0.000_661_413_224_079_251_6];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

xsref::tests! {
    test_itmodstruve0_f64: "itmodstruve0", "d-d", |x| xsf::itmodstruve0(x[0]);
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
        let rel = (xsf::itmodstruve0(x) - expected).abs() / expected.abs();
        assert!(rel < 4e-15, "itmodstruve0({x}): rel = {rel:e}");
    }
}

xsref::tests! {
    test_struve_h_f64: "struve_h", "d_d-d", |x| xsf::struve_h(x[0], x[1]);
    test_struve_l_f64: "struve_l", "d_d-d", |x| xsf::struve_l(x[0], x[1]);
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
        let rel = (xsf::struve_l(v, z) - expected).abs() / expected.abs();
        assert!(rel < 2e-12, "struve_l({v}, {z}): rel = {rel:e}");
    }
}

#[test]
fn test_struve_h_l_accurate() {
    // near the asymptotic expansion's threshold, a large leading-term exponent, one that
    // overflows (so the Bessel series wins), and L near v = -k - 1/2 (reference values from Arb)
    let actual = [
        xsf::struve_h(1.243_519_904_577_050_6, 24.304_034_546_767_898),
        xsf::struve_h(-4.625_056_973_249_563, 7.710_273_901_087_9e-85),
        xsf::struve_h(-26.903_775_197_924_794, 6.123_636_709_931_788e-9),
        xsf::struve_h(800.0, 320.0),
        xsf::struve_l(-11.499_999_924_029_11, 11.984_087_896_438_211),
        xsf::struve_l(-10.500_000_254_090_715, 12.134_779_880_436_286),
    ];
    let expected = [
        1.195_317_234_632_318_2,
        9.796_550_282_644_997e305,
        6.924_530_589_547_383e245,
        2.599_493_516_768_558e-215,
        88.475_249_644_313_31,
        250.252_585_832_560_87,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 4e-15);
}
