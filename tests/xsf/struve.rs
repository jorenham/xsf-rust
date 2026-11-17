xsref::tests! {
    test_itstruve0_f64: "itstruve0", "d-d", |x| xsf::itstruve0(x[0]);
    test_it2struve0_f64: "it2struve0", "d-d", |x| xsf::it2struve0(x[0]);
}

#[test]
fn test_itstruve0_it2struve0_large_x() {
    // the Neumann series for 0.5 <= x < 40, and the asymptotic expansion beyond (mpmath)
    let actual = [xsf::itstruve0(29.75), xsf::it2struve0(24.5)];
    let expected = [3.083_654_670_049_084_3, 0.026_538_829_242_429_812];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);

    let actual = [xsf::it2struve0(43.25), xsf::it2struve0(1000.0)];
    let expected = [0.014_760_727_000_063_827, 0.000_661_413_224_079_251_6];
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
        (21.25, 151_236_397.186_908_4),
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
        xsf::struve_h(1.25, 24.3125),
        xsf::struve_h(-4.625, 7.75e-85),
        xsf::struve_h(-26.875, 6e-9),
        xsf::struve_h(800.0, 320.0),
        xsf::struve_l(-11.499_999_924_029_11, 11.984_087_896_438_211),
        xsf::struve_l(-10.500_000_254_090_715, 12.125),
    ];
    let expected = [
        1.212_287_286_917_358_4,
        9.505_033_476_174_251e305,
        5.884_482_434_547_347e245,
        2.599_493_516_768_558e-215,
        88.475_249_644_313_31,
        247.092_250_919_801_02,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 4e-15);
}

#[test]
fn test_struve_h_bessel_series_no_underflow() {
    // large orders, where the coefficients (z/2)^n / n! of the Bessel series are large, but no
    // Bessel function underflows (reference values from Arb)
    let actual = [xsf::struve_h(780.375, 288.375), xsf::struve_h(824.75, 320.125)];
    let expected = [6.404_200_113_649_084e-237, 1.232_368_280_124_034_4e-232];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 4e-15);
}

#[test]
fn test_struve_h_l_power_series_cancellation() {
    // power series whose largest term is 4e4 to 2e8 times the sum (reference values from Arb)
    let actual =
        [xsf::struve_h(-7.555, 15.3), xsf::struve_h(0.06, 13.18), xsf::struve_l(-34.03, 23.55)];
    let expected =
        [6.617_629_020_047_206e-5, 0.000_323_590_515_247_530_25, -0.000_412_097_474_182_587_1];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_struve_h_l_large_order() {
    // leading terms that over- or underflow in f64 (reference values from Arb)
    let actual = [
        xsf::struve_h(992.0625, 393.3125),
        xsf::struve_l(-176.75, 135.8125),
        xsf::struve_l(-421.75, 238.5),
    ];
    let expected =
        [4.561_590_588_241_463e-270, 2_409_448_390_992.781_7, -4.365_343_490_813_207_5e45];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_struve_l_large_z_iv() {
    // the large-z expansion, where `iv` of large order was the largest error (Arb)
    let actual = [
        xsf::struve_l(-520.0625, 361.0),
        xsf::struve_l(-441.5625, 306.0),
        xsf::struve_l(-1010.375, 701.5),
    ];
    let expected = [33_718_042_048.289_43, 200_138_468.051_341_92, 3.838_949_295_084_25e22];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}
