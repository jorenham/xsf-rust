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

#[test]
fn test_struve_h_bessel_series_no_underflow() {
    // large orders, where the coefficients (z/2)^n / n! of the Bessel series are large, but no
    // Bessel function underflows (reference values from Arb)
    let actual = [
        xsf::struve_h(780.37, 288.403_150_312_66),
        xsf::struve_h(824.714_790_883_217_1, 320.092_396_316_967_1),
    ];
    let expected = [6.968_982_334_069_738e-237, 1.200_536_822_257_174_5e-232];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 4e-15);
}

#[test]
fn test_struve_h_l_power_series_cancellation() {
    // power series whose largest term is 4e4 to 3e8 times the sum (reference values from Arb)
    let actual = [
        xsf::struve_h(-7.556_101_727_808_929_5, 15.297_459_023_655_978),
        xsf::struve_h(0.049_505_315_039_433_12, 13.208_434_290_082_854),
        xsf::struve_l(-34.020_820_230_172_97, 23.545_083_303_272_627),
    ];
    let expected = [
        4.824_120_592_343_075e-5,
        0.007_575_369_029_986_219_7,
        0.000_419_673_272_353_247_93,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_struve_h_l_large_order() {
    // leading terms that over- or underflow in f64 (reference values from Arb)
    let actual = [
        xsf::struve_h(992.061_461_574_324_9, 393.320_394_318_406),
        xsf::struve_l(-176.739_829_968_973_7, 135.785_696_147_920_33),
        xsf::struve_l(-421.776_397_824_793_66, 238.513_623_635_582_8),
    ];
    let expected = [
        4.661_060_885_142_991e-270,
        2_331_323_077_507.244,
        -4.755_481_199_035_166_4e45,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_struve_l_large_z_iv() {
    // the large-z expansion, where `iv` of large order was the largest error (Arb)
    let actual = [
        xsf::struve_l(-520.094_865_353_031_1, 361.026_378_627_839_9),
        xsf::struve_l(-441.512_943_384_891_46, 305.847_377_211_526_63),
        xsf::struve_l(-1_010.337_195_294_990_9, 701.198_096_234_033_9),
    ];
    let expected = [
        35_323_083_951.276_2,
        261_643_267.355_977_83,
        5.708_755_380_857_371e21,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}
