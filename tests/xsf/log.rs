use num_complex::c64;

xsref::tests! {
    test_log1p_f64: "log1p", "d-d", |x| xsf::log1p(x[0]);
    test_log1p_c64: "log1p", "cd-cd", |x| xsf::log1p(c64(x[0], x[1]));
    test_log1pmx_f64: "log1pmx", "d-d", |x| xsf::log1pmx(x[0]);
    test_xlogy_f64: "xlogy", "d_d-d", |x| xsf::xlogy(x[0], x[1]);
    test_xlogy_c64: "xlogy", "cd_cd-cd", |x| xsf::xlogy(c64(x[0], x[1]), c64(x[2], x[3]));
    test_xlog1py_f64: "xlog1py", "d_d-d", |x| xsf::xlog1py(x[0], x[1]);
    test_xlog1py_c64: "xlog1py", "cd_cd-cd", |x| xsf::xlog1py(c64(x[0], x[1]), c64(x[2], x[3]));
}

#[test]
fn test_log1p_c64_double_double() {
    // |1 + z|^2 = 1 + 2^-50 exactly
    let w = xsf::log1p(c64(-2f64.powi(-25), 2f64.powi(-12)));
    xsref::np_assert_allclose!([w.re], [0.5 * 2f64.powi(-50).ln_1p()], rtol = 1e-15);
}

#[test]
fn test_xlogy_xlog1py_zero() {
    assert_eq!(xsf::xlogy(0.0, f64::INFINITY), 0.0);
    assert_eq!(xsf::xlog1py(0.0, -1.0), 0.0);
    assert_eq!(xsf::xlogy(c64(0.0, 0.0), c64(0.0, 0.0)), c64(0.0, 0.0));
    assert!(xsf::xlogy(0.0, f64::NAN).is_nan());
}

#[test]
fn test_log1p_f64_special() {
    for x in [0.0, -0.0, 5e-324, -5e-324, 1e-310, -f64::MIN_POSITIVE] {
        assert_eq!(xsf::log1p(x).to_bits(), x.to_bits());
    }
    assert_eq!(xsf::log1p(-1.0), f64::NEG_INFINITY);
    assert_eq!(xsf::log1p(f64::INFINITY), f64::INFINITY);
    for x in [-1.000_000_000_000_000_2, -1.5, f64::NEG_INFINITY, f64::NAN] {
        assert!(xsf::log1p(x).is_nan());
    }
}

#[test]
fn test_log1p_f64_hard() {
    // MPFR, rounded to nearest
    for (x, expected) in [
        // near 1 / sqrt(2) - 1 and sqrt(2) - 1
        (-0.292_893_470_340_391_5, -0.346_573_945_992_844_36),
        (-0.292_968_749_962_865_7, -0.346_680_413_161_215_3),
        (0.398_826_611_133_659_14, 0.335_633_750_285_246_3),
        (0.414_213_280_687_302_3, 0.346_573_391_098_018_6),
        // next to -2^-8 and 2^-7
        (-0.003_906_250_000_000_001, -0.003_913_899_321_136_329_5),
        (0.007_812_500_000_000_002, 0.007_782_140_442_054_951),
        // tiny
        (5.551_115_123_125_783e-17, 5.551_115_123_125_783e-17),
        (1e-10, 9.999_999_999_500_001e-11),
        // 1 + x is not exact
        (0.576_440_225_422_903_6, 0.455_169_283_282_184_6),
        (1.711_198_943_589_982_6, 0.997_390_951_663_971_3),
        // large |log1p(x)|
        (-0.999_999_999_999_999_9, -36.736_800_569_677_1),
        (-0.626_303_888_064_771_1, -0.984_312_346_785_580_9),
        (-0.377_447_143_737_160_1, -0.473_926_744_588_714_87),
        (9_007_199_254_740_992.0, 36.736_800_569_677_1),
        (f64::MAX, 709.782_712_893_384),
    ] {
        assert_eq!(xsf::log1p(x).to_bits(), f64::to_bits(expected), "x = {x:e}");
    }
}
