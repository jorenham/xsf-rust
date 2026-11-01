#[test]
fn test_pbwa() {
    xsref::test("pbwa", "d_d-d_d", |x| xsf::pbwa(x[0], x[1]));
}

#[test]
fn test_pbdv() {
    xsref::test("pbdv", "d_d-d_d", |x| xsf::pbdv(x[0], x[1]));
}

#[test]
fn test_pbdv_negative_x() {
    // non-integer v > 0 and -10 <= x < -5.8 (reference values from mpmath)
    for (v, x, d, dp) in [
        (
            7.821_919_754_766_831_5,
            -6.126_185_993_244_382,
            142.563_617_884_111_64,
            -29.698_107_853_800_717,
        ),
        (
            0.373_152_185_433_848_87,
            -5.814_209_871_190_306,
            -287.813_586_230_075_5,
            762.841_520_864_751_4,
        ),
        (9.5, -8.0, 9_081.000_904_297_447, -20_350.620_138_313_516),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let rel = ((y - d) / d).abs().max(((yp - dp) / dp).abs());
        assert!(rel < 1e-13, "v = {v}, x = {x}: rel = {rel:e}");
    }

    // near an integer v, or v >= x^2 / 3, where the series is not used
    for (v, x, d, dp, tol) in [
        (
            1e-20,
            -8.0,
            1.125_351_464_193_068_7e-7,
            4.501_408_120_768_454_4e-7,
            1e-6,
        ),
        (
            37.999_998_999_999,
            -9.0,
            -9.676_543_736_549_405e21,
            -7.936_223_618_237_496e21,
            1e-13,
        ),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let rel = ((y - d) / d).abs().max(((yp - dp) / dp).abs());
        assert!(rel < tol, "v = {v}, x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_pbdv_moderate_x() {
    // 1 < x < 12 (4.5 < x for v >= 0), and -30 <= v < 0 for -10 <= x < -5.8 (reference values
    // from mpmath)
    for (v, x, d, dp) in [
        (
            -2.5,
            5.9,
            1.754_156_195_473_642_3e-6,
            -5.856_502_185_207_843e-6,
        ),
        (
            -18.85,
            1.9,
            2.267_296_304_203_875_8e-12,
            -9.977_778_237_404_963e-12,
        ),
        (10.5, 7.0, 981.861_674_474_586_7, -1_483.395_230_884_999_4),
        (-0.01, -5.9, 26.989_634_433_251_59, -74.782_517_607_607_91),
        (-25.5, -7.0, 5_967.673_485_180_04, -36_282.507_255_654_57),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let rel = ((y - d) / d).abs().max(((yp - dp) / dp).abs());
        assert!(rel < 1e-13, "v = {v}, x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_pbdv_large_x() {
    // 36.5 <= x < 60 (53 <= x for v >= 0), near the underflow (reference values from mpmath)
    for (v, x, d, dp) in [
        (
            -19.5,
            38.0,
            2.262_826_187_826_714e-188,
            -4.413_886_212_120_611_4e-187,
        ),
        (
            -5.5,
            45.0,
            1.101_359_989_411_643_5e-229,
            -2.491_478_130_368_397e-228,
        ),
        (
            10.5,
            54.5,
            5.420_662_759_314_265e-305,
            -1.466_653_517_047_072_7e-303,
        ),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let rel = ((y - d) / d).abs().max(((yp - dp) / dp).abs());
        assert!(rel < 5e-13, "v = {v}, x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_pbdv_large_v() {
    // reference values from mpmath
    for (v, x, d, dp) in [
        (
            150.3,
            -30.0,
            -2.454_190_581_677_082_3e144,
            2.101_458_980_855_297_5e145,
        ),
        (
            400.3,
            80.0,
            1.361_794_629_135_228_6e61,
            -4.716_957_354_895_083e62,
        ),
        (
            13.1,
            -30.0,
            1.739_974_241_615_124_7e86,
            -2.526_760_712_361_696_6e87,
        ),
        (25.0, -12.0, -17_039_004_193.468_716, -57_453_241_630.516_9),
        (
            -200.3,
            0.5,
            5.771_946_734_121_321e-191,
            -8.160_152_171_862_542e-190,
        ),
        (
            -250.3,
            2.5,
            7.169_904_499_357_33e-264,
            -1.136_836_715_169_693_2e-262,
        ),
        (
            -500.7,
            -60.0,
            7.165_454_288_010_767e159,
            -2.680_875_133_938_43e161,
        ),
        // integer v
        (
            25.0,
            -45.0,
            -2.532_971_203_579_280_5e-179,
            -5.556_756_936_014_709e-178,
        ),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        for (y, r) in [(y, d), (yp, dp)] {
            let rel = ((y - r) / r).abs();
            assert!(rel < 5e-13, "v = {v}, x = {x}: {y:e} != {r:e}");
        }
    }

    // at and near x = 0, where D_v(0) = 0 for odd v, and D_v'(0) = 0 for even v
    assert_eq!(xsf::pbdv(300.0, 0.0).1, 0.0);
    assert_eq!(xsf::pbdv(301.0, 0.0).0, 0.0);
    for (v, x, i, r) in [
        (300.0, 0.0, 0, 3.753_274_111_571_926e306),
        (300.0, 5e-6, 1, -5.639_294_345_575_914e303),
        (301.0, 6e-6, 0, 6.778_413_033_236_749e303),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let y = [y, yp][i];
        assert!(
            (y / r - 1.0).abs() < 5e-13,
            "v = {v}, x = {x}: {y:e} != {r:e}"
        );
    }

    // overflow, e.g. D_v(-60) = -8.9e367 and D_v'(-60) = 2.6e369 for v = 30.3
    assert_eq!(xsf::pbdv(30.3, -60.0), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(xsf::pbdv(-30.3, -60.0), (f64::INFINITY, f64::NEG_INFINITY));
    assert_eq!(xsf::pbdv(1000.5, -1e-200), (f64::INFINITY, f64::INFINITY));
    assert_eq!(
        xsf::pbdv(350.5, 3.0),
        (f64::NEG_INFINITY, f64::NEG_INFINITY)
    );
}

#[test]
fn test_pbvv() {
    xsref::test("pbvv", "d_d-d_d", |x| xsf::pbvv(x[0], x[1]));
}

#[test]
fn test_pbdv_pbvv_v_rounding() {
    // v +- 1 rounding to an integer (reference values from mpmath)
    let dp = xsf::pbdv(-0.999_999_999_999_999_9, -3.0).1;
    let wp = xsf::pbvv(1.999_999_999_999_999_8, -1.0).1;
    for (y, r) in [
        (dp, -35.730_584_217_091_32),
        (wp, -0.025_839_571_037_084_53),
    ] {
        assert!((y / r - 1.0).abs() < 1e-13, "{y:e} != {r:e}");
    }
}

#[test]
fn test_pbvv_extreme() {
    // v >= 0 and 7.5 < |x| <= 10, and v near an integer (reference values from mpmath)
    for (v, x, w, wp) in [
        (
            8.463_751_277_865_178,
            7.641_451_046_470_024,
            0.023_227_491_821_967_712,
            0.050_155_659_621_339_22,
        ),
        (
            9.954_321_571_892_986,
            -7.634_719_796_831_582,
            -0.001_729_313_216_016_740_6,
            0.002_985_548_082_122_786,
        ),
        (6.0, -8.0, -5.663_320_431_230_510_5, 16.810_707_783_832_23),
        // near an integer v
        (
            2.000_000_01,
            7.025_168_798_105_803_5,
            601.282_319_691_077_8,
            1_829.270_702_453_187_4,
        ),
    ] {
        let (y, yp) = xsf::pbvv(v, x);
        let rel = ((y - w) / w).abs().max(((yp - wp) / wp).abs());
        assert!(rel < 1e-12, "v = {v}, x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_pbvv_negative_x() {
    // -7.5 <= x < 0 and v + 1 <= x^2 / 2, and x < -10 (reference values from mpmath)
    for (v, x, w, wp) in [
        (
            9.48,
            -7.46,
            0.000_247_241_678_108_458_4,
            -0.000_155_355_245_089_156_9,
        ),
        (
            6.5,
            -7.0,
            -0.000_528_639_039_615_62,
            -0.001_289_473_434_129_782_2,
        ),
        (
            4.25,
            -5.5,
            -0.319_862_046_579_462_3,
            0.397_396_669_652_463_9,
        ),
        (
            3.0,
            -20.0,
            1.374_924_131_582_966_8e38,
            -1.347_072_098_819_639_2e39,
        ),
        (
            9.25,
            -33.0,
            2.796_250_555_935_439e102,
            -4.526_043_015_443_377e103,
        ),
        (
            2.5,
            -15.0,
            -9.681_709_256_248_095e-23,
            -7.098_834_733_194_973e-22,
        ),
    ] {
        let (y, yp) = xsf::pbvv(v, x);
        let (rel, relp) = (((y - w) / w).abs(), ((yp - wp) / wp).abs());
        assert!(
            rel < 1e-12 && relp < 1e-12,
            "v = {v}, x = {x}: {rel:e}, {relp:e}"
        );
    }

    // half-integer v, where V_v(-x) overflows (subnormal results)
    let (y, yp) = xsf::pbvv(2.5, -54.0);
    let (rel, relp) = (
        (y / -1.616_037_873_6e-313 - 1.0).abs(),
        (yp / -4.355_816_750_767e-312 - 1.0).abs(),
    );
    assert!(rel < 1e-9 && relp < 1e-9, "{rel:e}, {relp:e}");
}
