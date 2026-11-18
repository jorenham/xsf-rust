use num_complex::c64;

xsref::tests! {
    // the table is called "hyperu" instead of "hypu"
    test_hypu: "hyperu", "d_d_d-d", |x| xsf::hypu(x[0], x[1], x[2]);
}

#[test]
fn test_hypu_tiny_a() {
    // mpmath (scipy/xsf#350)
    let abx = [(1e-17, 0.5, 100.0), (-1e-17, 0.5, 100.0), (-1.0 + f64::EPSILON / 2.0, 0.5, 200.0)];
    let expected = [1.0, 1.0, 199.499_999_999_999_9];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-16);

    // Arb (scipy/xsf#436)
    let abx = [(1e-20, 0.5, 50.0), (1e-15, 0.5, 50.0), (1e-20, 0.5, 40.0), (-1e-20, -0.5, 30.0)];
    let expected = [1.0, 0.999_999_999_999_996_1, 1.0, 1.0];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb (scipy/xsf#459)
    let abx = [(1e-20, 10.0, 0.01), (-1e-20, 15.0, 0.1), (1e-30, 30.0, 1.0), (1e-40, 20.0, 0.002)];
    let expected = [
        408.762_021_103_157_6,
        -6_934.305_358_080_271,
        1.859_511_646_545_143_5,
        1.223_736_600_747_381e27,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 1.2014, 0.99999999996, and 1.00028382 in xsf
    let abx = [(7.5e-16, 6.0, 0.0024), (-2.9e-14, 28.0, 4.55), (2.387e-14, 6.0, 0.01834)];
    let expected = [1.226_735_389_035_132, 0.997_697_971_141_589_2, 1.000_282_508_087_199];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb, near b = 2; 1.1332942, 160.62366, 24.321023, and 0.71354002 in xsf
    let abx = [
        (2e-14, 2.0 - 2e-12, 1.5e-13),
        (8e-14, 2.0 - 3e-14, 5e-16),
        (1e-13, 2.0 + 1e-15, 5e-15),
        (-3e-13, 2.0 - 1e-14, 1e-12),
    ];
    let expected = [
        1.133_333_333_326_205_4,
        160.999_999_999_844,
        21.000_000_000_005_023,
        0.699_999_999_991_843_7,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb, near integer b > 2; 100.97221, -8.7183826e10, 3.1884856e30, -4.2000327e45,
    // 2.4406707e33, and -3.6770736e16 in xsf
    let abx = [
        (1e-10, 3.0 - 1e-12, 1e-6),
        (-5e-16, 3.0 - 4e-15, 7.5e-14),
        (2e-13, 4.0 - 1e-13, 5e-15),
        (-7e-16, 5.0 + 1e-11, 1e-15),
        (0.3, 40.0 + 1e-10, 2.0),
        (-2.3e-7, 3.0 - 4e-13, 2.5e-12),
    ];
    let expected = [
        101.000_200_005_729_76,
        -88_888_888_887.891_33,
        3.199_999_999_989_569_5e30,
        -4.200_000_001_503_388_5e45,
        2.440_613_605_416_168_2e33,
        -3.679_999_511_422_975e16,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; U(a + 1, b, x) overflows, or a step of its recurrence in b (NaN and inf in xsf)
    let actual = [xsf::hypu(0.05, 15.5, 2.5e-21), xsf::hypu(0.001, 3.0, 1e-155)];
    let expected = [6.367_508_711_417_84e307, 1.000_576_559_744_994e307];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb: 1.2125e-13, next to a zero, where DLMF 13.3.9 and the recurrence in a both cancel
    let actual = xsf::hypu(-0.49, 1.999, 0.821_185_784_361_567_1);
    assert!(actual.is_nan() || (actual / 1.212_531_236_482_603_3e-13 - 1.0).abs() < 1e-6);
}

#[test]
fn test_hypu_tiny_a_x_below_b() {
    // Arb, for x < b, where the asymptotic series is 1; 1220.1739 in xsf
    let abx = [
        (1e-15, 35.0 + 1e-14, 4.0),
        (-1e-16, 54.0 - 1e-3, 9.3),
        (1e-12, 51.0 - 1e-5, 19.0),
        (4e-15, 28.0 + 1.4e-4, 5.9),
    ];
    let expected = [
        1_829.955_866_859_770_1,
        -50_259.716_622_619_83,
        1.000_020_655_722_755_4,
        1.000_001_176_181_804,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hypu_tiny_a_x_near_b() {
    // Arb, where the asymptotic series has 9 digits; 0.999995306 and 1.00000089 in xsf
    let abx = [(1e-6, 690.0, 795.0), (-2e-7, 878.0, 957.0)];
    let expected = [0.999_995_308_569_181_3, 1.000_000_887_955_178];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hypu_x0() {
    // DLMF 13.2.7; +inf in xsf
    let actual = [xsf::hypu(0.0, 3.42, 0.0), xsf::hypu(-2.0, 5.0, 0.0)];
    xsref::np_assert_allclose!(actual, [1.0, 30.0], rtol = 1e-15);
    // (-75 - 1.4e-14)_204 = 1.36e311
    assert_eq!(xsf::hypu(-204.0, -75.000_000_000_000_01, 0.0), f64::INFINITY);
    let actual = [-3.0, -2.0].map(|a| xsf::hypu(a, f64::INFINITY, -0.0));
    assert_eq!(actual, [f64::NEG_INFINITY, f64::INFINITY]);
    // exact; finite, but the product (b + m - 1) ... b over- or underflows (0 for 1e-300 in xsf)
    let abx = [
        (-172.0, 0.015_625, 0.0),
        (-200.0, 1e-300, 0.0),
        (-170.0, 1e-310, 0.0),
        (-172.0, 1e-310, 0.0),
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    let expected = [
        2.120_016_193_072_081_3e307,
        3.943_289_336_823_953e72,
        4.269_068_009_004_692e-6,
        0.124_101_807_021_766_4,
    ];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // DLMF 13.2.16-19, with the sign of 1 / Gamma(a); +inf in xsf
    assert_eq!(xsf::hypu(-0.5, 2.0, 0.0), f64::NEG_INFINITY);
    assert_eq!(xsf::hypu(-10.5, 1.0, 0.0), f64::NEG_INFINITY);
    assert_eq!(xsf::hypu(-1.5, 2.0, 0.0), f64::INFINITY);

    // mpmath, Gamma(1 - b) / Gamma(a - b + 1) (DLMF 13.2.22); 0 in xsf
    let actual = [
        xsf::hypu(-5.421_052_631_578_947, -3.421_052_631_578_947_3, 0.0),
        xsf::hypu(-1.5, -0.500_000_000_000_000_9, 0.0),
    ];
    let expected = [-4.632_697_274_741_578e-15, 7.871_276_301_443_317e-16];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hypu_huge_a() {
    // NaN without ~2^31 iterations, like in xsf
    assert!(xsf::hypu(-2_147_483_649.0, 1.0, 1.0).is_nan());
    assert!(xsf::hypu(1.0, f64::INFINITY, 1.0).is_nan());
    // 1 in xsf
    assert!(xsf::hypu(-2_147_483_648.0, 1.0, 1.0).is_nan());
}

#[test]
fn test_hypu_huge_b() {
    // Arb; took 5 s and 90 s, or was wrong (2 and -1.3e9), like in xsf (and 16 s for the inf)
    let abx = [
        (-2.0, 1e9, 2.0),
        (1.0, 2_147_483_648.0, 2_147_483_648.0),
        (-1.0, 2_147_483_648.0, 2.0),
        (0.0, 1e10, 2.0),
    ];
    let expected = [9.999_999_97e17, 2.704_487_856_590_993_5e-5, -2_147_483_646.0, 1.0];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-12);
    assert_eq!(xsf::hypu(1.5, 2_147_483_647.5, 2.0), f64::INFINITY);
    assert_eq!(xsf::hypu(1.0, 1001.0, 2.0), f64::INFINITY);
    // no method is accurate here (Arb: 1e-6), like in xsf
    assert!(xsf::hypu(1.0, -1e6, 0.5).is_nan());
}

#[test]
fn test_hypu_neg_a() {
    // mpmath; NaN in xsf
    let abx = [(-9.6, 8.22, 29.4), (-0.5, 0.0, 2.0), (-95.93, 0.2502, 322.25), (-150.5, 0.3, 50.0)];
    let expected = [
        26_083_742_718.484_07,
        1.568_794_902_417_112,
        -6.372_049_560_642_321e217,
        2.060_433_058_188_866_4e273,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-11);
}

#[test]
fn test_hypu_neg_a_large_b() {
    // mpmath; NaN in xsf, and with the sign wrong if `chgul` ignores that its sum is small
    let actual = [xsf::hypu(-8.5, 99.1, 81.0), xsf::hypu(-8.2, 94.2, 74.6)];
    let expected = [-9_737_194_064.088_581, -4_550_121_582.442_501];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-7);

    // Arb: -2.6e7458, instead of a wrong 1 from `chgul`
    let actual = xsf::hypu(-4.2e-66, 26.000_000_000_000_018, 1e-300);
    assert_eq!(actual, f64::NEG_INFINITY);
}

#[test]
fn test_hypu_large_neg_a_large_b() {
    // mpmath; NaN in xsf
    let abx = [(-70.1, 45.4, 11.4), (-60.2, 30.6, 10.0), (-150.3, 60.1, 30.0), (-2.4, 98.3, 29.8)];
    let expected = [
        -8.440_298_308_670_45e121,
        -1.875_196_314_005_285_6e96,
        4.672_086_732_836_285e291,
        -3.999_673_741_113_012e20,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-8);
}

#[test]
fn test_hypu_near_neg_int_a_large_b() {
    // NaN, as the recurrence in b is unstable past b = x (Arb: 4.67e33, 1.64e146, 4.06e71,
    // -1.63e93, and 1.15e56), instead of -6.94e34, -5.00e146, -8.76e71, 3.06e92, and 1.39e56 (or
    // -1.81e9 from `chgul`); NaN, or 15% to 17% off, in xsf
    let abx = [
        (-12.0 + 4e-15, 814.0 + 4e-7, 503.0),
        (-35.999_999_999_999_986, 839.121_528_298_220_1, 341.787_507_068_174_76),
        (-43.999_999_999_999_97, 19.879_646_794_531_492, 0.164_247_209_584_363_92),
        (-6.000_000_000_000_011_5, 682.626_481_969_116_4, 257.665_363_374_949_04),
        (-5.000_000_000_000_002, 72.000_168_440_354_41, 2.888_910_315_901_288_8),
    ];
    for (a, b, x) in abx {
        assert!(xsf::hypu(a, b, x).is_nan());
    }
}

#[test]
fn test_hypu_near_neg_int_a_huge_b() {
    // Arb; by `chgul`, as there is no recurrence in b for b >= 1000
    let actual = [xsf::hypu(-5.0 + 1e-13, 1200.1, 1050.0)];
    xsref::np_assert_allclose!(actual, [-46_724_726_946.010_12], rtol = 1e-11);
}

#[test]
fn test_hypu_near_neg_int_a_small_x() {
    // Arb; 9.8e-9 and 6.9e-12 off by the recurrence in a, whose step to a near 0 cancels
    // (-2159.9538 for the first in xsf)
    let abx =
        [(-2.000_000_1, 4.999_999_999_9, 0.004_85), (-2.000_000_1, 2.999_999_999_9, 0.004_85)];
    let expected = [-2_159.997_821_507_856, 11.952_552_198_173_894];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 3.6e-12 off by the recurrence in a, which underestimated its rounding errors (NaN in
    // xsf)
    let actual = xsf::hypu(-169.13, 0.986_37, 0.0116);
    xsref::np_assert_allclose!([actual], [7.407_921_790_935_49e302], rtol = 1e-13);
}

#[test]
fn test_hypu_zero_start_b() {
    // Arb: 8.08e52; NaN, not 0 by the recurrence in b, whose starting values underflow to 0
    assert!(xsf::hypu(180.0, 180.5, 0.5).is_nan());
}

#[test]
fn test_hypu_poly() {
    // mpmath; -6.90e95, 7.04e186, and NaN in xsf
    let abx = [(-60.0, 30.0, 10.0), (-100.0, -5.01, 46.6), (-95.0, 15.0, 0.004)];
    let expected =
        [-6.398_262_182_829_726e95, -8.235_589_936_550_95e165, -1.614_746_581_537_901e165];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // mpmath; 3.27e96 and 5.35e190 in xsf
    let actual = [xsf::hypu(-60.5, 30.5, 10.0), xsf::hypu(-100.5, 40.5, 20.0)];
    let expected = [8.123_216_159_896_29e95, -4.104_645_857_033_403e177];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-8);
}

#[test]
fn test_hypu_poly_overflow() {
    // Arb: 1.46e336, -1.38e338, and -2.24e367; -inf, inf, and NaN in xsf
    assert_eq!(xsf::hypu(-181.0, 1.0, 30.2), f64::INFINITY);
    assert_eq!(xsf::hypu(-182.0, 1.0, 30.0), f64::NEG_INFINITY);
    assert_eq!(xsf::hypu(-197.0, 1.0, 0.5), f64::NEG_INFINITY);

    // Arb; NaN and inf in xsf
    let actual = [xsf::hypu(-171.0, 1.0, 0.5), xsf::hypu(-168.0, 2.0, 30.0)];
    let expected = [-1.277_404_819_038_795_2e308, 9.572_203_868_570_836e307];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb: 2.88e468; -inf in xsf
    assert_eq!(xsf::hypu(-237.0, -107.25, 101.0), f64::INFINITY);

    // Arb; 0, -0, inf, and inf in xsf
    let abx = [
        (-176.75, -143.75, 0.003),
        (-277.625, -199.625, 0.0329),
        (-157.8125, 33.1875, 8.75),
        (5.125, 330.125, 14.9),
    ];
    let expected = [
        2.692_507_169_282_326e-295,
        -3.854_488_043_803_876e-115,
        -3.436_886_227_248_617_6e301,
        3.387_267_647_106_950_4e303,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb: 3.80e-83005 and 3.59e201603; NaN in xsf
    assert_eq!(xsf::hypu(100_000.0, 200_002.0, 1000.0), 0.0);
    assert_eq!(xsf::hypu(150_000.0, 300_002.0, 100.0), f64::INFINITY);
}

#[test]
fn test_hypu_small_a() {
    // mpmath; 0.045763206, 0.067619064, and 7.5410439102e-4 in xsf
    let abx = [
        (1.1, 0.5, 15.0),
        (1.088_244_511_969_767_5, -2.443_990_400_454_129_7, 7.704_068_190_231_311),
        (2.5, 0.5, 15.0),
    ];
    let expected =
        [0.045_763_173_934_103_04, 0.067_618_993_534_387_55, 0.000_754_104_391_072_057_4];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb; 0.0831444376 in xsf
    let actual = xsf::hypu(3.5, 0.0, 0.002);
    xsref::np_assert_allclose!([actual], [0.083_144_400_245_006_18], rtol = 1e-13);
}

#[test]
fn test_hypu_digits() {
    // Arb; 0.1098838166, 0.2518536383, and 0.9378950474 in xsf
    let abx = [(0.8, 0.95, 15.0), (0.5, 0.7, 15.0), (0.03, -0.001, 7.5)];
    let expected = [0.109_883_867_556_135_05, 0.251_853_637_950_820_1, 0.937_895_047_309_956_8];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb; -48.29100258 in xsf
    let actual = xsf::hypu(-8.36, 0.9956, 7.86);
    xsref::np_assert_allclose!([actual], [-48.291_009_037_510_32], rtol = 1e-10);

    // Arb; still finite if `chgubi` fails (NaN with 15 digits)
    let actual = xsf::hypu(0.3, 172.0, 100.0);
    xsref::np_assert_allclose!([actual], [122_292_314.484_337_23], rtol = 1e-13);

    // Arb; `chgul` has all digits here, but claimed only 10
    let actual = [xsf::hypu(-36.0, -32.2, 16.2), xsf::hypu(-39.0, 23.5, 10.0)];
    let expected = [-1.512_832_103_737_844e44, 1.358_838_881_564_021_4e54];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hypu_chgus() {
    // mpmath; 3.07e67 and 6.13e-197 (poles of Gamma(a) and Gamma(1 + a - b)) in xsf
    let actual = [xsf::hypu(-28.0, 97.6, 0.0055), xsf::hypu(-74.2, -71.2, 0.0017)];
    let expected = [1.766_091_302_349_042e57, 5.933_445_224_488_155e-197];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // mpmath; 501.62007 and 902602731 (near the poles of Gamma(b)) in xsf
    let actual = [xsf::hypu(-2.17, -1.9977, 16.7), xsf::hypu(-7.5, -6.999_996, 15.0)];
    let expected = [501.620_214_189_965_45, 902_731_434.181_382_9];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-10);

    // Arb; -8.59e8 and 1.1254 (near the poles of Gamma(b)) in xsf
    let actual = [xsf::hypu(3.0, 1e-12, 1.5), xsf::hypu(-1.9, -0.999_999_999_6, 0.9)];
    let expected = [0.012_613_223_506_698_312, 1.081_954_303_595_999_8];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb; the terms of the two series in `chgus` cancel here, so that their sum barely changes
    // long before it converges
    let actual = xsf::hypu(10.0, 1e-17, 0.5);
    xsref::np_assert_allclose!([actual], [1.141_388_465_789_265_9e-8], rtol = 1e-13);

    // mpmath: -4.41e397 and 3.70e388; NaN in xsf
    assert_eq!(xsf::hypu(-14.2, 95.7, 0.0028), f64::NEG_INFINITY);
    assert_eq!(xsf::hypu(0.41, 93.1, 0.002), f64::INFINITY);

    // mpmath; inf in xsf
    let actual = xsf::hypu(100.2, 95.7, 0.000_426_150_094_022_470_2);
    xsref::np_assert_allclose!([actual], [1.797_665_191_584_464_4e308], rtol = 1e-13);

    // Arb; here r2 is subnormal, but negligible (the terms of r1 overflow if scaled by 1/|r2|)
    let abx = [(1.0, -0.01, 1e-307), (1.2, -0.01, 1e-307), (2.0, -0.01, 1e-307)];
    let expected = [0.990_099_009_900_990_1, 0.897_529_488_896_726_7, 0.492_586_572_090_044_8];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);
}

#[test]
fn test_hypu_chgus_sinpi() {
    // Arb; 2.17835e279, 5.0010e199, and 7.52252e149 (sin(pi b) near an integer b) in xsf
    let abx = [(0.2, 2.0 - 1e-12, 1e-280), (3.0, 2.0 + 1e-12, 1e-200), (2.5, 2.0 - 1e-10, 1e-150)];
    let expected =
        [2.178_248_840_713_436_5e279, 5.000_000_002_299_904e199, 7.522_527_521_252_041e149];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hypu_chgus_gamma_poles() {
    // Arb; 0.850016, 0.729927, 3.86817017062829e-23, 8.2206319e32, and 1.3849818e92 in xsf
    // (Gamma(b) or Gamma(2 - b) near a pole)
    let abx = [
        (0.3, -1.0 + 1e-14, 1e-100),
        (-0.15, -1.0 + 2e-16, 1e-280),
        (22.0, -1.0 - 3e-15, 1.6e-14),
        (20.0, 3.0 - 1e-9, 1e-25),
        (24.4, 6.0 - 7e-10, 1.8e-23),
    ];
    let expected = [
        0.857_109_621_959_464_5,
        1.057_517_111_631_981_5,
        3.868_170_170_629_289_4e-23,
        8.220_634_769_930_943e32,
        1.384_981_943_556_528_2e92,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-15);
}

#[test]
fn test_hypu_chgus_rounding() {
    // Arb; 1.381478037678 and 3.443455701978e64 in xsf. Here chgus has fewer digits than it would
    // claim without the rounding of 1 + a - b (first) and of its terms (second).
    let abx = [
        (0.001_540_284_073_910_872_3, 0.999_999_968_212_571_8, 3.409_133_030_522_77e-108),
        (-50.3, 1.001, 0.999_999_999_999_999_9),
    ];
    let expected = [1.381_478_037_771_745_5, 3.443_455_703_062_540_5e64];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-14);
}

#[test]
fn test_hypu_chguit_unconverged() {
    // Arb; here the quadrature in chguit doesn't converge, so that chgu keeps chgus
    let abx = [(25.5, -1e-14, 1e-10), (25.5, 1e-10, 1e-6), (25.5, 2.000_000_000_1, 1e-6)];
    let expected =
        [1.270_443_779_292_559_8e-26, 1.270_106_821_887_141e-26, 3.238_802_988_847_028e-19];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-9);
}

#[test]
fn test_hypu_chguit_overflow() {
    // Arb; 2.11429e88, 2.80252e44, 5.5268225094e227, NaN, and NaN in xsf (inf before the
    // integral without x^(1 - b) in `chguit`, as Gamma(a) U, or t^(a-1) in it, overflows)
    let abx = [
        (20.0, 7.0 - 2e-11, 6e-18),
        (21.7, 5.0 + 2e-11, 1.8e-16),
        (13.2, 8.0 - 1.2e-7, 4e-34),
        (29.4, 85.76, 0.0032),
        (26.8, 16.0 - 4e-14, 3.3e-22),
    ];
    let expected = [
        2.114_360_915_592_612e88,
        2.802_512_228_456_701_7e44,
        5.526_822_512_984_652e227,
        2.852_274_462_640_627_7e307,
        6.940_743_519_297_277e306,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hypu_b_near_zero() {
    // mpmath (scipy/xsf#533)
    let abx = [(2.5, -1e-12, 15.0), (0.75, -1e-13, 3.0), (2.5, 1e-12, 15.0)];
    let expected =
        [0.000_707_019_879_343_321_6, 0.324_157_202_834_937_1, 0.000_707_019_879_343_501_9];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hypu_b_near_one() {
    // Arb; 45.474459, 258.47514, 231.697, -2.911784, -169.7466186, 124.033632, and 1.0000000005
    // in xsf
    let abx = [
        (1.0 - 1e-10, 1.0 - 1e-10, 1e-20),
        (1.5, 1.0 + 1e-10, 1e-100),
        (0.3, 1.0 + 1e-12, 1e-300),
        (-0.5, 1.0 - 1e-12, 1e-5),
        (-2.8125, 1.000_001, 3.5e-87),
        (3.5, 1.0 - 5e-9, 1e-180),
        (2.1e-11, 1.000_000_15, 3.2e-25),
    ];
    let expected = [
        45.474_486_089_040_32,
        258.475_097_127_408_1,
        231.692_340_668_562_74,
        -2.911_763_413_684_127_5,
        -169.746_618_603_492_12,
        124.033_634_780_699_4,
        1.000_000_001_184_435_9,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hypu_b_near_one_recurrence() {
    // Arb; through the recurrence in a, where the relative condition number of U in a is 6e6
    let actual = xsf::hypu(-5.1, 0.999_999_999_999, 0.329_999);
    xsref::np_assert_allclose!([actual], [-0.000_200_046_396_593_785_78], rtol = 1e-9);
}

#[test]
fn test_hypu_b_near_one_cancel() {
    // Arb; 0.00354666023, -0.00410532296, -0.192757095, and 7.119e172 in xsf (the two series of
    // DLMF 13.2.42 cancel, here also near zeros of U)
    let abx = [
        (-0.056_635_517_479_212_36, 1.000_000_039_777_816_4, 1.417_381_804_653_536_6e-8),
        (-0.063_925_708_762_279_1, 0.999_999_999_035_314_8, 9.429_167_746_765_604e-8),
        (-4.083_168_773_834_823, 0.999_987_773_001_151_6, 4.969_885_307_388_012e-7),
    ];
    let expected = [0.003_546_661_069_829_69, -0.004_105_300_291_470_05, -0.192_757_094_020_925_53];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    let actual =
        xsf::hypu(-108.048_671_397_557_14, 0.999_999_999_949_052_2, 9.822_176_387_052_842e-12);
    xsref::np_assert_allclose!([actual], [7.106_291_122_209_209e172], rtol = 3e-12);

    // Arb; near a pole of Gamma(a)
    let actual = xsf::hypu(-168.02, 0.9983, 1e-38);
    xsref::np_assert_allclose!([actual], [-1.435_814_008_821_082_3e302], rtol = 5e-15);

    // Arb; at a pole of Gamma(a), where a - (b - 1) rounds to a
    let actual = xsf::hypu(-168.0, 1.000_000_000_000_015, 1e-38);
    xsref::np_assert_allclose!([actual], [2.526_075_744_973_415_8e302], rtol = 5e-15);

    // Arb; at large |a| and |b - 1|
    let actual = xsf::hypu(-166.125, 1.0099, 1e-4);
    xsref::np_assert_allclose!([actual], [9.921_163_464_098_675e297], rtol = 5e-15);

    // Arb; 1e-2 <= |b - 1| < 2e-2, where the series cancel too (1.0e-13, 1.9e-13 and 1.5e-12 off
    // in xsf)
    let abx = [
        (-9.02, 0.9883, 1e-11),
        (-19.999_999_999_512_784, 0.982_396_189_761_965_6, 9.401_094_988_206_325e-8),
        (-152.22, 1.0118, 3.75e-4),
    ];
    let expected = [-224_960.382_251_122_97, 2.283_020_666_022_678e18, 1.734_174_831_273_375e267];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-15);

    // Arb; the other methods are more accurate here
    let actual =
        [xsf::hypu(165.0, 1.0099, 0.01), xsf::hypu(16.5, 0.994, 0.1), xsf::hypu(0.3, 0.9955, 6.7)];
    let expected =
        [3.703_508_884_277_675e-295, 2.397_573_509_052_701_8e-14, 0.558_249_828_304_120_8];
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-15);

    // Arb; the recurrence in a has 10 digits, the sum 13 (but estimates 11)
    let actual =
        xsf::hypu(-152.571_517_551_837_93, 0.997_562_922_686_568_3, 0.001_997_420_588_184_022);
    xsref::np_assert_allclose!([actual], [1.190_151_804_933_915_8e266], rtol = 1e-13);

    // Arb: 2.79e-307, but 1/Gamma(a) underflows (and the fallback returns 0)
    let actual = xsf::hypu(171.625, 1.009, 1e-20);
    assert!((0.0..=2.790_460_095_786_538_4e-307).contains(&actual));
}

#[test]
fn test_hypu_b_one_small_a() {
    // Arb; 0.0121547807106452, -0.147635975145256, and -3.5403980739e-5 in xsf (DLMF 13.3.7
    // cancels near the zeros of U)
    let ax = [
        (-0.040_066_081_206_398_19, 1.186_782_123_234_391_8e-11),
        (-0.031_365_424_999_681_31, 6.961_552_106_224_565e-17),
    ];
    let expected = [0.012_154_780_710_647_068, -0.147_635_975_145_240_85];
    let actual = ax.map(|(a, x)| xsf::hypu(a, 1.0, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 3e-14);

    // Arb; the relative condition numbers of U in a are 3e4 and 7e4
    let ax = [
        (-0.017_884_036_913_673_778, 3.002_462_249_449_102e-25),
        (-0.038_610_143_988_774_756, 3.383_438_104_338_344e-12),
    ];
    let expected = [-3.540_398_072_921_225e-5, -1.325_756_475_271_633_5e-5];
    let actual = ax.map(|(a, x)| xsf::hypu(a, 1.0, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-11);
}

#[test]
fn test_hypu_b_near_int() {
    // Arb; 0.084043687, 1.3983195, 0.0034193771, 3.4728082e-7, and 7.3899616e-6 in xsf
    let abx = [
        (11.5, 2.0 + 2e-13, 1e-6),
        (10.3, 2.0 - 1e-14, 1e-6),
        (5.5, 1e-14, 5e-4),
        (8.9, -1.0 - 3e-12, 4e-6),
        (6.9, -2.0 - 3e-15, 6e-6),
    ];
    let expected = [
        0.084_027_749_274_997_87,
        1.395_660_734_247_923_6,
        0.003_419_492_016_194_043_8,
        3.484_572_598_279_22e-7,
        6.899_556_502_797_434e-6,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 59304.064, 1.4934268e262, and 6.8120873e74 in xsf (NaN before the recurrence in b)
    let abx =
        [(32.0, 4.0 - 6e-14, 1.6e-13), (25.0, 4.0 + 5e-13, 6e-96), (90.0, 3.0 + 3e-14, 3e-106)];
    let expected = [59_381.105_544_359_43, 1.492_349_603_031_082_9e262, 6.730_761_625_464_286e74];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-14);
}

#[test]
fn test_hypu_large_a() {
    // mpmath; NaN in xsf
    let actual = xsf::hypu(79.2, 8.12, 3.28);
    xsref::np_assert_allclose!([actual], [2.837_981_244_377_152_2e-124], rtol = 1e-13);

    // Arb; by the asymptotic series, where chgus and chguit overflow, or chgus has too few digits
    let actual = [
        xsf::hypu(172.0, 173.000_000_000_000_03, 50.0),
        xsf::hypu(171.0, 172.000_000_000_000_03, 50.0),
    ];
    let expected = [5.986_310_706_507_632e-293, 2.993_155_353_253_815_5e-291];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-11);
}

#[test]
fn test_hypu_large_b() {
    // mpmath; 2.84e-8 and -614.4 in xsf
    let abx = [(3.69, -87.59, 20.53), (-1.28, 70.68, 126.76)];
    let expected = [2.845_728_132_149_186_6e-8, 171.372_404_205_867_78];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-11);

    // Arb; 6.526, -0.04274, and NaN in xsf, where x^(1 - b) underflows in `chgus`
    let abx = [(0.0013, 148.6, 205.0), (0.94, 165.43, 139.13), (-16.24, 143.86, 191.27)];
    let expected = [0.994_721_980_416_835_8, 1.829_378_809_443_741_5, -9.510_439_577_470_078e24];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 1.0000000026, 1.00000036, 1.000001, and 0.9999999976 in xsf, where the terms of the
    // second series in `chgus` grow until j ~ x
    let abx = [
        (-5e-10, 162.4, 94.0),
        (-1e-7, 167.5, 133.0),
        (-4e-7, 171.5, 163.6),
        (-1.24e-13, 143.69, 204.76),
    ];
    let expected = [
        0.898_935_414_508_470_5,
        0.999_994_633_555_940_3,
        1.000_000_405_055_020_2,
        1.000_000_000_000_514,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 6.52e-237, 2.6069004e-98, 1.58498e-174, 1.466e-196, and -1.207e-253 in xsf, where the
    // terms of the second series in `chgus` grow again after j = b - 1
    let abx = [
        (149.6, 146.67, 39.25),
        (73.61, 79.99, 23.51),
        (119.57, 116.57, 26.86),
        (123.77, 128.58, 42.47),
        (150.8, 157.76, 53.95),
    ];
    let expected = [
        7.905_260_961_079_413e-242,
        2.606_902_923_039_754_6e-98,
        1.584_879_429_039_851e-174,
        5.762_982_455_359_295e-200,
        1.954_035_627_593_264e-258,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 1.0e-282 and NaN in xsf, where the r2 terms of `chgus` underflow before they grow again
    let abx = [(164.9, 168.900_01, 55.0), (164.9, 168.91, 57.5)];
    let expected = [6.671_445_339_859_178e-286, 4.012_537_078_507_572_6e-289];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 1.029_625_838e-306, 7.489_818_350e-304, 4.094_197_594e-310, 5.499e-312, and 1.77e-314 in
    // xsf, where x^(1 - b) in `chgus` is subnormal (the last three correctly rounded)
    let actual = [xsf::hypu(145.5, 150.5, 129.43), xsf::hypu(124.5, 129.5, 275.4)];
    let expected = [1.029_625_833_883_920_1e-306, 7.489_818_349_462_037e-304];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
    let actual = [
        xsf::hypu(100.5, 103.5, 1200.0),
        xsf::hypu(136.5, 141.5, 193.68),
        xsf::hypu(165.5, 170.5, 81.07),
    ];
    assert_eq!(actual, [4.094_197_596_010_1e-310, 5.531_686_062_314e-312, 1.054_496_25e-314]);
}

#[test]
fn test_hypu_large_neg_b() {
    // mpmath; 25.65 and 2.203e6 in xsf
    let abx = [(0.2, -138.4, 40.4), (-0.4, -141.6, 43.9)];
    let expected = [0.354_135_821_455_708_52, 8.085_089_469_342_279];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hypu_large_neg_b_gamma_overflow() {
    // Arb; 1.498e-226, 3.47e-287, -2.36e-303, -4.65e-207, and -2.96e-314 in xsf, where
    // Gamma(1 + a - b) overflows in `chgus`
    let abx =
        [(115.7, -62.2, 4.86), (142.4, -28.5, 0.39), (115.2, -67.99, 0.51), (96.94, -85.08, 6.1)];
    let expected = [
        5.280_753_643_679_679e-243,
        3.384_055_764_610_558_4e-280,
        3.099_375_276_064_330_5e-241,
        8.380_400_275_246_349e-209,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
    // (correctly rounded)
    assert_eq!(xsf::hypu(163.0, -15.99, 0.18), 3.405_301_615e-315);
}

#[test]
fn test_hypu_kummer_overflow() {
    // Arb; NaN or inf in xsf, where U(a - b + 1, 2 - b, x) or x^(1 - b) over- or underflows
    let abx = [
        (-7.0, -7.0, 1e-187),
        (-1.25, -4.0, 1e-219),
        (0.5, -100.0, 0.01),
        (-107.0, -173.0, 0.006),
        (-0.5, -300.0, 2.0),
        (0.25, -250.0, 100.0),
        (-0.0015, -170.625, 1.5),
    ];
    let expected = [
        5040.0,
        5.426_195_543_094_943,
        0.099_621_961_951_682_32,
        6.809_106_301_661_743e220,
        17.385_389_278_273_84,
        0.231_064_916_158_941_04,
        1.007_756_643_938_103_8,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-12);

    // Arb; still from `chgus`, as the integral misses the near-singularity at s = x if 0 < b < 1
    let actual = xsf::hypu(0.999_999_99, 0.999_999_99, 1e-200);
    xsref::np_assert_allclose!([actual], [459.938_742_565_798_5], rtol = 1e-8);
}

#[test]
fn test_hypu_int_b() {
    // Arb; NaN, NaN, NaN, and 4.534890e-44 (from `chgubi`) in xsf
    let abx = [(13.5, 25.0, 0.005), (1.0, -26.0, 0.02), (25.0, 2.0, 5.0), (30.0, -10.0, 2.0)];
    let expected = [
        2.541_392_992_138_2e68,
        0.037_008_569_781_592_95,
        5.851_475_774_467_137e-33,
        4.534_945_524_706_779e-44,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb; 20.24, 1.0210079, -2.5019, -11.2, and 4.17144e16 in xsf (`a - n`, `gamma2`, and
    // `psi_spec` lose digits near the poles of Gamma(a))
    let abx = [
        (-2.0 + 1e-15, 4.0, 0.157),
        (2.1e-10, 6.0, 0.048),
        (-1.0 + 1e-13, 3.0, 0.5),
        (-1.0 + 1e-14, 12.0, 2.4),
        (-16.0 - 1e-14, -16.0, 7.6),
    ];
    let expected = [
        18.454_649_000_001_735,
        1.021_005_617_856_604_2,
        -2.500_000_000_000_773_5,
        -9.600_000_000_047_481,
        4.171_421_288_044_156e16,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; 1.00502, 1.00161, and 0.999999997 in xsf (x^-n underflows in `chgubi`)
    let abx = [(-1e-3, 160.0, 150.0), (-3e-4, 165.0, 212.0), (-2e-8, 152.0, 139.0)];
    let expected = [1.000_521_817_125_236_8, 1.001_171_575_012_627, 0.999_999_992_869_047_8];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-14);

    // Arb; 5.9e-4, 0.4382008983836371, and 0.5039275260563122 in xsf (`a - n` loses the low
    // digits of `a` in `chgubi`)
    let abx = [(0.936, 150.0, 170.0), (0.3, 140.0, 150.0), (0.43, 111.0, 106.9)];
    let expected = [0.046_627_861_862_793_184, 0.438_200_898_383_651_9, 0.503_927_526_056_324_3];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 8e-15);

    // Arb; 0, inf, 0, 0, and NaN in xsf (n! overflows in `chgubi` for n > 170)
    let abx = [
        (0.936, 172.0, 197.0),
        (1.5, 173.0, 3.0),
        (5.5, 172.0, 150.0),
        (1.5, -170.0, 2.0),
        (1.5, -175.0, 2.0),
    ];
    let expected = [
        0.039_699_178_704_971_555,
        2.401_138_892_408_052_4e228,
        3.952_873_928_825_523e-6,
        4.384_667_341_835_304_7e-4,
        4.201_507_804_107_285_5e-4,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-14);

    // Arb; NaN in xsf, where x^-n and hm3 in `chgubi` under- and overflow
    let actual = xsf::hypu(15.000_000_000_013_445, 864.0, 721.784_972_617_366_9);
    xsref::np_assert_allclose!([actual], [2.869_995_870_011_826_5e-17], rtol = 1e-14);

    // Arb; the over- and underflows aren't NaN
    assert_eq!(xsf::hypu(171.61, 1.0, 2.0), 0.0);
    assert_eq!(xsf::hypu(-3.0 + 1e-15, 999.0, 1e-300), f64::NEG_INFINITY);
    assert_eq!(xsf::hypu(-1e10, 10.0, 1e-300), f64::INFINITY);
    assert_eq!(xsf::hypu(-169.0 - 3e-14, 14.0, 4.0), f64::INFINITY);
    assert_eq!(xsf::hypu(-167.0 - 3e-14, -16.0, 385.0), f64::NEG_INFINITY);

    // Arb; still from `chgubi` (11 digits) if `chguit` overflows
    let actual = xsf::hypu(150.0, 9.0, 0.1);
    xsref::np_assert_allclose!([actual], [2.229_672_231_677_937e-250], rtol = 1e-12);

    // mpmath; near a zero of U, where `chgubi` has as many estimated digits as the recurrence in
    // `a` (in `b` for the last)
    let abx = [
        (-0.271_543_086_172_344_7, 1.0, 0.023_848_456_913_827_66),
        (-0.379_759_519_038_076_14, 1.0, 0.084_368_893_787_575_16),
        (-41.1575, 1.0, 2.154_434_690_031_882_3e-5),
        (-28.6575, 4.0, 0.215_443_469_003_188_67),
    ];
    let expected = [
        -7.824_866_488_540_172e-4,
        3.360_727_288_593_018e-4,
        3.386_947_244_049_369e47,
        1.089_830_835_818_625e30,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-12);
}

#[test]
fn test_hypu_int_b_overflow() {
    // Arb; where x^n, the sum that it multiplies in `chgubi`, or Gamma(a) (for the last four) over-
    // or underflow
    let abx = [
        (150.0, 3.0, 1e-200),
        (170.0, 2.0, 1e-310),
        (5.5, 999.0, 355.0),
        (170.5, 701.0, 350.0),
        (172.5, 2.0, 1e-310),
        (172.5, 701.0, 350.0),
        (-172.5, -150.0, 0.03),
        (-175.5, -140.0, 0.01),
    ];
    let expected = [
        2.625_414_310_389_023e139,
        234_243.164_524_601_7,
        1.187_305_823_643_144e165,
        7.685_813_590_383_2e-294,
        0.614_855_825_106_415_4,
        8.708_663_934_780_86e-299,
        4.485_540_746_959_017_7e282,
        -7.551_826_164_911_572e279,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb (5.559_687_150_9e-319); a subnormal result, of subnormal terms
    assert_eq!(xsf::hypu(170.5, -5.0, 0.1), 5.559_67e-319);

    // Arb; `chgubi` has 11 digits, and the log form of `chguit` 5
    let actual = xsf::hypu(179.975, 16.0, 0.2177);
    xsref::np_assert_allclose!([actual], [7.264_974_891_266_849e-308], rtol = 1e-12);
}

xsref::tests! {
    test_hyp1f1: "hyp1f1", "d_d_cd-cd", |x| xsf::hyp1f1(x[0], x[1], c64(x[2], x[3]));
    test_pmv: "pmv", "d_d_d-d", |x| xsref::with_int(x[0], |m| xsf::pmv(m, x[1], x[2]));
}

#[test]
fn test_pmv_large_m() {
    // NaN in xsf, which overflows `-m`
    assert!(xsf::pmv(i32::MIN, 0.5, 0.5).is_nan());
    // no `int` overflow of `j * j` in `lpmv0`
    let _ = xsf::pmv(50_000, 0.5, 0.5);
}

#[test]
fn test_hyp1f1_huge_a() {
    // `(int)a` is INT_MIN like in xsf, instead of ~2^31 iterations
    let _ = xsf::hyp1f1(3e9, 1.0, 1.0);
    let _ = xsf::hyp1f1(3e9, 1.0, c64(1.0, 0.0));
}

#[test]
fn test_hyp1f1_pole_b_zero() {
    // unlike in xsf: these took ~a steps (or were inf for the last two)
    assert!(xsf::hyp1f1(1.13e9, 0.0, 1.0).is_nan());
    assert!(xsf::hyp1f1(1.13e9, -0.0, 0.0).is_nan());
    assert!(xsf::hyp1f1(-1.13e9, 0.0, -1.0).is_nan());
    assert!(xsf::hyp1f1(2.0, 0.0, 31.0).is_nan());
    assert!(xsf::hyp1f1(3e9, 0.0, 1.0).is_nan());
}

#[test]
fn test_hyp1f1_huge_b_minus_a() {
    // Arb; these took ~(b - a) steps (up to 21 s), and were 0
    let abx = [(0.5, 1e9, -1e8), (100.3, 2e8, -2e6), (3.7, 2.1e9, -1e9), (-3.7, 1e7, -5e6)];
    let expected = [
        0.953_462_589_242_637_4,
        0.368_609_232_323_196_9,
        0.236_686_301_386_977_8,
        4.482_678_936_020_103,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb: 0.9999999990000010005 and 0.99999999996849996850; these were NaN and 1
    let actual = [
        xsf::hyp1f1(1_000_001.5, -1_000_002.5, 1e-9),
        xsf::hyp1f1(1_000_001.0, -4e15 - 0.5, 0.126),
    ];
    let expected = [0.999_999_999_000_001, 0.999_999_999_968_5];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);
    // as before
    assert_eq!(xsf::hyp1f1(1_000_001.0, -4e15 - 0.5, 2e-7), 1.0);
    assert_eq!(xsf::hyp1f1(1_000_001.0, -1.0, 0.0), 1.0);

    // DLMF 13.6.3 and 13.6.1; these took ~1e9 steps, and were NaN
    assert_eq!(xsf::hyp1f1(-0.0, 6.3e8, -4.6e10), 1.0);
    assert_eq!(xsf::hyp1f1(1e9, 1e9, f64::INFINITY), f64::INFINITY);
    // these took ~1e9 steps; mpmath: -3.7e30130, and this was NaN
    assert_eq!(xsf::hyp1f1(-1.2e9, -10.2, -1.0), f64::NEG_INFINITY);
    assert!(xsf::hyp1f1(1e9, 2.5, f64::NAN).is_nan());
}

#[test]
fn test_hyp1f1_neg_x_no_convergence() {
    // Arb; the series after DLMF 13.2.39 doesn't converge within 500 terms (b + 500 < 0), and these
    // were e^-x times too large
    let abx =
        [(2e7 + 0.3, -2e5 - 0.3, -1.0), (10.3, -1e6 - 0.5, -100.0), (-1000.3, -2e5 - 0.3, -100.0)];
    let expected = [2.756_493_631_860_608e43, 1.001_030_581_673_679_5, 0.606_364_673_033_219_6];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp1f1_beyond_500_terms() {
    // Arb; the series needs over 500 terms (b < 0), and these were 0.854 and NaN
    let abx = [(2.61, -2_880_000.5, -477.5), (2.5, -50_000.5, -750.0)];
    let expected = [1.000_432_863_835_773_4, 1.038_506_615_033_215_6];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hyp1f1_truncated() {
    // Arb: -1.64e183 and 9.31e-145; 500 terms are too few, and these were 1.02e266 and 1.39e126
    assert!(xsf::hyp1f1(-2e6 - 0.3, -600.5, 0.1).is_nan());
    // Arb: 2.7e84592; the terms shrink, but grow again near -b
    assert!(xsf::hyp1f1(-800.3, -100_000.5, 100_000.5).is_nan());
    assert!(xsf::hyp1f1(-2e6 - 0.3, 2_880_000.5, 477.5).is_nan());
}

#[test]
fn test_hyp1f1_large_a() {
    // mpmath; in xsf these are 1.2e87, 3.3e251, 7.9e303, inf, and NaN
    let abx = [
        (300.0, 400.0, 1e-10),
        (952.25, 952.24, 5.94),
        (57.05, 1.1e7, 1.86),
        (5.0, 2000.0, -650.0),
        (2.5, 3000.0, -600.0),
    ];
    let expected = [
        1.000_000_000_075,
        379.958_556_784_186_2,
        1.000_009_646_683_708_2,
        0.244_751_044_949_691_12,
        0.633_912_466_368_326_7,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp1f1_neg_int_a() {
    // mpmath; in xsf these are 5.4e29, -1.5e24, -1.3387308e-8, and 4.0e63
    let abx = [
        (-1000.0, 0.5, 3.0),
        (-100.0, 0.5, 30.0),
        (-7.0, 950.0, 850.0),
        (-200.0, 0.001, 50.0),
        // `p * x` overflows in the recurrence, but the series is fine
        (-3.0, 5e307, 1.7e308),
    ];
    let expected = [
        -4.108_819_740_932_209,
        288_862.419_837_060_44,
        -1.338_730_683_179_536_7e-8,
        -297_046_447_677_026.5,
        -13.824,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hyp1f1_complex_large_a() {
    // mpmath; in xsf these have a relative error of 3.5e3, 1.3e-6, 9.7e72, and 8.6e248
    let abz = [
        (30.0, 70.0, c64(0.0, 20.0)),
        (10.0, 50.0, c64(5.0, 30.0)),
        (0.5, 400.0, c64(-30.0, -1.0)),
        (952.25, 952.24, c64(5.94, 0.0)),
    ];
    let expected = [
        c64(-0.320_661_290_060_071_9, 0.381_802_776_457_525_86),
        c64(0.602_663_421_239_993_6, -0.092_796_424_643_608_8),
        c64(0.964_479_298_872_970_2, -0.001_121_758_543_499_657_6),
        c64(379.958_556_784_186_2, 0.0),
    ];
    let actual = abz.map(|(a, b, z)| xsf::hyp1f1(a, b, z));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-11);
}

#[test]
fn test_hyp1f1_complex_asymptotic() {
    // mpmath; scipy/xsf#398 (~1e-8 because the asymptotic expansion is truncated)
    let abz =
        [(0.5, 1.5, c64(5.0, -25.0)), (1.25, 2.5, c64(-11.0, 20.0)), (0.25, 3.5, c64(-2.0, 30.0))];
    let expected = [
        c64(0.242_995_298_734_033_47, 2.780_418_451_436_882),
        c64(0.007_180_913_061_977_942_4, 0.028_317_020_033_024_477),
        c64(0.520_732_899_622_106_7, 0.194_475_861_451_502_78),
    ];
    let actual = abz.map(|(a, b, z)| xsf::hyp1f1(a, b, z));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-8);
    let actual_conj = abz.map(|(a, b, z)| xsf::hyp1f1(a, b, z.conj()).conj());
    assert_eq!(actual_conj, actual);
}

#[test]
fn test_hyp1f1_neg_a() {
    // mpmath; in xsf these are 5.5e8, -9.7e28, -2.6, -0.037, -1.5e15, and NaN
    let abx = [
        (-1000.5, 1.0, 1.0),
        (-1000.5, 0.5, 3.0),
        (-1000.5, 0.5, 0.4),
        (1100.5, 100.0, -5.0),
        (-100.7, 10_000.0, 11_000.0),
        (-20_000.25, 3.0, 50.0),
    ];
    let expected = [
        0.155_652_734_227_881_56,
        -4.156_338_175_187_512,
        -0.827_579_317_701_538,
        1.667_059_039_325_649_8e-31,
        4.939_183_825_283_591e-113,
        -1_296.785_197_037_247_8,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-11);
}

#[test]
fn test_hyp1f1_neg_a_cancellation() {
    // mpmath; in xsf these have a relative error of 9.6e-10, 1.3e-9, 2.7e-9, and 2.9e-10
    let abx = [
        (-14.5636, 15.2965, 12.1738),
        (-10.4812, 28.7818, 52.9447),
        (-420.49, 2.1344, 0.20916),
        (38.8876, 25.9916, -17.909),
    ];
    let expected = [
        4.568_079_699_381_592e-5,
        5.599_261_569_881_319e-4,
        -0.014_790_768_773_875_595,
        -1.752_305_173_630_152e-13,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hyp1f1_huge_neg_a() {
    // Arb; DLMF 13.2.39 gives a polynomial of degree 10^5 < a - b <= 10^6
    let abx = [(3e5, 1.0, -0.01), (3e5, 30.0, -1.0), (1e6, 30.0, -0.1)];
    let expected =
        [-0.027_636_911_294_932_15, 4.867_432_831_661_024e-51, -6.454_678_759_269_63e-44];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-12);
}

#[test]
fn test_hyp1f1_poly_underflow() {
    // Arb, rounded to the nearest subnormal
    let abx = [(-1e5, 300.0, 15.0), (-1e5, 300.0, 17.0)];
    let expected = [5.488_687_915_453e-311, 8.469_77e-319];
    assert_eq!(abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x)), expected);
    // Arb: -1.3e-415 and -4.8e-459, which round to -0
    let abx = [(-1e5, 300.0, 100.0), (100_300.0, 300.0, -100.0)];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x).to_bits());
    assert_eq!(actual, [(-0.0_f64).to_bits(); 2]);
}

#[test]
fn test_hyp1f1_huge_neg_a_tricomi() {
    // Arb; the series cancels, and the recurrences in a lose digits or take over 10^6 steps
    // (condition numbers 1e3, 1e2, 4e2, 1e2, 2e5, and 3e3)
    let abx = [
        (-999_999.8, 0.5, 0.1),
        (1e6 + 0.3, 1.0, -0.1),
        (100_000.3, 0.5, -0.1),
        (-100_000.3, 300.0, 0.1),
        (1e6 + 0.3, 300.0, -0.1),
        (-2.3, 1e-20, 4.343_478_260_869_565e-21),
        // integer `a` past the 10^6 steps of the polynomial
        (-1e7, 1.0, 1.0),
        (-2_147_483_647.0, 1.0, 1e-6),
    ];
    let expected = [
        -0.572_048_140_852_985_2,
        -0.029_517_956_327_592_655,
        0.463_461_893_294_257_1,
        3.699_009_173_383_038e-16,
        -3.953_672_642_019_636_4e-140,
        0.001_000_000_000_000_050_2,
        -0.016_002_420_737_843_012,
        -0.058_394_171_687_540_5,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-10);
}

#[test]
fn test_hyp1f1_huge_neg_a_tricomi_scaling() {
    // Arb; `b` far above `t = sqrt(2 (b - 2a) x)`, where the `J_v` of Tricomi's expansion underflow
    // (condition numbers 1e2, 45, and 3e2)
    let abx =
        [(1e6 + 0.3, 3000.0, -0.1), (-299_989.44, 801.9, 0.039_59), (-647_118.44, 1839.0, 0.2452)];
    let expected =
        [2.767_870_543_528_383e-15, 3.212_622_595_478_475_4e-7, 3.858_006_441_225_035e-39];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb; the prefactor `Gamma(b) (t/2)^(1-b)` of Tricomi's expansion for large and tiny `b`
    // (condition numbers 12 and 29)
    let abx = [(-3000.3, 79.07, 0.1), (-79.5, 1e-300, 2.0)];
    let expected = [0.020_380_553_037_493_01, 3.559_294_019_842_940_4e300];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb; where only that prefactor (times `e^(x/2)`) overflows, for tiny `b` and large `x`
    // (condition numbers 29, 2e2, and 3e6)
    let abx = [(-79.5, 1e-307, 2.0), (10_000.3, 5e-308, -0.1)];
    let expected = [3.559_294_019_842_941e307, 2.153_758_520_051_220_8e307];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
    let actual = xsf::hyp1f1(-100_000_000.5, 0.5, 1420.0);
    xsref::np_assert_allclose!([actual], [-4.922_030_132_579_361e307], rtol = 1e-9);
}

#[test]
fn test_hyp1f1_large_x() {
    // mpmath; in xsf these are -1.0e286, -1.6e194, 1788, NaN, 6.4e214, -2.56741098e157,
    // -1.2e252, accurate, accurate, and 1 (stopping at the negligible first term)
    let abx = [
        (-0.5, 0.5, 700.0),
        (-0.25, 10.0, 500.0),
        (0.5, 100.0, 140.0),
        (0.5, 200.0, 545.0),
        (-0.5, -10.5, 450.0),
        (10.5, -50.25, 200.0),
        (1.2, 0.7, -600.0),
        (1e-20, 0.5, 100.0),
        (1e-20, -0.5, 100.0),
        (-1e-20, 1000.0, 1500.0),
    ];
    let expected = [
        -7.260_094_358_424_212e300,
        -2.309_342_355_757_666e194,
        1_839.818_976_462_603_8,
        1.682_375_001_297_950_4e63,
        6.642_084_214_977_351e214,
        -2.567_411_032_674_188e157,
        -1.702_996_271_829_136_5e-4,
        4.788_753_000_993_948e22,
        -9.481_239_629_384_337e24,
        -2.716_385_691_893_235e20,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp1f1_large_neg_x() {
    // Arb; by the asymptotic expansion (DLMF 13.7.2), and these were off by 1.5e-12, 2.2e-10, 9.4e-9,
    // NaN, NaN, and 2.9e-9 (`b < 0` with `b - a > 0`)
    let abx = [
        (2.0, 0.5, -200_000.0),
        (5.3, 2.7, -300_000.0),
        (3.56, 0.57, -960_000.0),
        (1.5, 1.0, -100_000.0),
        (0.7, 0.2, -500_000.0),
        (-0.031_25, -0.015_625, -220_344.75),
    ];
    let expected = [
        1.875_046_876_230_505_5e-11,
        -1.626_960_239_562_368e-29,
        -4.672_030_103_421_194_5e-23,
        -8.920_821_300_999_49e-9,
        -1.327_423_448_749_78e-4,
        -1.495_570_968_154_628,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);
}

#[test]
fn test_hyp1f1_neg_x_near_pole() {
    // Arb; b - a is within 1e-14 of -31, -14, and -82 (where it rounds to -82), 1e-20 and 6e-17
    // (above) of -1, and -1 - 8.3e-17 (rounding to -1) for the last two, i.e. near a zero of
    // 1/Gamma(b - a), and these were -3.9e-75, -8.0e-118, -1.3e-153, 0, 0, -0.0, and -0.0
    let abx = [
        (31.7, 0.699_999_999_999_988_6, -1000.0),
        (48.277_939_327_037_45, 34.277_939_327_037_46, -1_450.237_668_835_392_7),
        (21.553_847_421_355_318, -60.446_152_578_644_686, -646.529_513_537_224),
        (1.0, -1e-20, -1000.0),
        (0.3, -0.7, -50_000.0),
        (1.1, 0.1, -1e5),
        (1.1, 0.1, -1e6),
    ];
    let expected = [
        2.580_045_808_721_416_3e-75,
        5.612_220_190_937_725e-119,
        5.634_897_230_748_735e-34,
        -1.002_006_024_120_725_1e-3,
        9.236_244_627_873_464e-18,
        2.505_080_704_702_706_7e-21,
        1.989_816_934_354_528_4e-22,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hyp1f1_neg_b() {
    // mpmath; in xsf these are 1.03, -7844.06, -1.5e56, NaN, NaN, and NaN
    let abx = [
        (-0.0942, -85.783, 31.993),
        (20.0, -17.75, -36.75),
        (536.5, -66.75, -7.75),
        (874.25, -334.75, -138.5),
        (900.0, -830.5, 0.25),
        (737.5, -61.75, 518.0),
    ];
    let expected = [
        -830_724_968_298.684,
        -7_845.064_845_794_902,
        -1.334_256_059_470_499_8e26,
        -1.922_004_000_563_705_4e133,
        0.762_741_486_096_305,
        f64::INFINITY,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-12);
}

#[test]
fn test_hyp1f1_neg_b_cancellation() {
    // mpmath; in xsf these are 4.6e169, -0.0, and NaN
    let abx = [(500.5, -200.25, 5.5), (99.7, -0.3, -800.0), (340.5, -90.25, 140.5)];
    let expected =
        [1.952_103_586_824_771_7e-6, 9.487_619_910_475_397e-140, -6.941_454_309_937_287e307];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
}

#[test]
fn test_hyp1f1_neg_b_rounded_b_minus_a() {
    // Arb; `b - a` rounds to (or near) an integer, which changes the terms of DLMF 13.2.39 a lot
    let abx = [
        (0.037, -0.963, -25.2),
        (0.037, -1.963, -25.2),
        (0.109, -5.891_000_001, -30.0),
        (0.759, -0.241_000_000_001, -29.0),
        // (an integer far beyond the terms that are summed)
        (1e150, -1e-200, -1e-150),
    ];
    let expected = [
        3.089_167_835_910_344e-10,
        4.123_043_467_765_241e-9,
        1.430_320_106_089_673e-7,
        3.044_924_220_068_411_5e-11,
        5.767_248_077_568_734e199,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // here the series in double-double is only accurate to 1e-11, so it shouldn't precede
    // `hyp1f1_neg_b`
    let actual =
        xsf::hyp1f1(44.999_946_731_604_33, -5.326_409_422_435_637_5e-5, -52.614_841_190_696_32);
    xsref::np_assert_allclose!([actual], [-1.332_909_333_817_075_8e-7], rtol = 1e-13);

    // the sum overflows, and is redone with `e^(x - 1400)`, where `x - 1400` rounds
    let actual = xsf::hyp1f1(10_000.0, -1e-309, -1e-7);
    xsref::np_assert_allclose!([actual], [9.995_000_333_513_848e305], rtol = 1e-15);

    // and here that redone sum underflows to 0
    let actual = xsf::hyp1f1(1200.0, -60.3, -100.0);
    xsref::np_assert_allclose!([actual], [-6.069_130_008_803_468e50], rtol = 1e-14);

    // and here it is subnormal, but still accurate
    let actual = xsf::hyp1f1(200_000.0, -1e-298, -0.001);
    xsref::np_assert_allclose!([actual], [1.464_685_090_427_059_9e298], rtol = 1e-13);
}

#[test]
fn test_hyp1f1_rounded_b_minus_a() {
    // Arb; as above, but for `b > 0`
    let abx = [
        (1.037, 0.037, -25.2),
        (2.3, 0.3, -30.0),
        (1.963, 0.963, -25.2),
        (1.547_000_001, 0.547, -29.0),
    ];
    let expected = [
        -7.732_852_991_965_603e-9,
        1.973_234_762_445_072_6e-10,
        -2.861_749_604_034_37e-10,
        -3.253_371_482_745_043_5e-12,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp1f1_neg_b_large_a() {
    // Arb; DLMF 13.3.1 (forward in a) is unstable here, and these were 6.4e25, 0.99975004419684,
    // NaN, and NaN
    let abx = [
        (10.3, -2e5 - 0.3, 1.0),
        (2.5, -1e7 - 0.5, 1000.0),
        (1e6 + 0.5, -1e7 - 0.3, 30.0),
        (-3e5 - 0.5, -2e5 - 0.3, -1.0),
    ];
    let expected = [
        0.999_948_501_532_098,
        0.999_750_043_755_938_4,
        0.049_787_244_613_952_18,
        0.223_130_522_734_272_4,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
    // this was -1.4e41
    assert_eq!(xsf::hyp1f1(10.3, -1e7 - 0.5, 0.0), 1.0);
    // Arb: -5.4e5011; the terms shrink, but grow again near -b, and this was 0.59
    assert_eq!(xsf::hyp1f1(2.0, -1.2e5 - 0.3, 3.6e4), f64::NEG_INFINITY);
    // as before (here DLMF 13.3.1 is still used)
    assert_eq!(xsf::hyp1f1(2.0, -31.4, 1e5), f64::INFINITY);
}

#[test]
fn test_hyp1f1_neg_b_lossy() {
    // Arb: -165.6, -183.1, -505.3, and -2.1e715; the series cancels, or overflows while its terms
    // still change sign, within 500 terms
    assert!(xsf::hyp1f1(1e5 + 0.3, -0.5, -0.1).is_nan());
    assert!(xsf::hyp1f1(-1e5 - 0.3, -0.5, 0.1).is_nan());
    assert!(xsf::hyp1f1(1e6 + 0.3, -0.5, -0.1).is_nan());
    assert!(xsf::hyp1f1(-1e5 - 0.3, -600.5, 100.0).is_nan());
    // Arb: -7.8e1116 and -6.4e1081; the terms keep their sign after the overflow
    assert_eq!(xsf::hyp1f1(-0.3, -165.5, 2000.0), f64::NEG_INFINITY);
    assert_eq!(xsf::hyp1f1(-1.7, -0.5, 2500.0), f64::NEG_INFINITY);
    // Arb: 2.6e4157; a polynomial with positive terms
    assert_eq!(xsf::hyp1f1(-1e5, -1e6 - 0.5, 1e5), f64::INFINITY);
}

#[test]
fn test_hyp1f1_neg_b_unstable() {
    // Arb; DLMF 13.3.2 (backward in b from b + n > 0) is unstable here
    let abx = [
        (10.3, -99_999.5, 30.0),
        (100.5, -1000.3, 100.0),
        (1.7, -99_999.5, 10.0),
        (31.7, -99_999.5, 300.0),
    ];
    let expected = [
        0.996_915_215_768_620_1,
        7.233_114_758_583_054e-5,
        0.999_830_022_097_624_5,
        0.909_411_072_212_928,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp1f1_neg_b_minimal() {
    // Arb; the series cancels (its terms reach 3e13 for the first one), DLMF 13.3.2 is unstable in
    // both directions, and these were NaN
    let abx = [
        (316.3, -1000.3, 100.0),
        (1188.75, -2499.125, 54.75),
        (24_929.625, -2777.75, 32.625),
        (-21_502.4, -604.6, -1.25),
    ];
    let expected = [
        1.216_742_332_127_376e-13,
        7.373_079_150_971_805e-12,
        9.679_705_093_095_042e-104,
        2.113_210_012_489_118_2e-19,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb: 6.773662582844621e-53, but the continued fraction takes 30573 steps here, and K + I
    // amplifies its rounding by 5e6
    assert!(xsf::hyp1f1(1e6, -1000.3, 0.110_081_953_697_908_98).is_nan());
}

#[test]
fn test_hyp1f1_neg_b_underflow() {
    // Arb; the terms of the series underflow before they grow again near n = -b
    let abx = [(-0.6875, -2491.75, 832.5), (62.75, -3297.25, 960.75), (24.9375, -2301.25, 696.5)];
    let expected =
        [-1.732_239_273_271_290_8e252, 1.526_278_903_660_539_6e224, 2.788_870_728_214_917e168];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-13);
}

#[test]
fn test_hyp1f1_neg_b_underflow_recovered() {
    // Arb; the series cancels, and the terms of the sum that is tried last underflow before they
    // grow again near n = -b too
    let abx = [(16_989.625, -4777.625, 152.75), (-6787.4, -3640.6, -331.375)];
    let expected = [2.438_679_688_403_589e192, -1.353_909_816_060_362_7e-59];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // Arb: -1.4e1096; the terms are negligible after 500 terms, long before they grow again near
    // n = -b
    let (a, b, x) = (0.015_625, -1710.125, 1295.75);
    assert_eq!(xsf::hyp1f1(a, b, x), f64::NEG_INFINITY);
    // Arb: -6.1e2002; the terms underflow to 0 from n = 476 on, long before n = -b
    let (a, b, x) = (-300.375, -12_891.875, 5129.625);
    assert_eq!(xsf::hyp1f1(a, b, x), f64::NEG_INFINITY);
    // Arb: 1.4e3136; the alternating terms before n = -b overflow, even times e^-1400, so there
    // is no reliable sign
    let (a, b, x) = (-2147.625, -5933.125, 6410.875);
    assert!(xsf::hyp1f1(a, b, x).is_nan());
}

#[test]
fn test_hyp1f1_neg_b_poly_beyond_500_terms() {
    // Arb; the polynomials end after 500 terms, long before n = -b
    let abx = [(-500.0, -200_000.5, 100_000.0), (-600.0, -100_000.5, 200_001.0)];
    let expected = [1.190_070_170_014_810_2e88, 4.173_835_946_190_143_6e286];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp1f1_large_b_or_neg_x() {
    // mpmath; in xsf these are 0, 0, 14.511, and -1.3e197
    let abx = [
        (1.0, 1000.0, -800.0),
        (2.0, 3000.0, -3000.0),
        (0.5, 20_000.0, 20_000.0),
        (874.25, 0.25, -138.5),
    ];
    let expected = [
        0.555_445_795_880_786,
        0.249_937_491_318_142_4,
        14.512_932_899_154_44,
        6.546_132_637_255_607e-30,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-12);
}

#[test]
fn test_hyp1f1_tiny_a_or_b() {
    // Arb (scipy/xsf#454)
    let abx = [(1.0, 1e-17, 1.0), (0.5, 1e-10, 1.0), (1e-10, 0.5, 10.0), (0.5, 1e-10, -10.0)];
    let expected = [
        2.718_281_828_459_045e17,
        10_892_917_406.618_07,
        1.000_001_309_153_013,
        -978_427_283.172_141_4,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);

    // Arb (scipy/xsf#454)
    let abz = [
        (0.5, 1e-17, c64(1.0, 1.0)),
        (-2.0, 1e-17, c64(1.0, 1.0)),
        (0.5, 1e-10, c64(2.0, 1.0)),
        (1e-10, 0.5, c64(10.0, 1.0)),
    ];
    let expected = [
        c64(-3_278_454_977_900_773.0, 1.504_208_434_788_265e17),
        c64(-2e17, -2.0),
        c64(14_042_459_676.442_57, 52_867_238_789.940_98),
        c64(1.000_000_766_147_509_9, 1.056_275_226_895_348_6e-6),
    ];
    let actual = abz.map(|(a, b, z)| xsf::hyp1f1(a, b, z));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);
}

#[test]
fn test_hyp1f1_tiny_b_overflow() {
    // Arb; in xsf these are inf, inf, NaN, and NaN, since the terms (~ 1/b) or ln Gamma(b) overflow
    let abx = [
        (-1.0, 1e-300, -20.0),
        (-1.5, 1e-300, -20.0),
        (20.0, 1e-300, -30.0),
        (-50.0, 1e-300, 30.0),
    ];
    let expected = [
        1.999_999_999_999_999_9e301,
        1.046_860_343_479_010_6e302,
        -7.920_209_383_780_516e293,
        7.714_620_828_937_928e306,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 2e-15);

    // Arb; in xsf these are NaN, inf - inf i, and NaN
    let abz = [
        (1e-300, -1e-300, c64(20.0, 1.0)),
        (-10.5, 1e-300, c64(-5.0, 1.0)),
        (0.5, 1e-300, c64(-50.0, 1.0)),
    ];
    let expected = [
        c64(-262_135_871.806_876_27, -408_252_434.775_991_8),
        c64(7.907_968_079_087_868e304, -1.454_232_915_254_029_6e305),
        c64(-4.051_003_078_579_772e298, -4.180_074_037_583_19e296),
    ];
    let actual = abz.map(|(a, b, z)| xsf::hyp1f1(a, b, z));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-15);

    // Arb: 3.6808580837139879e-19; the series with the exact `b - a`, not rounded to -8
    let actual = xsf::hyp1f1(8.0, 1e-300, -780.0);
    xsref::np_assert_allclose!([actual], [3.680_858_083_713_988e-19], rtol = 1e-15);
    // NaN instead of the series with `b - a` rounded to -50, or with `|b| M` underflowing
    assert!(xsf::hyp1f1(50.0, 1e-300, c64(-1000.0, 0.0)).is_nan());
}

#[test]
fn test_hyp1f1_tiny_a_early_stop() {
    // Arb; in xsf the relative errors are 2.0e-9 and 2.3e-6
    let abx = [(1e-17, 0.5, 20.0), (-1e-17, 1.0, 29.5)];
    let expected = [1.000_000_001_975_100_3, 0.999_997_722_693_935_5];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);

    // Arb; in xsf the relative errors are 3.6e-11, 8.9e-9, and 2.2e-2; the last doesn't reach the
    // stopping criterion within 500 terms
    let abz = [
        (1e-17, 1.0, c64(18.0, 7.0)),
        (-1e-17, -2.5, c64(14.0, 6.0)),
        (1e-16, -7.5, c64(19.5, 0.0)),
        (1e-17, -300.25, c64(0.0, 110.0)),
    ];
    let expected = [
        c64(1.000_000_000_033_934, 1.138_606_472_512_429_6e-11),
        c64(1.000_000_006_236_292, 6.308_164_997_372_973_4e-9),
        c64(1.022_280_862_276_224_3, 0.0),
        c64(1.0, -3.409_453_691_196_780_6e-18),
    ];
    let actual = abz.map(|(a, b, z)| xsf::hyp1f1(a, b, z));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);

    // Arb; here the later terms (up to 1e9) cancel, and are below their rounding error
    let actual = xsf::hyp1f1(1e-17, -171.5, c64(12.0, 63.0));
    let expected = c64(1.000_000_000_000_773_4, -3.516_259_387_310_227e-13);
    xsref::np_assert_allclose!([actual], [expected], rtol = 1e-12);
}

#[test]
fn test_bernoulli_small() {
    // B_4, B_6, B_8 = -1/30, 1/42, -1/30
    let b = xsf::bernoulli::<9>();
    let expected = [-1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0];
    xsref::np_assert_allclose!([b[4], b[6], b[8]], expected, rtol = 1e-15);
}

#[test]
fn test_bernoulli_table_and_overflow() {
    let b = xsf::bernoulli::<263>();
    assert_eq!(b[258], 1.335_278_418_735_463_4e306);
    assert_eq!(b[256], -7.950_212_504_588_525e302);
    assert_eq!(b[259], 0.0);
    assert_eq!(b[260], f64::NEG_INFINITY);
    assert_eq!(b[261], 0.0);
    assert_eq!(b[262], f64::INFINITY);
}
