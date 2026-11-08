use num_complex::c64;
use num_traits::ToPrimitive;

#[test]
fn test_hypu() {
    // the table is called "hyperu" instead of "hypu"
    xsref::test("hyperu", "d_d_d-d", |x| xsf::hypu(x[0], x[1], x[2]));
}

#[test]
fn test_hypu_tiny_a() {
    // mpmath (scipy/xsf#350)
    let abx = [
        (1e-17, 0.5, 100.0),
        (-1e-17, 0.5, 100.0),
        (-1.0 + f64::EPSILON / 2.0, 0.5, 200.0),
    ];
    let expected = [1.0, 1.0, 199.499_999_999_999_9];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 5e-16);

    // Arb (scipy/xsf#436)
    let abx = [
        (1e-20, 0.5, 50.0),
        (1e-15, 0.5, 50.0),
        (1e-20, 0.5, 40.0),
        (-1e-20, -0.5, 30.0),
    ];
    let expected = [1.0, 0.999_999_999_999_996_1, 1.0, 1.0];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // Arb; 1 in xsf (`chgul` rounds `a + k - 1`)
    let abx = [
        (1e-20, 10.0, 0.01),
        (-1e-20, 15.0, 0.1),
        (1e-30, 30.0, 1.0),
        (1e-40, 20.0, 0.002),
    ];
    let expected = [
        408.762_021_103_157_6,
        -6_934.305_358_080_271,
        1.859_511_646_545_143_5,
        1.223_736_600_747_381e27,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
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
    let expected = [
        9.999_999_97e17,
        2.704_487_856_590_993_5e-5,
        -2_147_483_646.0,
        1.0,
    ];
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
    let abx = [
        (
            -9.596_157_126_249_977,
            8.220_449_743_267_231,
            29.368_858_516_382_556,
        ),
        (-0.5, 0.0, 2.0),
        (
            -95.934_412_853_269_4,
            0.250_167_345_651_739_6,
            322.254_120_856_906_47,
        ),
        (-150.5, 0.3, 50.0),
    ];
    let expected = [
        25_688_997_285.068_825,
        1.568_794_902_417_112,
        -6.715_069_973_718_165e217,
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

    // still NaN (mpmath: -inf), instead of a wrong 1 from `chgul`
    assert!(xsf::hypu(-4.2e-66, 26.000_000_000_000_018, 1e-300).is_nan());
}

#[test]
fn test_hypu_large_neg_a_large_b() {
    // mpmath; NaN in xsf
    let abx = [
        (
            -70.143_758_091_504_14,
            45.332_941_134_973_275,
            11.429_524_306_986_952,
        ),
        (-60.2, 30.6, 10.0),
        (-150.3, 60.1, 30.0),
        (-2.4, 98.3, 29.8),
    ];
    let expected = [
        -9.122_398_426_133_06e121,
        -1.875_196_314_005_285_6e96,
        4.672_086_732_836_285e291,
        -3.999_673_741_113_012e20,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-8);
}

#[test]
fn test_hypu_poly() {
    // mpmath; -6.90e95, -4.94e187, and NaN in xsf
    let abx = [
        (-60.0, 30.0, 10.0),
        (-100.0, -5.014_823_015_669_013, 46.614_071_395_839_73),
        (-95.0, 15.0, 0.004_221_448_686_339_405_5),
    ];
    let expected = [
        -6.398_262_182_829_726e95,
        -8.403_084_016_177_146e165,
        -1.612_479_198_337_786_4e165,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // mpmath; 3.27e96 and 5.35e190 in xsf
    let actual = [xsf::hypu(-60.5, 30.5, 10.0), xsf::hypu(-100.5, 40.5, 20.0)];
    let expected = [8.123_216_159_896_29e95, -4.104_645_857_033_403e177];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-8);
}

#[test]
fn test_hypu_small_a() {
    // mpmath; 0.045763206, 0.067619064, and 7.5410439102e-4 in xsf
    let abx = [
        (1.1, 0.5, 15.0),
        (
            1.088_244_511_969_767_5,
            -2.443_990_400_454_129_7,
            7.704_068_190_231_311,
        ),
        (2.5, 0.5, 15.0),
    ];
    let expected = [
        0.045_763_173_934_103_04,
        0.067_618_993_534_387_55,
        0.000_754_104_391_072_057_4,
    ];
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
    let expected = [
        0.109_883_867_556_135_05,
        0.251_853_637_950_820_1,
        0.937_895_047_309_956_8,
    ];
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
    let actual = [
        xsf::hypu(-28.0, 97.6, 0.0055),
        xsf::hypu(-74.2, -71.2, 0.0017),
    ];
    let expected = [1.766_091_302_349_042e57, 5.933_445_224_488_155e-197];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);

    // mpmath; 502.06005 and 902602731 (near the poles of Gamma(b)) in xsf
    let actual = [
        xsf::hypu(
            -2.169_910_971_917_65,
            -1.997_705_917_101_09,
            16.709_077_020_224_107,
        ),
        xsf::hypu(-7.5, -6.999_996, 15.0),
    ];
    let expected = [502.060_132_017_751_15, 902_731_434.181_382_9];
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-10);

    // Arb; -8.59e8 and 1.1254 (near the poles of Gamma(b)) in xsf
    let actual = [
        xsf::hypu(3.0, 1e-12, 1.5),
        xsf::hypu(-1.9, -0.999_999_999_6, 0.9),
    ];
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
}

#[test]
fn test_hypu_large_a() {
    // mpmath; NaN in xsf
    let actual = xsf::hypu(
        79.190_551_504_429_4,
        8.120_259_935_805_485,
        3.282_468_195_642_108,
    );
    xsref::np_assert_allclose!([actual], [2.923_069_295_358_514e-124], rtol = 1e-13);
}

#[test]
fn test_hypu_large_b() {
    // mpmath; 2.86e-8 and -666.1 in xsf
    let abx = [
        (
            3.688_232_632_414_102_6,
            -87.588_828_219_563_88,
            20.526_008_089_799_742,
        ),
        (
            -1.280_941_997_975_485_2,
            70.676_192_639_949_22,
            126.757_411_617_707_94,
        ),
    ];
    let expected = [2.870_107_296_557_804_4e-8, 172.021_142_044_546_17];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-11);
}

#[test]
fn test_hypu_large_neg_b() {
    // mpmath; 25.65 and 2.638e6 in xsf
    let abx = [
        (0.2, -138.4, 40.4),
        (
            -0.431_695_526_174_340_5,
            -141.653_296_128_945_6,
            43.936_589_915_444_64,
        ),
    ];
    let expected = [0.354_135_821_455_708_52, 9.542_883_153_865_268];
    let actual = abx.map(|(a, b, x)| xsf::hypu(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
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
        (
            -0.001_455_282_650_156_073,
            -170.624_101_975_403_85,
            1.490_954_765_747_425_4,
        ),
    ];
    let expected = [
        5040.0,
        5.426_195_543_094_943,
        0.099_621_961_951_682_32,
        6.809_106_301_661_743e220,
        17.385_389_278_273_84,
        0.231_064_916_158_941_04,
        1.007_524_453_974_522_3,
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
    let abx = [
        (13.5, 25.0, 0.005),
        (1.0, -26.0, 0.02),
        (25.0, 2.0, 5.0),
        (30.0, -10.0, 2.0),
    ];
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

    // Arb; the over- and underflows aren't NaN
    assert_eq!(xsf::hypu(171.61, 1.0, 2.0), 0.0);
    assert_eq!(xsf::hypu(-3.0 + 1e-15, 999.0, 1e-300), f64::NEG_INFINITY);
    assert_eq!(xsf::hypu(-1e10, 10.0, 1e-300), f64::INFINITY);
    assert_eq!(xsf::hypu(-169.0 - 3e-14, 14.0, 4.0), f64::INFINITY);
    assert_eq!(xsf::hypu(-167.0 - 3e-14, -16.0, 385.0), f64::NEG_INFINITY);

    // Arb; still from `chgubi` (11 digits) if `chguit` overflows
    let actual = xsf::hypu(150.0, 9.0, 0.1);
    xsref::np_assert_allclose!([actual], [2.229_672_231_677_937e-250], rtol = 1e-12);
}

#[test]
fn test_hyp1f1() {
    xsref::test("hyp1f1", "d_d_cd-cd", |x| {
        xsf::hyp1f1(x[0], x[1], c64(x[2], x[3]))
    });
}

#[test]
fn test_pmv() {
    xsref::test("pmv", "d_d_d-d", |x| {
        xsf::pmv(x[0].to_i32().unwrap(), x[1], x[2])
    });
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
    let abx = [
        (0.5, 1e9, -1e8),
        (100.3, 2e8, -2e6),
        (3.7, 2.1e9, -1e9),
        (-3.7, 1e7, -5e6),
    ];
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
    // these took ~1e9 steps
    assert!(xsf::hyp1f1(-1.2e9, -10.2, -1.0).is_nan());
    assert!(xsf::hyp1f1(1e9, 2.5, f64::NAN).is_nan());
}

#[test]
fn test_hyp1f1_neg_x_no_convergence() {
    // Arb; the series after DLMF 13.2.39 doesn't converge within 500 terms (b + 500 < 0), and these
    // were e^-x times too large
    let abx = [
        (2e7 + 0.3, -2e5 - 0.3, -1.0),
        (10.3, -1e6 - 0.5, -100.0),
        (-1000.3, -2e5 - 0.3, -100.0),
    ];
    let expected = [
        2.756_493_631_860_608e43,
        1.001_030_581_673_679_5,
        0.606_364_673_033_219_6,
    ];
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
    let abz = [
        (0.5, 1.5, c64(5.0, -25.0)),
        (1.25, 2.5, c64(-11.0, 20.0)),
        (0.25, 3.5, c64(-2.0, 30.0)),
    ];
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
    let abx = [
        (500.5, -200.25, 5.5),
        (99.7, -0.3, -800.0),
        (340.5, -90.25, 140.5),
    ];
    let expected = [
        1.952_103_586_824_771_7e-6,
        9.487_619_910_475_397e-140,
        -6.941_454_309_937_287e307,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-13);
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
    assert!(xsf::hyp1f1(2.0, -1.2e5 - 0.3, 3.6e4).is_nan());
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
    // Arb; in xsf the first is inf, and the others have a relative error of 8.3e-8, 1.1e-13, and
    // 8.3e-8
    let abx = [
        (1.0, 1e-17, 1.0),
        (0.5, 1e-10, 1.0),
        (1e-10, 0.5, 10.0),
        (0.5, 1e-10, -10.0),
    ];
    let expected = [
        2.718_281_828_459_045e17,
        10_892_917_406.618_07,
        1.000_001_309_153_013,
        -978_427_283.172_141_4,
    ];
    let actual = abx.map(|(a, b, x)| xsf::hyp1f1(a, b, x));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-15);

    // Arb; in xsf the first two are NaN, and the others have a relative error of 8.3e-8 and 1.1e-13
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
