use num_complex::c64;
use num_traits::ToPrimitive;

#[test]
fn test_hypu() {
    // the table is called "hyperu" instead of "hypu"
    xsref::test("hyperu", "d_d_d-d", |x| xsf::hypu(x[0], x[1], x[2]));
}

#[test]
fn test_hypu_tiny_a() {
    // NaN in xsf, because of an undefined `int` conversion
    assert_eq!(xsf::hypu(1e-17, 0.5, 100.0), 1.0);
}

#[test]
fn test_hypu_huge_a() {
    // NaN without ~2^31 iterations, like in xsf
    assert!(xsf::hypu(-2_147_483_649.0, 1.0, 1.0).is_nan());
    assert!(xsf::hypu(1.0, f64::INFINITY, 1.0).is_nan());
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

    // still NaN (mpmath: 1.0075 and -inf), instead of a wrong 0 if Gamma(a) overflows in
    // `chguit`, or a wrong 1 from `chgul`
    assert!(
        xsf::hypu(
            -0.001_455_282_650_156_073,
            -170.624_101_975_403_85,
            1.490_954_765_747_425_4
        )
        .is_nan()
    );
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
    // mpmath; in xsf these have a relative error of 0.13, 2.4e-5, and 4.9e-6 (and ~1e-8 here,
    // because the asymptotic expansion is truncated)
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
