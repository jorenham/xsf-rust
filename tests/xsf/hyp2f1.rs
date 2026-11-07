use num_complex::c64;
use xsref::np_assert_allclose;

#[test]
fn test_hyp2f1_c64_transform2_euler() {
    // unlike in xsf: the 1/z transformation of the Euler transform (mpmath references)
    for (a, b, c, z, expected) in [
        // a - b an integer: the limiting case
        (
            8.5,
            6.5,
            -4.25,
            c64(-1.0, 0.5),
            c64(-14.914_154_962_801_387, 2.508_617_290_210_916_7),
        ),
        (
            4.3,
            4.4,
            -4.2,
            c64(-0.9, 0.8),
            c64(-13.468_832_852_470_712, -27.983_525_735_799_486),
        ),
        // the limiting case with b within 1e-15 of an integer
        (
            8.881_784_197_001_252e-16,
            8.881_784_197_001_252e-16,
            -1.25,
            c64(-1.0, 0.5),
            c64(1.0, -4.661_160_101_737_811e-31),
        ),
        // (1 - z)^(c - a - b) is subnormal
        (
            9.48,
            9.48,
            -7.69,
            c64(2.97e11, 8.92e11),
            c64(3.884_057_170_000_241_6e-108, -6.451_235_141_304_403e-108),
        ),
        // c - (c - a) loses a
        (
            5e-6,
            -3.125,
            -9.0625,
            c64(133.5, -390.75),
            c64(3.466_873_805_327_974, -0.823_490_927_100_644_6),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        assert!(
            (actual - expected).norm() < 1e-13 * expected.norm(),
            "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
        );
    }
}

#[test]
fn test_hyp2f1_c64_transform1_maclaurin() {
    // unlike in xsf: the Maclaurin series if the 1 - z transformation cancels (mpmath references)
    for (a, b, c, z, expected) in [
        // c - a - b = 50 (the limiting case) and 50.5
        (
            4.5,
            2.5,
            57.0,
            c64(0.64, 0.65),
            c64(1.122_209_544_812_728_6, 0.158_257_844_232_734_11),
        ),
        (
            4.5,
            2.5,
            57.5,
            c64(0.64, 0.65),
            c64(1.121_218_275_094_999_2, 0.156_607_234_073_543_14),
        ),
        // |z| > 1
        (
            3.5,
            2.5,
            40.0,
            c64(0.9, 0.6),
            c64(1.211_419_901_187_38, 0.188_557_093_819_745_4),
        ),
        // c - a - b just misses an integer, so the 1 - z transformation is NaN
        (
            0.1,
            0.2,
            50.3,
            c64(0.64, 0.65),
            c64(1.000_254_319_428_871_4, 2.627_839_705_956_765e-4),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        assert!(
            (actual - expected).norm() < 1e-14 * expected.norm(),
            "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
        );
    }
}

#[test]
fn test_hyp2f1_c64_huge_num_terms() {
    // unlike in xsf: these summed up to 2^64 terms (mpmath references)
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        for (a, b, c, z, expected) in [
            (4.39, 2.14, 1e17, c64(0.64, 0.65), c64(1.0, 6.106_49e-17)),
            (4.39, 2.14, 1e300, c64(0.64, 0.65), c64(1.0, 6.106_49e-300)),
            (
                0.5,
                0.25,
                1e6 + 0.75,
                c64(0.64, 0.65),
                c64(1.000_000_079_999_938_5, 8.125_003_656_241_712e-8),
            ),
            // the most terms that are still summed: 2^20
            (
                -1_048_575.0,
                1.0,
                2.0,
                c64(1e-7, 0.0),
                c64(0.949_356_716_098_484_6, 0.0),
            ),
        ] {
            let actual = xsf::hyp2f1(a, b, c, z);
            assert!(
                (actual - expected).norm() < 1e-15 * expected.norm(),
                "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
            );
        }
        // these exceed the cap on the number of terms
        for (a, b, c, z) in [
            (-1_048_576.0, 1e-30, 1.0, c64(-1e-4, 0.0)),
            (-1e15, 1.0, 2.0, c64(0.5, 0.0)),
            (1e17, 2.0, 4.39, c64(3.0, 0.65)),
            (1e9, 2.0, 4.39, c64(3.0, 0.65)),
            (2_251_799_813_685_248.5, 0.25, 0.75, c64(0.95, 0.1)),
        ] {
            assert!(xsf::hyp2f1(a, b, c, z).re.is_nan());
        }
        tx.send(()).unwrap();
    });
    rx.recv_timeout(std::time::Duration::from_secs(10)).unwrap();
}

#[test]
fn test_hyp2f1_f64_c_minus_a_minus_b_le_minus_1() {
    // the Euler transformation cancels here (mpmath)
    let actual = [
        xsf::hyp2f1(-0.5, 125.0, 62.5, 0.4),
        xsf::hyp2f1(-0.8, 19.7, 9.4, 0.9),
    ];
    let desired = [0.453_216_679_866_897_07, -487.198_231_355_781_05];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but the recurrence of the power series is inaccurate here (mpmath)
    let actual = [
        xsf::hyp2f1(-1.000_000_000_000_11, 2.000_002, 1e-6, 0.999_999_7),
        xsf::hyp2f1(-49.5, 61.5, 6.5, 0.961_281_567_773_062_2),
    ];
    let desired = [-2_000_000.033_618_169_8, -6.640_002_805_058_352e-5];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // and here the power series without that recurrence is more accurate than both (mpmath)
    let actual = [
        xsf::hyp2f1(-3.1, 17.4, 12.3, 0.97),
        xsf::hyp2f1(18.2, -0.9, 12.3, 0.98),
        xsf::hyp2f1(-6.2, 7.7, -2.5, 0.99),
        xsf::hyp2f1(-0.0076, 19.42, 9.7355, 0.8761),
    ];
    let desired = [
        -0.010_434_102_923_645_322,
        -141.307_919_251_094_5,
        45_908_976.313_493_14,
        -30.905_387_588_773_838,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // (but not past the special cases of the Euler transformation, where c - (c - a) rounds to
    // -5 or -12) (mpmath)
    let actual = [
        xsf::hyp2f1(-5.000_000_000_000_101, 106.100_000_000_000_2, 100.1, 0.95),
        xsf::hyp2f1(-11.999_999_999_999_899, 122.100_000_000_000_2, 100.1, 0.9),
    ];
    let desired = [2.580_226_780_133_956_7e-8, 3.231_061_401_736_383e-14];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // (nor if the series doesn't converge: this is -6.33e35, not its partial sum 0.856)
    let actual = xsf::hyp2f1(-1e-13, 3.100_000_000_000_1, 0.1, 0.999_999_999_999_999_9);
    np_assert_allclose!([actual], [-6.326_846_914_851_421e35], rtol = 1e-4);
}

#[test]
fn test_hyp2f1_f64_recurrence_in_a() {
    // the recurrence in `a` can be unstable, while its error estimate `loss` is small (mpmath)
    let actual = [
        xsf::hyp2f1(29.0, -14.5, 0.1, 0.99),
        xsf::hyp2f1(14.0, -14.5, 0.7, 0.98),
        xsf::hyp2f1(-15.8, 25.0, 0.5, 0.98),
    ];
    let desired = [
        -5.131_087_011_432_369e21,
        -0.093_683_512_477_367_49,
        268_344_349.251_342_5,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_x_lt_minus_2() {
    // the 1/x transformation is inaccurate here, but the Pfaff transformations are not (mpmath)
    let actual = [
        xsf::hyp2f1(12.0, 3.1, 99.3, -2.8),
        xsf::hyp2f1(81.0, 6.2, -95.1, -8.4),
        xsf::hyp2f1(-58.5, 0.07, 50.9, -4.3),
    ];
    let desired = [
        0.416_874_471_029_815_73,
        6.875_723_600_176_229e-5,
        8_488_942_027_952.015,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but here the starting values of their recurrence need a recurrence as well (mpmath)
    let actual = xsf::hyp2f1(-20.000_000_000_2, -14.47, 2.93, -5.48);
    np_assert_allclose!([actual], [64_520_380_601.436_584], rtol = 1e-12);

    // or need the series without that recurrence (mpmath)
    let actual = [
        xsf::hyp2f1(-19.6, -4.1, 16.7, -9.9),
        xsf::hyp2f1(-4.2, -19.6, 14.6, -9.6),
    ];
    let desired = [-1_046_594.825_411_287_6, -4_908_149.039_599_091];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_c_lt_0_x_lt_0() {
    // the Euler transformation cancels here, but the Pfaff ones don't (mpmath)
    let actual = [
        xsf::hyp2f1(6.5, 8.9, -9.8, -0.08),
        xsf::hyp2f1(9.9, 4.5, -9.9, -0.01),
        xsf::hyp2f1(16.5, -19.2, -19.7, -0.3),
        xsf::hyp2f1(14.5, -15.5, -19.6, -0.63),
    ];
    let desired = [
        1.707_641_797_398_377,
        1.046_566_861_351_764,
        1.874_594_164_146_283_7e-4,
        -0.143_812_194_757_134_78,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_c_lt_0() {
    // the terms of the power series are tiny before k = -c, and huge after it (mpmath)
    let actual = [
        xsf::hyp2f1(0.01, 0.3, -181.2, 0.6),
        xsf::hyp2f1(0.1, 0.2, -1000.5, 0.55),
        xsf::hyp2f1(0.1, 0.2, -12_000.5, 0.02),
        xsf::hyp2f1(-7001.5, 12.2, -7000.5, 0.5),
    ];
    let desired = [
        3.196_698_793_921_328_4e28,
        -6.233_446_442_910_298e83,
        0.999_999_966_668_092_2,
        4_713.283_633_348_583,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // the ratio of its terms underflows here, but not for the backward recurrence in c (the direct
    // sum of the power series in high precision)
    let actual = [
        xsf::hyp2f1(1e-200, 1e-200, -1000.5, 0.8),
        xsf::hyp2f1(1e-160, 1e-160, -1000.5, 0.8),
    ];
    let desired = [-5.767_096_442_064_845e199, -5.767_096_442_064_845e279];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_parameters() {
    // unlike in xsf: the backward recurrence in c for large parameters (references: the direct sum
    // of the power series of the Pfaff transformation in high precision)
    let actual = [
        xsf::hyp2f1(35.17, 102.3, 0.0054, -9.74),
        xsf::hyp2f1(658.99, 0.0047, 51.38, -2.72),
        xsf::hyp2f1(-0.5, 692.8, 83.35, -3.98),
        xsf::hyp2f1(-21.2, -10.2, -1e-8, -0.46),
    ];
    let desired = [
        1.955_342_940_723_380_9e-60,
        0.983_273_621_397_356_7,
        5.845_197_863_408_547,
        -53_404_320_641.112_48,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // x << -1: the 1/x transformation (references: its two series summed in high precision)
    let actual = [
        xsf::hyp2f1(1.35, 87.4, -32.1, -8478.0),
        xsf::hyp2f1(-41.6, -16.3, 720.1, -204.7),
        xsf::hyp2f1(387.2, -2.01, 1.7, -6264.9),
        xsf::hyp2f1(849.1, -1.54, 55.0, -236.7),
    ];
    let desired = [
        -4.397_527_737_978_469e-6,
        -1.151_698_277_128_985_8e17,
        1_471_625_105_699.610_7,
        304_675.781_003_058_8,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // long recurrences, whose solutions oscillate in c (references: Arb), e.g. 1.6e617 overflows
    let actual = [
        xsf::hyp2f1(-522.3, -123.4, 533.9, -0.9),
        xsf::hyp2f1(-552.3, -953.1, -55.5, -4.13),
        xsf::hyp2f1(-679.9, -678.6, -28.7, -4.79),
    ];
    let desired = [5.957_219_702_847_165e-52, f64::INFINITY, f64::NEG_INFINITY];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_opposite_signs() {
    // unlike in xsf: the 1 - x transformation if 0 < x < 1, and the 1/(1 - x) transformation if
    // x < 0 (with a and b of the same sign then), where the Pfaff and Euler forms cancel
    // (references: Arb)
    let actual = [
        xsf::hyp2f1(-975.2, 241.4, 0.00135, 0.39),
        xsf::hyp2f1(-880.55, -75.65, 0.091, -0.472),
    ];
    let desired = [-1.469_500_580_590_125_3e-75, 8.556_925_324_233_141e83];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_c_lt_0() {
    // unlike in xsf: the backward recurrence in c from start values by other forms, if those of
    // its own form lose precision (references: Arb)
    let actual = [
        xsf::hyp2f1(758.1, 7.7, -151.2, -0.3),
        xsf::hyp2f1(16.1, 428.5, -230.7, -0.91),
        xsf::hyp2f1(-559.7, 42.4, -136.1, 0.27),
    ];
    let desired = [
        -8_553.508_491_012_49,
        -1_050_587.466_964_691_6,
        -2.046_272_877_671_090_4e24,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // from the start values at c + 1 and c + 2, which only help since the other forms nearly
    // succeeded (reference: Arb)
    let actual = xsf::hyp2f1(130.25, 125.25, 1.96, -145.46);
    np_assert_allclose!([actual], [-3.043_786_288_636_934_7e-280], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_near_1() {
    // unlike in xsf: the 1 - x transformation of `hyp2f1_recur` for x >= 0.9, e.g. these were
    // -1.83e6, inf and NaN (references: Arb)
    let actual = [
        xsf::hyp2f1(40.1, -133.3, -38.6, 0.901),
        xsf::hyp2f1(-18.5, 75.0, 172.4, 0.901),
        xsf::hyp2f1(-2.9, -667.5, 1.8, 1.0 - 1e-14),
    ];
    let desired = [
        43.700_421_521_573_67,
        1.614_794_417_233_753e-4,
        9_480_012.959_767_725,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_integer_b_minus_a() {
    // unlike in xsf: the log case of the 1/x transformation of `hyp2f1_recur` if b - a is an
    // integer, e.g. the first was NaN (references: Arb)
    let actual = [
        xsf::hyp2f1(12.5, 812.5, 300.25, -25_000.0),
        xsf::hyp2f1(25.0, 25.0, -0.5, -1e6),
        xsf::hyp2f1(1.5, -50.5, 2.25, -100.0),
        xsf::hyp2f1(-40.5, -20.5, 2.25, -100.0),
        xsf::hyp2f1(25.5, 26.5, 27.25, -100.0),
    ];
    let desired = [
        3.490_725_745_476_049e-61,
        -8.873_628_661_650_908e-148,
        1.101_295_428_271_828_5e100,
        -2.154_549_921_333_831e66,
        7.125_787_400_134_933e-51,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_polynomials() {
    // unlike in xsf: the recurrence in the degree, since the series cancel (references: Arb)
    let actual = [
        xsf::hyp2f1(-522.0, -123.5, 533.5, -0.9),
        xsf::hyp2f1(-522.0, -123.5, 533.5, -0.5),
        xsf::hyp2f1(-3.75, 10_675.75, -0.25, -5600.0),
    ];
    let desired = [
        6.485_127_939_853_472_5e-52,
        -4.200_865_766_102_724e-45,
        -2.143_761_498_506_179_6e29,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}
