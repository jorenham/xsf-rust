use num_complex::c64;
use xsref::np_assert_allclose;

xsref::tests! {
    test_hyp2f1_f64: "hyp2f1", "d_d_d_d-d", |x| xsf::hyp2f1(x[0], x[1], x[2], x[3]);
    test_hyp2f1_c64: "hyp2f1", "d_d_d_cd-cd", |x| xsf::hyp2f1(x[0], x[1], x[2], c64(x[3], x[4]));
}

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

    // but its error estimate shouldn't be too large if a starting value is tiny, so that the less
    // accurate alternatives aren't used instead (Arb)
    let actual = [
        xsf::hyp2f1(-9.5, 10.5, 0.5, 0.95),
        xsf::hyp2f1(10.000_000_001, 9.999_999_999_999_99, 1.5, -3.0),
        xsf::hyp2f1(9.999_999_999, 20.0, -10.5, -3.0),
    ];
    let desired = [
        -0.897_769_973_095_626_9,
        4.289_686_233_768_906e-8,
        0.001_570_408_561_128_409_7,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // and its starting values, by the power series with c = 1 - (c - a - b) = -8 - 1e-14 in AMS55
    // #15.3.6, shouldn't stop before the near pole of c (Arb)
    let actual = [
        xsf::hyp2f1(-12.5, 8.999_999_999_999_99, 5.5, 0.998),
        xsf::hyp2f1(-11.5, 8.999_999_999_999_99, 6.5, 0.998),
    ];
    let desired = [2.259_678_955_985_816e-8, -3.542_896_384_613_633e-8];
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
fn test_hyp2f1_f64_x_lt_minus_2_b_minus_a_near_int() {
    // unlike in xsf: these were -1.26e5, -2.06e69 and -5.53e3, the last of which is the second
    // function of the 1/x transformation of the one before it (Arb)
    let actual = [
        xsf::hyp2f1(-13.2, -29.200_000_000_001, 5.5, -2.5),
        xsf::hyp2f1(-13.2, -33.21, 50.7, -1000.0),
        xsf::hyp2f1(-33.21, -82.91, -19.01, -0.001),
    ];
    let desired = [
        -602_275.348_686_216_6,
        4.298_609_352_869_04e65,
        1.155_996_198_053_486_3,
    ];
    np_assert_allclose!(actual, desired, rtol = 5e-13);

    // and this was -2.74e6 (Arb)
    let actual = xsf::hyp2f1(-4.6, -20.600_000_000_1, 50.7, -10.0);
    np_assert_allclose!([actual], [-65.931_816_188_886_59], rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_x_lt_minus_2_pfaff_overflow() {
    // unlike in xsf: these were NaN for a < b, as one of the Pfaff series overflows there, and
    // the 1/x transformation is NaN (Arb)
    let args = [
        (1e-20, 3000.5, 1e7, -100.0),
        (3.1e-111, 344_352.06, 1.77e10, -621.8),
        (0.1, 3000.5, 3600.6, -100.0),
        (1.5, 344_352.06, 413_222.472, -3.0),
    ];
    let desired = [1.0, 1.0, 0.641_802_447_619_336_4, 0.152_720_780_375_729_64];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        // also with a and b swapped, which was already accurate
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-14);
    }
}

#[test]
fn test_hyp2f1_f64_tiny_first_term() {
    // unlike in xsf: these were 1, as the series stopped after a tiny first term (Arb)
    let actual = [
        xsf::hyp2f1(1e-20, -300.5, 300.5, -10.0),
        xsf::hyp2f1(1e-300, -500.5, 0.7, -10.0),
        xsf::hyp2f1(1e-20, 150.5, 100.7, 0.9),
        xsf::hyp2f1(-1e-20, 50.5, 0.7, 0.8),
    ];
    let desired = [
        6.528_217_172_881_034e123,
        2.950_367_578_022_261_5e219,
        1.000_000_046_617_229,
        -63_717_533_596_884.64,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // Arb: 1.88e1141
    assert_eq!(xsf::hyp2f1(1e-300, -3000.5, 3000.5, -10.0), f64::INFINITY);
    // the terms grow again here, but stay negligible
    assert_eq!(
        xsf::hyp2f1(1e-300, 4.000_000_000_000_008e30, 1e30, 0.25),
        1.0
    );
}

#[test]
fn test_hyp2f1_f64_euler_tiny_a() {
    // unlike in xsf: these were 1, 1, inf, 1 and -1.24e43, as `a` or `b` was lost in `c - a` or
    // `c - b` (Arb)
    let actual = [
        xsf::hyp2f1(1e-20, 300.5, 0.7, 0.5),
        xsf::hyp2f1(300.5, 1e-20, 0.7, 0.5),
        xsf::hyp2f1(1e-20, -800.5, 0.1, -1.5),
        xsf::hyp2f1(1e-20, 40.5, 20.5, 0.99),
        xsf::hyp2f1(1e-20, -370.5, 53.5, -1.0),
    ];
    let desired = [
        6.921_479_145_123_81e68,
        6.921_479_145_123_81e68,
        1.824_849_365_063_929_7e299,
        63_529_516.407_577_99,
        2.948_782_385_828_795e37,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // and these are 1 to rounding
    let actual = [
        xsf::hyp2f1(1e-20, 3.0, 0.7, -1.0),
        xsf::hyp2f1(1e-20, 300.5, 0.7, -0.999_999_999_999_99),
    ];
    np_assert_allclose!(actual, [1.0, 1.0], rtol = 1e-14);
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
fn test_hyp2f1_f64_a_eq_minus_b_neg_int() {
    // Arb
    let args = [
        (13.0, -13.0, 0.5, 0.9),
        (16.0, -16.0, 0.5, 0.95),
        (19.0, -19.0, 0.5, 0.99),
        (22.0, -22.0, -3.5, 0.95),
        (195.0, -195.0, 1.0, 0.9),
    ];
    let desired = [
        0.489_514_398_515_199_2,
        0.595_230_641_982_666_2,
        0.787_059_088_545_409,
        -6_540_705.712_015_875,
        -0.019_154_323_187_421_253,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        // also with a and b swapped, which was already accurate
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 2e-14);
    }

    // via the recurrence in the degree, accurate to ~1e-13 here
    let actual = [
        xsf::hyp2f1(195.0, -195.0, -200.5, 0.01),
        xsf::hyp2f1(-195.0, 195.0, -200.5, 0.01),
    ];
    np_assert_allclose!(actual, [6.722_964_015_601_578; 2], rtol = 1e-13);

    // and beyond the recurrence limit of the power series, in both orders by `hyp2f1_recur`
    let actual = [
        xsf::hyp2f1(10_001.0, -10_001.0, 9_999.5, 1e-12),
        xsf::hyp2f1(-10_001.0, 10_001.0, 9_999.5, 1e-12),
    ];
    np_assert_allclose!(actual, [0.999_999_989_997_499_8; 2], rtol = 1e-12);
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
    // also for a tiny c, whose `ln Gamma(c)` has no second-order error term (Arb)
    let actual = [xsf::hyp2f1(20.5, 10.25, 1e-200, -1e6)];
    np_assert_allclose!(actual, [-5.514_252_855_663_787e132], rtol = 1e-12);

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
fn test_hyp2f1_f64_huge_c() {
    // unlike in xsf: these were NaN, and the first and fourth hung (Arb, or the limit 1)
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        tx.send([
            xsf::hyp2f1(1.3, 6.5, f64::INFINITY, 0.999),
            xsf::hyp2f1(1.3, 6.5, f64::INFINITY, 1.0),
            xsf::hyp2f1(1.3, 6.5, f64::INFINITY, -3.0),
            xsf::hyp2f1(1.3, 6.5, 1e300, 0.999),
            xsf::hyp2f1(100_000_000.5, 1_000_000_000.5, 1e17, 0.95),
        ])
        .unwrap();
    });
    let actual = rx.recv_timeout(std::time::Duration::from_secs(10)).unwrap();
    let desired = [1.0, 1.0, 1.0, 1.0, 2.585_709_685_660_995_7];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_overflow() {
    // unlike in xsf: these were NaN, as the Gamma functions in the psi function expansion for
    // x > 0.9 and integer c - a - b over- or underflow (Arb)
    let actual = [
        xsf::hyp2f1(1.3, 6.5, 107.8, 0.95),
        xsf::hyp2f1(-1.5, -1.5, 110.0, 0.99),
        xsf::hyp2f1(0.25, 20.5, 283.75, 0.95),
        xsf::hyp2f1(20.5, 30.5, 10_051.0, 0.95),
    ];
    let desired = [
        1.080_572_289_607_784_4,
        1.020_272_592_726_024_6,
        1.017_973_915_497_389_6,
        1.061_033_031_294_106_8,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-15);

    // but a correct overflow (to 3.96e320) is kept
    assert_eq!(xsf::hyp2f1(6.5, 6.5, 1e-300, 0.95), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_ams55_15_3_6_gamma_range() {
    // unlike in xsf: these were 1, NaN, NaN and NaN, as a factor of AMS55 #15.3.6 left the normal
    // range, e.g. `Gamma(-d) / (Gamma(a) Gamma(b))` for a tiny `b`, or `Gamma(c)` (Arb)
    let actual = [
        xsf::hyp2f1(211.205, 4.1e-92, 141.983, 0.9992),
        xsf::hyp2f1(90.5, 1e-180, 1.25, 0.99999),
        xsf::hyp2f1(26.5, -508.5, 590.5, 0.999_999),
        xsf::hyp2f1(-1.5, -693.5, -176.25, 0.999_999_97),
    ];
    let desired = [
        3.803_424_824_536_675e64,
        5.865_541_318_634_591e263,
        5.295_559_383_862_064_7e-8,
        5.041_937_663_630_442,
    ];
    np_assert_allclose!(actual, desired, rtol = 5e-13);

    // but not at a pole of `Gamma(c - a - b)` (this gave -inf), where `hyp2f1_recur` overflows
    // to inf, as the true value does (Arb)
    assert_eq!(xsf::hyp2f1(1e5, 7.5, 21.5, 0.95), f64::INFINITY);
    // and here an inner series fails (this gave -3.96e59), but `hyp2f1_recur` does not (Arb)
    let actual = [xsf::hyp2f1(862.559, -310.305, 906.676, 0.999_993_238)];
    np_assert_allclose!(actual, [2.416_601_765_882_886e-243], rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_near_int() {
    // unlike in xsf: c - a - b (or b - a for x < -1) within 1e-13 of an integer cancelled in
    // AMS55 #15.3.6, which gave 573, -1.03, -1.33e14, 5.23345 and 3.49e-4 (Arb)
    let actual = [
        xsf::hyp2f1(7.5, 2.5 - 1e-14, 10.0, 0.999),
        xsf::hyp2f1(-12.5, 23.0 - 1e-14, 5.5, 0.95),
        xsf::hyp2f1(-2.5, 2.5, 1e-14, 0.99),
        xsf::hyp2f1(1.5, 4.5 - 1e-14, 8.0, 0.999),
        xsf::hyp2f1(3.0, 3.0 - 1e-14, 50.7, -1000.0),
    ];
    let desired = [
        465.329_511_494_182_9,
        0.066_364_450_275_830_67,
        -113_060_171_546_044.02,
        5.233_384_708_085_254,
        7.920_150_150_695_394e-5,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but the power series is tried first (the psi function expansion gives 6.746770e-4)
    let actual = [xsf::hyp2f1(7.5, 10.5 - 1e-13, 50.7, -10.0)];
    np_assert_allclose!(actual, [6.746_767_352_078_354e-4], rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_int() {
    // unlike in xsf: the power series is also tried first if c - a - b is an integer, here after
    // Pfaff, where the psi function expansion gave 6.7467671e-4, 2.0396740e-3, 0.1223477772 and
    // 0.358811271431 (Arb)
    let actual = [
        xsf::hyp2f1(7.5, 10.5, 50.7, -10.0),
        xsf::hyp2f1(7.5, 8.5, 50.7, -10.0),
        xsf::hyp2f1(3.0, 6.0, 50.7, -10.0),
        xsf::hyp2f1(1.5, 5.5, 50.7, -10.0),
    ];
    let desired = [
        6.746_767_352_078e-4,
        2.039_673_789_485_802e-3,
        0.122_347_777_162_347_68,
        0.358_811_271_434_322,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_huge_c_minus_1_lt_c_minus_a_minus_b_lt_0() {
    // unlike in xsf: the recurrence on c overflowed the stack for |c| >= 2^53, or was wrong
    // (6.78e137 for the last) (Arb)
    let actual = [
        xsf::hyp2f1(1e17 + 16.0, -15.5, 1e17, 0.5),
        xsf::hyp2f1(1e17 + 16.0, -15.5, 1e17, 0.95),
        xsf::hyp2f1(1e17 - 96.0, 96.5, 1e17, -0.1),
        xsf::hyp2f1(1e17 - 96.0, 96.5, 1e17, 0.7),
        xsf::hyp2f1(1e16 - 1000.0, 1000.9, 1e16, -0.3),
    ];
    let desired = [
        2.157_918_643_757_769e-5,
        6.823_937_919_615_83e-21,
        1.012_991_688_969_088_8e-4,
        2.869_451_704_379_998e50,
        8.997_019_796_944_917e-115,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    assert!(xsf::hyp2f1(1e17, 0.5, 1e17, f64::NAN).is_nan());
}

#[test]
fn test_hyp2f1_f64_tiny_c_or_c_minus_a() {
    // unlike in xsf: tiny c and c - a were treated as 0, which gave 2, inf, NaN, 1.0 and
    // 1.0000000000000023e-5 (Arb)
    let actual = [
        xsf::hyp2f1(-1e-30, 1.0, 1e-17, 0.5),
        xsf::hyp2f1(5.2, 12.5, -1e-17, 0.6),
        xsf::hyp2f1(12.5, 0.3, 1e-15, 0.2),
        xsf::hyp2f1(1e-14, -1.5, 1e-17, -1.5),
        xsf::hyp2f1(-0.4, -2.5, -0.400_000_000_000_09, 0.99),
    ];
    let desired = [
        0.999_999_999_999_9,
        -5.650_224_541_964_126e27,
        6.795_217_861_653_541e15,
        2_953.847_075_210_480_4,
        1.000_000_042_269_799_7e-5,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_euler_pow_underflow() {
    // unlike in xsf: (1 - x)^(c - a - b) under- or overflowed in the Euler transformation, which
    // gave 0, 0, 0, 0 and inf (Arb)
    let actual = [
        xsf::hyp2f1(9.7, 23.71, -20.7, -1e6),
        xsf::hyp2f1(0.5, 3.7, 0.7, -1e100),
        xsf::hyp2f1(0.5, 9.7, -5.5, -1e30),
        xsf::hyp2f1(2.5, 50.7, 20.37, -1e10),
        xsf::hyp2f1(3.0, 20.37, 4.000_000_000_000_001, 0.999_999_999_999_999_9),
    ];
    let desired = [
        2.339_902_740_854_614e-57,
        1.643_323_996_265_172_4e-51,
        2.444_777_663_146_785_3e-287,
        8.912_118_823_854_448e-27,
        1.699_790_832_401_936_6e308,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_euler_both_out_of_range() {
    // Arb; both factors of the Euler transformation leave the normal range here
    let actual = [
        xsf::hyp2f1(14.5, 1e-290, -368.7, -135.0),
        xsf::hyp2f1(97.0, 1e-240, -16.1, -7.7e5),
        xsf::hyp2f1(4.1, -6.5e-8, -68.9, -1.25e7),
        xsf::hyp2f1(-1.5e-9, 932.0, 8.5, -56.8),
    ];
    let desired = [1.0, 1.0, 1.000_001_498_528_907_5, 1.000_000_013_194_575_4];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // but not where x / (x - 1) rounds to 1 (Arb: 8.12e10)
    assert!(xsf::hyp2f1(1e-20, -0.1, -2.5, -1e300).is_nan());
}

#[test]
fn test_hyp2f1_f64_nan_if_lost() {
    // unlike in xsf: NaN if both error estimates say that no digit is correct, instead of 2.36e5
    // and 4.83e51 (Arb: 1 and 1)
    assert!(xsf::hyp2f1(-3.53e-26, -1592.2, -17.239, 0.04516).is_nan());
    assert!(xsf::hyp2f1(-3.25e-18, 491.04, -0.010_315, -0.328).is_nan());
    // but `hyp2f1_recur` bounds its error for these, which were 7.7e149, 2.4e7, -5.6e85 and 7.5e-4
    // (Arb)
    let actual = [
        xsf::hyp2f1(6.08e-26, -908.762, -19.967, 0.5),
        xsf::hyp2f1(-50.7, 0.1, 344.3, -1.0),
        xsf::hyp2f1(1e-14, -1500.5, 0.7, 0.2),
        xsf::hyp2f1(90.0, -90.0, 180.7, 0.8),
    ];
    let desired = [
        1.0,
        1.016_019_727_590_008,
        0.999_999_999_999_930_7,
        3.771_813_135_532_291e-17,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // a series that doesn't converge is estimated as 1, but this is accurate (Arb)
    let actual = [xsf::hyp2f1(-8.999_999_999_999_99, 4.0, 2.5, -1.0)];
    np_assert_allclose!(actual, [2_246.910_399_198_305_6], rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_nan_if_err_above_1() {
    // NaN if the error estimate with the recurrences in `a` is above 1, instead of 9.47e-62 and
    // 3.62e16 (Arb: 1.18e-64 and -7.42e9)
    assert!(xsf::hyp2f1(702.2, 161.2, 0.1188, -0.4165).is_nan());
    assert!(xsf::hyp2f1(198.4, 38.54, -105.4, -0.5366).is_nan());
    // but not for a polynomial, here T_10 near a root, which is accurate in absolute terms (Arb)
    let actual = xsf::eval_chebyshev_t(10.0, -core::f64::consts::FRAC_1_SQRT_2);
    assert!((actual - 6.835_808_657_661_923e-16).abs() < 1e-15);
    // and `hyp2f1_recur` bounds its error for these (Arb)
    let actual = [
        xsf::hyp2f1(29.0, -29.1, 0.18, 0.86),
        xsf::hyp2f1(-38.9, 34.5, -0.07, 0.885),
        xsf::hyp2f1(-0.15, -45.9, 51.0, -0.94),
    ];
    let desired = [
        9.595_688_178_599_744,
        0.802_925_766_046_597_9,
        0.783_576_573_101_434_1,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
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
    // and of the 1 - x transformation near x = 1 if c - a - b is an integer (Arb)
    let actual = xsf::hyp2f1(
        -151.265_625,
        -167.078_125,
        -314.343_75,
        0.996_171_893_522_054_3,
    );
    np_assert_allclose!([actual], [1.140_175_032_534_544_3e-91], rtol = 1e-12);
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
