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
        (8.5, 6.5, -4.25, c64(-1.0, 0.5), c64(-14.914_154_962_801_387, 2.508_617_290_210_916_7)),
        (4.3, 4.4, -4.2, c64(-0.9, 0.8), c64(-13.468_832_852_470_712, -27.983_525_735_799_486)),
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
        (4.5, 2.5, 57.0, c64(0.64, 0.65), c64(1.122_209_544_812_728_6, 0.158_257_844_232_734_11)),
        (4.5, 2.5, 57.5, c64(0.64, 0.65), c64(1.121_218_275_094_999_2, 0.156_607_234_073_543_14)),
        // |z| > 1
        (3.5, 2.5, 40.0, c64(0.9, 0.6), c64(1.211_419_901_187_38, 0.188_557_093_819_745_4)),
        // c - a - b just misses an integer, so the 1 - z transformation is NaN
        (0.1, 0.2, 50.3, c64(0.64, 0.65), c64(1.000_254_319_428_871_4, 2.627_839_705_956_765e-4)),
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
            (-1_048_575.0, 1.0, 2.0, c64(1e-7, 0.0), c64(0.949_356_716_098_484_6, 0.0)),
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
    let actual = [xsf::hyp2f1(-0.5, 125.0, 62.5, 0.4), xsf::hyp2f1(-0.8, 19.7, 9.4, 0.9)];
    let desired = [0.453_216_679_866_897_07, -487.198_231_355_781_05];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but the recurrence of the power series is inaccurate here (mpmath)
    let actual = [
        xsf::hyp2f1(-1.000_000_000_000_11, 2.000_002, 1e-6, 0.999_999_7),
        xsf::hyp2f1(-49.5, 61.5, 6.5, 0.96),
    ];
    let desired = [-2_000_000.033_618_169_8, -8.458_289_503_668_05e-5];
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
    let desired = [-5.131_087_011_432_369e21, -0.093_683_512_477_367_49, 268_344_349.251_342_5];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but its error estimate shouldn't be too large if a starting value is tiny, so that the less
    // accurate alternatives aren't used instead (Arb)
    let actual = [
        xsf::hyp2f1(-9.5, 10.5, 0.5, 0.95),
        xsf::hyp2f1(10.000_000_001, 9.999_999_999_999_99, 1.5, -3.0),
        xsf::hyp2f1(9.999_999_999, 20.0, -10.5, -3.0),
    ];
    let desired = [-0.897_769_973_095_626_9, 4.289_686_233_768_906e-8, 0.001_570_408_561_128_409_7];
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
    let desired = [0.416_874_471_029_815_73, 6.875_723_600_176_229e-5, 8_488_942_027_952.015];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but here the starting values of their recurrence need a recurrence as well (mpmath)
    let actual = xsf::hyp2f1(-20.000_000_000_2, -14.47, 2.93, -5.48);
    np_assert_allclose!([actual], [64_520_380_601.436_584], rtol = 1e-12);

    // or need the series without that recurrence (mpmath)
    let actual = [xsf::hyp2f1(-19.6, -4.1, 16.7, -9.9), xsf::hyp2f1(-4.2, -19.6, 14.6, -9.6)];
    let desired = [-1_046_594.825_411_287_6, -4_908_149.039_599_091];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but the 1/x transformation if its error estimate is smaller than that of the Pfaff series,
    // e.g. where it cancels, as xsf does. These were 1.37682208856282, -5.197163034404184e15,
    // 7.081003593530957e-7, -183294.3855496918 and 4100386506.2481966 (Arb).
    let args = [
        (-1.8, -2.5, 5.0, -10.0),
        (9.5, -4.9, -8.7, -1000.0),
        (9.8, 7.4, 0.3, -5.0),
        (-2.6, -6.1, 7.1, -50.0),
        (-7.8, -5.6, 9.0, -100.0),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [
        1.376_822_088_562_981_8,
        -5_197_163_034_410_069.0,
        7.081_003_593_523_789e-7,
        -183_294.385_549_550_82,
        4_100_386_506.245_826_7,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // unlike in xsf, the Pfaff series is kept where both error estimates are too large, or `c - b`
    // is near a pole of its gamma function (Arb)
    let actual = [
        xsf::hyp2f1(
            0.076_598_182_228_762_47,
            7.362_007_163_390_169,
            -3.923_401_672_797_962_4,
            -153.095_573_546_133_72,
        ),
        xsf::hyp2f1(-3.2, 1e-12, -3.000_000_04, -10.0),
    ];
    let desired = [1.258_875_876_921_877e-6, 1.003_134_633_107_021_5];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // or the terms of the 1/x transformation are subnormal: this is 9.37e-324 (Arb)
    assert_eq!(xsf::hyp2f1(1.375, 6.875, 7.875, -1e235).to_bits(), 2);
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
    let desired = [-602_275.348_686_216_6, 4.298_609_352_869_04e65, 1.155_996_198_053_486_3];
    np_assert_allclose!(actual, desired, rtol = 5e-13);

    // and this was -2.74e6 (Arb)
    let actual = xsf::hyp2f1(-4.6, -20.600_000_000_1, 50.7, -10.0);
    np_assert_allclose!([actual], [-65.931_816_188_886_59], rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_x_le_minus_9_b_minus_a_int() {
    // unlike in xsf: these lost about log10(-x) digits in the Pfaff transformation, as
    // 1.48876236e-33, 1.32938105e-33, -6.54168264e11, inf, -1297.62704740, 8.4284e-7 and
    // 14742952.50335 (Arb)
    let actual = [
        xsf::hyp2f1(3.0, 3.0, 4.0, -3.7e11),
        xsf::hyp2f1(3.0, 3.0, 6.0, -8.2e11),
        xsf::hyp2f1(-0.5, -1.5, 4.0, -742_448_874.85),
        xsf::hyp2f1(-1.97, -1.97, 5.84, -2.34e13),
        xsf::hyp2f1(1.25, -0.75, -3.45, -58_000.0),
        xsf::hyp2f1(0.375, 3.375, -1.5, -1e15),
        xsf::hyp2f1(-2.25, 0.75, 7.5, -1e4),
    ];
    let desired = [
        1.488_762_887_248_107_1e-33,
        1.329_379_117_616_510_5e-33,
        -654_168_241_771.272_2,
        2.019_700_022_161_330_8e25,
        -1_297.627_047_387_899,
        8.408_236_071_439_006e-7,
        14_742_952.503_354_19,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
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
    assert_eq!(xsf::hyp2f1(1e-300, 4.000_000_000_000_008e30, 1e30, 0.25), 1.0);
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
fn test_hyp2f1_f64_tiny_negative_a_or_b() {
    // unlike in xsf, where a tiny negative `a` or `b` counts as the integer 0: these were inf, 1,
    // 1, 1, 1, 1, 1, -4.27e6, 1, 1 and -inf (Arb)
    let actual = [
        xsf::hyp2f1(26.5, -3e-21, 0.001, -2e4),
        xsf::hyp2f1(-120.5, -2.5e-28, -9.94, -0.5005),
        xsf::hyp2f1(-1e-20, -30.5, 2.5, -100.0),
        xsf::hyp2f1(-1e-20, -50.5, 3.5, -1.5),
        xsf::hyp2f1(-1e-20, 150.5, 0.7, 0.99),
        // also for small parameters, and if the other is a negative integer
        xsf::hyp2f1(-4e-23, -9.75, -5.5, -5800.0),
        xsf::hyp2f1(-7.0, -1e-27, 9.5, -48000.0),
        xsf::hyp2f1(8.25, -3e-17, 1.25, 0.999_99),
        xsf::hyp2f1(-3e-23, 7.0, 1.625, 0.999_98),
        xsf::hyp2f1(-0.25, -1e-25, -9.25, 0.9985),
        xsf::hyp2f1(-9.9e-14, 5.25, -0.0112, -3.75),
    ];
    let desired = [
        1.0,
        -10_864.775_782_791_036,
        -3.386_320_822_599_271e37,
        0.999_973_883_705_280_6,
        -1.562_482_690_977_387e278,
        -5.396_467_458_796_281e16,
        -8.719_200_573_711_808,
        -2.337_318_708_006_997_3e17,
        -27.935_689_206_440_796,
        1.146_412_218_027_53,
        0.999_999_999_991_500_6,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // and these are 1 to rounding, also if the 1/x transformation would lose `a` in `c - (c - a)`
    let actual = [
        xsf::hyp2f1(-5e-324, -9.75, -5.5, -5800.0),
        xsf::hyp2f1(-1e-300, -120.5, -9.94, -0.5005),
        xsf::hyp2f1(-4.8e-28, -1.1, 117.5, -1.63),
    ];
    np_assert_allclose!(actual, [1.0, 1.0, 1.0], rtol = 1e-15);
    // also where the error estimates of the paths for polynomials say that no digit is correct
    // (Arb)
    let actual = [
        xsf::hyp2f1(-3.53e-26, -1592.2, -17.239, 0.04516),
        xsf::hyp2f1(-3.25e-18, 491.04, -0.010_315, -0.328),
    ];
    np_assert_allclose!(actual, [1.0, 1.0], rtol = 1e-13);
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
        (195.0, -195.0, -200.5, 0.01),
    ];
    let desired = [
        0.489_514_398_515_199_2,
        0.595_230_641_982_666_2,
        0.787_059_088_545_409,
        -6_540_705.712_015_875,
        -0.019_154_323_187_421_253,
        6.722_964_015_601_578,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        // also with a and b swapped, which was already accurate
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 2e-14);
    }

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
    let actual =
        [xsf::hyp2f1(1e-200, 1e-200, -1000.5, 0.8), xsf::hyp2f1(1e-160, 1e-160, -1000.5, 0.8)];
    let desired = [-5.767_096_442_064_845e199, -5.767_096_442_064_845e279];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_neg_int_c_minus_a_or_c_minus_b() {
    // unlike in xsf: by `hyp2f1_recur` if c - a or c - b is only near a negative integer and
    // -1 < c - a - b < 0, where the recurrence in `a` is estimated as accurate, but isn't. These
    // are 1.18e11, 24.58, 1.86e7, 61.73 and -1.73e13 in xsf (mpmath, Wolfram)
    let actual = [
        xsf::hyp2f1(-40.5, 40.9, -0.1, 0.8),
        xsf::hyp2f1(-25.25, 25.7, -0.3, 0.75),
        xsf::hyp2f1(-30.5, 31.2, 0.2, 0.85),
        xsf::hyp2f1(35.6, -35.1, -0.4, 0.7),
        xsf::hyp2f1(-45.5, 45.99, -0.01, 0.8),
    ];
    let desired = [
        116.439_612_499_602_16,
        24.543_627_215_364_04,
        -15.110_473_882_473_225,
        74.102_598_993_798_93,
        -4.645_596_140_460_149,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // also for x near 1 (mpmath, Wolfram)
    let actual = xsf::hyp2f1(-59.25, 70.1, 10.1, 0.999);
    np_assert_allclose!([actual], [-1.622_052_617_419_225_7e-11], rtol = 1e-12);
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
fn test_hyp2f1_f64_large_parameters_symmetric() {
    // unlike in xsf: the backward recurrence in c, for b near a negative integer. With a and b
    // swapped, these were 7.3575092e-14, NaN and -6.95925e-24 (Arb)
    let args = [
        (36.0, -104.0 - 9e-7, 0.0074, 0.66),
        (61.0, -193.0 - 2.8e-7, 0.000_107, 0.64),
        (59.0, -173.0 - 7.6e-9, 0.0274, 0.64),
    ];
    let desired =
        [7.357_509_219_073_461e-14, -1.626_123_558_941_129_7e-26, -6.960_243_881_349_578e-24];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-12);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }

    // and in both orders, but less accurate, where xsf gives -3.75e-33 and 2.99e12 (Arb)
    let args = [(-206.0 + 5.2e-7, 47.0, 3.6e-7, 0.72), (142.0, -372.0 + 9.1e-5, 2.1e-6, 0.86)];
    let desired = [-6.514_961_820_519_327e-38, 7.812_813_451_522_689e-90];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-10);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }
}

#[test]
fn test_hyp2f1_f64_large_integer_a_near_integer_b() {
    // unlike in xsf: the 1 - x transformation, for an integer a > 0, where one of its terms is a
    // polynomial. These were 1.7143275e-23, NaN, -1.87323e-22, NaN and NaN, where xsf gives
    // 3.28e21, 7.58e22, 8.07e91, -2.43e20 and 1.04e33 (Arb).
    let args = [
        (94.0, -172.000_005_2, 0.000_69, 0.796),
        (144.0, -307.0 - 3.2e-8, 1.07, 0.83),
        (107.0, -144.0 - 7.3e-6, 0.145, 0.945),
        (85.0, -171.999_998, 0.0369, 0.8847),
        (100.0, -185.000_000_01, 0.0018, 0.866),
    ];
    let desired = [
        1.714_327_734_183_920_2e-23,
        -2.798_543_910_807_295e-65,
        -1.871_761_033_753_452_6e-22,
        -4.884_956_947_532_335e-39,
        1.330_815_061_015_127_8e-34,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-13);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }

    // but not at the cost of the other forms, whose bounds can be smaller where the two terms of
    // the transformation cancel, e.g. for x < -1 (Arb; xsf gives NaN)
    let (a, b, c, x) = (194.001, 350.0, 55.065, -5.0);
    let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
    np_assert_allclose!(actual, [-9.278_721_967_709_765e-268; 2], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_parameters_near_log_case() {
    // unlike in xsf: the 1/x and 1 - x transformations for b - a or c - a - b within 1e-3 of an
    // integer, where their two terms cancel. These were NaN, NaN, NaN, 4688.54, 2.89638e-42,
    // -2.16715878e-48, -1.88005587e-30 and NaN, where xsf gives NaN, NaN, NaN, 4.84e3, 1.03e-31,
    // 5.88e-38, -1.88005587e-30 and NaN (Arb).
    let args = [
        (1.07, 3.07 + 2.5e-13, 251.94, -524_766.9),
        (79.25, 85.25 + 4e-12, -27.19, -411.2),
        (-312.005, 9.5, -302.505 - 6e-14, 0.999_97),
        (47.25, 1.5, 48.75 - 2.5e-11, 0.999_999_98),
        (53.75, 27.75 - 3.5e-11, 0.22, -10.4),
        (-68.75, 107.26, 50.51 + 3.6e-11, 0.9976),
        (115.25, -73.5, 40.75 + 1.8e-6, 0.9998),
        (46.25, 76.25 + 5e-5, -300.6, -145.7),
    ];
    let desired = [
        0.000_131_036_002_073_166_09,
        3.918_853_590_964_261e-188,
        -1.038_538_785_178_233_9e19,
        4_688.172_838_647_12,
        2.896_395_178_110_93e-42,
        -2.167_158_642_601_661_3e-48,
        -1.880_055_697_263_691_3e-30,
        5.591_540_566_735_008e-62,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-13);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }
}

#[test]
fn test_hyp2f1_f64_small_parameters_near_log_case() {
    // unlike in xsf: also for small parameters, with c - a - b (for x >= 0.9) or b - a (for
    // x < -2) within 1e-3 of an integer. The first four were NaN, where xsf gives 2.27784513,
    // 4247.41591, -1673017043.2 and 12.27095. The others were -9.186935e-4, 1.31533e19, -883812,
    // 1.3971776e-3 and -64173143.8372, as in xsf except for its 1.31530e19 and -64179200 (Arb).
    let args = [
        (0.5, 0.25, 0.750_000_001, 0.999),
        (0.5, 2.5, 1.999_999_999_99, 0.9999),
        (0.5, -1.5, -0.999_999_999, 0.999),
        (0.5, 7.5, 7.999_999_999_99, 0.9999),
        (-8.75, 8.875, 2.125_000_01, 0.985),
        (-7.5, -7.500_000_000_01, 5.5, -1000.0),
        (-2.5, -2.500_000_000_01, 5.5, -1000.0),
        (1.5, 1.499_999_999_99, 5.5, -1000.0),
        (-9.25, -9.249_999_999_99, 0.75, -6.5),
    ];
    let desired = [
        2.277_844_918_975_149_3,
        4_247.415_771_798_882,
        -1_673_017_343.073_581,
        12.271_558_655_659_248,
        -9.192_870_681_661_742e-4,
        1.314_998_615_206_617_7e19,
        -883_900.539_293_860_8,
        1.397_157_529_639_970_1e-3,
        -64_173_143.838_225_69,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-14);
    }
}

#[test]
fn test_hyp2f1_f64_near_log_case_subnormal() {
    // b - a or c - a - b within a subnormal of an integer, where that distance underflows in the
    // combined form of `hyp2f1_recur`. Unlike in xsf, the first two aren't NaN (Arb).
    let args = [
        (5.4e-323, 7.0, -210.3, -12_265.5),
        (2.0, 3.5e-323, 275.9, -114_255.8),
        (4.535e-320, 155.75, 158.75, 0.999_987_6),
        (1.859_931_117e-315, -79.875, -78.875, 0.999_998_1),
    ];
    for (a, b, c, x) in args {
        np_assert_allclose!(
            [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)],
            [1.0; 2],
            rtol = 1e-12
        );
    }
}

#[test]
fn test_hyp2f1_f64_euler_less_accurate() {
    // unlike in xsf: the Euler transformation for c - a - b <= -1, where the error bound of
    // `hyp2f1_recur` shows its result to be less accurate, or it is NaN. These were -1.74008243e44,
    // 9.2789205253e35, 1.41743332e52, 1.02457247e39, NaN and NaN, where xsf gives 2.11e75,
    // 9.27905840e35, 1.60e95, -3.41e97, NaN and NaN (Arb; the last is -2.8e-380).
    let actual = [
        xsf::hyp2f1(149.0, -58.0 - 1.1e-7, 7.3e-7, 0.87),
        xsf::hyp2f1(189.0, -34.0 - 1.1e-8, 1.4e-7, 0.603),
        xsf::hyp2f1(185.0, -72.0 - 1.8e-9, 7.1e-7, 0.862),
        xsf::hyp2f1(136.0, -67.0 - 4.1e-7, 9.55e-5, 0.925),
        xsf::hyp2f1(9.25, 18.25 - 7.5e-8, -200.97, -90.43),
        xsf::hyp2f1(241.75, 220.75 + 2.65e-9, -3.67, -46.87),
    ];
    let desired = [
        -1.740_082_505_618_863_8e44,
        9.278_920_509_129_896e35,
        1.417_433_100_506_847_6e52,
        1.024_572_396_591_105_2e39,
        -3.130_221_533_372_943e-5,
        -0.0,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // with its power factor split in two, and for a subnormal `a`. These were -1.7311832e-211 and
    // NaN, where xsf gives 0.0 and NaN (Arb).
    let actual = [
        xsf::hyp2f1(98.0, 124.0 - 1.1e-13, 3.06, -71.03),
        xsf::hyp2f1(1e-323, 1.0, -159.8, -104_073.0),
    ];
    np_assert_allclose!(actual, [-1.725_017_797_892_166_1e-211, 1.0], rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_pfaff_recip_less_accurate() {
    // unlike in xsf: the Pfaff transformation for x < -1 and the 1/x one for x < -2, with large a
    // and b of the same sign, where the error bound of `hyp2f1_recur` shows their result to be
    // less accurate, or it is not finite. These were -9.13e168, 1.16e146 and -inf, as in xsf
    // (Arb).
    let actual = [
        xsf::hyp2f1(-129.666, -113.666, 56.25, -52.03),
        xsf::hyp2f1(-122.125, -86.125 + 6.9e-5, 78.26, -76.27),
        xsf::hyp2f1(-176.9375, -130.9375 + 3.3e-6, -21.0 + 5.2e-8, -35.47),
    ];
    let desired =
        [4.791_283_530_597_28e162, -1.531_524_570_886_761_4e144, 2.979_197_152_672_975e274];
    np_assert_allclose!(actual, desired, rtol = 1e-9);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_error() {
    // unlike in xsf: the psi function expansion for x > 0.9 and an integer c - a - b has an error
    // estimate, so this is replaced by the result of `hyp2f1_recur`. These were 1.36200088831e-8
    // and 4.8456802331e-31, as in xsf (Arb).
    let actual = [
        xsf::hyp2f1(
            -61.245_538_520_640_43,
            52.213_688_621_304_07,
            0.968_150_100_663_642_4,
            0.957_677_717_016_758_5,
        ),
        xsf::hyp2f1(
            -29.096_463_013_348_03,
            81.222_637_166_243_06,
            71.126_174_152_895_03,
            0.912_150_259_055_036_1,
        ),
    ];
    let desired = [1.362_000_886_533_284_2e-8, 4.845_680_233_889_792e-31];
    np_assert_allclose!(actual, desired, rtol = 1e-11);
    // and for the Pfaff transformation for x < -1, where these were -1.94e102, -1.95e210 and 0 (as
    // its factors underflow), while xsf gives -1.96e102, inf and 0 (Arb)
    let actual = [
        xsf::hyp2f1(-138.75, -85.75, 31.625, -13.86),
        xsf::hyp2f1(-127.75, -164.75 - 2.8e-14, 1.0 + 8.9e-9, -25.89),
        xsf::hyp2f1(-161.75, -108.75, 63.125, -14.8),
    ];
    let desired =
        [1.476_521_363_124_945_4e100, 3.311_772_741_737_614_5e207, 5.121_708_481_485_257e107];
    np_assert_allclose!(actual, desired, rtol = 1e-9);
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
fn test_hyp2f1_f64_huge_c_at_1() {
    // unlike in xsf: these are NaN in xsf, as `Gamma(c)` overflows in `F(a, b; c; 1)`
    // (DLMF 15.4.20; Arb)
    let actual = [
        xsf::hyp2f1(1.3, 6.5, 1e17, 1.0),
        xsf::hyp2f1(1.3, 6.5, 1e300, 1.0),
        xsf::hyp2f1(200.0, 0.5, 300.0, 1.0),
        xsf::hyp2f1(17.5, 250.5, 312.0, 1.0),
        xsf::hyp2f1(6.9e16, -469.5, 1.65e17, 1.0),
        // a term with an overflowing denominator, and a slowly decaying series (d = 2.7e-14,
        // 0 in xsf)
        xsf::hyp2f1(6e153, 6e153, 1e308, 1.0),
        xsf::hyp2f1(1e-15, 179.999_999_999_999_97, 180.0, 1.0),
        // a sum close to `f64::MAX`, an overflowing `a b`, and an underflowing `a / c`
        xsf::hyp2f1(1_946.607_092_231_488_2, 5_848_583_033.383_596, 19_141_581_561.342_18, 1.0),
        xsf::hyp2f1(2e154, 2e154, 1e308, 1.0),
        xsf::hyp2f1(1e-16, 2e307, 1e308, 1.0),
    ];
    let desired = [
        1.0,
        1.0,
        1.736_405_723_697_015_8,
        24_384_542_401_255.008,
        3.695_150_263_060_27e-111,
        1.433_329_414_560_340_4,
        1.036_467_456_652_845,
        1.797_693_134_862_289e308,
        54.598_150_033_144_25,
        1.0,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // overflows (+inf for both in xsf)
    let c = 5.967_056_889_428_375e-309;
    assert_eq!(xsf::hyp2f1(-1.25, 1.125, c, 1.0), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(-1.25, 1.125, -c, 1.0), f64::NEG_INFINITY);
    // symmetric in a and b, also for equal |a| and |b|
    let (a, c) = (190.513_976_388_490_14, 19_977.707_966_310_958);
    assert_eq!(xsf::hyp2f1(a, -a, c, 1.0), xsf::hyp2f1(-a, a, c, 1.0));
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
fn test_hyp2f1_f64_tiny_negative_c_large_parameters() {
    // unlike in xsf, where a tiny negative c counts as the integer 0 (Arb)
    let actual = [
        xsf::hyp2f1(28.5, 32.25, -6.5e-267, -1.85),
        xsf::hyp2f1(-44.5, 46.25, -3.5e-103, 0.6),
        xsf::hyp2f1(41.25, 47.5, -9e-292, -8.75),
        xsf::hyp2f1(-17.5, 37.25, -6e-271, -9.75),
        xsf::hyp2f1(44.75, -32.5, -3e-234, 0.83),
        xsf::hyp2f1(-33.5, 23.25, -8e-16, 0.79),
        xsf::hyp2f1(-40.5, -47.25, -4.5e-193, -1.8),
        // and where the inner function of the Euler transformation overflows
        xsf::hyp2f1(60.0, 2.5e-6, -2e-266, -7.7e4),
        xsf::hyp2f1(43.0, 9.25, -2.5e-285, -600.0),
        xsf::hyp2f1(-15.5, 38.0, -3.5e-257, -23.5),
    ];
    let desired = [
        -4.619_620_997_681_311e252,
        -2.351_118_090_908_530_3e103,
        7.447_544_914_471_166e247,
        -6.278_839_855_668_874e302,
        -6.149_027_963_048_359e238,
        302_521_514_418.673_6,
        -3.702_964_233_045_192e212,
        1.249_950_267_344_792e260,
        -2.981_261_708_812_470_7e249,
        -7.762_342_418_664_52e291,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
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
fn test_hyp2f1_f64_euler_inner_overflow() {
    // Arb; the inner function of the Euler transformation overflows here, but its product with
    // (1 - x)^(c - a - b) doesn't
    let actual = [
        xsf::hyp2f1(21.0, 43.0, 1e-276, -24.0),
        xsf::hyp2f1(0.5, 20.0, 1e-233, -7.5e4),
        xsf::hyp2f1(-15.5, 48.0, 5e-261, -78.0),
        xsf::hyp2f1(34.0, -1.25, 2e-278, -1.7e4),
        // also where (1 - x)^(c - a - b) underflows
        xsf::hyp2f1(47.0, 43.0, 2e-237, -3000.0),
        // also where the rounding of 1 - x matters
        xsf::hyp2f1(1e14, 10.0, 1e-300, -1e-12),
        // also where that product is small, and the polynomial that the inner function is
        // rounded to lacks the dominant part of the result
        xsf::hyp2f1(22.0, 4098.5, 1e-189, -8e4),
    ];
    let desired = [
        -3.248_055_641_674_947_7e228,
        -2.347_637_291_530_230_8e229,
        1.281_349_226_923_448_9e305,
        8.831_616_558_293_24e284,
        -1.238_995_047_182_927e80,
        3.787_960_193_849_875e270,
        5.393_827_268_558_922e-167,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_nan_if_lost() {
    // unlike in xsf: NaN if both error estimates say that no digit is correct (no input is known to
    // reach this, as `hyp2f1_recur` bounds the errors of those below), where xsf gives 6.75e10 and
    // -2.81e238 (Arb)
    let actual = [
        xsf::hyp2f1(-7.000_000_000_000_021, -208.513_348_904_413_08, 0.5, 0.471_817_095_982_536_96),
        xsf::hyp2f1(
            -17.000_000_000_000_007,
            738.317_738_730_952_2,
            -128.5,
            -0.362_556_171_941_935_9,
        ),
    ];
    let desired = [119_399_279_643.698_67, -210.050_825_314_093_3];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // and these, which were 7.7e149, 2.4e7, -5.6e85 and 7.5e-4 (Arb)
    let actual = [
        xsf::hyp2f1(6.08e-26, -908.762, -19.967, 0.5),
        xsf::hyp2f1(-50.7, 0.1, 344.3, -1.0),
        xsf::hyp2f1(1e-14, -1500.5, 0.7, 0.2),
        xsf::hyp2f1(90.0, -90.0, 180.7, 0.8),
    ];
    let desired = [1.0, 1.016_019_727_590_008, 0.999_999_999_999_930_7, 3.771_813_135_532_291e-17];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // a series that doesn't converge is estimated as 1, but this is accurate (Arb)
    let actual = [xsf::hyp2f1(-8.999_999_999_999_99, 4.0, 2.5, -1.0)];
    np_assert_allclose!(actual, [2_246.910_399_198_305_6], rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_nan_if_err_above_1() {
    // NaN if the error estimate with the recurrences in `a` is above 1 (no input is known to reach
    // this, as `hyp2f1_recur` bounds the errors of those below), but not for a polynomial, here
    // T_10 near a root, which is accurate in absolute terms (Arb)
    let actual = xsf::eval_chebyshev_t(10.0, -core::f64::consts::FRAC_1_SQRT_2);
    assert!((actual - 6.835_808_657_661_923e-16).abs() < 1e-15);
    // and `hyp2f1_recur` bounds its error for these (Arb)
    let actual = [
        xsf::hyp2f1(29.0, -29.1, 0.18, 0.86),
        xsf::hyp2f1(-38.9, 34.5, -0.07, 0.885),
        xsf::hyp2f1(-0.15, -45.9, 51.0, -0.94),
    ];
    let desired = [9.595_688_178_599_744, 0.802_925_766_046_597_9, 0.783_576_573_101_434_1];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
    // and these, which were NaN in one order of a and b, where xsf gives -2.29e82 and 1.03e-13
    // (Arb)
    let args = [(640.0, -128.04, -9.5, 0.4248), (-688.5, 147.62, -28.5, 0.3898)];
    let desired = [1.621_532_650_837_629_3e79, -1.825_420_437_271_943_6e-18];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-12);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }
}

#[test]
fn test_hyp2f1_f64_large_opposite_signs() {
    // unlike in xsf: the 1 - x transformation if 0 < x < 1, and the 1/(1 - x) transformation if
    // x < 0 (with a and b of the same sign then), where the Pfaff and Euler forms cancel
    // (references: Arb)
    let actual =
        [xsf::hyp2f1(-975.2, 241.4, 0.00135, 0.39), xsf::hyp2f1(-880.55, -75.65, 0.091, -0.472)];
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
    let desired = [-8_553.508_491_012_49, -1_050_587.466_964_691_6, -2.046_272_877_671_090_4e24];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // from the start values at c + 1 and c + 2, which only help since the other forms nearly
    // succeeded (reference: Arb)
    let actual = xsf::hyp2f1(130.25, 125.25, 1.96, -145.46);
    np_assert_allclose!([actual], [-3.043_786_288_636_934_7e-280], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_c_lt_0_second_solution() {
    // unlike in xsf: for c < 0, from the second solution at 0 and one term of the 1/x
    // transformation (or of the 1 - x transformation if 0 < x < 1), since the other term loses
    // precision; xsf gives 8.29e102, 6.42e201, NaN and NaN (references: Arb)
    let args = [
        (-404.0, -124.75, -561.75, -8.25),
        (-390.0, -183.25, -671.25, -22.5),
        (-342.5, -186.25, -639.75, -36.25),
        (761.5, 211.5, -483.25, 0.04),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [
        -3.962_625_779_100_878e105,
        2.525_038_208_734_701_5e202,
        -1.322_092_032_643_801_3e232,
        3.140_497_898_387_16e-6,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
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
    let desired = [43.700_421_521_573_67, 1.614_794_417_233_753e-4, 9_480_012.959_767_725];
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
    let actual = xsf::hyp2f1(-151.265_625, -167.078_125, -314.343_75, 0.996);
    np_assert_allclose!([actual], [1.867_105_850_532_836e-91], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_polynomials() {
    // unlike in xsf: the recurrence in the degree, since the series cancel (references: Arb)
    let actual = [
        xsf::hyp2f1(-522.0, -123.5, 533.5, -0.9),
        xsf::hyp2f1(-522.0, -123.5, 533.5, -0.5),
        xsf::hyp2f1(-3.75, 10_675.75, -0.25, -5600.0),
    ];
    let desired =
        [6.485_127_939_853_472_5e-52, -4.200_865_766_102_724e-45, -2.143_761_498_506_179_6e29];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but their sums first, which are more accurate if their terms have one sign, e.g. this was
    // inf (Arb)
    let actual = xsf::hyp2f1(-193.25, -361.875, -784.25, -8.5);
    np_assert_allclose!([actual], [5.602_051_816_485_101e124], rtol = 1e-14);

    // and in double-double if those and the recurrences are inaccurate, e.g. near a root, where
    // these were -3.5662163414410e-9, -3.4478745961026 and -388.081847556425 (Arb)
    let args = [
        (-8.0, 18.016_500_331_508_873, 10.805_295_997_850_628, 0.908_756_475_07),
        (15.695_195_439_032_172, -283.0, -1.239_002_148_045_167_6, 0.022_551_819_176_619_81),
        (-4.999_999_999_790_708, -270.0, -0.418_468_928_086_676_85, -0.026_109_555_998_151_627),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [-3.566_216_341_442_087e-9, -3.447_874_596_102_943_6, -388.081_847_556_405_76];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    // but not where the low word of `1 / (1 - x)` is subnormal and so inexact, where that gave
    // 1.3e-4 and 1.3e-8 (exact)
    let actual = xsf::hyp2f1(-7.692_307_692_292_307e-308, -1.0, 10.0, -1.3e308);
    np_assert_allclose!([actual], [1.999_956_492_953_669_5e-12], rtol = 1e-6);
    let actual = xsf::hyp2f1(-6.666_666_653_333_333e-306, -1.0, 100.0, -1.5e307);
    np_assert_allclose!([actual], [2.000_000_054_253_873e-9], rtol = 1e-9);
}

#[test]
fn test_hyp2f1_f64_large_polynomials_in_1_over_x() {
    // unlike in xsf: the polynomial in 1/x, i.e. the term of the 1/x transformation that is left
    // for an integer a <= 0, e.g. these were 1.717382430e38, 6.50e50 and 7.6598e35 (references:
    // Arb)
    let args = [
        (-390.0, -60.3187, -547.357, -4.72131),
        (-315.0, -71.5, -426.5, -5.75),
        (-253.0, -50.75, -446.5, -7.5),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired =
        [1.717_382_032_883_699_5e38, 6.196_191_652_916_035_5e50, 7.661_996_786_019_798e35];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // and for an integer c - b <= 0, also with 0 < x < 1 (in 1 - 1/x), e.g. these were NaN and
    // 1.057e-11 (Arb)
    let actual = [
        xsf::hyp2f1(-177.375, -108.5, -208.5, -18.25),
        xsf::hyp2f1(-9.75, -958.25, -1022.25, 0.987),
    ];
    let desired = [6.812_936_622_801_884e131, 6.764_515_364_329_936e-12];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_polynomials_fallback() {
    // unlike in xsf: the result of `hyp2f1_recur` above its tolerance if the other methods are
    // worse, e.g. these were -3.97e161, 1.30e11, -4.37e70 and 1.3003887417e7 (references: Arb)
    let args = [
        (-578.0, 43.03, -211.37, -0.07),
        (-155.0, -13.11, -268.18, -5.78),
        (-405.0, 284.75, -356.78, -0.16),
        (-102.0, -11.3, -173.19, -5.63),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [
        5.727_126_898_619_048e-4,
        -6_321_179_103.061_271,
        3.088_344_341_495_295_5e-21,
        13_003_890.426_056_918,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_recur_fallback() {
    // unlike in xsf: also for non-polynomials, e.g. these were NaN, NaN and -5.62543958e-4, as the
    // recurrence in `a` loses most digits (references: Arb, Wolfram)
    let args = [
        (-43.53, 28.38, 0.850_66, 0.8976),
        (-58.8, 58.95, 29.1511, 0.883),
        (-39.4122, 29.003_22, -0.466_54, 0.506_55),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired =
        [-7.920_164_213_871_762e-10, 1.473_319_142_367_457_4e-34, -5.625_453_919_713_718e-4];
    np_assert_allclose!(actual, desired, rtol = 1e-11);

    // and these were 9.47e-62 and 3.62e16 in xsf, and NaN before (Arb, Wolfram)
    let actual =
        [xsf::hyp2f1(702.2, 161.2, 0.1188, -0.4165), xsf::hyp2f1(198.4, 38.54, -105.4, -0.5366)];
    let desired = [1.176_016_677_157_160_2e-64, -7_416_153_466.305_625];
    np_assert_allclose!(actual, desired, rtol = 1e-9);
}

#[test]
fn test_hyp2f1_f64_divergent_at_1() {
    // unlike in xsf: the sign of the limit x -> 1- for c - a - b <= 0 (DLMF 15.4.21, 15.4.23; Arb
    // at x = 1 - 2^-3000), where xsf has +inf for the first two and NaN for the last two
    let inf = f64::INFINITY;
    let actual = [
        xsf::hyp2f1(-2.109, 0.891, -1.723, 1.0),
        xsf::hyp2f1(1.5, -0.5, 1.0, 1.0),
        xsf::hyp2f1(2.5, -1.5, 1.0, 1.0),
        xsf::hyp2f1(-2.5, 7.5, 3.5, 1.0),
        xsf::hyp2f1(202.75, 102.25, 303.0, 1.0),
        xsf::hyp2f1(-99.25, -44.5, -231.875, 1.0),
    ];
    assert_eq!(actual, [-inf, -inf, inf, -inf, inf, -inf]);
    // +inf at a non-positive integer c, as for x < 1
    assert_eq!(xsf::hyp2f1(0.5, 0.5, -3.0, 1.0), inf);
    assert_eq!(xsf::hyp2f1(0.5, 0.5, -3.0, 0.5), inf);
}

#[test]
fn test_hyp2f1_f64_near_polynomial() {
    // unlike in xsf: not the polynomial if a, b, c - a or c - b is only near a negative integer,
    // e.g. c - a = -58 - 5e-15 for the first. These are 3.35e11, NaN, NaN, 4.6e49, -inf and
    // -3.95e-6 in xsf (Arb)
    let actual = [
        xsf::hyp2f1(65.9, -57.6, 7.9, 0.861),
        xsf::hyp2f1(8.8, 5.6, 2.8, -233.9),
        xsf::hyp2f1(-3.000_000_000_000_001, 9.8, 7.0, -1.45),
        xsf::hyp2f1(97.7, -11.000_000_000_000_002, -27.3, -0.649),
        xsf::hyp2f1(-6.2, 3.4, -8.2, -4.4),
        xsf::hyp2f1(-49.6, -83.999_999_999_999_94, -91.2, 0.773),
    ];
    let desired = [
        -5.460_690_518_245_093e-10,
        -7.659_905_486_275_666e-18,
        26.565_667_750_000_028,
        -11_568.209_528_386_858,
        0.005_964_933_841_297_133,
        1.981_179_670_464_318_1e-29,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_tiny_c_integer_b() {
    // unlike in xsf: c - b = -n + c is only near a negative integer, where the polynomial misses
    // the tail of the series. These are inf in xsf (Arb).
    let actual = [
        xsf::hyp2f1(41.5, 4.0, -1.6e-154, -35_557.0),
        xsf::hyp2f1(-42.25, 9.0, -1.5e-14, 0.9998),
        xsf::hyp2f1(15.0, -34.25, -1e-15, 0.97),
    ];
    let desired =
        [6.500_894_648_275_649_5e-24, -2.692_035_181_549_083_2e-9, 7.064_137_056_127_026e-5];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_just_below_1() {
    // unlike in xsf: x within 1e-13 below 1 is not x = 1, where xsf gives inf for c - a - b <= 0
    // (the first three and the last) and the value at 1 else, here 32.27, 440.5 and -7.106 (Arb)
    let x = 0.999_999_999_999_99;
    let actual = [
        xsf::hyp2f1(0.5, 0.7, 1.0, x),
        xsf::hyp2f1(0.5, 0.7, 1.0, 0.999_999_999_999_999_9),
        xsf::hyp2f1(-2.109, 0.891, -1.723, x),
        xsf::hyp2f1(0.5, 0.49, 1.0, x),
        xsf::hyp2f1(1.5, 2.25, 3.76, x),
        xsf::hyp2f1(0.81, -0.04, 0.775, 0.999_999_999_999_997_6),
        xsf::hyp2f1(-3.5, 1.25, -2.5, x),
    ];
    let desired = [
        1_258.099_570_916_464,
        3_095.910_120_747_267_7,
        -11_864_866.701_448_103,
        9.532_032_686_970_835,
        115.683_428_221_038_92,
        -0.256_813_810_821_179_77,
        -44_263.439_899_139_13,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // but x = 1 is: Gauss's sum (DLMF 15.4.20) for c - a - b > 0, and the limit else
    let actual = [xsf::hyp2f1(0.5, 0.2, 1.0, 1.0)];
    np_assert_allclose!(actual, [1.449_724_260_959_791_1], rtol = 1e-15);
    assert_eq!(xsf::hyp2f1(0.5, 0.7, 1.0, 1.0), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_underflow() {
    // unlike in xsf: a gamma factor of the larger sum of the psi function expansion (AMS55
    // 15.3.10-12) underflows to 0, where xsf returns the other sum alone, here 0 (Arb)
    let actual = xsf::hyp2f1(-57.593_75, 251.375, 134.781_25, 0.962);
    np_assert_allclose!([actual], [1.447_158_991_775_232_9e-19], rtol = 1e-8);
}

#[test]
fn test_hyp2f1_f64_log_case_cancels() {
    // unlike in xsf: an integer b - a for x < 0, or c - a - b for 0 < x < 1 (the last two), where
    // the limit of the 1/x or 1 - x transformation cancels: xsf gives 0, NaN, NaN, 1.732e114,
    // 3.148e12 and -6.611e-36 (Arb)
    let actual = [
        xsf::hyp2f1(-130.25, -106.25, 88.5, -56.5),
        xsf::hyp2f1(105.25, 137.25, -182.8, -19.875),
        // (b - a is 8.5e-14 from an integer)
        xsf::hyp2f1(239.1, 215.1 - 8.5e-14, -205.6, -74.1),
        xsf::hyp2f1(-90.75, -126.75, 74.95, -27.75),
        xsf::hyp2f1(-96.25, 186.5, 31.25, 0.93),
        // (a polynomial of degree 446, and c - a - b is 4.3e-14 from an integer)
        xsf::hyp2f1(-446.0, 51.25, -35.75 + 4e-14, 0.565),
    ];
    let desired = [
        -9.451_533_499_353_577e147,
        -8.900_607_975_593_249e-76,
        -2.530_210_061_769_191_5e-297,
        1.678_815_006_969_739_4e106,
        26.381_995_650_576_812,
        -7.385_094_653_794_124e-36,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_c_near_pole() {
    // unlike in xsf: c near (not at) a non-positive integer -n, where 2F1 is finite, and its power
    // series can jump by about 1 / (c + n) at the term n + 1; xsf gives inf (Arb)
    let actual = [
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, -0.4),
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, -0.99),
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, 0.4),
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, 0.9),
        // (where the Pfaff transformation rounds c - b near -n)
        xsf::hyp2f1(7.25, 1e-13, -2.0 + f64::EPSILON, -0.7),
        xsf::hyp2f1(-2.75, 1e-13, -1.0 - f64::EPSILON, -0.3),
    ];
    let desired = [
        1.022_933_013_832_101_3,
        349.569_306_386_504_9,
        3_671_424_765.509_292_6,
        2.052_009_450_569_774_3e63,
        -362.316_350_911_936_1,
        112.785_458_910_322_3,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // where the Pfaff power factor underflows, or a is near (not at) an integer that would
    // terminate the series
    let c = -1.0 + 2f64.powi(-45);
    let actual = [
        xsf::hyp2f1(c + 61.0, c + 31.0, c, -300_000.0),
        xsf::hyp2f1(-1.0 + f64::EPSILON, 2.0, -2.0 + f64::EPSILON, -0.4),
    ];
    np_assert_allclose!(
        actual,
        [1.795_587_581_895_763e-295, 0.655_976_676_384_839_7],
        rtol = 1e-13
    );
    // where the terms are near overflow: 2F1(c, b; c; x) = (1 - x)^-b
    let c = -300.0 + 2f64.powi(-44);
    let actual = [xsf::hyp2f1(c, 300.0, c, 0.904)];
    np_assert_allclose!(actual, [2.082_716_166_487_401_1e305], rtol = 1e-13);
    // NaN if the series cancels more than its rounding estimate allows
    let y = xsf::hyp2f1(0.1234, 7.1, -8.000_000_000_000_002, 0.010_557_580_327_416_858);
    assert!(y.is_nan());
    // at c = -n it is a pole, unless the series terminates before the term n + 1
    assert_eq!(xsf::hyp2f1(1.25, 2.0, -45.0, -0.4), f64::INFINITY);
    let actual = [xsf::hyp2f1(-3.0, 2.0, -45.0 - 4e-14, 0.5)];
    np_assert_allclose!(actual, [1.068_974_630_021_141_6], rtol = 1e-15);
}

#[test]
fn test_hyp2f1_f64_polynomial_neg_int_c() {
    // unlike in xsf: the polynomial for a negative integer c below the integer a (or b), whose
    // terms can cancel completely; xsf gives NaN or no correct digit (Arb, Wolfram)
    let actual = [
        xsf::hyp2f1(-186.0, -26.75, -333.0, -9.0),
        xsf::hyp2f1(-2.5, -123.0, -618.0, -12.181_556_567_107_146),
        xsf::hyp2f1(-282.0, 497.137_097_191_493_7, -301.0, -0.029_204_871_958_424_083),
        // (by its Pfaff transformation)
        xsf::hyp2f1(-192.346_044_530_068_72, -303.0, -310.0, 0.290_678_693_500_093),
        // (b is an integer below c, so the polynomial is that of a)
        xsf::hyp2f1(-359.0, -401.0, -400.0, 0.012_918_049_583_765_393),
    ];
    let desired = [
        -5.330_688_234_478_712e36,
        21.568_509_940_337_676,
        1.457_505_009_443_554_2e-6,
        9.959_364_424_242_371e-29,
        0.009_283_916_262_782_09,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    // the sign of an overflow, also for x > 1 (about 1.8e623 and -3.1e379)
    assert_eq!(
        xsf::hyp2f1(-150.0, 39.883_200_484_784_936, -155.0, -8_782.652_226_040_95),
        f64::INFINITY
    );
    assert_eq!(xsf::hyp2f1(-227.0, -223.0, -400.0, 163.795_069_436_680_82), f64::NEG_INFINITY);
    // a subnormal b, both a and b polynomial (of degree 2, not 1e6), a degree of 1e6 whose terms
    // vanish after about 60, and x = 1 (DLMF 15.4.24)
    let actual = [
        xsf::hyp2f1(-1100.0, f64::from_bits(1), -1101.0, 2.0),
        xsf::hyp2f1(-1e6, -2.0, -1_000_001.0, 0.25),
        xsf::hyp2f1(-1e6, 2.0, -1_000_001.0, 0.25),
        xsf::hyp2f1(-29.0, -27.0, -57.0, 1.0),
    ];
    let desired = [
        223.050_453_452_347_38,
        0.562_500_374_999_625,
        1.777_776_592_593_777_8,
        1.995_525_599_273_078_8e-15,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    assert_eq!(xsf::hyp2f1(-20.0, -20.0, -30.0, 1.0), 0.0);
    // its terms cancel by 1e80, so NaN (or the exact sum), not that of the terms up to where the
    // rest is negligible
    let y = xsf::hyp2f1(-2048.0, 512.0, -2049.0, -0.25);
    assert!(y.is_nan() || (y / 2.530_769_231_467_989e-50 - 1.0).abs() < 1e-13);
    // a pole if a is not an integer, also if c - a rounds to one (in the Euler transformation)
    assert_eq!(xsf::hyp2f1(-1.000_000_000_000_2, 0.25, -8192.0, 0.5), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_polynomial_neg_int_c_overflow() {
    // unlike in xsf: the sign of an overflowing polynomial for a negative integer c below the
    // integer a or b, e.g. of 1 - 4x + 38/9 x^2 and 1 + 4/3 x + x^2 for the first two, where the
    // power series stops at its first infinite term, or the recurrence in a gives NaN; xsf
    // gives -inf for the first two (references: exact)
    let inf = f64::INFINITY;
    for ((a, b, c, x), desired) in [
        ((-2.0, -20.0, -10.0, 1e308), inf),
        ((-2.0, 2.0, -3.0, -1e308), inf),
        ((5.0, -5.0, -7.0, -1.271_610_395_627_417_3e231), -inf),
        ((-4.0, -5.0, -7.0, 4.334_366_077_440_268e137), inf),
        // (by `hyp2f1_recur`, whose first term overflows at the scale of its sums)
        ((-11.0, 24.0, -20.0, 1.146_319_809_459_681_3e308), inf),
        ((-32.0, -136.199_875_679_093_6, -34.0, 8.858_494_428_306_513e307), inf),
        // (just above the cutoff, by 3.5e-16 and 2.5e-16, which the double-double sum of
        // `hyp2f1_recur` resolves)
        ((-3.0, 2.0_f64.powi(-960), -4.0, 2.759_926_699_133_471_5e199), inf),
        ((-2.0, 35.0, -40.0, -1.491_883_142_149_389e154), inf),
    ] {
        assert_eq!(xsf::hyp2f1(a, b, c, x), desired);
        assert_eq!(xsf::hyp2f1(b, a, c, x), desired);
    }
}

#[test]
fn test_hyp2f1_f64_degree_1() {
    // unlike in xsf: `1 - a x / c` also for b = -1, where xsf gives NaN (references: exact)
    let actual = [
        xsf::hyp2f1(10.0, -1.0, -10.0, -1.0),
        xsf::hyp2f1(1.0, -1.0, -10.0, -10.0),
        xsf::hyp2f1(-10.0, -1.0, -10.0, 1.0),
    ];
    assert_eq!(actual, [0.0; 3]);
    // unlike in xsf: not rounded to 0 below 1e-15; xsf gives 0 (in this order of a and b), 0 and
    // 2.2e-16 (references: exact)
    let x = 1.0 - f64::EPSILON / 2.0;
    let actual = [xsf::hyp2f1(-1.0, 1.0, 1.0, x), xsf::hyp2f1(1.0, -1.0, 1.0, x)];
    assert_eq!(actual, [f64::EPSILON / 2.0; 2]);
    // (where `b x` is subnormal)
    let t = f64::from_bits(1);
    assert_eq!([xsf::hyp2f1(-1.0, t, t, 0.5), xsf::hyp2f1(t, -1.0, t, 0.5)], [0.5; 2]);
    let actual = [
        xsf::hyp2f1(-1.0, -18.0, -6.361_004_257_813_726, 0.353_389_125_434_095_84),
        xsf::hyp2f1(-8.0, -1.0, 7.0, -0.874_999_999_999_999_8),
    ];
    np_assert_allclose!(
        actual,
        [1.745_358_090_684_158e-16, 2.537_652_627_714_643_4e-16],
        rtol = 1e-15
    );
}

#[test]
fn test_hyp2f1_f64_c_eq_a_or_b_neg_int() {
    // unlike in xsf: for c = a = -m or c = b = -m, the polynomial of degree m of DLMF 15.2.5, also
    // for x > 1 and m = 0; xsf gives inf, inf, 0.5, 2^-10, inf, (1 - x)^20 and inf (references:
    // exact)
    let actual = [
        xsf::hyp2f1(10.0, -10.0, -10.0, 1.0),
        xsf::hyp2f1(-10.0, -10.0, -10.0, -10.0),
        xsf::hyp2f1(-10.0, 1.0, -10.0, -1.0),
        xsf::hyp2f1(-10.0, 10.0, -10.0, -1.0),
        xsf::hyp2f1(2.5, -3.0, -3.0, 3.0),
        xsf::hyp2f1(-10.0, -20.0, -10.0, 1.0 / 3.0),
        xsf::hyp2f1(0.0, 2.5, 0.0, 3.0),
    ];
    let desired = [184_756.0, 25_937_424_601.0, 1.0, 60_172.0, 225.0625, 44_458.0 / 59_049.0, 1.0];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    // (whose terms cancel by 1e18 and 3e20; xsf gives inf and NaN)
    let actual = [
        xsf::hyp2f1(-50.0, -17.704_152_726_819_963, -50.0, 1.124_144_544_693_068),
        xsf::hyp2f1(-51.0, -50.0, -50.0, 0.656_800_541_976_558_8),
    ];
    np_assert_allclose!(
        actual,
        [-4.974_610_285_018_047e-13, 4.887_745_972_669_103e-10],
        rtol = 1e-14
    );
    // (exactly 0; xsf gives 0.5, 2^-10 and NaN)
    let actual = [
        xsf::hyp2f1(-9.0, 1.0, -9.0, -1.0),
        xsf::hyp2f1(-9.0, -10.0, -9.0, 0.5),
        xsf::hyp2f1(-10.0, -3.0, -10.0, 1.0),
    ];
    assert_eq!(actual, [0.0; 3]);
    // (1 - x)^400 for a = b = c = -400, which overflows
    assert_eq!(xsf::hyp2f1(-400.0, -400.0, -400.0, 6.984_779_578_287_557), f64::INFINITY);
    // (the sum of (-1)^k for k <= 10001)
    assert_eq!(xsf::hyp2f1(-10_001.0, 1.0, -10_001.0, -1.0), 0.0);
    // ((1 - x)^m, where 1 - x rounds to 1; reference: Arb)
    let actual = [xsf::hyp2f1(-1e20, -1e20, -1e20, 1e-20)];
    np_assert_allclose!(actual, [0.367_879_441_171_442_33], rtol = 1e-15);
    // (3^m and 2^-m, which overflow and underflow)
    let actual = [-2.0, 0.5].map(|x| xsf::hyp2f1(-1e20, -1e20, -1e20, x));
    assert_eq!(actual, [f64::INFINITY, 0.0]);
    // (a subnormal (1 - x)^m; reference: mpmath)
    let actual = [xsf::hyp2f1(-1_000_001.0, -1_000_001.0, -1_000_001.0, 0.000_708_5)];
    np_assert_allclose!(actual, [1.559_546_624_838_864e-308], rtol = 1e-15);
    // (the smallest subnormal, which exp(ln((1 - x)^m)) rounds up to)
    let actual = xsf::hyp2f1(-1e7, -1e7, -1e7, 7.451_054_586_157_437e-5);
    assert_eq!(actual, f64::from_bits(1));
    // (a tiny x and a huge m, where `ln(1 - x)` must be accurate relative to x; reference: mpmath)
    let a = -12_772_208_543_222_726_656.0;
    let actual = [xsf::hyp2f1(a, a, a, 2.0_f64.powi(-54))];
    np_assert_allclose!(actual, [1.216_780_750_623_399e-308], rtol = 1e-15);
    // (NaN for c = -inf)
    assert!(xsf::hyp2f1(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, 0.5).is_nan());
}

#[test]
fn test_hyp2f1_f64_c_pole() {
    // unlike in xsf: +inf at a pole of c (where no a or b ends the series before it) for any
    // finite x != 0; xsf gives -inf for the first three, and NaN for the others
    let actual = [
        xsf::hyp2f1(10.0, 1.0, -10.0, 10.0),
        xsf::hyp2f1(1.0, 1.0, -1.0, 10.0),
        xsf::hyp2f1(5.0, 1.0, -9.0, 1.5),
        xsf::hyp2f1(750.274_991_146_857_7, 22.0, -670.0, -7_048.169_400_335_604),
        xsf::hyp2f1(-22.992_639_872_801_874, 54.0, -752.0, -6_207.512_565_903_636),
        xsf::hyp2f1(
            -4.673_419_326_051_837,
            289.810_719_716_530_9,
            -995.0,
            -0.855_413_125_222_624_2,
        ),
    ];
    assert_eq!(actual, [f64::INFINITY; 6]);
    // (also for a non-integer a next to the integer above c, and for c = 0, where xsf gives 1.25
    // and -inf; but NaN for c = -inf)
    let actual =
        [xsf::hyp2f1(-1.0 - f64::EPSILON, 1.0, -2.0, 0.5), xsf::hyp2f1(10.0, 1.0, 0.0, 10.0)];
    assert_eq!(actual, [f64::INFINITY; 2]);
    assert!(xsf::hyp2f1(1.0, 1.0, f64::NEG_INFINITY, 0.5).is_nan());
}

#[test]
fn test_hyp2f1_f64_neg_int_c_x_1() {
    // unlike in xsf: the Chu-Vandermonde product at x = 1 for a negative integer c below the
    // integer a or b, also for small parameters; xsf gives NaN (references: exact)
    let actual = [xsf::hyp2f1(-2.0, -3.0, -4.0, 1.0), xsf::hyp2f1(-3.0, -2.0, -4.0, 1.0)];
    assert_eq!(actual, [0.0; 2]);
    // (a nonzero product, where xsf is an ulp below)
    let actual = [xsf::hyp2f1(-4.0, 2.0, -6.0, 1.0), xsf::hyp2f1(2.0, -4.0, -6.0, 1.0)];
    assert_eq!(actual, [14.0 / 3.0; 2]);
}

#[test]
fn test_hyp2f1_f64_x_gt_1() {
    // unlike in xsf: NaN where 2F1 is complex, in either order of a and b; xsf gives inf
    for (a, b, c, x) in [(32.000_000_000_000_01, 75.0, 97.0, 2.0), (32.0, 75.0, 97.0, 2.0)] {
        assert!(xsf::hyp2f1(a, b, c, x).is_nan());
        assert!(xsf::hyp2f1(b, a, c, x).is_nan());
    }
    // (but real for a polynomial, and for (1 - x)^-2 (1 - x / 2) by DLMF 15.8.1; exact)
    assert_eq!(xsf::hyp2f1(-2.0, 1.0, 3.0, 3.0), 0.5);
    assert_eq!(xsf::hyp2f1(3.0, 1.0, 2.0, 3.0), -0.125);
    // and for c = b, (1 - x)^-a (DLMF 15.4.6), here 2^1040, even where c - a rounds to c
    let c = 2.0_f64.powi(60);
    assert_eq!(xsf::hyp2f1(20.0, c, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(c, 20.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(-3.0, 2.5, 2.5, 3.0), -8.0);
    // (and inf otherwise, as for the other complex powers, also where c - a - b rounds to an
    // integer, and for x = inf)
    let c = 2.0_f64.powi(52);
    for (a, b, c, x) in [(0.5, 2.5, 2.5, 3.0), (1.25, c, c, 3.0), (0.5, 2.5, 2.5, f64::INFINITY)] {
        assert_eq!(xsf::hyp2f1(a, b, c, x), f64::INFINITY);
        assert_eq!(xsf::hyp2f1(b, a, c, x), f64::INFINITY);
    }
    // inf for (1 - x)^(c - a - b) times a polynomial (DLMF 15.8.1) where c - a - b rounds, here
    // to -25.000000000000114, as it overflows (about 2^1294.6)
    let c = -1024.0 + 2.0_f64.powi(-43);
    assert_eq!(xsf::hyp2f1(c + 1.0, 24.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(24.0, c + 1.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    // (but an exponent within `HYP2F1_EPS` of an integer counts as one: (-1/4)^2 and 1)
    let actual = [xsf::hyp2f1(2.5, -1.999_999_999_999_999_8, 2.5, 1.25)];
    np_assert_allclose!(actual, [0.0625], rtol = 1e-14);
    assert_eq!(xsf::hyp2f1(-2.0_f64.powi(-55), 2.5, 2.5, 1.25), 1.0);
    // inf as in xsf where c - a rounds to c, as the skipped Euler transformation overflows (its
    // polynomial of degree 128 is about 2^1223)
    let c = 2.0_f64.powi(59);
    assert_eq!(xsf::hyp2f1(24.0, c + 128.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(c + 128.0, 24.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
}
