use core::f64::consts::PI;
use xsf::{hyp0f0, hyp0f1, hyp0f2, hyp1f0, hyp1f1, hyp1f2, hyp2f0, hyp2f1, hyp2f2, hyp3f2, hyppfq};

/// Asserts that `hyp3f2` is within a relative error of `1e-13` of `expected`
fn assert_hyp3f2(a: [f64; 3], b: [f64; 2], z: f64, expected: f64) {
    let actual = hyp3f2(a[0], a[1], a[2], b[0], b[1], z);
    assert!(
        (actual - expected).abs() <= 1e-13 * expected.abs(),
        "hyp3f2({a:?}, {b:?}, {z}) = {actual} != {expected}"
    );
}

#[test]
fn test_hyp3f2_series() {
    // mpmath references, cross-checked with Arb
    for (a, b, z, expected) in [
        ([1.0, 1.0, 1.0], [2.0, 2.0], 0.5, 1.164_481_052_930_025),
        ([0.5, 1.5, 2.5], [3.25, 4.75], -0.75, 0.923_404_451_082_218),
        ([-0.3, 2.2, 12.5], [0.7, -3.4], 0.6, -11_935_539_783.064_94),
        ([30.5, 0.25, -1.75], [20.25, 8.5], 0.95, 0.931_616_739_208_828_7),
        ([1.5, 2.5, 3.5], [-0.5, 6.0], -0.125, 1.256_376_554_812_066),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
    // the parts of the error estimate would overflow (Arb)
    let y = hyp3f2(2.0, 1.0, 900.2, 4.5, 350.2, 0.8);
    assert!(y.is_nan() || (y / 3.717_681_528_889_451e152 - 1.0).abs() < 1e-13);
}

#[test]
fn test_hyp3f2_polynomial() {
    // mpmath references, cross-checked with Arb
    for (a, b, z, expected) in [
        ([-4.0, 2.5, 3.25], [0.75, -6.5], 12.5, 355_339.477_100_550_95),
        ([-25.0, 10.5, 0.1], [30.25, 1.5], -3.0, 5_267_612.739_117_058),
        // terminates before the zero of (b_1)_n
        ([-3.0, -7.0, 1.5], [-5.0, 2.5], -40.0, 79_088.466_666_666_66),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
    assert_eq!(hyp3f2(0.0, 2.5, 3.5, -1.0, 2.0, 0.75), 1.0);
}

#[test]
fn test_hyp3f2_cancels_to_1f0() {
    // 1F0(a;;2) = (1 - 2)^-a, also past the radius of convergence
    assert_eq!(hyp3f2(2.0, 2.0, 3.0, 2.0, 3.0, 2.0), 1.0);
    assert_eq!(hyp3f2(3.0, 2.0, 3.0, 2.0, 3.0, 2.0), -1.0);
}

#[test]
fn test_hyp3f2_cancels_to_2f1() {
    // 2F1(2, 2; 1; z) = (1 + z) / (1 - z)^3, also past the radius of convergence
    assert_eq!(hyp3f2(2.0, 2.0, 3.0, 1.0, 3.0, 2.0), -3.0);
    assert_eq!(hyppfq([2.0, 2.0, 3.0], [1.0, 3.0], 2.0), Some(-3.0));
    assert_eq!(hyppfq([2.0, 2.0], [1.0], 2.0), Some(-3.0));
    // 2F1(5/2, 1; 3/2; z) = (1 - z/3) / (1 - z)^2
    assert_hyp3f2([1.0, 2.5, 1.0], [1.0, 1.5], 10.0, -7.0 / 243.0);
    // complex, where `hyp2f1` is +inf
    assert!(hyp3f2(1.75, 3.0, 4.0, 2.0, 4.0, 2.0).is_nan());
    assert!(hyppfq([1.75, 3.0], [2.0], 2.0).is_some_and(f64::is_nan));
    // complex, with c - a = -5 + 2^-51 (Arb: imaginary part -1.7e-13), which `hyp2f1` rounds
    let (a, c, z) = (-3.645_644_290_993_088_6, -8.645_644_290_993_088, 5.231_771_887_127_298);
    assert!(hyppfq([a, 3.0], [c], z).is_some_and(f64::is_nan));
    assert!(hyp3f2(a, 3.0, 4.0, c, 4.0, z).is_nan());
}

#[test]
fn test_hyp3f2_special() {
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 0.0), 1.0);
    assert!(hyp3f2(f64::NAN, 2.5, 3.5, 0.5, 4.5, 0.25).is_nan());
    assert!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, f64::NAN).is_nan());
    // the limits of `hyppfq` for an infinite argument
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, f64::INFINITY, 0.25), 1.0);
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, f64::NEG_INFINITY), 0.0);
    assert!(hyp3f2(f64::INFINITY, 2.5, 3.5, 0.5, 4.5, 0.25).is_nan());

    // a zero of (b_2)_n, also if the polynomial terminates after it
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, -2.0, 0.25), f64::INFINITY);
    assert_eq!(hyp3f2(1.5, -3.0, 3.5, 0.5, -2.0, 0.25), f64::INFINITY);

    // complex for z > 1
    assert!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 1.5).is_nan());
    // where the series diverges
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 1.0), f64::INFINITY);
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 1.5, 1.5, 1.0), f64::INFINITY);
    // with the sign of the late terms (gh-1645)
    assert_eq!(hyp3f2(-0.5, 2.5, 3.5, 0.5, 3.25, 1.0), f64::NEG_INFINITY);
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, -0.5, 3.25, 1.0), f64::NEG_INFINITY);
    assert_eq!(hyp3f2(-0.5, 2.5, 3.5, -0.25, 3.25, 1.0), f64::INFINITY);
    // the 2F1 left by a cancellation, with a zero excess
    assert_eq!(hyp3f2(-2.0, -0.5, 1.25, -2.0, 0.75, 1.0), f64::NEG_INFINITY);
    let (a, b) = (
        [1.256_488_541_180_713_6, -0.452_965_067_141_963_83, 0.175_520_285_123_398_8],
        [0.968_097_853_412_922_8, -1.366_013_156_809_728_5],
    );
    assert_eq!(hyp3f2(a[0], a[1], a[2], b[0], b[1], 1.0), f64::NEG_INFINITY);
    assert_eq!(hyppfq([0.5, 1.5, 2.5, 3.5], [-0.5, 1.25, 2.25], 1.0), Some(f64::NEG_INFINITY));
    // b_1 + b_2 - a_1 - a_2 - a_3 = 2.8e-17 rounds to 0; the result is 2.5e15 (Arb)
    let y = hyp3f2(5.0, 0.1, 0.2, 5.0, 0.300_000_000_000_000_04, 1.0);
    assert!(y.is_nan() || (y / 2.467_832_059_538_944e15 - 1.0).abs() < 1e-13);
    // the excess 1e-300 is below the rounding of the double-double sums
    let y = hyp3f2(0.750_000_000_000_000_2, 2.5, -1e-300, 1.750_000_000_000_000_2, 1.5, 1.0);
    assert!(y.is_nan() || (y / 0.499_999_999_999_999_85 - 1.0).abs() < 1e-13);

    // a common parameter cancels, also a non-positive integer
    for z in [-30.0, 0.25, 0.95] {
        assert_hyp3f2([1.5, 2.5, 3.25], [2.5, 6.5], z, hyp2f1(1.5, 3.25, 6.5, z));
        assert_hyp3f2([-2.0, 1.5, -3.0], [4.5, -2.0], z, hyp2f1(1.5, -3.0, 4.5, z));
        // two of them
        assert_hyp3f2([1.5, 2.5, -3.0], [2.5, -3.0], z, (1.0 - z).powf(-1.5));
    }
    // z = 1 with Gauss's sum (Arb reference), also for a tiny excess c - a - b (exact
    // reference), and near a polynomial, where it cancels (exact and Arb references)
    assert_hyp3f2([1.5, 2.5, 3.25], [2.5, 6.5], 1.0, 4.324_537_409_748_778);
    let half_eps = f64::EPSILON / 2.0;
    for (a1, a2) in [(0.1, 1.0), (1.0, 0.1)] {
        let b = [5.0, 1.100_000_000_000_000_3];
        assert_hyp3f2([a1, a2, 5.0], b, 1.0, 3.275_345_183_542_189_4e14);
    }
    for (a1, a2) in [(-1.0, 1.0 - half_eps), (1.0 - half_eps, -1.0)] {
        assert_hyp3f2([a1, a2, 5.0], [1.0, 5.0], 1.0, half_eps);
    }
    let a = [-1.0 + half_eps, 1.0 + 2.0 * half_eps, 5.0];
    assert_hyp3f2(a, [1.0, 5.0], 1.0, -2.0 * half_eps);
    // just below f64::MAX (exact reference)
    let y = hyp3f2(-536.0, -537.0, 5.0, 7.875_869_841_387_024, 5.0, 1.0);
    assert!(y.is_nan() || (y / 1.797_693_134_862_315_3e308 - 1.0).abs() < 1e-13);
    assert!(hyp3f2(1.5, 2.5, 3.25, 2.5, 6.5, 1.5).is_nan());
    assert!(hyp3f2(1.5, 2.5, 3.25, 2.5, 1.5, 1.5).is_nan());

    // z < -1 with a common parameter; Arb references
    // the factors of a coefficient of the continuation underflow, unlike their product
    let (a, b) = ([80.25, 90.5, 3.5], [50.125, 3.5]);
    assert_hyp3f2(a, b, -1000.0, -4.705_675_953_632_075e-278);
    // a_1 - a_2 is nearly an integer
    let a = [0.999_999_999_999_999_9, -20.000_000_000_000_004, 3.5];
    assert_hyp3f2(a, [-60.75, 3.5], -2.0, 0.599_745_273_229_746_2);
}

#[test]
fn test_hyp3f2_continuation() {
    // mpmath references at 50 and 120 digits, cross-checked with Mathematica
    for (a, b, z, expected) in [
        ([0.5, 1.25, 2.75], [3.25, 4.75], -100.0, 0.267_996_937_273_824_95),
        ([10.5, 20.25, 0.125], [5.75, 30.5], -50.0, 0.594_745_542_710_990_4),
        ([0.25, 1.375, 2.875], [1.75, 3.875], -1e6, 0.038_052_289_002_429_94),
        // b_1 - a_1 = -1, so that the term of a_1 vanishes
        ([2.5, 0.25, 1.125], [1.5, 4.25], -3.0, 0.822_861_429_867_887_1),
        // the factors of each coefficient underflow or overflow, unlike their product
        ([80.25, 120.5, 160.75], [150.125, 140.5], -1000.0, 1.689_865_300_058_061_6e-236),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
    // integer a_i - a_j, where it is a limit; Arb references, as the sum of the residues of the
    // integrand of DLMF 16.5.1
    for (a, b, z, expected) in [
        // Li_2(z) / z
        ([1.0, 1.0, 1.0], [2.0, 2.0], -50.0, 0.185_539_903_706_652_44),
        ([0.5, 1.5, 2.25], [3.25, 4.75], -100.0, 0.263_996_241_269_326_5),
        ([2.5, 0.5, 1.5], [1.25, 4.75], -1e6, 0.001_256_920_615_740_740_4),
        // b_1 - a_1 and b_1 - a_2 are integers too
        ([1.0, 3.0, 2.5], [2.0, 5.5], -30.0, 0.048_293_858_423_504_636),
        // a huge |z|, where the sequence transformation does not converge
        ([2e-15, 1.0, 1.0], [2.0, 2.0], -1e150, 0.999_999_999_999_313_2),
        // a_1 - a_2 is so close to an integer that its distance to it is lost in an `f64`
        ([1e-10, 1.0, 1.0], [2.0, 2.0], -1e10, 0.999_999_997_897_414_9),
        // residues 1 and 2 nearly vanish at this z, unlike those after them (mpmath reference)
        (
            [0.125, 0.125, 0.375],
            [-1.5, -0.333_950_201_884_253],
            -13.364_862_417_958_138,
            1.037_547_848_133_466,
        ),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp3f2_weniger() {
    // mpmath references at 50 and 120 digits, cross-checked with Mathematica
    for (a, b, z, expected) in [
        // -Li_2(-1) = pi^2 / 12
        ([1.0, 1.0, 1.0], [2.0, 2.0], -1.0, PI * PI / 12.0),
        ([0.5, 1.5, 2.5], [3.25, 4.75], -1.0, 0.902_702_134_732_558_7),
        ([0.5, 1.25, 2.75], [3.25, 4.75], -0.99, 0.910_875_206_300_700_3),
        // z < -1, where the terms of the continuation cancel
        ([0.5, 1.25, 2.75], [3.25, 4.75], -2.5, 0.820_925_159_215_043_4),
        // z < -1 with integer a_i - a_j
        ([1.0, 1.0, 1.0], [2.0, 2.0], -2.0, 0.718_373_183_441_840_4),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp3f2_double_double() {
    // the terms cancel in `f64`; mpmath references at 50 and 120 digits, cross-checked with Arb
    for (a, b, z, expected) in [
        ([3.0, 0.5, -19.5], [0.25, 0.125], 0.25, -1.187_310_231_031_717_4),
        ([10.5, 12.5, 9.25], [1.5, 2.5], -0.25, -0.037_909_280_209_963_286),
        ([20.5, 0.75, 3.25], [1.5, 2.5], -0.5, 0.099_410_623_211_045_2),
        // near a zero of the polynomial
        ([-2.0, 1.0, 1.0], [0.5, 0.5], 0.143_237_254_218_789_44, -4.403_648_831_270_954e-17),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp3f2_conformal() {
    // Wolfram references
    for (a, b, z, expected) in [
        // more than 10^4 terms
        ([1.5, 2.5, 3.5], [0.5, 4.5], 0.999, 2.209_183_333_563_127e8),
        ([0.5, 0.5, 0.5], [1.5, 1.5], 0.999_999, 1.088_792_804_483_117),
        // the terms cancel too much
        ([1.5, 2.5, 3.5], [0.5, 4.5], -1.0, -0.263_538_393_682_344_94),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp3f2_long_series() {
    // near z = 1, the terms cancel in `f64`, and the series in double-double takes more than
    // 2500 terms; Arb references
    for (a, b, z, expected) in [
        ([30.0, -0.015_625, -23.125], [-1.75, 0.359_375], 0.98, -167_348_338.199_215),
        ([12.5, 0.6875, -10.125], [-0.015_625, 0.171_875], 0.995, 13_118_233.367_026_178),
        // more than 10^4 terms
        ([0.468_75, -0.031_25, -0.8125], [-24.25, -17.5], 0.99, -2.407_886_558_460_386e90),
        // terms below half an ulp of the sum, which `f64` drops
        ([6e-30, 12.0, 13.0], [7.0, 11.0], 0.999, 1.000_000_000_000_594_9),
        ([2.4e-26, 12.0, 13.0], [7.0, 11.0], 0.996, 1.000_000_000_000_148_7),
    ] {
        assert_hyp3f2(a, b, z, expected);
        assert_eq!(hyppfq(a, b, z), Some(hyp3f2(a[0], a[1], a[2], b[0], b[1], z)));
    }
    // the terms cancel too much, also in double-double (Arb reference)
    let (y, expected) =
        (hyp3f2(-84.75, 59.5, 5.75, 17.25, 39.25, 0.996), 1.206_571_594_222_440_5e-6);
    assert!(y.is_nan() || (y - expected).abs() <= 1e-13 * expected, "{y}");
}

#[test]
fn test_hyp3f2_recurrence() {
    // the terms cancel too much, also in double-double; Arb references
    for (a, b, z, expected) in [
        // from the smaller lower parameter
        ([50.5, 60.5, 70.5], [1.5, 2.5], -0.75, -0.001_405_626_005_107_606_3),
        // from the larger one, as the steps to the smaller one cancel
        (
            [20.0, 19.0, -191.883_718_479_187_68],
            [-0.550_007_120_797_645_4, 47.286_239_465_585_32],
            0.164_230_372_686_943_13,
            -0.003_420_261_981_933_175_4,
        ),
        // from the second start, as the errors of the series of the first grow too much; the
        // b_1 + m + k of the starts are not floats
        (
            [94.378_191_005_498_73, 173.300_865_140_052, 147.500_025_613_686_18],
            [-13.643_774_007_711_395, 17.503_608_568_778_155],
            -0.008_800_179_444_837_896,
            17_600_470_648.833_354,
        ),
        // b_1 is close to a pole
        (
            [55.702_960_750_015_784, 104.952_956_865_680_7, 26.757_585_999_689_19],
            [-16.000_000_359_538_15, -0.049_520_245_084_025_07],
            -0.363_135_541_088_158_7,
            -5.398_186_830_203_316e27,
        ),
        // z is tiny, so that the low part of c_0 underflows, which adds little to the error
        (
            [2.427_499_880_600_77, -14.741_567_093_881_988, 25.381_848_363_326_903],
            [-1.143_065_210_718_003_2e-299, -2.751_195_480_316_527_5e-250],
            1.097_287_567_318_436_7e-307,
            -3.169_232_623_501_818e244,
        ),
        // down in a_1, as the steps in a lower parameter are unstable for z > 1/2
        ([-150.5, 20.25, 0.5], [2.5, 1.5], 0.75, 0.025_047_672_505_746_516),
        // z < 0, from the third start, as the series of the first two cancel too much
        ([-220.25, -250.5, 30.5], [2.5, 1.5], -0.6, -2.331_472_191_451_277e74),
        // in 8993 steps
        ([-9000.5, 2.25, 0.75], [2.5, 1.5], 0.5, 0.001_489_134_401_576_868_4),
        // z < 0 and large a_i > 0, where the steps down in b_j make the rounding errors too
        // large: up in a_1
        (
            [97.956_824_227_648_03, -0.616_795_639_499_172_6, -0.616_794_877_033_292_7],
            [0.063_737_300_150_794_8, -5.409_944_020_397_493],
            -0.942_430_082_997_230_7,
            86.828_604_252_284_17,
        ),
        // up in a_2 for each a_1 of the starts of the steps up in a_1
        (
            [86.669_015_218_836_36, 86.669_015_214_123_63, -0.342_161_258_796_916_06],
            [0.006_099_327_489_741_903, -24.871_773_151_539_312],
            -0.536_003_254_536_937_3,
            -227_601.611_548_551_53,
        ),
        // down in a large negative b_2 from where it is about -8, ..., with the starts up in the
        // large a_i, as the steps from b_2 + m > 0 make the errors too large
        (
            [-0.172_201_469_004_036_18, 96.137_189_053_415_8, -0.172_816_733_809_348_1],
            [0.002_365_014_890_838_624_6, -36.366_935_171_618_344],
            -0.388_254_698_398_966_54,
            29.149_920_307_359_714,
        ),
        (
            [49.967_682_358_836_456, 49.967_682_358_834_17, -0.003_228_150_097_381_314_3],
            [0.031_508_936_960_521_56, -85.215_606_245_586_29],
            -0.345_785_155_904_633_43,
            821_447.247_158_456_2,
        ),
        // down in a large negative a_1 from where the rounding errors of the steps are within
        // half of the tolerance, as from where it is about -8, ..., the series cancels too much
        // or the errors grow too much
        (
            [-42.218_236_520_288_855, 82.937_075_502_875_77, -0.010_699_330_775_812_416],
            [-0.007_539_046_179_480_939, 2.163_311_561_090_144_8],
            0.756_325_475_197_407_3,
            622_095_232.427_604,
        ),
        // the same from series of more than 2500 terms
        (
            [-0.114_022_201_911_358_77, 16.983_936_871_431_336, -36.194_302_717_639_12],
            [0.871_800_762_335_526_7, 0.001_410_048_592_030_205_8],
            0.977_528_259_874_388,
            194.132_778_948_619_8,
        ),
        // down in b_2 from where it is about -8, ..., with the starts from the steps down in b_1
        // from where it is positive, as the series cancel where both are negative
        ([0.5, 0.25, -0.125], [-50.75, -80.25], -0.75, 1.000_002_876_951_059_2),
        // the same with a large a_1
        ([80.5, 0.5, 0.25], [-30.75, -60.25], -0.5, -223_263.106_822_880_22),
        // up in all three a_i
        ([300.5, 200.25, 100.75], [0.5, 1.25], -0.25, -6.193_192_091_840_847e24),
        // from F = 1 at a_1 = 0
        ([90.0, 80.5, 0.5], [0.75, -1.25], -0.7, 0.006_795_836_258_444_880_5),
        // in 256 steps, as the errors of the first ones of all 280 grow too much
        ([280.0, 282.0, 247.5], [-0.25, 6.5], -0.2, 1.539_240_808_399_805e39),
        // in 8998 steps
        ([9000.5, 0.25, 0.75], [0.5, 1.25], -0.5, 0.078_582_730_877_915_69),
        // in 64 steps, as a series that all 70 start from is not accurate
        (
            [72.0, 0.029_710_275_578_575_418, 1.162_868_243_060_262],
            [-0.970_289_724_421_424_6, 0.162_868_243_060_261_93],
            -0.5,
            -7.683_075_807_302_188e-10,
        ),
        // near z = -1, from series of more than 2500 terms
        ([80.5, 0.125, 0.5], [-0.25, 0.25], -0.98, 0.523_156_315_386_979_5),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
    // b_1 is subnormal, so that the last step divides by a product that underflows (Arb)
    let y = hyp3f2(0.25, 1.2, 0.3, 1e-316, 2.0, 7e-9);
    assert!(y.is_nan() || (y / 3.150_000_064_610_025e306 - 1.0).abs() < 1e-13);
    // the steps up make the errors too large where b_2 - a_3 = -2, as F is then tiny (Arb)
    let y = hyppfq([90.0, 80.0, 0.5], [0.75, -1.5], -0.7);
    assert!(y.is_none_or(|y| (y / -2.343_795_263_623_968e-17 - 1.0).abs() < 1e-13));
}

#[test]
fn test_hyp3f2_nan() {
    // z = 1, and too close to it
    assert!(hyp3f2(1.5, 2.5, 3.5, 0.5, 8.5, 1.0).is_nan());
    assert!(hyp3f2(0.5, 0.5, 0.5, 1.5, 1.5, 1.0 - 1e-8).is_nan());
    // a huge negative parameter
    assert!(hyp3f2(0.5 - 1e12, 1.0, 1.5, 2.0, 2.5, -0.9).is_nan());
    // a method with a NaN error estimate is not trusted, but another one is (Arb reference)
    let y = hyp3f2(28.5, 1.25, 80.625, 20.125, 5.5, -1.125);
    assert!((y / 0.014_103_262_754_441 - 1.0).abs() < 1e-3);
}

#[test]
fn test_hyp3f2_underflow() {
    // the product of the a_i underflows, but not the terms; Arb reference
    let (a, b) = ([1e-170, 1e-170, 100.0], [1e-150, 2e-150]);
    assert_hyp3f2(a, b, 0.75, 8.034_690_221_294_95e19);
    // the terms underflow, and stay negligible
    assert_eq!(hyp3f2(-30.0, 1.0, 1.0, 2.0, 2.0, 1e-20), 1.0);
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 1e200, 1e200, 0.5), 1.0);
    assert_eq!(hyp3f2(1e-320, 1.0, 1.0, 2.0, 2.0, 0.5), 1.0);
    // the first term underflows, and later ones are large
    assert!(hyp3f2(1e-200, 1e-200, 500.0, 1.0, 1.0, 0.9).is_nan());
    // the first terms are negligible, and later ones overflow, as does the sum (Arb: 5.0e432)
    assert_eq!(hyp3f2(1e-30, 1e-30, 500.0, 1.0, 1.0, 0.9), f64::INFINITY);
    // the same at z = -1, where the series is summed by the sequence transformation
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 2.0, 1e300, -1.0), 1.0);
    assert!(hyp3f2(1e-200, 1e-200, -1500.5, 1.0, 1.0, -1.0).is_nan());
}

#[test]
fn test_hyp3f2_overflow() {
    // the sum overflows near z = 1, with the sign of the terms from where b_2 + n > 0 on (Arb:
    // 2.8e361 and -3.5e358; the others are above partial sums of 1e331 in Arb)
    assert_eq!(hyp3f2(90.5, 90.5, 0.25, 13.5, -47.75, 0.95), f64::INFINITY);
    assert_eq!(hyp3f2(90.5, 90.5, 0.25, 13.5, -46.25, 0.95), f64::NEG_INFINITY);
    assert_eq!(hyp3f2(50.0, 50.0, 50.0, 1.0, 1.0, 0.999_999), f64::INFINITY);
    let (a, b) = ([30.5, 20.25, 10.75, 5.5, 1.25], [-40.5, -30.25, 2.5, 3.5]);
    assert_eq!(hyppfq(a, b, 0.9999), Some(f64::INFINITY));
    // the terms shrink fast after a large t_1, but the sum is about 1 (Arb: |F - 1| < 5e-250)
    for a2 in [1e-150, -1e-150] {
        let y = hyppfq([1e300, a2, 1.0], [1e100, 1e100], 2e-200);
        assert!(y.is_none_or(|y| (y - 1.0).abs() < 1e-13), "{y:?}");
    }
}

/// Asserts that `hyp2f2` is within a relative error of `1e-13` of `expected`
fn assert_hyp2f2(a: [f64; 2], b: [f64; 2], z: f64, expected: f64) {
    let actual = hyp2f2(a[0], a[1], b[0], b[1], z);
    assert!(
        (actual - expected).abs() <= 1e-13 * expected.abs(),
        "hyp2f2({a:?}, {b:?}, {z}) = {actual} != {expected}"
    );
}

#[test]
fn test_hyp2f2_series() {
    // Arb references
    for (a, b, z, expected) in [
        ([1.0, 1.0], [2.0, 2.0], 0.5, 1.140_302_841_043_172_1),
        ([0.5, 1.5], [2.5, 3.25], -0.75, 0.936_756_382_954_840_9),
        ([-0.3, 2.2], [0.7, -3.4], 6.0, -176_725.747_541_789_55),
        ([2.5, -1.5], [3.25, 0.75], 40.0, 4_335_768_632_723.463),
        // in double-double
        ([1.5, 2.25], [3.5, 4.75], -20.0, 0.135_925_339_819_463_36),
        ([-2.5, 0.5], [1.25, 3.5], -15.0, 14.521_375_665_855_382),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp2f2_asymptotic_errors() {
    // two one-term algebraic series, whose phases make them add (Arb reference)
    let (a, c) = (-2.0 + 49.0 * 2f64.powi(-52), -1.0 - 2f64.powi(-51));
    assert_hyp2f2([a, c], [a + 1.0, c + 1.0], 49.0, 2.252_357_221_940_723_2e18);
    // and where they do cancel (Arb reference)
    let (a, c) = (-159.875, -159.875 + 2f64.powi(-45));
    assert_hyp2f2([a, c], [a + 1.0, c + 1.0], 355.0, 1.032_653_437_407_810_2e154);
    // the estimate is just too large, so this falls back to the series (Arb reference)
    let (a, c) = (-199.5, -199.4375);
    assert_hyp2f2([a, c], [a + 1.0, c + 1.0], 352.0, 1.334_619_473_051_105_3e153);
}

#[test]
fn test_hyp2f2_polynomial() {
    // Arb references
    for (a, b, z, expected) in [
        ([-4.0, 2.5], [0.75, -6.5], 12.5, 1_424.990_151_262_878_5),
        ([-25.0, 10.5], [30.25, 1.5], -3.0, 1_431.131_238_039_312_8),
        // terminates before the zero of (b_1)_n
        ([-3.0, 1.5], [-5.0, 2.5], -40.0, -266.098_412_698_412_7),
        // cancelled terminating polynomial longer than hyp1f1's 500-term cap (mpmath 50-dps)
        ([-600.0, 1.25], [-700.0, 1.25], 1000.0, 1.084_726_577_449_080_7e303),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp2f2_special() {
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, 4.5, 0.0), 1.0);
    assert!(hyp2f2(f64::NAN, 2.5, 3.5, 4.5, 0.25).is_nan());
    assert!(hyp2f2(1.5, 2.5, 3.5, 4.5, f64::NAN).is_nan());
    // the limits of `hyppfq` for an infinite argument
    assert_eq!(hyp2f2(1.5, 2.5, f64::INFINITY, 4.5, 0.25), 1.0);
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, 4.5, f64::NEG_INFINITY), 0.0);
    assert_eq!(hyp2f2(1.5, -2.5, 3.5, 4.5, f64::INFINITY), f64::NEG_INFINITY);
    assert!(hyp2f2(f64::INFINITY, 2.5, 3.5, 4.5, 0.25).is_nan());

    // a zero of (b_2)_n, also if the polynomial terminates after it
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, -2.0, 0.25), f64::INFINITY);
    assert_eq!(hyp2f2(1.5, -3.0, 3.5, -2.0, 0.25), f64::INFINITY);

    // a common parameter cancels, also at a pole (Arb references of 1F1)
    assert_hyp2f2([1.5, 2.5], [3.5, 2.5], -2.5, 0.401_651_681_832_597_7);
    assert_hyp2f2([1.5, 2.5], [3.5, 2.5], 3.0, 4.636_813_090_326_332);
    assert_hyp2f2([-3.0, 2.5], [4.0, -3.0], -12.5, 0.010_921_207_435_195_95);
    assert_eq!(hyp2f2(1.5, 2.5, 2.5, 1.5, -3.0), (-3.0_f64).exp());
    // cancelled 1F1(20; 1; -750) = -4.0345e-289 must not underflow to a false zero
    let y = hyp2f2(20.0, 6.25, 1.0, 6.25, -750.0);
    assert!(y.is_nan() || (y / -4.034_504_366_728_602e-289 - 1.0).abs() <= 1e-13, "{y}");
}

#[test]
fn test_hyp2f2_asymptotic_gamma_factors() {
    // Arb references
    for (a, b, z, expected) in [
        // a_1 - a_2 is small, so that the terms of the algebraic part cancel
        ([0.125, -0.046_875], [0.0625, -0.0625], -500.0, 0.091_856_794_837_619_5),
        ([-0.046_875, 0.1875], [-0.0625, 0.015_625], -5000.0, -0.305_131_497_954_354),
        ([-0.0625, -0.031_25], [0.125, 0.0625], -1e4, -0.375_850_505_307_578_44),
        ([-0.008, -0.005], [0.01, 0.0125], -16_000.0, 0.554_541_814_716_236_9),
        (
            [-0.008_629_729_408_648_11, -0.005_077_050_360_388_232],
            [0.010_536_830_963_144_745, 0.013_783_784_427_297_33],
            -16_634.047_727_591_347,
            0.570_530_374_389_448_6,
        ),
        // the Gamma factors of the algebraic part are large
        (
            [2.280_379_798_760_627_4, 2.394_547_702_434_450_4],
            [-1.949_645_532_661_522_3, -1.018_504_185_809_934_8],
            -2_369.596_792_449_503_5,
            1.501_206_621_894_914_6e-4,
        ),
        ([1.5, 1.75], [4.75, -8.75], -200.0, -0.077_679_580_102_834_04),
        ([-5.5, -4.25], [2.75, 7.75], -1000.0, -73_846_233.289_640_3),
        // Gamma(b_1) is near the overflow threshold
        ([-0.5, 0.25], [171.0, 10.5], 700.0, -1.502_346_983_364_605_8e100),
        ([2.5, 3.75], [171.6, 0.25], -1e4, -1.548_924_264_912_872_6e-5),
        ([0.25, 0.75], [1e-290, 1.5], -1e4, -2.870_495_954_723_938_5e288),
        ([0.25, 0.75], [-1e-290, 1.5], -1e4, 2.870_495_954_723_938_5e288),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp2f2_asymptotic() {
    // Arb references
    for (a, b, z, expected) in [
        ([1.5, 2.25], [3.5, 4.75], 500.0, 5.523_108_556_223_365e206),
        ([-10.5, 0.25], [20.0, 3.5], 300.0, -1.057_385_809_228_637_7e71),
        ([1.5, 2.25], [3.5, 4.75], -500.0, 0.001_936_926_332_203_224_2),
        ([0.5, -1.25], [2.0, 3.0], -100.0, 17.512_557_901_365_43),
        ([3.0, 5.5], [3.5, -4.25], -1000.0, -1.124_031_587_749_003_6e-8),
        ([2.5, 3.25], [-0.5, 1.75], -60.0, 4.645_783_675_124_214e-5),
        // a common parameter cancels, also at a pole
        ([3.0, 5.5], [3.0, 7.0], -200.0, 1.770_351_282_323_549_3e-10),
        ([-3.0, 2.5], [-3.0, 4.0], -300.0, 4.324_983_868_207_274_5e-6),
        // a_1 - a_2 is close to an integer
        ([0.6, 2.600_000_1], [1.25, 2.75], -1000.0, 0.010_879_355_879_731_946),
        ([0.3, 2.300_000_01], [1.25, 2.75], -1200.0, 0.112_436_870_589_296_16),
        ([3.600_000_01, 0.6], [1.25, 2.75], -1000.0, 0.008_367_721_726_240_012),
        // so are 1 + a_1 - a_2 and 1 + a_1 - b_1
        (
            [-4.875_000_000_000_001, 0.125],
            [0.124_999_999_999_999_97, 15.5],
            -300.0,
            1_452_433.832_806_379,
        ),
        // 1 + a_1 - a_2 is not 0, but rounds to it
        ([0.1, 1.1], [-0.2, 0.3], 500.0, -2.514_071_867_080_611e220),
        // c_3 of the exponentially large part vanishes, but not c_4
        ([-1.5, -1.5], [-2.5, 1.5], 60.0, -4.921_886_930_234_788e21),
        ([-1.5, -1.5], [-2.5, 1.5], 100.0, -4.114_326_732_425_412e38),
        // a_1 - a_2 is an integer, where the algebraic part is a limit
        ([1.0, 2.0], [3.5, 4.75], -200.0, 0.043_106_214_672_287_14),
        ([1.0, 1.0], [2.0, 2.0], -50.0, 0.089_784_773_406_593_58),
        ([0.5, 2.5], [1.25, 3.0], -1e4, 0.008_371_495_152_126_434),
        ([1.0, 3.0], [2.5, 4.75], 300.0, 1.900_125_531_576_091_6e123),
        // so is b_1 - a_1
        ([2.0, 3.0], [1.0, 4.5], -100.0, -2.565_124_215_174_344_6e-5),
        // b_1 - a_1 is so close to an integer that its distance to it is lost in an `f64`
        ([1.0, 1.0], [1e-10, 2.0], -1e10, -1.000_000_002_160_306_6),
        // the terms of the algebraic part are below eps after a few, long before the smallest
        ([0.25, 0.5], [1.75, -0.5], -1e5, 0.087_476_240_625_041_52),
        ([2.5, 0.5], [0.25, -0.25], -1e5, -0.001_784_104_610_406_174_4),
        // the terms of the series of a_2 grow again after a small one, from where
        // 1 + a_2 - a_1 + n is small and positive, and then shrink to smaller ones
        (
            [80.649_718_544_819_85, -1.079_334_838_212_048_5],
            [-462.060_008_081_445_2, -466.754_127_058_247],
            -7_953.603_219_495_746,
            -4.858_950_728_706_407e16,
        ),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
    // residues 3 and 4 of the algebraic part nearly vanish at this z, unlike those after them:
    // accurate or NaN (Arb reference)
    let (b, z) = ([-1.5, -8.886_363_636_363_637], -110.305_786_588_507_15);
    let y = hyp2f2(0.75, 0.75, b[0], b[1], z);
    let expected = -0.205_817_644_249_452_2_f64;
    assert!(y.is_nan() || (y - expected).abs() <= 1e-13 * expected.abs());
    // equal a_i within sqrt(eps) of b_2, where the algebraic part changes much with a_1
    let t = -(2.0_f64.powi(-40));
    let b = [t - 1.0, t + f64::EPSILON.sqrt()];
    assert_hyp2f2([t, t], b, 30.0, 18_916_422_175.581_42);
    // overflow
    assert_eq!(hyp2f2(1.5, 2.25, 3.5, 4.75, 800.0), f64::INFINITY);
    assert_eq!(hyp2f2(2.0, 4.0, 4.0, 7.5, 900.0), f64::INFINITY);
    assert_eq!(hyp2f2(0.25, -0.5, 1.5, 3.5, 1e5), f64::NEG_INFINITY);
    // also for a huge z, with the sign of Gamma(b_1)
    assert_eq!(hyp2f2(1.5, 2.0, 2.5, 3.5, 1e300), f64::INFINITY);
    assert_eq!(hyp2f2(1.5, 2.0, -2.5, 3.5, 1e300), f64::NEG_INFINITY);
}

#[test]
fn test_hyp2f2_asymptotic_double_double() {
    // Arb references
    for (a, b, z, expected) in [
        // the terms of a series of the algebraic part cancel
        ([-32.5, -29.75], [8.0, -3.75], -1500.0, 3.955_046_741_892_129e63),
        ([-43.25, -62.125], [47.5, 30.75], -1640.0, -3.765_474_203_422_35e50),
        // a_1 - a_2 is small, so that the terms of the algebraic part cancel
        ([-0.0075, -0.0065], [0.01, 0.00625], -2800.0, 0.007_365_834_079_958_2),
        ([0.0425, 0.0475], [-0.04, 0.0055], -3e5, 0.577_434_031_258_155_7),
        ([2.5, 2.500_000_001], [-0.75, 1.25], -1e5, -7.313_762_124_125_236e-13),
        // both, with Gamma functions far below zero
        ([-18.25, -19.234_375], [1.25, 18.25], -2700.0, 1.319_031_511_167_807_6e36),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp2f2_inaccurate() {
    // the series cancels, and |z| is too small for the asymptotic expansion: about 8 digits, too
    // few for `hyppfq` (Arb reference)
    let y = hyp2f2(2.5, 3.25, -0.5, 1.75, -40.0);
    assert!((y / 2.190_480_405_532_719e-4 - 1.0).abs() < 1e-7);
    assert_eq!(hyppfq([2.5, 3.25], [-0.5, 1.75], -40.0), None);
    // e^z times a polynomial, where the algebraic part vanishes
    assert!(hyp2f2(2.0, 3.0, 1.0, 1.0, -100.0).is_nan());
    // an error of 0.499 of the result, which is twice the exact one (Arb: -503_991_821.928_347)
    assert!(hyp2f2(-13.5, 1e-4, -39.5, -26.999_999_999_999_9, -74.5).is_nan());
    // the result underflows to the smallest subnormal
    assert!(hyp2f2(1.0, 2.0, 1.0, 2.0, -745.0).is_nan());
    // subnormal results, accurate to their spacing (mpmath references)
    assert_eq!(hyp2f2(1.0, 2.0, 1.0, 2.0, -709.0), 1.216_780_750_623_423e-308);
    assert_eq!(hyp2f2(80.0, 80.5, 100.0, 100.0, -1e6), 9.244e-321);
    // overflows that it detects (Arb: 6.730_996_089e594 and 1.306_688_482e355)
    assert_eq!(hyp2f2(100.0, 100.5, 1.0, 2.0, 750.0), f64::INFINITY);
    assert_eq!(hyp2f2(-160.5, 0.25, 1.0, 1.0, -1e4), f64::INFINITY);
}

/// Asserts that `hyp1f2` is within a relative error of `1e-13` of `expected`
fn assert_hyp1f2(a: f64, b: [f64; 2], z: f64, expected: f64) {
    let actual = hyp1f2(a, b[0], b[1], z);
    assert!(
        (actual - expected).abs() <= 1e-13 * expected.abs(),
        "hyp1f2({a}, {b:?}, {z}) = {actual} != {expected}"
    );
}

#[test]
fn test_hyp1f2_series() {
    // Arb references
    for (a, b, z, expected) in [
        // Si(4) / 4
        (0.5, [1.5, 1.5], -4.0, 0.439_550_784_737_263_3),
        (-0.3, [0.7, -3.4], 6.0, -1.790_186_724_968_910_7),
        (2.5, [3.25, 0.75], 400.0, 7_318_534_180_885_220.0),
        // in double-double
        (1.0, [1.5, 2.5], -30.0, 0.031_954_931_160_083_2),
        (1.5, [2.25, 3.5], -200.0, 0.000_984_924_460_046_181_6),
        (0.25, [-2.5, 4.75], -60.0, 0.588_049_559_007_606_4),
    ] {
        assert_hyp1f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp1f2_polynomial() {
    // Arb references
    for (a, b, z, expected) in [
        (-4.0, [0.75, -6.5], 12.5, 47.894_455_552_464_21),
        (-25.0, [30.25, 1.5], -300.0, 2_757_805_550.158_651_4),
        // terminates before the zero of (b_1)_n
        (-3.0, [-5.0, 2.5], -40.0, -8.261_375_661_375_661),
    ] {
        assert_hyp1f2(a, b, z, expected);
    }
    // its estimated absolute error is exactly the tolerance
    let y = hyp1f2(-2.0, 6.160_584_865_114_734, 7.458_449_265_939_545, 30.228_957_868_979_32);
    assert_eq!(y, 0.012_571_498_510_638_38);
}

#[test]
fn test_hyp1f2_special() {
    assert_eq!(hyp1f2(1.5, 3.5, 4.5, 0.0), 1.0);
    assert!(hyp1f2(f64::NAN, 3.5, 4.5, 0.25).is_nan());
    assert!(hyp1f2(1.5, 3.5, 4.5, f64::NAN).is_nan());
    // the limits of `hyppfq` for an infinite argument
    assert_eq!(hyp1f2(1.5, f64::INFINITY, 4.5, 0.25), 1.0);
    assert_eq!(hyp1f2(1.5, -2.5, 3.5, f64::INFINITY), f64::NEG_INFINITY);
    assert_eq!(hyp1f2(-3.0, 2.5, 3.5, f64::NEG_INFINITY), f64::INFINITY);
    assert!(hyp1f2(1.5, 0.25, 0.5, f64::NEG_INFINITY).is_nan());

    // a zero of (b_2)_n, also if the polynomial terminates after it
    assert_eq!(hyp1f2(1.5, 3.5, -2.0, 0.25), f64::INFINITY);
    assert_eq!(hyp1f2(-3.0, 3.5, -2.0, 0.25), f64::INFINITY);

    // a common parameter cancels, also at a pole (Arb reference of 0F1)
    assert_hyp1f2(3.5, [3.5, 4.0], -12.5, -0.024_903_987_549_189_11);
    assert_hyp1f2(-3.0, [4.0, -3.0], -12.5, -0.024_903_987_549_189_11);
}

#[test]
fn test_hyp1f2_asymptotic() {
    // Arb references
    for (a, b, z, expected) in [
        (1.5, [2.25, 3.5], -2000.0, 3.386_681_684_937_679e-5),
        (1.0, [1.5, 2.5], -900.0, 0.000_953_821_660_161_545),
        (12.5, [20.25, 30.5], -5000.0, 4.908_658_768_179_874e-16),
        // a terminating algebraic part
        (2.5, [0.25, 4.5], -400.0, 0.092_109_577_006_595_11),
        // no algebraic part
        (2.5, [1.5, 4.0], -500.0, -0.000_951_630_992_123_725_2),
        (3.5, [3.5, 4.0], -1000.0, 5.674_432_067_916_389e-6),
        // mpmath and Wolfram reference
        (0.7, [1.3, 2.9], -1e6, 6.306_440_261_147_195e-5),
        // the Gamma factors of the algebraic part are large
        (2.25, [-1.5, 9.5], -5000.0, 1.100_224_978_096_844_3e-6),
        (3.5, [9.75, -8.75], -500.0, 6_638.401_626_292_566),
        // (and for small parameters)
        (0.03125, [-0.0625, 0.1875], -1e5, -0.116_965_526_511_898_26),
        // algebraic part has zero sum with negligible error next to the exponential part (mpmath)
        (25.5, [1.25, 0.05], -4000.0, 5.803_456_048_071_905e18),
    ] {
        assert_hyp1f2(a, b, z, expected);
    }
    // overflow
    assert_eq!(hyp1f2(1.5, 0.5, 3.5, 2e5), f64::INFINITY);
    assert_eq!(hyp1f2(1.5, -0.5, 3.5, 2e5), f64::NEG_INFINITY);
    assert_eq!(hyp1f2(1.5, 2.25, 3.5, 2e5), f64::INFINITY);
    assert_eq!(hyp1f2(-0.5, 2.25, 3.5, 2e5), f64::NEG_INFINITY);
    // also for a huge z
    assert_eq!(hyp1f2(1.5, 2.5, 3.5, 1e100), f64::INFINITY);
    assert_eq!(hyp1f2(1.5, -2.5, 3.5, 1e100), f64::NEG_INFINITY);
}

#[test]
fn test_hyp1f2_inaccurate() {
    // the series cancels, and |z| is too small for the asymptotic expansion
    assert!(hyp1f2(25.5, -50.25, 50.05, -4000.0).is_nan());
    // the phase 2 sqrt(-z) of the oscillating part is large: about 11 digits (mpmath:
    // `gamma(1/4) x^(3/8) besselj(-3/4, 2 sqrt(x))` for x = -z)
    let y = hyp1f2(0.5, 0.5, 0.25, -1e40);
    assert!((y / -98_243.753_444_667_18 - 1.0).abs() < 1e-11);
    // overflows that it detects (Arb: 2.617_687_386e341 and 1.026_573_294e324)
    assert_eq!(hyp1f2(100.0, 200.0, 200.0, 5e5), f64::INFINITY);
    assert_eq!(hyp1f2(-120.5, 1.0, 1.0, -1e6), f64::INFINITY);
}

#[test]
fn test_hyp0f2() {
    // the series, its cancellation for z < 0, and the asymptotic expansion (Arb)
    for (b, z, expected) in [
        ([1.5, 2.5], 0.5, 1.137_183_373_732_688_9),
        ([0.25, -1.75], 3.0, 28.452_593_016_473_912),
        ([1.5, 2.5], -50.0, 0.986_667_956_983_830_3),
        ([-0.5, 3.25], -12.5, -4.132_555_224_551_472),
        ([-3.5, 1.25], 200.0, 75_161_504.300_073_98),
        ([1.5, 2.5], 1e4, 1.243_554_196_262_664_3e23),
        ([2.5, 3.5], -1e6, -1.115_800_988_799_774_4e55),
    ] {
        let actual = hyp0f2(b[0], b[1], z);
        assert!(
            (actual - expected).abs() <= 1e-13 * expected.abs(),
            "hyp0f2({b:?}, {z}) = {actual} != {expected}"
        );
    }
    assert_eq!(hyp0f2(1.5, 2.5, 0.0), 1.0);
    assert!(hyp0f2(f64::NAN, 2.5, 0.5).is_nan());
    assert!(hyp0f2(1.5, 2.5, f64::NAN).is_nan());
    // a pole
    assert_eq!(hyp0f2(-2.0, 1.5, 0.5), f64::INFINITY);
    // they overflow (Arb: -2.448e647 and 1.093e552)
    assert_eq!(hyp0f2(1.0, 1.0, -1e9), f64::NEG_INFINITY);
    assert_eq!(hyp0f2(200.0, 200.0, 4e8), f64::INFINITY);
    assert_eq!(hyppfq([], [1.5, 2.5], -0.75), Some(hyp0f2(1.5, 2.5, -0.75)));
}

#[test]
fn test_hyp2f0() {
    // the Borel sum for z < 0 (Arb: `(-1/z)^a1 U(a1, a1 - a2 + 1, -1/z)`, and mpmath)
    for (a, z, expected) in [
        ([1.5, 2.5], -0.01, 0.964_059_468_995_989_7),
        ([0.5, 1.25], -2.0, 0.606_215_898_670_491_7),
        ([2.5, -1.5], -50.0, 1_611.766_578_245_452_8),
        ([-0.75, 3.25], -0.5, 2.041_085_184_015_574),
        ([4.5, 0.25], -1e4, 0.071_227_775_698_633),
        // a polynomial
        ([-3.0, 2.5], 0.5, -1.109_375),
        ([-3.0, 2.5], -4.0, 2971.0),
    ] {
        let actual = hyp2f0(a[0], a[1], z);
        assert!(
            (actual - expected).abs() <= 1e-13 * expected.abs(),
            "hyp2f0({a:?}, {z}) = {actual} != {expected}"
        );
    }
    assert_eq!(hyp2f0(1.5, 2.5, 0.0), 1.0);
    // the series diverges
    assert!(hyp2f0(1.5, 2.5, 0.5).is_nan());
    assert!(hyp2f0(f64::NAN, 2.5, -0.5).is_nan());
    assert!(hyp2f0(1.5, 2.5, f64::NAN).is_nan());
    assert_eq!(hyppfq([1.5, 2.5], [], -0.75), Some(hyp2f0(1.5, 2.5, -0.75)));
    // a1 is a spacing from -31, past which a1 - (a1 - a2 + 1) + 1 rounds (Arb: 3.457_054_166e108)
    let (a1, a2) = (-31.000_000_000_000_004, -41.5);
    assert!((hyp2f0(a1, a2, -100.0) / 3.457_054_166_467_876e108 - 1.0).abs() < 1e-6);
    for a in [[a1, a2], [a2, a1]] {
        let actual = hyppfq(a, [], -100.0).unwrap();
        assert!((actual / 3.457_054_166_467_876e108 - 1.0).abs() < 1e-14, "{actual}");
    }
    // an a_i a few spacings from a non-positive integer, which chgu's rounded a1 - b + 1 moves
    // by a spacing, and in the last an a2 at 2.4e-9 from -7 that it moves by one: U has no
    // digits, or 8 to 11, where the residues are accurate (Arb)
    for (a1, a2, z, expected) in [
        (-5.999_999_999_999_998, -22.999_999_999_999_996, -1e10, 2.675_976_546_400_021e231),
        (1.776_356_839_400_250_5e-15, -20.000_000_000_000_004, -1e13, 2.160_851_060_908_083e262),
        (-15.000_000_000_000_012, -30.999_999_999_999_993, -1e7, 2.126_314_811_847_468_7e227),
        (-3.000_000_000_000_001_3, -5.875, -336_536.0, 63_343_197_259_831.39),
        (
            -7.999_999_999_999_998,
            -7.000_000_002_437_216,
            -118_035_487.239_557_88,
            -1.240_805_934_615_225_3e61,
        ),
    ] {
        let actual = hyp2f0(a1, a2, z);
        assert!((actual / expected - 1.0).abs() < 1e-14, "{actual}");
        assert_eq!(hyppfq([a1, a2], [], z), Some(actual));
    }
    // the result is subnormal (Arb: 2.792_686_307e-324), or underflows (Arb: 2.002_611_525e-334)
    assert!(hyp2f0(20.5, 200.5, -3.2e13).is_nan());
    assert_eq!(hyppfq([20.5, 200.5], [], -3.2e13), None);
    assert!(hyp2f0(20.5, 200.5, -1e14).is_nan());
    // the polynomial 1 + 2 z overflows (Arb: -2.000_000_000_000_000_02e308)
    assert_eq!(hyp2f0(-1.0, -2.0, -1e308), f64::NEG_INFINITY);
    // the Borel sum with fewer digits than `hyppfq` needs, for a_i that differ by an integer
    // (Arb)
    let actual = hyp2f0(0.5, 90.5, -1e4);
    assert!((actual / 0.001_055_556_987_423_832_8 - 1.0).abs() < 1e-12, "{actual}");
    assert_eq!(hyppfq([0.5, 90.5], [], -1e4), None);
    // chgu estimates 12 digits of the Borel sum, which has 11.8, where the residues of a_i that
    // differ by an integer are accurate (Arb)
    let (a1, a2, z) = (-33.469_052_875_157_17, -30.469_052_875_157_17, -36.589_122_667_683_775);
    let actual = hyppfq([a1, a2], [], z).unwrap();
    assert!((actual / -2.953_959_436_326_337e85 - 1.0).abs() < 1e-14, "{actual}");
    // the Weniger sum, with 3.7 digits, where hyppfq has U (Arb: 2.968_695_239_976_95e209)
    let (a1, a2, z) = (-62.084_104_456_978_59, -49.999_999_999_996_75, -249.109_088_375_704_2);
    let actual = hyp2f0(a1, a2, z);
    assert!((actual / 2.968_695_239_976_95e209 - 1.0).abs() < 3e-4, "{actual}");
    let actual = hyppfq([a1, a2], [], z).unwrap();
    assert!((actual / 2.968_695_239_976_95e209 - 1.0).abs() < 1e-14, "{actual}");
    // the zero 1 - 0.5 * 2 of a polynomial, where the error estimate is not 0
    assert!(hyp2f0(-1.0, 0.5, 2.0).is_nan());
}

#[test]
fn test_pfq_extreme_products() {
    // exact, where a partial product of the parameters is about 10^-320
    let (tiny, huge) = (1e-160, 1e160);
    assert_eq!(hyppfq([-1.0, tiny, tiny, huge, huge], [], 2.0), Some(-1.0));
    // exact, where a partial product of the lower parameters is 2^-1023
    let (tiny, huge) = (f64::MIN_POSITIVE / 2.0, 2.0_f64.powi(1023));
    assert_eq!(hyp1f2(-1.0, tiny, -huge, 1.0), 2.0);
    assert_eq!(hyp2f2(-1.0, 2.0, tiny, -huge, 1.0), 3.0);
    assert_eq!(hyp3f2(-1.0, 1.0, 2.0, tiny, -huge, 1.0), 3.0);
    // a term ratio that over- or underflows before it is multiplied by z (the exact sums of the
    // terms, from Python `fractions`, rounded)
    let (a, b, z) = ([-1.0, 2.0_f64.powi(600)], [2.0, 2.0_f64.powi(-425)], f64::MIN_POSITIVE);
    assert_eq!(hyppfq(a, b, z), Some(-3.0));
    let (a, z) = ([-3.0, -f64::from_bits(1)], 2.0_f64.powi(500));
    assert_eq!(hyppfq(a, [1.5], z), Some(2.640_637_803_924_586e127));
    // the power of 2 of a product exceeds `i32`, and the second term overflows
    let mut a = vec![huge; 4_198_404];
    a.extend([-1.0, 16.0]);
    assert_eq!(hyppfq(&a, [], 2.0), None);
    // subnormal parameters (Arb)
    for (actual, expected) in [
        (hyp1f2(1.5e-309, 1.5e-319, -2.102, -4.284), 188_562_247_355.745_18),
        (hyp2f2(-3.0, -1.3e-315, -1.5e-45, 8e-317, 0.829), -7.693_166_057_368_33e45),
        (hyp3f2(-2.534, -1.4e-312, -3.0, -1.1e-249, -5e-147, 0.019), -3.783_943_130_613_643e82),
    ] {
        assert!((actual - expected).abs() <= 1e-15 * expected.abs(), "{actual} != {expected}");
    }
}

#[test]
fn test_hyppfq_routing() {
    let z = -0.75;
    assert_eq!(hyppfq([], [], z), Some(hyp0f0(z)));
    assert_eq!(hyppfq([2.5], [], z), Some(hyp1f0(2.5, z)));
    assert_eq!(hyppfq([], [1.5], z), Some(hyp0f1(1.5, z)));
    assert_eq!(hyppfq([0.5], [1.5], z), Some(hyp1f1(0.5, 1.5, z)));
    assert_eq!(hyppfq([0.3, 0.7], [1.5], z), Some(hyp2f1(0.3, 0.7, 1.5, z)));
    assert_eq!(hyppfq([0.5], [1.5, 2.5], z), Some(hyp1f2(0.5, 1.5, 2.5, z)));
    assert_eq!(hyppfq([0.5, -1.5], [1.5, 2.5], z), Some(hyp2f2(0.5, -1.5, 1.5, 2.5, z)));
    assert_eq!(hyppfq([0.5, 1.5, 3.5], [2.25, 4.5], z), Some(hyp3f2(0.5, 1.5, 3.5, 2.25, 4.5, z)));
}

#[test]
fn test_hyppfq_2f0() {
    // the Borel sum for z < 0 (Arb, as x^a1 U(a1, a1 - a2 + 1, x) with x = -1/z), symmetric in
    // a1 and a2, also where x^a1 overflows (series)
    for (a, z, expected) in [
        ([0.3, 1.7], -0.5, 0.848_525_611_949_669_8),
        ([-2.5, 0.75], -4.0, 88.108_287_408_568_78),
        ([40.0, 2.5], -1e-3, 0.906_748_253_693_675_3),
        ([1.5, 2.5], -1e-300, 1.0),
        ([1.5, 2.5], -1e-308, 1.0),
        ([1.5, 2.5], -1e-320, 1.0),
        // where `chguit` has more than 9 digits
        ([2.25, 3.25], -2.0, 0.025_935_294_832_104_445),
        ([1.25, 2.25], -1e3, 1.557_385_299_869_812_4e-4),
        ([2.5, 3.25], -2.0, 0.018_699_308_601_436_784),
        ([3.5, 7.5], -0.5, 0.008_273_723_913_148_375),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // with a1 and a2 swapped, as `chgu` estimates 12 digits for U(a1, a1 - a2 + 1, x), by `chgus`
    // and `chguit` (Arb)
    for (a, z, expected) in [
        ([7.3, -7.69], -2.7, 169_644_321_290.509_92),
        ([4.3, 7.0], -0.97, 5.585_769_788_936_731e-4),
    ] {
        let actual = hyppfq(a, [], z).unwrap();
        xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
    }
    // for a1, a2 << 0, where U is by the recurrence in a, with more than 12 digits (Arb)
    for (a, z, expected) in [
        ([-36.96, -50.1], -0.009_96, 2.252_440_982_724_057e-13),
        ([-37.14, -63.61], -0.0107, 4.832_372_813_708_457e-15),
        ([-68.36, -79.5], -0.009_69, 4.287_653_719_473_163_5e-21),
        ([-77.66, -116.73], -0.004_77, -4.990_584_018_414_777_5e-29),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // but not the Borel sum where U by the recurrence is 1.1e-13 off (Arb)
    let actual =
        hyppfq([-66.955_195_369_874_3, -86.323_704_933_975_42], [], -0.356_603_437_508_210_8);
    xsref::np_assert_allclose!([actual.unwrap()], [-3.721_878_069_495_089e77], rtol = 1e-14);
    // nor where it is 2.2e-13 off by the rounding of the coefficients of the recurrence (Arb)
    let a = [-121.661_194_083_177_43, -125.496_107_001_270_77];
    let actual = hyppfq(a, [], -0.011_817_257_637_100_627);
    assert!(actual.is_none_or(|y| (y / 6.013_041_362_636_676e-17 - 1.0).abs() < 1e-13));
    // nor where the rounding of x = -1/z or a1 - a2 + 1 changes it by more than 1e-13 (Arb)
    for (a, z, expected) in [
        (
            [-1.968_075_690_757_376e-12, -6.000_000_000_001_97],
            -734.029_490_840_056_6,
            -3.700_105_420_631_436e7,
        ),
        (
            [-129.715_870_814_732_23, -147.715_870_814_734_2],
            -0.009_864_846_651_846_624,
            1.152_931_589_568_623_8e-21,
        ),
        (
            [-12.000_000_000_000_002, -41.429_461_690_122_9],
            -9.964_632_067_709_93,
            -2.947_780_833_978_018_3e65,
        ),
        ([-1e-48, -20.5], -30.0, -0.069_303_167_698_671_19),
    ] {
        if let Some(actual) = hyppfq(a, [], z) {
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-13);
        }
    }
    // the Weniger fallback (Arb)
    for (a, z, expected) in [
        ([0.05, -0.06], -0.1, 1.000_286_818_712_551_8),
        ([0.25, 0.25], -0.1, 0.994_172_149_016_540_8),
        ([-0.25, 0.5], -0.075, 1.009_016_452_746_753),
        ([0.75, -0.5], -0.08, 1.029_048_187_106_747_8),
        ([-2.4, 2.75], -0.1875, 2.910_254_702_412_351),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // the Weniger fallback in double-double, for z < -1 and where the terms cancel (Arb)
    for (a, z, expected) in [
        ([1.1, 1.15], -1.09, 0.513_551_761_026_235_1),
        ([-1.47, -2.49], -1.76, -0.976_714_154_400_257_8),
        ([-17.57, -20.06], -1.04, 5.336_115_439_600_885_6e16),
        ([0.7, 0.8], -0.46, 0.835_046_888_976_750_3),
        ([-3.29, -13.69], -0.178, -1.070_282_071_194_967_9),
        ([-8.27, -33.08], -0.067, -0.114_157_126_503_264_89),
        ([-14.93, -14.85], -0.32, -31_736.919_105_697_503),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // where it stops with steps below its noise, but 1.2e-13 off by truncation (Arb)
    let y = hyppfq([-28.697, -30.146], [], -16.155);
    assert!(y.is_none_or(|y| (y / -4.139_259_483_009_185e66 - 1.0).abs() <= 1e-13));
    // and where only its last two steps are below its noise, 1.4e-13 off (Arb)
    let y = hyppfq([-13.0541, -13.0541], [], -8.2586);
    assert!(y.is_none_or(|y| (y / -3.712_522_542_395_385e19 - 1.0).abs() <= 1e-13));
    // a terminating sum whose drift bound would overflow unscaled (Arb)
    let actual = hyppfq([-125.0, -1000.0], [], -0.3).unwrap();
    xsref::np_assert_allclose!([actual], [-8.284_573_243_742_23e305], rtol = 1e-14);
    // the Weniger fallback, after 1549 terms (Arb)
    let y = hyppfq([-1548.5, 1.000_000_000_000_028_4], [], -0.001_953_124_999_999_999_8);
    assert!(y.is_none_or(|y| (y / 1.320_723_958_873_985_9e296 - 1.0).abs() <= 1e-13));
    // `chgu` estimates too few digits, and `(-1/z)^a1 U` gives 5.9e-19 (Arb:
    // -1.000_327_845_432_870_3e-19)
    assert_eq!(hyppfq([-10.5, -9.5], [], -0.033_448_280_861_346_91), None);
}

#[test]
fn test_hyppfq_2f0_borel_digits() {
    // where chgu claims the digits of `PFQ_TOL`, but the Borel sum is 1.5e-13 off (Arb)
    let y = hyppfq([-6.004_858_572_925_059, -21.079_297_738_677_543], [], -151_588_853.883_445_23);
    assert!(y.is_none_or(|y| (y / -1.040_990_391_775_226_8e184 - 1.0).abs() <= 1e-13));
}

#[test]
fn test_hyppfq_2f0_weniger_f64() {
    // where the noise that stops Weniger in f64 is below its error of 2.0e-13 and 1.4e-13 (Arb)
    let z = -59_302.810_486_050_94;
    for (a1, expected) in [
        (-1.368_297_746_698_496_2e-15, 0.429_248_028_314_945_7),
        (-1.137_055_427_506_450_3e-15, 0.525_705_111_529_719_7),
    ] {
        let y = hyppfq([a1, -2.999_999_999_999_947], [], z);
        assert!(y.is_none_or(|y| (y / expected - 1.0).abs() <= 1e-13), "{y:?}");
    }
}

#[test]
fn test_hyppfq_2f0_hard_cases() {
    // Cases where `chgu` lacks the extra digit margin; `hyppfq` still resolves these to 1e-14
    // via Weniger dd, `pfq_borel_dd`, or cached `margin = 0` U (Arb references).
    for (a, z, expected) in [
        ([2.2, 2.0], -1.5, 0.111_476_049_921_291_88),
        ([2.0, 3.2], -1.6, 0.051_620_474_992_384_33),
        ([-44.04, -60.33], -0.0903, 1.535_642_058_980_361_2e15),
        ([4.0, 4.5], -2.3, 5.149_428_916_838_061e-4),
        ([2.9, 2.8], -2.8, 0.009_916_196_734_627_417),
        ([7.6, 6.8], -3.3, 1.827_708_981_196_445e-8),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // Harder cases where `chgu` estimates < 13 digits: `hyppfq` stays conservative (`None` or
    // accurate to 1e-13), while `hyp2f0` still returns a finite accurate value.
    for (a, z, expected, rtol) in [
        ([-25.76, -92.98], -0.0111, 3.688_195_274_697_654_4e-13, 1e-13),
        ([-51.8, -36.5], -0.026, -6.120_387_741_382_293e-9, 1e-13),
        ([-81.6, -89.2], -0.002_72, 4.957_519_047_755_375e-13, 1e-13),
        ([-74.73, -90.73], -0.018_23, -1.280_977_214_916_313_9e-9, 1e-13),
        ([-111.52, -110.52], -0.0137, -2.235_055_638_338_494e-13, 1e-13),
        ([-122.89, -56.73], -0.0238, 519.670_814_874_989_4, 1e-11),
        ([-1606.98, 114.09], -6.78e-8, 1.012_508_691_661_776_5, 1e-13),
        ([445.47, -1302.36], -4.21e-5, 63_016_049_258.297_64, 1e-13),
        ([-100.0, -39.22], -0.0268, 0.312_116_603_209_507_8, 1e-13),
        ([5.5, 7.5], -1.0, 6.067_555_507_251_268_4e-5, 1e-12),
        ([7.0, 9.0], -1.8, 6.911_864_984_250_697e-8, 1e-13),
    ] {
        for a in [a, [a[1], a[0]]] {
            let p = hyppfq(a, [], z);
            assert!(p.is_none_or(|y| (y - expected).abs() <= 1e-13 * expected.abs()));
            let h = hyp2f0(a[0], a[1], z);
            xsref::np_assert_allclose!([h], [expected], rtol = rtol);
        }
    }
    let (a, z) = ([-50.999_999_999_699_924, -104.536_535_990_336_39], -0.056_613_881_920_597_09);
    let y = hyppfq(a, [], z);
    assert!(y.is_none_or(|y| (y / 2.794_502_416_895_517_4e22 - 1.0).abs() < 1e-13));
    // subnormal / extreme cancellation cases return `None` in `hyppfq` and `NaN` in `hyp2f0`
    assert_eq!(hyppfq([-2000.0, -2000.0], [], -0.000_122_549_019_607_843_14), None);
    for (a, z) in [([-238.57, -244.72], -0.0065), ([600.0, 4600.0], -0.000_514_2)] {
        assert_eq!(hyppfq(a, [], z), None);
        assert!(hyp2f0(a[0], a[1], z).is_nan());
    }
    // underflowing (zero) or nonzero subnormal t_1 in the asymptotic series still yields exact 1.0
    for (a, z) in [
        ([1e-300, 1e-300], -1e160),
        ([1.5e-308, -0.999_999_999_999_999_9], -1e-300),
        ([1.5e-308, -0.999_999_999_999_999_9], -0.125),
    ] {
        assert_eq!(hyppfq(a, [], z), Some(1.0));
        assert_eq!(hyp2f0(a[0], a[1], z), 1.0);
    }
}

#[test]
fn test_hyppfq_2f0_chgu() {
    // where U has the digits of the tolerance, but not those beyond (Arb)
    for (a, z, expected) in [
        (
            [0.426_888_586_508_577_86, 0.520_572_398_653_129_9],
            -4.231_846_256_386_015,
            0.723_332_978_613_73,
        ),
        (
            [-1.551_522_628_340_824, -1.383_752_833_722_752],
            -4.665_207_575_026_064,
            -4.942_797_210_602_93,
        ),
        (
            [0.791_211_332_690_890_1, 0.991_733_912_074_225_5],
            -4.008_904_464_785_383,
            0.407_420_222_646_740_3,
        ),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // where chgus would claim 13 digits without the rounding of pi z in gamma2 (Arb)
    for (a, z, expected) in [
        (
            [-52.086_010_359_946_73, -21.984_045_607_822_907],
            -13.677_918_500_493_34,
            3.597_641_568_386_878e109,
        ),
        (
            [-22.983_102_797_885_714, -35.543_635_003_749_316],
            -43.587_828_860_606_53,
            -1.313_003_920_498_504_4e87,
        ),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z);
            assert!(actual.is_none_or(|y| (y / expected - 1.0).abs() <= 1e-13), "{actual:?}");
        }
    }
    // U with 13 digits, but not 14, where the rounding of x adds to its error (Arb)
    let (a, z) = ([-41.063_445_135_539_524, -85.142_031_763_609_13], -244.069_262_850_980_4);
    for a in [a, [a[1], a[0]]] {
        let actual = hyppfq(a, [], z).unwrap();
        assert!((actual / 4.147_636_829_929_707_7e304 - 1.0).abs() <= 1.5e-13, "{actual}");
    }
    // where x^a1 U with 13 digits, but not 14, rounds to infinity, but the result is finite (Arb)
    let (a, z) = ([-41.529_323_181_743_06, -152.000_000_000_000_7], -3.414_969_771_369_697);
    for a in [a, [a[1], a[0]]] {
        let actual = hyppfq(a, [], z);
        let expected = 1.797_693_134_862_228_7e308;
        assert!(actual.is_none_or(|y| (y / expected - 1.0).abs() <= 1e-13), "{actual:?}");
    }
    // and where the result overflows by 1.3e-9 (Arb: 1.797_693_137_280_425e308)
    for a in [a, [a[1], a[0]]] {
        assert_eq!(hyppfq(a, [], -3.414_969_771_4), Some(f64::INFINITY));
    }
}

#[test]
fn test_hyppfq_2f0_residues() {
    // the residues in double-double, where U lacks a digit by the estimate of chgu and the other
    // methods fail: the two terms cancel, an a_i is near a negative integer, or the a_i are large
    // (Arb)
    for (a, z, expected) in [
        (
            [-0.630_819_077_138_073, -0.275_665_062_986_481_63],
            -7.558_864_736_580_486,
            0.040_883_462_433_760_266,
        ),
        (
            [-2.119_115_473_067_204, -2.859_944_142_412_956_6],
            -36.660_124_137_565_52,
            -329.568_040_645_147_56,
        ),
        (
            [-0.094_507_904_192_269_51, -0.325_624_596_942_675_3],
            -258.944_637_515_321_4,
            0.000_730_242_530_534_494_1,
        ),
        (
            [14.268_010_365_162_187, 5.147_607_270_817_833],
            -0.810_709_749_254_161_2,
            6.120_460_595_313_305e-6,
        ),
        (
            [-61.036_263_054_217_294, -59.649_519_065_312_376],
            -23.359_289_672_927_265,
            -2.744_638_996_882_791e163,
        ),
        (
            [-6.004_858_572_925_059, -21.079_297_738_677_543],
            -151_588_853.883_445_23,
            -1.040_990_391_775_226_8e184,
        ),
        ([0.5, 190.0], -100.0, 0.007_268_927_663_459_597),
    ] {
        for a in [a, [a[1], a[0]]] {
            let actual = hyppfq(a, [], z).unwrap();
            xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
        }
    }
    // where the result overflows (Arb: 2.119_933_336_826_352e381)
    let (a, z) = ([-66.319_864_448_545_77, -31.004_863_381_873_136], -48_506.922_888_667_235);
    assert_eq!(hyppfq(a, [], z), Some(f64::INFINITY));
    // but not within its error of the threshold, where it rounds to -f64::MAX or to -inf (Arb)
    let (a, z) = ([-7.000_000_000_000_008, -2.999_999_999_999_999_6], -1.016_898_363_184_092_2e46);
    let actual = hyppfq(a, [], z);
    assert!(actual.is_none_or(|y| (y / -f64::MAX - 1.0).abs() <= 1e-13), "{actual:?}");
    let (a, z) = ([-31.054_392_382_065_505, -30.054_392_339_243_95], -832_485_448.947_283_6);
    assert!([None, Some(f64::NEG_INFINITY)].contains(&hyppfq(a, [], z)));
    // nor where it underflows (Arb: 3.169e-395), as for a subnormal one
    let (a, z) = ([82.468_601_872_360_3, 58.426_626_520_618_98], -113_336.361_826_738_76);
    assert_eq!(hyppfq(a, [], z), None);
}

#[test]
fn test_hyppfq_borel() {
    // the Borel sum for p > q + 1 and z < 0, by the series up to its smallest term, the residues
    // of the Meijer G function, or its Weniger transformation (in f64 or double-double), also for
    // a_i that differ by integers or lie near a negative one (mpmath meijerg at 40 digits;
    // Wolfram agrees for a_i = 3/4, 7/4, 11/4 and for the a_i that differ by 2)
    for (a, b, z, expected) in [
        (&[1.0, 2.0, 3.0][..], &[][..], -0.01, 0.946_071_933_519_652_4),
        (&[1.5, -2.25, 0.75], &[], -0.05, 1.144_479_810_161_067_8),
        (&[0.3, 1.7, 0.6], &[], -2.5, 0.790_728_400_803_319_9),
        (&[0.5, 1.25, 2.5], &[], -1e4, 0.010_153_563_989_012_29),
        (&[0.75, 1.75, 2.75], &[], -30.0, 0.049_226_365_267_863_18),
        (&[1.5, 2.25, 0.75], &[3.5], -0.1, 0.936_956_810_640_734_1),
        (&[2.5, -3.5, 1.5], &[0.5], -1e-300, 1.0),
        (&[0.25, -1.75, 0.6, 1.1], &[], -40.0, 265.971_225_821_188_4),
        (&[-12.5, 3.25, -7.75, 1.5], &[2.5], -0.002, 0.698_147_990_588_820_4),
        (&[0.3, 1.7, 0.6, 2.2], &[1.5], -2.5, 0.743_181_240_510_237_2),
        (&[0.25, 2.25, 0.6, 1.1], &[1.5], -5.0, 0.763_430_979_979_843_5),
        (&[0.5, 1.25, -2.75, 3.5], &[1.5, 0.25], -1e6, 3.247_409_892_494_204_4e18),
        (&[1.2, 0.5, 1.75], &[], -0.06, 0.948_555_763_846_379_5),
        (&[1.3, 0.8, 0.7], &[0.3], -0.3, 0.574_543_997_541_148_2),
        (
            &[-3.0 + 2f64.powi(-33), 2.375, 0.9375, -1.8125, -4.998_046_875],
            &[],
            -0.8125,
            -761.523_562_944_142_1,
        ),
    ] {
        let actual = hyppfq(a, b, z).unwrap();
        xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
    }
    // complex for z > 0, as for 2F0
    assert_eq!(hyppfq([1.0, 2.0, 3.0], [], 0.01), None);
    // a polynomial whose terms cancel, with a_i that differ by an integer, has no residues (Arb)
    for (a, b, z, expected) in [
        (
            &[-200.0, -200.0, 0.5],
            &[][..],
            -0.001_225_490_196_078_431_4,
            1.184_410_924_689_460_4e165,
        ),
        (
            &[-2000.0, -2000.0, 0.5],
            &[0.3],
            -0.000_122_549_019_607_843_14,
            -0.023_231_370_699_859_81,
        ),
    ] {
        let actual = hyppfq(a, b, z);
        assert!(actual.is_none_or(|y| (y / expected - 1.0).abs() <= 1e-13), "{actual:?}");
    }
    // but the residues of other a_i that differ by an integer vanish (exact: 2.02e326, 5.64e324)
    let a = [-100.0, 1.5, 1.5, 2.25, 2.75];
    for b in [[0.5, 0.5], [0.7, 0.9]] {
        assert_eq!(hyppfq(a, b, -1.0), Some(f64::INFINITY));
    }
}

#[test]
fn test_hyppfq_borel_residues() {
    // the residues in double-double, where they cancel, e.g. for nearly equal a_i (mpmath
    // meijerg at 60 and 90 digits)
    for (a, b, z, expected) in [
        (
            &[1.0, 1.0 + 2f64.powi(-10), 1.0 - 2f64.powi(-10)][..],
            &[][..],
            -10.0,
            0.305_688_032_951_964_06,
        ),
        (
            &[-0.129_485_044_004_643_95, -0.205_791_231_073_163_58, -0.139_400_601_750_289_84],
            &[],
            -1.366_646_600_797_494_7,
            1.004_192_741_663_797_7,
        ),
        (
            &[7.386_572_466_708_953, 4.279_189_223_877_347_5, 5.834_117_507_796_039],
            &[],
            -1.276_126_675_055_175_9,
            1.989_747_179_487_340_7e-6,
        ),
        (
            &[8.314_471_273_616_888, 7.121_784_537_610_146_5, 11.597_238_385_168_016],
            &[5.099_214_049_145_988_5],
            -0.505_004_091_197_126_9,
            -6.420_708_860_397_208e-8,
        ),
        (
            &[-0.391_026_626_255_645_1, -0.296_214_915_415_556_1, -0.167_097_330_487_964_32],
            &[-0.088_446_353_843_811_5],
            -13.477_086_789_248_288,
            -0.651_536_282_336_249_4,
        ),
        (
            &[
                10.523_985_021_186_846,
                15.119_243_058_928_173,
                18.040_699_529_609_864,
                3.707_160_461_795_112_4,
            ],
            &[],
            -0.006_190_173_344_422_203,
            0.000_108_004_063_237_030_19,
        ),
        (
            &[
                0.891_919_465_010_852_6,
                0.751_648_993_887_893_2,
                2.016_692_712_135_52,
                1.018_510_355_057_419_5,
            ],
            &[-1.870_954_877_964_515_3],
            -35.210_729_664_195_696,
            -0.782_284_272_649_868_6,
        ),
    ] {
        let actual = hyppfq(a, b, z).unwrap();
        xsref::np_assert_allclose!([actual], [expected], rtol = 1e-14);
    }
    // but not within its error of the overflow threshold, where it rounds to f64::MAX (mpmath)
    let a = [-2.553_111_606_766_203, -2.553_111_499_272_675_6, -2.553_111_780_257_576_7];
    let actual = hyppfq(a, [], -8.199_439_242_455_9e118);
    assert!(actual.is_none_or(|y| (y / f64::MAX - 1.0).abs() <= 1e-13), "{actual:?}");
    // far above it (the logarithm of its largest residue is about 2.5e12)
    let (a, b) = ([-10_000_000_000.25, 0.25, 0.625], [-9_999_999_999.25]);
    assert_eq!(hyppfq(a, b, -1e100), Some(f64::INFINITY));
}

#[test]
fn test_hyppfq_cancel() {
    let z = 0.25;
    let f1f1 = Some(hyp1f1(0.5, 1.5, z));
    assert_eq!(hyppfq([0.5, 2.5], [2.5, 1.5], z), f1f1);
    // also a non-positive integer, which would make it a polynomial
    assert_eq!(hyppfq([-3.0, 0.5], [1.5, -3.0], z), f1f1);
    // each b_j cancels one a_i
    assert_eq!(hyppfq([2.5, 2.5], [2.5], z), Some(hyp1f0(2.5, z)));
    assert_eq!(hyppfq([1.0, 2.0], [2.0, 1.0], z), Some(hyp0f0(z)));
    assert_eq!(hyppfq([1.0, 2.0, 3.0, 4.0, 5.0], [5.0, 4.0, 3.0, 2.0], z), Some(hyp1f0(1.0, z)));
}

#[test]
fn test_hyppfq_special() {
    assert_eq!(hyppfq([0.5], [-2.0], 0.0), Some(1.0));
    // a pole, unless it is a polynomial of a lower degree
    assert_eq!(hyppfq([], [-2.0], 0.5), Some(f64::INFINITY));
    assert_eq!(hyppfq([-3.0], [-2.0], 0.5), Some(f64::INFINITY));
    assert_eq!(hyppfq([-1.0], [-2.0], 0.5), Some(1.25));
    // complex for p = q + 1 and z > 1, unless it is a polynomial or (1 - z)^-a
    assert!(hyppfq([0.3, 0.7], [1.5], 2.0).is_some_and(f64::is_nan));
    let (a, b) = ([0.5, 1.5, 2.5, 3.5], [1.25, 2.25, 3.25]);
    assert!(hyppfq(a, b, 2.0).is_some_and(f64::is_nan));
    assert_eq!(hyppfq([2.0], [], 3.0), Some(0.25));
    assert!(hyppfq([2.5], [], 3.0).is_some_and(f64::is_nan));
    assert_eq!(hyppfq([2.5], [], 1.0), Some(f64::INFINITY));
    assert_eq!(hyppfq([-2.5], [], 1.0), Some(0.0));
    assert!(hyppfq([0.5], [f64::NAN], 0.5).is_some_and(f64::is_nan));
}

#[test]
fn test_hyppfq_infinite() {
    let inf = f64::INFINITY;
    // equal infinities cancel, and z = 0 wins
    let f1f1 = Some(hyp1f1(0.5, 1.5, 0.25));
    assert_eq!(hyppfq([0.5, inf], [inf, 1.5], 0.25), f1f1);
    assert_eq!(hyppfq([inf], [-inf], 0.0), Some(1.0));
    // a polynomial, with the sign of its last term
    assert_eq!(hyppfq([-2.0], [2.5], -inf), Some(inf));
    assert_eq!(hyppfq([-3.0], [2.5], inf), Some(-inf));
    assert_eq!(hyppfq([-3.0], [-4.5], -inf), Some(-inf));
    assert_eq!(hyppfq([-3.0, 2.5], [-3.5, 4.5], inf), Some(inf));
    assert_eq!(hyppfq([-3.0], [], inf), Some(-inf));
    assert_eq!(hyppfq([0.0, 2.5], [1.5], -inf), Some(1.0));
    // for p <= q at z = +inf, with the sign of prod Gamma(b_j) / prod Gamma(a_i)
    assert_eq!(hyppfq([], [], inf), Some(inf));
    assert_eq!(hyppfq([0.5], [1.5], inf), Some(inf));
    assert_eq!(hyppfq([], [-2.5], inf), Some(-inf));
    assert_eq!(hyppfq([-2.5], [1.5], inf), Some(-inf));
    assert_eq!(hyppfq([-1.5, 2.5], [-0.5, -1.75, 2.25], inf), Some(-inf));
    // for q <= p <= q + 1 at z = -inf if all a_i > 0
    assert_eq!(hyppfq([], [], -inf), Some(0.0));
    assert_eq!(hyppfq([1.5], [-2.5], -inf), Some(0.0));
    assert_eq!(hyppfq([1.5, 2.5, 3.5], [-0.5, -2.25], -inf), Some(0.0));
    // the functions of p <= 2 and q <= 1, which have no limit for 0F1(0.25; z)
    assert_eq!(hyppfq([-2.5], [], -inf), Some(inf));
    assert_eq!(hyppfq([2.0], [], inf), Some(0.0));
    assert_eq!(hyppfq([-1.5], [-0.5], -inf), Some(-inf));
    assert_eq!(hyppfq([], [0.25], -inf), None);
    // complex for p = q + 1
    assert!(hyppfq([2.5], [], inf).is_some_and(f64::is_nan));
    assert!(hyppfq([1.5, 2.5, 3.5], [0.5, 4.5], inf).is_some_and(f64::is_nan));
    // the terms n >= 1 vanish for a b_j = +inf, unless it is complex
    assert_eq!(hyppfq([1.5, 2.5, 3.5], [0.5, inf], -10.0), Some(1.0));
    assert_eq!(hyppfq([1.5, 2.5, 3.5], [0.5, inf], 1.0), Some(1.0));
    assert_eq!(hyppfq([-2.0, 2.5, 3.5], [0.5, inf], 2.0), Some(1.0));
    assert_eq!(hyppfq([1.5, 2.5, 3.5], [0.5, inf], 2.0), None);
    // p > q + 1 diverges, unless it is a polynomial
    assert_eq!(hyppfq([0.5, 1.5, 2.5], [inf], 0.5), None);
    assert_eq!(hyppfq([-2.0, 1.5, 2.5], [inf], 0.5), Some(1.0));
    // a pole wins
    assert_eq!(hyppfq([1.5], [-2.0, inf], 0.25), Some(inf));
    // no limit
    assert_eq!(hyppfq([1.5], [-inf], 0.25), None);
    assert_eq!(hyppfq([1.5], [inf], inf), None);
    assert_eq!(hyppfq([inf], [2.5], 0.25), None);
}

#[test]
fn test_hyppfq_polynomial() {
    // the exact sums of the terms (Python `fractions`), rounded
    for (a, b, z, expected) in [
        (&[-5.0, 2.5][..], &[][..], -0.3, 35.447_635_937_499_996),
        (&[-10.0], &[2.5], 7.5, -0.230_642_806_801_874_54),
        (&[-3.0, 0.7], &[1.5], 2.0, 0.049_600_000_000_000_05),
        (&[-6.0, 1.5, 2.25, -0.5], &[3.5, 0.75, -8.0], 2.5, -2.188_432_770_889_141_6),
        (&[-20.0, 1.25, 1.5, 0.5], &[1.75, 1.5, 3.5], 5.0, 20_470_378.725_082_297),
        (&[-80.0, 2.5, 0.5, 0.5], &[1.5, 1.5, 3.5], 20.0, 1.998_050_134_345_888_1e96),
        (&[-7.0, 0.5, 1.5, 2.5, 3.5], &[0.25, 1.25, 2.25, -9.5], -1.75, -132.181_417_408_943_6),
    ] {
        let actual = hyppfq(a, b, z).unwrap();
        assert!(
            (actual - expected).abs() <= f64::EPSILON * expected.abs(),
            "hyppfq({a:?}, {b:?}, {z}) = {actual} != {expected}"
        );
    }
    // the sum cancels too much, but (1 - z)^70 does not
    assert_eq!(hyppfq([-70.0], [], 3.0), Some(2.0_f64.powi(70)));
}

/// Asserts that `hyppfq` is within a relative error of `1e-13` of `expected`
fn assert_hyppfq(a: &[f64], b: &[f64], z: f64, expected: f64) {
    let actual = hyppfq(a, b, z).unwrap();
    assert!(
        (actual - expected).abs() <= 1e-13 * expected.abs(),
        "hyppfq({a:?}, {b:?}, {z}) = {actual} != {expected}"
    );
}

#[test]
fn test_hyppfq_series() {
    // Arb references (the series with a bound on its rest)
    for (a, b, z, expected) in [
        (&[0.5, 1.5, 2.5, 3.5][..], &[1.25, 2.25, 3.25][..], 0.5, 1.701_434_243_680_788_3),
        (&[-0.5, 1.25, 2.75, -1.5], &[3.25, -2.75, 0.75], 0.3, 0.882_694_919_188_910_9),
        (
            &[0.5, 1.25, 2.75, 3.5, 0.75, 1.75],
            &[1.5, 2.25, 3.25, 4.75, 5.5],
            0.99,
            1.030_868_874_524_536_5,
        ),
        // p <= q
        (&[0.5, 1.5], &[1.25, 2.25, 3.25], 12.5, 3.369_349_110_002_69),
        (&[0.5, 1.25, 2.75], &[1.5, 2.25, 3.25], 25.0, 24_248_376.494_787_9),
        (&[-0.5], &[0.75, -2.75, 3.25], -28.0, -0.952_401_233_591_807_8),
        (&[], &[0.25, 1.75], -29.5, 12.856_431_709_673_524),
        (&[1.5, 2.5, -0.75], &[0.5, 3.25, 4.5, -1.25], 9.0, -3.884_219_564_556_629_3),
        // in double-double, as the terms cancel
        (&[6.5, 5.25, 4.75, 3.5], &[1.5, 1.25, 0.75], -0.55, -0.326_402_725_424_534),
        (&[3.5, 5.25, 4.75, 6.5, 2.25], &[1.5, 1.25, 0.75, 0.5], -0.3, -12.406_169_497_317_718),
        (&[1.5, 2.25, 0.75], &[3.5, 4.75, 1.25], -28.0, 0.252_024_878_533_297_47),
        // terms below half an ulp of the sum, whose roundings add up
        (&[2.17e-24, 2.0, 5.0, 12.0], &[1.0, 3.0, 10.0], 0.9978, 1.000_000_000_000_773_6),
        (&[3.35e-17, 2.5, 2.5], &[0.5, 3.5], 0.999, 1.000_000_000_000_111_3),
        // and where double-double underflows
        (&[2.3e-296, 45.0, 45.0, 45.0], &[1.0; 3], 0.98, 1.000_000_000_000_079_6),
        (&[2.89e-296, 45.0, 45.0, 45.0], &[1.0; 3], 0.98, 1.000_000_000_000_1),
        (&[2.8909e-296, 45.0, 45.0, 45.0], &[1.0; 3], 0.98, 1.000_000_000_000_1),
    ] {
        assert_hyppfq(a, b, z, expected);
    }
}

#[test]
fn test_hyppfq_series_flat() {
    // the sum of the |t_k| is just above 2^-20 (Arb reference)
    let y = hyppfq([-1.59065e-7, -1.000_000_001_4, 2.0, 3.0], [1.0; 3], 0.999_25).unwrap();
    let expected = 1.000_000_953_674_056;
    assert!((y - expected).abs() <= 1e-15 * expected, "{y} != {expected}");
    // the series actually grows past `f64::MAX`, while an unpadded ratio bound rounds below 1
    // (#1784): must reject rather than accept the partial sum near 1
    let (a, b) = (f64::from_bits(0x44f4_0000_0000_0003), f64::from_bits(0x44f4_0000_0000_0000));
    let c = -2_000.000_000_000_000_2;
    let y =
        hyppfq([5e-19, 1.0, a, a, c, c], [0.01, b, b, -1_999.95, -1_999.95], 0.999_999_999_999_999);
    assert_eq!(y, None);
}

#[test]
fn test_hyppfq_asymptotic() {
    // Arb references, where the series has too many terms
    let (a, b) = ([0.5, 1.25, 2.75], [1.5, 2.25, 3.25]);
    let a_neg = [-0.5, 1.25, 2.75];
    let (a4, b4) = ([0.5, 1.25, 2.75, -3.375], [1.5, 2.25, 3.25, -4.125]);
    assert_hyppfq(&a, &b, 500.0, 2.494_168_748_039_884e210);
    assert_hyppfq(&a_neg, &b, 400.0, -2.048_065_752_931_792_7e164);
    // p < q
    assert_hyppfq(&[0.5, 1.25], &b, 1e5, 2.827_653_514_543_852_6e262);
    assert_hyppfq(&[], &[0.25, 1.75], 1e7, 7.070_838_755_413_847_6e277);
    assert_hyppfq(&[0.5], &b4, 1e8, -1.021_654_858_110_642_5e170);
    // z < 0, where the series cancels
    assert_hyppfq(&a, &b, -200.0, 0.115_205_192_301_031_53);
    assert_hyppfq(&a, &b, -1e4, 0.016_481_346_737_557_533);
    assert_hyppfq(&a4, &b4, -80.0, 3_084.118_291_519_921_5);
    // (small a_i - a_j, so that the terms of the algebraic part cancel)
    let (a_close, b_close) = ([0.03125, 0.125, 0.015_625], [-0.03125, 0.1875, -0.0625]);
    assert_hyppfq(&a_close, &b_close, -3e4, 1.004_849_422_125_088_7);
    let (a_close, b_close) = ([-0.03, 0.08, -0.06], [0.015, -0.11, 0.03]);
    assert_hyppfq(&a_close, &b_close, -6e4, 1.719_455_843_480_630_1);
    // (large Gamma factors)
    let (a_large, b_large) = ([1.5, 1.25, -0.5], [-9.25, -11.25, 4.5]);
    assert_hyppfq(&a_large, &b_large, -500.0, 2.827_497_437_819_951);
    assert_hyppfq(&[], &[-8.5, 7.5], -1e5, 3.091_767_325_440_432_3e30);
    // (below eps after a few terms)
    let (a_few, b_few) = ([0.5, 3.25, 0.25], [-1.25, -0.25, 0.75]);
    assert_hyppfq(&a_few, &b_few, -1e5, 0.141_836_995_072_879_2);
    let (a_few, b_few) = ([-0.5, 2.75, 3.25], [-1.5, -1.25, 1.75]);
    assert_hyppfq(&a_few, &b_few, -3e4, 4.996_596_806_437_423e-13);
    // p = q - 1, where it oscillates
    assert_hyppfq(&[0.5, 1.25], &b, -1e5, 0.007_402_197_567_641_71);
    assert_hyppfq(&[0.5, 1.25], &b, -1e8, 0.000_234_112_089_078_822_63);
    assert_hyppfq(&[0.5, 1.25, -2.625], &b4, -3e4, -91_536_564.254_330_46);
    // p <= q - 2, where it oscillates and grows
    assert_hyppfq(&[], &[0.25, 1.75], -1e7, 5.532_241_318_399_256e137);
    assert_hyppfq(&[0.5], &b4, -1e8, -2.517_853_121_265_854e119);
    // (the algebraic part is more than half of these)
    let (a_alg, b_alg) = ([-0.375, 1.875], [3.25, 4.125, 5.5, 6.25]);
    assert_hyppfq(&[-1.4], &b_alg[..3], -4e4, 10_094.591_796_559_325);
    assert_hyppfq(&a_alg, &b_alg, -4e4, 9.716_270_976_402_212);
    // 1 + a_m - a_i and 1 + a_m - b_j are close to the same integer
    let a_near = [0.125, -4.875_000_000_000_001, 0.375];
    let b_near = [0.124_999_999_999_999_97, 15.5, -8.25];
    assert_hyppfq(&a_near, &b_near, -30.0, 0.739_525_258_049_542_4);
    let a_pair = [-0.125, -5.125_000_000_000_001];
    let b_pair = [-0.125_000_000_000_000_06, 15.5, -8.25, 1.5];
    assert_hyppfq(&a_pair, &b_pair[..3], -1000.0, -36_610.554_660_278_02);
    let a_pair4 = [-0.125, -5.124_999_999_999_999];
    assert_hyppfq(&a_pair4, &b_pair, -1e4, -14_340_142.864_468_437);
    // the terms of the series of a tiny a_1 grow after its second one, and then shrink to
    // smaller ones: accurate or None
    let a_tiny = [
        -1.177_387_855_629_240_8e-18,
        -15.152_718_139_177_594,
        -6.428_157_404_437_108,
        -7.171_717_841_041_328,
    ];
    let b_tiny = [
        0.523_592_411_670_334_3,
        -203.347_306_012_400_4,
        -241.850_859_418_556_47,
        -191.562_668_472_961_35,
    ];
    let y = hyppfq(a_tiny, b_tiny, -3_439.845_313_014_641);
    assert!(y.is_none_or(|y| (y / 2.050_874_042_107_11 - 1.0).abs() <= 1e-13));
    // (and here they are below eps of the sum after the second one, while their ratio is below
    // 1, before they grow)
    let a_tiny = [2.0_f64.powi(-65), -22.5, -27.375, -20.625];
    let b_tiny = [0.25, -475.5, -346.375, -392.625];
    assert_hyppfq(&a_tiny, &b_tiny, -7000.0, 2.755_031_164_893_358);
    // overflow, with the sign of the Gamma functions
    assert_eq!(hyppfq(a, b, 800.0), Some(f64::INFINITY));
    assert_eq!(hyppfq(a_neg, b, 800.0), Some(f64::NEG_INFINITY));
    assert_eq!(hyppfq([], [0.25, 1.75], 1e9), Some(f64::INFINITY));
    assert_eq!(hyppfq([-0.5, 1.25], b, 1e7), Some(f64::NEG_INFINITY));
    // integer a_i - a_j, where the expansion has a limit only
    let c = 1.0 + 2.0_f64.powi(-10);
    let (a_int, b_int) = ([1.0, 1.0, 10.5], [c, c, 11.25]);
    assert_hyppfq(&a_int, &b_int, 21.25, 739_410_485.211_035_5);
    // large |z| with integer a_i - a_j, where the algebraic part is the limit of DLMF 16.11.2,
    // of which only a bound is used for z > 0
    let (a_half, b_half) = ([0.5, 1.5, 2.5], [1.25, 2.25, 3.25, 0.75]);
    assert_hyppfq(&a_half, &b_half[..3], 500.0, 1.490_249_682_266_956_5e211);
    assert_hyppfq(&a_half, &b_half[..3], -200.0, 0.086_921_283_327_525_05);
    let b_half3 = [1.75, 2.25, 3.25];
    assert_hyppfq(&a_half[..2], &b_half3, 1e5, 2.998_583_212_151_733e262);
    assert_hyppfq(&a_half[..2], &b_half3, -1e5, 0.007_069_421_929_959_013);
    assert_hyppfq(&a_half[..2], &b_half, -1e6, -4.192_062_553_632_357e55);
    // (a pole of order 3)
    assert_hyppfq(&[1.0; 3], &[2.0; 3], -1e3, 0.028_834_862_048_815_51);
    // (an a_i and a b_j also differ by integers)
    let b_mixed = [0.5, 4.5, 2.25];
    assert_hyppfq(&[1.5, 2.5], &b_mixed, -2e4, -6.818_735_291_796_145e-6);
    assert_hyppfq(&[1.5, 2.5], &b_mixed, 2e4, 4.540_531_903_354_946e117);
}

#[test]
fn test_hyppfq_asymptotic_double_double() {
    // Arb references, for p = q and z < 0
    // (the terms of a series of the algebraic part cancel)
    let a = [-19.313_050_582_041_633, -13.850_622_854_249_043, -2.763_519_547_786_056];
    let b = [-31.346_275_528_683_137, -31.480_497_899_162_355, 17.040_681_912_798_888];
    assert_hyppfq(&a, &b, -1_340.224_629_533_191, -5_982_750.795_680_105);
    // (the terms of the series of a_3 grow again after its 1 + a_3 - a_i + n are positive)
    let (a, b) = ([55.24, 55.240_001, -0.87, 51.03], [17.79, -39.49, -50.97, -46.83]);
    assert_hyppfq(&a, &b, -3000.0, 100.085_095_585_652_68);
    // (the terms of the series of a tiny a_1 grow after the second one, to 1e41 times its size:
    // no bound shows that one of them bounds the rest)
    let (a, b) = ([2.0_f64.powi(-60), 0.125, 0.25], [-79.5, -144.25, -239.5]);
    assert_eq!(hyppfq(a, b, -1732.0), None);
    // (and here they grow from the fourth one on, long before n = |z|)
    let (a, b) = ([2.0_f64.powi(-70), -19.75, -29.625], [1.0, -219.75, -229.625]);
    assert_eq!(hyppfq(a, b, -100.0), None);
    // (the parameters of the series are positive from n = 2^60 on)
    let (a, b) = ([0.5, 2.0_f64.powi(60)], [2.0_f64.powi(54), 1.25]);
    assert_eq!(hyppfq(a, b, -1.75), None);
    // (small a_i - a_j, so that the terms of the algebraic part cancel)
    let (a, b) = ([0.5, 0.501, 0.502], [1.25, -0.75, 2.5]);
    assert_hyppfq(&a, &b, -2e5, -0.032_631_251_264_031_06);
    let (a, b) = ([0.02, 0.03, 0.01, 0.04], [-0.01, 0.05, 0.06, -0.02]);
    assert_hyppfq(&a, &b, -2e4, 0.711_858_233_663_966_2);
}

#[test]
fn test_hyppfq_continuation() {
    // mpmath references at 50 and 120 digits
    let (a, b) = ([0.5, 1.25, 2.75, 3.375], [3.25, 4.75, 1.75]);
    for (a, b, z, expected) in [
        (&a[..], &b[..], -10.0, 0.475_228_739_431_369_9),
        (&a, &b, -1e3, 0.059_103_932_783_523_536),
        (&a, &b, -1e6, 0.001_886_952_243_591_441_8),
        (
            &[0.5, 1.25, 2.75, 3.375, 0.125],
            &[3.25, 4.75, 1.75, 2.5],
            -50.0,
            0.865_198_007_619_403_1,
        ),
        (
            &[0.5, 1.25, 2.75, 3.375, 0.125, 1.0625],
            &[3.25, 4.75, 1.75, 2.5, 0.375],
            -7.5,
            0.888_755_695_438_054,
        ),
        // by `pfq_weniger`
        (&[0.5, 1.25, 2.75, 3.5], &[3.25, 4.75, 1.75], -1.0, 0.839_958_420_055_184_8),
        (&[0.5, 1.25, 2.75, 3.5, 0.125], &[3.25, 4.75, 1.75, 2.5], -1.25, 0.987_921_427_741_094),
    ] {
        assert_hyppfq(a, b, z, expected);
    }
    // integer a_i - a_j, where it is a limit; Arb references, as the sum of the residues of the
    // integrand of DLMF 16.5.1
    for (a, b, z, expected) in [
        // Li_3(z) / z and Li_4(z) / z
        (&[1.0; 4][..], &[2.0; 3][..], -40.0, 0.361_479_025_355_404_44),
        (&[1.0; 4], &[2.0; 3], -1e12, 3.561_375_815_418_521_6e-9),
        (&[1.0; 5], &[2.0; 4], -1e4, 0.037_150_546_634_823_31),
        // (1 + 11 z + 11 z^2 + z^3) / (1 - z)^5
        (&[2.0; 4], &[1.0; 3], -3.0, 0.039_062_5),
        // the b_j - a_i are integers too
        (&[5.0, 5.0, 3.0, 3.0], &[1.0, 4.0, 2.0], -8.5, -4.548_075_461_510_543e-6),
        (&[5.0, 2.0, 5.0, 5.0], &[1.0, 4.0, 4.0], -6e10, -3.215_020_574_792_095_3e-55),
        // two a_i of three differ from another by integers
        (&[0.25, 1.25, 3.25, 0.75], &[2.0, 1.5, 0.375], -1e3, 0.077_230_437_610_217_39),
        // boundary around z = -1, near/exact integer diffs, and genuine b_j = 1 with tiny a_j
        // (mpmath 60-dps)
        (&[0.5, 1.25, 2.75], &[2.25, 3.5], -1.0 - f64::EPSILON, 0.843_791_691_785_246_3),
        (&[0.5, 1.25, 2.75], &[2.25, 3.5], -1.0, 0.843_791_691_785_246_4),
        (&[0.5, 1.25, 2.75], &[2.25, 3.5], -1.0 + f64::EPSILON / 2.0, 0.843_791_691_785_246_4),
        (&[0.5, 1.25, 2.75], &[2.25, 3.5], -1.01, 0.842_648_061_384_124_5),
        (&[0.5, 1.5 + 1e-8, 2.75], &[2.25, 3.5], -2.0, 0.716_075_583_409_688_4),
        (&[0.5, 1.5, 2.75], &[2.25, 3.5], -2.0, 0.716_075_584_737_973_7),
        (&[2.0_f64.powi(-60), 0.5, 1.25], &[1.0, 2.5], -2.0, 1.0),
        // zero inner polynomial term and mixed exact-integer class + cancelling singletons
        (
            &[0.5, 1.5 - 2.0_f64.powi(-40), 1.75],
            &[2.5, 1.0],
            -2.0_f64.powi(40),
            7.959_621_827_055_401e-7,
        ),
        (
            &[0.5, 0.5 + 2.0_f64.powi(-32), 2.75, 3.75],
            &[2.25, 3.5, 4.125],
            -1000.0,
            0.191_031_737_123_466_5,
        ),
        // near-pole `b_j` resolved by double-double continuation (Arb references)
        (&[-1.5, 0.75, -0.125, -77.75], &[-3.75, -78.875, 21.5], -2.0, 1.441_597_903_532_934),
        (
            &[0.17549, -33.943, -0.1108, 0.011_587, 0.041_28],
            &[-48.912, 0.001_225_3, 29.137, 0.009_330_8],
            -6.5618,
            1.118_370_664_887_746,
        ),
    ] {
        assert_hyppfq(a, b, z, expected);
    }
}

#[test]
fn test_hyppfq_closed_forms() {
    // The formulas `section/number` of
    // https://functions.wolfram.com/HypergeometricFunctions/HypergeometricPFQ/03/01/, whose
    // special functions are evaluated with Arb
    for (a, b, z, expected) in [
        // 04/0004: Li_q(z) / z, also for z < -1 with integer a_i - a_j
        (&[1.0; 4][..], &[2.0; 3][..], 0.5, 1.074_426_387_216_080_4),
        (&[1.0; 4], &[2.0; 3], -0.75, 0.922_271_471_587_279_3),
        (&[1.0; 4], &[2.0; 3], -2.0, 0.834_141_681_983_285_6),
        (&[1.0; 5], &[2.0; 4], 0.9375, 1.075_162_383_357_896_7),
        // 06/0008: 3 zeta(3) / 4, and 06/0007 for q = 4 and 6
        (&[1.0; 4], &[2.0; 3], -1.0, 0.901_542_677_369_695_7),
        (&[1.0; 5], &[2.0; 4], -1.0, 7.0 * PI.powi(4) / 720.0),
        (&[1.0; 7], &[2.0; 6], -1.0, 31.0 * PI.powi(6) / 30240.0),
        // 04/0003: a^q Phi(z, q, a) with the Lerch transcendent
        (&[1.0, 0.25, 0.25, 0.25], &[1.25; 3], -1.0, 0.993_051_520_249_976_6),
        (&[1.0, 0.25, 0.25, 0.25], &[1.25; 3], 0.875, 1.008_602_569_162_220_2),
        (&[1.0, 0.375, 0.375, 0.375, 0.375], &[1.375; 4], -0.5, 0.997_372_909_080_830_9),
        // 04/0002: a sum of 2F1(a, a_k; a_k + 1; z), also for z < -1
        (&[0.375, 0.25, 0.625, 0.875], &[1.25, 1.625, 1.875], 0.75, 1.011_649_399_153_953_5),
        (&[0.375, 0.25, 0.625, 0.875], &[1.25, 1.625, 1.875], -1.25, 0.985_603_765_207_305_6),
        (&[0.375, 0.25, 0.625, 0.875], &[1.25, 1.625, 1.875], -3.0, 0.970_328_335_697_056_8),
        (
            &[2.375, 0.25, 0.625, 0.875, 1.5],
            &[1.25, 1.625, 1.875, 2.5],
            -0.9,
            0.962_230_307_834_818_8,
        ),
        (
            &[2.375, 0.25, 0.625, 0.875, 1.5],
            &[1.25, 1.625, 1.875, 2.5],
            -30.0,
            0.701_332_543_330_772_8,
        ),
        // 04/0010: a rational function of z
        (&[2.0; 4], &[1.0; 3], -0.5, -20.0 / 81.0),
        (&[2.0; 4], &[1.0; 3], 0.5, 300.0),
        (&[2.0; 6], &[1.0; 5], -0.25, 0.236_544),
        // 06/0001: digamma functions of the a_k
        (&[1.0, 0.25, 0.5, 0.75], &[1.25, 1.5, 1.75], -1.0, 0.975_967_713_426_429_7),
        (&[1.0, 0.25, 0.5, 0.75, 2.375], &[1.25, 1.5, 1.75, 3.375], -1.0, 0.982_477_625_546_169_4),
        // 06/0009: zeta(3) and ln 2
        (&[1.0; 4], &[3.0; 3], -1.0, 0.966_812_171_669_881_1),
        (&[1.0; 5], &[3.0; 4], -1.0, 0.988_341_738_313_950_6),
        // 06/0010: pi^3 / 32, and Hurwitz zeta functions for q = 4
        (&[1.0, 0.5, 0.5, 0.5], &[1.5; 3], -1.0, PI.powi(3) / 32.0),
        (&[1.0, 0.5, 0.5, 0.5, 0.5], &[1.5; 4], -1.0, 0.988_944_551_741_105_3),
        // 30/0012 and 30/0013 for d = 1/4
        (
            &[1.0, 1.0, 1.0, 0.25, 0.25, 1.75, 1.75],
            &[2.0, 2.0, 1.25, 1.25, 2.75, 2.75],
            -1.0,
            0.996_208_699_633_151_4,
        ),
        (
            &[1.0, 1.5, 1.5, 0.25, 0.25, 0.75, 0.75],
            &[0.5, 0.5, 1.25, 1.25, 1.75, 1.75],
            -1.0,
            0.949_477_168_790_850_3,
        ),
    ] {
        assert_hyppfq(a, b, z, expected);
    }
}

#[test]
fn test_hyppfq_at_one() {
    // +inf where the series diverges, for sum(b) <= sum(a)
    let a = [0.5, 1.5, 2.5, 3.5];
    assert_eq!(hyppfq(a, [1.25, 2.25, 3.25], 1.0), Some(f64::INFINITY));
    assert_eq!(hyppfq(a, [1.25, 2.25, 4.5], 1.0), Some(f64::INFINITY));
    // no method where it converges, as to zeta(6) (30/0011 of the formulas above)
    assert_eq!(hyppfq(a, [1.25, 2.25, 4.75], 1.0), None);
    assert_eq!(hyppfq([1.0; 7], [2.0; 6], 1.0), None);
    // sum(b) - sum(a) = 2^-99 - 254 2^-200 > 0 is within the rounding errors of the sums
    let pow2 = |n| 2.0_f64.powi(n);
    let mut a = vec![0.5, 0.5, pow2(-54), pow2(-54), pow2(-99)];
    a.extend([pow2(-200); 254]);
    let mut b = vec![1.0, pow2(-53)];
    b.extend([pow2(-106); 256]);
    assert_eq!(hyppfq(&a, &b, 1.0), None);
}

#[test]
fn test_hyppfq_none() {
    // the sum cancels too much
    assert_eq!(hyppfq([-60.0, 0.5, 0.5], [1.5, 1.5, 1.5], 200.0), None);
    assert_eq!(hyppfq([0.5, 1.5], [], 0.5), None);
}

#[test]
fn test_hyppfq_conformal() {
    // the sum cancels too much; Wolfram reference
    let (a, b) = ([10.5, 8.25, 9.75, 7.5], [0.5, 1.25, 0.75]);
    assert_hyppfq(&a, &b, -0.5, 891.432_204_916_963_3);
}

#[test]
fn test_hyp_macro() {
    use xsf::hyp;

    assert_eq!(hyp![1.0, 2.0; 3.0; -0.5], hyppfq([1.0, 2.0], [3.0], -0.5));
    assert_eq!(hyp![1.0, 2.0, 3.0; 4.0, 5.0; 0.25,], hyppfq([1.0, 2.0, 3.0], [4.0, 5.0], 0.25));
    assert_eq!(hyp![0.5; 1.5,; 0.3], hyppfq([0.5], [1.5], 0.3));
    // empty lists
    assert_eq!(hyp![; ; 0.7], hyppfq([], [], 0.7));
    assert_eq!(hyp![0.5; ; -0.3], hyppfq([0.5], [], -0.3));
    assert_eq!(hyp![; 1.5; 0.3], hyppfq([], [1.5], 0.3));
    // each argument is evaluated once, in order
    let mut n = 0.0;
    let mut next = || {
        n += 1.0;
        n
    };
    assert_eq!(hyp![next(), next(); next(); 0.1 * next()], hyppfq([1.0, 2.0], [3.0], 0.4));
    assert_eq!(n, 4.0);
}

#[test]
fn test_hyp2f2_series_sum_above_1e154() {
    // the squared partial sums of the error estimate overflow
    assert_hyp2f2([76.5, 4.25], [-0.25, -18.75], 200.0, 7.879_700_512_996_357e201);
}

#[test]
fn test_hyppfq_series_sum_above_1e154_rounding() {
    // the rounding of a + n adds up, so the estimate of the f64 series was too small
    let a = 255.999_999_999_999_97;
    let y = hyppfq([a; 3], [256.0; 3], 700.0).unwrap();
    let expected = 1.014_232_054_734_890_5e304;
    assert!((y - expected).abs() <= 5e-14 * expected, "{y} != {expected}");
}

#[test]
fn test_hyp1f2_series_sum_near_max() {
    // the error estimate of the series sum 1.2e308 overflowed
    assert_hyp1f2(-2.0, [5e-308, 1.0], -2.0, 1.2e308);
}

#[test]
fn test_hyppfq_series_sum_above_1e154_rounding_below() {
    // the roundings of a + n and b + n go the other way
    let (a, b) = (255.000_000_000_000_03, 254.999_999_999_999_97);
    let y = hyppfq([a; 3], [b; 3], 355.0).unwrap();
    let expected = 1.494_655_400_472_756_6e154;
    assert!((y - expected).abs() <= 5e-14 * expected, "{y} != {expected}");
}
