use core::f64::consts::PI;
use xsf::{hyp0f0, hyp0f1, hyp1f0, hyp1f1, hyp1f2, hyp2f1, hyp2f2, hyp3f2, hyppfq};

/// Asserts that `hyp3f2` is within a relative error of `1e-13` of `expected`
fn assert_hyp3f2(a: [f64; 3], b: [f64; 2], z: f64, expected: f64) {
    let actual = hyp3f2(a[0], a[1], a[2], b[0], b[1], z).unwrap();
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
        (
            [30.5, 0.25, -1.75],
            [20.25, 8.5],
            0.95,
            0.931_616_739_208_828_7,
        ),
        ([1.5, 2.5, 3.5], [-0.5, 6.0], -0.125, 1.256_376_554_812_066),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
    // the parts of the error estimate would overflow (Arb)
    let y = hyp3f2(2.0, 1.0, 900.2, 4.5, 350.2, 0.8);
    assert!(y.is_none_or(|y| (y / 3.717_681_528_889_451e152 - 1.0).abs() < 1e-13));
}

#[test]
fn test_hyp3f2_polynomial() {
    // mpmath references, cross-checked with Arb
    for (a, b, z, expected) in [
        (
            [-4.0, 2.5, 3.25],
            [0.75, -6.5],
            12.5,
            355_339.477_100_550_95,
        ),
        (
            [-25.0, 10.5, 0.1],
            [30.25, 1.5],
            -3.0,
            5_267_612.739_117_058,
        ),
        // terminates before the zero of (b_1)_n
        ([-3.0, -7.0, 1.5], [-5.0, 2.5], -40.0, 79_088.466_666_666_66),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
    assert_eq!(hyp3f2(0.0, 2.5, 3.5, -1.0, 2.0, 0.75), Some(1.0));
}

#[test]
fn test_hyp3f2_special() {
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 0.0), Some(1.0));
    assert!(hyp3f2(f64::NAN, 2.5, 3.5, 0.5, 4.5, 0.25).unwrap().is_nan());
    assert!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, f64::NAN).unwrap().is_nan());
    // the limits of `hyppfq` for an infinite argument
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, f64::INFINITY, 0.25), Some(1.0));
    assert_eq!(
        hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, f64::NEG_INFINITY),
        Some(0.0)
    );
    assert_eq!(hyp3f2(f64::INFINITY, 2.5, 3.5, 0.5, 4.5, 0.25), None);

    // a zero of (b_2)_n, also if the polynomial terminates after it
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, -2.0, 0.25), Some(f64::INFINITY));
    assert_eq!(hyp3f2(1.5, -3.0, 3.5, 0.5, -2.0, 0.25), Some(f64::INFINITY));

    // complex for z > 1
    assert!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 1.5).unwrap().is_nan());
    // where the series diverges
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 1.0), Some(f64::INFINITY));
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 1.5, 1.5, 1.0), Some(f64::INFINITY));
    // b_1 + b_2 - a_1 - a_2 - a_3 = 2.8e-17 rounds to 0; the result is 2.5e15 (Arb)
    let y = hyp3f2(5.0, 0.1, 0.2, 5.0, 0.300_000_000_000_000_04, 1.0);
    assert!(y.is_none_or(|y| (y / 2.467_832_059_538_944e15 - 1.0).abs() < 1e-13));
    // the excess 1e-300 is below the rounding of the double-double sums
    let y = hyp3f2(
        0.750_000_000_000_000_2,
        2.5,
        -1e-300,
        1.750_000_000_000_000_2,
        1.5,
        1.0,
    );
    assert!(y.is_none_or(|y| (y / 0.499_999_999_999_999_85 - 1.0).abs() < 1e-13));

    // a common parameter cancels, also a non-positive integer
    for z in [-30.0, 0.25, 0.95] {
        assert_hyp3f2([1.5, 2.5, 3.25], [2.5, 6.5], z, hyp2f1(1.5, 3.25, 6.5, z));
        assert_hyp3f2([-2.0, 1.5, -3.0], [4.5, -2.0], z, hyp2f1(1.5, -3.0, 4.5, z));
        // two of them
        assert_hyp3f2([1.5, 2.5, -3.0], [2.5, -3.0], z, (1.0 - z).powf(-1.5));
    }
    assert!(hyp3f2(1.5, 2.5, 3.25, 2.5, 6.5, 1.5).unwrap().is_nan());
    assert!(hyp3f2(1.5, 2.5, 3.25, 2.5, 1.5, 1.5).unwrap().is_nan());

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
        (
            [0.5, 1.25, 2.75],
            [3.25, 4.75],
            -100.0,
            0.267_996_937_273_824_95,
        ),
        (
            [10.5, 20.25, 0.125],
            [5.75, 30.5],
            -50.0,
            0.594_745_542_710_990_4,
        ),
        (
            [0.25, 1.375, 2.875],
            [1.75, 3.875],
            -1e6,
            0.038_052_289_002_429_94,
        ),
        // b_1 - a_1 = -1, so that the term of a_1 vanishes
        (
            [2.5, 0.25, 1.125],
            [1.5, 4.25],
            -3.0,
            0.822_861_429_867_887_1,
        ),
        // the factors of each coefficient underflow or overflow, unlike their product
        (
            [80.25, 120.5, 160.75],
            [150.125, 140.5],
            -1000.0,
            1.689_865_300_058_061_6e-236,
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
        (
            [0.5, 1.25, 2.75],
            [3.25, 4.75],
            -0.99,
            0.910_875_206_300_700_3,
        ),
        // z < -1, where the terms of the continuation cancel
        (
            [0.5, 1.25, 2.75],
            [3.25, 4.75],
            -2.5,
            0.820_925_159_215_043_4,
        ),
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
        (
            [3.0, 0.5, -19.5],
            [0.25, 0.125],
            0.25,
            -1.187_310_231_031_717_4,
        ),
        (
            [10.5, 12.5, 9.25],
            [1.5, 2.5],
            -0.25,
            -0.037_909_280_209_963_286,
        ),
        (
            [20.5, 0.75, 3.25],
            [1.5, 2.5],
            -0.5,
            0.099_410_623_211_045_2,
        ),
        // near a zero of the polynomial
        (
            [-2.0, 1.0, 1.0],
            [0.5, 0.5],
            0.143_237_254_218_789_44,
            -4.403_648_831_270_954e-17,
        ),
    ] {
        assert_hyp3f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp3f2_none() {
    // z = 1
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 8.5, 1.0), None);
    // more than 10^4 terms
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 0.999), None);
    // z < -1 with integer a_i - a_j, where the sequence transformation does not converge
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 2.0, 2.0, -50.0), None);
    // the same for a huge |z|, where its orders hardly change before they converge
    assert_eq!(hyp3f2(2e-15, 1.0, 1.0, 2.0, 2.0, -1e150), None);
    // the terms cancel too much
    assert_eq!(hyp3f2(50.5, 60.5, 70.5, 1.5, 2.5, -0.75), None);
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, -1.0), None);
    // a huge negative parameter
    assert_eq!(hyp3f2(0.5 - 1e12, 1.0, 1.5, 2.0, 2.5, -0.9), None);
}

#[test]
fn test_hyp3f2_underflow() {
    // the product of the a_i underflows, but not the terms; Arb reference
    let (a, b) = ([1e-170, 1e-170, 100.0], [1e-150, 2e-150]);
    assert_hyp3f2(a, b, 0.75, 8.034_690_221_294_95e19);
    // the terms underflow, and stay negligible
    assert_eq!(hyp3f2(-30.0, 1.0, 1.0, 2.0, 2.0, 1e-20), Some(1.0));
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 1e200, 1e200, 0.5), Some(1.0));
    assert_eq!(hyp3f2(1e-320, 1.0, 1.0, 2.0, 2.0, 0.5), Some(1.0));
    // the first term underflows, and later ones are large
    assert_eq!(hyp3f2(1e-200, 1e-200, 500.0, 1.0, 1.0, 0.9), None);
    // the same at z = -1, where the series is summed by the sequence transformation
    assert_eq!(hyp3f2(1.0, 1.0, 1.0, 2.0, 1e300, -1.0), Some(1.0));
    assert_eq!(hyp3f2(1e-200, 1e-200, -1500.5, 1.0, 1.0, -1.0), None);
}

/// Asserts that `hyp2f2` is within a relative error of `1e-13` of `expected`
fn assert_hyp2f2(a: [f64; 2], b: [f64; 2], z: f64, expected: f64) {
    let actual = hyp2f2(a[0], a[1], b[0], b[1], z).unwrap();
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
fn test_hyp2f2_polynomial() {
    // Arb references
    for (a, b, z, expected) in [
        ([-4.0, 2.5], [0.75, -6.5], 12.5, 1_424.990_151_262_878_5),
        ([-25.0, 10.5], [30.25, 1.5], -3.0, 1_431.131_238_039_312_8),
        // terminates before the zero of (b_1)_n
        ([-3.0, 1.5], [-5.0, 2.5], -40.0, -266.098_412_698_412_7),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
}

#[test]
fn test_hyp2f2_special() {
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, 4.5, 0.0), Some(1.0));
    assert!(hyp2f2(f64::NAN, 2.5, 3.5, 4.5, 0.25).unwrap().is_nan());
    assert!(hyp2f2(1.5, 2.5, 3.5, 4.5, f64::NAN).unwrap().is_nan());
    // the limits of `hyppfq` for an infinite argument
    assert_eq!(hyp2f2(1.5, 2.5, f64::INFINITY, 4.5, 0.25), Some(1.0));
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, 4.5, f64::NEG_INFINITY), Some(0.0));
    assert_eq!(
        hyp2f2(1.5, -2.5, 3.5, 4.5, f64::INFINITY),
        Some(f64::NEG_INFINITY)
    );
    assert_eq!(hyp2f2(f64::INFINITY, 2.5, 3.5, 4.5, 0.25), None);

    // a zero of (b_2)_n, also if the polynomial terminates after it
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, -2.0, 0.25), Some(f64::INFINITY));
    assert_eq!(hyp2f2(1.5, -3.0, 3.5, -2.0, 0.25), Some(f64::INFINITY));

    // a common parameter cancels, also at a pole (Arb references of 1F1)
    assert_hyp2f2([1.5, 2.5], [3.5, 2.5], -2.5, 0.401_651_681_832_597_7);
    assert_hyp2f2([1.5, 2.5], [3.5, 2.5], 3.0, 4.636_813_090_326_332);
    assert_hyp2f2([-3.0, 2.5], [4.0, -3.0], -12.5, 0.010_921_207_435_195_95);
    assert_eq!(hyp2f2(1.5, 2.5, 2.5, 1.5, -3.0), Some((-3.0_f64).exp()));
}

#[test]
fn test_hyp2f2_asymptotic() {
    // Arb references
    for (a, b, z, expected) in [
        ([1.5, 2.25], [3.5, 4.75], 500.0, 5.523_108_556_223_365e206),
        (
            [-10.5, 0.25],
            [20.0, 3.5],
            300.0,
            -1.057_385_809_228_637_7e71,
        ),
        (
            [1.5, 2.25],
            [3.5, 4.75],
            -500.0,
            0.001_936_926_332_203_224_2,
        ),
        ([0.5, -1.25], [2.0, 3.0], -100.0, 17.512_557_901_365_43),
        (
            [3.0, 5.5],
            [3.5, -4.25],
            -1000.0,
            -1.124_031_587_749_003_6e-8,
        ),
        ([2.5, 3.25], [-0.5, 1.75], -60.0, 4.645_783_675_124_214e-5),
        // a common parameter cancels, also at a pole
        ([3.0, 5.5], [3.0, 7.0], -200.0, 1.770_351_282_323_549_3e-10),
        ([-3.0, 2.5], [-3.0, 4.0], -300.0, 4.324_983_868_207_274_5e-6),
        // c_3 of the exponentially large part vanishes, but not c_4
        ([-1.5, -1.5], [-2.5, 1.5], 60.0, -4.921_886_930_234_788e21),
        ([-1.5, -1.5], [-2.5, 1.5], 100.0, -4.114_326_732_425_412e38),
    ] {
        assert_hyp2f2(a, b, z, expected);
    }
    // overflow
    assert_eq!(hyp2f2(1.5, 2.25, 3.5, 4.75, 800.0), Some(f64::INFINITY));
    assert_eq!(hyp2f2(2.0, 4.0, 4.0, 7.5, 900.0), Some(f64::INFINITY));
}

#[test]
fn test_hyp2f2_none() {
    // the series cancels, and |z| is too small for the asymptotic expansion
    assert_eq!(hyp2f2(2.5, 3.25, -0.5, 1.75, -40.0), None);
    // an integer a_1 - a_2
    assert_eq!(hyp2f2(1.0, 2.0, 3.5, 4.75, -200.0), None);
    // the result underflows
    assert_eq!(hyp2f2(80.0, 80.5, 100.0, 100.0, -1e6), None);
    assert_eq!(hyp2f2(1.0, 2.0, 1.0, 2.0, -745.0), None);
}

/// Asserts that `hyp1f2` is within a relative error of `1e-13` of `expected`
fn assert_hyp1f2(a: f64, b: [f64; 2], z: f64, expected: f64) {
    let actual = hyp1f2(a, b[0], b[1], z).unwrap();
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
}

#[test]
fn test_hyp1f2_special() {
    assert_eq!(hyp1f2(1.5, 3.5, 4.5, 0.0), Some(1.0));
    assert!(hyp1f2(f64::NAN, 3.5, 4.5, 0.25).unwrap().is_nan());
    assert!(hyp1f2(1.5, 3.5, 4.5, f64::NAN).unwrap().is_nan());
    // the limits of `hyppfq` for an infinite argument
    assert_eq!(hyp1f2(1.5, f64::INFINITY, 4.5, 0.25), Some(1.0));
    assert_eq!(
        hyp1f2(1.5, -2.5, 3.5, f64::INFINITY),
        Some(f64::NEG_INFINITY)
    );
    assert_eq!(
        hyp1f2(-3.0, 2.5, 3.5, f64::NEG_INFINITY),
        Some(f64::INFINITY)
    );
    assert_eq!(hyp1f2(1.5, 0.25, 0.5, f64::NEG_INFINITY), None);

    // a zero of (b_2)_n, also if the polynomial terminates after it
    assert_eq!(hyp1f2(1.5, 3.5, -2.0, 0.25), Some(f64::INFINITY));
    assert_eq!(hyp1f2(-3.0, 3.5, -2.0, 0.25), Some(f64::INFINITY));

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
    ] {
        assert_hyp1f2(a, b, z, expected);
    }
    // overflow
    assert_eq!(hyp1f2(1.5, 0.5, 3.5, 2e5), Some(f64::INFINITY));
    assert_eq!(hyp1f2(1.5, -0.5, 3.5, 2e5), Some(f64::NEG_INFINITY));
    assert_eq!(hyp1f2(1.5, 2.25, 3.5, 2e5), Some(f64::INFINITY));
    assert_eq!(hyp1f2(-0.5, 2.25, 3.5, 2e5), Some(f64::NEG_INFINITY));
}

#[test]
fn test_hyp1f2_none() {
    // the series cancels, and |z| is too small for the asymptotic expansion
    assert_eq!(hyp1f2(25.5, 1.25, 0.05, -4000.0), None);
    // the phase 2 sqrt(-z) of the oscillating part is too large
    assert_eq!(hyp1f2(0.5, 0.5, 0.25, -1e40), None);
}

#[test]
fn test_pfq_extreme_products() {
    // exact, where a partial product of the parameters is about 10^-320
    let (tiny, huge) = (1e-160, 1e160);
    assert_eq!(
        hyppfq(&[-1.0, tiny, tiny, huge, huge], &[], 2.0),
        Some(-1.0)
    );
    // exact, where a partial product of the lower parameters is 2^-1023
    let (tiny, huge) = (f64::MIN_POSITIVE / 2.0, 2.0_f64.powi(1023));
    assert_eq!(hyp1f2(-1.0, tiny, -huge, 1.0), Some(2.0));
    assert_eq!(hyp2f2(-1.0, 2.0, tiny, -huge, 1.0), Some(3.0));
    assert_eq!(hyp3f2(-1.0, 1.0, 2.0, tiny, -huge, 1.0), Some(3.0));
    // the power of 2 of a product exceeds `i32`, and the second term overflows
    let mut a = vec![huge; 4_198_404];
    a.extend([-1.0, 16.0]);
    assert_eq!(hyppfq(&a, &[], 2.0), None);
    // subnormal parameters (mpmath)
    for (actual, expected) in [
        (
            hyp1f2(1.390_671_161_567e-309, 1.618_95e-319, -2.102, -4.284),
            161_971_933_906.752_96,
        ),
        (
            hyp2f2(
                -3.0,
                -1.326_247_37e-315,
                -1.536_324_334_650_837e-45,
                8.289_046e-317,
                0.829,
            ),
            -7.395_713_298_118_337e45,
        ),
        (
            hyp3f2(
                -2.534,
                -1.358_077_306_22e-312,
                -3.0,
                -1.055_757_816_013_944_4e-249,
                -5.005_207_737_957_752_3e-147,
                0.019,
            ),
            -3.820_474_752_858_302e82,
        ),
    ] {
        let actual = actual.unwrap();
        assert!(
            (actual - expected).abs() <= 1e-15 * expected.abs(),
            "{actual} != {expected}"
        );
    }
}

#[test]
fn test_hyppfq_routing() {
    let z = -0.75;
    assert_eq!(hyppfq(&[], &[], z), Some(hyp0f0(z)));
    assert_eq!(hyppfq(&[2.5], &[], z), Some(hyp1f0(2.5, z)));
    assert_eq!(hyppfq(&[], &[1.5], z), Some(hyp0f1(1.5, z)));
    assert_eq!(hyppfq(&[0.5], &[1.5], z), Some(hyp1f1(0.5, 1.5, z)));
    assert_eq!(
        hyppfq(&[0.3, 0.7], &[1.5], z),
        Some(hyp2f1(0.3, 0.7, 1.5, z))
    );
    assert_eq!(hyppfq(&[0.5], &[1.5, 2.5], z), hyp1f2(0.5, 1.5, 2.5, z));
    assert_eq!(
        hyppfq(&[0.5, -1.5], &[1.5, 2.5], z),
        hyp2f2(0.5, -1.5, 1.5, 2.5, z)
    );
    assert_eq!(
        hyppfq(&[0.5, 1.5, 3.5], &[2.25, 4.5], z),
        hyp3f2(0.5, 1.5, 3.5, 2.25, 4.5, z)
    );
}

#[test]
fn test_hyppfq_cancel() {
    let z = 0.25;
    let f1f1 = Some(hyp1f1(0.5, 1.5, z));
    assert_eq!(hyppfq(&[0.5, 2.5], &[2.5, 1.5], z), f1f1);
    // also a non-positive integer, which would make it a polynomial
    assert_eq!(hyppfq(&[-3.0, 0.5], &[1.5, -3.0], z), f1f1);
    // each b_j cancels one a_i
    assert_eq!(hyppfq(&[2.5, 2.5], &[2.5], z), Some(hyp1f0(2.5, z)));
    assert_eq!(hyppfq(&[1.0, 2.0], &[2.0, 1.0], z), Some(hyp0f0(z)));
    assert_eq!(
        hyppfq(&[1.0, 2.0, 3.0, 4.0, 5.0], &[5.0, 4.0, 3.0, 2.0], z),
        Some(hyp1f0(1.0, z))
    );
}

#[test]
fn test_hyppfq_special() {
    assert_eq!(hyppfq(&[0.5], &[-2.0], 0.0), Some(1.0));
    // a pole, unless it is a polynomial of a lower degree
    assert_eq!(hyppfq(&[], &[-2.0], 0.5), Some(f64::INFINITY));
    assert_eq!(hyppfq(&[-3.0], &[-2.0], 0.5), Some(f64::INFINITY));
    assert_eq!(hyppfq(&[-1.0], &[-2.0], 0.5), Some(1.25));
    // complex for p = q + 1 and z > 1, unless it is a polynomial or (1 - z)^-a
    assert!(hyppfq(&[0.3, 0.7], &[1.5], 2.0).is_some_and(f64::is_nan));
    let (a, b) = ([0.5, 1.5, 2.5, 3.5], [1.25, 2.25, 3.25]);
    assert!(hyppfq(&a, &b, 2.0).is_some_and(f64::is_nan));
    assert_eq!(hyppfq(&[2.0], &[], 3.0), Some(0.25));
    assert!(hyppfq(&[2.5], &[], 3.0).is_some_and(f64::is_nan));
    assert_eq!(hyppfq(&[2.5], &[], 1.0), Some(f64::INFINITY));
    assert_eq!(hyppfq(&[-2.5], &[], 1.0), Some(0.0));
    assert!(hyppfq(&[0.5], &[f64::NAN], 0.5).is_some_and(f64::is_nan));
}

#[test]
fn test_hyppfq_infinite() {
    let inf = f64::INFINITY;
    // equal infinities cancel, and z = 0 wins
    let f1f1 = Some(hyp1f1(0.5, 1.5, 0.25));
    assert_eq!(hyppfq(&[0.5, inf], &[inf, 1.5], 0.25), f1f1);
    assert_eq!(hyppfq(&[inf], &[-inf], 0.0), Some(1.0));
    // a polynomial, with the sign of its last term
    assert_eq!(hyppfq(&[-2.0], &[2.5], -inf), Some(inf));
    assert_eq!(hyppfq(&[-3.0], &[2.5], inf), Some(-inf));
    assert_eq!(hyppfq(&[-3.0], &[-4.5], -inf), Some(-inf));
    assert_eq!(hyppfq(&[-3.0, 2.5], &[-3.5, 4.5], inf), Some(inf));
    assert_eq!(hyppfq(&[-3.0], &[], inf), Some(-inf));
    assert_eq!(hyppfq(&[0.0, 2.5], &[1.5], -inf), Some(1.0));
    // for p <= q at z = +inf, with the sign of prod Gamma(b_j) / prod Gamma(a_i)
    assert_eq!(hyppfq(&[], &[], inf), Some(inf));
    assert_eq!(hyppfq(&[0.5], &[1.5], inf), Some(inf));
    assert_eq!(hyppfq(&[], &[-2.5], inf), Some(-inf));
    assert_eq!(hyppfq(&[-2.5], &[1.5], inf), Some(-inf));
    assert_eq!(hyppfq(&[-1.5, 2.5], &[-0.5, -1.75, 2.25], inf), Some(-inf));
    // for q <= p <= q + 1 at z = -inf if all a_i > 0
    assert_eq!(hyppfq(&[], &[], -inf), Some(0.0));
    assert_eq!(hyppfq(&[1.5], &[-2.5], -inf), Some(0.0));
    assert_eq!(hyppfq(&[1.5, 2.5, 3.5], &[-0.5, -2.25], -inf), Some(0.0));
    // the functions of p <= 2 and q <= 1, which have no limit for 0F1(0.25; z)
    assert_eq!(hyppfq(&[-2.5], &[], -inf), Some(inf));
    assert_eq!(hyppfq(&[2.0], &[], inf), Some(0.0));
    assert_eq!(hyppfq(&[-1.5], &[-0.5], -inf), Some(-inf));
    assert_eq!(hyppfq(&[], &[0.25], -inf), None);
    // complex for p = q + 1
    assert!(hyppfq(&[2.5], &[], inf).is_some_and(f64::is_nan));
    assert!(hyppfq(&[1.5, 2.5, 3.5], &[0.5, 4.5], inf).is_some_and(f64::is_nan));
    // the terms n >= 1 vanish for a b_j = +inf, unless it is complex
    assert_eq!(hyppfq(&[1.5, 2.5, 3.5], &[0.5, inf], -10.0), Some(1.0));
    assert_eq!(hyppfq(&[1.5, 2.5, 3.5], &[0.5, inf], 1.0), Some(1.0));
    assert_eq!(hyppfq(&[-2.0, 2.5, 3.5], &[0.5, inf], 2.0), Some(1.0));
    assert_eq!(hyppfq(&[1.5, 2.5, 3.5], &[0.5, inf], 2.0), None);
    // p > q + 1 diverges, unless it is a polynomial
    assert_eq!(hyppfq(&[0.5, 1.5, 2.5], &[inf], 0.5), None);
    assert_eq!(hyppfq(&[-2.0, 1.5, 2.5], &[inf], 0.5), Some(1.0));
    // a pole wins
    assert_eq!(hyppfq(&[1.5], &[-2.0, inf], 0.25), Some(inf));
    // no limit
    assert_eq!(hyppfq(&[1.5], &[-inf], 0.25), None);
    assert_eq!(hyppfq(&[1.5], &[inf], inf), None);
    assert_eq!(hyppfq(&[inf], &[2.5], 0.25), None);
}

#[test]
fn test_hyppfq_polynomial() {
    // the exact sums of the terms (Python `fractions`), rounded
    for (a, b, z, expected) in [
        (&[-5.0, 2.5][..], &[][..], -0.3, 35.447_635_937_499_996),
        (&[-10.0], &[2.5], 7.5, -0.230_642_806_801_874_54),
        (&[-3.0, 0.7], &[1.5], 2.0, 0.049_600_000_000_000_05),
        (
            &[-6.0, 1.5, 2.25, -0.5],
            &[3.5, 0.75, -8.0],
            2.5,
            -2.188_432_770_889_141_6,
        ),
        (
            &[-20.0, 1.25, 1.5, 0.5],
            &[1.75, 1.5, 3.5],
            5.0,
            20_470_378.725_082_297,
        ),
        (
            &[-80.0, 2.5, 0.5, 0.5],
            &[1.5, 1.5, 3.5],
            20.0,
            1.998_050_134_345_888_1e96,
        ),
        (
            &[-7.0, 0.5, 1.5, 2.5, 3.5],
            &[0.25, 1.25, 2.25, -9.5],
            -1.75,
            -132.181_417_408_943_6,
        ),
    ] {
        let actual = hyppfq(a, b, z).unwrap();
        assert!(
            (actual - expected).abs() <= f64::EPSILON * expected.abs(),
            "hyppfq({a:?}, {b:?}, {z}) = {actual} != {expected}"
        );
    }
    // the sum cancels too much, but (1 - z)^70 does not
    assert_eq!(hyppfq(&[-70.0], &[], 3.0), Some(2.0_f64.powi(70)));
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
        (
            &[0.5, 1.5, 2.5, 3.5][..],
            &[1.25, 2.25, 3.25][..],
            0.5,
            1.701_434_243_680_788_3,
        ),
        (
            &[-0.5, 1.25, 2.75, -1.5],
            &[3.25, -2.75, 0.75],
            0.3,
            0.882_694_919_188_910_9,
        ),
        (
            &[0.5, 1.25, 2.75, 3.5, 0.75, 1.75],
            &[1.5, 2.25, 3.25, 4.75, 5.5],
            0.99,
            1.030_868_874_524_536_5,
        ),
        // in double-double, as the terms cancel
        (
            &[6.5, 5.25, 4.75, 3.5],
            &[1.5, 1.25, 0.75],
            -0.55,
            -0.326_402_725_424_534,
        ),
        (
            &[3.5, 5.25, 4.75, 6.5, 2.25],
            &[1.5, 1.25, 0.75, 0.5],
            -0.3,
            -12.406_169_497_317_718,
        ),
    ] {
        assert_hyppfq(a, b, z, expected);
    }
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
        (
            &[0.5, 1.25, 2.75, 3.5],
            &[3.25, 4.75, 1.75],
            -1.0,
            0.839_958_420_055_184_8,
        ),
        (
            &[0.5, 1.25, 2.75, 3.5, 0.125],
            &[3.25, 4.75, 1.75, 2.5],
            -1.25,
            0.987_921_427_741_094,
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
        (
            &[1.0, 0.25, 0.25, 0.25],
            &[1.25; 3],
            -1.0,
            0.993_051_520_249_976_6,
        ),
        (
            &[1.0, 0.25, 0.25, 0.25],
            &[1.25; 3],
            0.875,
            1.008_602_569_162_220_2,
        ),
        (
            &[1.0, 0.375, 0.375, 0.375, 0.375],
            &[1.375; 4],
            -0.5,
            0.997_372_909_080_830_9,
        ),
        // 04/0002: a sum of 2F1(a, a_k; a_k + 1; z), also for z < -1
        (
            &[0.375, 0.25, 0.625, 0.875],
            &[1.25, 1.625, 1.875],
            0.75,
            1.011_649_399_153_953_5,
        ),
        (
            &[0.375, 0.25, 0.625, 0.875],
            &[1.25, 1.625, 1.875],
            -1.25,
            0.985_603_765_207_305_6,
        ),
        (
            &[0.375, 0.25, 0.625, 0.875],
            &[1.25, 1.625, 1.875],
            -3.0,
            0.970_328_335_697_056_8,
        ),
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
        (
            &[1.0, 0.25, 0.5, 0.75],
            &[1.25, 1.5, 1.75],
            -1.0,
            0.975_967_713_426_429_7,
        ),
        (
            &[1.0, 0.25, 0.5, 0.75, 2.375],
            &[1.25, 1.5, 1.75, 3.375],
            -1.0,
            0.982_477_625_546_169_4,
        ),
        // 06/0009: zeta(3) and ln 2
        (&[1.0; 4], &[3.0; 3], -1.0, 0.966_812_171_669_881_1),
        (&[1.0; 5], &[3.0; 4], -1.0, 0.988_341_738_313_950_6),
        // 06/0010: pi^3 / 32, and Hurwitz zeta functions for q = 4
        (&[1.0, 0.5, 0.5, 0.5], &[1.5; 3], -1.0, PI.powi(3) / 32.0),
        (
            &[1.0, 0.5, 0.5, 0.5, 0.5],
            &[1.5; 4],
            -1.0,
            0.988_944_551_741_105_3,
        ),
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
    assert_eq!(hyppfq(&a, &[1.25, 2.25, 3.25], 1.0), Some(f64::INFINITY));
    assert_eq!(hyppfq(&a, &[1.25, 2.25, 4.5], 1.0), Some(f64::INFINITY));
    // no method where it converges, as to zeta(6) (30/0011 of the formulas above)
    assert_eq!(hyppfq(&a, &[1.25, 2.25, 4.75], 1.0), None);
    assert_eq!(hyppfq(&[1.0; 7], &[2.0; 6], 1.0), None);
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
    assert_eq!(hyppfq(&[-60.0, 0.5, 0.5], &[1.5, 1.5, 1.5], 200.0), None);
    assert_eq!(hyppfq(&[0.5, 1.5], &[], 0.5), None);
    let (a, b) = ([10.5, 8.25, 9.75, 7.5], [0.5, 1.25, 0.75]);
    assert_eq!(hyppfq(&a, &b, -0.5), None);
    // z < -1 with integer a_i - a_j
    assert_eq!(hyppfq(&[1.0; 4], &[2.0; 3], -40.0), None);
    assert_eq!(hyppfq(&[2.0; 4], &[1.0; 3], -3.0), None);
}
