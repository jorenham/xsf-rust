use core::f64::consts::PI;
use xsf::{hyp2f1, hyp2f2, hyp3f2};

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
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, f64::INFINITY, 0.25), None);
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, f64::NEG_INFINITY), None);

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
    assert_eq!(hyp2f2(1.5, 2.5, f64::INFINITY, 4.5, 0.25), None);
    assert_eq!(hyp2f2(1.5, 2.5, 3.5, 4.5, f64::NEG_INFINITY), None);

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
fn test_hyp2f2_none() {
    // the series cancels too much, also in double-double precision
    assert_eq!(hyp2f2(2.5, 3.25, -0.5, 1.75, -40.0), None);
    assert_eq!(hyp2f2(1.5, 2.25, 3.5, 4.75, -500.0), None);
    // the terms overflow
    assert_eq!(hyp2f2(1.5, 2.25, 3.5, 4.75, 800.0), None);
}
