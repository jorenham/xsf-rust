use xsf::{hyp2f1, hyp3f2};

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

    // a common parameter cancels
    for z in [-2.5, 0.25, 0.95, 1.0] {
        assert_eq!(
            hyp3f2(1.5, 2.5, 3.5, 2.5, 6.5, z),
            Some(hyp2f1(1.5, 3.5, 6.5, z))
        );
        assert_eq!(
            hyp3f2(-2.0, 1.5, -3.0, 4.5, -2.0, z),
            Some(hyp2f1(1.5, -3.0, 4.5, z))
        );
    }
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
fn test_hyp3f2_none() {
    // |z| = 1
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, -1.0), None);
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 8.5, 1.0), None);
    // z < -1 with an integer a_2 - a_1
    assert_eq!(hyp3f2(1.5, 2.5, 3.25, 0.5, 4.5, -2.5), None);
    // more than 10^4 terms
    assert_eq!(hyp3f2(1.5, 2.5, 3.5, 0.5, 4.5, 0.999), None);
    assert_eq!(hyp3f2(1.5, 2.25, 3.125, 0.5, 4.5, -1.001), None);
    // the terms cancel
    assert_eq!(hyp3f2(50.5, 60.5, 70.5, 1.5, 2.5, -0.75), None);
    // near a zero of the polynomial
    assert_eq!(
        hyp3f2(-2.0, 1.0, 1.0, 0.5, 0.5, 0.143_237_254_218_789_44),
        None
    );
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
}
