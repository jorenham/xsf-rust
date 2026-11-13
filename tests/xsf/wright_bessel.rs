xsref::tests! {
    test_wright_bessel_f64: "wright_bessel", "d_d_d-d", |x| xsf::wright_bessel(x[0], x[1], x[2]);
    test_log_wright_bessel_f64: "log_wright_bessel", "d_d_d-d", |x| {
        xsf::log_wright_bessel(x[0], x[1], x[2])
    };
}

#[test]
fn test_log_wright_bessel_near_zero() {
    // unlike in xsf, ln(Phi) does not lose accuracy for Phi ~ 1 (mpmath)
    for (a, b, x, expected) in [
        (10.0, 1.0, 1e-3, 2.755_731_922_018_888_3e-10), // was 2.6e-7 off
        (10.0, 2.0, 1e-3, 2.505_210_838_512_792_5e-11), // was 3.0e-6 off
        (20.0, 1.0, 2.5, 1.027_579_405_828_041_2e-18),  // was 0
        (1e-4, 1.0, 1e-12, 1.000_057_715_007_667_5e-12), // was 5.8e-5 off
    ] {
        let actual = xsf::log_wright_bessel(a, b, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-13,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_wright_bessel_large_a() {
    // unlike in xsf, the terms around the maximum term account for b (Arb)
    let actual = xsf::wright_bessel(5.0, 100.0, 1e12);
    let expected = 4.066_520_826_827_190_5e-147; // was 3.3e-2 off
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-13);
    let abx = [(5.0, 100.0, 1e12), (5.0, 100.0, 1e20), (177.0, 177.0, 2.5)];
    let actual = abx.map(|(a, b, x)| xsf::log_wright_bessel(a, b, x));
    let expected = [
        -337.077_220_869_952_3,  // was 9.9e-5 off
        2_588.247_122_998_684_5, // was 1.2e-6 off
        -737.509_837_141_777_4,  // was 3.5e-7 off
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_wright_bessel_huge_x() {
    // unlike in xsf, no `int` overflow of k_max = 2147367898 and 2.3e9 (Arb)
    let actual = [3.064e59, 1e60].map(|x| xsf::log_wright_bessel(5.0, 0.0, x));
    let expected = [
        12_884_207_403.503_17,  // was -inf
        15_691_925_841.974_194, // was -inf
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
    assert_eq!(xsf::wright_bessel(5.0, 0.0, 3.064e59), f64::INFINITY);
    assert_eq!(xsf::wright_bessel(5.0, 0.0, 1e60), f64::INFINITY);
}

#[test]
fn test_wright_bessel_peak() {
    // unlike in xsf, the terms around the maximum term instead of the integral (Arb)
    let abx = [(0.5, 100.0, 1e3), (1.0, 100.0, 1e5)];
    let actual = abx.map(|(a, b, x)| xsf::wright_bessel(a, b, x));
    let expected = [
        3.411_236_758_044_524_3e-117, // was 5.1e-8 off
        1.026_933_459_623_076_3e22,   // was NaN
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);
    let abx = [(0.001, 100.0, 1e5), (0.001, 50.0, 1e5), (1.0, 100.0, 1e5)];
    let actual = abx.map(|(a, b, x)| xsf::log_wright_bessel(a, b, x));
    let expected = [
        99_143.465_191_656_53, // was NaN
        99_401.240_922_855_65, // was 4.4e-4 off
        50.683_449_183_697,    // was NaN
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_wright_bessel_asymptotic() {
    // unlike in xsf, the series around the maximum term just above Z = 14 + b^2 / (2 (1 + a)) (Arb)
    let abx = [(0.5, 10.0, 709.782_712_893_384_1), (4.0, 1.0, 224_388.0)];
    let actual = abx.map(|(a, b, x)| xsf::wright_bessel(a, b, x));
    let expected = [
        2.680_788_404_494_657e48, // was 8.1e-8 off
        11_866_812.621_349_995,   // was 1.5e-6 off
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);
}
