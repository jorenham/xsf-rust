xsref::tests! {
    test_wright_bessel_f64: "wright_bessel", "d_d_d-d", |x| xsf::wright_bessel(x[0], x[1], x[2]);
    test_log_wright_bessel_f64: "log_wright_bessel", "d_d_d-d", |x| {
        xsf::log_wright_bessel(x[0], x[1], x[2])
    };
}

#[test]
fn test_log_wright_bessel_near_zero() {
    // Keep the precision when Phi is close to 1, where its logarithm is close to 0.
    // In the last three cases b is one ulp below 1, so log(Gamma(b)) is also close to 0.
    let b1 = 0.999_999_999_999_999_9;
    for (a, b, x, expected) in [
        (10.0, 1.0, 1e-3, 2.755_731_922_018_888_3e-10),
        (10.0, 2.0, 1e-3, 2.505_210_838_512_792_5e-11),
        (20.0, 1.0, 2.5, 1.027_579_405_828_041_2e-18),
        (1e-4, 1.0, 1e-12, 1.000_057_715_007_667_5e-12),
        (1e-4, 1.0, 1e-10, 1.000_057_715_007_667_5e-10),
        (100.0, 1.0, 1.5, 1.607_265_432_188_200_3e-158),
        (10.0, b1, 3.63e-10, 3.594_925_664_826_874e-17),
        (100.0, b1, 2.5, -6.408_381_213_480_008e-17),
        (1e-4, b1, 1e-20, -6.407_381_155_765e-17),
    ] {
        let actual = xsf::log_wright_bessel(a, b, x);
        assert!((actual / expected - 1.0).abs() < 1e-13, "{actual:e} != {expected:e}");
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

#[test]
fn test_wright_bessel_asymptotic_overflow() {
    // Regression test: in the large-x asymptotic expansion, pow(Z, 0.5 - b) * exp(.../a * Z)
    // overflowed in exp (or underflowed in pow) before the product was formed, so the
    // finite results below were returned as inf. The last case overflows: it was NaN.
    // The asymptotic expansion has a relative error of up to ~1e-7 here.
    //
    // Unlike in xsf, these cases use the series around the maximum term.
    let abx = [
        (0.75, 20.0, 10.0_f64.powf(4.5)),
        (2.0, 30.0, 1e8),
        (2.0, 30.0, 1.077e8),
        (4.9, 80.0, 5e15),
    ];
    let actual = abx.map(|(a, b, x)| xsf::wright_bessel(a, b, x));
    let expected = [
        1.153_548_523_932_893_3e271,
        3.980_356_527_302_694_7e298,
        6.652_307_057_991_812e307,
        3.334_148_696_466_592_6e91,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-6);
    assert_eq!(xsf::wright_bessel(4.9, 150.0, 4.9e18), f64::INFINITY);

    let actual = xsf::wright_bessel(4.0, 10.0, 2e13);
    let expected = 2.557_441_567_905_184e300;
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-12);
    assert_eq!(xsf::wright_bessel(1.0, 120.0, 1e9), f64::INFINITY);
    assert_eq!(xsf::wright_bessel(2.0, 30.0, 1e308), f64::INFINITY);
}

#[test]
fn test_log_wright_bessel_ax_overflow() {
    // a x overflows, but Z = (a x)^(1/(1+a)) of the asymptotic expansion does not (#2264)
    let abx = [(2.0, 30.0, 1e308), (4.9, 1.0, f64::MAX), (1.5, 0.0, 1e308)];
    let actual = abx.map(|(a, b, x)| xsf::log_wright_bessel(a, b, x));
    let expected = [8.772_053_214_638_599e102, 2.780_995_443_939_785e52, 3.106_599_394_326_126e123];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);
}
