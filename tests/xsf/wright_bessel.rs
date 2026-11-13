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
