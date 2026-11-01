#[test]
fn test_beta() {
    xsref::test("beta", "d_d-d", |x| xsf::beta(x[0], x[1]));
}

#[test]
fn test_betaln() {
    xsref::test("betaln", "d_d-d", |x| xsf::betaln(x[0], x[1]));
}

#[test]
fn test_beta_lgam_ratio() {
    // unlike in xsf (scipy/scipy#26423): no cancellation in lgam(a + b) - lgam(a) for |b| << a
    for (a, b, expected) in [
        (999_999.0, 0.5, 0.001_772_454_958_690_183_8),
        (300_000.0, -2.5, -46_598_342_334_957.34),
    ] {
        let actual = xsf::beta(a, b);
        assert!(
            (actual / expected - 1.0).abs() < 5e-15,
            "{actual:e} != {expected:e}"
        );
    }
    for (a, b, expected) in [
        (1e100, 1e94, -1.481_551_105_796_410_7e95),
        (10_000.0, 0.25, -1.014_553_193_217_846_2),
    ] {
        let actual = xsf::betaln(a, b);
        assert!(
            (actual / expected - 1.0).abs() < 5e-15,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_beta_reflected() {
    // unlike in xsf (scipy/xsf#357): no cancellation in lgam(a) + lgam(b) - lgam(a + b) for
    // |b| << -a or 0 <= a + b << a
    for (a, b, expected) in [
        (-1_000_000.25, 0.5, -0.001_772_453_407_792_191_7),
        (5_000.25, -4_999.5, 0.304_878_142_880_374_1),
    ] {
        let actual = xsf::beta(a, b);
        assert!(
            (actual / expected - 1.0).abs() < 5e-15,
            "{actual:e} != {expected:e}"
        );
    }
    let actual = xsf::betaln(-1_000_000.25, 0.5);
    assert!((actual / -6.335_390_586_057_390_5 - 1.0).abs() < 5e-15);
}

#[test]
fn test_beta_near_171() {
    // unlike in xsf (scipy/xsf#357): no under- or overflow of the Gamma functions for arguments
    // below -170.6
    for (a, b, expected) in [
        (-171.1, 0.5, -0.416_731_771_605_593_1),
        (
            -171.331_359_822_697_1,
            0.004_863_556_345_095_009,
            198.157_298_160_572_34,
        ),
        (171.5, -171.25, 0.025_871_582_394_479_65),
    ] {
        let actual = xsf::beta(a, b);
        assert!(
            (actual / expected - 1.0).abs() < 5e-14,
            "{actual:e} != {expected:e}"
        );
        let actual = xsf::betaln(a, b);
        let expected = expected.abs().ln();
        assert!(
            (actual / expected - 1.0).abs() < 5e-14,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_beta_negint_int_overflow() {
    // unlike in xsf, 1 - a doesn't overflow the C++ `int` for a = -2^31 and -2^31 + 1
    for a in [-2_147_483_648.0, -2_147_483_647.0] {
        let expected = 1.0 / a; // (-1)^1 B(1 - a - 1, 1)
        assert!((xsf::beta(a, 1.0) / expected - 1.0).abs() < 1e-14);
        assert!((xsf::betaln(a, 1.0) / expected.abs().ln() - 1.0).abs() < 1e-14);
    }
}
