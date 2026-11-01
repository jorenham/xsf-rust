#[test]
fn test_pow_rising() {
    xsref::test("poch", "d_d-d", |x| xsf::pow_rising(x[0], x[1]));
}

#[test]
fn test_pow_rising_lgam_ratio() {
    // unlike in xsf (scipy/xsf#357): `lgam_ratio` for 10 < a <= 1e4
    for (a, m, expected) in [
        (
            5_878.861_847_196_475_5,
            -0.029_080_395_984_806_877,
            0.776_941_938_775_515_5,
        ),
        (150.5, 0.75, 42.941_972_832_542_55),
        (50.25, -0.5, 0.142_132_904_591_272_65),
    ] {
        let actual = xsf::pow_rising(a, m);
        assert!(
            (actual / expected - 1.0).abs() < 5e-15,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_pow_rising_reflection() {
    // unlike in xsf (scipy/xsf#357): the reflection formula for a < -10
    for (a, m, expected) in [
        (
            -96_483_701_341.528_73,
            -0.001_068_826_099_618_102_1,
            0.973_630_092_687_070_2,
        ),
        (-10_000.25, 0.5, -100.002_499_984_374_61),
        (-50.75, 2.5, 17_673.802_508_037_614),
    ] {
        let actual = xsf::pow_rising(a, m);
        assert!(
            (actual / expected - 1.0).abs() < 5e-15,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_pow_falling_near_pole() {
    for (x, m, expected) in [
        (
            -59.000_000_000_1,
            -7.654_166_228_426_843,
            0.000_141_662_338_550_097_54,
        ),
        (-0.999_999_999_9, 20.5, -3.527_763_913_031_247e28),
        (-1_000.25, -104.5, -9.343_131_111_265e-312), // subnormal
    ] {
        let actual = xsf::pow_falling(x, m);
        assert!(
            (actual / expected - 1.0).abs() < 1e-12,
            "{actual:e} != {expected:e}"
        );
    }
}
