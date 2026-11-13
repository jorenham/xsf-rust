xsref::tests! {
    test_pow_rising: "poch", "d_d-d", |x| xsf::pow_rising(x[0], x[1]);
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
fn test_pow_rising_poles() {
    // unlike in xsf (scipy/xsf#357), the limit where both a and a + m are poles, which the loops
    // can't reach when they overflow or can't reduce m (Arb); the poch table covers the others
    assert_eq!(xsf::pow_rising(-300.0, 250.0), f64::INFINITY);
    let am = [
        (-1e17, -10.0),
        (-1e17, 3.0),
        (-1_152_921_504_606_846_976.0, 1.0),
    ];
    let actual = am.map(|(a, m)| xsf::pow_rising(a, m));
    let expected = [
        9.999_999_999_999_994e-171,
        -1e51,
        -1_152_921_504_606_846_976.0,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_pow_rising_sum_near_pole() {
    // unlike in xsf (scipy/xsf#357), without the rounding error of a + m near (or onto) a pole
    // of Gamma(a + m) (Arb)
    let am = [
        (-0.033_978_567_503_241_18, -0.966_021_432_496_758_9),
        (-0.013_566_220_790_266_651, 0.013_566_220_790_266_646),
    ];
    let actual = am.map(|(a, m)| xsf::pow_rising(a, m));
    let expected = [-685_298_206_918_409.0, 2_586_070_687_183_794.2];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
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
        // x - m + 1 near (or rounded onto) a pole
        (
            -6.010_755_014_317_816,
            1.989_244_985_682_278_7,
            3.627_414_616_553_739_6e-10,
        ),
        (
            0.862_338_700_445_937_9,
            47.862_338_700_445_94,
            -5.220_110_314_612_468e42,
        ),
        (2.5e-17, 1.0, 2.5e-17),
        // sin(pi x) is subnormal, and the product of the first two factors would be
        (5e-324, 2.5, 0.423_142_187_660_817_2),
        (-999.999_999_999_999_9, -102.0, 2.309_688_563_031_087e-304),
    ] {
        let actual = xsf::pow_falling(x, m);
        assert!(
            (actual / expected - 1.0).abs() < 1e-12,
            "{actual:e} != {expected:e}"
        );
    }
}
