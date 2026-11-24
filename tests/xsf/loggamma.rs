use num_complex::c64;

xsref::tests! {
    test_loggamma_f64: "loggamma", "d-d", |x| xsf::loggamma(x[0]);
    test_loggamma_c64: "loggamma", "cd-cd", |x| xsf::loggamma(c64(x[0], x[1]));
}

#[test]
fn test_loggamma_c64_taylor_cutoff() {
    for (z, expected) in [
        (
            c64(0.919_660_915_069_406_1, 0.183_154_665_331_011_46),
            c64(0.021_089_888_719_317_225, -0.128_428_682_915_316_14),
        ),
        (
            c64(1.937_555_037_215_635_2, 0.190_001_649_000_368_64),
            c64(-0.037_211_529_359_037_95, 0.073_021_395_401_304_9),
        ),
    ] {
        let w = xsf::loggamma(z);
        assert!((w - expected).norm() < 2e-15 * expected.norm());
        assert_eq!(xsf::loggamma(z.conj()), w.conj());
    }
}

#[test]
fn test_loggamma_c64_stirling_overflow() {
    let x = 8.988_465_674_311_05e307;
    for y in [x, -x] {
        assert_eq!(xsf::loggamma(c64(-x, y)), c64(f64::NEG_INFINITY, f64::INFINITY.copysign(y)));
    }
    let z = c64(-6e307, 2.659_791_818_217_255_6e305);
    let w = xsf::loggamma(z);
    assert_eq!(w.re, f64::NEG_INFINITY);
    assert!((w.im / 8.224_587_161_873_372e290 - 1.0).abs() < 1e-13);
    assert_eq!(xsf::loggamma(z.conj()), w.conj());

    let z = c64(3.978_414_405_427_773e305, f64::MAX);
    let w = xsf::loggamma(z);
    assert!((w.re / -2.901_477_737_953_569e292 - 1.0).abs() < 1e-13);
    assert_eq!(w.im, f64::INFINITY);
    assert_eq!(xsf::loggamma(z.conj()), w.conj());
}

#[test]
fn test_loggamma_c64_subnormal_angle() {
    let x = 8.988_465_674_311_05e307;
    let z = c64(x, 1e-30);
    let w = xsf::loggamma(z);
    assert_eq!(w.re, f64::INFINITY);
    xsref::np_assert_allclose!([w.im], [7.090_895_657_128_241e-28], rtol = 4.0 * f64::EPSILON);
    assert_eq!(xsf::loggamma(z.conj()), w.conj());

    for y in [f64::from_bits(1), -f64::from_bits(1)] {
        let w = xsf::loggamma(c64(x, y));
        assert_eq!(w.re, f64::INFINITY);
        assert_eq!(w.im.to_bits(), f64::from_bits(709).copysign(y).to_bits());
    }
}

#[test]
fn test_loggamma_c64_negative_axis() {
    // the sign of the zero imaginary part selects the side of the branch cut
    let w = xsf::loggamma(c64(-2.5, 0.0));
    assert_eq!(xsf::loggamma(c64(-2.5, -0.0)), w.conj());
    // mpmath: loggamma(-2.5 + 1e-30j)
    assert!((w.re / -0.056_243_716_497_674_05 - 1.0).abs() < 1e-13);
    assert_eq!(w.im, -3.0 * core::f64::consts::PI);

    // poles of gamma, zeros of rgamma
    let w = xsf::loggamma(c64(-2.0, 0.0));
    assert!(w.re.is_nan() && w.im.is_nan());
    assert_eq!(xsf::rgamma(c64(-2.0, 0.0)), c64(0.0, 0.0));
}

xsref::tests! {
    test_rgamma_f64: "rgamma", "d-d", |x| xsf::rgamma(x[0]);
}

#[test]
fn test_rgamma_f64_subnormal() {
    // 1 / Gamma(x) = x for subnormal x, also for x < 0 (scipy/xsf#417)
    for x in [5e-324, -5e-324, 1e-310, -1e-310, -2e-308] {
        assert_eq!(xsf::rgamma(x), x);
    }
}

#[test]
fn test_rgamma_f64_near_overflow() {
    // no spurious overflow for `-178 < x < -170` (mpmath)
    for (x, expected) in [
        (-170.917_162_917_922_44, -6.637_711_907_120_738e307),
        (-173.999_999_999_999_97, 1.826_215_811_702_844e302),
        (-176.000_000_000_000_03, -5.624_744_700_046_412e306),
    ] {
        let actual = xsf::rgamma(x);
        assert!((actual / expected - 1.0).abs() < 1e-15, "{actual:e} != {expected:e}");
    }
    assert_eq!(xsf::rgamma(-172.5), f64::NEG_INFINITY);
    assert_eq!(xsf::rgamma(-173.5), f64::INFINITY);
}

xsref::tests! {
    test_rgamma_c64: "rgamma", "cd-cd", |x| xsf::rgamma(c64(x[0], x[1]));
}

#[test]
fn test_rgamma_c64_negative_real_axis() {
    for im in [0.0_f64, -0.0] {
        for (x, expected, im_neg) in [
            (-0.5, -0.282_094_791_773_878_14, false),
            (-0.9, -0.094_602_330_550_059_98, true),
            (-169.1, 7.016_094_122_732_882e303, true),
        ] {
            let z = xsf::rgamma(c64(x, im));
            assert!((z.re / expected - 1.0).abs() <= 4.0 * f64::EPSILON);
            assert_eq!(z.im.to_bits(), if im_neg { -im } else { im }.to_bits());
        }

        let min_subnormal = f64::from_bits(1);
        let z = xsf::rgamma(c64(-min_subnormal, im));
        assert_eq!(z.re, -min_subnormal);
        assert_eq!(z.im.to_bits(), im.to_bits());
    }
}

#[test]
fn test_rgamma_c64_overflow_at_large_even_real_part() {
    let inf = f64::INFINITY;
    for (x, y, expected) in [
        (-9_007_199_254_740_992.0, f64::from_bits(1), c64(inf, inf)),
        (-2.998_076_996_061_15e154, 1.268_434_288_203_711_7e-10, c64(inf, inf)),
        (-2.998_076_996_061_15e154, 8.333_333_333_333_334, c64(-inf, inf)),
        (-8.988_465_674_311_05e307, 10.0, c64(-inf, -inf)),
    ] {
        assert_eq!(xsf::rgamma(c64(x, y)), expected);
        assert_eq!(xsf::rgamma(c64(x, -y)), expected.conj());
        assert_eq!(xsf::rgamma(c64(x, 0.0)), c64(0.0, 0.0));
    }
}

#[test]
fn test_rgamma_c64_near_nonpositive_integer() {
    for (x, re, im) in [
        (0.0, -5.772_156_649_015_329_4e-61, 1.0e-30),
        (-5.0, -2.047_341_202_118_160_8e-58, -1.2e-28),
        (-10.0, 8.534_039_795_205_319e-54, 3.6288e-24),
    ] {
        let result = xsf::rgamma(c64(x, 1e-30));
        assert!((result.re / re - 1.0).abs() <= 4.0 * f64::EPSILON);
        assert!((result.im / im - 1.0).abs() <= 4.0 * f64::EPSILON);
        assert_eq!(xsf::rgamma(c64(x, -1e-30)), result.conj());
    }
}
