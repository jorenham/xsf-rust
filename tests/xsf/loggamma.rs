use num_complex::c64;

#[test]
fn test_loggamma_f64() {
    xsref::test("loggamma", "d-d", |x| xsf::loggamma(x[0]));
}

#[test]
fn test_loggamma_c64() {
    xsref::test("loggamma", "cd-cd", |x| xsf::loggamma(c64(x[0], x[1])));
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

#[test]
fn test_rgamma_f64() {
    xsref::test("rgamma", "d-d", |x| xsf::rgamma(x[0]));
}

#[test]
fn test_rgamma_f64_near_overflow() {
    // unlike in xsf: no spurious overflow for `-178 < x < -170` (mpmath)
    for (x, expected) in [
        (-170.917_162_917_922_44, -6.637_711_907_120_738e307),
        (-173.999_999_999_999_97, 1.826_215_811_702_844e302),
        (-176.000_000_000_000_03, -5.624_744_700_046_412e306),
    ] {
        let actual = xsf::rgamma(x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-15,
            "{actual:e} != {expected:e}"
        );
    }
    assert_eq!(xsf::rgamma(-172.5), f64::NEG_INFINITY);
    assert_eq!(xsf::rgamma(-173.5), f64::INFINITY);
}

#[test]
fn test_rgamma_c64() {
    xsref::test("rgamma", "cd-cd", |x| xsf::rgamma(c64(x[0], x[1])));
}
