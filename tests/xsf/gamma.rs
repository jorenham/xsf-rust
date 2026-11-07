use num_complex::c64;

#[test]
fn test_gamma_f64() {
    xsref::test("gamma", "d-d", |x| xsf::gamma(x[0]));
}

#[test]
fn test_gamma_c64() {
    xsref::test("gamma", "cd-cd", |x| xsf::gamma(c64(x[0], x[1])));
}

#[test]
fn test_gamma_c64_real_axis() {
    let w = xsf::gamma(c64(-2.0, 0.0));
    assert!(w.re.is_nan() && w.im.is_nan());
    // overflow: inf + 0i instead of inf + NaN i
    let w = xsf::gamma(c64(200.0, 0.0));
    assert_eq!(w.re, f64::INFINITY);
    assert_eq!(w.im.to_bits(), 0.0_f64.to_bits());
}

#[test]
fn test_gammainc() {
    xsref::test("gammainc", "d_d-d", |x| xsf::gammainc(x[0], x[1]));
}

#[test]
fn test_gammaincc() {
    xsref::test("gammaincc", "d_d-d", |x| xsf::gammaincc(x[0], x[1]));
}

#[test]
fn test_gammainc_large_a() {
    // unlike in xsf (scipy/xsf#339): the power series is no longer truncated for large `a`
    // (mpmath)
    for (a, x, expected) in [
        (300_000.0, 297_500.0, 2.364_252_205_564_727e-6),
        (1e6, 995_000.0, 2.749_580_359_270_071e-7),
        (1e8, 99_950_000.0, 2.854_642_139_958_626e-7),
        (1e10, 9_999_500_000.0, 2.865_326_545_108_891e-7),
        (1e12, 999_970_000_000.0, 4.862_750_805_536_816e-198),
    ] {
        let actual = xsf::gammainc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-13,
            "{actual:e} != {expected:e}"
        );
    }
    let q = xsf::gammaincc(1e8, 99_950_000.0);
    assert!((q / 0.999_999_714_535_786 - 1.0).abs() < 1e-15, "{q}");
    let x = xsf::gammaincinv(1e10, 2.865_326_545_108_891e-7);
    assert!((x / 9_999_500_000.0 - 1.0).abs() < 1e-15, "{x}");
}

#[test]
fn test_gammainc_huge_a() {
    // unlike in xsf, these are not NaN (Arb: min(P, Q) < 10^(-6e306))
    assert_eq!(xsf::gammainc(f64::MAX, 3.0), 0.0);
    assert_eq!(xsf::gammaincc(4.4e306, 3.7e180), 1.0);
    assert_eq!(xsf::gammaincc(2.532_737_276_080_075_8e305, f64::MAX), 0.0); // MAX / ln(MAX)
    assert_eq!(xsf::chdtr(1e308, 3.0), 0.0);
    assert_eq!(xsf::pdtr(1e307, 1e306), 1.0);
}

#[test]
fn test_gammainc_tiny_a() {
    // unlike in xsf, these are accurate for tiny `a` (Arb)
    assert_eq!(xsf::gammainc(1e-310, 0.5), 1.0); // was 0
    assert_eq!(xsf::gammainc(1e-300, 1.0), 1.0); // was 1.0000000000000238
    for (a, x, expected) in [
        (3e-309, 1e-4, 2.589_967_411_372_412e-308),
        (1e-300, 1.0, 2.193_839_343_955_203e-301),
        (1e-200, 2.0, 4.890_051_070_806_111_5e-202),
    ] {
        let actual = xsf::gammaincc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-14,
            "{actual:e} != {expected:e}"
        );
    }
    let x = xsf::gammainccinv(1e-200, 1e-210);
    assert!((x / 19.984_172_765_540_563 - 1.0).abs() < 1e-15, "{x}");
}

#[test]
fn test_gammaincinv() {
    xsref::test("gammaincinv", "d_d-d", |x| xsf::gammaincinv(x[0], x[1]));
}

#[test]
fn test_gammainccinv() {
    xsref::test("gammainccinv", "d_d-d", |x| xsf::gammainccinv(x[0], x[1]));
}

#[test]
fn test_gammaln() {
    xsref::test("gammaln", "d-d", |x| xsf::gammaln(x[0]));
}

#[test]
fn test_gammasgn() {
    xsref::test("gammasgn", "d-d", |x| xsf::gammasgn(x[0]));
}
