use num_complex::c64;

#[test]
fn test_riemann_zeta_f64() {
    xsref::test("riemann_zeta", "d-d", |x| xsf::riemann_zeta(x[0]));
}

#[test]
fn test_riemann_zeta_c64() {
    xsref::test("riemann_zeta", "cd-cd", |x| {
        xsf::riemann_zeta(c64(x[0], x[1]))
    });
}

#[test]
fn test_riemann_zeta_c64_extreme() {
    // (k + 1)^z overflows in the Borwein sum, and 1 / inf = 0 (C99 Annex G)
    let w = xsf::riemann_zeta(c64(1000.0, 0.5));
    assert_eq!((w.re, w.im), (1.0, -3.170_083_028_960_87e-302));
    // the reflection factor overflows, and zeta(1 - z) is finite
    assert_eq!(xsf::riemann_zeta(c64(-1000.0, 0.5)).re, f64::INFINITY);
    // too many Euler-Maclaurin terms (undefined behavior in xsf)
    for im in [7.4e19, 1e300, f64::INFINITY] {
        let w = xsf::riemann_zeta(c64(3.0, im));
        assert!(w.re.is_nan() && w.im.is_nan());
    }
}

#[test]
fn test_zeta_f64() {
    xsref::test("zeta", "d_d-d", |x| xsf::zeta(x[0], x[1]));
}

#[test]
fn test_zeta_large_q() {
    // asymptotic expansion for q > 1e8 (reference value from mpmath)
    let (x, q) = (36.0, 1.1e8);
    let expected = 1.016_688_811_413_570_4e-283;
    let rel = (xsf::zeta(x, q) - expected).abs() / expected.abs();
    assert!(rel < 1e-15, "rel = {rel:e}");

    assert_eq!(xsf::zeta(f64::INFINITY, q), 0.0);
    assert_eq!(xsf::zeta(f64::INFINITY, f64::INFINITY), 0.0);
}

#[test]
fn test_zeta_near_1_tiny_q() {
    // the tail sum ~ 1 / (x - 1) isn't negligible (reference values from mpmath)
    for (x, q, expected) in [
        (1.0 + f64::EPSILON, 1e-15, 5.503_599_627_370_504e15),
        (1.000_000_000_001, 1e-17, 1.000_009_999_150_220_6e17),
    ] {
        let rel = (xsf::zeta(x, q) - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "zeta({x}, {q}): rel = {rel:e}");
    }
}

#[test]
fn test_zeta_negative_q() {
    // reference values from mpmath, by summing the terms with q + k < 0 directly
    for (x, q, expected) in [
        (4.0, -431.636_697_197_915_6, 63.992_140_081_094_475),
        (2.0, -1e10 - 0.5, 9.869_604_400_989_358),
        (31.0, -7.5, 1.593_135_539_215_446_3e-29),
        (31.0, -0.500_000_000_000_000_1, 3.304_063_430_367_368e-5),
    ] {
        let rel = (xsf::zeta(x, q) - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "zeta({x}, {q}): rel = {rel:e}");
    }
    // all terms underflow
    assert_eq!(xsf::zeta(1e300, 1.5), 0.0);
    assert_eq!(xsf::zeta(f64::INFINITY, 1.5), 0.0);
    assert_eq!(xsf::zeta(f64::INFINITY, -0.5), f64::INFINITY);
}

#[test]
fn test_zetac_f64() {
    xsref::test("zetac", "d-d", |x| xsf::zetac(x[0]));
}

#[test]
fn test_zetac_large_x() {
    // zetac(x) = 2^-x to double precision, which underflows to 0 for x >= 1075
    for x in [127.0, 130.0, 1000.0, 1074.0, 1076.0, 1e300] {
        assert_eq!(xsf::zetac(x), (-x).exp2());
    }
}
