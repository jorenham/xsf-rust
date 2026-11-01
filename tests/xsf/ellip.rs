#[test]
fn test_ellipk() {
    xsref::test("ellipk", "d-d", |x| xsf::ellipk(x[0]));
}

#[test]
fn test_ellipkm1() {
    xsref::test("ellipkm1", "d-d", |x| xsf::ellipkm1(x[0]));
}

#[test]
fn test_ellipkinc() {
    xsref::test("ellipkinc", "d_d-d", |x| xsf::ellipkinc(x[0], x[1]));
}

#[test]
fn test_ellipkinc_ellipeinc_huge_phi() {
    // phi - npio2 * pi/2 isn't in [-pi/2, pi/2] (reference values from mpmath)
    let (phi, m) = (-1.788_575_894_471_164_1e19, -1.244_349_504_000_301_7e35);
    for (got, expected) in [
        (xsf::ellipkinc(phi, m), -1_348.954_827_564_235_2),
        (xsf::ellipeinc(phi, m), -4.016_600_941_581_635_3e36),
    ] {
        let rel = (got - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "rel = {rel:e}");
    }
    // overflows (was NaN, since the residual's term overflowed with the opposite sign)
    let (phi, m) = (-7.497_579_196_689_92e261, -2.657_819_711_198_385_3e169);
    assert_eq!(xsf::ellipeinc(phi, m), f64::NEG_INFINITY);
    assert_eq!(xsf::ellipeinc(-phi, m), f64::INFINITY);
}

#[test]
fn test_ellipe() {
    xsref::test("ellipe", "d-d", |x| xsf::ellipe(x[0]));
}

#[test]
fn test_ellipeinc() {
    xsref::test("ellipeinc", "d_d-d", |x| xsf::ellipeinc(x[0], x[1]));
}

#[test]
fn test_ellipj() {
    // Workaround for https://github.com/scipy/xsref/issues/11:
    // The xsref table contains incorrect values for the phase φ (am),
    // so we reconstruct it from sn using sn.asin()
    xsref::test("ellipj", "d_d-d_d_d_d", |x| {
        let (sn, cn, dn, am) = xsf::ellipj(x[0], x[1]);
        (sn, cn, dn, am.sin().asin())
    });
}
