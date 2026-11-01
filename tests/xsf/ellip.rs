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

#[test]
fn test_ellipj_near_1() {
    // m close to 1 with large |u| (reference values from mpmath)
    for (u, m, expected, rtol) in [
        (
            -607.281_543_046_733_3,
            0.999_999_999_999_475_6,
            [
                0.999_999_999_997_620_9,
                2.181_313_477_196_956_8e-6,
                2.298_366_119_731_372e-6,
                -61.261_058_926_314_45,
            ],
            1e-8,
        ),
        (
            30.0,
            0.999_999_999_999_9,
            [
                0.991_112_190_818_886_4,
                -0.133_028_663_115_087_02,
                0.133_028_663_115_456_34,
                1.704_220_507_506_328_7,
            ],
            1e-14,
        ),
        (
            // just above the cutoff of the expansion in 1 - m, which has a rel. error of 4.5e-8 here
            21.5988,
            0.999_999_999_999,
            [
                0.999_999_954_996_353_5,
                -3.000_121_513_706_456e-4,
                3.000_138_179_281_605_5e-4,
                1.571_096_338_950_767_8,
            ],
            1e-9,
        ),
        (
            -400.0,
            1.0,
            [
                -1.0,
                3.830_339_193_428_011_4e-174,
                3.830_339_193_428_011_4e-174,
                -core::f64::consts::FRAC_PI_2,
            ],
            1e-14,
        ),
        (
            720.0,
            1.0,
            [
                1.0,
                4.064_461_604_84e-313,
                4.064_461_604_84e-313,
                core::f64::consts::FRAC_PI_2,
            ],
            1e-10,
        ),
    ] {
        let (sn, cn, dn, ph) = xsf::ellipj(u, m);
        for (got, expected) in [sn, cn, dn, ph].into_iter().zip(expected) {
            let rel = (got - expected).abs() / expected.abs();
            assert!(
                rel < rtol,
                "ellipj({u}, {m}): {got} != {expected} (rel = {rel:e})"
            );
        }
    }
    // the smallest subnormal sech(u)
    let (_, cn, dn, _) = xsf::ellipj(-745.5, 1.0);
    assert_eq!((cn, dn), (f64::from_bits(1), f64::from_bits(1)));
}
