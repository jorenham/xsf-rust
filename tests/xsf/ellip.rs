xsref::tests! {
    test_ellipk: "ellipk", "d-d", |x| xsf::ellipk(x[0]);
    test_ellipkm1: "ellipkm1", "d-d", |x| xsf::ellipkm1(x[0]);
    test_ellipkinc: "ellipkinc", "d_d-d", |x| xsf::ellipkinc(x[0], x[1]);
}

#[test]
fn test_ellipkinc_ellipeinc_huge_phi() {
    // the large-phi limit (reference values from mpmath)
    let (phi, m) = (-1.788_575_894_471_164_1e19, -1.244_349_504_000_301_7e35);
    let (phi2, m2) = (3.560_108_682_609_496e16, -1.24e35);
    for (got, expected) in [
        (xsf::ellipkinc(phi, m), -1_348.954_827_564_235_2),
        (xsf::ellipeinc(phi, m), -4.016_600_941_581_635_3e36),
        (xsf::ellipkinc(phi2, m2), 2.689_647_972_067_577),
        (xsf::ellipeinc(phi2, m2), 7.980_943_234_086_402e33),
    ] {
        let rel = (got - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "rel = {rel:e}");
    }
    // overflows (was NaN)
    let (phi, m) = (-7.497_579_196_689_92e261, -2.657_819_711_198_385_3e169);
    assert_eq!(xsf::ellipeinc(phi, m), f64::NEG_INFINITY);
    assert_eq!(xsf::ellipeinc(-phi, m), f64::INFINITY);
}

#[test]
fn test_ellipkinc_m1_pi2() {
    // fl(pi/2) < pi/2, so F(fl(pi/2) | 1) is finite (reference value from mpmath and Arb)
    let expected = 38.025_003_373_828_866;
    for s in [1.0, -1.0] {
        let got = xsf::ellipkinc(s * core::f64::consts::FRAC_PI_2, 1.0);
        let rel = (got - s * expected).abs() / expected;
        assert!(rel < 1e-15, "{got} != {}", s * expected);
    }
    // but infinite just above
    let above = f64::from_bits(core::f64::consts::FRAC_PI_2.to_bits() + 1);
    assert_eq!(xsf::ellipkinc(above, 1.0), f64::INFINITY);
}

xsref::tests! {
    test_ellipe: "ellipe", "d-d", |x| xsf::ellipe(x[0]);
    test_ellipeinc: "ellipeinc", "d_d-d", |x| xsf::ellipeinc(x[0], x[1]);

    // Workaround for https://github.com/scipy/xsref/issues/11:
    // The xsref table contains incorrect values for the phase φ (am),
    // so we reconstruct it from sn using sn.asin()
    test_ellipj: "ellipj", "d_d-d_d_d_d", |x| {
        let (sn, cn, dn, am) = xsf::ellipj(x[0], x[1]);
        (sn, cn, dn, am.sin().asin())
    };
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
            1e-12,
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
            // between K and 2K
            21.5988,
            0.999_999_999_999,
            [
                0.999_999_954_996_353_5,
                -3.000_121_513_706_456e-4,
                3.000_138_179_281_605_5e-4,
                1.571_096_338_950_767_8,
            ],
            1e-14,
        ),
        (
            // near -21 K
            -407.508_951_116_008_1,
            0.999_999_999_999_999_8,
            [
                -1.0,
                9.183_585_813_193_36e-10,
                1.492_943_359_303_742e-8,
                -32.986_722_861_774_47,
            ],
            1e-11,
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

    // dn near the quarter period, where 1 - m sn^2 cancels (reference value from mpmath)
    let (_, _, dn, _) = xsf::ellipj(12.706_487_231_104_472, 0.999_999_999_85);
    let expected = 1.224_806_159_815_916_4e-5;
    assert!(
        (dn - expected).abs() < 1e-13 * expected,
        "{dn} != {expected}"
    );

    // too large |u| to reduce by quarter periods
    for u in [-1e18, 1e17, 1e18] {
        let (sn, cn, dn, _) = xsf::ellipj(u, 0.999_999_999_9);
        assert!(
            sn.abs() <= 1.0 && cn.abs() <= 1.0 && dn.abs() <= 1.0,
            "ellipj({u})"
        );
    }
}

#[test]
fn test_ellipj_am_small_u() {
    // m = 1 or close to 1, where 2 atan(e^u) - pi/2 cancels (reference values from mpmath)
    for (u, m, expected) in [
        (1e-30, 1.0, 1e-30),
        (1.29e-15, 1.0, 1.29e-15),
        (-1.27e-10, 1.0, -1.27e-10),
        (1e-30, 0.999_999_999_95, 1e-30),
        (1e-5, 0.999_999_999_95, 9.999_999_999_833_334e-6),
    ] {
        let (_, _, _, ph) = xsf::ellipj(u, m);
        let rel = (ph - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "ellipj({u}, {m}): {ph} != {expected}");
    }
}
