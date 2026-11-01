#[test]
fn test_iv_ratio() {
    xsref::test("iv_ratio", "d_d-d", |x| xsf::iv_ratio(x[0], x[1]));
    assert_eq!(xsf::iv_ratio(0.5, 3.0), 3.0_f64.tanh());
}

#[test]
fn test_iv_ratio_c() {
    xsref::test("iv_ratio_c", "d_d-d", |x| xsf::iv_ratio_c(x[0], x[1]));
}

#[test]
fn test_iv_ratioinv() {
    // mpmath: x such that besseli(v, x) / besseli(v - 1, x) == r
    let vr = [
        (0.5, 0.5),
        (0.75, 1e-10),
        (1.0, 0.5),
        (2.5, 0.1),
        (10.0, 0.9),
        (100.0, 0.999),
    ];
    let xs = [
        0.549_306_144_334_054_9,
        1.5e-10,
        1.159_319_920_750_138_3,
        0.503_608_676_589_03,
        90.499_984_217_183_92,
        99_450.725_114_819_78,
    ];
    xsref::np_assert_allclose!(vr.map(|(v, r)| xsf::iv_ratioinv(v, r)), xs, rtol = 1e-15);

    assert_eq!(xsf::iv_ratioinv(3.0, 0.0), 0.0);
    assert_eq!(xsf::iv_ratioinv(3.0, 1.0), f64::INFINITY);
    assert!(xsf::iv_ratioinv(0.25, 0.5).is_nan());
    assert!(xsf::iv_ratioinv(f64::INFINITY, 0.5).is_nan());
    assert!(xsf::iv_ratioinv(3.0, 1.5).is_nan());
    assert!(xsf::iv_ratioinv(3.0, f64::NAN).is_nan());
}
