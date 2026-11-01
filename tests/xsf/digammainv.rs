#[test]
fn test_digammainv() {
    // mpmath: x such that digamma(x) == y
    let ys = [0.0, 1.0, -1.0, -10.0, 5.0, 700.0];
    let xs = [
        1.461_632_144_968_362_2,
        3.203_171_468_376_931_2,
        0.785_003_325_374_520_7,
        0.104_357_198_770_116_51,
        148.912_878_356_218_87,
        1.014_232_054_735_004_5e304,
    ];
    xsref::np_assert_allclose!(ys.map(xsf::digammainv), xs, rtol = 1e-14);

    assert!(xsf::digammainv(f64::NAN).is_nan());
    assert_eq!(xsf::digammainv(f64::NEG_INFINITY), 0.0);
    assert_eq!(xsf::digammainv(710.0), f64::INFINITY);
}

#[test]
fn test_digammainv_newton_cycle() {
    // unlike in xsf (scipy/xsf#361): Newton cycles here, which used to return NaN
    let x = xsf::digammainv(9.220_027_678_716_601);
    xsref::np_assert_allclose!([x], [10_097.843_801_671_605], rtol = 1e-14);
}
