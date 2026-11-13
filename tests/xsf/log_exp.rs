use xsref::np_assert_allclose;

xsref::tests! {
    test_expit_f32: "expit", "f-f", |x| xsf::expit(x[0] as f32);
    test_expit_f64: "expit", "d-d", |x| xsf::expit(x[0]);
    test_logit_f32: "logit", "f-f", |x| xsf::logit(x[0] as f32);
    test_logit_f64: "logit", "d-d", |x| xsf::logit(x[0]);
    test_log_expit_f32: "log_expit", "f-f", |x| xsf::log_expit(x[0] as f32);
    test_log_expit_f64: "log_expit", "d-d", |x| xsf::log_expit(x[0]);
    test_exprel: "exprel", "d-d", |x| xsf::exprel(x[0]);
}

#[test]
fn test_expit_subnormal() {
    // mpmath
    let x = [-709.0, -710.0, -730.0, -745.0, -746.0];
    let desired = [
        1.216_780_750_623_423e-308,
        4.476_286_225_675_13e-309,
        9.226_315e-318,
        5e-324,
        0.0,
    ];
    np_assert_allclose!(x.map(xsf::expit), desired, rtol = 1e-15, atol = 5e-324);
    assert_eq!(xsf::expit(f64::NEG_INFINITY), 0.0);
    assert!((xsf::expit(-89.0_f32) / 2.227_364e-39 - 1.0).abs() < 1e-6);
}

#[test]
fn test_exprel_near_overflow() {
    // mpmath
    let x = [709.0, 710.0, 715.0, 716.0];
    let desired = [
        1.159_154_790_064_171e305,
        3.146_471_501_636_212_5e305,
        4.637_121_967_338_202e307,
        1.258_739_962_544_279_3e308,
    ];
    np_assert_allclose!(x.map(xsf::exprel), desired, rtol = 1e-15);
    assert_eq!(
        [716.4, 717.0, f64::INFINITY].map(xsf::exprel),
        [f64::INFINITY; 3]
    );
    assert!((xsf::exprel(90.0_f32) / 1.356_003_7e37 - 1.0).abs() < 1e-6);
    assert_eq!([94.0, f32::INFINITY].map(xsf::exprel), [f32::INFINITY; 2]);
}

#[test]
fn test_log1mexp() {
    let xs: [f64; 5] = [-20.0, -5.0, -1.5, -1.0, -0.75];
    let expected = xs.map(|x| (-x.exp_m1()).ln());
    np_assert_allclose!(xs.map(xsf::log1mexp), expected, atol = 1e-15);

    assert!(xsf::log1mexp(0.5_f64).is_nan());

    let zero_limit = xsf::log1mexp(0.0_f64);
    assert!(zero_limit.is_infinite() && zero_limit.is_sign_negative());
}
