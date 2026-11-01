#![allow(clippy::cast_possible_truncation)]

use xsref::np_assert_allclose;

#[test]
fn test_expit_f32() {
    xsref::test("expit", "f-f", |x| xsf::expit(x[0] as f32));
}

#[test]
fn test_expit_f64() {
    xsref::test("expit", "d-d", |x| xsf::expit(x[0]));
}

#[test]
fn test_logit_f32() {
    xsref::test("logit", "f-f", |x| xsf::logit(x[0] as f32));
}

#[test]
fn test_logit_f64() {
    xsref::test("logit", "d-d", |x| xsf::logit(x[0]));
}

#[test]
fn test_log_expit_f32() {
    xsref::test("log_expit", "f-f", |x| xsf::log_expit(x[0] as f32));
}

#[test]
fn test_log_expit_f64() {
    xsref::test("log_expit", "d-d", |x| xsf::log_expit(x[0]));
}

#[test]
fn test_exprel() {
    xsref::test("exprel", "d-d", |x| xsf::exprel(x[0]));
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
