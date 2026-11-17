use xsref::np_assert_allclose;

/// NumPy's `TestLogAddExp(2).test_logaddexp(2)_values`: the log of 1 + 5, 2 + 4, ... = 6
#[track_caller]
fn check_logaddexp_values(
    log: fn(f64) -> f64,
    f_f64: fn(f64, f64) -> f64,
    f_f32: fn(f32, f32) -> f32,
) {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0].map(log);
    let y = [5.0, 4.0, 3.0, 2.0, 1.0].map(log);
    let z = [log(6.0); 5];

    let z64: Vec<f64> = x.iter().zip(&y).map(|(&a, &b)| f_f64(a, b)).collect();
    np_assert_allclose!(z64, z, atol = 1e-15);

    let z32: Vec<f64> =
        x.iter().zip(&y).map(|(&a, &b)| f64::from(f_f32(a as f32, b as f32))).collect();
    np_assert_allclose!(z32, z, atol = 1.5e-7);
}

#[test]
fn test_logaddexp_values() {
    check_logaddexp_values(f64::ln, xsf::logaddexp, xsf::logaddexp);
}

#[test]
fn test_logaddexp2_values() {
    check_logaddexp_values(f64::log2, xsf::logaddexp2, xsf::logaddexp2);
}
