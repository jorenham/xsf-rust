use num_complex::c64;

#[test]
fn test_log1p_f64() {
    xsref::test("log1p", "d-d", |x| xsf::log1p(x[0]));
}

#[test]
fn test_log1p_c64() {
    xsref::test("log1p", "cd-cd", |x| xsf::log1p(c64(x[0], x[1])));
}

#[test]
fn test_log1pmx_f64() {
    xsref::test("log1pmx", "d-d", |x| xsf::log1pmx(x[0]));
}

#[test]
fn test_xlogy_f64() {
    xsref::test("xlogy", "d_d-d", |x| xsf::xlogy(x[0], x[1]));
}

#[test]
fn test_xlogy_c64() {
    xsref::test("xlogy", "cd_cd-cd", |x| {
        xsf::xlogy(c64(x[0], x[1]), c64(x[2], x[3]))
    });
}

#[test]
fn test_xlog1py_f64() {
    xsref::test("xlog1py", "d_d-d", |x| xsf::xlog1py(x[0], x[1]));
}

#[test]
fn test_xlog1py_c64() {
    xsref::test("xlog1py", "cd_cd-cd", |x| {
        xsf::xlog1py(c64(x[0], x[1]), c64(x[2], x[3]))
    });
}

#[test]
fn test_log1p_c64_double_double() {
    // |1 + z|^2 = 1 + 2^-50 exactly
    let w = xsf::log1p(c64(-2f64.powi(-25), 2f64.powi(-12)));
    xsref::np_assert_allclose!([w.re], [0.5 * 2f64.powi(-50).ln_1p()], rtol = 1e-15);
}

#[test]
fn test_xlogy_xlog1py_zero() {
    assert_eq!(xsf::xlogy(0.0, f64::INFINITY), 0.0);
    assert_eq!(xsf::xlog1py(0.0, -1.0), 0.0);
    assert_eq!(xsf::xlogy(c64(0.0, 0.0), c64(0.0, 0.0)), c64(0.0, 0.0));
    assert!(xsf::xlogy(0.0, f64::NAN).is_nan());
}
