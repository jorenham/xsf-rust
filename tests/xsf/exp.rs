use num_complex::c64;

xsref::tests! {
    test_expm1_f64: "expm1", "d-d", |x| xsf::expm1(x[0]);
    test_expm1_c64: "expm1", "cd-cd", |x| xsf::expm1(c64(x[0], x[1]));
    test_exp2_f64: "exp2", "d-d", |x| xsf::exp2(x[0]);
    test_exp10_f64: "exp10", "d-d", |x| xsf::exp10(x[0]);
}

#[test]
fn test_expm1_c64_edges() {
    // small |z|, where exp(z) - 1 would lose precision
    let w = xsf::expm1(c64(1e-10, 1e-10));
    xsref::np_assert_allclose!([w.re, w.im], [1e-10, 1e-10 + 1e-20], rtol = 1e-15);

    // each branch: re <= -40, -40 < re <= -1, and re > -1
    for re in [-41.0_f64, -40.0, -1.5, -1.0, -0.25, 2.0] {
        let w = xsf::expm1(c64(re, 2.0));
        let expected = [re.exp() * 2.0_f64.cos() - 1.0, re.exp() * 2.0_f64.sin()];
        xsref::np_assert_allclose!([w.re, w.im], expected, rtol = 1e-14);
    }
    assert_eq!(xsf::expm1(c64(-40.0, 1.0)).re, -1.0);
}
