use xsref::np_assert_allclose;

#[test]
fn test_softplus_f64() {
    let x = [-100.0, -1.0, 0.0, 1.0, 100.0];
    let y = x.map(|xi: f64| xi.exp().ln_1p());
    np_assert_allclose!(x.map(xsf::softplus), y, atol = 1e-15);
}
