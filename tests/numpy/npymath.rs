use xsref::np_assert_allclose;

#[test]
fn test_logaddexp_f32() {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [5.0, 4.0, 3.0, 2.0, 1.0];
    let z = [6.0, 6.0, 6.0, 6.0, 6.0];

    let xf = x.map(f64::log2);
    let yf = y.map(f64::log2);
    let zf = z.map(f64::log2);

    #[allow(clippy::cast_possible_truncation)]
    let zr: [f32; 5] = std::array::from_fn(|i| xsf::logaddexp2(xf[i] as f32, yf[i] as f32));
    np_assert_allclose!(zr.map(f64::from), zf, atol = 1.5e-7);
}

#[test]
fn test_logaddexp_f64() {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [5.0, 4.0, 3.0, 2.0, 1.0];
    let z = [6.0, 6.0, 6.0, 6.0, 6.0];

    let xf = x.map(f64::log2);
    let yf = y.map(f64::log2);
    let zf = z.map(f64::log2);

    let zr: [f64; 5] = std::array::from_fn(|i| xsf::logaddexp2(xf[i], yf[i]));
    np_assert_allclose!(zr, zf, atol = 1e-15);
}
