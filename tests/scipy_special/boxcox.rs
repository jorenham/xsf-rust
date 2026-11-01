fn map2<X1, X2, Y, const N: usize>(f: fn(X1, X2) -> Y, x1: &[X1; N], x2: &[X2; N]) -> [Y; N]
where
    X1: Copy,
    X2: Copy,
{
    core::array::from_fn(|i| f(x1[i], x2[i]))
}

// Ported from `scipy.special.tests.test_boxcox`

#[test]
fn test_boxcox_basic() {
    let x = [0.5, 1.0, 2.0, 4.0];

    // lambda = 0  =>  y = x.ln()
    let y = map2(xsf::boxcox, &x, &[0.0; 4]);
    let expected = x.map(f64::ln);
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);

    // lambda = 1  =>  y = x - 1
    let y = map2(xsf::boxcox, &x, &[1.0; 4]);
    let expected = x.map(|xi| xi - 1.0);
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);

    // lambda = 2  =>  y = 0.5*(x**2 - 1)
    let y = map2(xsf::boxcox, &x, &[2.0; 4]);
    let expected = x.map(|xi| 0.5 * (xi * xi - 1.0));
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);

    // x = 0 and lambda > 0  =>  y = -1 / lambda
    let lam = [0.5, 1.0, 2.0];
    let y = map2(xsf::boxcox, &[0.0; 3], &lam);
    let expected = lam.map(|l| -1.0 / l);
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);
}

#[test]
fn test_boxcox_underflow() {
    let x = 1.0 + 1e-15;
    let lmbda = 1e-306;
    let y = xsf::boxcox(x, lmbda);
    xsref::np_assert_allclose!(&[y], &[x.ln()], rtol = 1e-14);
}

#[test]
fn test_boxcox_nonfinite() {
    // x < 0  =>  y = nan
    let x = [-1.0, -1.0, -0.5];
    let y = map2(xsf::boxcox, &x, &[0.5, 2.0, -1.5]);
    xsref::np_assert_equal!(&y, &[f64::NAN; 3]);

    // x = 0 and lambda <= 0  =>  y = -inf
    let x = 0.0;
    let y = map2(xsf::boxcox, &[x; 2], &[-2.5, 0.0]);
    xsref::np_assert_equal!(&y, &[f64::NEG_INFINITY; 2]);
}

#[test]
fn test_boxcox1p_basic() {
    let x = [-0.25, -1e-20, 0.0, 1e-20, 0.25, 1.0, 3.0];

    // lambda = 0  =>  y = log(1+x)
    let y = map2(xsf::boxcox1p, &x, &[0.0; 7]);
    let expected = x.map(|xi| (1.0 + xi).ln());
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);

    // lambda = 1  =>  y = x
    let y = map2(xsf::boxcox1p, &x, &[1.0; 7]);
    xsref::np_assert_allclose!(&y, &x, atol = 1.5e-7);

    // lambda = 2  =>  y = 0.5*((1+x)**2 - 1) = 0.5*x*(2 + x)
    let y = map2(xsf::boxcox1p, &x, &[2.0; 7]);
    let expected = x.map(|xi| f64::midpoint(xi * xi, 2.0 * xi));
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);

    // x = -1 and lambda > 0  =>  y = -1 / lambda
    let lam = [0.5, 1.0, 2.0];
    let y = map2(xsf::boxcox1p, &[-1.0; 3], &lam);
    let expected = lam.map(|l| -1.0 / l);
    xsref::np_assert_allclose!(&y, &expected, atol = 1.5e-7);
}

#[test]
fn test_boxcox1p_underflow() {
    let x = [1e-15, 1e-306];
    let lmbda = [1e-306, 1e-18];
    let y = map2(xsf::boxcox1p, &x, &lmbda);
    xsref::np_assert_allclose!(&y, &x, rtol = 1e-14);
}

#[test]
fn test_boxcox1p_nonfinite() {
    // x < -1  =>  y = nan
    let x = [-2.0, -2.0, -1.5];
    let y = map2(xsf::boxcox1p, &x, &[0.5, 2.0, -1.5]);
    xsref::np_assert_equal!(&y, &[f64::NAN; 3]);

    // x = -1 and lambda <= 0  =>  y = -inf
    let x = -1.0;
    let y = map2(xsf::boxcox1p, &[x; 2], &[-2.5, 0.0]);
    xsref::np_assert_equal!(&y, &[f64::NEG_INFINITY; 2]);
}

#[test]
fn test_inv_boxcox() {
    let x = [0.0, 1.0, 2.0];
    let lam = [0.0, 1.0, 2.0];
    let y = map2(xsf::boxcox, &x, &lam);
    let x2 = map2(xsf::inv_boxcox, &y, &lam);
    xsref::np_assert_allclose!(&x, &x2, atol = 1.5e-7);

    let x = [0.0, 1.0, 2.0];
    let lam = [0.0, 1.0, 2.0];
    let y = map2(xsf::boxcox1p, &x, &lam);
    let x2 = map2(xsf::inv_boxcox1p, &y, &lam);
    xsref::np_assert_allclose!(&x, &x2, atol = 1.5e-7);
}

#[test]
fn test_inv_boxcox1p_underflow() {
    let x = 1e-15;
    let lam = 1e-306;
    let y = xsf::inv_boxcox1p(x, lam);
    xsref::np_assert_allclose!(&[y], &[x], rtol = 1e-14);

    // the limit for lambda * x -> 0 is expm1(x), not x (scipy/xsf#332)
    let y = xsf::inv_boxcox1p(1.0, 1e-160);
    xsref::np_assert_allclose!(&[y], &[1.0_f64.exp_m1()], rtol = 1e-15);
}

#[test]
fn test_boxcox_premature_overflow() {
    let x = [100.0, 155.0];
    let lmb = [0.01, -155.0];

    // test boxcox & inv_boxcox
    let y = map2(xsf::boxcox, &x, &lmb);
    assert!(y.iter().all(|&v| v.is_finite()));
    let x_inv = map2(xsf::inv_boxcox, &y, &lmb);
    xsref::np_assert_allclose!(&x, &x_inv, rtol = 1e-7, atol = f64::EPSILON);

    // test boxcox1p & inv_boxcox1p
    let x1m = x.map(|v| v - 1.0);
    let y1p = map2(xsf::boxcox1p, &x1m, &lmb);
    assert!(y1p.iter().all(|&v| v.is_finite()));
    let x1p_inv = map2(xsf::inv_boxcox1p, &y1p, &lmb);
    xsref::np_assert_allclose!(&x1m, &x1p_inv, rtol = 1e-7, atol = f64::EPSILON);
}
