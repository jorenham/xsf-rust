use core::f64::consts::PI;

#[test]
fn test_bivariate_normal_cdf() {
    // mpmath: the integral of npdf(x) * ncdf((k - r*x) / sqrt(1 - r^2)) over x <= h
    let hkr = [
        (1.0, 2.0, 0.2),
        (-0.5, 0.3, -0.6),
        (2.0, -1.0, 0.8),
        (0.4, 0.7, -0.95),
        (-1.5, -2.0, 0.99),
        (3.0, 3.0, -0.999),
    ];
    let expected = [
        0.825_336_189_341_017_8,
        0.101_592_099_747_496_54,
        0.158_655_206_571_701_87,
        0.413_466_121_994_546_5,
        0.022_749_515_662_487_264,
        0.997_300_203_936_739_8,
    ];
    let actual = hkr.map(|(h, k, r)| xsf::bivariate_normal_cdf(h, k, r));
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);

    // orthant probability
    for r in [-1.0, -0.8, -0.5, 0.0, 0.5, 0.95, 1.0] {
        let p = xsf::bivariate_normal_cdf(0.0, 0.0, r);
        xsref::np_assert_allclose!([p], [0.25 + r.asin() / (2.0 * PI)], atol = 1e-15);
    }

    let inf = f64::INFINITY;
    assert_eq!(xsf::bivariate_normal_cdf(-inf, 0.0, 0.5), 0.0);
    assert_eq!(xsf::bivariate_normal_cdf(inf, inf, 0.5), 1.0);
    assert_eq!(xsf::bivariate_normal_cdf(inf, 1.0, 0.5), xsf::ndtr(1.0));
    assert!(xsf::bivariate_normal_cdf(f64::NAN, 0.0, 0.5).is_nan());
    assert!(xsf::bivariate_normal_cdf(1.0, -1.0, f64::NAN).is_nan());
    assert!(xsf::bivariate_normal_cdf(f64::NAN, -41.0, 0.5).is_nan());
    // NaN in xsf
    assert_eq!(xsf::bivariate_normal_cdf(1e70, 38.0, 0.97), 1.0);
    assert_eq!(xsf::bivariate_normal_cdf(-1e200, -1e200, 0.5), 0.0);
}

#[test]
fn test_trivariate_normal_cdf() {
    // P(X < 0, Y < 0, Z < 0) = 1 / 4 for the equicorrelation matrix with rho = 1/2
    let p = xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1e-7);
    xsref::np_assert_allclose!([p], [0.25], atol = 1e-15);

    // mpmath: the integral of npdf(x) * bivariate_normal_cdf(...) over x <= h1
    let args = [
        [0.5, -0.3, 1.2, 0.3, -0.2, 0.5],
        [-1.0, 0.0, 1.0, 0.9, 0.8, 0.7],
        [1.0, 2.0, -0.5, -0.4, 0.6, -0.1],
    ];
    let expected = [
        0.292_929_637_622_110_7,
        0.157_894_486_614_603_12,
        0.292_319_806_202_785_64,
    ];
    let actual = args.map(|[h1, h2, h3, r12, r13, r23]| {
        xsf::trivariate_normal_cdf(h1, h2, h3, r12, r13, r23, 1e-14)
    });
    xsref::np_assert_allclose!(actual, expected, rtol = 1e-14);
    // permuting the variables describes the same event
    let [h1, h2, h3, r12, r13, r23] = args[0];
    let p = xsf::trivariate_normal_cdf(h3, h1, h2, r13, r23, r12, 1e-14);
    xsref::np_assert_allclose!([p], [expected[0]], rtol = 1e-14);

    // independent variables
    let p = xsf::trivariate_normal_cdf(1.0, -0.5, 0.25, 0.0, 0.0, 0.0, 1e-7);
    let p_ref = xsf::ndtr(1.0) * xsf::ndtr(-0.5) * xsf::ndtr(0.25);
    xsref::np_assert_allclose!([p], [p_ref], rtol = 1e-15);

    let inf = f64::INFINITY;
    assert_eq!(
        xsf::trivariate_normal_cdf(-inf, 0.0, 0.0, 0.2, -0.3, 0.4, 1e-7),
        0.0
    );
    assert_eq!(
        xsf::trivariate_normal_cdf(inf, inf, inf, 0.2, -0.3, 0.4, 1e-7),
        1.0
    );
    assert_eq!(
        xsf::trivariate_normal_cdf(1.0, inf, 2.0, 0.2, -0.3, 0.4, 1e-7),
        xsf::bivariate_normal_cdf(1.0, 2.0, -0.3)
    );
    // invalid correlations
    assert!(xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, 1.1, 0.0, 0.0, 1e-7).is_nan());
    assert!(xsf::trivariate_normal_cdf(0.0, 0.0, 0.0, 0.9, 0.9, -0.9, 1e-7).is_nan());
    assert!(xsf::trivariate_normal_cdf(f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-7).is_nan());
}
