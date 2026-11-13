use xsref::np_assert_allclose;

#[test]
fn test_stdtridf() {
    // `test_students_t_functions_spot_checks` in scipy/special/tests/test_cdflib.py
    let cases = [
        (10.0, 10.0, 0.999_999_205_223_412_2),
        (10.0, -10.0, 7.947_765_877_982_061e-7),
        (1.0, -2000.0, 0.000_159_154_929_828_985_4),
        (1000.0, -2.0, 0.022_885_173_246_625_867),
        (1000.0, 5.0, 0.999_999_661_637_181_8),
        (1000.0, -30.0, 7.687_343_722_021_527e-142),
    ];
    for (df, t, p) in cases {
        np_assert_allclose!([xsf::stdtridf(p, t)], [df], rtol = 1e-9);
    }
    assert!(xsf::stdtridf(0.5, 0.0).is_nan());
    assert!(xsf::stdtridf(1.5, 1.0).is_nan());
    // p is below stdtr(0+, 1) = 1/2, or above stdtr(inf, 1) = ndtr(1)
    assert_eq!(xsf::stdtridf(0.1, 1.0), 0.0);
    assert_eq!(xsf::stdtridf(0.9, 1.0), f64::INFINITY);
    // the limits where stdtr is smallest or largest
    assert_eq!(xsf::stdtridf(0.0, 1.0), 0.0);
    assert_eq!(xsf::stdtridf(1.0, -1.0), 0.0);
    assert_eq!(xsf::stdtridf(0.0, -40.0), f64::INFINITY);
}

#[test]
fn test_fdtridfd() {
    for dfd in [0.5, 1.0, 2.0, 3.0, 100.0] {
        let p = xsf::fdtr(2.0, dfd, 15.0);
        np_assert_allclose!([xsf::fdtridfd(2.0, p, 15.0)], [dfd], rtol = 1e-10);
    }
    // fdtr(30, dfd, 0.3) peaks at dfd ~ 0.3, and has two roots for this p
    for dfd in [0.01, 1.0] {
        let p = xsf::fdtr(30.0, dfd, 0.3);
        let root = xsf::fdtridfd(30.0, p, 0.3);
        np_assert_allclose!([xsf::fdtr(30.0, root, 0.3)], [p], rtol = 1e-12);
    }
    // both roots, ~0.307 and 0.5, lie between the search points 1/8 and 1
    np_assert_allclose!(
        [xsf::fdtridfd(30.0, 0.093_702_302_794_215_2, 0.3)],
        [0.306_842_477_303_175_2],
        rtol = 1e-13
    );
    // not where fdtr underflows to 0
    assert_eq!(xsf::fdtridfd(1.0, 0.0, 1e-300), 0.0);
    assert!(xsf::fdtridfd(1.0, 0.5, 0.0).is_nan());
    assert!(xsf::fdtridfd(0.0, 0.5, 1.0).is_nan());
}
