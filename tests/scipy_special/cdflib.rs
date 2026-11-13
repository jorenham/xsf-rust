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

#[test]
fn test_ncfdtrinc() {
    // `test_ncfdtr_ncfdtri` in scipy/special/tests/test_cdflib.py: `(dfn, dfd, nc, f, ncfdtr)`
    let cases = [
        (100.0, 0.1, 0.1, 100.0, 0.297_873_964_100_926_76),
        (100.0, 100.0, 0.01, 0.1, 4.434_473_759_869_042_4e-26),
        (10.0, 0.01, 1.0, 0.1, 0.012_339_557_729_057_956),
        (1.0, 100.0, 100.0, 0.1, 1.792_594_052_682_130_4e-22),
        (1.0, 0.01, 10.0, 0.01, 0.000_219_445_252_902_99),
        (0.1, 0.1, 1.0, 1.0, 0.313_633_581_342_323_9),
    ];
    for (dfn, dfd, nc, f, p) in cases {
        np_assert_allclose!([xsf::ncfdtrinc(dfn, dfd, p, f)], [nc], rtol = 1e-10);
    }
    // `test_ncfdtrinc_zero_nc`
    for (dfn, dfd, f) in [(1.0, 10.0, 2.0), (3.0, 5.0, 1.0), (10.0, 10.0, 4.0)] {
        let p = xsf::fdtr(dfn, dfd, f);
        np_assert_allclose!([xsf::ncfdtrinc(dfn, dfd, p, f)], [0.0], atol = 1e-15);
    }
    // ncfdtr(2, 2, nc, 1) = exp(-nc / 4) / 2, so that d ln(p) / d ln(nc) ~ 1e-11 here
    np_assert_allclose!(
        [xsf::ncfdtrinc(2.0, 2.0, 0.499_999_999_993_75, 1.0)],
        [5.000_000_413_733_105e-11],
        rtol = 1e-5
    );
    // the central incomplete beta function of the series underflows; the root is from mpmath
    np_assert_allclose!(
        [xsf::ncfdtrinc(0.001, 0.001, 0.03, 1e-300)],
        [4.936_046_727_488_337],
        rtol = 1e-12
    );
    assert_eq!(xsf::ncfdtrinc(1.0, 1.0, 0.0, 1.0), f64::INFINITY);
    assert_eq!(xsf::ncfdtrinc(1.0, 1.0, 1.0, 1.0), 0.0);
    assert!(xsf::ncfdtrinc(1.0, 1.0, 0.5, 0.0).is_nan());
}

#[test]
fn test_ncfdtridfd_ncfdtridfn() {
    // `test_ncfdtridfd` and `test_ncfdtridfn` in scipy/special/tests/test_basic.py, with
    // `ncfdtr(2, dfd, 0.25, 15)` and `ncfdtr(dfn, 2, 0.25, 15)` from SciPy
    let cases = [
        (1.0, 0.809_713_927_482_644_6),
        (2.0, 0.930_204_317_118_978_3),
        (3.0, 0.967_878_689_542_765_6),
    ];
    for (dfd, p) in cases {
        np_assert_allclose!([xsf::ncfdtridfd(2.0, p, 0.25, 15.0)], [dfd], rtol = 1e-12);
    }
    let cases = [
        (0.1, 0.892_442_257_903_718_6),
        (1.0, 0.925_623_741_271_307_1),
        (2.0, 0.930_204_317_118_978_3),
        (3.0, 0.931_884_110_629_833_7),
    ];
    for (dfn, p) in cases {
        np_assert_allclose!([xsf::ncfdtridfn(p, 2.0, 0.25, 15.0)], [dfn], rtol = 1e-12);
    }
    // `dfn f` overflows, and the subnormal `1 - x` limits the accuracy; the root is from mpmath
    np_assert_allclose!(
        [xsf::ncfdtridfd(1e6, 0.5, 1.0, 1e308)],
        [0.001_937_357_577_091_185_3],
        rtol = 1e-9
    );
    assert!(xsf::ncfdtridfd(2.0, 0.5, f64::NAN, 15.0).is_nan());
    assert!(xsf::ncfdtridfn(0.5, 2.0, -1.0, 15.0).is_nan());
}
