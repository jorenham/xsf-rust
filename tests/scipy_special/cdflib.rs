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
    // the root lies far below `1e-84`, where an absolute tolerance of `1e-100` would limit its
    // precision; the root is from mpmath
    np_assert_allclose!(
        [xsf::fdtridfd(2.0, 1.102_686_716_518_312e-93, 3.0)],
        [1e-95],
        rtol = 1e-13
    );
    // not where fdtr underflows to 0
    assert_eq!(xsf::fdtridfd(1.0, 0.0, 1e-300), 0.0);
    assert!(xsf::fdtridfd(1.0, 0.5, 0.0).is_nan());
    assert!(xsf::fdtridfd(0.0, 0.5, 1.0).is_nan());
}

#[test]
fn test_chndtr() {
    // `test_chndtr_and_inverses_against_wolfram_alpha` and `test_tails` in
    // scipy/special/tests/test_cdflib.py: `(x, df, nc, chndtr)`
    let cases = [
        (25.0, 20.0, 400.0, 4.121_065_511_239_62e-57),
        (0.001, 8.0, 40.0, 5.376_180_620_136_604e-24),
        (20.0, 2.0, 107.0, 1.393_907_435_558_196e-9),
        (3.0, 2.0, 1.0, 0.620_643_653_219_543_6),
        (700.0, 20.0, 400.0, 0.999_999_999_256_806_5),
        (150.0, 13.5, 10.0, 0.999_999_999_999_999_8),
        (0.1, 200.0, 50.0, 1.131_122_486_720_548_1e-299),
        (1e-12, 20.0, 50.0, 3.737_446_313_006_551e-141),
        (9e4, 1e5, 1e3, 1.689_553_370_421_756_6e-141),
    ];
    for (x, df, nc, p) in cases {
        np_assert_allclose!([xsf::chndtr(x, df, nc)], [p], rtol = 1e-13);
    }
    // The central term of the series underflows, where SciPy returns 0, and its exponent
    // `~-1150` limits the accuracy to `~1150 eps`. From Wolfram.
    np_assert_allclose!(
        [xsf::chndtr(20.0, 63.5, 774.0)],
        [1.858_113_891_312_916e-147],
        rtol = 5e-13
    );
    // `test_chndtr_with_nc_zero_equals_chdtr` and `test_chndtr_edge_cases`
    assert_eq!(xsf::chndtr(1.0, 3.0, 0.0), xsf::chdtr(3.0, 1.0));
    assert_eq!(xsf::chndtr(1.0, 0.0, 0.0), 1.0);
    assert_eq!(xsf::chndtr(1.0, f64::INFINITY, f64::INFINITY), 0.0);
    assert_eq!(xsf::chndtr(f64::INFINITY, 1.0, 1.0), 1.0);
    assert_eq!(xsf::chndtr(1e3, 10.0, 10.0), 1.0);
    assert!(xsf::chndtr(1.0, 0.0, 1.0).is_nan());
    assert!(xsf::chndtr(1.0, f64::INFINITY, 1.0).is_nan());
    assert!(xsf::chndtr(f64::INFINITY, 1.0, f64::INFINITY).is_nan());
    assert!(xsf::chndtr(1.0, 1.0, -1.0).is_nan());
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

#[test]
fn test_nctdtridf() {
    // `test_nctdtridf_accuracy` in scipy/special/tests/test_cdflib.py: `(df, nc, t, nctdtr)`
    let cases = [
        (3000.0, 3.0, 0.1, 0.001_865_778_082_632_332_8),
        (0.98, -3.8, 15.0, 0.999_999_026_459_194_5),
        (98.0, 0.38, 0.0015, 0.352_528_178_485_963_13),
        (980.0, 0.38, 1.5, 0.868_424_706_851_729_3),
    ];
    for (df, nc, t, p) in cases {
        // the CDF barely depends on df here, e.g. d ln(p) / d ln(df) ~ 2e-5 for df = 3000
        np_assert_allclose!([xsf::nctdtridf(p, nc, t)], [df], rtol = 1e-8);
    }
    // The CDF peaks near df = 1e-4, and is lost to the same `~17 eps` at df = 1/8, 1 and 8, which
    // also limits the accuracy. The root is from mpmath.
    np_assert_allclose!(
        [xsf::nctdtridf(1e-7, 100.0, 1.0)],
        [0.001_359_433_745_890_397_4],
        rtol = 1e-8
    );
    // the central beta terms of the series underflow; the root is from mpmath
    np_assert_allclose!(
        [xsf::nctdtridf(3.389_564_203_107_224_7e-262, -20.0, -60.0)],
        [5000.0],
        rtol = 1e-12
    );
    assert!(xsf::nctdtridf(0.5, 1.0, 0.0).is_nan());
}

#[test]
fn test_nctdtrinc() {
    // `test_nctdtrinc_gh19896` in scipy/special/tests/test_cdflib.py
    let df = [0.001, 0.98, 9.8, 98.0, 10000.0, 98.0, 0.98, 0.001];
    let p = [
        0.001_000_212_411_985_872_6,
        0.099_999_999_527_969_73,
        0.299_999_996_819_771_9,
        0.799_999_958_571_485_6,
        0.001_005_062_230_783_708_5,
        0.099_999_999_324_025_53,
        0.799_999_998_436_441_5,
        0.998_999_999_551_634_8,
    ];
    let t = [0.0015, 0.15, 1.5, 15.0, 0.0015, 0.15, 15.0, 300.0];
    let nc = [
        3.090_232_306_168_629,
        1.406_141_304_556_198,
        2.014_225_177_124_157,
        13.727_067_118_283_456,
        3.090_232_306_168_629,
        1.431_242_787_793_622_2,
        3.712_743_137_978_295,
        -3.086_951_096_691_082,
    ];
    for i in 0..df.len() {
        np_assert_allclose!([xsf::nctdtrinc(df[i], p[i], t[i])], [nc[i]], rtol = 5e-12);
    }
    // `test_nctdtrinc_zero_nc`
    for (df, t) in [(10.0, 2.0), (100.0, -3.0), (1000.0, 4.0)] {
        let p = xsf::stdtr(df, t);
        np_assert_allclose!([xsf::nctdtrinc(df, p, t)], [0.0], atol = 1e-15);
    }
    assert_eq!(xsf::nctdtrinc(1.0, 0.0, 1.0), f64::INFINITY);
    assert_eq!(xsf::nctdtrinc(1.0, 1.0, 1.0), f64::NEG_INFINITY);
    assert!(xsf::nctdtrinc(0.0, 0.5, 1.0).is_nan());
}
