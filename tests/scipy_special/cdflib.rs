use xsref::np_assert_allclose;

#[test]
fn test_pdtrik() {
    // `TestPdtrik` in scipy/special/tests/test_cdflib.py
    let cases = [
        (1.897_610_755_368_228_5e-40, 100.0, 2.0),
        (0.486_701_201_720_851_35, 100.0, 99.0),
        (8.303_834_066_990_52e-69, 1000.0, 500.0),
        (2.252_837_804_125_894e-227, 100_000.0, 90_000.0),
    ];
    for (p, m, k) in cases {
        np_assert_allclose!([xsf::pdtrik(p, m)], [k], rtol = 1e-15);
    }
    for m in [10.0, 50.0, 500.0] {
        np_assert_allclose!([xsf::pdtrik(xsf::pdtr(5.0, m), m)], [5.0], rtol = 1e-15);
    }
    assert_eq!(xsf::pdtrik(0.0, 0.5), 0.0);
    assert_eq!(xsf::pdtrik(0.0, f64::INFINITY), 0.0);
    assert!(xsf::pdtrik(0.5, f64::INFINITY).is_nan());
    assert_eq!(xsf::pdtrik(0.5, 0.0), 0.0);
    // p lies below pdtr(0, m) = exp(-m)
    assert_eq!(xsf::pdtrik(0.95, 1e-6), 0.0);
    assert!(xsf::pdtrik(1.0, 1.0).is_nan());
}

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
fn test_bdtrik() {
    // `test_bdtrik` in scipy/special/tests/test_cdflib.py
    let cases = [
        (6.110_236_824_190_975e-78, 100.0, 0.9, 10.0),
        (4.018_825_542_119_935_7e-224, 1000.0, 0.9, 500.0),
        (0.999_999_831_905_892_8, 1000.0, 0.5, 580.0),
        (0.999_973_438_601_112_4, 100.0, 0.9, 99.0),
    ];
    for (y, n, p, k) in cases {
        np_assert_allclose!([xsf::bdtrik(y, n, p)], [k], rtol = 1e-11);
    }
    // y lies below bdtr(0, 10, 0.5) = 2^-10
    assert_eq!(xsf::bdtrik(1e-4, 10.0, 0.5), 0.0);
    assert_eq!(xsf::bdtrik(0.5, 10.0, 0.0), 0.0);
    assert_eq!(xsf::bdtrik(1.0, 10.0, 0.5), 10.0);
    assert_eq!(xsf::bdtrik(0.5, 10.0, 1.0), 10.0);
    assert!(xsf::bdtrik(0.5, f64::INFINITY, 0.5).is_nan());
}

#[test]
fn test_bdtrin() {
    for (k, n, p) in [(0.0, 1.0, 0.5), (3.0, 10.0, 0.25), (10.0, 3000.0, 0.01)] {
        let y = xsf::bdtr(k, n as i32, p);
        np_assert_allclose!([xsf::bdtrin(k, y, p)], [n], rtol = 1e-12);
    }
    assert_eq!(xsf::bdtrin(3.0, 0.0, 0.5), f64::INFINITY);
    assert_eq!(xsf::bdtrin(3.0, 1.0, 0.5), 3.0);
    assert_eq!(xsf::bdtrin(3.0, 0.0, 0.0), f64::INFINITY);
    assert_eq!(xsf::bdtrin(3.0, 1.0, 1.0), 3.0);
    assert!(xsf::bdtrin(3.0, 0.5, 0.0).is_nan());
    assert!(xsf::bdtrin(3.0, 0.5, 1.0).is_nan());
}

#[test]
fn test_btdtria_btdtrib() {
    // `test_betainc_and_inverses` in scipy/special/tests/test_basic.py: `(a, b, x, p)`
    let cases = [
        (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
        (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
        (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
        (
            0.000_274_279_474_979_266_5,
            289_206.031_25,
            1.639_984_034_231_756e-56,
            0.968_870_878_219_604_5,
        ),
        (4.0, 99997.0, 0.000_194_784_157_889_212_1, 0.999_995),
    ];
    for (a, b, x, p) in cases {
        np_assert_allclose!([xsf::btdtria(p, b, x)], [a], rtol = 1e-13);
        np_assert_allclose!([xsf::btdtrib(a, p, x)], [b], rtol = 1e-13);
    }
    assert_eq!(xsf::btdtria(0.0, 2.0, 0.5), f64::INFINITY);
    assert_eq!(xsf::btdtria(1.0, 2.0, 0.5), 0.0);
    assert_eq!(xsf::btdtrib(2.0, 0.0, 0.5), 0.0);
    assert_eq!(xsf::btdtrib(2.0, 1.0, 0.5), f64::INFINITY);
    // the limits where the CDF is 1 (or 0 for `btdtrib` with `a = inf`) for all `a` or `b`
    assert_eq!(xsf::btdtria(0.0, 2.0, 1.0), f64::INFINITY);
    assert_eq!(xsf::btdtria(1.0, f64::INFINITY, 0.5), 0.0);
    assert_eq!(xsf::btdtrib(2.0, 1.0, 1.0), f64::INFINITY);
    assert_eq!(xsf::btdtrib(f64::INFINITY, 0.0, 0.5), 0.0);
    assert!(xsf::btdtria(0.5, 2.0, 1.0).is_nan());
    assert!(xsf::btdtria(0.5, f64::INFINITY, 0.5).is_nan());
    assert!(xsf::btdtria(0.0, 2.0, 1.5).is_nan());
    assert!(xsf::btdtrib(2.0, 0.5, 1.0).is_nan());
    assert!(xsf::btdtrib(0.0, 0.5, 0.5).is_nan());
    // a tiny root, by the symmetry `betainc(a, a, 1/2) = 1/2`
    np_assert_allclose!([xsf::btdtria(0.5, 1e-200, 0.5)], [1e-200], rtol = 1e-13);
    np_assert_allclose!([xsf::btdtrib(1e-200, 0.5, 0.5)], [1e-200], rtol = 1e-13);
}

#[test]
fn test_chdtriv() {
    for (v, x) in [(0.5, 0.1), (1.0, 1.0), (3.0, 10.0), (100.0, 80.0)] {
        let p = xsf::chdtr(v, x);
        np_assert_allclose!([xsf::chdtriv(p, x)], [v], rtol = 1e-13);
    }
    // `test_chdtriv_p_equals_1_returns_0` in scipy/special/tests/test_cdflib.py
    assert_eq!(xsf::chdtriv(1.0, 0.1), 0.0);
    assert_eq!(xsf::chdtriv(1.0, 100.0), 0.0);
    assert_eq!(xsf::chdtriv(0.0, 1.0), f64::INFINITY);
    assert!(xsf::chdtriv(0.5, f64::INFINITY).is_nan());
    assert!(xsf::chdtriv(0.5, 0.0).is_nan());
    // the smallest subnormal, which `x/2` loses, and the next one, with a root from mpmath
    assert!(xsf::chdtriv(0.5, f64::from_bits(1)).is_nan());
    np_assert_allclose!(
        [xsf::chdtriv(0.5, f64::from_bits(2))],
        [0.001_863_640_484_035_382_4],
        rtol = 1e-13
    );
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
fn test_nbdtrik_nbdtrin() {
    // `test_inverses` in scipy/special/tests/test_cdflib.py: `(k, y, p, n)`
    let cases = [
        (5.0, 0.321_456_909_179_687_5, 0.25, 3.0),
        (10.0, 3.820_379_007_878_081e-70, 0.15, 100.0),
        (100.0, 0.975_574_722_927_043, 0.15, 10.0),
        (2.0, 3.503_253_250_000_002e-73, 0.001, 25.0),
        (2.0, 0.999_999_327_102_995_4, 0.999, 15.0),
        (500.0, 0.999_999_999_860_317_8, 0.1, 15.0),
        (0.0, 1.000_000_000_000_000_6e-10, 0.1, 10.0),
    ];
    for (k, y, p, n) in cases {
        np_assert_allclose!([xsf::nbdtrin(k, y, p)], [n], rtol = 1e-8);
        np_assert_allclose!([xsf::nbdtrik(y, n, p)], [k], rtol = 1e-8, atol = 1e-15);
    }
    // y lies below nbdtr(0, 10, 0.5) = 2^-10
    assert_eq!(xsf::nbdtrik(1e-4, 10.0, 0.5), 0.0);
    assert_eq!(xsf::nbdtrik(1.0, 10.0, 0.5), f64::INFINITY);
    assert_eq!(xsf::nbdtrik(0.5, 10.0, 0.0), f64::INFINITY);
    assert!(xsf::nbdtrik(0.5, 0.0, 0.5).is_nan());
    assert_eq!(xsf::nbdtrin(3.0, 0.0, 0.5), f64::INFINITY);
    assert_eq!(xsf::nbdtrin(3.0, 1.0, 0.5), 0.0);
    assert!(xsf::nbdtrin(3.0, 0.5, 1.0).is_nan());
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
