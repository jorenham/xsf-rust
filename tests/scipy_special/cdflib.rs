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
    // subnormal `x`, whose half can round; from mpmath
    np_assert_allclose!(
        [
            xsf::chndtr(1e-310, 1.0, 1.0),
            xsf::chndtr(f64::from_bits(1), 1.0, 1.0),
            xsf::chndtr(f64::from_bits(3), 1.0, 1.0),
            xsf::chndtr(f64::from_bits(3), 0.5, 1.0),
            xsf::chndtr(f64::from_bits(5), 1e-3, 800.0),
        ],
        [
            4.839_414_490_382_86e-156,
            1.075_685_090_088_338_5e-162,
            1.863_141_228_977_307_2e-162,
            1.104_080_283_333_011_4e-81,
            1.320_925_586_328_177_9e-174,
        ],
        rtol = 1e-15
    );
    // the power series of the central term converges too slowly
    np_assert_allclose!(
        [xsf::chndtr(99_448_457.707_916_78, 1.0, 1e8)],
        [3.650_311_188_203_897e-168],
        rtol = 1e-12
    );
    assert!(xsf::chndtr(1.0, 0.0, 1.0).is_nan());
    assert!(xsf::chndtr(1.0, f64::INFINITY, 1.0).is_nan());
    assert!(xsf::chndtr(f64::INFINITY, 1.0, f64::INFINITY).is_nan());
    assert!(xsf::chndtr(1.0, 1.0, -1.0).is_nan());
}

#[test]
fn test_chndtri_chndtridf_chndtrinc() {
    // `test_chndtr_and_inverses_against_wolfram_alpha` and `test_tails` in
    // scipy/special/tests/test_cdflib.py: `(x, df, nc, chndtr)`
    let cases = [
        (25.0, 20.0, 400.0, 4.121_065_511_239_62e-57),
        (0.001, 8.0, 40.0, 5.376_180_620_136_604e-24),
        (20.0, 2.0, 107.0, 1.393_907_435_558_196e-9),
        (3.0, 2.0, 1.0, 0.620_643_653_219_543_6),
        (350.0, 300.0, 10.0, 0.938_801_280_062_764_1),
        (1e-12, 20.0, 50.0, 3.737_446_313_006_551e-141),
        (9e4, 1e5, 1e3, 1.689_553_370_421_756_6e-141),
    ];
    for (x, df, nc, p) in cases {
        np_assert_allclose!([xsf::chndtri(p, df, nc)], [x], rtol = 1e-13);
        np_assert_allclose!([xsf::chndtridf(x, p, nc)], [df], rtol = 1e-13);
        np_assert_allclose!([xsf::chndtrinc(x, df, p)], [nc], rtol = 1e-13);
    }
    // `test_chndtrix_gh2158`, from R
    np_assert_allclose!(
        [xsf::chndtri(0.999_999, 2.0, 1.0 + 1e-6)],
        [35.257_285_899_505_4],
        rtol = 1e-10
    );
    // a root far below 1e-100, with `d ln(x) / d ln(p) = 2 / df ~ 500`; from mpmath
    np_assert_allclose!(
        [xsf::chndtri(
            0.193_123_859_806_063_84,
            0.004_028_403_030_778_61,
            0.503_188_864_805_273_5
        )],
        [5.426_817_558_374_058e-301],
        rtol = 1e-9
    );
    // a subnormal root, with `d ln(x) / d ln(p) ~ 1000`; from mpmath
    np_assert_allclose!(
        [xsf::chndtri(0.3, 0.001_98, 1.0)],
        [1.702_160_569_678_88e-309],
        rtol = 1e-12
    );
    // the limits where the CDF is smallest or largest
    assert_eq!(xsf::chndtri(0.0, 1.0, 1.0), 0.0);
    assert_eq!(xsf::chndtri(1.0, 1.0, 1.0), f64::INFINITY);
    assert_eq!(xsf::chndtridf(1.0, 0.0, 1.0), f64::INFINITY);
    assert_eq!(xsf::chndtrinc(1.0, 1.0, 1.0), 0.0);
    // p lies above the CDF at df = 0, or at nc = 0
    assert_eq!(xsf::chndtridf(1.0, 0.99, 1.0), 0.0);
    assert_eq!(xsf::chndtrinc(1.0, 1.0, 0.9), 0.0);
    assert!(xsf::chndtri(0.5, 0.0, 1.0).is_nan());
    assert!(xsf::chndtridf(0.0, 0.5, 1.0).is_nan());
    assert!(xsf::chndtrinc(1.0, 1.0, f64::NAN).is_nan());
}

#[test]
fn test_ncfdtr_ncfdtri() {
    // `test_ncfdtr_ncfdtri` in scipy/special/tests/test_cdflib.py: `(dfn, dfd, nc, f, ncfdtr)`
    let cases = [
        (100.0, 0.1, 0.1, 100.0, 0.297_873_964_100_926_76),
        (100.0, 100.0, 0.01, 0.1, 4.434_473_759_869_042_4e-26),
        (100.0, 0.01, 0.1, 0.01, 0.002_848_616_633_080_384),
        (10.0, 0.01, 1.0, 0.1, 0.012_339_557_729_057_956),
        (100.0, 100.0, 0.01, 0.01, 1.892_647_742_096_493_6e-72),
        (1.0, 100.0, 100.0, 0.1, 1.792_594_052_682_130_4e-22),
        (1.0, 0.01, 100.0, 10.0, 0.012_334_711_965_024_968),
        (1.0, 0.01, 10.0, 0.01, 0.000_219_445_252_902_99),
        (10.0, 1.0, 0.1, 100.0, 0.921_934_555_507_070_5),
        (0.1, 0.1, 1.0, 1.0, 0.313_633_581_342_323_9),
        (1.0, 0.1, 100.0, 10.0, 0.029_260_642_796_808_97),
        (1e-100, 3.0, 1.5, 1e100, 0.611_815_287_345_399),
    ];
    for (dfn, dfd, nc, f, p) in cases {
        np_assert_allclose!([xsf::ncfdtr(dfn, dfd, nc, f)], [p], rtol = 1e-13);
        np_assert_allclose!([xsf::ncfdtri(dfn, dfd, nc, p)], [f], rtol = 5e-11);
    }
    assert_eq!(xsf::ncfdtr(100.0, 100.0, 0.1, 10.0), 1.0);
    // the term j = 1 of the series, where it starts, underflows, while the term j = 0, i.e.
    // `(2 / pi) atan(sqrt(f))`, does not; SciPy returns NaN
    np_assert_allclose!(
        [xsf::ncfdtr(1.0, 1.0, 1e-200, 1e-100)],
        [6.366_197_723_675_814e-51],
        rtol = 1e-13
    );
    // e^(-nc / 2) (2 / pi) atan(sqrt(f)) + O(f^(3 / 2)), from Wolfram
    np_assert_allclose!(
        [xsf::ncfdtr(1.0, 1.0, 1.0, 1e-310)],
        [3.861_294_105_202_156e-156],
        rtol = 2e-12
    );
    // `dfn f` and `x` underflow, where the series is its first term; from mpmath
    np_assert_allclose!(
        [xsf::ncfdtr(1.0, 2.0, 1.0, f64::from_bits(1))],
        [9.533_010_901_443_612e-163],
        rtol = 1e-15
    );
    np_assert_allclose!(
        [xsf::ncfdtr(1e-300, 2.0, 1.0, 1e-100)],
        [0.606_530_659_712_633_4],
        rtol = 1e-15
    );
    // `dfn / 2` is subnormal, or `dfn f` is but `x` is not; from mpmath
    np_assert_allclose!(
        [
            xsf::ncfdtr(1e-308, 2.0, 1.0, 1.0),
            xsf::ncfdtr(1e-308, 2.0, 1.0, 1e308),
            xsf::ncfdtr(0.3, 1e-13, 1.0, 1e-320),
        ],
        [
            0.606_530_659_712_633_4,
            0.716_531_310_573_789_3,
            1.504_177_081_709_522e-59
        ],
        rtol = 1e-14
    );
    // `x` underflows, also if `dfn f` or `dfn / (dfd + dfn f)` does not; from mpmath
    np_assert_allclose!(
        [
            xsf::ncfdtr(1.0, 1e100, 1.0, 1e-300),
            xsf::ncfdtr(1e-200, 1e124, 1.0, 1e-120),
        ],
        [4.839_414_490_382_867e-151, 0.606_530_659_712_633_4],
        rtol = 1e-14
    );
    // `x` is subnormal and the noncentral chi-square limit is exact; from mpmath
    np_assert_allclose!(
        [
            xsf::ncfdtr(1.0, 1e308, 1.0, 1e-6),
            xsf::ncfdtr(1.0, 1.7e308, 1.0, 3e-16),
        ],
        [4.839_414_490_382_786e-4, 8.382_111_776_228_172e-9],
        rtol = 1e-14
    );
    // a subnormal root, the exact `tan(pi p / 2)^2`
    np_assert_allclose!(
        [xsf::ncfdtri(1.0, 1.0, 0.0, 1e-155)],
        [2.467_401_100_272_34e-310],
        rtol = 1e-12
    );
    assert_eq!(xsf::ncfdtr(1.0, 1.0, 1.0, 0.0), 0.0);
    assert_eq!(xsf::ncfdtr(1.0, 1.0, 1.0, f64::INFINITY), 1.0);
    assert_eq!(xsf::ncfdtri(1.0, 1.0, 1.0, 0.0), 0.0);
    assert_eq!(xsf::ncfdtri(1.0, 1.0, 1.0, 1.0), f64::INFINITY);
    // `test_ncfdtri_domain_error`
    for (dfn, dfd, nc, p) in [
        (-1.0, 0.1, 0.1, 0.5),
        (1.0, -1.0, 0.1, 0.5),
        (1.0, 1.0, -1.0, 0.5),
        (1.0, 1.0, 1.0, 100.0),
        (1.0, 1.0, 1.0, -1.0),
    ] {
        assert!(xsf::ncfdtri(dfn, dfd, nc, p).is_nan());
    }
    assert!(xsf::ncfdtr(0.0, 1.0, 1.0, 1.0).is_nan());
    assert!(xsf::ncfdtr(1.0, 1.0, -1.0, 1.0).is_nan());
}

#[test]
fn test_ncfdtr_tiny_dfn() {
    // `dfn / 2` rounds to 0, or is lost in `dfn / 2 + icent`, while `I_x(dfn / 2, dfd / 2)` depends
    // on it relatively, e.g. `~dfd / dfn` if both are tiny; from Arb
    np_assert_allclose!(
        [
            xsf::ncfdtr(5e-324, 1e-300, 1.0, 1.0),
            xsf::ncfdtr(5e-324, 2.0, 1.0, f64::MAX),
            xsf::ncfdtr(1e-310, 1e-300, 3.0, 1.0),
            xsf::ncfdtr(1e-80, 1e-290, 2.0, 1e-20),
            xsf::ncfdtr(1e-180, 1e-184, 0.75, 1e12),
            xsf::ncfdtr(0.003, 0.002, 300.0, 0.008),
        ],
        [
            0.606_530_659_712_633_4,
            0.606_530_659_712_633_5,
            0.223_130_160_126_116_8,
            3.678_794_411_714_423e-211,
            6.872_205_567_352_987e-5,
            2.872_768_058_013_942_6e-66,
        ],
        rtol = 1e-13
    );
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
fn test_nctdtr() {
    // `test_gh19896` in scipy/special/tests/test_cdflib.py: `(df, nc, t, nctdtr)`
    let cases = [
        (0.98, -3.8, 0.0015, 0.999_927_998_751_481_5),
        (0.98, -3.8, 0.15, 0.999_952_836_170_050_5),
        (0.98, -3.8, 1.5, 0.999_990_882_301_694_2),
        (0.98, -3.8, 15.0, 0.999_999_026_459_194_5),
        (0.98, 0.38, 0.0015, 0.352_415_331_226_93),
        (0.98, 0.38, 0.15, 0.397_496_972_671_469_83),
        (0.98, 0.38, 1.5, 0.716_862_963_488_558),
        (0.98, 0.38, 15.0, 0.965_624_644_925_749_4),
        (0.98, 3.8, 0.0015, 7.269_733_549_422_93e-05),
        (0.98, 3.8, 0.15, 0.000_124_164_811_475_891_05),
        (0.98, 3.8, 1.5, 0.035_388_035_775_454_095),
        (0.98, 3.8, 15.0, 0.795_482_697_543_058_3),
        (0.98, 38.0, 0.0015, 3.021_069_43e-316),
        (0.98, 38.0, 0.15, 6.069_970_616_996_603e-309),
        (0.98, 38.0, 1.5, 2.591_995_360_483_094e-97),
        (0.98, 38.0, 15.0, 0.011_927_265_886_910_935),
        (9.8, -3.8, 0.0015, 0.999_928_077_619_278_6),
        (9.8, -3.8, 0.15, 0.999_959_941_068_544_2),
        (9.8, -3.8, 1.5, 0.999_999_743_239_478_8),
        (9.8, -3.8, 15.0, 0.999_999_999_999_998_4),
        (9.8, 0.38, 0.0015, 0.352_515_597_910_749_1),
        (9.8, 0.38, 0.15, 0.407_631_201_403_791_94),
        (9.8, 0.38, 1.5, 0.847_679_401_702_465_1),
        (9.8, 0.38, 15.0, 0.999_999_929_711_626_8),
        (9.8, 3.8, 0.0015, 7.277_620_328_149_153e-05),
        (9.8, 3.8, 0.15, 0.000_130_248_022_209_006_52),
        (9.8, 3.8, 1.5, 0.013_477_432_800_072_933),
        (9.8, 3.8, 15.0, 0.999_850_151_230_648),
        (9.8, 38.0, 0.0015, 3.050_660_95e-316),
        (9.8, 38.0, 0.15, 1.790_655_146_76e-313),
        (9.8, 38.0, 1.5, 2.093_594_016_590_074_6e-249),
        (9.8, 38.0, 15.0, 2.252_076_291_604_796e-09),
        (98.0, -3.8, 0.0015, 0.999_928_087_514_910_9),
        (98.0, -3.8, 0.15, 0.999_960_825_017_045_2),
        (98.0, -3.8, 1.5, 0.999_999_930_475_768_2),
        (98.0, -3.8, 15.0, 1.0),
        (98.0, 0.38, 0.0015, 0.352_528_178_485_963_13),
        (98.0, 0.38, 0.15, 0.408_902_530_017_948_46),
        (98.0, 0.38, 1.5, 0.866_467_283_000_655_2),
        (98.0, 0.38, 15.0, 1.0),
        (98.0, 3.8, 0.0015, 7.278_609_891_281_275e-05),
        (98.0, 3.8, 0.15, 0.000_131_031_867_482_700_4),
        (98.0, 3.8, 1.5, 0.010_990_879_189_991_727),
        (98.0, 3.8, 15.0, 0.999_999_999_999_998_9),
        (98.0, 38.0, 0.0015, 3.054_373_85e-316),
        (98.0, 38.0, 0.15, 9.166_833_616_6e-314),
        (98.0, 38.0, 1.5, 1.808_588_423_656_392_6e-288),
        (98.0, 38.0, 15.0, 2.774_053_279_203_590_7e-50),
        (980.0, -3.8, 0.0015, 0.999_928_088_518_896_5),
        (980.0, -3.8, 0.15, 0.999_960_914_455_927_3),
        (980.0, -3.8, 1.5, 0.999_999_941_005_097_9),
        (980.0, -3.8, 15.0, 1.0),
        (980.0, 0.38, 0.0015, 0.352_529_454_879_281_2),
        (980.0, 0.38, 0.15, 0.409_031_532_465_738_2),
        (980.0, 0.38, 1.5, 0.868_424_706_851_729_3),
        (980.0, 0.38, 15.0, 1.0),
        (980.0, 3.8, 0.0015, 7.278_710_289_828_983e-05),
        (980.0, 3.8, 0.15, 0.000_131_111_316_679_065_73),
        (980.0, 3.8, 1.5, 0.010_750_678_886_113_882),
        (980.0, 3.8, 15.0, 1.0),
        (980.0, 38.0, 0.0015, 3.054_750_6e-316),
        (980.0, 38.0, 0.15, 8.619_164_631_3e-314),
        (980.0, 38.0, 15.0, 5.407_535_300_713_606e-105),
    ];
    for (df, nc, t, p) in cases {
        np_assert_allclose!([xsf::nctdtr(df, nc, t)], [p], rtol = 1e-13, atol = 1e-303);
    }
}

#[test]
fn test_nctdtr_accuracy() {
    // `test_gh19896_loose_tolerance`, `test_nctdtr_gh8344` and `test_nctdtr_accuracy`, with one
    // tolerance for the quadrature
    let cases = [
        (980.0, 38.0, 1.5, 1.182_445_411_141_349_3e-291),
        (3000.0, 3.0, 0.1, 0.001_865_778_082_632_332_8),
        (3.0, 5.0, -2.0, 1.564_537_399_914_962_2e-9),
        (1000.0, 10.0, 1.0, 1.149_355_213_382_662_3e-19),
        (1e-5, -6.0, 2.0, 0.999_999_999_013_500_3),
        (10.0, 20.0, 0.15, 6.426_530_505_957_303e-88),
    ];
    for (df, nc, t, p) in cases {
        np_assert_allclose!([xsf::nctdtr(df, nc, t)], [p], rtol = 1e-13);
    }
    // from Wolfram Alpha, where `t < 0 < nc` cancels in the series
    np_assert_allclose!(
        [xsf::nctdtr(10.0, 10.0, -1.0), xsf::nctdtr(1.0, 1.0, -1.0)],
        [2.635_102_038_752_395e-27, 0.057_480_091_794_325_825],
        rtol = 1e-14
    );
    // from mpmath, with extreme `t` or `df`
    let cases = [
        (0.01, 0.0, -3.44e198, 0.005_018_938_223_303_685),
        (0.01, 1.0, -1e100, 0.015_322_457_680_105_562),
        (0.1, 0.0, -1e30, 4.173_803_137_173_216e-4),
        (1e-300, 2.0, -1e300, 0.022_750_131_948_179_22),
        (1e300, 1.0, -1.0, 0.022_750_131_948_179_21),
        (1e300, -3.0, 2.0, 0.999_999_713_348_428_1),
        (1e300, 0.0, 1e300, 1.0),
        (1.0, 0.0, -1e300, 3.183_098_861_837_907e-301),
    ];
    for (df, nc, t, p) in cases {
        np_assert_allclose!([xsf::nctdtr(df, nc, t)], [p], rtol = 1e-13);
    }
    // from Wolfram, with extreme `df` or large `nc`
    let cases = [
        (1e-20, 40.0, 1.0, 1.939_525_002_723_684e-19),
        (1e-20, 1000.0, 1.0, 1.617_606_190_878_827_6e-19),
        (1e-300, 2e6, 1.0, 3.309_370_719_684_12e-298),
        (0.01, 1e4, 1e4, 0.023_373_375_328_696_513),
        (1.0, 1e4, 1e4, 0.317_310_510_282_621_3),
        // `erfc(1 / sqrt(2))`, the limit for `nc = t -> inf`
        (1.0, 1e20, 1e20, 0.317_310_507_862_914_1),
        (1e30, 1e20, 1e20, 0.499_999_999_999_999_8),
        (1.0, -1e308, 1e308, 1.0),
        (f64::MAX, 0.0, 1.0, 0.841_344_746_068_542_9),
    ];
    for (df, nc, t, p) in cases {
        np_assert_allclose!([xsf::nctdtr(df, nc, t)], [p], rtol = 1e-13);
    }
}

#[test]
fn test_nctdtr_edge_cases() {
    // `test_nctdtr` in scipy/special/tests/test_basic.py
    assert_eq!(xsf::nctdtr(1.0, 0.0, 0.0), 0.5);
    assert_eq!(xsf::nctdtr(9.0, 65536.0, 45.0), 0.0);
    np_assert_allclose!([xsf::nctdtr(f64::INFINITY, 1.0, 1.0)], [0.5], atol = 1e-4);
    assert!(xsf::nctdtr(2.0, f64::INFINITY, 10.0).is_nan());
    assert_eq!(xsf::nctdtr(2.0, 1.0, f64::INFINITY), 1.0);
    assert_eq!(xsf::nctdtr(1.0, 1.0, f64::NEG_INFINITY), 0.0);
    assert!(xsf::nctdtr(f64::NAN, 1.0, 1.0).is_nan());
    assert!(xsf::nctdtr(2.0, f64::NAN, 1.0).is_nan());
    assert!(xsf::nctdtr(2.0, 1.0, f64::NAN).is_nan());
    // Unlike in SciPy, not NaN for `t != 0` if `df` is tiny, but close to `Phi(-nc)`
    np_assert_allclose!(
        [
            xsf::nctdtr(1e-300, 1.0, -1.0),
            xsf::nctdtr(1e-300, 1.0, 1.0),
            xsf::nctdtr(5e-324, 1.0, 1.0)
        ],
        [0.158_655_253_931_457_07; 3],
        rtol = 1e-15
    );
}

#[test]
fn test_nctdtrit() {
    // `test_nctdtrit` in scipy/special/tests/test_cdflib.py: `(df, nc, t, nctdtr)`
    let cases = [
        (0.98, 38.0, 1.5, 2.591_995_360_483_094e-97),
        (3000.0, 3.0, 0.1, 0.001_865_778_082_632_332_8),
        (0.98, -3.8, 15.0, 0.999_999_026_459_194_5),
        (9.8, 38.0, 15.0, 2.252_076_291_604_796e-9),
    ];
    for (df, nc, t, p) in cases {
        np_assert_allclose!([xsf::nctdtrit(df, nc, p)], [t], rtol = 1e-10);
    }
    // from mpmath
    np_assert_allclose!(
        [xsf::nctdtrit(1.0, 2.0, 0.5), xsf::nctdtrit(0.01, 0.0, 0.01)],
        [2.809_050_333_659_89, -3.960_440_137_152_441_5e168],
        rtol = 1e-13
    );
}

#[test]
fn test_nctdtrit_edge_cases() {
    let nan = f64::NAN;
    for (df, nc, p) in [(nan, 1.0, 0.5), (1.0, nan, 0.5), (1.0, 1.0, nan)] {
        assert!(xsf::nctdtrit(df, nc, p).is_nan());
    }
    for (df, nc, p) in [(0.0, 1.0, 0.5), (1.0, f64::INFINITY, 0.5), (1.0, 1.0, 1.5)] {
        assert!(xsf::nctdtrit(df, nc, p).is_nan());
    }
    // `nc + ndtri(p)` if `df = inf` or `p = 1`
    assert_eq!(
        xsf::nctdtrit(f64::INFINITY, 2.0, 0.3),
        2.0 + xsf::ndtri(0.3)
    );
    assert_eq!(xsf::nctdtrit(1.0, 1.0, 1.0), f64::INFINITY);
    // Unlike in SciPy, -inf if `p = 0`, instead of inf
    assert_eq!(xsf::nctdtrit(1.0, 1.0, 0.0), f64::NEG_INFINITY);
    // near `t = 0`, where the rounding of `nctdtr` can jump
    for (df, nc, p) in [
        (100.0, 1.0, 0.158_655_253_931_457_1),
        (1.0, 30.0, xsf::ndtr(-30.0)),
    ] {
        assert!(xsf::nctdtrit(df, nc, p).abs() <= 1e-15);
    }
    // `stdtri` if `nc = 0`, also in the tails
    for p in [1e-300, 0.3, 1.0 - f64::EPSILON / 2.0] {
        np_assert_allclose!(
            [xsf::nctdtrit(5.0, 0.0, p)],
            [xsf::stdtri(5.0, p)],
            rtol = 1e-13
        );
    }
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
