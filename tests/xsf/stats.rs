#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

use num_complex::c64;

#[test]
fn test_stdtr() {
    const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
    const X: [f64; 9] = [-8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0];
    // values from scipy.special.stdtr
    const P_REF: [[f64; 9]; 6] = [
        [
            0.113_252_546_403_789,
            0.159_610_041_494_336,
            0.222_757_445_091_566,
            0.301_121_610_841_322,
            0.500_000_000_000_000,
            0.698_878_389_158_678,
            0.777_242_554_908_434,
            0.840_389_958_505_664,
            0.886_747_453_596_211,
        ],
        [
            0.039_583_424_160_566,
            0.077_979_130_377_369,
            0.147_583_617_650_433,
            0.250_000_000_000_000,
            0.500_000_000_000_000,
            0.750_000_000_000_000,
            0.852_416_382_349_567,
            0.922_020_869_622_631,
            0.960_416_575_839_434,
        ],
        [
            0.007_634_036_082_669,
            0.028_595_479_208_968,
            0.091_751_709_536_137,
            0.211_324_865_405_187,
            0.500_000_000_000_000,
            0.788_675_134_594_813,
            0.908_248_290_463_863,
            0.971_404_520_791_032,
            0.992_365_963_917_331,
        ],
        [
            0.002_038_288_793_893,
            0.014_004_228_005_073,
            0.069_662_984_279_422,
            0.195_501_109_477_885,
            0.500_000_000_000_000,
            0.804_498_890_522_115,
            0.930_337_015_720_578,
            0.985_995_771_994_927,
            0.997_961_711_206_107,
        ],
        [
            0.000_661_948_454_609,
            0.008_065_044_950_046,
            0.058_058_261_758_408,
            0.186_950_483_150_030,
            0.500_000_000_000_000,
            0.813_049_516_849_970,
            0.941_941_738_241_592,
            0.991_934_955_049_954,
            0.999_338_051_545_391,
        ],
        [
            0.000_005_887_471_395,
            0.001_259_166_312_368,
            0.036_694_017_385_370,
            0.170_446_566_151_030,
            0.500_000_000_000_000,
            0.829_553_433_848_970,
            0.963_305_982_614_630,
            0.998_740_833_687_632,
            0.999_994_112_528_605,
        ],
    ];
    for (i, &nu) in NU.iter().enumerate() {
        let p = X.map(|x| xsf::stdtr(nu, x));
        xsref::np_assert_allclose!(p, P_REF[i], atol = 1e-15);
    }
}

#[test]
fn test_stdtr_extreme() {
    // reference values from mpmath
    for &(nu, x, p) in &[
        // x^2 << nu
        (7.8e5, 5.9e-8, 0.500_000_023_537_587),
        (7.8e5, -5.9e-8, 0.499_999_976_462_413),
        // large nu and 1 < x^2 <= 2
        (1e9, 1.1, 0.864_333_938_921_217_7),
        (1e8, -1.3, 0.096_800_486_083_800_24),
        // huge nu (ndtr)
        (1e200, -2.33, 0.009_903_075_559_164_25),
        // huge |x|
        (0.5, -1e200, 3.207_009_754_142_229e-101),
        (0.0016, -3.8e274, 0.180_685_979_862_586_23),
        (0.0016, 3.8e274, 0.819_314_020_137_413_8),
        (1e306, -1e200, 0.0),
        (1e306, 1e200, 1.0),
        (5e-324, f64::NEG_INFINITY, 0.0),
        (5e-324, f64::INFINITY, 1.0),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtr(nu, x)], &[p], rtol = 1e-13);
    }
}

#[test]
fn test_stdtri() {
    const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
    const P: [f64; 7] = [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999];
    // values from Wolfram Alpha
    const X_REF: [[f64; 7]; 6] = [
        [
            -102_849.115_630_175_55,
            -1_028.491_010_471_621_9,
            -10.270_324_410_234_506,
            0.0,
            10.270_324_410_234_506,
            1_028.491_010_471_621_9,
            102_849.115_630_175_55,
        ],
        [
            -318.308_838_985_550_45,
            -31.820_515_953_773_96,
            -3.077_683_537_175_253_4,
            0.0,
            3.077_683_537_175_253_4,
            31.820_515_953_773_96,
            318.308_838_985_550_45,
        ],
        [
            -22.327_124_770_119_875,
            -6.964_556_734_283_274,
            -1.885_618_083_164_126_7,
            0.0,
            1.885_618_083_164_126_7,
            6.964_556_734_283_274,
            22.327_124_770_119_875,
        ],
        [
            -10.214_531_852_407_387,
            -4.540_702_858_568_133_6,
            -1.637_744_353_696_210_1,
            0.0,
            1.637_744_353_696_210_1,
            4.540_702_858_568_133_6,
            10.214_531_852_407_387,
        ],
        [
            -7.173_182_219_782_308_5,
            -3.746_947_387_979_197,
            -1.533_206_274_058_943_9,
            0.0,
            1.533_206_274_058_943_9,
            3.746_947_387_979_197,
            7.173_182_219_782_308_5,
        ],
        [
            -4.143_700_494_046_59,
            -2.763_769_458_112_696,
            -1.372_183_641_110_335_6,
            0.0,
            1.372_183_641_110_335_6,
            2.763_769_458_112_696,
            4.143_700_494_046_59,
        ],
    ];
    for (i, &nu) in NU.iter().enumerate() {
        let x = P.map(|p| xsf::stdtri(nu, p));
        xsref::np_assert_allclose!(x, X_REF[i], rtol = 1e-13);
    }
}

#[test]
fn test_stdtri_huge_t() {
    // reference values from mpmath
    for &(nu, p, t) in &[
        (0.0017, 0.79, 8.561_742_739_129_646e219),
        (0.0017, 0.21, -8.561_742_739_128_98e219),
        (1.08, 5.3e-315, -3.489_628_880_599_662e290),
        (2.02, 1e-320, -1.860_555_975_007_462_2e158),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtri(nu, p)], &[t], rtol = 2e-13);
    }
}

#[test]
fn test_stdtri_large_nu() {
    // reference values from mpmath
    for &(nu, p, t) in &[
        (10.0, 0.500_000_000_1, 2.569_978_247_571_428_5e-10),
        (3e19, 0.499_999_8, -5.013_256_549_406_37e-7),
        (1e20, 0.01, -2.326_347_874_040_841),
        (1e6, 0.01, -2.326_351_603_120_805_4),
        (1e10, 1e-100, -21.273_453_802_184_924),
        (1001.0, 1e-310, -55.953_590_999_699_58),
        (1e200, 0.01, -2.326_347_874_040_841),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtri(nu, p)], &[t], rtol = 1e-15);
    }
}

#[test]
fn test_ndtr_f64() {
    xsref::test("ndtr", "d-d", |x| xsf::ndtr(x[0]));
}

#[test]
fn test_ndtr_c64() {
    xsref::test("ndtr", "cd-cd", |x| xsf::ndtr(c64(x[0], x[1])));
}

// Based on `scipy.special.tests.test_ndtr.TestLogNdtr.test_log_ndtr_moderate_le8`
#[test]
#[allow(clippy::approx_constant)]
fn test_log_ndtr() {
    let x = [-0.75, -0.25, 0.0, 0.5, 1.5, 2.5, 3.0, 4.0, 5.0, 7.0, 8.0];
    let expected = [
        -1.484_448_229_919_656_2,
        -0.913_061_764_811_135_1,
        -0.693_147_180_559_945_3,
        -0.368_946_415_288_656_4,
        -0.069_143_455_612_233_98,
        -0.006_229_025_485_860_002,
        -0.001_350_809_964_748_193_8,
        -3.167_174_337_748_927e-05,
        -2.866_516_129_637_636e-07,
        -1.279_812_543_886_654e-12,
        -6.220_960_574_271_786e-16,
    ];
    let y = x.map(xsf::log_ndtr);
    xsref::np_assert_allclose!(&y, &expected, rtol = 1e-14);
}

#[test]
fn test_ndtri() {
    xsref::test("ndtri", "d-d", |x| xsf::ndtri(x[0]));
}

#[test]
fn test_von_mises_cdf() {
    // mpmath: the integral of the density from -pi to x, plus the period count
    for (k, x, expected) in [
        (1.0, 0.5, 0.664_076_474_604_993_1),
        (10.0, -1.0, 0.001_410_801_706_404_327_5),
        (0.1, 3.0, 0.979_653_423_965_240_7),
        (5.0, 7.0, 1.935_915_236_709_023),
        (49.0, 0.1, 0.757_381_643_056_756),
        (100.0, 0.2, 0.976_932_103_211_657),
        // left tail
        (20.0, -3.1, 3.150_549_424_635_017e-19),
        (25.0, -3.0, 5.908_127_636_104_764e-23),
        (100.0, -1.0, 5.101_678_232_611_33e-22),
        // `2k cos(x/2)^2` just above a Gauss-Laguerre node, and huge `k`
        (30.0, -1.201_347_072_474_133_7, 3.642_830_244_123_969e-10),
        (1e308, -5e-154, 2.866_515_718_791_936e-7),
        (f64::INFINITY, -1.0, 0.0),
        (f64::INFINITY, 1.0, 1.0),
    ] {
        xsref::np_assert_allclose!([xsf::von_mises_cdf(k, x)], [expected], rtol = 1e-14);
    }

    assert!(xsf::von_mises_cdf(-1.0, 0.5).is_nan());
    assert!(xsf::von_mises_cdf(f64::NAN, -2.0).is_nan());
    assert!(xsf::von_mises_cdf(f64::NAN, 2.0).is_nan());
    assert_eq!(xsf::von_mises_cdf(1.0, f64::INFINITY), f64::INFINITY);
}

#[test]
fn test_geninvgauss() {
    // (x, p, b), with the reference values from mpmath
    let args = [
        (0.01, 0.5, 1.0),
        (1.0, 0.5, 1.0),
        (5.0, 0.5, 1.0),
        (0.3, -2.5, 10.0),
        (7.0, 3.0, 0.25),
    ];
    let logpdf = [
        -47.621_353_440_210_626,
        -0.918_938_533_204_672_8,
        -3.323_657_489_421_723,
        -4.005_586_780_467_203,
        -3.924_725_429_203_324,
    ];
    let pdf = [
        2.081_176_820_202_832e-21,
        0.398_942_280_401_432_7,
        0.036_020_844_672_153_66,
        0.018_213_598_738_512_592,
        0.019_747_558_229_980_97,
    ];
    let actual = args.map(|(x, p, b)| xsf::geninvgauss_logpdf(x, p, b));
    xsref::np_assert_allclose!(&actual, &logpdf, rtol = 1e-14);
    let actual = args.map(|(x, p, b)| xsf::geninvgauss_pdf(x, p, b));
    xsref::np_assert_allclose!(&actual, &pdf, rtol = 1e-14);

    assert_eq!(xsf::geninvgauss_logpdf(0.0, 0.5, 0.5), f64::NEG_INFINITY);
    assert_eq!(xsf::geninvgauss_pdf(0.0, 0.5, 0.5), 0.0);
    assert_eq!(xsf::geninvgauss_pdf(2e6, 50.0, 2.0), 0.0);
    assert!(xsf::geninvgauss_logpdf(1.0, 0.5, 0.0).is_nan());
    assert!(xsf::geninvgauss_logpdf(1.0, 0.5, -1.0).is_nan());
}

#[test]
fn test_genhyperbolic() {
    // (x, p, a, b), with the reference values from mpmath
    let args = [
        (-3.0, 2.0, 3.0, 1.5),
        (0.0, 2.0, 3.0, 1.5),
        (1.2, 2.0, 3.0, 1.5),
        (4.0, -0.7, 0.5, -0.2),
        (-10.0, 1.0, 2.0, 0.0),
    ];
    let logpdf = [
        -11.470_885_563_926_888,
        -1.447_878_690_690_421_3,
        -0.982_289_498_087_753,
        -5.972_599_781_992_024,
        -18.825_827_120_241_21,
    ];
    let pdf = [
        1.042_936_092_158_763_5e-5,
        0.235_068_412_380_730_5,
        0.374_452_807_711_144_86,
        0.002_547_609_572_561_395,
        6.668_791_884_131_342e-9,
    ];
    let actual = args.map(|(x, p, a, b)| xsf::genhyperbolic_logpdf(x, p, a, b));
    xsref::np_assert_allclose!(&actual, &logpdf, rtol = 1e-14);
    let actual = args.map(|(x, p, a, b)| xsf::genhyperbolic_pdf(x, p, a, b));
    xsref::np_assert_allclose!(&actual, &pdf, rtol = 1e-14);

    // Cauchy limit: p = -1/2, a -> 0, b = 0
    let x = [-30.0, -1.0, 0.0, 0.5, 20.0];
    let actual = x.map(|x| xsf::genhyperbolic_pdf(x, -0.5, f64::from(f32::EPSILON), 0.0));
    let expected = x.map(|x| 1.0 / (core::f64::consts::PI * (1.0 + x * x)));
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-6);
}

#[test]
fn test_nrdtrimn_nrdtrisd() {
    // reference values from scipy.special
    xsref::np_assert_allclose!([xsf::nrdtrimn(0.3, 2.0, 1.0)], [2.048_801_025_416_082]);
    xsref::np_assert_allclose!([xsf::nrdtrisd(1.0, 0.7, 3.0)], [3.813_878_803_572_981]);
    assert!(xsf::nrdtrimn(0.5, 0.0, 1.0).is_nan());
    assert!(xsf::nrdtrisd(0.0, 1.0, 1.0).is_nan());
}

#[test]
fn test_kolmogorov() {
    xsref::test("kolmogorov", "d-d", |x| xsf::kolmogorov(x[0]));
}

#[test]
fn test_kolmogc() {
    xsref::test("kolmogc", "d-d", |x| xsf::kolmogc(x[0]));
}

#[test]
fn test_kolmogi() {
    xsref::test("kolmogi", "d-d", |x| xsf::kolmogi(x[0]));
}

#[test]
fn test_kolmogci() {
    xsref::test("kolmogci", "d-d", |x| xsf::kolmogci(x[0]));
}

#[test]
fn test_kolmogp() {
    xsref::test("kolmogp", "d-d", |x| xsf::kolmogp(x[0]));
}

#[test]
fn test_smirnov() {
    xsref::test("smirnov", "p_d-d", |x| xsf::smirnov(x[0] as i32, x[1]));
}

#[test]
fn test_smirnovc() {
    xsref::test("smirnovc", "p_d-d", |x| xsf::smirnovc(x[0] as i32, x[1]));
}

#[test]
fn test_smirnovi() {
    xsref::test("smirnovi", "p_d-d", |x| xsf::smirnovi(x[0] as i32, x[1]));
}

#[test]
fn test_smirnovci() {
    xsref::test("smirnovci", "p_d-d", |x| xsf::smirnovci(x[0] as i32, x[1]));
}

#[test]
fn test_smirnovp_large_n() {
    // in xsf this is -0.0, due to `int` overflow
    let (n, x, h) = (50_000, 2e-3, 1e-8);
    let expected = (xsf::smirnov(n, x + h) - xsf::smirnov(n, x - h)) / (2.0 * h);
    let actual = xsf::smirnovp(n, x);
    assert!(
        ((actual - expected) / expected).abs() < 1e-6,
        "{actual} != {expected}"
    );
}

#[test]
fn test_smirnovp() {
    xsref::test("smirnovp", "p_d-d", |x| xsf::smirnovp(x[0] as i32, x[1]));
}

#[test]
fn test_tukeylambdacdf_small_p_or_lambda() {
    // mpmath; the first is ill-conditioned (condition number ~2600)
    for (x, lambda, expected, tol) in [
        (
            -4.687_066_523_927_919_4,
            0.212_969_479_099_442_08,
            1.287_389_454_410_983_2e-13,
            1e-12,
        ),
        (-20.0, 9e-5, 2.024_340_980_204_265e-9, 1e-14),
        (-20.0, -9e-5, 2.098_545_034_667_178_4e-9, 1e-14),
        (3.0, 5e-5, 0.952_584_620_711_382_7, 1e-15),
    ] {
        let rel = (xsf::tukeylambdacdf(x, lambda) - expected).abs() / expected;
        assert!(rel < tol, "tukeylambdacdf({x}, {lambda}): rel = {rel:e}");
    }
}

#[test]
fn test_owens_t() {
    xsref::test("owens_t", "d_d-d", |x| xsf::owens_t(x[0], x[1]));
}

#[test]
fn test_chdtr() {
    xsref::test("chdtr", "d_d-d", |x| xsf::chdtr(x[0], x[1]));
}

#[test]
fn test_chdtrc() {
    xsref::test("chdtrc", "d_d-d", |x| xsf::chdtrc(x[0], x[1]));
}

#[test]
fn test_chdtri() {
    xsref::test("chdtri", "d_d-d", |x| xsf::chdtri(x[0], x[1]));
}

#[test]
fn test_fdtr() {
    xsref::test("fdtr", "d_d_d-d", |x| xsf::fdtr(x[0], x[1], x[2]));
}

#[test]
fn test_fdtrc() {
    xsref::test("fdtrc", "d_d_d-d", |x| xsf::fdtrc(x[0], x[1], x[2]));
}

#[test]
fn test_fdtr_fdtrc_w_near_1() {
    // reference values from mpmath
    let p = xsf::fdtr(1.38, 0.112, 5.77e62);
    xsref::np_assert_allclose!(&[p], &[0.999_743_687_449_994_5], rtol = 1e-13);
    let q = xsf::fdtrc(0.00118, 1.46e11, 0.00286);
    xsref::np_assert_allclose!(&[q], &[0.007_474_128_501_559_321], rtol = 1e-12);
    // `incbet` slightly exceeds 1 here (the true value is 4.6245e-10)
    let p = xsf::fdtr(
        458_884.656_463_958_03,
        2.951_258_965_120_171e-11,
        1_072.009_775_175_712_8,
    );
    assert!((0.0..1e-9).contains(&p));
    let q = xsf::fdtrc(
        2.951_258_965_120_171e-11,
        458_884.656_463_958_03,
        0.000_932_827_314_784_597_3,
    );
    assert!((0.0..1e-9).contains(&q));
}

#[test]
fn test_fdtr_fdtrc_tiny_df() {
    // reference values from mpmath; one df is tiny and the other large, so that `w` is near 1
    // but below the mean
    for &(a, b, x, p) in &[
        (1e7, 0.001, 0.5, 0.003_160_968_634_795_953_3),
        (1e7, 0.002, 1.0, 0.006_312_353_390_409_258),
    ] {
        xsref::np_assert_allclose!(&[xsf::fdtr(a, b, x)], &[p], rtol = 1e-12);
    }
    // inaccurate either way, but the complement would be far worse
    let p = xsf::fdtr(1e16, 1e-13, 6e-14);
    xsref::np_assert_allclose!(&[p], &[1.462_771_840_746_950_6e-14], rtol = 0.02);
    for &(a, b, x, q) in &[
        (0.002, 1e7, 1.0, 0.006_312_353_390_409_258),
        (
            0.001_153_561_542_308_523_6,
            731_611.160_931_477_6,
            3.267_656_625_181_253,
            0.003_281_744_887_778_894_7,
        ),
    ] {
        xsref::np_assert_allclose!(&[xsf::fdtrc(a, b, x)], &[q], rtol = 1e-12);
    }
}

#[test]
fn test_fdtr_bdtr_large_df() {
    // reference values from mpmath (see `incbet_lnpow`)
    let p = xsf::fdtr(100.0, 1e20, 0.5);
    xsref::np_assert_allclose!(&[p], &[6.953_305_247_616_099e-6], rtol = 1e-12);
    let p = xsf::bdtr(3e8, 1_000_000_000, 0.3);
    xsref::np_assert_allclose!(&[p], &[0.500_015_600_125_545], rtol = 1e-11);
}

#[test]
fn test_bdtr_nbdtrc_large_n() {
    // reference values from mpmath
    let p = [
        xsf::bdtr(0.0, 1_000_000_000, 1e-10),
        xsf::bdtr(2.0, 2_000_000_000, 1e-9),
        xsf::nbdtrc(2_000_000_000, 1, 1e-13),
        xsf::nbdtrc(2_000_000_000, 3, 1e-9),
    ];
    let p_ref = [
        0.904_837_418_031_435_4,
        0.676_676_416_183_063_4,
        0.999_800_019_998_566_7,
        0.676_676_415_371_051_7,
    ];
    xsref::np_assert_allclose!(&p, &p_ref, rtol = 1e-14);
}

#[test]
fn test_fdtr_stdtr_huge_df() {
    // reference values from mpmath (see `betainc_xy`)
    let abx = [
        (100.0, 1e20, 0.5),
        (30.0, 1e20, 2.0),
        (1e20, 30.0, 2.0),
        (3.0, 1e14, 2.0),
    ];
    let p = [
        6.953_305_247_616_099e-6,
        0.999_079_317_603_851_4,
        0.989_739_572_087_657_4,
        0.888_389_774_905_280_1,
    ];
    let q = [
        0.999_993_046_694_752_4,
        9.206_823_961_486_663e-4,
        0.010_260_427_912_342_617,
        0.111_610_225_094_719_86,
    ];
    let rtol = 1e-13;
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtr(a, b, x)), &p, rtol = rtol);
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtrc(a, b, x)), &q, rtol = rtol);

    let p = [xsf::stdtr(1e20, -3.0), xsf::stdtr(1e20, 3.0)];
    let q = 0.001_349_898_031_630_094_6;
    xsref::np_assert_allclose!(&p, &[q, 1.0 - q], rtol = 1e-14);
}

#[test]
fn test_fdtri() {
    xsref::test("fdtri", "d_d_d-d", |x| xsf::fdtri(x[0], x[1], x[2]));
}

#[test]
fn test_gdtr() {
    xsref::test("gdtr", "d_d_d-d", |x| xsf::gdtr(x[0], x[1], x[2]));
}

#[test]
fn test_gdtrc() {
    xsref::test("gdtrc", "d_d_d-d", |x| xsf::gdtrc(x[0], x[1], x[2]));
}

#[test]
fn test_pdtr() {
    xsref::test("pdtr", "d_d-d", |x| xsf::pdtr(x[0], x[1]));
}

#[test]
fn test_pdtrc() {
    xsref::test("pdtrc", "d_d-d", |x| xsf::pdtrc(x[0], x[1]));
}

#[test]
fn test_pdtri() {
    xsref::test("pdtri", "p_d-d", |x| xsf::pdtri(x[0] as i32, x[1]));
}

#[test]
fn test_bdtr() {
    xsref::test("bdtr", "d_p_d-d", |x| xsf::bdtr(x[0], x[1] as i32, x[2]));
}

#[test]
fn test_bdtri_small_y() {
    // mpmath: the p for which bdtr(96, 191, p) == 2.1012712563041036e-22
    let p = xsf::bdtri(96.0, 191, 2.101_271_256_304_103_6e-22);
    assert!((p / 0.814_405_035_443_584 - 1.0).abs() < 1e-14, "p = {p}");
}

#[test]
fn test_bdtrc() {
    xsref::test("bdtrc", "d_p_d-d", |x| xsf::bdtrc(x[0], x[1] as i32, x[2]));
}

#[test]
fn test_bdtri() {
    xsref::test("bdtri", "d_p_d-d", |x| xsf::bdtri(x[0], x[1] as i32, x[2]));
}

#[test]
fn test_poisson_binom_smoke() {
    // Subset of xsf's `test_poisson_binom.cpp`
    const K: [usize; 8] = [0, 7, 14, 21, 25, 26, 44, 51];
    const PMF_REF: [f64; 8] = [
        3.428_603_61e-26,
        2.575_893_01e-11,
        4.150_992_19e-5,
        4.122_831_29e-2,
        1.357_810_59e-1,
        1.357_810_59e-1,
        2.575_893_01e-11,
        3.428_603_61e-26,
    ];

    // np.linspace(1e-5, 1 - 1e-5, 51)
    let p: Vec<f64> = (0..51)
        .map(|i| 1e-5 + f64::from(i) * (1.0 - 2e-5) / 50.0)
        .collect();

    let pmf = xsf::poisson_binom_pmf_all(&p);
    let cdf = xsf::poisson_binom_cdf_all(&p);

    xsref::np_assert_allclose!(K.map(|k| pmf[k]), PMF_REF, rtol = 1e-8);

    xsref::np_assert_allclose!([cdf[25], cdf[51]], [0.5, 1.0], atol = 1e-15);
}

#[test]
fn test_wilcoxon() {
    // Subset of xsf's `test_wilcoxon.cpp`
    assert_eq!(xsf::wilcoxon_pmf_all(0), [1.0]);
    assert_eq!(xsf::wilcoxon_pmf_all(2), [0.25; 4]);
    assert_eq!(
        xsf::wilcoxon_pmf_all(5),
        [
            0.03125, 0.03125, 0.03125, 0.0625, 0.0625, 0.09375, 0.09375, 0.09375, 0.09375, 0.09375,
            0.09375, 0.0625, 0.0625, 0.03125, 0.03125, 0.03125,
        ]
    );
    assert_eq!(
        xsf::wilcoxon_cdf_all(3),
        [0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 1.0]
    );

    let pmf = xsf::wilcoxon_pmf_all(3);
    assert_eq!(xsf::wilcoxon_pmf(3.0, &pmf), 0.25);
    assert_eq!(xsf::wilcoxon_pmf(-1.0, &pmf), 0.0);
    assert_eq!(xsf::wilcoxon_pmf(7.0, &pmf), 0.0);
    assert!(xsf::wilcoxon_pmf(f64::NAN, &pmf).is_nan());

    let cdf = xsf::wilcoxon_cdf_all(5);
    let k = [-1.0, -0.5, 5.0, 7.9, 9.0, 15.0, f64::INFINITY];
    assert_eq!(
        k.map(|k| xsf::wilcoxon_cdf(k, &cdf)),
        [0.0, 0.03125, 0.3125, 0.5, 0.6875, 1.0, 1.0]
    );
    assert_eq!(
        k.map(|k| xsf::wilcoxon_sf(k, &cdf)),
        [1.0, 1.0, 0.78125, 0.59375, 0.40625, 0.03125, 0.0]
    );
    assert_eq!(xsf::wilcoxon_sf(f64::NEG_INFINITY, &cdf), 1.0);
    assert!(xsf::wilcoxon_cdf(f64::NAN, &cdf).is_nan());
    assert!(xsf::wilcoxon_sf(f64::NAN, &cdf).is_nan());
}

#[test]
fn test_nbdtr() {
    xsref::test("nbdtr", "p_p_d-d", |x| {
        xsf::nbdtr(x[0] as i32, x[1] as i32, x[2])
    });
}

#[test]
fn test_nbdtrc() {
    xsref::test("nbdtrc", "p_p_d-d", |x| {
        xsf::nbdtrc(x[0] as i32, x[1] as i32, x[2])
    });
}

#[test]
fn test_cdf_cvm_inf_smoke() {
    assert!(xsf::cdf_cvm_inf(f64::NAN).is_nan());
    assert!(xsf::cdf_cvm_inf(-1.0).abs() <= f64::EPSILON);

    let result = xsf::cdf_cvm_inf(0.5);
    assert!(result.is_finite());
    assert!(result > 0.0 && result < 1.0);
}

#[test]
fn test_cdf_cvm_smoke() {
    const N: i32 = 10;
    let lower = 1.0 / (12.0 * f64::from(N));
    let upper = f64::from(N) / 3.0;

    assert!(xsf::cdf_cvm(lower, N).abs() <= f64::EPSILON);
    assert!((xsf::cdf_cvm(upper, N) - 1.0).abs() <= f64::EPSILON);

    let result = xsf::cdf_cvm(0.5, N);
    assert!(result.is_finite());
    assert!((0.0..=1.0).contains(&result));
}

#[test]
fn test_cdf_cvm_nonfinite() {
    assert_eq!(xsf::cdf_cvm_inf(f64::INFINITY), 1.0);
    assert_eq!(xsf::cdf_cvm_inf(f64::MAX), 1.0);
    assert!(xsf::cdf_cvm_inf(f64::NAN).is_nan());
    for n in [-1, 5] {
        assert_eq!(xsf::cdf_cvm(f64::INFINITY, n), 1.0);
        assert!(xsf::cdf_cvm(f64::NAN, n).is_nan());
    }
    assert!(xsf::cdf_cvm(1e4, 1_000_000).is_nan());
}
