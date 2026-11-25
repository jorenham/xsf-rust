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
        // nu = 2, F(x) = 1/2 + x / (2 sqrt(2 + x^2))
        (2.0, -1e5, 4.999_999_999_25e-11),
        (2.0, -1e150, 5e-301),
        (2.0, -1e154, 5e-309),
        (2.0, 1e200, 1.0),
        (2.0, f64::NEG_INFINITY, 0.0),
        (2.0, f64::INFINITY, 1.0),
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
fn test_stdtr_invalid_degree() {
    for x in [f64::NEG_INFINITY, -1.0, -0.0, 0.0, 1.0, f64::INFINITY, f64::NAN] {
        assert!(xsf::stdtr(f64::NAN, x).is_nan());
    }
}

#[test]
fn test_stdtr_cauchy_tail() {
    // Cauchy CDF at exact binary64 inputs, evaluated in high precision.
    for (x, reference) in [
        (-1.000_000_000_000_000_2, 0.249_999_999_999_999_97),
        (-1.0, 0.25),
        (-0.999_999_999_999_999_9, 0.25),
        (-1e8, 3.183_098_861_837_906_6e-9),
        (-1e20, 3.183_098_861_837_906_7e-21),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtr(1.0, x)], &[reference], rtol = 6e-16);
    }
    for (x, reference) in
        [(-1e308, 3.183_098_861_837_906e-309_f64), (-f64::MAX, 1.770_657_516_629_887e-309)]
    {
        assert!(xsf::stdtr(1.0, x).to_bits().abs_diff(reference.to_bits()) <= 1);
        assert_eq!(xsf::stdtr(1.0, -x), 1.0);
    }
    assert_eq!(xsf::stdtr(1.0, f64::NEG_INFINITY).to_bits(), 0);
    assert_eq!(xsf::stdtr(1.0, -0.0), 0.5);
    assert_eq!(xsf::stdtr(1.0, 0.0), 0.5);
    assert_eq!(xsf::stdtr(1.0, f64::INFINITY), 1.0);
    assert!(xsf::stdtr(1.0, f64::NAN).is_nan());
}

#[test]
fn test_stdtr_degree_two_subnormal_tail() {
    // Arb evaluation of 1 / (sqrt(2 + x²) * (sqrt(2 + x²) - x)).
    assert_eq!(xsf::stdtr(2.0, -4.740_375_954_054_59e153).to_bits(), 0x000f_ffff_ffff_fffd);

    let x = -2.0_f64.powi(537);
    let above = f64::from_bits(x.to_bits() - 1);
    assert_eq!(xsf::stdtr(2.0, above).to_bits(), 1);
    assert_eq!(xsf::stdtr(2.0, x).to_bits(), 0);
}

#[test]
fn test_stdtr_tail() {
    // Arb
    for &(nu, x, p) in &[
        (30_003.0, -28.375, 4.205_739_883_671_263_5e-175),
        (5e4, -37.0, 5.762_069_766_213_904e-296),
        (1e5, -30.0, 3.689_268_436_111_116_6e-197),
        (2e5, -2.5, 0.006_210_062_456_181_062),
        (1e6, -20.0, 2.866_543_523_695_186e-89),
        (3e7, -10.0, 7.620_500_678_968_407e-24),
        (5e8, -35.0, 1.125_756_434_871_158_7e-268),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtr(nu, x)], &[p], rtol = 2e-15);
    }
}

#[test]
fn test_stdtr_huge_nu() {
    // reference values from the 1/nu series up to nu^-6 (mpmath), which agrees with Arb to 2e-27
    // at nu = 1e9
    for &(nu, x, p) in &[
        (2e9, -37.0, 5.726_914_668_651_248e-300),
        (3e10, -30.0, 4.906_747_121_138_99e-198),
        (5e12, -20.0, 2.753_624_140_745_236_8e-89),
        (7e13, 2.5, 0.993_790_334_674_222_7),
        (2e14, -37.0, 5.725_571_235_957_468e-300),
        (1e16, -35.0, 1.124_910_706_514_676_8e-268),
        (1e18, -30.0, 4.906_713_927_149_183e-198),
        (1e20, -8.0, 6.220_960_574_271_784e-16),
        (1e21, 1.5, 0.933_192_798_731_141_9),
        (1e30, -10.0, 7.619_853_024_160_525e-24),
        (1e50, -36.85, 1.462_274_937_250_674_6e-297),
        (1e100, 0.9, 0.815_939_874_653_240_5),
        (1e120, -3.0, 0.001_349_898_031_630_094_6),
        (1e150, -30.0, 4.906_713_927_148_187e-198),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtr(nu, x)], &[p], rtol = 1e-15);
    }
}

#[test]
fn test_stdtr_huge_nu_half_minsub() {
    // Adjacent x values on either side of 2^-1075, from the normalized density integral.
    for (nu, below, above) in [
        (1_000_000_000.000_000_1, -38.485_422_595_634_134, -38.485_422_595_634_13),
        (1e14, -38.485_408_335_709_95, -38.485_408_335_709_94),
        (9.894_668_870_662_003e17, -38.485_408_335_567_36, -38.485_408_335_567_35),
    ] {
        assert_eq!(xsf::stdtr(nu, below).to_bits(), 0);
        assert_eq!(xsf::stdtr(nu, above).to_bits(), 1);
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
fn test_stdtri_tiny_p() {
    // reference values from mpmath
    for &(nu, p, t) in &[
        (1.0, 1e-20, -3.183_098_861_837_907e19),
        (2.0, 5e-324, -3.181_212_452_095_196_2e161),
        (100.0, 5e-324, -16_559.999_468_166_71),
        (1e4, 5e-324, -39.936_370_979_369_5),
        (100_001.0, 5e-324, -38.610_245_498_234_78),
        (1e200, 5e-324, -38.467_405_617_144_346),
    ] {
        xsref::np_assert_allclose!(&[xsf::stdtri(nu, p)], &[t], rtol = 2e-15);
    }
}

#[test]
fn test_stdtri_tiny_nu() {
    // |t| > f64::MAX for every p != 1/2 if nu < 1.5e-19 (Arb)
    for nu in [5e-324, 1e-310, 1e-300, 1e-19] {
        assert_eq!(xsf::stdtri(nu, 1e-300), f64::NEG_INFINITY);
        assert_eq!(xsf::stdtri(nu, 0.5 - f64::EPSILON / 4.0), f64::NEG_INFINITY);
        assert_eq!(xsf::stdtri(nu, 0.5), 0.0);
        assert_eq!(xsf::stdtri(nu, 0.5 + f64::EPSILON / 2.0), f64::INFINITY);
        assert_eq!(xsf::stdtri(nu, 0.7), f64::INFINITY);
    }
}

xsref::tests! {
    test_ndtr_f64: "ndtr", "d-d", |x| xsf::ndtr(x[0]);
    test_ndtr_c64: "ndtr", "cd-cd", |x| xsf::ndtr(c64(x[0], x[1]));
}

#[test]
fn test_ndtr_complex_small_phase_overflow() {
    let min_subnormal = f64::from_bits(1);
    let phase_limit = 0.5 / 64.0;
    for (x, y) in [
        (1e-30, 316.227_766_016_837_96),
        (min_subnormal, f64::MAX),
        (min_subnormal, 64.0),
        (f64::from_bits(f64::to_bits(phase_limit) - 1), 64.0),
    ] {
        for x in [-x, x] {
            for y in [-y, y] {
                let expected = c64(f64::INFINITY.copysign(x), f64::INFINITY.copysign(y));
                assert_eq!(xsf::ndtr(c64(x, y)), expected);
            }
        }
    }
    for x in [-0.0, 0.0] {
        assert_eq!(xsf::ndtr(c64(x, 64.0)), c64(0.5, f64::INFINITY));
    }
}

#[test]
fn test_ndtr_complex_underflow_signs() {
    for (y, negative_re) in [(-0.0625, true), (-f64::from_bits(1), false)] {
        let z = c64(-40.0, y);
        let p = xsf::ndtr(z);
        assert_eq!(p, c64(0.0, 0.0));
        assert_eq!(p.re.is_sign_negative(), negative_re);
        assert!(p.im.is_sign_negative());

        let conjugate = xsf::ndtr(z.conj());
        assert_eq!(conjugate, c64(0.0, 0.0));
        assert_eq!(conjugate.re.to_bits(), p.re.to_bits());
        assert_eq!(conjugate.im.to_bits(), 0);

        let reflected = xsf::ndtr(-z);
        assert_eq!(reflected.re.to_bits(), 1.0_f64.to_bits());
        assert_eq!(reflected.im.to_bits(), 0);
        assert_eq!(xsf::ndtr(-z.conj()).im.to_bits(), p.im.to_bits());
    }
}

#[test]
fn test_ndtr_tails() {
    // Reference values from mpmath (`ncdf`, `owens_t` in tools/accuracy/specs.py, 50 digits),
    // confirmed with Wolfram (`CDF[NormalDistribution[], x]`, NIntegrate of Owen's T)
    let y = [-20.0, -30.0].map(xsf::ndtr);
    xsref::np_assert_allclose!(
        &y,
        &[2.753_624_118_606_233_7e-89, 4.906_713_927_148_187e-198],
        rtol = 1e-15
    );
    // subnormal, to a multiple of 2^-1074
    let y = [-38.0, -38.4].map(xsf::ndtr);
    xsref::np_assert_allclose!(&y, &[2.885_428_35e-316, 6.4e-323], atol = 5e-324);
    // ~0.50016 and ~0.419 times 2^-1074
    assert_eq!(xsf::ndtr(-38.4854), f64::from_bits(1));
    assert_eq!(xsf::ndtr(-38.49), 0.0);
    let y = xsf::log_ndtr(30.0);
    xsref::np_assert_allclose!(&[y], &[-4.906_713_927_148_187e-198], rtol = 1e-15);
    let t = xsf::owens_t(30.0, 2.0);
    xsref::np_assert_allclose!(&[t], &[2.453_356_963_574_093_6e-198], rtol = 1e-15);
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
fn test_ndtr_c64_bounded_overflow_phase() {
    for (z, expected) in [
        (c64(10.0, 316.227_766_016_837_96), c64(f64::INFINITY, f64::NEG_INFINITY)),
        (c64(-316.227_766_016_837_96, 10_000.0), c64(f64::NEG_INFINITY, f64::INFINITY)),
        // The initial quotient gives an adjacent quadrant; signs from the exact-product phase.
        (c64(5_000_003.0, 1_000_000_002.0), c64(f64::INFINITY, f64::NEG_INFINITY)),
    ] {
        assert_eq!(xsf::ndtr(z), expected);
        assert_eq!(xsf::ndtr(z.conj()), expected.conj());
        assert_eq!(xsf::ndtr(-z), -expected);
        assert_eq!(xsf::ndtr(-z.conj()), -expected.conj());
    }
}

#[test]
fn test_ndtr_c64_overflow_like_erfc() {
    // Arb at the exact binary64 inputs: a wide phase, and both components just above MAX.
    for (z, expected) in [
        (c64(-1e5, -1e12), c64(f64::INFINITY, f64::INFINITY)),
        (c64(-1e5, 1e50), c64(f64::INFINITY, f64::INFINITY)),
        (c64(-100.0, -106.962_610_289_764_34), c64(f64::NEG_INFINITY, f64::NEG_INFINITY)),
    ] {
        assert_eq!(xsf::ndtr(z), expected);
        assert_eq!(xsf::ndtr(z.conj()), expected.conj());
        assert_eq!(xsf::ndtr(-z), -expected);
        assert_eq!(xsf::ndtr(-z.conj()), -expected.conj());
    }

    // An overflowing squared diagonal, and a finite part next to an overflowing one.
    for (z, re, im) in [
        (c64(-1e300, -1e300), 2.703_582_075_983_466_4e-301, -8.052_273_672_793_895e-302),
        (
            c64(-250.0, -252.835_915_170_293_8),
            3.027_401_123_304_904e306,
            -4.022_582_231_100_544e306,
        ),
    ] {
        let actual = xsf::ndtr(z);
        xsref::np_assert_allclose!([actual.re, actual.im], [re, im], rtol = 2e-15);
    }
    let actual = xsf::ndtr(c64(-250.0, -252.863_599_594_722_22));
    // The real part is below 1% of the modulus.
    xsref::np_assert_allclose!([actual.re], [3.968_472_849_220_208_6e307], rtol = 5e-14);
    assert_eq!(actual.im, f64::NEG_INFINITY);
}

#[test]
fn test_log_ndtr_c64_near_0() {
    // mpmath `log1p(-erfc(z / sqrt(2)) / 2)`, where `ndtr(-z)` is small
    let z = [
        c64(5.901_265_878_562_874, 0.416_992_545_962_521_9),
        c64(14.9375, -14.0),
        c64(4.5, -1.25),
        c64(6.25, -0.75),
    ];
    let expected = [
        c64(1.604_864_238_438_167_5e-9, 1.130_791_271_554_178e-9),
        c64(2.051_146_923_142_219_4e-8, -1.436_868_598_019_861_2e-8),
        c64(-6.615_497_872_514_756e-6, 2.853_073_261_677_222_5e-6),
        c64(-2.406_718_333_043_747_7e-11, 2.691_130_913_270_477e-10),
    ];
    xsref::np_assert_allclose!(z.map(xsf::log_ndtr), expected, rtol = 5e-14);
}

xsref::tests! {
    test_ndtri: "ndtri", "d-d", |x| xsf::ndtri(x[0]);
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
    let args =
        [(0.01, 0.5, 1.0), (1.0, 0.5, 1.0), (5.0, 0.5, 1.0), (0.3, -2.5, 10.0), (7.0, 3.0, 0.25)];
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

xsref::tests! {
    test_kolmogorov: "kolmogorov", "d-d", |x| xsf::kolmogorov(x[0]);
    test_kolmogc: "kolmogc", "d-d", |x| xsf::kolmogc(x[0]);
    test_kolmogi: "kolmogi", "d-d", |x| xsf::kolmogi(x[0]);
    test_kolmogci: "kolmogci", "d-d", |x| xsf::kolmogci(x[0]);
    test_kolmogp: "kolmogp", "d-d", |x| xsf::kolmogp(x[0]);
}

#[test]
fn test_kolmogp_near_cutover() {
    // -8x sum_k (-1)^(k-1) k^2 exp(-2k^2 x^2), mpmath (dps = 40)
    let x = [0.83, 0.85, 0.9, 1.0];
    let expected = [
        -1.567_071_014_748_016_2,
        -1.519_197_893_798_649,
        -1.380_727_054_237_732_8,
        -1.071_948_558_356_941_7,
    ];
    xsref::np_assert_allclose!(&x.map(xsf::kolmogp), &expected, rtol = 1e-15);
    // 2 sum_k (-1)^(k-1) exp(-2k^2 x^2), mpmath (dps = 40)
    let sf = xsf::kolmogorov(0.83);
    xsref::np_assert_allclose!(&[sf], &[0.496_190_995_350_505_8], rtol = 1e-15);
}

xsref::tests! {
    test_smirnov: "smirnov", "p_d-d", |x| xsf::smirnov(x[0] as i32, x[1]);
    test_smirnovc: "smirnovc", "p_d-d", |x| xsf::smirnovc(x[0] as i32, x[1]);
    test_smirnovi: "smirnovi", "p_d-d", |x| xsf::smirnovi(x[0] as i32, x[1]);
    test_smirnovci: "smirnovci", "p_d-d", |x| xsf::smirnovci(x[0] as i32, x[1]);
}

#[test]
fn test_smirnovp_large_n() {
    // scipy/xsf#315: `int` overflow made this -0.0
    let (n, x, h) = (50_000, 2e-3, 1e-8);
    let expected = (xsf::smirnov(n, x + h) - xsf::smirnov(n, x - h)) / (2.0 * h);
    let actual = xsf::smirnovp(n, x);
    assert!(((actual - expected) / expected).abs() < 1e-6, "{actual} != {expected}");
}

#[test]
fn test_smirnovp_signed_zero() {
    // the derivative underflows; `f64::max(0.0, -0.0)` gave +0 in debug builds
    let p = xsf::smirnovp(248, 0.992_187_5);
    assert_eq!(p.to_bits(), (-0.0_f64).to_bits(), "{p:?}");
}

xsref::tests! {
    test_smirnovp: "smirnovp", "p_d-d", |x| xsf::smirnovp(x[0] as i32, x[1]);
}

#[test]
fn test_tukeylambdacdf_small_p_or_lambda() {
    // mpmath; the first is ill-conditioned (condition number ~2600)
    for (x, lambda, expected, tol) in [
        (-4.687_066_523_927_919_4, 0.212_969_479_099_442_08, 1.287_389_454_410_983_2e-13, 1e-12),
        (-20.0, 9e-5, 2.024_340_980_204_265e-9, 1e-14),
        (-20.0, -9e-5, 2.098_545_034_667_178_4e-9, 1e-14),
        (3.0, 5e-5, 0.952_584_620_711_382_7, 1e-15),
    ] {
        let rel = (xsf::tukeylambdacdf(x, lambda) - expected).abs() / expected.abs();
        assert!(rel < tol, "tukeylambdacdf({x}, {lambda}): rel = {rel:e}");
    }
}

xsref::tests! {
    test_owens_t: "owens_t", "d_d-d", |x| xsf::owens_t(x[0], x[1]);
}

#[test]
fn test_owens_t_large_h_small_a() {
    // scipy/xsf#326, with the reference values from mpmath
    for (h, a, expected) in [
        (5.0, 1e-10, 5.931_152_735_254_122e-17),
        (6.0, 1e-14, 2.423_926_559_560_457e-23),
        (8.0, 1e-3, 2.015_542_376_282_099_7e-18),
    ] {
        let t = xsf::owens_t(h, a);
        xsref::np_assert_allclose!(&[t], &[expected], rtol = 1e-13);
    }
}

xsref::tests! {
    test_chdtr: "chdtr", "d_d-d", |x| xsf::chdtr(x[0], x[1]);
    test_chdtrc: "chdtrc", "d_d-d", |x| xsf::chdtrc(x[0], x[1]);
    test_chdtri: "chdtri", "d_d-d", |x| xsf::chdtri(x[0], x[1]);
    test_fdtr: "fdtr", "d_d_d-d", |x| xsf::fdtr(x[0], x[1], x[2]);
    test_fdtrc: "fdtrc", "d_d_d-d", |x| xsf::fdtrc(x[0], x[1], x[2]);
}

#[test]
fn test_fdtr_fdtrc_inf() {
    for &(a, b) in &[(1.0, 1.0), (1e-300, 1e300), (1e300, 1e-300)] {
        assert_eq!(xsf::fdtr(a, b, f64::INFINITY), 1.0);
        assert_eq!(xsf::fdtrc(a, b, f64::INFINITY), 0.0);
    }
    assert!(xsf::fdtr(f64::NAN, 1.0, f64::INFINITY).is_nan());
    assert!(xsf::fdtr(1.0, f64::NAN, f64::INFINITY).is_nan());
    // huge `a` or subnormal `b`
    assert_eq!(xsf::fdtrc(1e10, 5e-324, f64::INFINITY), 0.0);
    assert_eq!(xsf::fdtr(5e-324, 5e-324, f64::INFINITY), 1.0);
    assert_eq!(xsf::fdtr(1e10, 1.0, -0.0), 0.0);
    assert_eq!(xsf::fdtr(f64::INFINITY, 1.0, 1e-310), 0.0);
    assert_eq!(xsf::fdtrc(f64::INFINITY, 1.0, 1e-310), 1.0);
}

#[test]
fn test_fdtr_fdtrc_w_near_1() {
    // reference values from mpmath
    let p = xsf::fdtr(1.38, 0.112, 5.77e62);
    xsref::np_assert_allclose!(&[p], &[0.999_743_687_449_994_5], rtol = 1e-13);
    let q = xsf::fdtrc(0.00118, 1.46e11, 0.00286);
    xsref::np_assert_allclose!(&[q], &[0.007_474_128_501_559_321], rtol = 1e-12);
    // `incbet` slightly exceeds 1 here (the true value is 4.6245e-10)
    let p = xsf::fdtr(458_884.656_463_958_03, 2.951_258_965_120_171e-11, 1_072.009_775_175_712_8);
    assert!((0.0..1e-9).contains(&p));
    let q =
        xsf::fdtrc(2.951_258_965_120_171e-11, 458_884.656_463_958_03, 0.000_932_827_314_784_597_3);
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
    // `incbet` with tiny b and x near 1, where the continued fractions are cut off (mpmath and
    // Wolfram)
    let p = xsf::fdtr(1000.0, 6.2e-10, 0.3);
    xsref::np_assert_allclose!(&[p], &[6.235_420_792_551_117e-9], rtol = 1e-14);
    // inaccurate either way, but the complement would be far worse
    let p = xsf::fdtr(1e16, 1e-13, 6e-14);
    xsref::np_assert_allclose!(&[p], &[1.462_771_840_746_950_6e-14], rtol = 0.02);
    for &(a, b, x, q) in &[
        (0.002, 1e7, 1.0, 0.006_312_353_390_409_258),
        (0.00115, 731_611.0, 3.25, 0.003_276_496_748_104_418_5),
    ] {
        xsref::np_assert_allclose!(&[xsf::fdtrc(a, b, x)], &[q], rtol = 1e-12);
    }
}

#[test]
fn test_fdtr_fdtrc_df_below_eps() {
    // reference values from mpmath and Arb; one df is below eps times the other (see `incbet_xy`)
    let p = [
        xsf::fdtrc(1e-300, 1.0, 1.0),
        xsf::fdtr(1.0, 1e-300, 1.0),
        xsf::fdtr(0.1, 1e-20, 1e10),
        xsf::fdtr(3.0, 1e-20, 1e-3),
        // both df are below eps / 1000
        xsf::fdtr(1e-200, 1e-210, 1e100),
        // the small result is not computed as `1 - (1 - I)`
        xsf::fdtrc(2.3e-16, 1e6, 172.0),
        xsf::fdtrc(1e-15, 100.0, 0.093_75),
    ];
    let p_ref = [
        3.460_809_111_296_668e-298,
        3.460_809_111_296_668e-298,
        4.334_779_851_161_283e-19,
        1.981_442_661_534_338_8e-19,
        9.999_999_999_000_001e-11,
        3.562_342_396_438_712e-15,
        1.851_593_242_835_031_3e-14,
    ];
    xsref::np_assert_allclose!(&p, &p_ref, rtol = 1e-14);
}

#[test]
fn test_fdtr_fdtrc_ax_overflow() {
    // reference values from Arb; `a x` overflows and `1 - w = b / (b + a x)` is tiny
    let abx = [(4e178, 0.006, 8e182), (1e308, 10.0, 2.0), (1e200, 1e-250, 1e200)];
    let p = [0.721_727_090_722_264_3, 0.891_178_018_914_151_3, 5.181_396_116_814_895_5e-248];
    let q = [0.278_272_909_277_735_74, 0.108_821_981_085_848_76, 1.0];
    let rtol = 1e-13;
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtr(a, b, x)), &p, rtol = rtol);
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtrc(a, b, x)), &q, rtol = rtol);
    let q = xsf::fdtrc(1.4e13, 1.0, 1.4e306);
    xsref::np_assert_allclose!(&[q], &[6.743_355_313_447_234e-154], rtol = rtol);
    // fdtr(inf, b, x) = Q(b/2, b / (2 x))
    let p = xsf::fdtr(f64::INFINITY, 1.0, 1.0);
    xsref::np_assert_allclose!(&[p], &[0.317_310_507_862_914_1], rtol = rtol);
}

#[test]
fn test_fdtr_fdtrc_subnormal_1mw() {
    // reference values from Arb; `1 - w = b / (b + a x)` is subnormal or underflows
    let abx = [(2.0, 1e-278, 1e307), (0.1, 1e-300, 1e300), (1e8, 0.002, 1e308), (1.0, 1.0, 1e308)];
    let p = [6.738_527_132_910_383e-276, 6.995_845_500_149_159e-298, 0.511_065_903_382_357_5, 1.0];
    let q = [1.0, 1.0, 0.488_934_096_617_642_5, 6.366_197_723_675_813_5e-155];
    let rtol = 1e-14;
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtr(a, b, x)), &p, rtol = rtol);
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtrc(a, b, x)), &q, rtol = rtol);
    // `b / a` underflows, but `1 - w = 1/3`
    assert_eq!(xsf::fdtr(4.0, 1e-323, 5e-324), 0.0);
}

#[test]
fn test_fdtrc_subnormal_w() {
    // reference values from Arb; `w = a x / (b + a x)` is subnormal or underflows, or `b = inf`,
    // where `fdtrc` is the chi-square survival function of `a x`
    let abx = [
        (3.0, f64::INFINITY, 0.5),
        (1e-20, f64::INFINITY, 1.0),
        (0.01, 1e300, 1e-30),
        (0.5, 1e308, 1e-8),
        (1e-10, 1e300, 1e-20),
        // `a x` is subnormal or underflows
        (0.001, f64::INFINITY, 5e-324),
        (2.7e-185, f64::INFINITY, 3.96e-318),
        (1e-20, f64::INFINITY, 1e-310),
    ];
    let q = [
        0.682_270_330_336_212_5,
        2.308_381_668_776_966_3e-19,
        0.308_584_120_241_167_6,
        0.992_198_754_982_112_5,
        3.459_674_209_291_373e-9,
        0.313_213_691_198_101_83,
        1.560_528_095_486_543e-182,
        3.799_845_061_018_467e-18,
    ];
    xsref::np_assert_allclose!(&abx.map(|(a, b, x)| xsf::fdtrc(a, b, x)), &q, rtol = 1e-14);
    let p = xsf::fdtr(0.001, f64::INFINITY, 5e-324);
    xsref::np_assert_allclose!(&[p], &[0.686_786_308_801_898_2], rtol = 1e-14);
    // `a / 2` rounds; Arb gives 1.1026785e-320, i.e. 2232 times the smallest subnormal
    assert_eq!(xsf::fdtrc(1.5e-323, f64::INFINITY, 5e-324), 1.103e-320);
    // `a x` is just below `2 MIN_POSITIVE`; Arb gives 8.7427288e-321
    assert_eq!(xsf::fdtrc(2.5e-323, f64::INFINITY, 1_801_439_850_948_198.2), 8.745e-321);
    assert_eq!(xsf::fdtr(5e-324, f64::INFINITY, 1.0), 1.0);
    assert_eq!(xsf::fdtr(5e-324, f64::INFINITY, 0.0), 0.0);
    assert_eq!(xsf::fdtrc(5e-324, f64::INFINITY, 0.0), 1.0);
    // `a / 2` rounds, while `a x` is normal, or `fdtrc` is the gamma limit `Q(a / 2, t)` as `b` is
    // large; Arb gives 1.7482875e-321 and 8.7427288e-321
    assert_eq!(xsf::fdtrc(5e-324, f64::INFINITY, 1e16), 1.75e-321);
    assert_eq!(xsf::fdtrc(2.5e-323, 1e20, 1_801_439_850_948_198.2), 8.745e-321);
    for (x, q) in [(0.0, 1.0), (1e-310, 1.0), (f64::INFINITY, 0.0)] {
        assert_eq!(xsf::fdtrc(3.0, f64::INFINITY, x), q);
    }
}

#[test]
fn test_fdtr_tiny_w() {
    // `w` underflows or `a x` is subnormal, with a finite `b`, where `fdtrc = O(a)` for a tiny (or
    // subnormal) `a`, or `fdtr = O(b)` for a subnormal `b`; from Arb
    let p = [
        xsf::fdtrc(1e-307, 1.0, 5e-324),
        xsf::fdtrc(1e-308, 1.0, 1e-320),
        xsf::fdtrc(1e-310, 1e308, 1e-310),
        xsf::fdtrc(4.4e-308, 1e308, 1e-310),
        xsf::fdtrc(1e-20, 1.0, 1e-320),
        xsf::fdtrc(1e-10, 1.0, 1e-320),
        xsf::fdtrc(6e-5, 1.0, 1e-315),
        xsf::fdtrc(3e-4, 0.01, 1e-315),
        xsf::fdtrc(0.01, 1e45, 1e-313),
        xsf::fdtr(1e-10, 1.0, 1e-320),
        xsf::fdtr(6e-307, 2.5e-323, 5e-324),
    ];
    let desired = [
        7.263_599_949_158_365e-305,
        7.237_048_719_471_299e-306,
        7.138_593_445_859_812e-308,
        3.127_590_211_779_119_6e-305,
        3.921_326_185_559_873e-18,
        3.806_196_858_474_085_4e-8,
        0.021_850_396_135_364_984,
        0.129_668_149_313_479_72,
        0.973_408_713_466_907_7,
        0.999_999_961_938_031_4,
        4.117_213_715_343_721e-17,
    ];
    xsref::np_assert_allclose!(&p, &desired, rtol = 1e-14);
}

#[test]
fn test_chdtr_fdtr_subnormal() {
    // reference values from mpmath; `x / 2`, `a x` or `w = a x / (b + a x)` is subnormal and rounds
    let p = [
        xsf::chdtr(1.0, f64::from_bits(1)),
        xsf::chdtr(1.0, f64::from_bits(3)),
        xsf::chdtr(0.5, 2.5e-322),
        xsf::chdtr(1.0, 1e-310),
        xsf::chdtr(0.5, 1e-310),
        xsf::fdtr(0.3, 1.0, 1e-320),
        xsf::fdtr(0.3, 1e-13, 1e-320),
        xsf::fdtr(0.7, 1e-10, 3e-319),
    ];
    let p_ref = [
        1.773_504_888_603_627_3e-162,
        3.071_800_574_533_264_4e-162,
        3.696_239_685_446_701e-81,
        7.978_845_608_028_641e-156,
        2.933_738_615_886_127e-78,
        6.990_002_794_838_8e-49,
        2.479_968_749_514_133_6e-59,
        1.311_231_953_921_511_3e-118,
    ];
    xsref::np_assert_allclose!(&p, &p_ref, rtol = 1e-14);
    assert_eq!(xsf::chdtr(2.002, f64::from_bits(3)), f64::from_bits(1));
    // the chi-square limit if `b = inf`
    assert_eq!(xsf::fdtr(1.0, f64::INFINITY, 1e-300), xsf::chdtr(1.0, 1e-300));
    // `x^(a / 2)` underflows, but the subnormal result does not; from mpmath
    let p = [xsf::fdtr(2.2, 1e-290, 1e-310), xsf::ncfdtr(2.2, 1e-290, 0.0, 1e-310)];
    xsref::np_assert_allclose!(&p, &[1.082_037_389_82e-312; 2], rtol = 1e-11);
    // `df / 2` rounds to 0
    let tiny = f64::from_bits(1);
    assert_eq!(xsf::chdtr(tiny, 0.0), 0.0);
    assert_eq!(xsf::chdtrc(tiny, 0.0), 1.0);
    assert_eq!(xsf::fdtr(tiny, 2.0, 0.0), 0.0);
    assert_eq!(xsf::fdtr(tiny, 2.0, 1.0), 1.0);
    assert_eq!(xsf::fdtr(tiny, tiny, 1e-300), 0.5);
    // the limit `B / (A + B)` as `A, B -> 0`, also if neither `a x` nor `w` is subnormal, which
    // mpmath confirms for `w` near 0, 1/2 and 1
    for x in [1e-300, 1.0, 1e300] {
        assert_eq!(xsf::fdtr(3.0 * tiny, 2.0 * tiny, x), 0.4);
        assert_eq!(xsf::fdtrc(3.0 * tiny, 2.0 * tiny, x), 0.6);
    }
    assert_eq!(xsf::fdtr(3.0 * tiny, tiny, 1e308), 0.25);
    assert_eq!(xsf::fdtr(5.0 * tiny, tiny, 1e10), 1.0 / 6.0);
    assert_eq!(xsf::fdtrc(tiny, 7.0 * tiny, 0.5), 0.125);
    assert_eq!(xsf::ncfdtr(3.0 * tiny, 2.0 * tiny, 0.0, 1.0), 0.4);
    let q = xsf::fdtrc(16.0 * tiny, f64::MIN_POSITIVE, 1e-127);
    xsref::np_assert_allclose!(&[q], &[3.552_713_678_800_488_3e-15], rtol = 1e-15);
    // from Arb; the small `fdtr` or `fdtrc` that is `O(df)` for the subnormal `df`
    let pq = [
        xsf::fdtr(1e-300, tiny, 1.0),
        xsf::fdtr(1e-300, tiny, 1e300),
        xsf::fdtrc(tiny, 1e-300, 1.0),
        xsf::fdtrc(tiny, 1e-300, 1e-10),
        xsf::fdtrc(3.0 * tiny, 1e-300, 1.0),
        xsf::fdtrc(tiny, 1e-18, 1.0),
        xsf::fdtrc(1e-310, 1.0, 1.0),
    ];
    let pq_ref = [
        4.940_656_458_412_465e-24,
        4.940_656_458_412_465e-24,
        4.940_656_458_412_465e-24,
        4.940_656_458_412_465e-24,
        1.482_196_937_523_739_7e-23,
        4.940_656_458_412_467e-306,
        3.575_938_365_946_359e-308,
    ];
    xsref::np_assert_allclose!(&pq, &pq_ref, rtol = 1e-15);
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
    // reference values from mpmath (see `incbet_xy`)
    let abx = [(100.0, 1e20, 0.5), (30.0, 1e20, 2.0), (1e20, 30.0, 2.0), (3.0, 1e14, 2.0)];
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
fn test_fdtr_fdtrc_huge_dfn() {
    // `dfn = 1e300`, where `-ln(B(dfn / 2, dfd / 2))` cancels `dfd / 2 ln(1 - w)`; from Arb
    let bx = [(5.0, 0.5), (5.0, 2.0), (0.5, 0.5)];
    let p = [0.075_235_246_146_512_18, 0.776_495_071_123_322_7, 0.153_513_595_808_322_47];
    let q = [0.924_764_753_853_487_8, 0.223_504_928_876_677_3, 0.846_486_404_191_677_5];
    let rtol = 1e-15;
    xsref::np_assert_allclose!(&bx.map(|(b, x)| xsf::fdtr(1e300, b, x)), &p, rtol = rtol);
    xsref::np_assert_allclose!(&bx.map(|(b, x)| xsf::fdtrc(1e300, b, x)), &q, rtol = rtol);
}

#[test]
fn test_fdtr_huge_dfn_rounded_w() {
    // the power of the rounded `w` underflows, unlike that of `1 - (1 - w)`; from Arb
    let p = xsf::fdtr(1.440_000_000_000_000_2e19, 2.0, 0.001_666_666_666_666_666_8);
    xsref::np_assert_allclose!(&[p], &[2.650_396_553_004_479e-261], rtol = 1e-13);
}

#[test]
fn test_fdtr_fdtrc_exact_1mw() {
    // reference values from Arb; the Beta argument (`w` for fdtr, `1 - w` for fdtrc) is in
    // (1/2, 0.95], where its complement is inexact
    let p = [xsf::fdtr(1.76e5, 1e4, 1.031), xsf::fdtr(1.2e7, 1.1e6, 1.0042)];
    let q = [xsf::fdtrc(5.165e4, 2.308e5, 0.987), xsf::fdtrc(2.5e6, 3.3e7, 0.9976)];
    let p_ref = [0.981_513_729_308_534_3, 0.998_524_516_014_640_7];
    let q_ref = [0.971_070_538_485_072_5, 0.995_187_634_483_848_4];
    xsref::np_assert_allclose!(&p, &p_ref, rtol = 1e-15);
    xsref::np_assert_allclose!(&q, &q_ref, rtol = 1e-15);
}

#[test]
fn test_fdtr_fdtrc_near_mean() {
    // reference values from Arb, within ~3 standard deviations of the mean, where the rounding of
    // `w` would be amplified ~`dfn + dfd` times
    let p = [xsf::fdtr(3.4e6, 8.9e6, 0.9987), xsf::fdtr(6.1e6, 7.9e5, 0.99666)];
    let q = [xsf::fdtrc(6.7e6, 3.4e6, 1.00176)];
    let p_ref = [0.074_583_358_906_805_76, 0.023_866_291_181_763_397];
    let q_ref = [0.030_944_149_108_297_337];
    xsref::np_assert_allclose!(&p, &p_ref, rtol = 2e-15);
    xsref::np_assert_allclose!(&q, &q_ref, rtol = 2e-15);
}

xsref::tests! {
    test_fdtri: "fdtri", "d_d_d-d", |x| xsf::fdtri(x[0], x[1], x[2]);
}

#[test]
fn test_fdtri_tiny_p() {
    // reference values from Arb
    let x = [
        xsf::fdtri(10.0, 20.0, 1e-20),
        xsf::fdtri(3.0, 4.0, 1e-15),
        xsf::fdtri(100.0, 50.0, 1e-30),
    ];
    let x_ref = [4.372_813_139_522_666e-5, 7.238_446_978_136_566e-11, 0.061_376_206_282_273_02];
    xsref::np_assert_allclose!(&x, &x_ref, rtol = 1e-14);
}

#[test]
fn test_fdtri_small_p_large_dfn() {
    // reference values from Arb; here `incbet(dfd/2, dfn/2, 1/2) > 1 - p`
    let x = [xsf::fdtri(1e4, 1e3, 1e-14), xsf::fdtri(1e4, 1.0, 1e-14), xsf::fdtri(1e3, 0.5, 1e-9)];
    let x_ref = [0.710_423_069_118_838, 0.016_644_796_655_940_685, 0.014_226_162_949_357_642];
    xsref::np_assert_allclose!(&x, &x_ref, rtol = 1e-14);
}

#[test]
fn test_fdtri_huge_dfn() {
    // reference values from Arb; here `w = dfn x / (dfd + dfn x)` is close to 1
    let x = [
        xsf::fdtri(1e10, 4.0, 1e-14),
        xsf::fdtri(1e8, 10.0, 1e-20),
        xsf::fdtri(1_318_190.8, 1133.5, 4.9e-33),
        xsf::fdtri(1e6, 1e3, 1e-100),
        xsf::fdtri(1e9, 4.0, 1e-8),
        xsf::fdtri(1e9, 2.0, 1e-7),
    ];
    let x_ref = [
        0.055_799_129_486_423_39,
        0.084_358_636_560_782_81,
        0.630_147_551_837_172_9,
        0.439_640_524_172_038_86,
        0.092_868_680_485_509_71,
        0.062_042_067_843_321_7,
    ];
    xsref::np_assert_allclose!(&x, &x_ref, rtol = 1e-15);
}

#[test]
fn test_fdtri_tiny_dfd_p_eps() {
    // reference values from Arb; here `1 - p` rounds to `1 - 2^-53` and
    // `incbet(dfd/2, dfn/2, 1/2)` to 1
    let x = [
        xsf::fdtri(3.0, 1e-16, 1e-16),
        xsf::fdtri(10.0, 1e-16, 1e-16),
        xsf::fdtri(100.0, 1e-16, 1e-16),
    ];
    let x_ref = [4.044_941_625_584_34e-16, 5.425_247_889_451_073e-16, 6.004_433_381_771_332e-16];
    xsref::np_assert_allclose!(&x, &x_ref, rtol = 1e-14);
}

#[test]
fn test_fdtri_tiny_dfd_small_p() {
    // reference values from Arb; here the rounding of `1 - p` changes p by far more than the
    // ~eps / (1 - w) lost in the other branch
    let x = [
        xsf::fdtri(30.0, 1e-12, 1e-13),
        xsf::fdtri(100.0, 1e-14, 1e-15),
        xsf::fdtri(1e3, 1e-10, 1e-13),
        xsf::fdtri(1e4, 1e-9, 1e-12),
    ];
    let x_ref = [
        4.415_753_830_305_983e-13,
        4.639_393_680_922_096e-15,
        1.097_587_864_249_657_5e-11,
        1.103_072_688_087_557_6e-10,
    ];
    xsref::np_assert_allclose!(&x, &x_ref, rtol = 1e-14);
}

xsref::tests! {
    test_gdtr: "gdtr", "d_d_d-d", |x| xsf::gdtr(x[0], x[1], x[2]);
    test_gdtrc: "gdtrc", "d_d_d-d", |x| xsf::gdtrc(x[0], x[1], x[2]);
    test_pdtr: "pdtr", "d_d-d", |x| xsf::pdtr(x[0], x[1]);
    test_pdtrc: "pdtrc", "d_d-d", |x| xsf::pdtrc(x[0], x[1]);
    test_pdtri: "pdtri", "p_d-d", |x| xsf::pdtri(x[0] as i32, x[1]);
    test_bdtr: "bdtr", "d_p_d-d", |x| xsf::bdtr(x[0], x[1] as i32, x[2]);
}

#[test]
fn test_domain_ends() {
    // NaN only outside the domain, like the other inverses
    assert_eq!(xsf::fdtri(3.0, 7.5, 0.0), 0.0);
    assert_eq!(xsf::pdtri(3, 1.0), 0.0);
}

#[test]
fn test_poisson_zero_count() {
    let minsub = f64::from_bits(1);
    let p = 0.367_879_441_171_442_33;
    assert!((xsf::pdtr(0.0, 1.0) - p).abs() <= p * f64::EPSILON);
    assert_eq!(xsf::pdtrc(0.0, minsub), minsub);
    assert_eq!(xsf::pdtrc(0.5, minsub), minsub);

    let m = xsf::pdtri(0, 0.05);
    assert!((xsf::pdtr(0.0, m) - 0.05).abs() <= 0.05 * f64::EPSILON);
}

#[test]
fn test_pdtri_k_max() {
    // `k + 1` overflows `int` in xsf; the median of Gamma(k + 1), from mpmath
    let m = xsf::pdtri(i32::MAX, 0.5);
    xsref::np_assert_allclose!(&[m], &[2_147_483_647.666_666_7], rtol = 1e-15);
}

#[test]
fn test_bdtri_small_y() {
    // mpmath: the p for which bdtr(96, 191, p) == 2.1e-22
    let p = xsf::bdtri(96.0, 191, 2.1e-22);
    assert!((p / 0.814_406_563_996_110_7 - 1.0).abs() < 1e-14, "p = {p}");
}

xsref::tests! {
    test_bdtrc: "bdtrc", "d_p_d-d", |x| xsf::bdtrc(x[0], x[1] as i32, x[2]);
    test_bdtri: "bdtri", "d_p_d-d", |x| xsf::bdtri(x[0], x[1] as i32, x[2]);
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
fn test_poisson_binom_cdf() {
    // SciPy's `poisson_binom_cdf` docstring example
    let p = [0.2, 0.4, 0.6];
    let actual = [-1, 0, 1, 2, 3, 4].map(|k| xsf::poisson_binom_cdf(k, p));
    xsref::np_assert_allclose!(actual, [0.0, 0.192, 0.656, 0.952, 1.0, 1.0], rtol = 1e-15);
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
    assert_eq!(xsf::wilcoxon_cdf_all(3), [0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 1.0]);

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

xsref::tests! {
    test_nbdtr: "nbdtr", "p_p_d-d", |x| xsf::nbdtr(x[0] as i32, x[1] as i32, x[2]);
    test_nbdtrc: "nbdtrc", "p_p_d-d", |x| xsf::nbdtrc(x[0] as i32, x[1] as i32, x[2]);
}

#[test]
fn test_cdf_cvm_inf_smoke() {
    assert!(xsf::cdf_cvm_inf(f64::NAN).is_nan());
    assert!(xsf::cdf_cvm_inf(-1.0).abs() <= f64::EPSILON);

    // reference value from Arb (python-flint) and mpmath
    xsref::np_assert_allclose!(&[xsf::cdf_cvm_inf(0.5)], &[0.960_166_782_434_392_4], rtol = 4e-15);
}

#[test]
fn test_cdf_cvm_inf_near_1() {
    // reference values from mpmath (the series summed to convergence)
    for &(x, p) in &[
        (1.49, 0.999_817_995_788_864_9),
        (3.0, 0.999_999_924_322_565_5),
        (4.0, 0.999_999_999_526_554_7),
        (6.0, 0.999_999_999_999_979_9),
        (7.0, 0.999_999_999_999_999_9),
    ] {
        xsref::np_assert_allclose!(&[xsf::cdf_cvm_inf(x)], &[p], rtol = 4e-15);
    }
}

#[test]
fn test_cdf_cvm_smoke() {
    const N: i32 = 10;
    let lower = 1.0 / (12.0 * f64::from(N));
    let upper = f64::from(N) / 3.0;

    assert!(xsf::cdf_cvm(lower, N).abs() <= f64::EPSILON);
    assert!((xsf::cdf_cvm(upper, N) - 1.0).abs() <= f64::EPSILON);

    // reference value from mpmath (the series of Csorgo and Faraway (1996) summed to convergence)
    xsref::np_assert_allclose!(&[xsf::cdf_cvm(0.5, N)], &[0.962_304_211_246_718_8], rtol = 4e-15);
}

#[test]
fn test_cdf_cvm_finite_n() {
    // reference values from mpmath (the series of Csorgo and Faraway (1996) summed to convergence)
    for &(x, n, p) in &[
        (0.3, 50, 0.864_968_825_999_295_4),
        (1.0, 100, 0.997_617_780_218_287_6),
        (1.5, 10, 0.999_963_160_218_789_8),
        (2.0, 30, 0.999_993_433_827_002_8),
    ] {
        xsref::np_assert_allclose!(&[xsf::cdf_cvm(x, n)], &[p], rtol = 4e-15);
    }
}

#[test]
fn test_cdf_cvm_nonfinite() {
    assert_eq!(xsf::cdf_cvm_inf(f64::INFINITY), 1.0);
    assert_eq!(xsf::cdf_cvm_inf(1e308), 1.0);
    for n in [-1, 10] {
        assert_eq!(xsf::cdf_cvm(f64::INFINITY, n), 1.0);
        assert!(xsf::cdf_cvm(f64::NAN, n).is_nan());
    }
    assert_eq!(xsf::cdf_cvm(1e4, 1_000_000), 1.0);
    // around where `cdf_cvm_inf` reaches 1; mpmath gives 1 - 3e-18 and 1 - 9e-19
    xsref::np_assert_allclose!(&[xsf::cdf_cvm(7.75, 10_000)], &[1.0], rtol = 4e-15);
    assert_eq!(xsf::cdf_cvm(8.0, 10_000), 1.0);
    assert_eq!(xsf::cdf_cvm(1e8, 1_000_000_000), 1.0);
}
