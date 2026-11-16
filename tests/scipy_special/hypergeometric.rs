use core::f64::consts::{E, PI};
use num_complex::c64;
use xsf::hyp0f1;
use xsref::np_assert_allclose;

#[test]
fn test_hyp0f0_f64() {
    np_assert_allclose!([xsf::hyp0f0(0.0)], [1.0], atol = 1e-15);
    np_assert_allclose!([xsf::hyp0f0(1.0)], [E], atol = 1e-15);
    np_assert_allclose!([xsf::hyp0f0(-1.0)], [1.0 / E], atol = 1e-15);
}

#[test]
fn test_hyp0f0_c64() {
    np_assert_allclose!([xsf::hyp0f0(c64(0.0, 0.0))], [c64(1.0, 0.0)], atol = 1e-15);
    np_assert_allclose!([xsf::hyp0f0(c64(1.0, 0.0))], [c64(E, 0.0)], atol = 1e-15);
    np_assert_allclose!([xsf::hyp0f0(c64(0.0, PI))], [c64(-1.0, 0.0)], atol = 1e-15);
}

#[test]
fn test_hyp1f0_f64() {
    assert_eq!(xsf::hyp1f0(-1.0, -1.0), 2.0);
    assert_eq!(xsf::hyp1f0(0.0, -1.0), 1.0);
    assert_eq!(xsf::hyp1f0(1.0, -1.0), 0.5);
    assert_eq!(xsf::hyp1f0(2.0, -1.0), 0.25);

    assert_eq!(xsf::hyp1f0(-1.0, 0.0), 1.0);
    assert_eq!(xsf::hyp1f0(0.0, 0.0), 1.0);
    assert_eq!(xsf::hyp1f0(1.0, 0.0), 1.0);
    assert_eq!(xsf::hyp1f0(2.0, 0.0), 1.0);

    assert_eq!(xsf::hyp1f0(-1.0, 1.0), 0.0);
    assert_eq!(xsf::hyp1f0(0.0, 1.0), 1.0);
    assert!(xsf::hyp1f0(1.0, 1.0).is_nan());
    assert!(xsf::hyp1f0(2.0, 1.0).is_nan());

    assert_eq!(xsf::hyp1f0(-1.0, 2.0), -1.0);
    assert_eq!(xsf::hyp1f0(0.0, 2.0), 1.0);
    assert_eq!(xsf::hyp1f0(1.0, 2.0), -1.0);
    assert_eq!(xsf::hyp1f0(2.0, 2.0), 1.0);
}

#[test]
fn test_hyp1f0_c64() {
    assert_eq!(xsf::hyp1f0(-1.0, c64(-1.0, 0.0)), c64(2.0, 0.0));
    assert_eq!(xsf::hyp1f0(0.0, c64(-1.0, 0.0)), c64(1.0, 0.0));
    assert_eq!(xsf::hyp1f0(1.0, c64(-1.0, 0.0)), c64(0.5, 0.0));
    assert_eq!(xsf::hyp1f0(2.0, c64(-1.0, 0.0)), c64(0.25, 0.0));

    assert_eq!(xsf::hyp1f0(-1.0, c64(0.0, 0.0)), c64(1.0, 0.0));
    assert_eq!(xsf::hyp1f0(0.0, c64(0.0, 0.0)), c64(1.0, 0.0));
    assert_eq!(xsf::hyp1f0(1.0, c64(0.0, 0.0)), c64(1.0, 0.0));
    assert_eq!(xsf::hyp1f0(2.0, c64(0.0, 0.0)), c64(1.0, 0.0));

    assert_eq!(xsf::hyp1f0(-1.0, c64(1.0, 0.0)), c64(0.0, 0.0));
    assert_eq!(xsf::hyp1f0(0.0, c64(1.0, 0.0)), c64(1.0, 0.0));
    assert!(xsf::hyp1f0(1.0, c64(1.0, 0.0)).is_nan());
    assert!(xsf::hyp1f0(2.0, c64(1.0, 0.0)).is_nan());

    np_assert_allclose!(
        [xsf::hyp1f0(-1.0, c64(2.0, 0.0))],
        [c64(-1.0, 0.0)],
        atol = 1e-15
    );
    np_assert_allclose!(
        [xsf::hyp1f0(0.0, c64(2.0, 0.0))],
        [c64(1.0, 0.0)],
        atol = 1e-15
    );
    np_assert_allclose!(
        [xsf::hyp1f0(1.0, c64(2.0, 0.0))],
        [c64(-1.0, 0.0)],
        atol = 1e-15
    );
    np_assert_allclose!(
        [xsf::hyp1f0(2.0, c64(2.0, 0.0))],
        [c64(1.0, 0.0)],
        atol = 1e-15
    );
}

/// Corresponds to `scipy.special.tests.test_basic.TestHyper.test_hyp0f1`
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/tests/test_basic.py#L3144-L3168>
#[test]
fn test_hyp0f1() {
    // scalar input
    np_assert_allclose!([hyp0f1(2.5, 0.5)], [1.214_827_026_899_997], rtol = 1e-12);
    np_assert_allclose!([hyp0f1(2.5, 0.0)], [1.0], rtol = 1e-15);

    // float input, expected values match mpmath
    let x = [-1.5, -1.0, 0.0, 1.0, 1.5].map(|z| hyp0f1(3.0, z));
    let expected = [
        0.584_936_592_291_43,
        0.705_668_057_231_27,
        1.0,
        1.377_896_895_397_47,
        1.603_736_852_884_80,
    ];
    np_assert_allclose!(x, expected, rtol = 1e-12);

    // complex input
    let x = [-1.5, -1.0, 0.0, 1.0, 1.5].map(|re| hyp0f1(3.0, c64(re, 0.0)));
    np_assert_allclose!(x, expected.map(|v| c64(v, 0.0)), rtol = 1e-12);
}

/// Corresponds to `scipy.special.tests.test_basic.TestHyper.test_hyp0f1_gh5764`
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/tests/test_basic.py#L3170-L3176>
/// which looks like a regression test for scipy/scipy#5764
#[test]
fn test_hyp0f1_gh5764() {
    // Just checks the point that failed; there's a more systematic test in test_mpmath
    let res = hyp0f1(0.8, c64(0.5, 0.5));
    // The expected value was generated using mpmath
    np_assert_allclose!(
        [res],
        [c64(1.613_971_977_644_111_5, 0.808_930_540_617_906_7)],
        atol = 1.5e-7
    );
}

#[test]
fn test_hyp0f1_b_one() {
    // hyp0f1(1, z) = I_0(2 sqrt(z)), which overflows for z above about 1.2744e5
    np_assert_allclose!(
        [hyp0f1(1.0, 1.2e5)],
        [1.171_443_671_889_55e299],
        rtol = 1e-13
    );
    assert_eq!(hyp0f1(1.0, 2e5), f64::INFINITY);
    assert_eq!(hyp0f1(1.0, 1e300), f64::INFINITY);
}

#[test]
fn test_hyp0f1_large_negative_z() {
    // where `sqrt(-z)` is not exact, which shifts the phase of `J`; Arb references
    np_assert_allclose!(
        [hyp0f1(2.5, -123_456.789), hyp0f1(0.25, -1e30)],
        [-3.348_689_981_141_865_5e-6, -923.938_185_812_732],
        rtol = 1e-14
    );
}

#[test]
fn test_hyp0f1_complex_large_z() {
    // where `sqrt(-z)` is not exact, which shifts the argument of `J`; Arb references
    np_assert_allclose!(
        [hyp0f1(0.5, c64(-1e8, 10.0)), hyp0f1(2.3, c64(-1e6, 0.5))],
        [
            c64(0.813_200_097_194_249_9, 0.000_581_984_859_012_089_5),
            c64(1.669_670_181_234_794_6e-6, -1.009_120_641_436_902_8e-9),
        ],
        rtol = 1e-14
    );
}

#[test]
fn test_hyp0f1_near_poles() {
    // where `b - 1` is not exact, which is amplified near the poles of Gamma(b); Arb references
    np_assert_allclose!(
        [
            hyp0f1(-31.999_996_483_927_095, -48.410_101_072_009_944),
            hyp0f1(-31.000_010_454_454_72, 84.215_513_283_229_77)
        ],
        [4.723_194_235_861_527, 0.095_074_906_605_175_92],
        rtol = 1e-14
    );
    np_assert_allclose!(
        [
            hyp0f1(
                -15.003_121_251_499_058,
                c64(0.828_424_759_977_544_1, 18.026_879_250_310_117)
            ),
            hyp0f1(
                -31.999_988_706_823_49,
                c64(25.359_285_755_766_876, -33.027_826_633_854_495)
            ),
        ],
        [
            c64(0.335_420_362_037_963_1, -0.835_851_607_783_394_4),
            c64(0.240_787_837_165_309_3, 0.380_332_237_847_745_34),
        ],
        rtol = 1e-14
    );

    // where the recurrence from the orders `b` and `b + 1` would be less accurate: `b + 1` rounds
    // to 1, its two terms cancel, `2 sqrt(|z|)` is past the turning point `1 - b`, or the `sin(pi b)`
    // term of `I_{b-1}` is negligible, or the gain is too small to be reliable; Arb references
    // (here 2 sqrt(-z) is about the first zero of J_1, so rounding `b - 1` to -1 costs ~14%)
    let z = -12.304_614_080_423_65;
    for actual in [hyp0f1(1e-16, z), hyp0f1(1e-16, c64(z, 0.0)).re] {
        np_assert_allclose!([actual], [2.197_705_466_158_317_3], rtol = 0.2);
    }
    np_assert_allclose!(
        [hyp0f1(
            -127.086_052_907_759_45,
            c64(1_179.687_644_608_810_6, -1_589.740_742_745_845_2)
        )],
        [c64(0.004_507_218_203_138_257, 0.007_064_843_216_650_686)],
        rtol = 1e-13
    );
    let (b, z) = (-255.878_361_452_585_15, -700_206.011_134_528_7);
    let expected = -1.786_070_161_422_895_2e243;
    np_assert_allclose!([hyp0f1(b, z)], [expected], rtol = 1e-13);
    np_assert_allclose!([hyp0f1(b, c64(z, 0.0))], [c64(expected, 0.0)], rtol = 5e-13);
    np_assert_allclose!(
        [hyp0f1(
            -511.999_956_923_477_9,
            c64(35_623.457_087_526_804, 0.0)
        )],
        [c64(9.548_881_977_270_142e31, 0.0)],
        rtol = 1e-12
    );
    np_assert_allclose!(
        [hyp0f1(
            -127.000_000_000_005_58,
            c64(1_958.239_202_491_229_4, -383.322_084_776_425_24)
        )],
        [c64(-47_965_934_161.306_46, -54_008_764_700.608_65)],
        rtol = 1e-13
    );
    np_assert_allclose!(
        [hyp0f1(
            -511.895_066_714_241_5,
            c64(7_597.308_951_474_139, 0.0)
        )],
        [c64(4.411_420_808_676_181_4e-7, 0.0)],
        rtol = 2e-13
    );
    let (b, z) = (-511.025_334_625_666_97, -65_088.852_030_151_46);
    let expected = 5.575_669_507_071_876e68;
    np_assert_allclose!([hyp0f1(b, z)], [expected], rtol = 2e-13);
    np_assert_allclose!([hyp0f1(b, c64(z, 0.0))], [c64(expected, 0.0)], rtol = 2e-13);
}

#[test]
fn test_hyp0f1_complex_prefactor_overflow() {
    // where Gamma(b) or |z|^((1-b)/2) over- or underflows; Arb references
    np_assert_allclose!(
        [
            hyp0f1(200.0, c64(-1e6, 0.0)),
            hyp0f1(60.5, c64(-1.2e12, 0.0)),
            hyp0f1(-200.5, c64(-1e6, 10.0)),
            // the prefactor underflows, but not the result
            hyp0f1(100.0, c64(-1e12, 6e8)),
            // the prefactor overflows, but not the result (J_-201.5 itself has 3.7e-12 here)
            hyp0f1(-200.5, c64(-6_419_513.856_144_281, 1_013.470_386_041_489_2)),
        ],
        [
            c64(7.033_286_824_840_678e-227, 0.0),
            c64(2.542_862_730_728_592_6e-282, 0.0),
            c64(-7.063_051_849_311_629e226, -1.346_609_895_131_63e225),
            c64(9.749_203_779_137_917e-182, 1.900_960_397_697_128_3e-182),
            c64(1.791_080_341_627_737_6e308, -6.562_898_362_614_055e307),
        ],
        rtol = 1e-11
    );
}

#[test]
fn test_hyp0f1_large_b() {
    let bz = [
        // large-v asymptotic expansion
        (300.0, 1.0),
        (-661.112_997_069_224_4, 0.062_385_007_293_382_27),
        // b near a non-positive integer
        (-9.000_000_000_000_306, -0.019_499_612_460_819_94),
        // over- or underflow of Gamma(b), |z|^((1-b)/2), or J_{b-1}(2 sqrt(-z))
        (200.0, -1.0),
        (233.063_303_394_596_46, -511.589_207_216_921_75),
        (873.979_261_659_223_8, -1_963.183_700_462_148_6),
        (-628.010_908_208_124_6, -2_515.795_365_646_08),
        (130_356_679_484_338.69, 229_453_612.840_189_75),
        (-127.999_999_999_999_96, -0.023),
    ];
    // computed with mpmath
    let expected = [
        1.003_338_876_548_527_4,
        0.999_905_640_851_462_9,
        1.002_169_266_580_656_6,
        0.995_012_417_311_756_4,
        0.110_196_664_254_785_33,
        0.105_488_627_754_553_3,
        55.638_844_554_215_72,
        1.000_001_760_199_904_5,
        1.000_179_703_771_905_4,
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-12);

    // b < 0 and large z > 0, where I_{b-1}(2 sqrt(z)) or Gamma(b) z^((1-b)/2) over- or underflows
    let bz = [(-8_691.5, 8.6e6), (-22_268.5, 6.1e7), (-30_123.5, 1.1e8)];
    // computed with mpmath via DLMF 10.27.2, with K_{1-b} from its integral (DLMF 10.32.9)
    let expected = [
        1.830_168_255_735_138e-246,
        -4.861_992_754_365_575e201,
        2.616_354_312_299_041e60,
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-10);
    assert_eq!(hyp0f1(-85.5, 5.5e5), f64::INFINITY);
    assert_eq!(hyp0f1(-300.5, 2e5), f64::NEG_INFINITY);

    // |z| > |b| >> 1, where the asymptotic expansion lost ~eps |b| ln|b| in its exponent
    let bz = [
        (1e12, 2e12),
        (1e12, -2e12),
        (-1e12 - 0.5, 2e12),
        (-1e12 - 0.5, -2e12),
        (100_000.5, 1e7),
    ];
    // computed with mpmath from the Maclaurin series
    let expected = [
        7.389_056_098_915_872,
        0.135_335_283_236_342_03,
        0.135_335_283_237_018_7,
        7.389_056_098_938_039,
        2.555_911_194_624_412_4e43,
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-14);
    // close to overflow (computed with Arb), and with Gamma(b) = inf
    np_assert_allclose!(
        [hyp0f1(-10_001.5, 12_872_883.657_095_07)],
        [1.116_996_244_411_357_3e308],
        rtol = 1e-11
    );
    assert_eq!(hyp0f1(1e-20, 2e5), f64::INFINITY);

    // near zeros of J_{b-1}(2 sqrt(-z)), where it rounds to 0 (instead of 5.2e-17 and 6.1e-15)
    assert!(hyp0f1(1.674_445_551_227_833_6, -2.866_969_344_179_029).abs() < 1e-14);
    assert!(hyp0f1(-3.854_390_408_608_444_4, -5.550_754_014_064_338).abs() < 1e-14);

    // where the Bessel functions lost digits, through their prefactor or near poles
    let bz = [
        (-7.000_001, 0.1),
        (200.0, 100.0),
        (300.0, 1000.0),
        (300.0, -1000.0),
        (-300.5, -1000.0),
    ];
    // computed with Arb
    let expected = [
        0.985_832_546_664_187_5,
        1.647_699_631_540_315_4,
        27.526_315_983_675_59,
        0.035_011_885_022_349_28,
        28.404_724_822_364_255,
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-14);
    // Debye's expansion with b near an integer (with the J term), and with u_1(p) = 0
    np_assert_allclose!(
        [hyp0f1(-999.999_999_999_999_9, -216_000.0)],
        [3.099_721_235_145_35e109],
        rtol = 1e-11
    );
    np_assert_allclose!(
        [hyp0f1(3001.0, 1_500_000.0)],
        [1.242_442_611_337_837e202],
        rtol = 1e-12
    );
}

#[test]
fn test_hyp0f1_c64_large_b() {
    // Arb; the Bessel functions over- or underflow here
    let bz = [
        (1000.0, c64(30.0, 40.0)),
        (-500.5, c64(-20.0, 10.0)),
        (150.25, c64(-60.0, -45.0)),
    ];
    let expected = [
        c64(1.029_630_689_563_530_5, 0.041_205_970_834_947_474),
        c64(1.040_562_667_171_689, -0.020_794_893_995_743_39),
        c64(0.640_602_115_953_318_3, -0.198_368_510_709_189_82),
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-15);

    // Arb; |z| > b, where the Bessel function underflows (Debye's expansion, scipy/xsf#557)
    let bz = [
        (1000.0, c64(2000.0, 100.0)),
        (1000.0, c64(-2000.0, 100.0)),
        (1e12, c64(2e12, 1e11)),
        (1e14, c64(2e14, 1e13)),
    ];
    let expected = [
        c64(7.337_688_264_623_645, 0.734_749_622_069_595),
        c64(0.134_387_641_620_190_2, 0.013_510_970_206_071_45),
        c64(7.352_141_595_894_477, 0.737_674_716_148_388_2),
        c64(7.352_141_595_908_852, 0.737_674_716_151_300_9),
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-14);

    // Arb; on the branch cut, where `J_{b-1}` rounds to 0 near a zero
    np_assert_allclose!(
        [hyp0f1(1.400_800_854_970_417_8, c64(-2.25, 0.0))],
        [c64(5.471_961_797_816_518e-17, 0.0)],
        atol = 1e-16
    );

    // Arb; the norm of the partial sums overflows, or the partial sums themselves
    let actual = [c64(0.0, 1.0), c64(0.2, 0.9)].map(|z| hyp0f1(6e-309, z));
    let expected = [
        c64(-8.217_785_411_823_984e307, 1.528_356_022_306_726_8e308),
        c64(-3.690_208_731_438_365e307, 1.707_521_906_170_721_5e308),
    ];
    np_assert_allclose!(actual, expected, rtol = 1e-14);
}

#[test]
fn test_hyp0f1_nan() {
    assert!(hyp0f1(f64::NAN, 0.0).is_nan());
    assert!(hyp0f1(1.0, c64(f64::NAN, 0.0)).is_nan());
    // pole at v = -inf
    assert!(hyp0f1(f64::NEG_INFINITY, 0.0).is_nan());
    assert!(hyp0f1(f64::NEG_INFINITY, c64(0.0, 0.0)).is_nan());
}
