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
fn test_hyp1f0_f64_inexact_1mz() {
    // `1 - z` rounds to 1 (gh-1613)
    let (a, z) = (2f64.powi(60), 2f64.powi(-60));
    np_assert_allclose!([xsf::hyp1f0(a, z)], [E], rtol = 1e-15);
    np_assert_allclose!([xsf::hyp1f0(-a, z)], [1.0 / E], rtol = 1e-15);
    np_assert_allclose!([xsf::hyp1f0(a, -z)], [1.0 / E], rtol = 1e-15);

    // mpmath: (1 - z)^-a
    let (a, z) = (-4.694_146_575_744_287_4e17, -1.114_329_277_760_014e-16);
    np_assert_allclose!([xsf::hyp1f0(a, z)], [5.214_157_846_371_395e22], rtol = 1e-14);
    let (a, z) = (2.161_946_778_091_621e17, -1.949_623_411_343_651_5e-15);
    np_assert_allclose!([xsf::hyp1f0(a, z)], [8.823_793_205_924_179e-184], rtol = 1e-14);

    // `w^-a` overflows, but `(1 - z)^-a` does not
    let a = -3.624_940_246_312_303_4e17;
    np_assert_allclose!([xsf::hyp1f0(a, z)], [8.463_801_622_536_752e306], rtol = 1e-12);

    // at the over- and underflow thresholds
    assert_eq!(xsf::hyp1f0(-16.0, -2f64.powi(64)), f64::INFINITY);
    assert_eq!(xsf::hyp1f0(19.545_454_545_454_543, -2f64.powi(55)), f64::from_bits(1));
    assert_eq!(xsf::hyp1f0(5.0, 2f64.powi(215)), -f64::from_bits(1));
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

    np_assert_allclose!([xsf::hyp1f0(-1.0, c64(2.0, 0.0))], [c64(-1.0, 0.0)], atol = 1e-15);
    np_assert_allclose!([xsf::hyp1f0(0.0, c64(2.0, 0.0))], [c64(1.0, 0.0)], atol = 1e-15);
    np_assert_allclose!([xsf::hyp1f0(1.0, c64(2.0, 0.0))], [c64(-1.0, 0.0)], atol = 1e-15);
    np_assert_allclose!([xsf::hyp1f0(2.0, c64(2.0, 0.0))], [c64(1.0, 0.0)], atol = 1e-15);
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
    np_assert_allclose!([hyp0f1(1.0, 1.2e5)], [1.171_443_671_889_55e299], rtol = 1e-13);
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
fn test_hyp0f1_asymptotic_exponent() {
    // Debye's expansion, whose exponent `~ln(F)` amplifies the rounding of `sqrt(z)` or of `b - 1`;
    // Arb references
    np_assert_allclose!(
        [hyp0f1(939.25, 20_339.5), hyp0f1(-63.1, 69_347.25)],
        [1_993_204_536.773_062, 9.211_594_514_131_52e293],
        rtol = 5e-16
    );
    // where `v1 (z / v1)` overflows in the double-double division
    np_assert_allclose!([hyp0f1(4.2e307, f64::MAX)], [72.256_461_077_483_01], rtol = 1e-15);
    // ... without changing the rounding where it doesn't (Arb: 7.7303732707275472)
    assert_eq!(hyp0f1(8.79e307, f64::MAX), 7.730_373_270_727_547);
}

#[test]
fn test_hyp0f1_near_poles() {
    // where `b - 1` is not exact, which is amplified near the poles of Gamma(b); Arb references
    np_assert_allclose!(
        [hyp0f1(-31.999_996, -48.41), hyp0f1(-31.000_01, 84.5)],
        [4.723_178_622_128_975, 0.098_158_316_709_055_02],
        rtol = 1e-14
    );
    np_assert_allclose!(
        [hyp0f1(-15.0035, c64(0.83, 18.03)), hyp0f1(-31.999_988_7, c64(25.36, -33.03))],
        [
            c64(0.335_082_562_694_837, -0.835_865_956_742_366_4),
            c64(0.240_757_212_662_643_36, 0.380_338_960_913_804_4),
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
        [hyp0f1(-127.086, c64(1179.69, -1589.74))],
        [c64(0.004_509_268_583_556_574, 0.007_071_343_041_396_766)],
        rtol = 1e-13
    );
    let (b, z) = (-255.878_361_452_585_15, -700_206.011_134_528_7);
    let expected = -1.786_070_161_422_895_2e243;
    np_assert_allclose!([hyp0f1(b, z)], [expected], rtol = 1e-13);
    np_assert_allclose!([hyp0f1(b, c64(z, 0.0))], [c64(expected, 0.0)], rtol = 5e-13);
    np_assert_allclose!(
        [hyp0f1(-511.999_956_9, c64(35_623.46, 0.0))],
        [c64(9.544_129_436_495_515e31, 0.0)],
        rtol = 1e-12
    );
    np_assert_allclose!(
        [hyp0f1(-127.000_000_000_005_58, c64(1958.24, -383.32))],
        [c64(-47_958_952_193.679_64, -54_017_848_410.105_41)],
        rtol = 1e-13
    );
    np_assert_allclose!(
        [hyp0f1(-511.896, c64(7597.31, 0.0))],
        [c64(4.411_526_498_988_145e-7, 0.0)],
        rtol = 2e-13
    );
    let (b, z) = (-511.026, -65_088.5);
    let expected = 5.445_990_179_103_475e68;
    // (rounding `b - 1` costs ~4e-13 here)
    np_assert_allclose!([hyp0f1(b, z)], [expected], rtol = 5e-13);
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
            // the prefactor overflows, but not the result (J_-201.5 itself has 2.7e-12 here)
            hyp0f1(-200.5, c64(-6_419_513.86, 1013.47)),
        ],
        [
            c64(7.033_286_824_840_678e-227, 0.0),
            c64(2.542_862_730_728_592_6e-282, 0.0),
            c64(-7.063_051_849_311_629e226, -1.346_609_895_131_63e225),
            c64(9.749_203_779_137_917e-182, 1.900_960_397_697_128_3e-182),
            c64(1.791_082_883_940_286_4e308, -6.562_886_115_800_78e307),
        ],
        rtol = 1e-11
    );
}

#[test]
fn test_hyp0f1_prefactor_overflow() {
    // `Gamma(b) (-z)^((1-b)/2)` overflows, but not the result; Arb references
    np_assert_allclose!([hyp0f1(-421.55, -688_130.7)], [4.798_975_451_790_942_5e306], rtol = 1e-14);
    np_assert_allclose!([hyp0f1(-49.35, -5.77e14)], [-2.612_137_718_883_196e304], rtol = 1e-12);
}

#[test]
fn test_hyp0f1_complex_huge_b() {
    // past the turning point `4 |z| = (b-1)^2`, where AMOS fails and the result underflows, as
    // `Gamma(b) (-z)^((1-b)/2)` does and `|J_{b-1}| <= 1` (DLMF 10.14.1)
    for (b, z) in [(1e18, -5e35), (3e10, -1e21)] {
        for im in [0.0, -0.0] {
            assert_eq!(hyp0f1(b, c64(z, im)), c64(0.0, 0.0));
        }
    }
    // off the axis, where AMOS fails too, but the result need not underflow (mpmath contour
    // quadrature: 1.45 + 0.198i, 1.18e26 + 1.02e26i, about -5.2e-318 + 3.1e-317i), no wrong 0
    for (b, z) in [
        (3.334_286_344_827_101_4e17, c64(-1.000_514_555_786_581_5e33, 1.501_255_051_533_844_6e34)),
        (3.876_701_454_231_036e19, c64(-1.352_515_771_166_857_8e37, 2.029_426_880_398_685_5e38)),
        (3e15, c64(-8.099_535_869_576_17e28, 1.215_321_563_150_820_3e30)),
    ] {
        let f = hyp0f1(b, z);
        assert!(f.re != 0.0 || f.im != 0.0, "{f}");
    }
    // on the axis for small `b`, where the real path loses the phase (mpmath: 0.314, -5.35e-17,
    // about -1.48e-323), no wrong value
    for (b, z) in [(0.5, -1e32), (1.0, -1e64), (2.65, -1e300)] {
        assert!(hyp0f1(b, c64(z, 0.0)).re.is_nan());
    }
}

#[test]
fn test_hyp0f1_large_b() {
    let bz = [
        // large-v asymptotic expansion
        (300.0, 1.0),
        (-661.1, 0.0624),
        // b near a non-positive integer
        (-9.000_000_000_000_306, -0.0195),
        // over- or underflow of Gamma(b), |z|^((1-b)/2), or J_{b-1}(2 sqrt(-z))
        (200.0, -1.0),
        (233.06, -511.59),
        (873.98, -1963.18),
        (-628.01, -2515.8),
        (1.3e14, 2.3e8),
        (-127.999_999_999_999_96, -0.023),
    ];
    // computed with mpmath
    let expected = [
        1.003_338_876_548_527_4,
        0.999_905_616_320_190_9,
        1.002_169_309_745_669_4,
        0.995_012_417_311_756_4,
        0.110_192_808_171_306_21,
        0.105_489_276_506_703_78,
        55.639_583_267_914_766,
        1.000_001_769_232_334_4,
        1.000_179_703_771_905_4,
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-12);

    // b < 0 and large z > 0, where I_{b-1}(2 sqrt(z)) or Gamma(b) z^((1-b)/2) over- or underflows
    let bz = [(-8_691.5, 8.6e6), (-22_268.5, 6.1e7), (-30_123.5, 1.1e8)];
    // computed with mpmath via DLMF 10.27.2, with K_{1-b} from its integral (DLMF 10.32.9)
    let expected =
        [1.830_168_255_735_138e-246, -4.861_992_754_365_575e201, 2.616_354_312_299_041e60];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-10);
    assert_eq!(hyp0f1(-85.5, 5.5e5), f64::INFINITY);
    assert_eq!(hyp0f1(-300.5, 2e5), f64::NEG_INFINITY);

    // |z| > |b| >> 1, where the asymptotic expansion lost ~eps |b| ln|b| in its exponent
    let bz =
        [(1e12, 2e12), (1e12, -2e12), (-1e12 - 0.5, 2e12), (-1e12 - 0.5, -2e12), (100_000.5, 1e7)];
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
    let bz =
        [(-7.000_001, 0.1), (200.0, 100.0), (300.0, 1000.0), (300.0, -1000.0), (-300.5, -1000.0)];
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
    np_assert_allclose!([hyp0f1(3001.0, 1_500_000.0)], [1.242_442_611_337_837e202], rtol = 1e-12);
}

#[test]
fn test_hyp0f1_turning_point() {
    // below the turning point `z = -(b-1)^2/4`, where the Bessel functions lost ~eps |b| ln|b|
    // through their prefactor
    let bz = [(1501.0, -512_293.0), (-1499.5, -521_411.0)];
    // computed with Arb
    let expected = [4.951_925_098_315_959e-176, 7.295_375_461_642_578e179];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-13);

    // just above `|b-1| p^3 = 50`, with `p = sqrt(1 + 4z/(b-1)^2)`, where the 11 terms of Debye's
    // expansion lost up to ~1e-11
    let bz = [(1001.0, -215_619.0), (-999.5, -215_846.0), (-1999.5, -913_836.0)];
    // computed with Arb
    let expected =
        [1.039_046_359_128_989_3e-109, 2.992_759_281_182_759_5e109, 1.978_895_938_068_044_4e235];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-13);

    // past the turning point, where the Bessel functions lost ~eps |b| ln|b| through their prefactor
    let bz =
        [(732.0, -163_204.0), (1508.0, -683_246.0), (-520.5, -149_822.0), (-696.5, -307_912.0)];
    // computed with Arb
    let expected = [
        -1.775_972_146_185_652_8e-129,
        -1.014_882_538_628_246e-261,
        5.006_641_664_901_066_5e158,
        -3.421_118_017_417_493_8e233,
    ];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-13);

    // further past it, where AMOS `J` lost ~eps sqrt(-z) through its phase
    let bz = [(862.6, -268_538.7), (120.25, -1e7)];
    // computed with Arb
    let expected = [-1.497_621_087_257_290_3e-185, 3.570_732_130_117_078e-223];
    np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 2e-15);
}

#[test]
fn test_hyp0f1_c64_large_b() {
    // Arb; the Bessel functions over- or underflow here
    let bz = [(1000.0, c64(30.0, 40.0)), (-500.5, c64(-20.0, 10.0)), (150.25, c64(-60.0, -45.0))];
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

    // Arb; where 4 terms of Debye's expansion (DLMF 10.41.10) lose 1.7e-13
    np_assert_allclose!(
        [hyp0f1(488.06, c64(-1393.07, 630.22))],
        [c64(0.015_365_699_816_803_597, 0.055_112_429_428_133_845)],
        rtol = 1e-15
    );

    // Arb; b < 0 and |b| < |z| < 6 |b|, where the `K` (or `Y`) term of `I_{b-1}` (or `J_{b-1}`)
    // overflows and Gamma(b) underflows
    let bz = [
        (-438.28, c64(561.95, 38.63)),
        (-454.8, c64(-1117.57, 158.38)),
        (-969.22, c64(-6471.05, 0.0)),
    ];
    let expected = [
        c64(0.276_879_120_063_980_63, -0.024_395_825_299_763_675),
        c64(11.036_922_053_149_937, -4.030_633_238_805_921),
        c64(812.234_227_253_51, 0.0),
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

#[test]
fn test_hyp0f1_small_z() {
    // the terms after a near-zero `v + n` matter
    let v = -1.999_999_999_999_999_8;
    let z = 2f64.powi(-19);
    np_assert_allclose!([hyp0f1(v, z)], [1.002_603_214_235_023_4], rtol = 1e-14);
    np_assert_allclose!(
        [hyp0f1(v, c64(z, 0.0))],
        [c64(1.002_603_214_235_023_4, 0.0)],
        rtol = 1e-14
    );

    // `z * z` underflows
    let v = 2f64.powi(-600);
    np_assert_allclose!([hyp0f1(v, -v)], [1.204_959_932_551_442e-181], rtol = 1e-14);
    np_assert_allclose!(
        [hyp0f1(v, c64(-v, 0.0))],
        [c64(1.204_959_932_551_442e-181, 0.0)],
        rtol = 1e-14
    );
}

#[test]
fn test_hyp0f1_tiny_v_overflow() {
    let v = f64::from_bits(1);

    // only the real part overflows
    let res = hyp0f1(v, c64(-1e-8, 1e-100));
    assert_eq!(res.re, f64::NEG_INFINITY);
    np_assert_allclose!([res.im], [2.024_022_512_832_881e223], rtol = 1e-14);

    // `1 / v` overflows
    assert_eq!(hyp0f1(v, 1e-8), f64::INFINITY);
    assert_eq!(hyp0f1(v, -1e-8), f64::NEG_INFINITY);
}

#[test]
fn test_hyp0f1_tiny_v_unscaled() {
    // `z / v` is exact, unlike `1 / v * z`
    let v = 49.0 * 2f64.powi(-600);
    np_assert_allclose!([hyp0f1(v, -v)], [5.904_303_669_502_066e-180], rtol = 1e-14);
    np_assert_allclose!(
        [hyp0f1(v, c64(-v, 0.0))],
        [c64(5.904_303_669_502_066e-180, 0.0)],
        rtol = 1e-14
    );

    // `1 / v` overflows
    // (`2^-1024` and `2^-1025` are subnormal)
    let (v, expected) = (f64::from_bits(1 << 51), f64::from_bits(1 << 50));
    assert_eq!(hyp0f1(v, -v), expected);
    assert_eq!(hyp0f1(v, c64(-v, 0.0)), c64(expected, 0.0));

    // the first three terms cancel, to about 1e-9 of their size (Arb)
    let (v, z) = (2.828_005_925_933_391_5e-8, -2.828_005_965_921_478_8e-8);
    np_assert_allclose!([hyp0f1(v, z)], [3.842_043_098_418_636e-17], rtol = 1e-7);
}

#[test]
fn test_hyp0f1_scaled_retry() {
    // `v + z` overflows
    np_assert_allclose!([hyp0f1(1e308, 1e308)], [std::f64::consts::E], rtol = 1e-14);
    assert_eq!(hyp0f1(f64::INFINITY, 1.0), 1.0);

    // the scaled `z.im` underflows
    let v = f64::from_bits(1);
    let res = hyp0f1(v, c64(-1e-8, v));
    np_assert_allclose!([res.im], [0.999_999_99], rtol = 1e-14);

    // overflows, although the Bessel result is 2.4e-14 below the maximum
    assert_eq!(hyp0f1(v, c64(2f64.powi(-50), 0.0)).re, f64::INFINITY);
}
