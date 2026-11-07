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
fn test_hyp0f1_nan() {
    assert!(hyp0f1(f64::NAN, 0.0).is_nan());
    assert!(hyp0f1(1.0, c64(f64::NAN, 0.0)).is_nan());
    // pole at v = -inf
    assert!(hyp0f1(f64::NEG_INFINITY, 0.0).is_nan());
    assert!(hyp0f1(f64::NEG_INFINITY, c64(0.0, 0.0)).is_nan());
}
