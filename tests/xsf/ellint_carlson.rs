use core::f64::consts::{LN_2, PI};

use num_complex::{Complex64, c64};

const INF: f64 = f64::INFINITY;
const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

// The cases are from `scipy.special.tests.test_basic.TestCarlson` (also used in the xsf tests),
// with the reference values recomputed with mpmath.

#[test]
fn test_elliprc() {
    let actual = [
        xsf::elliprc(0.0, 0.25),
        xsf::elliprc(2.25, 2.0),
        xsf::elliprc(0.25, -2.0), // Cauchy principal value
    ];
    let expected = [PI, LN_2, 0.231_049_060_186_648_45];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    let actual = [
        xsf::elliprc(ZERO, I),
        xsf::elliprc(-I, I),
        xsf::elliprc(I, -ONE),
    ];
    let expected = [
        c64(1.110_720_734_539_591_5, -1.110_720_734_539_591_5),
        c64(1.226_084_956_907_219_9, -0.344_711_369_887_676_8),
        c64(0.777_785_969_204_473_8, 0.198_324_849_934_287_73),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    assert_eq!(xsf::elliprc(1.0, INF), 0.0);
    assert!(xsf::elliprc(1.0, 0.0).is_nan());
    assert_eq!(xsf::elliprc(ONE, c64(1.0, INF)), ZERO);
}

#[test]
fn test_elliprd() {
    let actual = [xsf::elliprd(0.0, 2.0, 1.0), xsf::elliprd(2.0, 3.0, 4.0)];
    let expected = [1.797_210_352_103_388_4, 0.165_105_272_942_610_54];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    let actual = [
        xsf::elliprd(I, -I, c64(2.0, 0.0)),
        xsf::elliprd(ZERO, I, -I),
        xsf::elliprd(ZERO, c64(-1.0, 1.0), I),
        xsf::elliprd(c64(-2.0, -1.0), -I, c64(-1.0, 1.0)),
    ];
    let expected = [
        c64(0.659_338_541_542_197_7, 0.0),
        c64(1.270_819_627_190_968_6, 2.781_112_015_952_057_7),
        c64(-1.857_723_543_923_906, -0.961_934_508_888_385_6),
        c64(1.824_902_739_370_380_5, -1.221_847_578_482_703_6),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    assert_eq!(xsf::elliprd(1.0, 1.0, INF), 0.0);
    assert_eq!(xsf::elliprd(1.0, 1.0, 0.0), INF);
    assert!(xsf::elliprd(1.0, 1.0, -1.0).is_nan());
    assert!(xsf::elliprd(1.0, 1.0, -f64::MIN_POSITIVE / 2.0).is_nan());
    assert_eq!(xsf::elliprd(ZERO, ONE, ZERO).re, INF);
    assert!(xsf::elliprd(ONE, ONE, -ONE).re.is_nan());
}

#[test]
fn test_elliprf() {
    let actual = [xsf::elliprf(0.0, 1.0, 2.0), xsf::elliprf(2.0, 3.0, 4.0)];
    let expected = [1.311_028_777_146_059_8, 0.584_082_841_677_151_7];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    let actual = [
        xsf::elliprf(I, -I, ZERO),
        xsf::elliprf(c64(-1.0, 1.0), I, ZERO),
        xsf::elliprf(I, -I, c64(2.0, 0.0)),
        xsf::elliprf(c64(-1.0, 1.0), I, c64(1.0, -1.0)),
    ];
    let expected = [
        c64(1.854_074_677_301_371_9, 0.0),
        c64(0.796_125_865_842_339_1, -1.213_856_669_836_496),
        c64(1.044_144_565_406_436, 0.0),
        c64(0.939_120_502_186_193_7, -0.532_962_520_186_352_7),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    assert_eq!(xsf::elliprf(1.0, INF, 1.0), 0.0);
    assert_eq!(xsf::elliprf(0.0, 1.0, 0.0), INF);
    assert!(xsf::elliprf(1.0, 1.0, -1.0).is_nan());
    assert_eq!(xsf::elliprf(c64(INF, 0.0), ZERO, ONE), ZERO);
    assert!(xsf::elliprf(ONE, ONE, c64(-INF, 1.0)).re.is_nan());

    // K(m) = R_F(0, 1 - m, 1)
    for m in [-1e300, -1e10, -1.0, 0.0, 0.5, 0.99, 1.0 - 1e-12] {
        xsref::np_assert_allclose!(
            [xsf::elliprf(0.0, 1.0 - m, 1.0)],
            [xsf::ellipk(m)],
            rtol = 1e-15
        );
    }
}

#[test]
fn test_elliprg() {
    let actual = [
        xsf::elliprg(0.0, 16.0, 16.0),
        xsf::elliprg(2.0, 3.0, 4.0),
        xsf::elliprg(0.0, 0.0796, 4.0),
    ];
    let expected = [PI, 1.725_503_028_069_227_8, 1.028_475_809_028_804];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    let actual = [
        xsf::elliprg(ZERO, I, -I),
        xsf::elliprg(c64(-1.0, 1.0), I, ZERO),
        xsf::elliprg(-I, c64(-1.0, 1.0), I),
    ];
    let expected = [
        c64(0.423_606_542_396_989_5, 0.0),
        c64(0.446_605_916_770_183_74, 0.707_683_523_575_153_8),
        c64(0.360_233_921_844_733_1, 0.403_486_234_017_221_1),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    assert_eq!(xsf::elliprg(0.0, 0.0, 0.0), 0.0);
    assert_eq!(xsf::elliprg(1.0, INF, 1.0), INF);
    assert_eq!(xsf::elliprg(c64(INF, 0.0), ONE, ONE).re, INF);

    // E(m) = 2 R_G(0, 1 - m, 1)
    for m in [-1e300, -1e10, -1.0, 0.0, 0.5, 0.99, 1.0] {
        xsref::np_assert_allclose!(
            [2.0 * xsf::elliprg(0.0, 1.0 - m, 1.0)],
            [xsf::ellipe(m)],
            rtol = 1e-15
        );
    }
}

#[test]
fn test_elliprj() {
    let actual = [
        xsf::elliprj(0.0, 1.0, 2.0, 3.0),
        xsf::elliprj(2.0, 3.0, 4.0, 5.0),
        // Cauchy principal values
        xsf::elliprj(2.0, 3.0, 4.0, -0.5),
        xsf::elliprj(2.0, 3.0, 4.0, -5.0),
        // "difficult arguments"
        xsf::elliprj(
            6.483_625_725_195_452e-8,
            1.164_913_652_819_688_6e-27,
            3.676_734_016_716_8e13,
            0.493_704_617_023_468,
        ),
        xsf::elliprj(
            14.375_105_857_849_121,
            9.993_988_969_725_365e-11,
            1.728_442_622_699_44e-26,
            5.898_871_222_598_245e-6,
        ),
    ];
    let expected = [
        0.776_886_237_785_823_3,
        0.142_975_796_671_567_54,
        0.247_238_197_030_515_64,
        -0.127_112_300_429_639_1,
        8.634_269_206_442_419e-6,
        829_774.142_480_162_7,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 5e-15);

    let actual = [
        xsf::elliprj(c64(2.0, 0.0), c64(3.0, 0.0), c64(4.0, 0.0), c64(-1.0, 1.0)),
        xsf::elliprj(I, -I, ZERO, c64(2.0, 0.0)),
        xsf::elliprj(c64(-1.0, 1.0), c64(-1.0, -1.0), ONE, c64(2.0, 0.0)),
        xsf::elliprj(I, -I, ZERO, c64(1.0, -1.0)),
        xsf::elliprj(c64(-1.0, 1.0), c64(-1.0, -1.0), ONE, c64(-3.0, 1.0)),
    ];
    let expected = [
        c64(0.136_139_458_277_705_35, -0.382_075_616_244_271_67),
        c64(1.649_001_166_271_088_4, 0.0),
        c64(0.941_483_588_412_202_4, 0.0),
        c64(1.826_011_522_900_931_7, 1.229_066_190_864_347_2),
        c64(-0.611_279_708_120_281_8, -1.068_403_839_000_680_7),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    assert_eq!(xsf::elliprj(1.0, 1.0, 1.0, INF), 0.0);
    assert_eq!(xsf::elliprj(1.0, 1.0, INF, 1.0), 0.0);
    assert!(xsf::elliprj(-1.0, 1.0, 1.0, 1.0).is_nan());
    assert!(xsf::elliprj(1.0, 0.0, 0.0, 0.0).is_nan());
}

#[test]
fn test_carlson_huge_and_tiny_args() {
    let p = |k| 2.0_f64.powi(k);
    // the tiny imaginary part must not be flushed to 0
    let rd_c = xsf::elliprd(c64(p(400), 0.0), c64(p(400), 0.0), c64(-p(400), p(-1000)));
    let actual = [
        xsf::elliprg(5.17e139, 2.97e11, 1.55e225),
        xsf::elliprj(1e-200, 2e-200, 3e-200, 4e-200),
        // the smallest argument must not become subnormal
        xsf::elliprd(p(300), p(300), 9.332_636_185_032_187e-302),
        // `x * y` underflows
        xsf::elliprj(p(-540), p(-540), p(250), p(-550)),
        rd_c.re,
        rd_c.im,
    ];
    let expected = [
        1.968_501_968_502_952_7e112,
        2.398_480_997_495_677_7e299,
        4.820_814_132_776_971_4e60,
        8.806_214_172_946_52e125,
        -4.015_121_944_122_943e-181,
        -1.361_995_467_683_141e-181,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}
