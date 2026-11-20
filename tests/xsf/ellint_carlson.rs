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

    let actual = [xsf::elliprc(ZERO, I), xsf::elliprc(-I, I), xsf::elliprc(I, -ONE)];
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
fn test_elliprc_subnormal_huge() {
    let actual = [
        xsf::elliprc(1.0, 5e-324),
        xsf::elliprc(1.0, -5e-324), // Cauchy principal value
        xsf::elliprc(0.0, 1e-320),
        xsf::elliprc(1e-320, 1e-320),
        xsf::elliprc(1e300, 1e308),
    ];
    // Arb
    let expected = [
        372.913_183_141_250_6,        // 372.91318314125057647
        372.913_183_141_250_6,        // 372.91318314125057647
        1.570_805_070_562_178e160,    // 1.5708050705621779667e160
        1.000_005_566_455_136_3e160,  // 1.0000055664551362865e160
        1.570_696_334_648_211_6e-154, // 1.5706963346482116368e-154
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // Arb: 286.75424458785220432 - 118.77749717634007236i
    let actual = xsf::elliprc(c64(1.0, 1.0), c64(1e-320, 1e-320));
    let expected = c64(286.754_244_587_852_2, -118.777_497_176_340_07);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);

    // Cauchy principal value, as the imaginary part of y is subnormal
    // Arb: 6.8720984451488491504e149 + 4.6016926624554902025e148i
    let actual = xsf::elliprc(c64(1e-300, 1e-300), c64(-1e-300, 1e-310));
    let expected = c64(6.872_098_445_148_849e149, 4.601_692_662_455_49e148);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
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
fn test_elliprd_subnormal() {
    let actual = [
        xsf::elliprd(1.0, 1.0, 5e-324),
        xsf::elliprd(1e10, 1e20, 1e-320),
        xsf::elliprd(5e-324, 0.0, 1.0),
        xsf::elliprd(0.0, 1e-310, 1.0),
    ];
    // Arb
    let expected = [
        1.349_674_138_362_959e162, // 1.3496741383629589148e162
        3.000_016_699_365_409e145, // 3.0000166993654088595e145
        1_117.818_990_965_431_6,   // 1117.8189909654315653
        1_071.860_951_325_591,     // 1071.8609513255909195
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // Arb: 9.5436373562855885717e161 - 9.5436373562855885717e161i
    let actual = xsf::elliprd(ONE, I, c64(5e-324, 0.0));
    let expected = c64(9.543_637_356_285_59e161, -9.543_637_356_285_59e161);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);

    // Arb: 6.0720675992193185506e323
    assert_eq!(xsf::elliprd(0.0, 1.0, 5e-324), INF);
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
fn test_elliprf_yz_underflow() {
    let actual = [
        xsf::elliprf(0.0, 1e-300, 1e-70),
        // with a correction of -4e-13 relative
        xsf::elliprf(1e-320, 1e-300, 1e-70),
        xsf::elliprf(0.0, 5e-324, 1e-20),
    ];
    // Arb
    let expected = [
        2.661_835_800_554_351_6e37, // 2.6618358005543514484e37
        2.661_835_800_553_351_4e37, // 2.6618358005533514540e37
        3_505_804_793_918.700_7,    // 3505804793918.7007452
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // Arb: 2.6601029326029515852e37 - 3.9269908169872415566e34i
    let actual = xsf::elliprf(ZERO, c64(1e-300, 1e-300), c64(1e-70, 0.0));
    let expected = c64(2.660_102_932_602_951_7e37, -3.926_990_816_987_241_3e34);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
}

#[test]
fn test_elliprf_tiny_x_branch() {
    // the correction for x is sqrt(x) / (sqrt(y) sqrt(z)), not sqrt(x / (y z)), if
    // arg(y) + arg(z) > pi; mpmath (dps = 50): 109.73848232833850925 - 264.93213234591670474i
    let actual = xsf::elliprf(c64(1e-312, 0.0), c64(-1e-295, 1e-295), c64(-1.0, 1.0));
    let expected = c64(109.738_482_328_338_5, -264.932_132_345_916_7);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
}

#[test]
fn test_elliprf_subnormal_huge() {
    let actual = [
        xsf::elliprf(5e-324, 5e-324, 1.0),
        xsf::elliprf(0.0, 5e-324, 1.0),
        xsf::elliprf(1e-310, 1e-309, 1.0),
        xsf::elliprf(1e-310, 1e-306, 1.0),
        xsf::elliprf(1e308, 1e308, 1e308),
    ];
    // Arb
    let expected = [
        372.913_183_141_250_6, // 372.91318314125057647
        373.606_330_321_810_5, // 373.60633032181052178
        356.860_921_336_291_6, // 356.86092133629160316
        353.671_863_258_355_7, // 353.67186325835571219
        1e-154,                // 9.9999999999999999451e-155
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // Arb: 356.91746341227886450 - 0.60464862014200434676i
    let actual = xsf::elliprf(c64(1e-310, 0.0), c64(0.0, 1e-309), ONE);
    let expected = c64(356.917_463_412_278_9, -0.604_648_620_142_004_3);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
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
fn test_carlson_agm_opposite_sides() {
    // `sqrt(y) + sqrt(z)` cancels, as y and z are on opposite sides of the negative real axis
    let actual = [
        xsf::elliprg(ZERO, c64(-0.043_001, -4.6e-11), c64(-0.043, 4.6e-11)),
        xsf::elliprf(ZERO, c64(-32.000_003, -1.3e-9), c64(-32.0, 1.3e-9)),
    ];
    // Arb
    let expected = [
        // -1.1862926450000208018 - 0.16285503760200789828i
        c64(-1.186_292_645_000_020_9, -0.162_855_037_602_007_9),
        // 3.3508415252610557203 + 0.27752697069943702297i
        c64(3.350_841_525_261_056, 0.277_526_970_699_437_05),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_elliprg_negligible_imag_huge() {
    // a negligible imaginary part doesn't prevent downscaling the huge arguments
    let actual = xsf::elliprg(
        c64(1.735_778_648_751_103e186, -3.309_532_202_948_146e-247),
        c64(1.106_380_492_215_430_7e198, 2.695_124_182_909_526_7e-261),
        c64(3.913_736_070_965_839e206, -6.593_607_292_907_775e-23),
    );
    let expected = c64(9.891_582_519_574_452e102, -8.332_346_327_267_998e-127);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
}

#[test]
fn test_carlson_huge_near_cut() {
    // tiny imaginary parts next to a negative real part are not negligible
    let p = |k| 2.0_f64.powi(k);
    let [x, y, z] = [c64(-p(450), p(-797)), c64(-p(450), -p(-797)), c64(1.0, f64::MIN_POSITIVE)];
    let actual = [xsf::elliprf(x, y, z), xsf::elliprd(x, y, z), xsf::elliprg(x, y, z)];
    let expected = [
        c64(1.606_891_115_915_371_2e-65, 0.0),
        c64(1.031_865_770_122_304_8e-135, 0.0),
        c64(-2.330_509_388_214_315e70, 0.0),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);

    let x = c64(p(260), f64::from_bits(1));
    let z = c64(p(260), 0.0);
    let actual = [
        xsf::elliprj(x, ZERO, z, c64(-p(260), 1.0)),
        xsf::elliprj(x, ZERO, z, c64(-p(260), -1.0)),
    ];
    let expected = [
        c64(-9.343_569_646_112_14e-118, -9.343_569_646_112_14e-118),
        c64(-9.343_569_646_112_14e-118, 9.343_569_646_112_14e-118),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_carlson_huge_negligible_imag_finite() {
    // negligible imaginary parts are flushed only if the result is not finite otherwise (Arb)
    let p = |k| 2.0_f64.powi(k);
    let z = c64(p(-826), p(-880));
    let h = c64(p(450), 0.0);
    let actual = [xsf::elliprc(h, z), xsf::elliprf(z, z, h), xsf::elliprd(z, h, z)];
    let expected = [
        c64(8.214_427_384_559_377e-66, -5.147_557_589_468_029e-85),
        c64(8.214_427_384_559_377e-66, -5.147_557_589_468_029e-85),
        c64(1.244_854_670_664_297_9e181, -6.910_331_588_418_35e164),
    ];
    for (a, e) in actual.iter().zip(expected) {
        xsref::np_assert_allclose!([a.re, a.im], [e.re, e.im], rtol = 1e-15);
    }

    // just below 1.5 times the smallest subnormal
    let z = c64(p(550), p(-550));
    assert_eq!(xsf::elliprj(z, z, z, c64(p(800), 0.0)).re, f64::from_bits(1));

    // NaN unscaled, so retried with the imaginary part flushed (Arb)
    let h = c64(p(400), 0.0);
    let res = xsf::elliprj(h + c64(0.0, f64::from_bits(1)), h, h, h * 2.0);
    xsref::np_assert_allclose!([res.re], [1.551_519_687_348_157_6e-181], rtol = 1e-15);
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
        // `x = y` are both real and conjugate (Arb)
        xsf::elliprj(ONE, ONE, c64(2.0, 0.0), c64(-1.0, 1.0)),
    ];
    let expected = [
        c64(0.136_139_458_277_705_35, -0.382_075_616_244_271_67),
        c64(1.649_001_166_271_088_4, 0.0),
        c64(0.941_483_588_412_202_4, 0.0),
        c64(1.826_011_522_900_931_7, 1.229_066_190_864_347_2),
        c64(-0.611_279_708_120_281_8, -1.068_403_839_000_680_7),
        c64(0.104_049_921_394_109_93, -0.870_513_482_495_226_4),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // based on scipy/xsf#568: equal positive real arguments (mpmath)
    let half = c64(0.5, 0.0);
    let actual = [
        xsf::elliprj(half, half, half, c64(-2.0, 0.1)),
        xsf::elliprj(ONE, ONE, ONE, c64(-1.0, -1e-3)),
    ];
    let expected = [
        c64(-1.259_664_965_962_783_9, -1.226_063_044_216_320_5),
        c64(-0.563_913_158_248_595, 1.665_754_193_026_62),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);

    assert_eq!(xsf::elliprj(1.0, 1.0, 1.0, INF), 0.0);
    assert_eq!(xsf::elliprj(1.0, 1.0, INF, 1.0), 0.0);
    assert!(xsf::elliprj(-1.0, 1.0, 1.0, 1.0).is_nan());
    assert!(xsf::elliprj(1.0, 0.0, 0.0, 0.0).is_nan());
}

#[test]
fn test_elliprj_nearly_imaginary_duplicated() {
    // the duplication step keeps the real parts nonnegative where they cancel
    let actual = xsf::elliprj(
        c64(9.916_661_834_316_465e-77, -5.966_600_187_505_507e238),
        c64(2.848_126_752_962_884_6e-29, 4.000_597_101_020_285_5e-234),
        c64(4.782_444_515_529_445e-116, -5.405_461_544_153_664e21),
        c64(9.584_517_345_332_541e93, 4.512_323_928_209_932_5e31),
    );
    let expected = c64(7.528_467_325_907_695e-212, 7.670_796_004_575_693e-212);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
}

#[test]
fn test_carlson_huge_negligible_real() {
    // a negligible real part doesn't prevent downscaling the huge arguments (Arb)
    let (s, e) = (2.0_f64.powi(620), 2.0_f64.powi(-80));
    let actual = [
        xsf::elliprg(
            c64(3.160_931_853_621_509e-222, 2.546_737_660_160_885_7e232),
            c64(2.770_788_258_627_483_4e-250, 8.876_654_259_025_239e205),
            c64(3.515_741_291_918_735e211, 9.916_737_044_917_146e-127),
        ),
        xsf::elliprg(c64(1e-308, s), c64(0.0, e * s), c64(s, -e * s)) / s.sqrt(),
    ];
    let expected = [
        c64(5.642_182_268_591_743_5e115, 5.642_182_268_591_743_5e115),
        c64(0.641_920_478_949_122_3, 0.265_892_168_345_759_3),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // but not if a real part is negative, where the terms of the arguments on both sides of the
    // cut cancel (Arb: 2^310 (-21.4746 + 9.20807i))
    let (y, z) = (c64(-s, e * s), c64(-s, -e * s));
    assert!(xsf::elliprg(c64(1e-308, s), y, z).is_nan());
}

#[test]
fn test_elliprj_nearly_imaginary_p() {
    // NaN for imaginary `p`, where the first `p + lambda` is 0 here (the integral is
    // 0.10806392013989438 - 0.0030081533140647775i)
    let i = |b| c64(0.0, b);
    assert!(xsf::elliprj(i(1.0), i(4.0), i(9.0), i(-11.0)).is_nan());

    // also with a tiny `Re p < 0`, which blocks rescaling
    let s = 2.0_f64.powi(300);
    let (z, p) = (c64(1e-15, 9.0 * s), c64(-1e-320, -11.0 * s));
    assert!(xsf::elliprj(i(s), i(4.0 * s), z, p).is_nan());

    // NaN for `x, y` that are not conjugate, even if their real parts are negligible
    let s = 2.0_f64.powi(450);
    let (x, y) = (c64(1e-308, s), c64(4e-308, -s));
    assert!(xsf::elliprj(x, y, c64(s, 0.0), c64(-s, s)).is_nan());

    // also if flushing them would make them conjugate (Arb: 2.42484e-163 - 3.23129e-163i)
    let s = 2.0_f64.powi(360);
    let (x, y) = (c64(2.0_f64.powi(-1024), s), c64(2.0_f64.powi(-1018), -s));
    assert!(xsf::elliprj(x, y, c64(s, 0.0), i(s)).is_nan());

    // NaN if `p + lambda` cancels in a duplication step for a wide span, where the `R_C` argument
    // underflows (Arb: 7.33688e-218 - 5.51996e-218i)
    let actual = xsf::elliprj(
        c64(1e-308, -4.361_399_083_161_572_3e133),
        c64(1e-170, -1e-33),
        c64(0.0, -5.554_010_425_978_977e151),
        c64(2.0_f64.powi(-582), 4.921_712_707_964_019e142),
    );
    assert!(actual.is_nan());

    // but not if it underflows without cancellation (Arb: 4.17416e-438)
    let a = 2.0_f64.powi(976);
    assert_eq!(xsf::elliprj(a, a, a, f64::from_bits(1)), 0.0);

    // a negligible `Re x` that blocks rescaling, also if `p + lambda` cancels in the first step
    // (with `Im p` 1 ULP toward 0 too), or in a later one (Arb)
    let s = 2.0_f64.powi(360);
    let x = c64(1e-308, s);
    let im = -11.0 * s;
    let actual = [
        xsf::elliprj(x, i(4.0 * s), i(9.0 * s), c64(1.0, im)),
        xsf::elliprj(x, i(4.0 * s), i(9.0 * s), c64(1.0, f64::from_bits(im.to_bits() - 1))),
        xsf::elliprj(x, i(s), i(s), c64(1.0, -15.0 * s)),
        // a huge `Im x` next to tiny real parts (#1655)
        xsf::elliprj(
            c64(3.252_281_333_249_279e-86, 9.628_653_385_766_011e269),
            c64(1.804_605_791_678_615_8e-256, 1.821_374_505_004_775_8e-95),
            c64(5.704_000_764_113_375e-244, -2.547_345_192_149_236_5e-279),
            c64(5.398_519_435_655_003e34, -3.613_203_968_954_653e-58),
        ),
    ];
    let expected = [
        c64(3.002_500_299_932_586_4e-164, -8.357_998_873_287_521e-166),
        c64(3.002_500_299_932_587e-164, -8.357_998_873_287_536e-166),
        c64(4.895_118_724_650_793e-164, 2.001_918_454_026_923_3e-164),
        c64(5.965_420_762_787_378e-168, -6.028_323_376_004_016e-168),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);

    // a wide span with a tiny `Re p > 0` (Wolfram, Arb, mpmath, and the integral of
    // https://dlmf.nist.gov/19.16.E2)
    let actual = xsf::elliprj(
        c64(1.043_155_444_981_761_3e-88, -5.565_490_407_805_559_5e135),
        c64(1.142_710_787_031_183_8e-145, 3.370_965_858_974_056e-27),
        c64(7.007_098_407_015_658e141, 1.559_315_032_887_089_7e-132),
        c64(1e-300, 2.694_388_702_830_613e-45),
    );
    let expected = c64(1.781_293_768_892_485_5e-124, 0.0);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-15);
}

#[test]
fn test_elliprj_nearly_imaginary_cancel() {
    // `p + lambda` cancels to rounding noise in the first step (Arb)
    let rj = |e: f64| xsf::elliprj(c64(e, 1.0), c64(e, 4.0), c64(e, 9.0), c64(e, -11.0));
    let x = c64(1e-20, 1.0);
    let actual = [
        rj(1e-5),
        rj(1e-10),
        rj(1e-16),
        rj(1e-20),
        rj(1e-300),
        xsf::elliprj(c64(1e-20, -1.0), c64(1e-20, -4.0), c64(1e-20, -9.0), c64(1e-20, 11.0)),
        xsf::elliprj(ZERO, c64(1e-30, 4.0), c64(1e-30, 9.0), c64(1e-30, -6.0)),
        // in the second and third step
        xsf::elliprj(x, x, x, c64(1e-20, -15.0)),
        xsf::elliprj(x, x, x, c64(1e-20, -63.0)),
        // to the left of the imaginary axis
        xsf::elliprj(
            c64(4.474_347_940_814_568e-33, -24.022_291_057_907_644),
            c64(4.474_347_940_814_568e-33, -3.185_473_969_853_751),
            c64(4.474_347_940_814_568e-33, -24.0),
            c64(4.474_347_940_814_568e-33, 911.661_966_189_504_5),
        ),
        // in the eleventh step, where it is not negligible
        xsf::elliprj(
            c64(1.270_044_622_016_335_7e-213, 8.846_463_798_237_172),
            c64(5e-324, 0.002_158_110_463_198_337_3),
            c64(0.0, 0.001_077_764_749_758_734_6),
            c64(1.136_426_608_358_368_5e-30, -1_477_429.516_167_255_8),
        ),
    ];
    let expected = [
        c64(0.108_063_759_433_333_95, -0.003_007_992_608_451_945_6),
        c64(0.108_063_920_138_287_32, -0.003_008_153_312_457_716_8),
        c64(0.108_063_920_139_894_38, -0.003_008_153_314_064_776),
        c64(0.108_063_920_139_894_38, -0.003_008_153_314_064_777_6),
        c64(0.108_063_920_139_894_38, -0.003_008_153_314_064_777_6),
        c64(0.108_063_920_139_894_38, 0.003_008_153_314_064_777_6),
        c64(0.223_280_575_032_133_25, 0.001_136_428_124_214_924),
        c64(0.176_181_737_250_065_22, 0.072_051_668_386_978_52),
        c64(0.039_133_136_225_438_95, 0.026_116_877_617_553_11),
        c64(0.000_708_473_680_389_505_8, -0.000_472_991_720_576_317_77),
        c64(2.421_074_041_680_476e-6, 2.417_363_009_911_755_5e-6),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);

    // where `s beta` is subnormal (Arb)
    let a = 2.0_f64.powi(-28);
    let p = c64(f64::MIN_POSITIVE, -11.0 * a);
    let actual = xsf::elliprj(c64(0.0, a), c64(0.0, 4.0 * a), c64(0.0, 9.0 * a), p);
    let expected = c64(475_270_146_947.483_76, -13_229_998_187.788_53);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-12);

    // NaN if the lower bound `4^-m Re p` of `Re(p + lambda)` underflows (Arb: 2.26497 - 0.0723348i)
    let c = |im| c64(5e-324, im);
    let (y, z) = (c(0.548_967_078_054_124_9), c(0.172_925_512_765_675_9));
    assert!(xsf::elliprj(c(1.0), y, z, c(-1.464_874_055_611_607_8)).is_nan());

    // a conjugate pair with a tiny `Im p` (Arb)
    let x = c64(-1.0, 2.0_f64.powi(-24));
    let z = c64(2.0_f64.powi(64), 0.0);
    let p = c64(-256.999_999_999_999_94, -(2.0_f64.powi(-1020)));
    let actual = xsf::elliprj(x, x.conj(), z, p);
    let expected = c64(-5.484_067_219_012_265e-11, 4.285_892_810_353_007_5e-12);
    xsref::np_assert_allclose!(&[actual], &[expected], rtol = 1e-13);
}

#[test]
fn test_elliprj_conjugate_cut() {
    // a conjugate pair with `p` near `-lambda` and a tiny `Im p` (Arb)
    let x = c64(-0.269_822_709_532_446_47, 0.103_570_618_267_107_2);
    let (y, z) = (x.conj(), ZERO);
    let expected = c64(-29.989_358_762_276_39, 83.216_442_883_721_38);
    for pi in [-7.837_813_607_397_922e-21, -1e-100, -1e-300] {
        let p = c64(-0.289_017_590_377_695_1, pi);
        for [a, b, c] in [[x, y, z], [y, x, z], [z, x, y], [z, y, x], [x, z, y], [y, z, x]] {
            xsref::np_assert_allclose!(
                &[xsf::elliprj(a, b, c, p), xsf::elliprj(a, b, c, p.conj())],
                &[expected, expected.conj()],
                rtol = 1e-13
            );
        }
    }

    // and with `p` near `-lambda` of the steps 0, 1 and 2 (Arb)
    let x = c64(2.0_f64.powi(-200), 0.0);
    let y = c64(1.0, 2.0);
    for (pr, expected) in [
        (-2.236_067_977_499_79, c64(-0.762_778_196_568_029_6, 0.828_385_166_894_746_5)),
        (-11.651_457_212_829_742, c64(-0.254_346_363_442_386, 0.107_783_055_278_282_39)),
        (-48.817_624_918_411_6, c64(-0.067_612_796_497_439_43, 0.013_527_579_279_401_574)),
    ] {
        let p = c64(pr, -1e-100);
        for [a, b, c] in [[x, y, y.conj()], [y, y.conj(), x], [y.conj(), x, y]] {
            xsref::np_assert_allclose!(
                &[xsf::elliprj(a, b, c, p), xsf::elliprj(a, b, c, p.conj())],
                &[expected, expected.conj()],
                rtol = 2e-13
            );
        }
    }

    // an exact pair where rounding gives `lambda` a tiny imaginary part, near the cut and with a
    // tiny `Im p` (Arb)
    let x = c64(-2.583_169_022_335_466, 1.251_316_297_465_678_4e-70);
    let y = c64(8.112_302_455_109_168e-7, 0.0);
    let p = c64(-3.525_393_241_933_433, 2.880_059_528_591_252e-140);
    let actual = [xsf::elliprj(x, y, x.conj(), p)];
    let expected = [c64(-320.523_995_035_498_9, -2.663_686_039_152_461_5)];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);

    // and with a subnormal `Im p`, which `4^-m Im p` rounds to 0 in the steps (Arb)
    for (x, z, p, expected) in [
        (
            c64(391.945_997_773_030_4, 155.956_030_228_474_33),
            c64(7.968_794_498_639_412e33, 0.0),
            c64(-3.846_870_561_155_529e35, 2.293_954_062_989_63e-309),
            c64(-3.206_047_660_640_072_3e-51, -1.954_913_392_456_703_8e-53),
        ),
        (
            c64(1.0, 1.0),
            ZERO,
            c64(-1e32, -5e-324),
            c64(-4.121_030_609_847_326_3e-32, 4.712_388_980_384_689_5e-48),
        ),
        (
            c64(-1.848_896_916_592_894e24, 1.501_031_164_846_787_4e-53),
            ZERO,
            c64(-2.416_826_398_198_505e26, -4e-323),
            c64(-1.651_932_324_420_299_7e-36, 1.263_886_626_718_720_3e-39),
        ),
    ] {
        let actual = [xsf::elliprj(x, x.conj(), z, p)];
        xsref::np_assert_allclose!(&actual, &[expected], rtol = 1e-13);
    }

    // not a pair: conjugate only up to a subnormal imaginary part (Arb)
    let x = c64(-1.0, 5e-324);
    let actual = [xsf::elliprj(x, x, c64(0.25, 0.0), c64(-1.0, -4.0))];
    let expected = [c64(1.431_004_976_908_286_7, -0.043_252_977_798_958_4)];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_elliprj_real_cut() {
    // real `x, y, z` with `p` near `-lambda` of the steps 0, 1 and 2 and a tiny `Im p` (Arb)
    let at_3 = c64(-0.544_010_195_874_729_4, 0.589_048_622_548_086_2);
    for (pr, expected) in [
        (-3.0_f64.next_up(), at_3),
        (-3.0, at_3),
        (-3.0_f64.next_down(), at_3),
        (-15.0, c64(-0.175_527_524_442_984_58, 0.073_631_077_818_510_78)),
        (-63.0, c64(-0.046_138_727_260_895_78, 0.009_203_884_727_313_847)),
    ] {
        let p = c64(pr, -1e-20);
        xsref::np_assert_allclose!(
            &[xsf::elliprj(ONE, ONE, ONE, p), xsf::elliprj(ONE, ONE, ONE, p.conj())],
            &[expected, expected.conj()],
            rtol = 1e-13
        );
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "overflows a debug stack (#1765)")]
fn test_elliprj_real_cut_huge_p() {
    // real `x, y, z` with a huge `p`, whose distance to them overflows (Arb)
    let x = c64(2.0_f64.powi(-256), 0.0);
    let actual = [xsf::elliprj(x, x, x, c64(-f64::MAX, f64::MAX))];
    let expected = [c64(-2.839_325_246_800_259_7e-270, -2.839_325_246_800_259_7e-270)];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
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

#[test]
fn test_elliprg_wide_span() {
    let actual = [
        xsf::elliprg(1e-300, 1e150, 1e300),
        xsf::elliprg(1e-300, 1e300, 1e300),
        xsf::elliprg(2.23e-308, 6.76e-279, 1.03e-74),
        xsf::elliprg(c64(1.0, 1.0), c64(1e300, 0.0), c64(1e300, 1e-300)).re,
    ];
    // Arb
    let expected = [
        5e149,                       // 5.0000000000000001313e149
        7.853_981_633_974_483_5e149, // 7.8539816339744833023e149
        5.074_445_782_546_11e-38,    // 5.0744457825461098929e-38
        7.853_981_633_974_483_5e149, // 7.8539816339744833023e149
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // 1e-40 isn't treated as 0, since -1 is invalid
    let invalid = xsf::elliprg(c64(1e-40, 0.0), ONE, c64(-1.0, 0.0));
    assert!(invalid.re.is_nan());
}

#[test]
fn test_elliprj_wide_span() {
    let actual = [
        xsf::elliprj(1e-200, 1e-200, 1.0, 1.0),
        xsf::elliprj(1e-100, 1.0, 1e200, 1.0),
        // the direct computation of the "huge z" case is NaN
        xsf::elliprj(1e-320, 1e75, 1e294, 1e4),
        // `8 z / (a + g)` of the "tiny y" case overflows
        xsf::elliprj(1e-200, 1e-200, 1e300, 1.0),
    ];
    // Arb
    let expected = [
        689.854_969_439_893_6,        // 689.85496943989354116
        3e-100,                       // 3.0000000000000000454e-100
        1.490_188_239_869_415_3e-186, // 1.4901882398694151904e-186
        6.907_755_278_982_136e-148,   // 6.9077552789821368710e-148
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // a tiny argument prevents the downscaling of the huge ones
    let actual = [
        xsf::elliprj(1e-320, 1e180, 1e208, 1e180),
        xsf::elliprj(1e-200, 1e120, 1e120, 1e20),
        xsf::elliprj(1e-300, 1e150, 1e150, 1e150),
        xsf::elliprj(0.5, 1e156, 1e158, 2.5e-323),
    ];
    // Arb
    let expected = [
        3e-284,                     // 2.9999999999999999995e-284
        4.712_388_980_384_69e-130,  // 4.7123889803846899519e-130
        2.356_194_490_192_345e-225, // 2.3561944901923449966e-225
        1.577_252_122_957_69e-154,  // 1.5772521229576897772e-154
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_elliprj_wide_span_complex() {
    let x = c64(-1e70, 1e-100);
    let p2 = |k: i32| 2.0_f64.powi(k);
    let t = ONE * p2(-300);
    let w = c64(-p2(900), 1.0);
    let [x1, x2, x3] = [(300, -100), (300, -900), (-150, -460)].map(|(j, k)| c64(-p2(j), p2(k)));
    let x4 = c64(-1.75 * p2(-207), 1.125 * p2(-412));
    let x5 = c64(-p2(-28), p2(-519));
    let x6 = c64(-p2(24), 1.916_433_886_452_932_6e-261);
    let x7 = c64(-1.871_947_571_772_804e-21, 7.334_545_347_385_581e-291);
    let x8 = c64(-3.404_076_565_287_425e255, 3.569_679_440_543_941e-182);
    let x9 = c64(-p2(-60), 1.25 * p2(-848));
    let v = c64(-1e215, 1e-99);
    let s = f64::from_bits(1);
    let actual = [
        xsf::elliprj(c64(1e-60, 1e-60), c64(3e142, 3.5e142), c64(1e120, 3e120), c64(3e112, 9e112)),
        xsf::elliprj(
            c64(4e124, 2e125),
            c64(2e-125, -1e-125),
            c64(2.5e116, 2.6e116),
            c64(3e140, -4e140),
        ),
        // the duplicated `x` and `y` are equal
        xsf::elliprj(ZERO, ONE * 1e-130, ONE * 1e80, c64(-1e75, 1e74)),
        // the duplicated `p` has a tiny imaginary part, which mustn't prevent the downscaling
        xsf::elliprj(ONE * 1e-130, ONE * 1e150, ONE * 1e72, c64(-5e-82, -5e-81)),
        // tiny arguments next to a larger one, where `delta` underflows (#1655)
        xsf::elliprj(
            c64(4.111_978_644_178_181e-155, -5.905_609_635_520_994e72),
            c64(6.526_168_150_887_142e-140, 3.950_623_117_497_389_6e-196),
            c64(1.056_940_585_316_440_9e-284, -1.029_633_047_706_145_6e-256),
            c64(8.319_008_270_198_852e-254, -1.881_357_037_524_383_9e-261),
        ),
        xsf::elliprj(
            c64(7.880_448_859_174_146e-106, 2.019_480_815_864_376_2e-212),
            c64(2.416_079_673_010_958_7e-291, 2.343_972_336_601_560_3e-244),
            c64(1.003_671_132_383_083_1e-160, -1.199_478_819_847_264_8e-15),
            c64(1.504_115_561_441_218e-232, -8.063_479_214_635_844e-270),
        ),
        xsf::elliprj(
            c64(8.845_099_372_497_506e-180, -5.343_768_257_914_506e-286),
            c64(1.858_782_969_335_321_4e-173, 3.985_589_097_049_306e-6),
            c64(7.453_035_834_548_298e-260, -6.874_725_959_934_287e-181),
            c64(2.258_262_926_985_609e-167, 1.507_805_556_876_459e-293),
        ),
        // `p = x` near the negative real axis, where `p + lam` of the duplication step cancels
        xsf::elliprj(x, x.conj(), ONE * 1e-170, x),
        // also with a huge `sqrt(p)` over a tiny square root sum, a quotient that overflows for `v`
        xsf::elliprj(w, w.conj(), t, w),
        xsf::elliprj(v, v.conj(), ONE * 1e-296, v),
        // a slightly cancelling `p + lam` keeps the tiny `Im p`
        xsf::elliprj(c64(-3e-214, 7e-51), c64(-3e-214, -7e-51), ONE * 2e237, c64(-4e95, -2e-90)),
        xsf::elliprj(x1, x1.conj(), t, c64(-1.25 * p2(300), p2(-200))),
        // the `R_C` term of the duplication step overflows if divided by the tiniest factor first
        xsf::elliprj(x2, x2.conj(), t, x2),
        // `p = x`, where `delta = 0` at any scale
        xsf::elliprj(x3, x3.conj(), ONE * p2(-750), x3),
        // also with a subnormal `Im z`, which is in the domain as long as it isn't upscaled
        xsf::elliprj(x3, x3.conj(), c64(p2(-750), s), x3),
        // and with an `s^2` of the loop that underflows at any scale
        xsf::elliprj(x7, x7.conj(), ONE * 1.547_841_309_340_748_8e-15, x7),
        // `delta` underflows to 0 or is subnormal, but is normal once the arguments are upscaled
        xsf::elliprj(x3, x3.conj(), ONE * p2(-750), c64(x3.re, x3.im.next_up())),
        xsf::elliprj(x4, x4.conj(), ONE * p2(-850), c64(x4.re, x4.im / 8.0)),
        xsf::elliprj(x5, x5.conj(), ONE * p2(-628), c64(x5.re, -x5.im / 8.0)),
        // the product of two tiny factors of `delta` underflows before the huge one
        xsf::elliprj(
            x6,
            x6.conj(),
            ONE * 1.020_561_412_397_825_7e-44,
            c64(x6.re, 5.218_315_697_306_068e-261),
        ),
        // or is subnormal, also if the last factor makes `delta` normal
        xsf::elliprj(x9, x9.conj(), ZERO, c64(x9.re, 0.75 * x9.im)),
        // a subnormal square root sum of the duplication step
        xsf::elliprj(x8, x8.conj(), ONE * 2.726_537_368_598_458_6e-132, x8),
        // `delta` and `s^2` of the loop overflow, in the right half-plane
        xsf::elliprj(
            c64(1.0, -1.0),
            c64(1e-308, -3.0 * p2(358)),
            c64(1.0, -p2(359)),
            c64(1.0, p2(314)),
        ),
        // subnormal arguments, where the product of their square root sums would be subnormal
        xsf::elliprj(I * 1e65, ONE * s, ONE * 2.0 * s, ONE * 3.0 * s),
    ];
    // Arb
    let expected = [
        c64(-4.372_042_713_197_493e-189, -3.982_831_329_114_256e-188),
        c64(1.470_038_697_995_964_6e-202, 3.995_469_583_856_571_3e-203),
        c64(-7.026_491_305_778_003e-113, -7.482_776_732_461_581e-114),
        c64(4.460_686_995_705_570_6e-71, 4.929_003_648_729_797_4e-71),
        c64(1.803_132_905_404_062_6e160, 1.859_820_912_247_241e160),
        c64(2.794_562_084_032_514e176, 2.794_558_943_206_208e176),
        c64(7.033_083_592_633_553e170, -6.879_496_782_700_67e170),
        c64(-5.872_783_610_260_013e-103, -2.999_999_999_999_999_7e65),
        c64(0.0, -1.031_865_770_122_304_8e-135),
        c64(-3.43e-320, -9.486_832_980_505_139e-9),
        c64(-2.820_053_088_338_412_6e-212, 2.634_305_524_140_272_6e-214),
        c64(-1.142_299_594_953_818_3e-132, -5.798_932_689_929_481e-135),
        c64(1.234_608_741_757_442e-43, -1.776_715_956_459_856_7e226),
        c64(2.285_460_988_448_679e87, -3.374_185_345_907_397_3e161),
        c64(2.285_460_988_448_679e87, -3.374_185_345_907_397_3e161),
        c64(2.035_188_486_092_850_5e28, -1.039_643_822_894_742_4e298),
        c64(3.489_340_983_474_825e87, -3.374_185_345_907_397_3e161),
        c64(-2.741_580_245_462_780_4e95, -4.454_929_293_580_444e155),
        c64(-2.251_135_237_744_843e15, 1.228_958_694_522_29e161),
        c64(-1.345_154_158_350_593_6e-8, -2.503_882_255_550_643e257),
        c64(-1.014_312_557_716_746_1e30, -5.284_903_246_867_195e264),
        c64(4.413_282_569_619_264e-140, -1.440_430_726_332_907_9e54),
        c64(1.268_243_976_859_445_2e-155, 1.268_243_587_482_952_2e-155),
        c64(3.188_107_080_038_646_5e290, -3.188_107_080_038_646_5e290),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_elliprj_tiny_square_root_sum() {
    let x = c64(-8.492_010_046_825_438e99, 4.925_013_126_945_802e-253);
    let x1 = c64(-0.5, f64::from_bits(4));
    let x2 = c64(-2.500_350_457_769_907_7e58, 2.344_608_882_961_892_2e-300);
    let actual = [
        // `sqrt(p) / 2` over the tiny `sqrt(p) + sqrt(x)` of the duplication step overflows
        xsf::elliprj(x, x.conj(), ONE * 7.564_518_093_696_574e-224, c64(x.re, -x.im / 2.0)),
        // `sqrt(x) + sqrt(y)` is subnormal at any scale, in the duplicated `x` and `y` too
        xsf::elliprj(x1, x1.conj(), ONE * 2.0_f64.powi(800), x1),
        xsf::elliprj(x2, x2.conj(), ONE * 1.842_999_700_315_177_4e-290, x2),
    ];
    // Arb
    let expected = [
        c64(9.383_207_441_725_736e-76, 7.992_940_847_747_7e202),
        c64(-1.986_267_058_634_141_7e-96, -5.878_659_972_788_811e202),
        c64(0.058_606_435_074_860_28, -8.091_898_117_828_573e270),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_elliprj_subnormal_difference() {
    let p2 = |k: i32| 2.0_f64.powi(k);
    let s = f64::from_bits(1);
    let x = c64(-0.5, f64::from_bits(5));
    let actual = [
        // quartering the arguments rounds the subnormal `p - y` of a square root sum that cancels
        xsf::elliprj(x, x.conj(), ONE * p2(800), x),
        xsf::elliprj(x, x.conj(), ONE * p2(800), c64(x.re, f64::from_bits(3))),
        // also next to a huge `z`, as the arguments aren't scaled up
        xsf::elliprj(x, x.conj(), ONE * p2(1022), x),
        // a rounded `p - x` that nothing cancelling is formed from
        xsf::elliprj(ONE, ZERO, ONE * p2(1022), c64(1.0, s)),
        xsf::elliprj(ONE * p2(717), ZERO, ONE * p2(717), c64(p2(717), s)),
    ];
    // Arb
    let expected = [
        c64(-6.356_054_587_629_254e-97, -4.702_927_978_231_049e202),
        c64(9.294_220_595_638_364e-121, -5.451_253_281_024_553e202),
        c64(-4.896_526_343_390_018e-130, -1.811_501_963_914_339_8e169),
        c64(4.475_004_438_720_124e-154, 0.0),
        c64(s, 0.0),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // 0 by underflow, where the unquartered duplicated arguments would overflow
    let h = p2(1023);
    assert_eq!(xsf::elliprj(ONE * h, ZERO, ONE * h, c64(h, s)), ZERO);
    // and where the quartered ones do
    let x = c64(1.75 * h, 1.75 * h);
    assert_eq!(xsf::elliprj(x, x, x.conj(), ONE * s), ZERO);
}

#[test]
fn test_elliprj_subnormal_imag() {
    let p2 = |k: i32| 2.0_f64.powi(k);
    let (m, s) = (f64::MIN_POSITIVE, f64::from_bits(1));
    let x = c64(-p2(-30), 4.342_393_763_554_186_5e-272);
    let x1 = c64(-1.231_626_472_742_281_6e226, 5.530_266_466_854_912e-223);
    let x2 = c64(-1.892_873_431_416_007_7e206, 6.975_985_885_499_649e-178);
    let actual = [
        // quartered, the duplicated arguments get subnormal imaginary parts, which count as 0
        xsf::elliprj(x, x.conj(), ONE * 7.834_207_360_637_343e168, c64(x.re.next_down(), x.im)),
        xsf::elliprj(x1, x1.conj(), ZERO, c64(x1.re, -2.018_394_559_566_016e-223)),
        xsf::elliprj(x2, x2.conj(), ONE * 4.124_929_518_853_05e207, x2),
    ];
    // Arb
    let expected = [
        c64(-2.947_533_875_963_715_5e-57, -8.141_477_836_328_409e-60),
        c64(4.243_660_222_025_578_4e-75, 6.285_526_858_619_49e109),
        c64(-8.074_964_604_501_736e-165, -6.547_336_703_589_425e73),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // not if nothing cancels, where taking them as real is as accurate as before
    let actual = [xsf::elliprj(c64(m, -4.0 * s), ZERO, ONE * m, c64(p2(549), m))];
    let expected = [c64(1.714_357_124_141_202_7e-11, 3.806_637_503_303_462e-27)];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 2e-15);
}

#[test]
fn test_elliprj_subnormal_delta() {
    let p2 = |k: i32| 2.0_f64.powi(k);
    let x = c64(-1.167_889_111_935_278e45, 1.471_454_003_441_240_8e-172);
    let x1 = c64(-4.434_065_827_113_934e34, 2.186_661_079_437_059e-161);
    let x2 = c64(-3.014_284_785_580_670_7e-77, 3.022_044_004_073_942e-88);
    let y = c64(p2(-900), p2(-901));
    let actual = [
        // `delta` of the loop underflows, or a partial product of it, unless the arguments are
        // upscaled
        xsf::elliprj(
            x,
            x.conj(),
            ONE * 1.210_487_892_549_367e43,
            c64(x.re, -3.205_624_184_419_254e-172),
        ),
        xsf::elliprj(
            x1,
            x1.conj(),
            ONE * 8.214_861_731_877_567e-204,
            c64(x1.re, -5.662_168_100_451_648_5e-161),
        ),
        // `delta` is normal, but `delta 4^-2m` in the loop isn't
        xsf::elliprj(ZERO, y, ONE * p2(78), ONE * p2(-500)),
        // nor the product of its first two factors, or `s^2` (in the first one)
        xsf::elliprj(ZERO, y, ONE * p2(78), ONE * p2(-550)),
        xsf::elliprj(ZERO, y, ONE * p2(81), ONE * p2(-550)),
        // also if the arguments can't be upscaled
        xsf::elliprj(ZERO, y, ONE * p2(256), ONE * p2(-589)),
        // the step and the upscaling alternate
        xsf::elliprj(ZERO, ONE * f64::from_bits(1), ONE * p2(253), c64(p2(-1020), p2(-1020))),
        // not by the step if the subnormal `Im z` takes the shifted arguments out of the domain
        xsf::elliprj(
            x2,
            x2.conj(),
            c64(1.763_600_756_787_408_2e-249, 1e-310),
            c64(1.5 * x2.re, x2.im),
        ),
    ];
    // Arb (the last one with `Im z = 0`)
    let expected = [
        c64(-1.845_574_348_790_100_6e-65, 4.337_305_767_725_274e149),
        c64(-7.214_117_479_787_963e-50, 4.378_123_669_686_43e143),
        c64(2.487_693_345_481_894e141, -4.141_019_586_341_121_4e138),
        c64(2.452_384_296_026_243e156, -4.662_373_566_494_95e153),
        c64(8.670_487_828_977_77e155, -1.648_397_982_646_743_8e153),
        c64(1.936_714_749_728_869_7e141, -4.141_019_586_341_121_4e138),
        c64(2.798_077_604_098_084e270, -2.688_054_912_283_786_6e270),
        c64(-9.257_333_890_993_706e116, -4.649_958_488_504_821e115),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_elliprj_tiny_y() {
    // `y <= 1e-26` but not small enough for the "tiny y" asymptotic case
    let actual = [xsf::elliprj(0.0, 1e-21, 6e-29, 3e-9), xsf::elliprj(5e-37, 1e-30, 1e-19, 2e-9)];
    // Arb
    let expected = [
        3.067_647_836_331_202e20, // 3.0676478363312025864e20
        6.664_402_648_223_819e19, // 6.6644026482238192011e19
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_elliprj_huge_p() {
    // complex `x, y, z << p`, where the `s` of the main loop cancels
    let actual = [
        xsf::elliprj(
            c64(2.4e-65, 1.8e-64),
            c64(1e-142, -5.2e-142),
            c64(4.2e-148, -1.35e-147),
            c64(4e16, -2e15),
        ),
        xsf::elliprj(
            c64(4.16e-18, 9.22e-18),
            c64(4.2e-113, -2.86e-112),
            c64(3.6e-32, 2.6e-33),
            c64(1.32e56, -4e56),
        ),
        xsf::elliprj(
            c64(1.5e-130, 2.7e-130),
            c64(2.4e-150, -8.6e-150),
            c64(3.4e-133, -3.4e-132),
            c64(2.2e-54, 2.5e-55),
        ),
    ];
    // Arb
    let expected = [
        // 400084498862560409.57 - 306435843021533679.75i
        c64(4.000_844_988_625_604e17, -3.064_358_430_215_337e17),
        // 3.0658237668252603484e-47 + 2.6265878247629335519e-47i
        c64(3.065_823_766_825_260_6e-47, 2.626_587_824_762_933_7e-47),
        // 2.8232070639983970999e119 - 9.0800781779622300614e118i
        c64(2.823_207_063_998_397e119, -9.080_078_177_962_23e118),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // the asymptotic form for `z / p <= 5e-14`
    let actual = [xsf::elliprj(1.0, 2.0, 3.0, 1e14), xsf::elliprj(0.0, 0.5, 7.0, 2e14)];
    // Arb
    let expected = [
        2.180_837_335_167_910_8e-14, // 2.1808373351679106680e-14
        1.551_984_603_581_836e-14,   // 1.5519846035818361032e-14
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_carlson_near_overflow_and_tiny() {
    let max = f64::MAX;
    let actual = [
        xsf::elliprd(1e308, 1e308, 1e-310),
        xsf::elliprd(1e300, max, 1e-320),
        xsf::elliprf(max, max, 2.3e-308),
        xsf::elliprc(1e-310, 1e308),
        xsf::elliprc(5e-324, max),
    ];
    // Arb
    let expected = [
        3.000_000_000_000_005e-153,   // 3.0000000000000045497e-153
        2.237_514_674_315_783_5e-144, // 2.2375146743157834977e-144
        1.171_553_422_455_405e-154,   // 1.1715534224554048805e-154
        1.570_796_326_794_896_6e-154, // 1.5707963267948966106e-154
        1.171_553_422_455_405e-154,   // 1.1715534224554048805e-154
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // Arb: 1.0649531960466087604e-459
    assert_eq!(xsf::elliprd(1e-310, 1.0, 1e308), 0.0);
}

#[test]
fn test_elliprd_near_underflow() {
    // Arb: 1.4583051862293731565e-308
    let actual = xsf::elliprd(6.5e-314, 2e211, 4.6e202);
    xsref::np_assert_allclose!(&[actual], &[1.458_305_186_229_373_3e-308], rtol = 1e-15);

    // subnormal results
    let actual = [xsf::elliprd(2.1e199, 3.5e-316, 4.4e206), xsf::elliprd(3.4e240, 1.5e-323, 1e190)];
    // Arb
    let expected = [
        2.865_320_665_163_82e-309, // 2.8653206651638199877e-309
        1.626_978_433_639_93e-310, // 1.6269784336399211517e-310
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-13);
}

#[test]
fn test_elliprd_overflow() {
    // just below the overflow threshold, where `3 * adt` overflows
    let actual = xsf::elliprd(0.0, 1.965_996_116_174_601_3e-8, 1.190_184_710_917_690_7e-304);
    xsref::np_assert_allclose!(&[actual], &[1.797_693_134_862_315_3e308], rtol = 1e-15);

    // the first term of the sum overflows, and the sum itself (2e308), respectively
    assert_eq!(xsf::elliprd(1e-290, 1e-70, 1e-290), INF);
    assert_eq!(xsf::elliprd(0.0, 1e-16, 1.5e-300), INF);

    // mpmath: 1.5e325 - 3.75e314 i
    let actual = xsf::elliprd(c64(1e-290, 1e-300), c64(1e-70, 0.0), c64(1e-290, 0.0));
    assert_eq!(actual, c64(INF, -INF));
}

#[test]
fn test_elliprj_cauchy_tiny() {
    // `x y` of DLMF 19.20.14 underflows, and `p q` as well in the first case
    let actual = [
        xsf::elliprj(1e-200, 1e-200, 2.0, -1e-210),
        xsf::elliprj(2e-89, 1e-289, 0.06, -1.5e-218),
        xsf::elliprj(1e-200, 1e-200, 2.0, -1.0),
        // `sqrt(z / (pn q))` overflows
        xsf::elliprj(1.0, 1.0, 1e300, -5e-324),
        xsf::elliprj(0.0, 1.0, 1e300, -5e-324),
    ];
    // Arb
    let expected = [
        2.442_260_300_028_494_3e201,  // 2.4422603000284942741e201
        5.773_502_691_896_257_5e118,  // 5.7735026918962579769e118
        -488.672_336_815_068_3,       // -488.67233681506833085
        1.116_660_107_882_071_9e-147, // 1.116660107882071864e-147
        -3e-150,                      // -2.999999999999999921e-150
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // NaN (or else correct), as the terms of the transform cancel; Arb: -122.39924517639211
    let actual = xsf::elliprj(1e-200, 1e-200, 2.0, -1e-200);
    assert!(actual.is_nan() || (actual + 122.399_245_176_392_11).abs() < 1e-13);
}

#[test]
fn test_elliprj_tiny_p() {
    let p = |k| 2.0_f64.powi(k);
    let actual = [
        xsf::elliprj(1e-276, 1e-48, 1e-276, 1e-276),
        xsf::elliprj(1.0, 1.0, 1.0, 1e-310),
        // `x * y` underflows
        xsf::elliprj(1e-200, 1e-120, 1e100, 1e-220),
        // `sqrt(x * y * z)` overflows
        xsf::elliprj(p(683), p(683), p(683), p(-300)),
        // `p / (p + lambda)` underflows
        xsf::elliprj(0.0, p(-1020), p(1022), 5e-324),
    ];
    let expected = [
        1.499_999_999_999_999_9e300,
        1_069.781_509_783_911,
        7.115_699_433_150_121e111,
        4.016_505_204_038_232_4e-306,
        1.060_031_584_427_373_4e162,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // complex input with `p << x, y, z`
    let actual = [
        xsf::elliprj(c64(9.37, 0.0), c64(2.45, 0.0), c64(1.56, 0.0), c64(8.8e-201, 4.8e-201)),
        xsf::elliprj(c64(3.0, 4.0), c64(2.0, -1.0), c64(0.5, 7.0), c64(8.8e-151, 4.8e-151)),
        xsf::elliprj(c64(4.0, 3.0), c64(2.0, -5.0), c64(6.0, 1.0), c64(1e-280, 0.0)),
        // `Re p < 0` with conjugate `y` and `z`
        xsf::elliprj(c64(6.5, 0.0), c64(4.3, 2.7), c64(4.3, -2.7), c64(-3e-194, 1.7e-194)),
    ];
    // Arb
    let expected = [
        // 115.53541545358995929 - 0.12516368047983641353i
        c64(115.535_415_453_589_96, -0.125_163_680_479_836_42),
        // 32.582194210797886059 - 48.772720866609911523i
        c64(32.582_194_210_797_89, -48.772_720_866_609_91),
        // 74.344070551451420321 + 14.347796406788596231i
        c64(74.344_070_551_451_42, 14.347_796_406_788_596),
        // 51.752360261895413592 - 0.30429517368485992370i
        c64(51.752_360_261_895_41, -0.304_295_173_684_859_9),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    // the duplication step rounds `Re x` below 0, or `Re p` to 0
    let actual = [
        xsf::elliprj(c64(1e-20, 1.0), c64(1e-20, 2.0), c64(1e-20, 3.0), c64(1e-20, 0.0)),
        xsf::elliprj(c64(0.0, 1.0), c64(0.0, 2.0), c64(0.0, 3.0), c64(5e-324, 1e-14)),
    ];
    // Arb
    let expected = [
        c64(-19.226_711_474_154_413, -20.587_060_997_330_077),
        c64(-13.924_594_681_017_782, -13.924_594_681_017_782),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_elliprj_subnormal() {
    let actual = [
        // R_D(0, 5e-324, 1)
        xsf::elliprj(0.0, 5e-324, 1.0, 1.0),
        // `z / p` underflows
        xsf::elliprj(0.0, 5e-324, 5e-324, 3.0),
        // `(x + y) / 2 + sqrt(x y)` is subnormal, and 0, respectively
        xsf::elliprj(5e-324, 1e-320, 1e-10, 1.0),
        xsf::elliprj(0.0, 5e-324, 1e20, 0.5),
    ];
    // Arb
    let expected = [
        1_117.818_990_965_431_6,    // 1117.8189909654315653
        7.066_877_263_035_343e161,  // 7.0668772630353430919e161
        107_479_496.807_218_07,     // 107479496.80721808114
        2.235_399_657_305_823_6e-7, // 2.2353996573058236229e-7
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    assert_eq!(xsf::elliprj(0.0, 0.0, 1.0, 5e-324), INF);
    // Arb: 8.2588819190784932822e308
    assert_eq!(xsf::elliprj(0.0, 5e-324, 0.5, 1e-307), INF);
}

#[test]
fn test_elliprj_negative_delta() {
    // `(p - x) (p - y) (p - z) < 0`, e.g. for `x, y << p << z`
    let actual = [
        xsf::elliprj(3.8e-231, 3.5e-20, 3.5, 3.7e-8),
        xsf::elliprj(1e-200, 1e-20, 1.0, 1e-8),
        xsf::elliprj(7e5, 1.4e5, 4.7e8, 4e-8),
    ];
    // Arb
    let expected = [
        630_004_492.889_770_3,    // 630004492.88977032757
        4_352_597_307.454_136,    // 4352597307.4541356094
        6.526_307_207_698_169e-9, // 6.5263072076981690311e-9
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

    let actual = [
        xsf::elliprj(c64(1.8e11, -3.5e11), c64(6e8, 3e7), c64(7e8, -8e8), c64(9e-12, 5e-12)),
        xsf::elliprj(c64(2.7e6, -1.6e6), c64(1.4e5, -8e4), c64(1.2e5, -2.2e5), c64(1e-11, 1.7e-11)),
        // `d + r 4^-m` cancels in a step, after which `w2` is inf
        xsf::elliprj(ONE * 1e-100, ONE * 1e-6, ONE * 1e30, c64(-3e102, 1e95)),
    ];
    // Arb
    let expected = [
        c64(8.193_228_841_310_918e-14, 1.099_562_954_582_950_9e-13),
        c64(8.189_511_559_078_331e-8, 1.316_376_299_672_930_8e-7),
        c64(-4.283_282_603_501_267e-116, -1.427_760_867_833_755_7e-123),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_carlson_near_conjugate() {
    // a conjugate pair close to the negative real axis
    let x = c64(-1e-4, 1e-30);
    let actual = [
        xsf::elliprf(x, c64(0.5, 0.0), x.conj()),
        xsf::elliprd(x, x.conj(), c64(0.5, 0.0)),
        // the real part is a tiny remainder of the cancelling terms, so it's inaccurate
        xsf::elliprd(x, c64(0.5, 0.0), x.conj()),
        xsf::elliprj(x, c64(0.5, 0.0), x.conj(), c64(0.3, 0.0)),
        // R_J(x, y, z, z) = R_D(x, y, z)
        xsf::elliprj(x, c64(0.5, 0.0), x.conj(), x.conj()),
        // nearly conjugate
        xsf::elliprf(x, c64(0.5, 0.0), c64(-1.000_001e-4, -1e-30)),
        // `x - conj(x)` overflows
        xsf::elliprf(c64(-6e307, 1e308), c64(-6e307, -1e308), c64(1.0, 0.0)),
    ];
    let expected = [
        c64(92.639_053_758_880_68, 0.0),
        c64(547.241_289_722_859_8, 0.0),
        c64(20_939.582_790_734_995, 4.242_216_486_679_578_4e30),
        c64(909.422_264_352_536_2, 0.0),
        c64(20_939.582_790_734_995, 4.242_216_486_679_578_4e30),
        c64(28.498_784_620_315_163, 2.221_219_358_137_306),
        c64(2.009_168_978_782_164_5e-154, 0.0),
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}
