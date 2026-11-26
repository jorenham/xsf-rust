use std::f64::consts::FRAC_PI_4;

use num_complex::c64;
use xsref::np_assert_allclose;

use crate::DENORM_MIN;

xsref::tests! {
    test_hyp2f1_f64: "hyp2f1", "d_d_d_d-d", |x| xsf::hyp2f1(x[0], x[1], x[2], x[3]);
    test_hyp2f1_c64: "hyp2f1", "d_d_d_cd-cd", |x| xsf::hyp2f1(x[0], x[1], x[2], c64(x[3], x[4]));
}

#[test]
fn test_hyp2f1_symmetric_real() {
    // Arb at the two ends of the transformed parameter range.
    for (h, expected) in [(0.5, 1.372_880_500_618_350_2), (32.0, 105_012_102.535_169_36)] {
        let actual = xsf::hyp2f1(h, h, 2.0 * h, 0.75);
        np_assert_allclose!([actual], [expected], rtol = 8e-15);
        for y in [0.0, -0.0] {
            let complex = xsf::hyp2f1(h, h, 2.0 * h, c64(0.75, y));
            assert_eq!(complex.re, actual);
            assert_eq!(complex.im.to_bits(), 0.0_f64.to_bits());
        }
    }
}

#[test]
fn test_hyp2f1_c64_unit_exponent_reciprocal() {
    for (a, b, c) in [
        (1.0, 1.0, 1.0),
        (2f64.powi(60), 1.0, 2f64.powi(60)),
        (1.0, 2f64.powi(60), 2f64.powi(60)),
    ] {
        for sign in [-1.0, 1.0] {
            let actual = xsf::hyp2f1(a, b, c, c64(10.0, sign * 10.0));
            np_assert_allclose!(
                [actual.re, actual.im],
                [-9.0 / 181.0, sign * 10.0 / 181.0],
                rtol = 4e-15
            );
            let actual = xsf::hyp2f1(a, b, c, c64(1.0, sign * f64::MAX));
            assert_eq!(actual.re.to_bits(), 0.0_f64.to_bits());
            assert_eq!(actual.im, sign / f64::MAX);
        }
    }
    let below = f64::from_bits(1.0_f64.to_bits() - 1);
    let above = f64::from_bits(1.0_f64.to_bits() + 1);
    for x in [below, above] {
        let actual = xsf::hyp2f1(1.0, 1.0, 1.0, c64(x, -1.0));
        assert_eq!(actual, c64(1.0 - x, -1.0));
    }
}

#[test]
fn test_hyp2f1_c64_common_parameter_polynomial() {
    for (a, b) in [(500.0, -500.0), (-500.0, 500.0)] {
        for y in [0.0, -0.0] {
            let actual = xsf::hyp2f1(a, b, 500.0, c64(0.75, y));
            assert_eq!(actual.re, 2.0_f64.powi(-1000));
            assert_eq!(actual.im.to_bits(), 0);
        }
    }
    // The rounded subtraction would give 1 before taking the power.
    let actual = xsf::hyp2f1(500.0, -1_048_575.0, 500.0, c64(2.0_f64.powi(-55), 0.0));
    assert_eq!(actual, c64(0.999_999_999_970_896_2, 0.0));
    // The real-axis shortcut must preserve small components off the axis.
    let actual = xsf::hyp2f1(2.0, -2.0, 2.0, c64(2.0_f64.powi(-54), 1.0));
    assert_eq!(actual, c64(-2.0_f64.powi(-53), -2.0));
    let y = DENORM_MIN;
    let actual = xsf::hyp2f1(4.0, -3.0, 4.0, c64(1.5, y));
    assert_eq!(actual, c64(-0.125, -y));
    // With c negative, retain termination at c rather than canceling it.
    assert_eq!(xsf::hyp2f1(-2.0, 3.0, -2.0, c64(0.5, 0.0)), c64(4.0, 0.0));
}

#[test]
fn test_hyp2f1_c64_linear_subnormal_denominator() {
    // For a=b=-1 the exact polynomial is 1 + z/c.
    let c = -f64::from_bits(1_u64 << 50);
    let large = 2.0_f64.powi(1022);
    for zero in [0.0, -0.0] {
        let value = xsf::hyp2f1(-1.0, -1.0, c, c64(0.25, zero));
        assert_eq!(value.re, -large);
        assert_eq!(value.im.to_bits(), 0);
    }
    for sign in [1.0, -1.0] {
        let value = xsf::hyp2f1(-1.0, -1.0, c * 0.5, c64(0.125, sign * 0.125));
        assert_eq!(value, c64(-large, -sign * large));
    }
    let value = xsf::hyp2f1(-1.0, -1.0, c, c64(-c, DENORM_MIN));
    assert_eq!(value, c64(0.0, -2.0_f64.powi(-50)));
    assert!(xsf::hyp2f1(-1.0, -1.0, c, c64(1.0, 0.0)).re.is_infinite());
    assert!(
        xsf::hyp2f1(-1.0, -1.0, 0.0, c64(0.25, 0.0))
            .re
            .is_infinite()
    );
}

#[test]
fn test_hyp2f1_c64_pfaff_accuracy() {
    // Arb at the exact binary64 inputs.
    for sign in [-1.0, 1.0] {
        let actual = xsf::hyp2f1(
            8.095_813_935_368_371,
            4.001_376_844_959_068_5,
            4.078_873_014_294_075,
            c64(-0.947_368_421_052_631_6, sign * 0.526_315_789_473_684_1),
        );
        let expected = c64(-0.001_809_357_394_137_878_3, sign * 0.003_481_887_377_423_739);
        assert!((actual - expected).norm() <= 4.0 * f64::EPSILON * expected.norm());
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_symmetric_beta_weight() {
    for (a, b, c, z, expected) in [
        (
            2.25,
            8.25,
            16.5,
            c64(0.113_793_103_448_275_98, 0.948_275_862_068_965_7),
            c64(0.446_860_075_021_192_6, 0.731_321_493_403_688_5),
        ),
        (
            16.25,
            2.25,
            4.5,
            c64(0.113_793_103_448_275_98, -1.024_137_931_034_483),
            c64(-0.021_965_227_124_574_663, 0.009_908_300_237_809_064),
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            let actual = xsf::hyp2f1(a, b, c, z);
            np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 8e-15);
        }
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_small_odd_terms() {
    let b = 16.3125_f64.to_bits();
    let z = c64(0.5, 0.866_025_403_784_438_6);
    // Arb at the exact binary64 inputs.
    for (b, expected) in [
        (b - 1, c64(-1.145_189_436_074_619, 2.764_731_868_057_758_6)),
        (b, c64(-1.145_189_436_074_624_2, 2.764_731_868_057_756_6)),
        (b + 1, c64(-1.145_189_436_074_629_2, 2.764_731_868_057_754_5)),
    ] {
        let actual = xsf::hyp2f1(15.75, f64::from_bits(b), 32.625, z);
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 1e-14);
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_parameter_order() {
    // Arb at the exact binary64 inputs, also at conj(z).
    for (a, b, c, z, expected) in [
        (
            9.854,
            -14.732,
            14.16,
            c64(0.638, -0.8683),
            c64(-0.001_838_354_430_051_382_3, 0.014_966_128_267_305_905),
        ),
        (
            8.849,
            -14.946,
            -1.694,
            c64(0.5993, -0.8826),
            c64(527_042_171.148_339_2, -558_878_784.684_578_3),
        ),
        // The odd terms vanish in one order, where c = 2 (c - a), or are small one ulp away.
        (
            -7.75,
            4.25,
            -15.5,
            c64(0.341_379_310_344_827_7, -0.948_275_862_068_965_5),
            c64(-1.050_983_485_011_692_1, -1.114_552_232_548_607_5),
        ),
        (
            -7.749_999_999_999_999,
            4.25,
            -15.5,
            c64(0.341_379_310_344_827_7, -0.948_275_862_068_965_5),
            c64(-1.050_983_485_011_692, -1.114_552_232_548_608),
        ),
        (
            -7.750_000_000_000_001,
            4.25,
            -15.5,
            c64(0.341_379_310_344_827_7, -0.948_275_862_068_965_5),
            c64(-1.050_983_485_011_692_6, -1.114_552_232_548_607_2),
        ),
        // c = 2b, with a term ratio dividing by c/2 + 7/2 = 2^-51 after a small term.
        (
            -3.999_999_999_999_999_6,
            -3.499_999_999_999_999_6,
            -6.999_999_999_999_999,
            c64(0.5, 0.866_025_403_784_438_6),
            c64(-0.378_896_182_077_845_76, -0.656_267_438_152_697_5),
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            for actual in [xsf::hyp2f1(a, b, c, z), xsf::hyp2f1(b, a, c, z)] {
                np_assert_allclose!(
                    [actual.re, actual.im],
                    [expected.re, expected.im],
                    rtol = 1e-13
                );
            }
        }
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_dip() {
    // Arb at the exact binary64 inputs, also at conj(z).
    for (a, b, c, z, expected) in [
        // The terms dip, and grow again past c + n = -1/4.
        (
            0.5,
            -95.125,
            -105.25,
            c64(0.375, 0.906_25),
            c64(0.879_448_710_243_232_6, 0.420_597_622_480_321),
        ),
        // c + n stays negative within the term limit.
        (
            0.5,
            -3000.5,
            -3000.25,
            c64(0.0, 0.9375),
            c64(0.794_258_086_626_283_2, 0.314_106_447_021_066),
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            let actual = xsf::hyp2f1(a, b, c, z);
            np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 1e-13);
        }
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_growing_terms() {
    // The terms overflow, and those with the upper parameters swapped start tiny but grow.
    let actual = xsf::hyp2f1(-2250.5, 1e-25, 50.25, c64(0.0, 1.0));
    // Arb at the exact binary64 inputs, or NaN.
    if !actual.is_nan() {
        let expected = c64(6.450_899_072_505_266_4e215, -1.134_025_706_455_491_9e215);
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 1e-13);
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_near_c() {
    // Arb at the exact binary64 inputs, also at conj(z).
    // With c near b, the terms grow far past the sum and almost cancel again.
    let z = c64(0.5625, -0.8125);
    for (c, expected) in [
        (-56.874_999, c64(-0.039_422_469_447_030_41, 0.885_589_629_206_983_8)),
        (-56.874_999_9, c64(-0.039_422_448_799_181_074, 0.885_589_619_897_664_7)),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            for actual in [xsf::hyp2f1(-1.5, -56.875, c, z), xsf::hyp2f1(-56.875, -1.5, c, z)] {
                np_assert_allclose!(
                    [actual.re, actual.im],
                    [expected.re, expected.im],
                    rtol = 1e-10,
                    atol = 1e-10
                );
            }
        }
    }
}

#[test]
fn test_hyp2f1_c64_lopez_temme_paired_sum_fallback() {
    let b = f64::from_bits(5e11_f64.to_bits() - 1);
    let actual = xsf::hyp2f1(1000.0, b, 1e12, c64(0.5, 0.866_025_403_784_438_6));
    // Lopez-Temme moments evaluated in Arb at the exact binary64 inputs.
    let expected = c64(-1.473_459_611_433_148e62, 2.552_106_909_903_471e62);
    np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 2e-13);
}

#[test]
fn test_hyp2f1_c64_transform2_euler() {
    // unlike in xsf: the 1/z transformation of the Euler transform (mpmath references)
    for (a, b, c, z, expected) in [
        // a - b an integer: the limiting case
        (8.5, 6.5, -4.25, c64(-1.0, 0.5), c64(-14.914_154_962_801_387, 2.508_617_290_210_916_7)),
        (4.3, 4.4, -4.2, c64(-0.9, 0.8), c64(-13.468_832_852_470_712, -27.983_525_735_799_486)),
        // the limiting case with b within 1e-15 of an integer
        (
            8.881_784_197_001_252e-16,
            8.881_784_197_001_252e-16,
            -1.25,
            c64(-1.0, 0.5),
            c64(1.0, -4.661_160_101_737_811e-31),
        ),
        // (1 - z)^(c - a - b) is subnormal
        (
            9.48,
            9.48,
            -7.69,
            c64(2.97e11, 8.92e11),
            c64(3.884_057_170_000_241_6e-108, -6.451_235_141_304_403e-108),
        ),
        // c - (c - a) loses a
        (
            5e-6,
            -3.125,
            -9.0625,
            c64(133.5, -390.75),
            c64(3.466_873_805_327_974, -0.823_490_927_100_644_6),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        assert!(
            (actual - expected).norm() < 1e-13 * expected.norm(),
            "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
        );
    }
}

#[test]
fn test_hyp2f1_c64_inverse_lopez_temme() {
    // The 1/z terms cancel, but the Lopez-Temme series does not. Arb at the binary64 inputs.
    for (a, b, c, z, expected) in [
        (
            15.75,
            16.5,
            32.625,
            c64(0.75, 1.0),
            c64(-1.587_295_690_150_468_2, -4.321_281_223_476_627),
        ),
        (
            27.048,
            29.901,
            46.704,
            c64(-1.483_752_897_983_188, 0.369_995_905_609_408_26),
            c64(-3.999_435_224_208_415e-8, 9.445_129_048_405_783e-9),
        ),
        (
            4.551,
            9.229,
            40.233,
            c64(0.487_808_055_590_822_96, -1.306_392_206_115_161_6),
            c64(0.109_329_995_590_975_23, -1.193_178_024_137_373_5),
        ),
        // With a < 0 the paired Lopez-Temme sum does not apply, so keep the 1/z sum.
        (
            -10.5,
            -32.75,
            -60.25,
            c64(0.0, 12.0),
            c64(-74_172_362.262_856_14, 60_548_017.841_367_66),
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            let actual = xsf::hyp2f1(a, b, c, z);
            assert!(
                (actual - expected).norm() <= 2e-14 * expected.norm(),
                "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
            );
        }
    }
}

#[test]
fn test_hyp2f1_c64_inverse_tiny_imag() {
    // The imaginary part of the Lopez-Temme prefactor underflows, so keep the 1/z sum.
    // Arb at the binary64 inputs, compared componentwise.
    let z = c64(-2.0, 2.0_f64.powi(-970));
    let expected = c64(1.415_783_188_282_483_4e-5, 4.484_815_861_380_78e-297);
    for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
        let actual = xsf::hyp2f1(128.25, 8.75, 100.125, z);
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 1e-13);
    }
}

#[test]
fn test_hyp2f1_c64_inverse_near_equal() {
    // a - b is near an integer, so the poles of the 1/z terms nearly cancel (the last two after
    // the Euler transform). Arb at the binary64 inputs, compared componentwise.
    for (a, b, c, z, expected) in [
        (
            1.0 + 2.0_f64.powi(-47),
            1.0 + f64::EPSILON,
            2.0,
            c64(2.0, 0.25),
            c64(0.163_330_184_812_119_58, 1.427_890_722_129_949_3),
        ),
        (
            0.618_776_231_689_070_1,
            0.618_776_231_689_068_9,
            2.982_633_969_189_627_5,
            c64(-3.459_733_432_750_732, -0.853_452_349_090_242_8),
            c64(0.764_767_445_893_113_2, -0.033_055_478_840_512_145),
        ),
        (
            0.997_388_517_972_448_7,
            -1.002_611_482_027_599_8,
            3.150_968_818_948_326_7,
            c64(2.045_922_291_999_997_5, 0.333_194_302_263_917),
            c64(0.351_688_010_827_500_9, -0.104_787_053_409_834_48),
        ),
        (
            -6.198_721_292_721_535,
            -3.198_721_292_721_532_3,
            -2.126_916_436_678_337,
            c64(-1.018_479_056_347_564_7, 1.446_640_991_662_744_5),
            c64(-2_615.177_010_022_977_4, -2_193.034_090_860_991),
        ),
        (
            2.891_850_628_025_062,
            2.891_850_628_025_058,
            -1.673_453_207_224_207,
            c64(-0.656_577_054_389_945_3, 0.973_480_163_347_360_1),
            c64(-0.486_955_198_463_085_9, -6.854_834_441_578_813),
        ),
        (
            3.0,
            2.000_000_000_000_001_3,
            4.0,
            c64(3.514_055_857_750_915, -0.331_333_151_666_662_44),
            c64(0.154_015_399_276_666_3, 0.445_105_599_187_927_8),
        ),
    ] {
        for (a, b) in [(a, b), (b, a)] {
            for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
                let actual = xsf::hyp2f1(a, b, c, z);
                // (2e-14, as the real part of the fifth is 1/14 of the norm)
                np_assert_allclose!(
                    [actual.re, actual.im],
                    [expected.re, expected.im],
                    rtol = 2e-14
                );
            }
        }
    }
}

#[test]
fn test_hyp2f1_c64_inverse_near_equal_edges() {
    // a - b is near an integer: with large c, with c - b an integer, near 1/8, or with an
    // overflowing (1 + d)_m. Arb at the binary64 inputs, compared componentwise.
    for (a, b, c, z, expected) in [
        (2.125, 1.0, 101.0, c64(100.0, 1.0), c64(-0.382_351_731_216_135_9, 1.089_877_562_791_618)),
        (
            28.125_000_000_000_227,
            3.125,
            6.125,
            c64(2.0, 0.25),
            c64(-5.932_929_836_608_319e-4, -3.525_838_182_852_035e-5),
        ),
        (
            1.125,
            1.0,
            2.0625,
            c64(2.0, 0.25),
            c64(-0.027_024_445_301_198_41, 1.453_291_172_228_273_7),
        ),
        (
            3_000.125 + 2.0_f64.powi(-20),
            0.125,
            1.25,
            c64(100.0, 1.0),
            c64(0.183_905_066_744_820_52, 0.075_906_798_189_663_85),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        // (the second and third lose about 3e-14 in one component to cancellation)
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 5e-14);
    }
    // On the real axis below -1, the imaginary part stays +0.
    for z in [c64(-1.125, 0.0), c64(-1.125, -0.0)] {
        let actual = xsf::hyp2f1(-2.875 + 2.0_f64.powi(-14), -2.875, 1.5, z);
        np_assert_allclose!([actual.re], [-0.697_502_355_699_381_9], rtol = 4e-15);
        assert_eq!(actual.im.to_bits(), 0.0_f64.to_bits());
    }
}

#[test]
fn test_hyp2f1_c64_inverse_integer_large_a() {
    // a - b is an integer and Gamma(a) overflows (gh-2244); in the last two, (1 - z)^(c - a - b)
    // underflows. Arb at the binary64 inputs, compared componentwise.
    for (a, b, c, z, expected) in [
        (
            180.0,
            100.0,
            125.5,
            c64(-2.0, 0.5),
            c64(3.020_509_328_819_671_4e-66, 1.708_047_496_988_71e-65),
        ),
        (
            187.0,
            105.0,
            125.5,
            c64(-2.0, 0.5),
            c64(-2.482_296_999_002_554_6e-71, 2.983_532_921_504_48e-72),
        ),
        (
            187.0,
            90.0,
            125.5,
            c64(-2.0, 0.5),
            c64(7.206_274_521_205_293e-60, -6.691_889_189_747_599e-60),
        ),
        (
            197.46,
            -0.539_999_999_999_992,
            24.096,
            c64(4.319_960_734_025_822, 9.391_402_568_842_924),
            c64(5.236_383_181_499_285, -9.689_051_624_988_759),
        ),
        (
            187.0,
            105.0,
            124.999_999_661_270_98,
            c64(-22.0, 0.999_999_999_767_169_4),
            c64(1.338_564_001_322_857_7e-176, -2.992_835_631_306_290_3e-174),
        ),
        (
            231.0,
            80.0,
            6.5,
            c64(-14.817_493_273_898_183, 23.847_506_713_102_423),
            c64(5.414_467_502_777_059e-190, 4.912_107_071_492_364e-190),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 2e-14);
    }
}

#[test]
fn test_hyp2f1_c64_inverse_near_integer_cancelling() {
    // a - b is near an integer, but the paired 1/z series cancel more than the unpaired ones.
    // Arb at the binary64 inputs, compared componentwise.
    let actual = xsf::hyp2f1(75.25 + 2.0_f64.powi(-14), -4.75, 5.0, c64(-1.25, 0.125));
    let expected = [422_907.130_552_140_6, -197_134.634_126_622_1];
    np_assert_allclose!([actual.re, actual.im], expected, rtol = 4e-15);
}

#[test]
fn test_hyp2f1_c64_inverse_near_integer_unpaired() {
    // a - b is near an integer, but pairing the 1/z series gains little, or loses a tiny a or b.
    // Arb at the binary64 inputs, relative to the norm.
    for (a, b, c, z, expected, rtol) in [
        (
            5.024_997_972_554_356e-7,
            -1.361_721_274_401_763_5e-6,
            -30.717_500_717_391_21,
            c64(-1.201_105_588_470_700_7, -0.669_078_036_954_691_5),
            c64(0.999_999_999_999_972_9, -1.555_693_488_344_063e-14),
            1e-14,
        ),
        (
            2.889_925_064_146_597e-10,
            5.492_583_293_852_574e-13,
            -108.719_964_352_302_91,
            c64(-0.193_675_008_963_609_06, -3.649_948_169_840_155_5),
            c64(1.0, 5.377_948_821_171_109_6e-24),
            1e-13,
        ),
        (
            1.333_146_694_014_018_5e-10,
            1.376_623_639_782_201_2e-10,
            -225.100_664_606_829_72,
            c64(-0.930_006_774_382_229_3, -3.541_234_746_238_994_4),
            c64(1.0, 2.898_743_782_923_748e-22),
            2e-11,
        ),
        (
            9.294_168_821_359_962,
            4.294_168_821_359_3,
            11.1875,
            c64(0.667_471_980_816_785, 1.352_831_044_772_488),
            c64(0.082_811_595_041_503_58, -0.395_157_581_779_192_1),
            1e-13,
        ),
        (
            9.193_099_974_568_938,
            37.193_099_985_524_036,
            44.8125,
            c64(-2.848_635_621_977_817_5, 1.422_743_258_795_358_8),
            c64(-9.501_738_164_549_64e-6, 7.011_859_829_330_347e-7),
            1e-13,
        ),
        (
            -5.999_999_999_999_815,
            -5.999_999_998_420_666,
            0.163_827_140_645_315_34,
            c64(-4.219_983_348_289_068, 0.0),
            c64(-2_937.138_296_482_582, 0.0),
            1e-13,
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            let actual = xsf::hyp2f1(a, b, c, z);
            assert!((actual - expected).norm() <= rtol * expected.norm(), "{z}: {actual}");
        }
    }
}

/// Compares each component, so that a small one is checked too.
fn assert_components_close(actual: num_complex::Complex64, expected: num_complex::Complex64) {
    np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 3e-15);
}

#[test]
fn test_hyp2f1_c64_binomial() {
    // F(a, b; a; z) = F(b, a; a; z) = (1 - z)^(-b) (DLMF 15.4.6), Arb references, also at conj(z).
    for (a, b, z, expected) in [
        // Pearson, Olver and Porter case 13 continued off the real axis.
        (
            500.0,
            -500.0,
            c64(0.75, 0.01),
            c64(5.814_852_388_330_934e-302, -1.264_531_485_912_276_4e-301),
        ),
        // A large degree with a representable power.
        (2.5, -5e5, c64(1e-3, 1e-3), c64(-3.927_742_331_616_963e-218, 5.946_959_811_064_571e-218)),
        (2.5, 5e5, c64(1e-3, 1e-3), c64(-7.732_767_069_881_334e216, -1.170_811_400_298_151_8e217)),
        // Outside the unit disk, for equal and nonintegral parameters.
        (2.5, 2.5, c64(0.83, -1.11), c64(-0.687_674_549_201_301_9, 0.295_186_369_637_935_7)),
        (0.5, 0.5, c64(0.39, 1.12), c64(0.761_295_329_201_635_7, 0.452_252_310_326_859_64)),
        (
            300.0,
            299.0,
            c64(-2.18, 2.83),
            c64(-6.230_327_053_389_383e-189, -4.947_746_322_762_214e-189),
        ),
        (1.5, 2.25, c64(1e10, -1e10), c64(8.055_272_254_422_263e-24, 1.205_556_687_874_377_5e-23)),
        // Exact rational reference from (2 - i)^(-10).
        (10.0, 10.0, c64(-1.0, 1.0), c64(-2.426_88e-5, -3.190_784e-4)),
        // Integer exponents beyond 2^52, for a tiny z and on the unit circle around z = 1.
        (2.5, 1e17, c64(1e-18, 1e-18), c64(1.099_649_666_829_409_1, 0.110_332_988_730_203_72)),
        (2.5, 1e300, c64(1e-300, 1e-300), c64(1.468_693_939_915_885_1, 2.287_355_287_178_842_7)),
        (
            2.5,
            1e17,
            c64(0.5, 0.866_025_403_784_438_6),
            c64(-8.406_158_824_930_644, 76.657_099_442_722_17),
        ),
        (
            1.5,
            -9_007_199_254_740_992.0,
            c64(1.000_000_000_000_1, 1.0),
            c64(0.066_246_702_206_136_83, -0.997_803_274_466_835_7),
        ),
        (
            1.5,
            -9_007_199_254_740_992.0,
            c64(2.000_000_000_000_01, 0.0),
            c64(9.017_628_405_030_062e39, 0.0),
        ),
        // At a pole of c, where b terminates the series first.
        (-5.0, -2.0, c64(0.3, 0.2), c64(0.45, -0.28)),
        // On the cut, and next to it for a tiny ordinate.
        (1.5, 0.25, c64(3.0, 0.0), c64(0.594_603_557_501_360_5, 0.594_603_557_501_360_5)),
        (1.5, 0.25, c64(3.0, 1e-300), c64(0.594_603_557_501_360_5, 0.594_603_557_501_360_5)),
        // Next to the cut, the real part of a half-integer power.
        (
            2.5,
            -1000.5,
            c64(3.0, 2.0_f64.powi(-1000)),
            c64(707.460_334_577_140_8, -1.515_342_004_482_324_3e301),
        ),
        // The power amplifies the smallest subnormal ordinate to a normal component.
        (
            1.0,
            -512.0,
            c64(-0.5, DENORM_MIN),
            c64(1.441_201_271_173_911e90, -2.430_457_299_045_594e-231),
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            assert_components_close(xsf::hyp2f1(a, b, a, z), expected);
            assert_components_close(xsf::hyp2f1(b, a, a, z), expected);
        }
    }
}

#[test]
fn test_hyp2f1_c64_binomial_exact() {
    for sign in [-1.0, 1.0] {
        // Exact powers on and next to the line Re z = 1.
        for (y, expected) in [(1.0, -1.0), (10.0, -1e-10)] {
            let value = xsf::hyp2f1(10.0, 10.0, 10.0, c64(1.0, sign * y));
            assert_eq!(value.re, expected);
            assert_eq!(value.im.to_bits(), 0);
        }
        let value = xsf::hyp2f1(10.0, 10.0, 10.0, c64(1.0 + f64::EPSILON, sign));
        assert_eq!(value.im, -sign * 10.0 * f64::EPSILON);
        // A half-integer power further from the cut.
        let actual = xsf::hyp2f1(2.5, -0.5, 2.5, c64(2.0, sign * 2.0_f64.powi(-100)));
        assert_eq!(actual, c64(2.0_f64.powi(-101), -sign));
    }
    // A component of 1 - z below the rounding of the other one reaches the power.
    let z = c64(1.0 + f64::EPSILON, 2.0_f64.powi(1023));
    assert_eq!(xsf::hyp2f1(2.5, -1.0, 2.5, z), c64(-f64::EPSILON, -z.im));
    let square = c64(f64::NEG_INFINITY, 2.0 * f64::EPSILON * z.im);
    assert_eq!(xsf::hyp2f1(2.5, -2.0, 2.5, z), square);
}

#[test]
fn test_hyp2f1_c64_binomial_overflow() {
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    for (b, z, expected) in [
        (-3000.0, c64(3.0, 1.0), c64(-inf, inf)),
        (3000.0, c64(3.0, 1.0), c64(0.0, 0.0)),
        // A power that is linear in an ordinate below its own reciprocal.
        (-(2.0_f64.powi(564)), c64(-(2.0_f64.powi(1000)), 2.0_f64.powi(-1000)), c64(inf, -inf)),
        // Beyond the accepted error, only a real power or an underflow is certain.
        (1e300, c64(0.5, 0.0), c64(inf, 0.0)),
        (-1e300, c64(0.3, 0.2), c64(0.0, 0.0)),
        (1e300, c64(0.3, 0.2), c64(nan, nan)),
        (1e300, c64(0.0, 1e-290), c64(nan, nan)),
        // An infinite parameter does not cancel.
        (inf, c64(3.0, 0.0), c64(nan, nan)),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            let actual = xsf::hyp2f1(2.5, b, 2.5, z);
            let same = actual == expected || actual.is_nan() && expected.is_nan();
            assert!(same, "b = {b}, z = {z}: {actual}");
        }
    }
}

#[test]
fn test_hyp2f1_c64_inverse_near_integer_index() {
    // Arb references at the exact binary64 arguments.
    for (a, b, c, z, expected) in [
        (
            1.0,
            2.0,
            f64::from_bits(3.0_f64.to_bits() - 1),
            c64(0.85, 0.71),
            c64(0.893_382_297_136_081_9, 1.038_774_367_667_433_7),
        ),
        (
            1.5,
            2.5,
            f64::from_bits(3.5_f64.to_bits() + 1),
            c64(1.25, 0.5),
            c64(-1.037_930_983_193_957, 2.216_972_649_752_144_3),
        ),
    ] {
        for (a, b) in [(a, b), (b, a)] {
            for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
                let actual = xsf::hyp2f1(a, b, c, z);
                np_assert_allclose!(
                    [actual.re, actual.im],
                    [expected.re, expected.im],
                    rtol = 2e-14
                );
            }
        }
    }

    // The integer power at the switch must not add a phase on the negative real axis.
    let c = f64::from_bits(6.5_f64.to_bits() - 1);
    assert_eq!(xsf::hyp2f1(3.5, 2.5, c, c64(-1.2, 0.0)).im, 0.0);
    let actual = xsf::hyp2f1(3.5, 2.5, c, c64(-1.2, 1e-30));
    np_assert_allclose!(
        [actual.re, actual.im],
        [0.312_448_405_332_954_3, 2.275_890_463_150_229_3e-31],
        rtol = 1e-11,
        atol = 0.0
    );
}

#[test]
fn test_hyp2f1_c64_triangle_map() {
    // Schwarz–Christoffel maps for thin triangles. Arb references at the binary64 inputs.
    for (a, d, x, y, re, im) in [
        (0.5, 0.01 / 180.0, 1.001, 0.001, 3.972_049_553_641_628, 1.175_344_862_076_758_1),
        (0.375, 0.1 / 180.0, 1.05, 0.0, 1.915_689_474_120_120_7, 1.154_794_618_703_662),
        (0.375, 0.1 / 180.0, 1.05, DENORM_MIN, 1.915_689_474_120_120_7, 1.154_794_618_703_662),
        (
            0.5,
            1e-7,
            0.875,
            0.484_122_918_275_927_2,
            1.195_306_715_757_189_3,
            0.502_528_800_257_228_6,
        ),
        (0.125, 1.0 / 64.0, 0.9001, 0.001, 1.262_194_663_684_517_4, 0.001_165_653_720_072_697_2),
    ] {
        let b = 1.0 - d;
        let c = 1.0 + a;
        for (a, b) in [(a, b), (b, a)] {
            for sign in [-1.0, 1.0] {
                let value = xsf::hyp2f1(a, b, c, c64(x, sign * y));
                np_assert_allclose!([value.re, value.im], [re, sign * im], rtol = 3e-15);
            }
        }
    }

    // A subnormal phase below the cut survives through the beta derivative.
    for y in [-8.3695e-320, 8.3695e-320] {
        let value = xsf::hyp2f1(
            0.964_843_75,
            1.0 - f64::EPSILON,
            1.964_843_75,
            c64(0.999_999_971_191_316_1, y),
        );
        assert!((value.im - 2.803_053_809_78e-312_f64.copysign(y)).abs() <= 2.0 * DENORM_MIN);
    }
}

#[test]
fn test_hyp2f1_c64_rounded_integer_balance() {
    for (a, b, c, x, y, re, im) in [
        (
            8.000_000_000_000_002,
            8.0,
            16.0,
            1.0001,
            1e-6,
            205_297.257_891_107_5,
            160_165.783_236_770_77,
        ),
        (
            32.0,
            31.999_999_999_999_996,
            64.0,
            1.0,
            DENORM_MIN,
            2.159_228_790_004_916_6e22,
            4.605_886_830_369_123e19,
        ),
        // Inside the unit disk, and with a rounded balance of -2.
        (23.9, 13.35, 37.25, 0.99949, 0.0011, 20_550_478_649.543_236, 24_851_828_396.636_16),
        (23.9, 13.35, 35.25, 0.99949, 0.0011, -15_283_485_606_553.012, 11_862_031_344_500.123),
        // Outside it, two ulps from a balance of 5.
        (
            -4.3,
            3.44,
            4.140_000_000_000_001,
            1.02,
            0.05,
            0.009_737_815_581_243_74,
            -0.001_776_337_350_623_636,
        ),
        // A balance of 1e-13, beyond the rounding of a + b.
        (
            1.3,
            2.1,
            3.400_000_000_000_1,
            0.95,
            0.1,
            3.785_068_120_859_939_3,
            2.306_708_495_240_381_6,
        ),
        // Gamma functions that overflow or underflow, also for an exact balance.
        (
            100.3,
            80.4,
            180.7,
            0.999_999,
            1e-7,
            8.216_080_747_133_156e54,
            2.119_177_392_195_293_8e53,
        ),
        (150.5, 60.25, 210.75, 0.9999, 1e-5, 4.019_727_082_272_486e54, 4.690_374_590_829_056e53),
        (
            181.719,
            -17.3476,
            162.3714,
            0.998,
            0.004,
            3.770_952_191_815_918e-22,
            -9.784_596_732_535_692e-23,
        ),
        // A rising factorial that overflows in one order of a and b.
        (
            -169.25,
            173.5,
            170.25,
            1.0001,
            1e-40,
            1.879_389_312_714_306e-110,
            -3.344_377_011_037_522e-148,
        ),
        // The Maclaurin series is closer than the limit.
        (
            -0.515_978_988_132_843_6,
            13.328_394_588_502_42,
            11.812_415_600_390_873,
            0.952_516_937_091_452_1,
            0.197_088_629_712_879_52,
            0.251_208_055_225_539_1,
            -0.388_357_387_394_512_25,
        ),
    ] {
        for (a, b) in [(a, b), (b, a)] {
            for sign in [-1.0, 1.0] {
                let actual = xsf::hyp2f1(a, b, c, c64(x, sign * y));
                np_assert_allclose!([actual.re, actual.im], [re, sign * im], rtol = 3e-12);
            }
        }
    }
    // The limit cancels eight digits, more than a balance of 1e-7 leaves.
    let actual = xsf::hyp2f1(40.0, 60.0, 112.000_000_1, c64(0.97, 0.0));
    np_assert_allclose!([actual.re], [3.906_774_039_774_269e17], rtol = 3e-12);
    // A term of the limit series vanishes before the series has converged.
    let actual = xsf::hyp2f1(8.0, 8.0, 16.000_000_1, c64(0.967_797_599_699_437_9, 0.0));
    np_assert_allclose!([actual.re], [6_343.004_188_994_46], rtol = 3e-12);
    // The same, at an exact balance where the gamma ratio overflows.
    let actual = xsf::hyp2f1(2.0, 170.0, 172.0, c64(0.989_486_905_267_847_8, 0.0));
    np_assert_allclose!([actual.re], [2_842.713_082_597_145_3], rtol = 3e-12);
}

#[test]
fn test_hyp2f1_c64_near_zero_balance() {
    // Balances from 1e-15 to 1e-3, where the two 1 - z series cancel, also near negative
    // integers a and b.
    for (a, b, c, x, y, re, im) in [
        (
            35.906_427_470_934_794,
            35.918_109_483_734_87,
            71.824_536_954_569_67,
            0.977_438_975_071_172,
            -0.005_209_794_246_290_756,
            1.083_548_523_330_061_6e17,
            -2.213_056_030_691_095e17,
        ),
        (
            13.979_526_584_792_934,
            31.593_534_152_054_417,
            45.573_060_736_847_346,
            0.947_394_922_390_871_2,
            -0.185_303_445_597_178_76,
            -394_690.911_800_178_9,
            -2_700_337.709_037_638_7,
        ),
        (
            19.751_104_662_473_47,
            32.354_129_416_194_446,
            52.105_244_078_667_92,
            0.999_809_355_316_824_6,
            -8.586_293_610_241_832e-5,
            1.848_063_021_882_126_2e15,
            -4.442_273_256_723_710_6e14,
        ),
        (
            1.849_560_119_568_044,
            2.709_219_916_618_709_4,
            4.558_680_036_186_754,
            1.009_122_957_134_769_3,
            0.434_559_770_297_138_1,
            0.366_264_253_368_914_9,
            3.030_751_344_512_619_7,
        ),
        (
            -19.273_999_226_725_604,
            -13.346_483_832_781_338,
            -32.619_483_059_506_94,
            0.998_137_980_533_285_1,
            0.000_140_539_004_778_811_05,
            2.629_387_207_616_451_7e-9,
            -9.556_512_838_089_924e-11,
        ),
        (
            8.000_000_000_001,
            8.0,
            16.0,
            1.0001,
            1e-6,
            205_297.257_891_800_95,
            160_165.783_237_828_03,
        ),
        (
            -31.999_999_978_277_355,
            -30.264_032_172_974_95,
            -62.264_032_151_236_54,
            0.952_624_109_802_563,
            0.127_226_516_469_02,
            -5.802_408_139_152_598_5e-12,
            2.360_913_625_077_687_5e-12,
        ),
        (
            -32.999_999_998_376_8,
            -25.836_224_093_355_025,
            -58.836_224_091_731_815,
            1.042_184_402_821_436_2,
            0.293_670_970_687_779_53,
            -5.307_870_591_538_909e-10,
            1.983_824_573_521_540_5e-9,
        ),
    ] {
        for (a, b) in [(a, b), (b, a)] {
            for sign in [-1.0, 1.0] {
                let actual = xsf::hyp2f1(a, b, c, c64(x, sign * y));
                np_assert_allclose!([actual.re, actual.im], [re, sign * im], rtol = 3e-12);
            }
        }
    }
}

#[test]
fn test_hyp2f1_c64_near_zero_balance_edges() {
    // A vanishing term at a balance of 1/8.
    let actual = xsf::hyp2f1(4.375, 4.375, 8.875, c64(0.904_973_484_375_869_2, 0.0));
    np_assert_allclose!([actual.re], [43.253_829_995_862_97], rtol = 1e-13);
    // No hang: a far below the term limit returns, but as NaN (a known limitation).
    let actual = xsf::hyp2f1(-16_777_216.5, 16_777_217.5, 1.0 + 2_f64.powi(-26), c64(0.95, 0.0));
    assert!(actual.is_nan());
    // An overflowing gamma ratio, at a balance of 2^-43.
    let actual = xsf::hyp2f1(500.0, 500.0, 1000_f64.next_up(), c64(1.0, 1e-8));
    let desired = [3.272_480_137_841_555e302, 1.049_871_344_510_691_8e302];
    np_assert_allclose!([actual.re, actual.im], desired, rtol = 3e-12);
    // Near a negative integer a, beyond the direct gamma ratio, at a balance of 1/8.
    let actual =
        xsf::hyp2f1(-100.999_999_999_999_99, -0.5, -101.374_999_999_999_99, c64(1.0, 1e-6));
    np_assert_allclose!([actual.re], [0.080_922_750_731_539_07], rtol = 1e-13);
    np_assert_allclose!([actual.im], [-4.670_398_751_111_399e-6], rtol = 1e-10);
    // A subnormal imaginary part of z.
    let actual = xsf::hyp2f1(20.0, 20.0, 40.000_000_000_000_01, c64(0.999, 1e-313));
    np_assert_allclose!([actual.re], [855_282_174_438.874_3], rtol = 1e-13);
    np_assert_allclose!([actual.im], [6.987_929_090_888_118e-299], rtol = 1e-13);
    // Near a negative integer a, beyond the direct gamma ratio, at a balance of 1e-14.
    let actual = xsf::hyp2f1(-98.999_999_99, -1.5, -100.499_999_989_999_99, c64(0.95, 0.0));
    np_assert_allclose!([actual.re], [0.016_051_788_054_881_427], rtol = 1e-14);
    // A small imaginary part, which the rotation of the paired sum cancels.
    let actual = xsf::hyp2f1(1.5, 2.5, 4.125, c64(1.0, 1e-100));
    let desired = [40.577_947_190_182_54, 3.211_635_692_333_511_7e-12];
    np_assert_allclose!([actual.re, actual.im], desired, rtol = 1e-14);
    // An underflowing gamma ratio of -5.1e-257 - 7.2e-257i returns NaN (a known limitation).
    let actual = xsf::hyp2f1(-600.25, -600.25, -1_200.499_999, c64(0.95, 0.001));
    assert!(actual.is_nan());
    // The Maclaurin series replaces the paired sum, at a balance of 2^-49.
    let actual = xsf::hyp2f1(6.0, 6.0, 12_f64.next_up(), c64(0.95, 0.0));
    np_assert_allclose!([actual.re], [459.334_127_269_937_5], rtol = 3e-14);
    // A log ratio of NaN, at a balance of 1/8.
    let actual = xsf::hyp2f1(-0.0625, 1.5, 1.5625, c64(0.99, 0.001));
    let desired = [0.778_651_057_800_017_5, -0.003_664_795_346_346_866_8];
    np_assert_allclose!([actual.re, actual.im], desired, rtol = 1e-14);
}

#[test]
fn test_hyp2f1_c64_large_near_zero_balance() {
    // Large a and b at a distance from 1, where the 1 - z and Maclaurin series cancel, at an
    // exact balance
    for (a, b, c, x, y, re, im) in [
        (
            39.358_066_203_017_785,
            37.510_684_884_703_18,
            76.868_751_087_720_97,
            0.989_591_664_228_343_8,
            -0.268_876_794_115_152_03,
            152_914_339_631.575_23,
            -213_148_183_992.845_9,
        ),
        (
            26.812_186_340_530_424,
            35.466_637_734_956_11,
            62.278_824_075_486_53,
            0.952_293_724_093_895_6,
            0.280_160_557_509_086_25,
            -231_308_370.703_926_53,
            282_156_881.884_880_07,
        ),
    ] {
        for sign in [1.0, -1.0] {
            let actual = xsf::hyp2f1(a, b, c, c64(x, sign * y));
            np_assert_allclose!([actual.re, actual.im], [re, sign * im], rtol = 1e-14);
        }
    }
    // Both lips of the cut, at an exact balance.
    let (a, b, c) = (32.442_996_112_826_41, 39.907_764_779_431_275, 72.350_760_892_257_68);
    for y in [0.0, -0.0] {
        let actual = xsf::hyp2f1(a, b, c, c64(1.055_851_490_374_344_6, y));
        np_assert_allclose!([actual.re], [-2.088_749_453_446_372e19], rtol = 1e-13);
        let im = 6.059_741_666_724_565e20_f64.copysign(-y);
        np_assert_allclose!([actual.im], [im], rtol = 1e-14);
    }
    // Small a with large b and c, where a continuation from z = 1/2 loses too much.
    for y in [0.3125, -0.3125] {
        let actual = xsf::hyp2f1(0.5, 96.25, 35.625, c64(0.925, y));
        let desired = [1_161.620_547_299_975_7, 1_680.282_191_320_187_6_f64.copysign(y)];
        np_assert_allclose!([actual.re, actual.im], desired, rtol = 1e-12);
    }
    // An integer balance close to 1, where the Maclaurin sums at 1/2 lose about 40 eps
    let (a, c) = (123.415_688_700_992_7, 245.832_377_401_985_4);
    let actual = xsf::hyp2f1(a, a, c, c64(0.997_013_321_126_109_9, 0.000_193_262_214_104_612));
    let desired = [7.879_161_399_246_598e69, 3.898_993_388_042_795e69];
    np_assert_allclose!([actual.re, actual.im], desired, rtol = 1e-13);
    // A continuation near 1 that does not beat the limit for a nearly integer balance
    let actual = xsf::hyp2f1(
        62.377_025_825_687_795,
        13.886_599_054_008_652,
        47.263_624_879_696_45,
        c64(0.901_941_236_573_627_1, 0.278_212_790_916_568_1),
    );
    let desired = [4_795_717_777.401_794, -4_480_686_255.456_943_5];
    np_assert_allclose!([actual.re, actual.im], desired, rtol = 5e-14);
}

#[test]
fn test_hyp2f1_c64_large_near_zero_balance_inverse() {
    // Large a and b outside the unit circle, where the 1/z and Maclaurin series cancel: an exact
    // balance and an integer a - b.
    for (a, b, c, x, y, re, im) in [
        (
            35.672_428_368_194_73,
            39.25,
            74.922_428_320_253_62,
            1.122_345_052_724_406_5,
            -0.217_285_801_179_134_7,
            22_818_029_254_577.574,
            3_152_305_726_487.289,
        ),
        (
            34.25,
            33.25,
            67.5,
            1.120_760_142_922_645_3,
            0.225_001_041_654_064_87,
            -646_162_001_297.730_7,
            -415_616_178_762.879_1,
        ),
    ] {
        for sign in [1.0, -1.0] {
            let actual = xsf::hyp2f1(a, b, c, c64(x, sign * y));
            np_assert_allclose!([actual.re, actual.im], [re, sign * im], rtol = 1e-14);
        }
    }
    // An integer a - b with a tiny a, and an imaginary part below the error bound
    let actual =
        xsf::hyp2f1(1e-20, 150.0, 112.5, c64(1.156_055_063_705_178_3, 0.340_986_535_059_630_6));
    np_assert_allclose!([actual.re, actual.im], [1.0, 2.447_334_736_479_153_5e-20], rtol = 1e-14);
    // If finite, then accurate (these give NaN, as on master, so this guards against wrong values).
    let tiny = 2.0_f64.powi(-538);
    for (a, b, c, x, y, re, im) in [
        (1.25, 331.25, 128.125, 1.0625, 0.3125, -7.248_754_802_291_216_5, -184.468_986_796_544_4),
        (
            1e-20,
            700.0,
            240.125,
            1.125,
            0.125,
            1.887_566_949_173_823_6e116,
            8.835_528_861_246_3e116,
        ),
        (0.5, 0.5, 10_000.5, 1.125, 0.125, 1.000_028_127_109_545_5, 3.125_634_892_208_265e-6),
        (0.5, 0.5, 1e100, 1.125, 0.125, 1.0, 3.125e-102),
        // (a b underflows)
        (tiny, tiny, 5e-324, 1.125, 0.125, 1.433_216_987_849_965_9, 0.589_048_622_548_086_2),
        (
            2.0 * tiny,
            3.0 * tiny,
            5e-324,
            1.125,
            0.125,
            3.599_301_927_099_795,
            3.534_291_735_288_517_3,
        ),
        (tiny, tiny, 5e-324, 1.125, 0.0, 1.519_860_385_419_959, FRAC_PI_4),
        (tiny, tiny, 5e-324, 1.125, -0.0, 1.519_860_385_419_959, -FRAC_PI_4),
    ] {
        let actual = xsf::hyp2f1(a, b, c, c64(x, y));
        if actual.is_finite() {
            np_assert_allclose!([actual.re, actual.im], [re, im], rtol = 1e-13);
        }
    }
    // A weak imaginary part on both lips of the cut, and just off it, keeps its sign.
    for y in [0.0_f64, -0.0, 1e-100, -1e-100] {
        let actual = xsf::hyp2f1(0.5, 0.5, 20.5, c64(1.125, y));
        np_assert_allclose!([actual.re], [1.014_626_580_089_177_6], rtol = 1e-11);
        assert_eq!(actual.im.is_sign_negative(), y.is_sign_negative(), "{actual}");
    }
}

#[test]
fn test_hyp2f1_c64_euler_tiny_cut_lips() {
    let a = 32.000_000_000_000_01;
    for (a, b) in [(a, 106.0), (106.0, a)] {
        for y in [0.0, DENORM_MIN, 1e-40] {
            for sign in [-1.0, 1.0] {
                let actual = xsf::hyp2f1(a, b, 97.0, c64(2.0, sign * y));
                np_assert_allclose!(
                    [actual.re, actual.im],
                    [0.000_370_569_924_328_136, sign * 8.271_994_658_390_709e-18],
                    rtol = 5e-15
                );
            }
        }
    }
}

#[test]
fn test_hyp2f1_c64_integral_euler_power() {
    // (1 - 3z/5) / (1 - z)^2.
    for y in [0.0, -0.0, 1e-30, -1e-30] {
        let actual = xsf::hyp2f1(3.5, 1.0, 2.5, c64(2.0, y));
        np_assert_allclose!([actual.re, actual.im], [-0.2, -0.2 * y], rtol = 5e-15);
    }

    let actual = xsf::hyp2f1(4097.0, 1024.0, 4096.0, c64(2.0, 0.0));
    assert_eq!(actual, c64(0.5, 0.0));
    let actual = xsf::hyp2f1(4097.0, 1024.0, 4096.0, c64(1.0, 1.0));
    assert_eq!(actual, c64(0.75, 0.25));

    let z = c64(1.5, 0.866_025_403_784_438_6);
    let actual = xsf::hyp2f1(8_589_934_593.0, 2_147_483_648.0, 8_589_934_592.0, z);
    np_assert_allclose!(
        [actual.re, actual.im],
        [-0.124_999_976_671_243_7, -0.649_519_120_182_641_7],
        rtol = 5e-15
    );
}

#[test]
fn test_hyp2f1_c64_euler_power_overflow() {
    // Arb at the binary64 inputs, where a half power near the cut overflows.
    for (a, b, c, z, expected) in [
        (
            -3.019_339_683_446_933_3,
            -27.5,
            -12.019_339_683_446_933,
            c64(6.804_317_682_837_257e29, 8.095e-320),
            c64(f64::INFINITY, f64::NEG_INFINITY),
        ),
        (
            -97.5,
            -289.5,
            -301.5,
            c64(5.886_203_089_524_244e19, 2.4e-181),
            c64(f64::NEG_INFINITY, f64::INFINITY),
        ),
        (
            -206.522_940_928_326_3,
            825.5,
            -440.522_940_928_326_3,
            c64(1.656_692_523_841_373_4, -1e-40),
            c64(f64::NEG_INFINITY, f64::INFINITY),
        ),
    ] {
        assert_eq!(xsf::hyp2f1(a, b, c, z), expected);
    }

    // Arb, where only the real part overflows.
    for im in [-5e-324_f64, 5e-324] {
        let z = c64(1.5, im);
        let actual = xsf::hyp2f1(-60.282_220_454_807_145, 697.0, -304.282_220_454_807_15, z);
        assert_eq!(actual.re, f64::NEG_INFINITY);
        np_assert_allclose!([actual.im], [1.754_357_958_579_083_2e176 * im.signum()], rtol = 1e-5);
    }
}

#[test]
fn test_hyp2f1_c64_euler_power_underflow() {
    // (1 - 15z/16) / (1 - z)^2, where (1 - z)^-2 underflows.
    let z = c64(-2.0_f64.powi(1021), 0.0);
    for (a, b) in [(1.0, 17.0), (17.0, 1.0)] {
        assert_eq!(xsf::hyp2f1(a, b, 16.0, z), c64(4.172_013_484_701_002_6e-308, 0.0));
    }

    // (1 - z)^(-A - 1) (1 + (A/16 - 1) z) with A = 720 2^52 (Arb), where -A - 1 rounds.
    let big = 720.0 * 2.0_f64.powi(52);
    let z = c64(1.0, 1.0 + f64::EPSILON);
    for (a, b) in [(big, 17.0), (17.0, big)] {
        let actual = xsf::hyp2f1(a, b, 16.0, z);
        np_assert_allclose!(
            [actual.re, actual.im],
            [-4.118_559_248_038_33e-296, 4.118_559_248_038_329e-296],
            rtol = 1e-12
        );
    }

    assert!(xsf::hyp2f1(f64::INFINITY, 17.0, 16.0, c64(0.5, 0.0)).is_nan());
}

#[test]
fn test_hyp2f1_c64_euler_polynomial_cancellation() {
    // Arb at the exact binary64 inputs. The first is ((1 - z)^-37 - 1) / (37 z).
    for (a, b, c, z, expected) in [
        (
            38.0,
            1.0,
            2.0,
            c64(2.0, 0.25),
            c64(-0.009_058_578_399_532_28, 0.002_685_488_481_771_803_6),
        ),
        (
            35.0,
            0.5,
            2.0,
            c64(1.875, 0.125),
            c64(0.073_108_585_745_817_23, 0.131_363_053_533_884_8),
        ),
        (
            144.0,
            2.500_000_000_000_000_4,
            5.0,
            c64(0.53, 0.8875),
            c64(7.398_391_777_232_239e-5, -6.000_132_525_033_908e-5),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        assert!((actual - expected).norm() <= 2e-14 * expected.norm());
    }
}

#[test]
fn test_hyp2f1_c64_euler_polynomial_tiny_a() {
    // Arb at the exact binary64 inputs. Here c - b is a negative integer, and c - a rounds to c
    // (#2247).
    for (a, b, c, x, expected, rtol) in [
        (1e-20, 300.5, 20.5, 0.5, 5.702_567_468_803_259e37, 2e-16),
        (-1e-20, 300.5, 20.5, 0.5, -5.702_567_468_803_259e37, 2e-16),
        (300.5, 1e-20, 20.5, 0.5, 5.702_567_468_803_259e37, 2e-16),
        (1e-20, 300.5, 20.5, 0.9, 1.602_722_697_015_226_8e228, 2e-16),
        // (here the exponent c - a - b of the Euler factor rounds)
        (1e-10, 300.5, 20.5, 0.5, 5.702_567_472_298_491_4e47, 2e-14),
    ] {
        let actual = xsf::hyp2f1(a, b, c, c64(x, 0.0));
        np_assert_allclose!([actual.re, actual.im], [expected, 0.0], rtol = rtol, atol = 0.0);
    }
    for (a, expected, rtol) in [
        (1e-20, c64(9.366_954_165_364_594e34, 1.251_049_345_038_304_3e35), 2e-14),
        (1e-10, c64(9.366_954_170_578_047e44, 1.251_049_345_846_006_5e45), 1e-13),
    ] {
        let actual = xsf::hyp2f1(a, 300.5, 20.5, c64(0.5, 0.1));
        np_assert_allclose!(
            [actual.re, actual.im],
            [expected.re, expected.im],
            rtol = rtol,
            atol = 0.0
        );
    }
}

#[test]
fn test_hyp2f1_c64_pfaff_polynomial_arb() {
    // Arb at the exact binary64 inputs, for both orders of (a, b): (a, b, c, z, expected, rtols).
    for (a, b, c, z, expected, rtol) in [
        // The terms in z cancel.
        (-26.0, 29.5, 27.75, c64(0.781_25, 0.0), c64(3.890_826_082_506_454e-16, 0.0), [1e-14; 2]),
        (
            -19.0,
            -35.5,
            -36.75,
            c64(0.95, 0.25),
            c64(1.236_696_750_629_711_3e-11, 6.406_230_812_945_963e-12),
            [1e-14; 2],
        ),
        (
            -40.0,
            33.0,
            9.5,
            c64(0.95, 0.125),
            c64(-1.030_992_516_789_373_6e-11, -7.086_575_071_894_045e-13),
            [1e-14, 1e-13],
        ),
        (-40.0, 33.0, 9.5, c64(0.95, 0.0), c64(6.340_337_088_462_662e-20, 0.0), [1e-14; 2]),
        // Only the real part of the Pfaff polynomial is resolved.
        (
            -18.0,
            0.75,
            0.5,
            c64(1.0 - f64::EPSILON, 2.0_f64.powi(-450)),
            c64(-0.042_035_401_442_955_29, 1.165_295_602_809_008_4e-137),
            [5e-16, 1e-8],
        ),
        // 1 - 6.4x + 9.6x^2 near x = 1/4, resolved to a quarter of its value, 0.8 * 2^-53.
        (
            -2.0,
            32.0,
            10.0,
            c64(0.25 - f64::EPSILON / 4.0, 0.0),
            c64(8.881_784_197_001_255e-17, 0.0),
            [0.3, 1e-14],
        ),
        // c - a rounds.
        (
            38.5,
            -160.0,
            -21.999_999_999,
            c64(0.728_030_478_673_954_3, 0.867_436_485_373_291_8),
            c64(-1.031_291_498_685_774_2e75, -1.378_021_763_441_624_9e75),
            [2e-14; 2],
        ),
        // The origin sum keeps the imaginary part.
        (
            -17.0,
            11.25,
            -3.078_209_853_523_439_2,
            c64(1.000_000_000_000_005_8, 2.780_885_689_443_515e-15),
            c64(-33.311_434_158_034_42, 1.059_732_248_933_062e-11),
            [1e-14, 1e-4],
        ),
        (
            -21.5,
            -16.0,
            -18.5,
            c64(1.000_000_000_718_648_3, -4.301_370_674_414_646e-9),
            c64(0.383_730_473_252_956_1, -3.154_416_939_670_981e-8),
            [1e-10; 2],
        ),
    ] {
        for (a, b) in [(a, b), (b, a)] {
            let actual = xsf::hyp2f1(a, b, c, z);
            np_assert_allclose!([actual.re], [expected.re], rtol = rtol[0]);
            np_assert_allclose!([actual.im], [expected.im], rtol = rtol[1]);
        }
    }

    // (1 - z)^(-1/2) F(-256, 1/2; c; z/(z - 1)) (DLMF 15.8.1), Arb at the binary64 inputs, in this
    // order only (#2268).
    let c = 2.0_f64.powi(60);
    let actual = xsf::hyp2f1(c + 256.0, 0.5, c, c64(2.0, 1e-10));
    np_assert_allclose!([actual.re], [4.999_999_999_999_998e-11], rtol = 1e-4);
    np_assert_allclose!([actual.im], [0.999_999_999_999_999_8], rtol = 1e-15);
}

#[test]
fn test_hyp2f1_c64_pfaff_polynomial_exact() {
    let bits = |z: num_complex::Complex64| [z.re.to_bits(), z.im.to_bits()];

    // Positive zeros: 1 + (-32/3)(3/32), and 1 - 6x + 8x^2 at x = 1/4.
    assert_eq!(bits(xsf::hyp2f1(-1.0, 32.0, 3.0, c64(0.093_75, 0.0))), [0, 0]);
    for (a, b) in [(-2.0, 15.0), (15.0, -2.0)] {
        assert_eq!(bits(xsf::hyp2f1(a, b, 5.0, c64(0.25, 0.0))), [0, 0]);
    }

    // A positive zero imaginary part on both lips of the cut.
    for (a, b) in [(-19.0, -4.5), (-4.5, -19.0)] {
        for im in [0.0, -0.0] {
            let actual = xsf::hyp2f1(a, b, -7.5, c64(2.0, im));
            np_assert_allclose!([actual.re], [103.486_713_286_713_28], rtol = 5e-16);
            assert_eq!(actual.im.to_bits(), 0);
        }
    }

    // (1 - z)(1 - 3z), where (1 - z)^2 underflows in the imaginary part.
    let tiny = 2.0_f64.powi(-600);
    assert_eq!(xsf::hyp2f1(-2.0, 2.0, 1.0, c64(1.0, tiny)), c64(0.0, 2.0 * tiny));

    // The same, with 1 - z the double just above zero in the real part and the subnormal 2^-1074
    // in the imaginary part, so that the imaginary part is the subnormal 2^-1073.
    let z = c64(1.0 - f64::EPSILON / 2.0, DENORM_MIN);
    let actual = xsf::hyp2f1(-2.0, 2.0, 1.0, z);
    np_assert_allclose!([actual.re], [-f64::EPSILON], rtol = 5e-16);
    assert_eq!(actual.im, f64::from_bits(2));
}

#[test]
fn test_hyp2f1_c64_pfaff_polynomial_limits() {
    // Rational arithmetic gives -0.008 + 0.0192i at the binary64 inputs, but 2(c + 1)
    // overflows: NaN.
    assert!(xsf::hyp2f1(-2.0, 1.2e308, 1e308, c64(0.9, 0.1)).is_nan());

    // Both components overflow.
    let z = c64(9.324_507_004_751_248e30, -1.997_223_793_753_443_4e31);
    let actual = xsf::hyp2f1(-20.5, 21_895.476_365_024_74, 21_885.476_365_024_74, z);
    assert_eq!(actual, c64(f64::NEG_INFINITY, f64::NEG_INFINITY));
}

#[test]
fn test_hyp2f1_c64_transform1_maclaurin() {
    // unlike in xsf: the Maclaurin series if the 1 - z transformation cancels (mpmath references)
    for (a, b, c, z, expected) in [
        // c - a - b = 50 (the limiting case) and 50.5
        (4.5, 2.5, 57.0, c64(0.64, 0.65), c64(1.122_209_544_812_728_6, 0.158_257_844_232_734_11)),
        (4.5, 2.5, 57.5, c64(0.64, 0.65), c64(1.121_218_275_094_999_2, 0.156_607_234_073_543_14)),
        // |z| > 1
        (3.5, 2.5, 40.0, c64(0.9, 0.6), c64(1.211_419_901_187_38, 0.188_557_093_819_745_4)),
        // c - a - b just misses an integer, so the 1 - z transformation is NaN
        (0.1, 0.2, 50.3, c64(0.64, 0.65), c64(1.000_254_319_428_871_4, 2.627_839_705_956_765e-4)),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        assert!(
            (actual - expected).norm() < 1e-14 * expected.norm(),
            "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
        );
    }
}

#[test]
fn test_hyp2f1_c64_huge_num_terms() {
    // unlike in xsf: these summed up to 2^64 terms (mpmath references)
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        for (a, b, c, z, expected) in [
            (4.39, 2.14, 1e17, c64(0.64, 0.65), c64(1.0, 6.106_49e-17)),
            (4.39, 2.14, 1e300, c64(0.64, 0.65), c64(1.0, 6.106_49e-300)),
            (
                0.5,
                0.25,
                1e6 + 0.75,
                c64(0.64, 0.65),
                c64(1.000_000_079_999_938_5, 8.125_003_656_241_712e-8),
            ),
            // the most terms that are still summed: 2^20
            (-1_048_575.0, 1.0, 2.0, c64(1e-7, 0.0), c64(0.949_356_716_098_484_6, 0.0)),
        ] {
            let actual = xsf::hyp2f1(a, b, c, z);
            assert!(
                (actual - expected).norm() < 1e-15 * expected.norm(),
                "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
            );
        }
        // these exceed the cap on the number of terms
        for (a, b, c, z) in [
            (-1_048_576.0, 1e-30, 1.0, c64(-1e-4, 0.0)),
            (-1e15, 1.0, 2.0, c64(0.5, 0.0)),
            (1e17, 2.0, 4.39, c64(3.0, 0.65)),
            (1e9, 2.0, 4.39, c64(3.0, 0.65)),
            (2_251_799_813_685_248.5, 0.25, 0.75, c64(0.95, 0.1)),
        ] {
            assert!(xsf::hyp2f1(a, b, c, z).re.is_nan());
        }
        tx.send(()).unwrap();
    });
    rx.recv_timeout(std::time::Duration::from_secs(10)).unwrap();
}

#[test]
fn test_hyp2f1_f64_c_minus_a_minus_b_le_minus_1() {
    // the Euler transformation cancels here (mpmath)
    let actual = [xsf::hyp2f1(-0.5, 125.0, 62.5, 0.4), xsf::hyp2f1(-0.8, 19.7, 9.4, 0.9)];
    let desired = [0.453_216_679_866_897_07, -487.198_231_355_781_05];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but the recurrence of the power series is inaccurate here (mpmath)
    let actual = [
        xsf::hyp2f1(-1.000_000_000_000_11, 2.000_002, 1e-6, 0.999_999_7),
        xsf::hyp2f1(-49.5, 61.5, 6.5, 0.96),
    ];
    let desired = [-2_000_000.033_618_169_8, -8.458_289_503_668_05e-5];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // and here the power series without that recurrence is more accurate than both (mpmath)
    let actual = [
        xsf::hyp2f1(-3.1, 17.4, 12.3, 0.97),
        xsf::hyp2f1(18.2, -0.9, 12.3, 0.98),
        xsf::hyp2f1(-6.2, 7.7, -2.5, 0.99),
        xsf::hyp2f1(-0.0076, 19.42, 9.7355, 0.8761),
    ];
    let desired = [
        -0.010_434_102_923_645_322,
        -141.307_919_251_094_5,
        45_908_976.313_493_14,
        -30.905_387_588_773_838,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // (but not past the special cases of the Euler transformation, where c - (c - a) rounds to
    // -5 or -12) (mpmath)
    let actual = [
        xsf::hyp2f1(-5.000_000_000_000_101, 106.100_000_000_000_2, 100.1, 0.95),
        xsf::hyp2f1(-11.999_999_999_999_899, 122.100_000_000_000_2, 100.1, 0.9),
    ];
    let desired = [2.580_226_780_133_956_7e-8, 3.231_061_401_736_383e-14];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // (nor if the series doesn't converge: this is -6.33e35, not its partial sum 0.856)
    let actual = xsf::hyp2f1(-1e-13, 3.100_000_000_000_1, 0.1, 0.999_999_999_999_999_9);
    np_assert_allclose!([actual], [-6.326_846_914_851_421e35], rtol = 1e-4);
}

#[test]
fn test_hyp2f1_f64_recurrence_in_a() {
    // the recurrence in `a` can be unstable, while its error estimate `loss` is small (mpmath)
    let actual = [
        xsf::hyp2f1(29.0, -14.5, 0.1, 0.99),
        xsf::hyp2f1(14.0, -14.5, 0.7, 0.98),
        xsf::hyp2f1(-15.8, 25.0, 0.5, 0.98),
    ];
    let desired = [-5.131_087_011_432_369e21, -0.093_683_512_477_367_49, 268_344_349.251_342_5];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but its error estimate shouldn't be too large if a starting value is tiny, so that the less
    // accurate alternatives aren't used instead (Arb)
    let actual = [
        xsf::hyp2f1(-9.5, 10.5, 0.5, 0.95),
        xsf::hyp2f1(10.000_000_001, 9.999_999_999_999_99, 1.5, -3.0),
        xsf::hyp2f1(9.999_999_999, 20.0, -10.5, -3.0),
    ];
    let desired = [-0.897_769_973_095_626_9, 4.289_686_233_768_906e-8, 0.001_570_408_561_128_409_7];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // and its starting values, by the power series with c = 1 - (c - a - b) = -8 - 1e-14 in AMS55
    // #15.3.6, shouldn't stop before the near pole of c (Arb)
    let actual = [
        xsf::hyp2f1(-12.5, 8.999_999_999_999_99, 5.5, 0.998),
        xsf::hyp2f1(-11.5, 8.999_999_999_999_99, 6.5, 0.998),
    ];
    let desired = [2.259_678_955_985_816e-8, -3.542_896_384_613_633e-8];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_x_lt_minus_2() {
    // the 1/x transformation is inaccurate here, but the Pfaff transformations are not (mpmath)
    let actual = [
        xsf::hyp2f1(12.0, 3.1, 99.3, -2.8),
        xsf::hyp2f1(81.0, 6.2, -95.1, -8.4),
        xsf::hyp2f1(-58.5, 0.07, 50.9, -4.3),
    ];
    let desired = [0.416_874_471_029_815_73, 6.875_723_600_176_229e-5, 8_488_942_027_952.015];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but here the starting values of their recurrence need a recurrence as well (mpmath)
    let actual = xsf::hyp2f1(-20.000_000_000_2, -14.47, 2.93, -5.48);
    np_assert_allclose!([actual], [64_520_380_601.436_584], rtol = 1e-12);

    // or need the series without that recurrence (mpmath)
    let actual = [xsf::hyp2f1(-19.6, -4.1, 16.7, -9.9), xsf::hyp2f1(-4.2, -19.6, 14.6, -9.6)];
    let desired = [-1_046_594.825_411_287_6, -4_908_149.039_599_091];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but the 1/x transformation if its error estimate is smaller than that of the Pfaff series,
    // e.g. where it cancels, as xsf does. These were 1.37682208856282, -5.197163034404184e15,
    // 7.081003593530957e-7, -183294.3855496918 and 4100386506.2481966 (Arb).
    let args = [
        (-1.8, -2.5, 5.0, -10.0),
        (9.5, -4.9, -8.7, -1000.0),
        (9.8, 7.4, 0.3, -5.0),
        (-2.6, -6.1, 7.1, -50.0),
        (-7.8, -5.6, 9.0, -100.0),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [
        1.376_822_088_562_981_8,
        -5_197_163_034_410_069.0,
        7.081_003_593_523_789e-7,
        -183_294.385_549_550_82,
        4_100_386_506.245_826_7,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // unlike in xsf, the Pfaff series is kept where both error estimates are too large, or `c - b`
    // is near a pole of its gamma function (Arb)
    let actual = [
        xsf::hyp2f1(
            0.076_598_182_228_762_47,
            7.362_007_163_390_169,
            -3.923_401_672_797_962_4,
            -153.095_573_546_133_72,
        ),
        xsf::hyp2f1(-3.2, 1e-12, -3.000_000_04, -10.0),
    ];
    let desired = [1.258_875_876_921_877e-6, 1.003_134_633_107_021_5];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // or the terms of the 1/x transformation are subnormal: this is 9.37e-324 (Arb)
    assert_eq!(xsf::hyp2f1(1.375, 6.875, 7.875, -1e235).to_bits(), 2);
}

#[test]
fn test_hyp2f1_f64_x_lt_minus_2_b_minus_a_near_int() {
    // unlike in xsf: these were -1.26e5, -2.06e69 and -5.53e3, the last of which is the second
    // function of the 1/x transformation of the one before it (Arb)
    let actual = [
        xsf::hyp2f1(-13.2, -29.200_000_000_001, 5.5, -2.5),
        xsf::hyp2f1(-13.2, -33.21, 50.7, -1000.0),
        xsf::hyp2f1(-33.21, -82.91, -19.01, -0.001),
    ];
    let desired = [-602_275.348_686_216_6, 4.298_609_352_869_04e65, 1.155_996_198_053_486_3];
    np_assert_allclose!(actual, desired, rtol = 5e-13);

    // and this was -2.74e6 (Arb)
    let actual = xsf::hyp2f1(-4.6, -20.600_000_000_1, 50.7, -10.0);
    np_assert_allclose!([actual], [-65.931_816_188_886_59], rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_x_le_minus_9_b_minus_a_int() {
    // unlike in xsf: these lost about log10(-x) digits in the Pfaff transformation, as
    // 1.48876236e-33, 1.32938105e-33, -6.54168264e11, inf, -1297.62704740, 8.4284e-7 and
    // 14742952.50335 (Arb)
    let actual = [
        xsf::hyp2f1(3.0, 3.0, 4.0, -3.7e11),
        xsf::hyp2f1(3.0, 3.0, 6.0, -8.2e11),
        xsf::hyp2f1(-0.5, -1.5, 4.0, -742_448_874.85),
        xsf::hyp2f1(-1.97, -1.97, 5.84, -2.34e13),
        xsf::hyp2f1(1.25, -0.75, -3.45, -58_000.0),
        xsf::hyp2f1(0.375, 3.375, -1.5, -1e15),
        xsf::hyp2f1(-2.25, 0.75, 7.5, -1e4),
    ];
    let desired = [
        1.488_762_887_248_107_1e-33,
        1.329_379_117_616_510_5e-33,
        -654_168_241_771.272_2,
        2.019_700_022_161_330_8e25,
        -1_297.627_047_387_899,
        8.408_236_071_439_006e-7,
        14_742_952.503_354_19,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_x_lt_minus_2_pfaff_overflow() {
    // unlike in xsf: these were NaN for a < b, as one of the Pfaff series overflows there, and
    // the 1/x transformation is NaN (Arb)
    let args = [
        (1e-20, 3000.5, 1e7, -100.0),
        (3.1e-111, 344_352.06, 1.77e10, -621.8),
        (0.1, 3000.5, 3600.6, -100.0),
        (1.5, 344_352.06, 413_222.472, -3.0),
    ];
    let desired = [1.0, 1.0, 0.641_802_447_619_336_4, 0.152_720_780_375_729_64];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        // also with a and b swapped, which was already accurate
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-14);
    }
}

#[test]
fn test_hyp2f1_f64_tiny_first_term() {
    // unlike in xsf: these were 1, as the series stopped after a tiny first term (Arb)
    let actual = [
        xsf::hyp2f1(1e-20, -300.5, 300.5, -10.0),
        xsf::hyp2f1(1e-300, -500.5, 0.7, -10.0),
        xsf::hyp2f1(1e-20, 150.5, 100.7, 0.9),
        xsf::hyp2f1(-1e-20, 50.5, 0.7, 0.8),
    ];
    let desired = [
        6.528_217_172_881_034e123,
        2.950_367_578_022_261_5e219,
        1.000_000_046_617_229,
        -63_717_533_596_884.64,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // Arb: 1.88e1141
    assert_eq!(xsf::hyp2f1(1e-300, -3000.5, 3000.5, -10.0), f64::INFINITY);
    // the terms grow again here, but stay negligible
    assert_eq!(xsf::hyp2f1(1e-300, 4.000_000_000_000_008e30, 1e30, 0.25), 1.0);
}

#[test]
fn test_hyp2f1_f64_euler_tiny_a() {
    // unlike in xsf: these were 1, 1, inf, 1 and -1.24e43, as `a` or `b` was lost in `c - a` or
    // `c - b` (Arb)
    let actual = [
        xsf::hyp2f1(1e-20, 300.5, 0.7, 0.5),
        xsf::hyp2f1(300.5, 1e-20, 0.7, 0.5),
        xsf::hyp2f1(1e-20, -800.5, 0.1, -1.5),
        xsf::hyp2f1(1e-20, 40.5, 20.5, 0.99),
        xsf::hyp2f1(1e-20, -370.5, 53.5, -1.0),
    ];
    let desired = [
        6.921_479_145_123_81e68,
        6.921_479_145_123_81e68,
        1.824_849_365_063_929_7e299,
        63_529_516.407_577_99,
        2.948_782_385_828_795e37,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // and these are 1 to rounding
    let actual = [
        xsf::hyp2f1(1e-20, 3.0, 0.7, -1.0),
        xsf::hyp2f1(1e-20, 300.5, 0.7, -0.999_999_999_999_99),
    ];
    np_assert_allclose!(actual, [1.0, 1.0], rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_tiny_negative_a_or_b() {
    // unlike in xsf, where a tiny negative `a` or `b` counts as the integer 0: these were inf, 1,
    // 1, 1, 1, 1, 1, -4.27e6, 1, 1 and -inf (Arb)
    let actual = [
        xsf::hyp2f1(26.5, -3e-21, 0.001, -2e4),
        xsf::hyp2f1(-120.5, -2.5e-28, -9.94, -0.5005),
        xsf::hyp2f1(-1e-20, -30.5, 2.5, -100.0),
        xsf::hyp2f1(-1e-20, -50.5, 3.5, -1.5),
        xsf::hyp2f1(-1e-20, 150.5, 0.7, 0.99),
        // also for small parameters, and if the other is a negative integer
        xsf::hyp2f1(-4e-23, -9.75, -5.5, -5800.0),
        xsf::hyp2f1(-7.0, -1e-27, 9.5, -48000.0),
        xsf::hyp2f1(8.25, -3e-17, 1.25, 0.999_99),
        xsf::hyp2f1(-3e-23, 7.0, 1.625, 0.999_98),
        xsf::hyp2f1(-0.25, -1e-25, -9.25, 0.9985),
        xsf::hyp2f1(-9.9e-14, 5.25, -0.0112, -3.75),
    ];
    let desired = [
        1.0,
        -10_864.775_782_791_036,
        -3.386_320_822_599_271e37,
        0.999_973_883_705_280_6,
        -1.562_482_690_977_387e278,
        -5.396_467_458_796_281e16,
        -8.719_200_573_711_808,
        -2.337_318_708_006_997_3e17,
        -27.935_689_206_440_796,
        1.146_412_218_027_53,
        0.999_999_999_991_500_6,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // and these are 1 to rounding, also if the 1/x transformation would lose `a` in `c - (c - a)`
    let actual = [
        xsf::hyp2f1(-5e-324, -9.75, -5.5, -5800.0),
        xsf::hyp2f1(-1e-300, -120.5, -9.94, -0.5005),
        xsf::hyp2f1(-4.8e-28, -1.1, 117.5, -1.63),
    ];
    np_assert_allclose!(actual, [1.0, 1.0, 1.0], rtol = 1e-15);
    // also where the error estimates of the paths for polynomials say that no digit is correct
    // (Arb)
    let actual = [
        xsf::hyp2f1(-3.53e-26, -1592.2, -17.239, 0.04516),
        xsf::hyp2f1(-3.25e-18, 491.04, -0.010_315, -0.328),
    ];
    np_assert_allclose!(actual, [1.0, 1.0], rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_c_lt_0_x_lt_0() {
    // the Euler transformation cancels here, but the Pfaff ones don't (mpmath)
    let actual = [
        xsf::hyp2f1(6.5, 8.9, -9.8, -0.08),
        xsf::hyp2f1(9.9, 4.5, -9.9, -0.01),
        xsf::hyp2f1(16.5, -19.2, -19.7, -0.3),
        xsf::hyp2f1(14.5, -15.5, -19.6, -0.63),
    ];
    let desired = [
        1.707_641_797_398_377,
        1.046_566_861_351_764,
        1.874_594_164_146_283_7e-4,
        -0.143_812_194_757_134_78,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_a_eq_minus_b_neg_int() {
    // Arb
    let args = [
        (13.0, -13.0, 0.5, 0.9),
        (16.0, -16.0, 0.5, 0.95),
        (19.0, -19.0, 0.5, 0.99),
        (22.0, -22.0, -3.5, 0.95),
        (195.0, -195.0, 1.0, 0.9),
        (195.0, -195.0, -200.5, 0.01),
    ];
    let desired = [
        0.489_514_398_515_199_2,
        0.595_230_641_982_666_2,
        0.787_059_088_545_409,
        -6_540_705.712_015_875,
        -0.019_154_323_187_421_253,
        6.722_964_015_601_578,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        // also with a and b swapped, which was already accurate
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 2e-14);
    }

    // and beyond the recurrence limit of the power series, in both orders by `hyp2f1_recur`
    let actual = [
        xsf::hyp2f1(10_001.0, -10_001.0, 9_999.5, 1e-12),
        xsf::hyp2f1(-10_001.0, 10_001.0, 9_999.5, 1e-12),
    ];
    np_assert_allclose!(actual, [0.999_999_989_997_499_8; 2], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_c_lt_0() {
    // the terms of the power series are tiny before k = -c, and huge after it (mpmath)
    let actual = [
        xsf::hyp2f1(0.01, 0.3, -181.2, 0.6),
        xsf::hyp2f1(0.1, 0.2, -1000.5, 0.55),
        xsf::hyp2f1(0.1, 0.2, -12_000.5, 0.02),
        xsf::hyp2f1(-7001.5, 12.2, -7000.5, 0.5),
    ];
    let desired = [
        3.196_698_793_921_328_4e28,
        -6.233_446_442_910_298e83,
        0.999_999_966_668_092_2,
        4_713.283_633_348_583,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // the ratio of its terms underflows here, but not for the backward recurrence in c (the direct
    // sum of the power series in high precision)
    let actual =
        [xsf::hyp2f1(1e-200, 1e-200, -1000.5, 0.8), xsf::hyp2f1(1e-160, 1e-160, -1000.5, 0.8)];
    let desired = [-5.767_096_442_064_845e199, -5.767_096_442_064_845e279];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_neg_int_c_minus_a_or_c_minus_b() {
    // unlike in xsf: by `hyp2f1_recur` if c - a or c - b is only near a negative integer and
    // -1 < c - a - b < 0, where the recurrence in `a` is estimated as accurate, but isn't. These
    // are 1.18e11, 24.58, 1.86e7, 61.73 and -1.73e13 in xsf (mpmath, Wolfram)
    let actual = [
        xsf::hyp2f1(-40.5, 40.9, -0.1, 0.8),
        xsf::hyp2f1(-25.25, 25.7, -0.3, 0.75),
        xsf::hyp2f1(-30.5, 31.2, 0.2, 0.85),
        xsf::hyp2f1(35.6, -35.1, -0.4, 0.7),
        xsf::hyp2f1(-45.5, 45.99, -0.01, 0.8),
    ];
    let desired = [
        116.439_612_499_602_16,
        24.543_627_215_364_04,
        -15.110_473_882_473_225,
        74.102_598_993_798_93,
        -4.645_596_140_460_149,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // also for x near 1 (mpmath, Wolfram)
    let actual = xsf::hyp2f1(-59.25, 70.1, 10.1, 0.999);
    np_assert_allclose!([actual], [-1.622_052_617_419_225_7e-11], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_parameters() {
    // unlike in xsf: the backward recurrence in c for large parameters (references: the direct sum
    // of the power series of the Pfaff transformation in high precision)
    let actual = [
        xsf::hyp2f1(35.17, 102.3, 0.0054, -9.74),
        xsf::hyp2f1(658.99, 0.0047, 51.38, -2.72),
        xsf::hyp2f1(-0.5, 692.8, 83.35, -3.98),
        xsf::hyp2f1(-21.2, -10.2, -1e-8, -0.46),
    ];
    let desired = [
        1.955_342_940_723_380_9e-60,
        0.983_273_621_397_356_7,
        5.845_197_863_408_547,
        -53_404_320_641.112_48,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // x << -1: the 1/x transformation (references: its two series summed in high precision)
    let actual = [
        xsf::hyp2f1(1.35, 87.4, -32.1, -8478.0),
        xsf::hyp2f1(-41.6, -16.3, 720.1, -204.7),
        xsf::hyp2f1(387.2, -2.01, 1.7, -6264.9),
        xsf::hyp2f1(849.1, -1.54, 55.0, -236.7),
    ];
    let desired = [
        -4.397_527_737_978_469e-6,
        -1.151_698_277_128_985_8e17,
        1_471_625_105_699.610_7,
        304_675.781_003_058_8,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // also for a tiny c, whose `ln Gamma(c)` has no second-order error term (Arb)
    let actual = [xsf::hyp2f1(20.5, 10.25, 1e-200, -1e6)];
    np_assert_allclose!(actual, [-5.514_252_855_663_787e132], rtol = 1e-12);

    // long recurrences, whose solutions oscillate in c (references: Arb), e.g. 1.6e617 overflows
    let actual = [
        xsf::hyp2f1(-522.3, -123.4, 533.9, -0.9),
        xsf::hyp2f1(-552.3, -953.1, -55.5, -4.13),
        xsf::hyp2f1(-679.9, -678.6, -28.7, -4.79),
    ];
    let desired = [5.957_219_702_847_165e-52, f64::INFINITY, f64::NEG_INFINITY];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_parameters_symmetric() {
    // unlike in xsf: the backward recurrence in c, for b near a negative integer. With a and b
    // swapped, these were 7.3575092e-14, NaN and -6.95925e-24 (Arb)
    let args = [
        (36.0, -104.0 - 9e-7, 0.0074, 0.66),
        (61.0, -193.0 - 2.8e-7, 0.000_107, 0.64),
        (59.0, -173.0 - 7.6e-9, 0.0274, 0.64),
    ];
    let desired =
        [7.357_509_219_073_461e-14, -1.626_123_558_941_129_7e-26, -6.960_243_881_349_578e-24];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-12);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }

    // and in both orders, but less accurate, where xsf gives -3.75e-33 and 2.99e12 (Arb)
    let args = [(-206.0 + 5.2e-7, 47.0, 3.6e-7, 0.72), (142.0, -372.0 + 9.1e-5, 2.1e-6, 0.86)];
    let desired = [-6.514_961_820_519_327e-38, 7.812_813_451_522_689e-90];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-10);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }
}

#[test]
fn test_hyp2f1_f64_large_integer_a_near_integer_b() {
    // unlike in xsf: the 1 - x transformation, for an integer a > 0, where one of its terms is a
    // polynomial. These were 1.7143275e-23, NaN, -1.87323e-22, NaN and NaN, where xsf gives
    // 3.28e21, 7.58e22, 8.07e91, -2.43e20 and 1.04e33 (Arb).
    let args = [
        (94.0, -172.000_005_2, 0.000_69, 0.796),
        (144.0, -307.0 - 3.2e-8, 1.07, 0.83),
        (107.0, -144.0 - 7.3e-6, 0.145, 0.945),
        (85.0, -171.999_998, 0.0369, 0.8847),
        (100.0, -185.000_000_01, 0.0018, 0.866),
    ];
    let desired = [
        1.714_327_734_183_920_2e-23,
        -2.798_543_910_807_295e-65,
        -1.871_761_033_753_452_6e-22,
        -4.884_956_947_532_335e-39,
        1.330_815_061_015_127_8e-34,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-13);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }

    // but not at the cost of the other forms, whose bounds can be smaller where the two terms of
    // the transformation cancel, e.g. for x < -1 (Arb; xsf gives NaN)
    let (a, b, c, x) = (194.001, 350.0, 55.065, -5.0);
    let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
    np_assert_allclose!(actual, [-9.278_721_967_709_765e-268; 2], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_parameters_near_log_case() {
    // unlike in xsf: the 1/x and 1 - x transformations for b - a or c - a - b within 1e-3 of an
    // integer, where their two terms cancel. These were NaN, NaN, NaN, 4688.54, 2.89638e-42,
    // -2.16715878e-48, -1.88005587e-30 and NaN, where xsf gives NaN, NaN, NaN, 4.84e3, 1.03e-31,
    // 5.88e-38, -1.88005587e-30 and NaN (Arb).
    let args = [
        (1.07, 3.07 + 2.5e-13, 251.94, -524_766.9),
        (79.25, 85.25 + 4e-12, -27.19, -411.2),
        (-312.005, 9.5, -302.505 - 6e-14, 0.999_97),
        (47.25, 1.5, 48.75 - 2.5e-11, 0.999_999_98),
        (53.75, 27.75 - 3.5e-11, 0.22, -10.4),
        (-68.75, 107.26, 50.51 + 3.6e-11, 0.9976),
        (115.25, -73.5, 40.75 + 1.8e-6, 0.9998),
        (46.25, 76.25 + 5e-5, -300.6, -145.7),
    ];
    let desired = [
        0.000_131_036_002_073_166_09,
        3.918_853_590_964_261e-188,
        -1.038_538_785_178_233_9e19,
        4_688.172_838_647_12,
        2.896_395_178_110_93e-42,
        -2.167_158_642_601_661_3e-48,
        -1.880_055_697_263_691_3e-30,
        5.591_540_566_735_008e-62,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-13);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }
}

#[test]
fn test_hyp2f1_f64_small_parameters_near_log_case() {
    // unlike in xsf: also for small parameters, with c - a - b (for x >= 0.9) or b - a (for
    // x < -2) within 1e-3 of an integer. The first four were NaN, where xsf gives 2.27784513,
    // 4247.41591, -1673017043.2 and 12.27095. The others were -9.186935e-4, 1.31533e19, -883812,
    // 1.3971776e-3 and -64173143.8372, as in xsf except for its 1.31530e19 and -64179200 (Arb).
    let args = [
        (0.5, 0.25, 0.750_000_001, 0.999),
        (0.5, 2.5, 1.999_999_999_99, 0.9999),
        (0.5, -1.5, -0.999_999_999, 0.999),
        (0.5, 7.5, 7.999_999_999_99, 0.9999),
        (-8.75, 8.875, 2.125_000_01, 0.985),
        (-7.5, -7.500_000_000_01, 5.5, -1000.0),
        (-2.5, -2.500_000_000_01, 5.5, -1000.0),
        (1.5, 1.499_999_999_99, 5.5, -1000.0),
        (-9.25, -9.249_999_999_99, 0.75, -6.5),
    ];
    let desired = [
        2.277_844_918_975_149_3,
        4_247.415_771_798_882,
        -1_673_017_343.073_581,
        12.271_558_655_659_248,
        -9.192_870_681_661_742e-4,
        1.314_998_615_206_617_7e19,
        -883_900.539_293_860_8,
        1.397_157_529_639_970_1e-3,
        -64_173_143.838_225_69,
    ];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-14);
    }
}

#[test]
fn test_hyp2f1_f64_near_log_case_subnormal() {
    // b - a or c - a - b within a subnormal of an integer, where that distance underflows in the
    // combined form of `hyp2f1_recur`. Unlike in xsf, the first two aren't NaN (Arb).
    let args = [
        (5.4e-323, 7.0, -210.3, -12_265.5),
        (2.0, 3.5e-323, 275.9, -114_255.8),
        (4.535e-320, 155.75, 158.75, 0.999_987_6),
        (1.859_931_117e-315, -79.875, -78.875, 0.999_998_1),
    ];
    for (a, b, c, x) in args {
        np_assert_allclose!(
            [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)],
            [1.0; 2],
            rtol = 1e-12
        );
    }
}

#[test]
fn test_hyp2f1_f64_euler_less_accurate() {
    // unlike in xsf: the Euler transformation for c - a - b <= -1, where the error bound of
    // `hyp2f1_recur` shows its result to be less accurate, or it is NaN. These were -1.74008243e44,
    // 9.2789205253e35, 1.41743332e52, 1.02457247e39, NaN and NaN, where xsf gives 2.11e75,
    // 9.27905840e35, 1.60e95, -3.41e97, NaN and NaN (Arb; the last is -2.8e-380).
    let actual = [
        xsf::hyp2f1(149.0, -58.0 - 1.1e-7, 7.3e-7, 0.87),
        xsf::hyp2f1(189.0, -34.0 - 1.1e-8, 1.4e-7, 0.603),
        xsf::hyp2f1(185.0, -72.0 - 1.8e-9, 7.1e-7, 0.862),
        xsf::hyp2f1(136.0, -67.0 - 4.1e-7, 9.55e-5, 0.925),
        xsf::hyp2f1(9.25, 18.25 - 7.5e-8, -200.97, -90.43),
        xsf::hyp2f1(241.75, 220.75 + 2.65e-9, -3.67, -46.87),
    ];
    let desired = [
        -1.740_082_505_618_863_8e44,
        9.278_920_509_129_896e35,
        1.417_433_100_506_847_6e52,
        1.024_572_396_591_105_2e39,
        -3.130_221_533_372_943e-5,
        -0.0,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // with its power factor split in two, and for a subnormal `a`. These were -1.7311832e-211 and
    // NaN, where xsf gives 0.0 and NaN (Arb).
    let actual = [
        xsf::hyp2f1(98.0, 124.0 - 1.1e-13, 3.06, -71.03),
        xsf::hyp2f1(1e-323, 1.0, -159.8, -104_073.0),
    ];
    np_assert_allclose!(actual, [-1.725_017_797_892_166_1e-211, 1.0], rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_pfaff_recip_less_accurate() {
    // unlike in xsf: the Pfaff transformation for x < -1 and the 1/x one for x < -2, with large a
    // and b of the same sign, where the error bound of `hyp2f1_recur` shows their result to be
    // less accurate, or it is not finite. These were -9.13e168, 1.16e146 and -inf, as in xsf
    // (Arb).
    let actual = [
        xsf::hyp2f1(-129.666, -113.666, 56.25, -52.03),
        xsf::hyp2f1(-122.125, -86.125 + 6.9e-5, 78.26, -76.27),
        xsf::hyp2f1(-176.9375, -130.9375 + 3.3e-6, -21.0 + 5.2e-8, -35.47),
    ];
    let desired =
        [4.791_283_530_597_28e162, -1.531_524_570_886_761_4e144, 2.979_197_152_672_975e274];
    np_assert_allclose!(actual, desired, rtol = 1e-9);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_error() {
    // unlike in xsf: the psi function expansion for x > 0.9 and an integer c - a - b has an error
    // estimate, so this is replaced by the result of `hyp2f1_recur`. These were 1.36200088831e-8
    // and 4.8456802331e-31, as in xsf (Arb).
    let actual = [
        xsf::hyp2f1(
            -61.245_538_520_640_43,
            52.213_688_621_304_07,
            0.968_150_100_663_642_4,
            0.957_677_717_016_758_5,
        ),
        xsf::hyp2f1(
            -29.096_463_013_348_03,
            81.222_637_166_243_06,
            71.126_174_152_895_03,
            0.912_150_259_055_036_1,
        ),
    ];
    let desired = [1.362_000_886_533_284_2e-8, 4.845_680_233_889_792e-31];
    np_assert_allclose!(actual, desired, rtol = 1e-11);
    // and for the Pfaff transformation for x < -1, where these were -1.94e102, -1.95e210 and 0 (as
    // its factors underflow), while xsf gives -1.96e102, inf and 0 (Arb)
    let actual = [
        xsf::hyp2f1(-138.75, -85.75, 31.625, -13.86),
        xsf::hyp2f1(-127.75, -164.75 - 2.8e-14, 1.0 + 8.9e-9, -25.89),
        xsf::hyp2f1(-161.75, -108.75, 63.125, -14.8),
    ];
    let desired =
        [1.476_521_363_124_945_4e100, 3.311_772_741_737_614_5e207, 5.121_708_481_485_257e107];
    np_assert_allclose!(actual, desired, rtol = 1e-9);
}

#[test]
fn test_hyp2f1_f64_huge_c() {
    // unlike in xsf: these were NaN, and the first and fourth hung (Arb, or the limit 1)
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        tx.send([
            xsf::hyp2f1(1.3, 6.5, f64::INFINITY, 0.999),
            xsf::hyp2f1(1.3, 6.5, f64::INFINITY, 1.0),
            xsf::hyp2f1(1.3, 6.5, f64::INFINITY, -3.0),
            xsf::hyp2f1(1.3, 6.5, 1e300, 0.999),
            xsf::hyp2f1(100_000_000.5, 1_000_000_000.5, 1e17, 0.95),
        ])
        .unwrap();
    });
    let actual = rx.recv_timeout(std::time::Duration::from_secs(10)).unwrap();
    let desired = [1.0, 1.0, 1.0, 1.0, 2.585_709_685_660_995_7];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
}

#[test]
fn test_hyp2f1_f64_huge_c_at_1() {
    // unlike in xsf: these are NaN in xsf, as `Gamma(c)` overflows in `F(a, b; c; 1)`
    // (DLMF 15.4.20; Arb)
    let actual = [
        xsf::hyp2f1(1.3, 6.5, 1e17, 1.0),
        xsf::hyp2f1(1.3, 6.5, 1e300, 1.0),
        xsf::hyp2f1(200.0, 0.5, 300.0, 1.0),
        xsf::hyp2f1(17.5, 250.5, 312.0, 1.0),
        xsf::hyp2f1(6.9e16, -469.5, 1.65e17, 1.0),
        // a term with an overflowing denominator, and a slowly decaying series (d = 2.7e-14,
        // 0 in xsf)
        xsf::hyp2f1(6e153, 6e153, 1e308, 1.0),
        xsf::hyp2f1(1e-15, 179.999_999_999_999_97, 180.0, 1.0),
        // a sum close to `f64::MAX`, an overflowing `a b`, and an underflowing `a / c`
        xsf::hyp2f1(1_946.607_092_231_488_2, 5_848_583_033.383_596, 19_141_581_561.342_18, 1.0),
        xsf::hyp2f1(2e154, 2e154, 1e308, 1.0),
        xsf::hyp2f1(1e-16, 2e307, 1e308, 1.0),
    ];
    let desired = [
        1.0,
        1.0,
        1.736_405_723_697_015_8,
        24_384_542_401_255.008,
        3.695_150_263_060_27e-111,
        1.433_329_414_560_340_4,
        1.036_467_456_652_845,
        1.797_693_134_862_289e308,
        54.598_150_033_144_25,
        1.0,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // overflows (+inf for both in xsf)
    let c = 5.967_056_889_428_375e-309;
    assert_eq!(xsf::hyp2f1(-1.25, 1.125, c, 1.0), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(-1.25, 1.125, -c, 1.0), f64::NEG_INFINITY);
    // symmetric in a and b, also for equal |a| and |b|
    let (a, c) = (190.513_976_388_490_14, 19_977.707_966_310_958);
    assert_eq!(xsf::hyp2f1(a, -a, c, 1.0), xsf::hyp2f1(-a, a, c, 1.0));
}

#[test]
fn test_hyp2f1_f64_psi_expansion_overflow() {
    // unlike in xsf: these were NaN, as the Gamma functions in the psi function expansion for
    // x > 0.9 and integer c - a - b over- or underflow (Arb)
    let actual = [
        xsf::hyp2f1(1.3, 6.5, 107.8, 0.95),
        xsf::hyp2f1(-1.5, -1.5, 110.0, 0.99),
        xsf::hyp2f1(0.25, 20.5, 283.75, 0.95),
        xsf::hyp2f1(20.5, 30.5, 10_051.0, 0.95),
    ];
    let desired = [
        1.080_572_289_607_784_4,
        1.020_272_592_726_024_6,
        1.017_973_915_497_389_6,
        1.061_033_031_294_106_8,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-15);

    // but a correct overflow (to 3.96e320) is kept
    assert_eq!(xsf::hyp2f1(6.5, 6.5, 1e-300, 0.95), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_ams55_15_3_6_gamma_range() {
    // unlike in xsf: these were 1, NaN, NaN and NaN, as a factor of AMS55 #15.3.6 left the normal
    // range, e.g. `Gamma(-d) / (Gamma(a) Gamma(b))` for a tiny `b`, or `Gamma(c)` (Arb)
    let actual = [
        xsf::hyp2f1(211.205, 4.1e-92, 141.983, 0.9992),
        xsf::hyp2f1(90.5, 1e-180, 1.25, 0.99999),
        xsf::hyp2f1(26.5, -508.5, 590.5, 0.999_999),
        xsf::hyp2f1(-1.5, -693.5, -176.25, 0.999_999_97),
    ];
    let desired = [
        3.803_424_824_536_675e64,
        5.865_541_318_634_591e263,
        5.295_559_383_862_064_7e-8,
        5.041_937_663_630_442,
    ];
    np_assert_allclose!(actual, desired, rtol = 5e-13);

    // but not at a pole of `Gamma(c - a - b)` (this gave -inf), where `hyp2f1_recur` overflows
    // to inf, as the true value does (Arb)
    assert_eq!(xsf::hyp2f1(1e5, 7.5, 21.5, 0.95), f64::INFINITY);
    // and here an inner series fails (this gave -3.96e59), but `hyp2f1_recur` does not (Arb)
    let actual = [xsf::hyp2f1(862.559, -310.305, 906.676, 0.999_993_238)];
    np_assert_allclose!(actual, [2.416_601_765_882_886e-243], rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_near_int() {
    // unlike in xsf: c - a - b (or b - a for x < -1) within 1e-13 of an integer cancelled in
    // AMS55 #15.3.6, which gave 573, -1.03, -1.33e14, 5.23345 and 3.49e-4 (Arb)
    let actual = [
        xsf::hyp2f1(7.5, 2.5 - 1e-14, 10.0, 0.999),
        xsf::hyp2f1(-12.5, 23.0 - 1e-14, 5.5, 0.95),
        xsf::hyp2f1(-2.5, 2.5, 1e-14, 0.99),
        xsf::hyp2f1(1.5, 4.5 - 1e-14, 8.0, 0.999),
        xsf::hyp2f1(3.0, 3.0 - 1e-14, 50.7, -1000.0),
    ];
    let desired = [
        465.329_511_494_182_9,
        0.066_364_450_275_830_67,
        -113_060_171_546_044.02,
        5.233_384_708_085_254,
        7.920_150_150_695_394e-5,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but the power series is tried first (the psi function expansion gives 6.746770e-4)
    let actual = [xsf::hyp2f1(7.5, 10.5 - 1e-13, 50.7, -10.0)];
    np_assert_allclose!(actual, [6.746_767_352_078_354e-4], rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_int() {
    // unlike in xsf: the power series is also tried first if c - a - b is an integer, here after
    // Pfaff, where the psi function expansion gave 6.7467671e-4, 2.0396740e-3, 0.1223477772 and
    // 0.358811271431 (Arb)
    let actual = [
        xsf::hyp2f1(7.5, 10.5, 50.7, -10.0),
        xsf::hyp2f1(7.5, 8.5, 50.7, -10.0),
        xsf::hyp2f1(3.0, 6.0, 50.7, -10.0),
        xsf::hyp2f1(1.5, 5.5, 50.7, -10.0),
    ];
    let desired = [
        6.746_767_352_078e-4,
        2.039_673_789_485_802e-3,
        0.122_347_777_162_347_68,
        0.358_811_271_434_322,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_huge_c_minus_1_lt_c_minus_a_minus_b_lt_0() {
    // unlike in xsf: the recurrence on c overflowed the stack for |c| >= 2^53, or was wrong
    // (6.78e137 for the last) (Arb)
    let actual = [
        xsf::hyp2f1(1e17 + 16.0, -15.5, 1e17, 0.5),
        xsf::hyp2f1(1e17 + 16.0, -15.5, 1e17, 0.95),
        xsf::hyp2f1(1e17 - 96.0, 96.5, 1e17, -0.1),
        xsf::hyp2f1(1e17 - 96.0, 96.5, 1e17, 0.7),
        xsf::hyp2f1(1e16 - 1000.0, 1000.9, 1e16, -0.3),
    ];
    let desired = [
        2.157_918_643_757_769e-5,
        6.823_937_919_615_83e-21,
        1.012_991_688_969_088_8e-4,
        2.869_451_704_379_998e50,
        8.997_019_796_944_917e-115,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    assert!(xsf::hyp2f1(1e17, 0.5, 1e17, f64::NAN).is_nan());
}

#[test]
fn test_hyp2f1_f64_tiny_c_or_c_minus_a() {
    // unlike in xsf: tiny c and c - a were treated as 0, which gave 2, inf, NaN, 1.0 and
    // 1.0000000000000023e-5 (Arb)
    let actual = [
        xsf::hyp2f1(-1e-30, 1.0, 1e-17, 0.5),
        xsf::hyp2f1(5.2, 12.5, -1e-17, 0.6),
        xsf::hyp2f1(12.5, 0.3, 1e-15, 0.2),
        xsf::hyp2f1(1e-14, -1.5, 1e-17, -1.5),
        xsf::hyp2f1(-0.4, -2.5, -0.400_000_000_000_09, 0.99),
    ];
    let desired = [
        0.999_999_999_999_9,
        -5.650_224_541_964_126e27,
        6.795_217_861_653_541e15,
        2_953.847_075_210_480_4,
        1.000_000_042_269_799_7e-5,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
}

#[test]
fn test_hyp2f1_f64_tiny_negative_c_large_parameters() {
    // unlike in xsf, where a tiny negative c counts as the integer 0 (Arb)
    let actual = [
        xsf::hyp2f1(28.5, 32.25, -6.5e-267, -1.85),
        xsf::hyp2f1(-44.5, 46.25, -3.5e-103, 0.6),
        xsf::hyp2f1(41.25, 47.5, -9e-292, -8.75),
        xsf::hyp2f1(-17.5, 37.25, -6e-271, -9.75),
        xsf::hyp2f1(44.75, -32.5, -3e-234, 0.83),
        xsf::hyp2f1(-33.5, 23.25, -8e-16, 0.79),
        xsf::hyp2f1(-40.5, -47.25, -4.5e-193, -1.8),
        // and where the inner function of the Euler transformation overflows
        xsf::hyp2f1(60.0, 2.5e-6, -2e-266, -7.7e4),
        xsf::hyp2f1(43.0, 9.25, -2.5e-285, -600.0),
        xsf::hyp2f1(-15.5, 38.0, -3.5e-257, -23.5),
    ];
    let desired = [
        -4.619_620_997_681_311e252,
        -2.351_118_090_908_530_3e103,
        7.447_544_914_471_166e247,
        -6.278_839_855_668_874e302,
        -6.149_027_963_048_359e238,
        302_521_514_418.673_6,
        -3.702_964_233_045_192e212,
        1.249_950_267_344_792e260,
        -2.981_261_708_812_470_7e249,
        -7.762_342_418_664_52e291,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_euler_pow_underflow() {
    // unlike in xsf: (1 - x)^(c - a - b) under- or overflowed in the Euler transformation, which
    // gave 0, 0, 0, 0 and inf (Arb)
    let actual = [
        xsf::hyp2f1(9.7, 23.71, -20.7, -1e6),
        xsf::hyp2f1(0.5, 3.7, 0.7, -1e100),
        xsf::hyp2f1(0.5, 9.7, -5.5, -1e30),
        xsf::hyp2f1(2.5, 50.7, 20.37, -1e10),
        xsf::hyp2f1(3.0, 20.37, 4.000_000_000_000_001, 0.999_999_999_999_999_9),
    ];
    let desired = [
        2.339_902_740_854_614e-57,
        1.643_323_996_265_172_4e-51,
        2.444_777_663_146_785_3e-287,
        8.912_118_823_854_448e-27,
        1.699_790_832_401_936_6e308,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_euler_both_out_of_range() {
    // Arb; both factors of the Euler transformation leave the normal range here
    let actual = [
        xsf::hyp2f1(14.5, 1e-290, -368.7, -135.0),
        xsf::hyp2f1(97.0, 1e-240, -16.1, -7.7e5),
        xsf::hyp2f1(4.1, -6.5e-8, -68.9, -1.25e7),
        xsf::hyp2f1(-1.5e-9, 932.0, 8.5, -56.8),
    ];
    let desired = [1.0, 1.0, 1.000_001_498_528_907_5, 1.000_000_013_194_575_4];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // but not where x / (x - 1) rounds to 1 (Arb: 8.12e10)
    assert!(xsf::hyp2f1(1e-20, -0.1, -2.5, -1e300).is_nan());
}

#[test]
fn test_hyp2f1_f64_euler_inner_overflow() {
    // Arb; the inner function of the Euler transformation overflows here, but its product with
    // (1 - x)^(c - a - b) doesn't
    let actual = [
        xsf::hyp2f1(21.0, 43.0, 1e-276, -24.0),
        xsf::hyp2f1(0.5, 20.0, 1e-233, -7.5e4),
        xsf::hyp2f1(-15.5, 48.0, 5e-261, -78.0),
        xsf::hyp2f1(34.0, -1.25, 2e-278, -1.7e4),
        // also where (1 - x)^(c - a - b) underflows
        xsf::hyp2f1(47.0, 43.0, 2e-237, -3000.0),
        // also where the rounding of 1 - x matters
        xsf::hyp2f1(1e14, 10.0, 1e-300, -1e-12),
        // also where that product is small, and the polynomial that the inner function is
        // rounded to lacks the dominant part of the result
        xsf::hyp2f1(22.0, 4098.5, 1e-189, -8e4),
    ];
    let desired = [
        -3.248_055_641_674_947_7e228,
        -2.347_637_291_530_230_8e229,
        1.281_349_226_923_448_9e305,
        8.831_616_558_293_24e284,
        -1.238_995_047_182_927e80,
        3.787_960_193_849_875e270,
        5.393_827_268_558_922e-167,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_nan_if_lost() {
    // unlike in xsf: NaN if both error estimates say that no digit is correct (no input is known to
    // reach this, as `hyp2f1_recur` bounds the errors of those below), where xsf gives 6.75e10 and
    // -2.81e238 (Arb)
    let actual = [
        xsf::hyp2f1(-7.000_000_000_000_021, -208.513_348_904_413_08, 0.5, 0.471_817_095_982_536_96),
        xsf::hyp2f1(
            -17.000_000_000_000_007,
            738.317_738_730_952_2,
            -128.5,
            -0.362_556_171_941_935_9,
        ),
    ];
    let desired = [119_399_279_643.698_67, -210.050_825_314_093_3];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // and these, which were 7.7e149, 2.4e7, -5.6e85 and 7.5e-4 (Arb)
    let actual = [
        xsf::hyp2f1(6.08e-26, -908.762, -19.967, 0.5),
        xsf::hyp2f1(-50.7, 0.1, 344.3, -1.0),
        xsf::hyp2f1(1e-14, -1500.5, 0.7, 0.2),
        xsf::hyp2f1(90.0, -90.0, 180.7, 0.8),
    ];
    let desired = [1.0, 1.016_019_727_590_008, 0.999_999_999_999_930_7, 3.771_813_135_532_291e-17];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // a series that doesn't converge is estimated as 1, but this is accurate (Arb)
    let actual = [xsf::hyp2f1(-8.999_999_999_999_99, 4.0, 2.5, -1.0)];
    np_assert_allclose!(actual, [2_246.910_399_198_305_6], rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_nan_if_err_above_1() {
    // NaN if the error estimate with the recurrences in `a` is above 1 (no input is known to reach
    // this, as `hyp2f1_recur` bounds the errors of those below), but not for a polynomial, here
    // T_10 near a root, which is accurate in absolute terms (Arb)
    let actual = xsf::eval_chebyshev_t(10.0, -core::f64::consts::FRAC_1_SQRT_2);
    assert!((actual - 6.835_808_657_661_923e-16).abs() < 1e-15);
    // and `hyp2f1_recur` bounds its error for these (Arb)
    let actual = [
        xsf::hyp2f1(29.0, -29.1, 0.18, 0.86),
        xsf::hyp2f1(-38.9, 34.5, -0.07, 0.885),
        xsf::hyp2f1(-0.15, -45.9, 51.0, -0.94),
    ];
    let desired = [9.595_688_178_599_744, 0.802_925_766_046_597_9, 0.783_576_573_101_434_1];
    np_assert_allclose!(actual, desired, rtol = 1e-14);
    // and these, which were NaN in one order of a and b, where xsf gives -2.29e82 and 1.03e-13
    // (Arb)
    let args = [(640.0, -128.04, -9.5, 0.4248), (-688.5, 147.62, -28.5, 0.3898)];
    let desired = [1.621_532_650_837_629_3e79, -1.825_420_437_271_943_6e-18];
    for ((a, b, c, x), y) in args.into_iter().zip(desired) {
        let actual = [xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x)];
        np_assert_allclose!(actual, [y, y], rtol = 1e-12);
        assert_eq!(actual[0].to_bits(), actual[1].to_bits());
    }
}

#[test]
fn test_hyp2f1_f64_large_opposite_signs() {
    // unlike in xsf: the 1 - x transformation if 0 < x < 1, and the 1/(1 - x) transformation if
    // x < 0 (with a and b of the same sign then), where the Pfaff and Euler forms cancel
    // (references: Arb)
    let actual =
        [xsf::hyp2f1(-975.2, 241.4, 0.00135, 0.39), xsf::hyp2f1(-880.55, -75.65, 0.091, -0.472)];
    let desired = [-1.469_500_580_590_125_3e-75, 8.556_925_324_233_141e83];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_c_lt_0() {
    // unlike in xsf: the backward recurrence in c from start values by other forms, if those of
    // its own form lose precision (references: Arb)
    let actual = [
        xsf::hyp2f1(758.1, 7.7, -151.2, -0.3),
        xsf::hyp2f1(16.1, 428.5, -230.7, -0.91),
        xsf::hyp2f1(-559.7, 42.4, -136.1, 0.27),
    ];
    let desired = [-8_553.508_491_012_49, -1_050_587.466_964_691_6, -2.046_272_877_671_090_4e24];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // from the start values at c + 1 and c + 2, which only help since the other forms nearly
    // succeeded (reference: Arb)
    let actual = xsf::hyp2f1(130.25, 125.25, 1.96, -145.46);
    np_assert_allclose!([actual], [-3.043_786_288_636_934_7e-280], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_c_lt_0_second_solution() {
    // unlike in xsf: for c < 0, from the second solution at 0 and one term of the 1/x
    // transformation (or of the 1 - x transformation if 0 < x < 1), since the other term loses
    // precision; xsf gives 8.29e102, 6.42e201, NaN and NaN (references: Arb)
    let args = [
        (-404.0, -124.75, -561.75, -8.25),
        (-390.0, -183.25, -671.25, -22.5),
        (-342.5, -186.25, -639.75, -36.25),
        (761.5, 211.5, -483.25, 0.04),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [
        -3.962_625_779_100_878e105,
        2.525_038_208_734_701_5e202,
        -1.322_092_032_643_801_3e232,
        3.140_497_898_387_16e-6,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_near_1() {
    // unlike in xsf: the 1 - x transformation of `hyp2f1_recur` for x >= 0.9, e.g. these were
    // -1.83e6, inf and NaN (references: Arb)
    let actual = [
        xsf::hyp2f1(40.1, -133.3, -38.6, 0.901),
        xsf::hyp2f1(-18.5, 75.0, 172.4, 0.901),
        xsf::hyp2f1(-2.9, -667.5, 1.8, 1.0 - 1e-14),
    ];
    let desired = [43.700_421_521_573_67, 1.614_794_417_233_753e-4, 9_480_012.959_767_725];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_integer_b_minus_a() {
    // unlike in xsf: the log case of the 1/x transformation of `hyp2f1_recur` if b - a is an
    // integer, e.g. the first was NaN (references: Arb)
    let actual = [
        xsf::hyp2f1(12.5, 812.5, 300.25, -25_000.0),
        xsf::hyp2f1(25.0, 25.0, -0.5, -1e6),
        xsf::hyp2f1(1.5, -50.5, 2.25, -100.0),
        xsf::hyp2f1(-40.5, -20.5, 2.25, -100.0),
        xsf::hyp2f1(25.5, 26.5, 27.25, -100.0),
    ];
    let desired = [
        3.490_725_745_476_049e-61,
        -8.873_628_661_650_908e-148,
        1.101_295_428_271_828_5e100,
        -2.154_549_921_333_831e66,
        7.125_787_400_134_933e-51,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
    // and of the 1 - x transformation near x = 1 if c - a - b is an integer (Arb)
    let actual = xsf::hyp2f1(-151.265_625, -167.078_125, -314.343_75, 0.996);
    np_assert_allclose!([actual], [1.867_105_850_532_836e-91], rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_polynomials() {
    // unlike in xsf: the recurrence in the degree, since the series cancel (references: Arb)
    let actual = [
        xsf::hyp2f1(-522.0, -123.5, 533.5, -0.9),
        xsf::hyp2f1(-522.0, -123.5, 533.5, -0.5),
        xsf::hyp2f1(-3.75, 10_675.75, -0.25, -5600.0),
    ];
    let desired =
        [6.485_127_939_853_472_5e-52, -4.200_865_766_102_724e-45, -2.143_761_498_506_179_6e29];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but their sums first, which are more accurate if their terms have one sign, e.g. this was
    // inf (Arb)
    let actual = xsf::hyp2f1(-193.25, -361.875, -784.25, -8.5);
    np_assert_allclose!([actual], [5.602_051_816_485_101e124], rtol = 1e-14);

    // and in double-double if those and the recurrences are inaccurate, e.g. near a root, where
    // these were -3.5662163414410e-9, -3.4478745961026 and -388.081847556425 (Arb)
    let args = [
        (-8.0, 18.016_500_331_508_873, 10.805_295_997_850_628, 0.908_756_475_07),
        (15.695_195_439_032_172, -283.0, -1.239_002_148_045_167_6, 0.022_551_819_176_619_81),
        (-4.999_999_999_790_708, -270.0, -0.418_468_928_086_676_85, -0.026_109_555_998_151_627),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [-3.566_216_341_442_087e-9, -3.447_874_596_102_943_6, -388.081_847_556_405_76];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    // but not where the low word of `1 / (1 - x)` is subnormal and so inexact, where that gave
    // 1.3e-4 and 1.3e-8 (exact)
    let actual = xsf::hyp2f1(-7.692_307_692_292_307e-308, -1.0, 10.0, -1.3e308);
    np_assert_allclose!([actual], [1.999_956_492_953_669_5e-12], rtol = 1e-6);
    let actual = xsf::hyp2f1(-6.666_666_653_333_333e-306, -1.0, 100.0, -1.5e307);
    np_assert_allclose!([actual], [2.000_000_054_253_873e-9], rtol = 1e-9);
}

#[test]
fn test_hyp2f1_f64_large_polynomials_in_1_over_x() {
    // unlike in xsf: the polynomial in 1/x, i.e. the term of the 1/x transformation that is left
    // for an integer a <= 0, e.g. these were 1.717382430e38, 6.50e50 and 7.6598e35 (references:
    // Arb)
    let args = [
        (-390.0, -60.3187, -547.357, -4.72131),
        (-315.0, -71.5, -426.5, -5.75),
        (-253.0, -50.75, -446.5, -7.5),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired =
        [1.717_382_032_883_699_5e38, 6.196_191_652_916_035_5e50, 7.661_996_786_019_798e35];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // and for an integer c - b <= 0, also with 0 < x < 1 (in 1 - 1/x), e.g. these were NaN and
    // 1.057e-11 (Arb)
    let actual = [
        xsf::hyp2f1(-177.375, -108.5, -208.5, -18.25),
        xsf::hyp2f1(-9.75, -958.25, -1022.25, 0.987),
    ];
    let desired = [6.812_936_622_801_884e131, 6.764_515_364_329_936e-12];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_large_polynomials_fallback() {
    // unlike in xsf: the result of `hyp2f1_recur` above its tolerance if the other methods are
    // worse, e.g. these were -3.97e161, 1.30e11, -4.37e70 and 1.3003887417e7 (references: Arb)
    let args = [
        (-578.0, 43.03, -211.37, -0.07),
        (-155.0, -13.11, -268.18, -5.78),
        (-405.0, 284.75, -356.78, -0.16),
        (-102.0, -11.3, -173.19, -5.63),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired = [
        5.727_126_898_619_048e-4,
        -6_321_179_103.061_271,
        3.088_344_341_495_295_5e-21,
        13_003_890.426_056_918,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_recur_fallback() {
    // unlike in xsf: also for non-polynomials, e.g. these were NaN, NaN and -5.62543958e-4, as the
    // recurrence in `a` loses most digits (references: Arb, Wolfram)
    let args = [
        (-43.53, 28.38, 0.850_66, 0.8976),
        (-58.8, 58.95, 29.1511, 0.883),
        (-39.4122, 29.003_22, -0.466_54, 0.506_55),
    ];
    let actual = args.map(|(a, b, c, x)| xsf::hyp2f1(a, b, c, x));
    let desired =
        [-7.920_164_213_871_762e-10, 1.473_319_142_367_457_4e-34, -5.625_453_919_713_718e-4];
    np_assert_allclose!(actual, desired, rtol = 1e-11);

    // and these were 9.47e-62 and 3.62e16 in xsf, and NaN before (Arb, Wolfram)
    let actual =
        [xsf::hyp2f1(702.2, 161.2, 0.1188, -0.4165), xsf::hyp2f1(198.4, 38.54, -105.4, -0.5366)];
    let desired = [1.176_016_677_157_160_2e-64, -7_416_153_466.305_625];
    np_assert_allclose!(actual, desired, rtol = 1e-9);
}

#[test]
fn test_hyp2f1_gauss_sum_at_1() {
    // The excess c - a - b is positive in every case. The fourth c - a rounds to -1.
    // References: Arb on the exact binary64 arguments.
    let tiny = f64::from_bits(1 << 44); // 2^-1030
    let args = [
        (0.1, 1.0, 1.1, 1.200_959_900_632_133_2e15),
        (0.1, 1.0, 1.100_000_000_000_000_3, 3.275_345_183_542_189_4e14),
        (-1.0 + f64::EPSILON / 2.0, 1.0 + f64::EPSILON, 1.0, -2.220_446_049_250_313_6e-16),
        (1.75, -1.25, 0.75 - f64::EPSILON / 2.0, 4.932_590_933_657_687e-16),
        (1.750_000_000_002, -1.25, 0.75 - f64::EPSILON / 2.0, 8.886_062_567_066_966e-12),
        (171.125, -170.999_999_999_999_97, 0.125_000_000_000_031_97, -6_099.234_743_153_951),
        (2.0_f64.powi(-600), -2.0_f64.powi(-600), 2.0_f64.powi(-1000), -2.0_f64.powi(800)),
        (tiny, tiny, 3.0 * tiny, 4.0 / 3.0),
    ];
    for (a, b, c, expected) in args {
        for (a, b) in [(a, b), (b, a)] {
            let real = xsf::hyp2f1(a, b, c, 1.0);
            np_assert_allclose!([real], [expected], rtol = 5e-15);
            for im in [0.0, -0.0] {
                let complex = xsf::hyp2f1(a, b, c, c64(1.0, im));
                assert_eq!(complex.re, real);
                assert_eq!(complex.im.to_bits(), 0);
            }
            assert_eq!(xsf::hyppfq([a, b], [c], 1.0), Some(real));
        }
    }

    // A true pole of a denominator gamma function gives zero.
    assert_eq!(xsf::hyp2f1(1.75, -1.25, 0.75, 1.0).to_bits(), (-0.0_f64).to_bits());
    for im in [0.0_f64, -0.0] {
        let y = xsf::hyp2f1(1.75, -1.25, 0.75, c64(1.0, im));
        assert_eq!(y.re.to_bits(), im.to_bits());
        assert_eq!(y.im.to_bits(), 0);
    }

    // A nonzero negative result can also round to -0.
    let y = xsf::hyp2f1(1_200.000_01, -401.25, 800.0, c64(1.0, 0.0));
    assert_eq!(y.re.to_bits(), (-0.0_f64).to_bits());
}

#[test]
fn test_hyp2f1_f64_divergent_at_1() {
    // unlike in xsf: the sign of the limit x -> 1- for c - a - b <= 0 (DLMF 15.4.21, 15.4.23; Arb
    // at x = 1 - 2^-3000), where xsf has +inf for the first two and NaN for the last two
    let inf = f64::INFINITY;
    let actual = [
        xsf::hyp2f1(-2.109, 0.891, -1.723, 1.0),
        xsf::hyp2f1(1.5, -0.5, 1.0, 1.0),
        xsf::hyp2f1(2.5, -1.5, 1.0, 1.0),
        xsf::hyp2f1(-2.5, 7.5, 3.5, 1.0),
        xsf::hyp2f1(202.75, 102.25, 303.0, 1.0),
        xsf::hyp2f1(-99.25, -44.5, -231.875, 1.0),
    ];
    assert_eq!(actual, [-inf, -inf, inf, -inf, inf, -inf]);
    // The excess is tiny and negative for the first case, and overflows for the second.
    for (a, b, c) in [(-1.0 + f64::EPSILON / 2.0, 2.0, 1.0), (1e308, 1e308, -0.5)] {
        for (a, b) in [(a, b), (b, a)] {
            assert_eq!(xsf::hyp2f1(a, b, c, 1.0), -inf);
            assert_eq!(xsf::hyp2f1(a, b, c, c64(1.0, 0.0)), c64(inf, 0.0));
        }
    }
    // +inf at a non-positive integer c, as for x < 1
    assert_eq!(xsf::hyp2f1(0.5, 0.5, -3.0, 1.0), inf);
    assert_eq!(xsf::hyp2f1(0.5, 0.5, -3.0, 0.5), inf);
}

#[test]
fn test_hyp2f1_f64_near_polynomial() {
    // unlike in xsf: not the polynomial if a, b, c - a or c - b is only near a negative integer,
    // e.g. c - a = -58 - 5e-15 for the first. These are 3.35e11, NaN, NaN, 4.6e49, -inf and
    // -3.95e-6 in xsf (Arb)
    let actual = [
        xsf::hyp2f1(65.9, -57.6, 7.9, 0.861),
        xsf::hyp2f1(8.8, 5.6, 2.8, -233.9),
        xsf::hyp2f1(-3.000_000_000_000_001, 9.8, 7.0, -1.45),
        xsf::hyp2f1(97.7, -11.000_000_000_000_002, -27.3, -0.649),
        xsf::hyp2f1(-6.2, 3.4, -8.2, -4.4),
        xsf::hyp2f1(-49.6, -83.999_999_999_999_94, -91.2, 0.773),
    ];
    let desired = [
        -5.460_690_518_245_093e-10,
        -7.659_905_486_275_666e-18,
        26.565_667_750_000_028,
        -11_568.209_528_386_858,
        0.005_964_933_841_297_133,
        1.981_179_670_464_318_1e-29,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hyp2f1_f64_tiny_c_integer_b() {
    // unlike in xsf: c - b = -n + c is only near a negative integer, where the polynomial misses
    // the tail of the series. These are inf in xsf (Arb).
    let actual = [
        xsf::hyp2f1(41.5, 4.0, -1.6e-154, -35_557.0),
        xsf::hyp2f1(-42.25, 9.0, -1.5e-14, 0.9998),
        xsf::hyp2f1(15.0, -34.25, -1e-15, 0.97),
    ];
    let desired =
        [6.500_894_648_275_649_5e-24, -2.692_035_181_549_083_2e-9, 7.064_137_056_127_026e-5];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_hyp2f1_f64_just_below_1() {
    // unlike in xsf: x within 1e-13 below 1 is not x = 1, where xsf gives inf for c - a - b <= 0
    // (the first three and the last) and the value at 1 else, here 32.27, 440.5 and -7.106 (Arb)
    let x = 0.999_999_999_999_99;
    let actual = [
        xsf::hyp2f1(0.5, 0.7, 1.0, x),
        xsf::hyp2f1(0.5, 0.7, 1.0, 0.999_999_999_999_999_9),
        xsf::hyp2f1(-2.109, 0.891, -1.723, x),
        xsf::hyp2f1(0.5, 0.49, 1.0, x),
        xsf::hyp2f1(1.5, 2.25, 3.76, x),
        xsf::hyp2f1(0.81, -0.04, 0.775, 0.999_999_999_999_997_6),
        xsf::hyp2f1(-3.5, 1.25, -2.5, x),
    ];
    let desired = [
        1_258.099_570_916_464,
        3_095.910_120_747_267_7,
        -11_864_866.701_448_103,
        9.532_032_686_970_835,
        115.683_428_221_038_92,
        -0.256_813_810_821_179_77,
        -44_263.439_899_139_13,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // but x = 1 is: Gauss's sum (DLMF 15.4.20) for c - a - b > 0, and the limit else
    let actual = [xsf::hyp2f1(0.5, 0.2, 1.0, 1.0)];
    np_assert_allclose!(actual, [1.449_724_260_959_791_1], rtol = 1e-15);
    assert_eq!(xsf::hyp2f1(0.5, 0.7, 1.0, 1.0), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_psi_expansion_underflow() {
    // unlike in xsf: a gamma factor of the larger sum of the psi function expansion (AMS55
    // 15.3.10-12) underflows to 0, where xsf returns the other sum alone, here 0 (Arb)
    let actual = xsf::hyp2f1(-57.593_75, 251.375, 134.781_25, 0.962);
    np_assert_allclose!([actual], [1.447_158_991_775_232_9e-19], rtol = 1e-8);
}

#[test]
fn test_hyp2f1_f64_log_case_cancels() {
    // unlike in xsf: an integer b - a for x < 0, or c - a - b for 0 < x < 1 (the last two), where
    // the limit of the 1/x or 1 - x transformation cancels: xsf gives 0, NaN, NaN, 1.732e114,
    // 3.148e12 and -6.611e-36 (Arb)
    let actual = [
        xsf::hyp2f1(-130.25, -106.25, 88.5, -56.5),
        xsf::hyp2f1(105.25, 137.25, -182.8, -19.875),
        // (b - a is 8.5e-14 from an integer)
        xsf::hyp2f1(239.1, 215.1 - 8.5e-14, -205.6, -74.1),
        xsf::hyp2f1(-90.75, -126.75, 74.95, -27.75),
        xsf::hyp2f1(-96.25, 186.5, 31.25, 0.93),
        // (a polynomial of degree 446, and c - a - b is 4.3e-14 from an integer)
        xsf::hyp2f1(-446.0, 51.25, -35.75 + 4e-14, 0.565),
    ];
    let desired = [
        -9.451_533_499_353_577e147,
        -8.900_607_975_593_249e-76,
        -2.530_210_061_769_191_5e-297,
        1.678_815_006_969_739_4e106,
        26.381_995_650_576_812,
        -7.385_094_653_794_124e-36,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-11);
}

#[test]
fn test_hyp2f1_f64_c_near_pole() {
    // unlike in xsf: c near (not at) a non-positive integer -n, where 2F1 is finite, and its power
    // series can jump by about 1 / (c + n) at the term n + 1; xsf gives inf (Arb)
    let actual = [
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, -0.4),
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, -0.99),
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, 0.4),
        xsf::hyp2f1(1.25, 2.0, -45.0 - 4e-14, 0.9),
        // (where the Pfaff transformation rounds c - b near -n)
        xsf::hyp2f1(7.25, 1e-13, -2.0 + f64::EPSILON, -0.7),
        xsf::hyp2f1(-2.75, 1e-13, -1.0 - f64::EPSILON, -0.3),
    ];
    let desired = [
        1.022_933_013_832_101_3,
        349.569_306_386_504_9,
        3_671_424_765.509_292_6,
        2.052_009_450_569_774_3e63,
        -362.316_350_911_936_1,
        112.785_458_910_322_3,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
    // where the Pfaff power factor underflows, or a is near (not at) an integer that would
    // terminate the series
    let c = -1.0 + 2f64.powi(-45);
    let actual = [
        xsf::hyp2f1(c + 61.0, c + 31.0, c, -300_000.0),
        xsf::hyp2f1(-1.0 + f64::EPSILON, 2.0, -2.0 + f64::EPSILON, -0.4),
    ];
    np_assert_allclose!(
        actual,
        [1.795_587_581_895_763e-295, 0.655_976_676_384_839_7],
        rtol = 1e-13
    );
    // where the terms are near overflow: 2F1(c, b; c; x) = (1 - x)^-b
    let c = -300.0 + 2f64.powi(-44);
    let actual = [xsf::hyp2f1(c, 300.0, c, 0.904)];
    np_assert_allclose!(actual, [2.082_716_166_487_401_1e305], rtol = 1e-13);
    // NaN if the series cancels more than its rounding estimate allows
    let y = xsf::hyp2f1(0.1234, 7.1, -8.000_000_000_000_002, 0.010_557_580_327_416_858);
    assert!(y.is_nan());
    // at c = -n it is a pole, unless the series terminates before the term n + 1
    assert_eq!(xsf::hyp2f1(1.25, 2.0, -45.0, -0.4), f64::INFINITY);
    let actual = [xsf::hyp2f1(-3.0, 2.0, -45.0 - 4e-14, 0.5)];
    np_assert_allclose!(actual, [1.068_974_630_021_141_6], rtol = 1e-15);
}

#[test]
fn test_hyp2f1_f64_polynomial_neg_int_c() {
    // unlike in xsf: the polynomial for a negative integer c below the integer a (or b), whose
    // terms can cancel completely; xsf gives NaN or no correct digit (Arb, Wolfram)
    let actual = [
        xsf::hyp2f1(-186.0, -26.75, -333.0, -9.0),
        xsf::hyp2f1(-2.5, -123.0, -618.0, -12.181_556_567_107_146),
        xsf::hyp2f1(-282.0, 497.137_097_191_493_7, -301.0, -0.029_204_871_958_424_083),
        // (by its Pfaff transformation)
        xsf::hyp2f1(-192.346_044_530_068_72, -303.0, -310.0, 0.290_678_693_500_093),
        // (b is an integer below c, so the polynomial is that of a)
        xsf::hyp2f1(-359.0, -401.0, -400.0, 0.012_918_049_583_765_393),
        // (by its reversal about 1, for x > 1)
        xsf::hyp2f1(-43.0, -22.5, -47.0, 1.023_180_005_600_207_6),
        xsf::hyp2f1(-45.0, -35.157_074_913_997_75, -53.0, 1.000_000_000_003_527),
    ];
    let desired = [
        -5.330_688_234_478_712e36,
        21.568_509_940_337_676,
        1.457_505_009_443_554_2e-6,
        9.959_364_424_242_371e-29,
        0.009_283_916_262_782_09,
        4.883_556_412_840_876e-19,
        -3.889_513_609_716_145e-24,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    // the sign of an overflow, also for x > 1 (about 1.8e623 and -3.1e379)
    assert_eq!(
        xsf::hyp2f1(-150.0, 39.883_200_484_784_936, -155.0, -8_782.652_226_040_95),
        f64::INFINITY
    );
    assert_eq!(xsf::hyp2f1(-227.0, -223.0, -400.0, 163.795_069_436_680_82), f64::NEG_INFINITY);
    // a subnormal b, both a and b polynomial (of degree 2, not 1e6), a degree of 1e6 whose terms
    // vanish after about 60, and x = 1 (DLMF 15.4.24)
    let actual = [
        xsf::hyp2f1(-1100.0, DENORM_MIN, -1101.0, 2.0),
        xsf::hyp2f1(-1e6, -2.0, -1_000_001.0, 0.25),
        xsf::hyp2f1(-1e6, 2.0, -1_000_001.0, 0.25),
        xsf::hyp2f1(-29.0, -27.0, -57.0, 1.0),
    ];
    let desired = [
        223.050_453_452_347_38,
        0.562_500_374_999_625,
        1.777_776_592_593_777_8,
        1.995_525_599_273_078_8e-15,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    assert_eq!(xsf::hyp2f1(-20.0, -20.0, -30.0, 1.0), 0.0);
    // its terms cancel by 1e80, so NaN (or the exact sum), not that of the terms up to where the
    // rest is negligible
    let y = xsf::hyp2f1(-2048.0, 512.0, -2049.0, -0.25);
    assert!(y.is_nan() || (y / 2.530_769_231_467_989e-50 - 1.0).abs() < 1e-13);
    // a pole if a is not an integer, also if c - a rounds to one (in the Euler transformation)
    assert_eq!(xsf::hyp2f1(-1.000_000_000_000_2, 0.25, -8192.0, 0.5), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_polynomial_neg_int_c_overflow() {
    // unlike in xsf: the sign of an overflowing polynomial for a negative integer c below the
    // integer a or b, e.g. of 1 - 4x + 38/9 x^2 and 1 + 4/3 x + x^2 for the first two, where the
    // power series stops at its first infinite term, or the recurrence in a gives NaN; xsf
    // gives -inf for the first two (references: exact)
    let inf = f64::INFINITY;
    for ((a, b, c, x), desired) in [
        ((-2.0, -20.0, -10.0, 1e308), inf),
        ((-2.0, 2.0, -3.0, -1e308), inf),
        ((5.0, -5.0, -7.0, -1.271_610_395_627_417_3e231), -inf),
        ((-4.0, -5.0, -7.0, 4.334_366_077_440_268e137), inf),
        // (by `hyp2f1_recur`, whose first term overflows at the scale of its sums)
        ((-11.0, 24.0, -20.0, 1.146_319_809_459_681_3e308), inf),
        ((-32.0, -136.199_875_679_093_6, -34.0, 8.858_494_428_306_513e307), inf),
        // (just above the cutoff, by 3.5e-16 and 2.5e-16, which the double-double sum of
        // `hyp2f1_recur` resolves)
        ((-3.0, 2.0_f64.powi(-960), -4.0, 2.759_926_699_133_471_5e199), inf),
        ((-2.0, 35.0, -40.0, -1.491_883_142_149_389e154), inf),
    ] {
        assert_eq!(xsf::hyp2f1(a, b, c, x), desired);
        assert_eq!(xsf::hyp2f1(b, a, c, x), desired);
    }
}

#[test]
fn test_hyp2f1_f64_polynomial_overflow() {
    // unlike in xsf: polynomials for c > 0 or a non-integer c whose power series is not finite,
    // e.g. of 1 - 4x/3 + x^2/2 and 1 - 8x + 8x^2 for the first two, where it stops at its first
    // infinite term, or the recurrence in a gives NaN; xsf gives -inf for the first two
    // (references: exact)
    let (inf, max) = (f64::INFINITY, f64::MAX);
    // (next to the pole c = -1, which the degree 2 passes)
    let c_1 = -1.0 + f64::EPSILON / 2.0;
    for ((a, b, c, x), desired) in [
        ((-2.0, 2.0, 3.0, 1e308), inf),
        ((-2.0, 2.0, 0.5, 1e308), inf),
        ((-22.0, 18.049_542_264_449_997, 24.0, 6.296_222_190_809_398e148), inf),
        ((-13.0, 2.362_904_706_447_58, 27.076_332_302_285_26, 1.248_008_371_891_642_5e96), -inf),
        ((-3.0, -3.0, -28.791_853_612_315_794, 3.675_816_243_990_002e122), -inf),
        ((-25.0, 36.527_592_225_483_26, -4.159_605_496_329_195_5, 4.708_584_498_582_971e71), inf),
        (
            (-6.0, -8.269_095_529_373_416, 2.511_868_208_292_315_6, -1.188_302_376_460_491_2e151),
            inf,
        ),
        (
            (-8.0, -4.555_226_504_260_197, -2.488_187_984_357_836, -2.120_778_378_045_006_4e159),
            inf,
        ),
        // (just below the cutoff, which the double-double sum resolves)
        ((-7.0, -6.219_584_855_883_66, 20.955_424_871_235_71, 1.171_551_334_021_426_5e45), max),
        (
            (-2.0, -0.976_556_258_795_454_8, -5.143_678_704_031_934_5, -4.090_963_831_756_783e155),
            -max,
        ),
        // (where a term ratio over- or underflows, but the terms don't, and one with a finite
        // value, which was -inf)
        ((-2.0, -2.0, 1e-310, 1.0), inf),
        ((-2.0, 1.0, 1e-310, 2.0), inf),
        ((-2.0, -2.0, 1e-310, -3.0), inf),
        ((-2.0, 2e200, 1e308, 1e300), inf),
        ((-2.0, 1e308, 1e307, 2.0), 361.0),
        ((-2.0, 1.0, c_1, 1e200), -inf),
        ((-2.0, 1.0, c_1, -1e200), -inf),
        // (Euler transformations whose c - a rounds to a negative integer, where the polynomial
        // overflows, but not its product with (1 - x)^(c - a - b); xsf gives NaN and inf;
        // references: Arb)
        (
            (20.0, 19.578_718_423_577_527, 3.761_908_130_024e-312, -1_240.829_428_855_755_4),
            1.302_688_042_625_831_2e251,
        ),
        (
            (44.0, 70.443_210_106_362_3, -3.5e-323, -1.226_182_118_621_034_2),
            -6.120_084_936_549_461_6e302,
        ),
        // (finite by the other paths, e.g. the Pfaff series where (1 - x)^(c - a - b)
        // underflows, which this leaves correctly rounded)
        ((60.0, 1e-12, -1e-315, -300_000.0), 1.000_000_001_501_041_5e303),
        ((2.0, 1.25, 2e-309, -8.000_976_562_5), 9.665_443_073_997_591e302),
        (
            (
                -43.0,
                -103.894_138_368_212_03,
                4.176_330_477_074_113e-309,
                -1.639_378_690_091_260_6e-18,
            ),
            -1.753_654_077_006_681_6e294,
        ),
    ] {
        assert_eq!(xsf::hyp2f1(a, b, c, x), desired);
        assert_eq!(xsf::hyp2f1(b, a, c, x), desired);
    }
}

#[test]
fn test_hyp2f1_f64_degree_1() {
    // unlike in xsf: `1 - a x / c` also for b = -1, where xsf gives NaN (references: exact)
    let actual = [
        xsf::hyp2f1(10.0, -1.0, -10.0, -1.0),
        xsf::hyp2f1(1.0, -1.0, -10.0, -10.0),
        xsf::hyp2f1(-10.0, -1.0, -10.0, 1.0),
    ];
    assert_eq!(actual, [0.0; 3]);
    // unlike in xsf: not rounded to 0 below 1e-15; xsf gives 0 (in this order of a and b), 0 and
    // 2.2e-16 (references: exact)
    let x = 1.0 - f64::EPSILON / 2.0;
    let actual = [xsf::hyp2f1(-1.0, 1.0, 1.0, x), xsf::hyp2f1(1.0, -1.0, 1.0, x)];
    assert_eq!(actual, [f64::EPSILON / 2.0; 2]);
    // (where `b x` is subnormal)
    let t = DENORM_MIN;
    assert_eq!([xsf::hyp2f1(-1.0, t, t, 0.5), xsf::hyp2f1(t, -1.0, t, 0.5)], [0.5; 2]);
    let actual = [
        xsf::hyp2f1(-1.0, -18.0, -6.361_004_257_813_726, 0.353_389_125_434_095_84),
        xsf::hyp2f1(-8.0, -1.0, 7.0, -0.874_999_999_999_999_8),
    ];
    np_assert_allclose!(
        actual,
        [1.745_358_090_684_158e-16, 2.537_652_627_714_643_4e-16],
        rtol = 1e-15
    );
}

#[test]
fn test_hyp2f1_f64_c_eq_a_or_b_neg_int() {
    // unlike in xsf: for c = a = -m or c = b = -m, the polynomial of degree m of DLMF 15.2.5, also
    // for x > 1 and m = 0; xsf gives inf, inf, 0.5, 2^-10, inf, (1 - x)^20 and inf (references:
    // exact)
    let actual = [
        xsf::hyp2f1(10.0, -10.0, -10.0, 1.0),
        xsf::hyp2f1(-10.0, -10.0, -10.0, -10.0),
        xsf::hyp2f1(-10.0, 1.0, -10.0, -1.0),
        xsf::hyp2f1(-10.0, 10.0, -10.0, -1.0),
        xsf::hyp2f1(2.5, -3.0, -3.0, 3.0),
        xsf::hyp2f1(-10.0, -20.0, -10.0, 1.0 / 3.0),
        xsf::hyp2f1(0.0, 2.5, 0.0, 3.0),
    ];
    let desired = [184_756.0, 25_937_424_601.0, 1.0, 60_172.0, 225.0625, 44_458.0 / 59_049.0, 1.0];
    np_assert_allclose!(actual, desired, rtol = 1e-15);
    // (whose terms cancel by 1e18 and 3e20; xsf gives inf and NaN)
    let actual = [
        xsf::hyp2f1(-50.0, -17.704_152_726_819_963, -50.0, 1.124_144_544_693_068),
        xsf::hyp2f1(-51.0, -50.0, -50.0, 0.656_800_541_976_558_8),
    ];
    np_assert_allclose!(
        actual,
        [-4.974_610_285_018_047e-13, 4.887_745_972_669_103e-10],
        rtol = 1e-14
    );
    // (exactly 0; xsf gives 0.5, 2^-10 and NaN)
    let actual = [
        xsf::hyp2f1(-9.0, 1.0, -9.0, -1.0),
        xsf::hyp2f1(-9.0, -10.0, -9.0, 0.5),
        xsf::hyp2f1(-10.0, -3.0, -10.0, 1.0),
    ];
    assert_eq!(actual, [0.0; 3]);
    // (1 - x)^400 for a = b = c = -400, which overflows
    assert_eq!(xsf::hyp2f1(-400.0, -400.0, -400.0, 6.984_779_578_287_557), f64::INFINITY);
    // (the sum of (-1)^k for k <= 10001)
    assert_eq!(xsf::hyp2f1(-10_001.0, 1.0, -10_001.0, -1.0), 0.0);
    // ((1 - x)^m, where 1 - x rounds to 1; reference: Arb)
    let actual = [xsf::hyp2f1(-1e20, -1e20, -1e20, 1e-20)];
    np_assert_allclose!(actual, [0.367_879_441_171_442_33], rtol = 1e-15);
    // (3^m and 2^-m, which overflow and underflow)
    let actual = [-2.0, 0.5].map(|x| xsf::hyp2f1(-1e20, -1e20, -1e20, x));
    assert_eq!(actual, [f64::INFINITY, 0.0]);
    // (a subnormal (1 - x)^m; reference: mpmath)
    let actual = [xsf::hyp2f1(-1_000_001.0, -1_000_001.0, -1_000_001.0, 0.000_708_5)];
    np_assert_allclose!(actual, [1.559_546_624_838_864e-308], rtol = 1e-15);
    // (the smallest subnormal, which exp(ln((1 - x)^m)) rounds up to)
    let actual = xsf::hyp2f1(-1e7, -1e7, -1e7, 7.451_054_586_157_437e-5);
    assert_eq!(actual, DENORM_MIN);
    // (a tiny x and a huge m, where `ln(1 - x)` must be accurate relative to x; reference: mpmath)
    let a = -12_772_208_543_222_726_656.0;
    let actual = [xsf::hyp2f1(a, a, a, 2.0_f64.powi(-54))];
    np_assert_allclose!(actual, [1.216_780_750_623_399e-308], rtol = 1e-15);
    // (NaN for c = -inf)
    assert!(xsf::hyp2f1(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, 0.5).is_nan());
}

#[test]
fn test_hyp2f1_f64_c_pole() {
    // unlike in xsf: +inf at a pole of c (where no a or b ends the series before it) for any
    // finite x != 0; xsf gives -inf for the first three, and NaN for the others
    let actual = [
        xsf::hyp2f1(10.0, 1.0, -10.0, 10.0),
        xsf::hyp2f1(1.0, 1.0, -1.0, 10.0),
        xsf::hyp2f1(5.0, 1.0, -9.0, 1.5),
        xsf::hyp2f1(750.274_991_146_857_7, 22.0, -670.0, -7_048.169_400_335_604),
        xsf::hyp2f1(-22.992_639_872_801_874, 54.0, -752.0, -6_207.512_565_903_636),
        xsf::hyp2f1(
            -4.673_419_326_051_837,
            289.810_719_716_530_9,
            -995.0,
            -0.855_413_125_222_624_2,
        ),
    ];
    assert_eq!(actual, [f64::INFINITY; 6]);
    // (also for a non-integer a next to the integer above c, and for c = 0, where xsf gives 1.25
    // and -inf; but NaN for c = -inf)
    let actual =
        [xsf::hyp2f1(-1.0 - f64::EPSILON, 1.0, -2.0, 0.5), xsf::hyp2f1(10.0, 1.0, 0.0, 10.0)];
    assert_eq!(actual, [f64::INFINITY; 2]);
    assert!(xsf::hyp2f1(1.0, 1.0, f64::NEG_INFINITY, 0.5).is_nan());
}

#[test]
fn test_hyp2f1_f64_neg_int_c_x_1() {
    // unlike in xsf: the Chu-Vandermonde product at x = 1 for a negative integer c below the
    // integer a or b, also for small parameters; xsf gives NaN (references: exact)
    let actual = [xsf::hyp2f1(-2.0, -3.0, -4.0, 1.0), xsf::hyp2f1(-3.0, -2.0, -4.0, 1.0)];
    assert_eq!(actual, [0.0; 2]);
    // (a nonzero product, where xsf is an ulp below)
    let actual = [xsf::hyp2f1(-4.0, 2.0, -6.0, 1.0), xsf::hyp2f1(2.0, -4.0, -6.0, 1.0)];
    assert_eq!(actual, [14.0 / 3.0; 2]);
}

#[test]
fn test_hyp2f1_f64_polynomial_x_1() {
    // unlike in xsf: the Chu-Vandermonde product (DLMF 15.4.24) at x = 1 for a non-integer c,
    // where xsf loses most digits, and gives 1.3365133773992304e17 for the fifth (references:
    // Arb)
    let args = [
        (-35.0, -30.83, -37.75, 1.015_046_708_249_615_6e-13),
        (-35.0, -30.5, -37.75, -1.889_066_610_240_864e-13),
        (-35.0, -28.3, -37.75, -2.423_744_082_535_354_7e-14),
        (-34.0, -30.83, -36.25, 1.605_183_033_836_263_2e-11),
        (-151.0, 102.4, 163.4, 3.857_040_397_885_821e-40),
        (-64.0, 44.8, 145.4, 5.748_680_614_749_967e-9),
        (-141.0, -48.1, -194.9, 1.296_718_313_473_943e-39),
    ];
    for (a, b, c, expected) in args {
        let actual = [xsf::hyp2f1(a, b, c, 1.0), xsf::hyp2f1(b, a, c, 1.0)];
        np_assert_allclose!(actual, [expected; 2], rtol = 1e-15);
    }
    // (also after a common parameter cancels)
    let actual = xsf::hyppfq([-62.0, -35.0, -30.83], [-37.75, -62.0], 1.0);
    np_assert_allclose!([actual.unwrap()], [args[0].3], rtol = 1e-15);
    // (correctly rounded to a subnormal; reference: exact rational)
    let actual = xsf::hyp2f1(-20.0, 29_389_293_564_442_344.0, 29_389_293_564_442_348.0, 1.0);
    assert_eq!(actual, 1.864_490_772_489_121e-308);
}

#[test]
fn test_hyp2f1_f64_x_gt_1() {
    // unlike in xsf: NaN where 2F1 is complex, in either order of a and b; xsf gives inf
    for (a, b, c, x) in [(32.000_000_000_000_01, 75.0, 97.0, 2.0), (32.0, 75.0, 97.0, 2.0)] {
        assert!(xsf::hyp2f1(a, b, c, x).is_nan());
        assert!(xsf::hyp2f1(b, a, c, x).is_nan());
    }
    // (but real for a polynomial, and for (1 - x)^-2 (1 - x / 2) by DLMF 15.8.1; exact)
    assert_eq!(xsf::hyp2f1(-2.0, 1.0, 3.0, 3.0), 0.5);
    assert_eq!(xsf::hyp2f1(3.0, 1.0, 2.0, 3.0), -0.125);
    // and for c = b, (1 - x)^-a (DLMF 15.4.6), here 2^1040, even where c - a rounds to c
    let c = 2.0_f64.powi(60);
    assert_eq!(xsf::hyp2f1(20.0, c, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(c, 20.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(-3.0, 2.5, 2.5, 3.0), -8.0);
    // (and inf otherwise, as for the other complex powers, also where c - a - b rounds to an
    // integer, and for x = inf)
    let c = 2.0_f64.powi(52);
    for (a, b, c, x) in [(0.5, 2.5, 2.5, 3.0), (1.25, c, c, 3.0), (0.5, 2.5, 2.5, f64::INFINITY)] {
        assert_eq!(xsf::hyp2f1(a, b, c, x), f64::INFINITY);
        assert_eq!(xsf::hyp2f1(b, a, c, x), f64::INFINITY);
    }
    // inf for (1 - x)^(c - a - b) times a polynomial (DLMF 15.8.1) where c - a - b rounds, here
    // to -25.000000000000114, as it overflows (about 2^1294.6)
    let c = -1024.0 + 2.0_f64.powi(-43);
    assert_eq!(xsf::hyp2f1(c + 1.0, 24.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(24.0, c + 1.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    // (but an exponent within `HYP2F1_EPS` of an integer counts as one: (-1/4)^2 and 1)
    let actual = [xsf::hyp2f1(2.5, -1.999_999_999_999_999_8, 2.5, 1.25)];
    np_assert_allclose!(actual, [0.0625], rtol = 1e-14);
    assert_eq!(xsf::hyp2f1(-2.0_f64.powi(-55), 2.5, 2.5, 1.25), 1.0);
    // inf as in xsf where c - a rounds to c, as the skipped Euler transformation overflows (its
    // polynomial of degree 128 is about 2^1223)
    let c = 2.0_f64.powi(59);
    assert_eq!(xsf::hyp2f1(24.0, c + 128.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
    assert_eq!(xsf::hyp2f1(c + 128.0, 24.0, c, 1.0 + f64::EPSILON), f64::INFINITY);
}

#[test]
fn test_hyp2f1_f64_rounded_euler_exponent_on_cut() {
    let a = 32.000_000_000_000_01;
    for (a, b) in [(a, 106.0), (106.0, a)] {
        assert!(xsf::hyp2f1(a, b, 97.0, 2.0).is_nan());
    }
    for (a, b) in [(32.0, 106.0), (106.0, 32.0)] {
        let y = xsf::hyp2f1(a, b, 97.0, 2.0);
        np_assert_allclose!([y], [3.705_699_243_281_37e-4], rtol = 1e-13);
    }

    // The Euler polynomial is exactly 1 - x/2, so a nonintegral power still gives zero at x=2.
    let a = f64::EPSILON / 2.0;
    let b = 1.0 + f64::EPSILON;
    let c = f64::EPSILON;
    for (p, q) in [(a, b), (b, a)] {
        assert_eq!(xsf::hyp2f1(p, q, c, 2.0).to_bits(), (-0.0_f64).to_bits());
        for x in [2.0 - f64::EPSILON, 2.0 + 2.0 * f64::EPSILON] {
            assert!(xsf::hyp2f1(p, q, c, x).is_nan());
        }
    }
}

#[test]
fn test_hyp2f1_c64_centered_population_polynomial() {
    // Pearson, Olver and Porter, Table A.2 case 16, with an off-axis ordinate.
    let z = c64(core::f64::consts::FRAC_1_SQRT_2, 0.01);
    for (a, b) in [(-100.0, -200.0), (-200.0, -100.0)] {
        let y = xsf::hyp2f1(a, b, -299.999_999_999, z);
        np_assert_allclose!(
            [y.re, y.im],
            [-1.832_446_212_804_294_4e-32, -2.699_124_054_158_986e-31],
            rtol = 1e-13
        );
    }

    // Arb at both binary64 inputs around a near-integer population root.
    for (x, expected) in [
        (0.955_144_145_039_975_6, 4.465_417_203_032_561e-124),
        (0.955_144_145_039_975_7, -1.856_100_292_981_34e-123),
    ] {
        let y = xsf::hyp2f1(-100.0, -200.0, -224.999_999_999_999_97, c64(x, -0.0));
        assert!((y.re - expected).abs() < 1e-131);
        assert_eq!(y.im.to_bits(), 0.0_f64.to_bits());
    }
}

#[test]
fn test_hyp2f1_c64_unresolved_population_polynomial() {
    let z = c64(1.0, 1.0);
    let y = xsf::hyp2f1(-1200.0, -2000.0, -3200.0, z);
    // Arb finite sum. A future solution can return a finite result at this scale.
    let expected = c64(-1.982_714_984_309_758e-98, 2.089_452_355_842_588e-98);
    assert!(y.is_nan() || (y - expected).norm() <= expected.norm());
}

#[test]
fn test_hyp2f1_c64_nonnegative_population_polynomial() {
    for imag in [0.0, -0.0] {
        // Pearson, Olver and Porter, Table A.2 case 16.
        let z = c64(core::f64::consts::FRAC_1_SQRT_2, imag);
        let y = xsf::hyp2f1(-100.0, -200.0, -299.999_999_999, z);
        np_assert_allclose!([y.re], [2.653_635_302_903_685_6e-31], rtol = 1e-13);
        assert_eq!(y.im.to_bits(), 0.0_f64.to_bits());

        // N = K + n - 1 gives a zero constant coefficient in 1 - z.
        let y = xsf::hyp2f1(-100.0, -200.0, -299.0, c64(1.0, imag));
        assert_eq!(y.re, 0.0);
        assert_eq!(y.im.to_bits(), 0.0_f64.to_bits());
    }
}

#[test]
fn test_hyp2f1_c64_real_polynomial_near_one() {
    // Arb; the terms of the series cancel near z = 1 (jorenham/xsf-rust#2211).
    for (a, b, c, x, expected) in [
        (-44.0, 1.06, -1.856_984_912_833_326, 1.0, -0.035_524_249_588_294_71),
        (-30.0, 20.25, -9.75, 1.0, 2_276_618_356.587_293_6),
        (
            -70.0,
            34.547_259_441_692_82,
            -23.998_630_037_442_19,
            1.000_000_000_011_891_6,
            -716_883_681.526_757_1,
        ),
        (
            0.914_864_001_003_731_6,
            -75.0,
            -24.971_086_673_375_95,
            0.999_999_999_013_385,
            -2.156_742_200_728_316_6,
        ),
    ] {
        for imag in [0.0, -0.0] {
            let y = xsf::hyp2f1(a, b, c, c64(x, imag));
            np_assert_allclose!([y.re], [expected], rtol = 4e-16);
            assert_eq!(y.im.to_bits(), 0.0_f64.to_bits());
        }
    }

    // Arb: 0.5000000000000000001 and -0.5000000000000000004 times the smallest subnormal, just
    // past the ties that the rounded products sit on.
    let (p80, p81) = (2.0_f64.powi(80), 2.0_f64.powi(81));
    for imag in [0.0, -0.0] {
        let y = xsf::hyp2f1(-1075.0, p80, p81, c64(1.0, imag));
        assert_eq!(y.re.to_bits(), DENORM_MIN.to_bits());
        let y = xsf::hyp2f1(-1075.0, -3.0 * p80, -p81, c64(1.0, imag));
        assert_eq!(y.re.to_bits(), (-DENORM_MIN).to_bits());
    }

    // Exact rationals (Chu-Vandermonde); the rounding of 2^-1022 and below is open to the
    // double-double bound, and the series is far off.
    for (a, bits) in [
        (-1022.0, 0x0010_0000_0000_0000),
        (-1023.0, 0x0008_0000_0000_0000),
        (-1024.0, 0x0004_0000_0000_0000),
    ] {
        for imag in [0.0, -0.0] {
            let y = xsf::hyp2f1(a, -p80, -p81, c64(1.0, imag));
            assert_eq!(y.re.to_bits(), bits);
        }
    }

    // Arb; the real overload keeps its rounding near the subnormal range.
    for (a, b, bits) in [
        (-1022.0, -549_755_813_888.5, 0x000f_ffff_805f_f1fb),
        (-1023.0, -549_755_813_888.5, 0x0007_ffff_c010_00ff),
        (-1024.0, -549_755_813_888.5, 0x0003_ffff_dff8_0080),
        (-1025.0, -549_755_813_888.75, 0x0001_ffff_efef_fd40),
    ] {
        let y = xsf::hyp2f1(a, b, -1_099_511_627_776.0, 1.0);
        assert_eq!(y.to_bits(), bits);
    }

    // Arb: 70368752623712.5049 times the smallest subnormal, which rounds up.
    for imag in [0.0, -0.0] {
        let y = xsf::hyp2f1(-1028.0, 2.0_f64.powi(41), 2.0_f64.powi(42), c64(1.0, imag));
        assert_eq!(y.re.to_bits(), 0x0000_4000_0080_e061);
    }
}

#[test]
fn test_hyp2f1_c64_beta_binomial_polynomial() {
    // Pearson, Olver and Porter, Table A.2 case 18, and its off-axis continuation.
    for (a, b) in [(-300.0, 5.0), (5.0, -300.0)] {
        for y in [0.0, -0.0, DENORM_MIN] {
            let value = xsf::hyp2f1(a, b, 10.0, c64(0.5, y));
            np_assert_allclose!([value.re], [1.661_006_238_211_309e-7], rtol = 1e-14);
            assert_eq!(value.im.to_bits(), if y == 0.0 { 0 } else { (-0.0_f64).to_bits() });
        }
        let value = xsf::hyp2f1(a, b, 10.0, c64(0.5, 0.01));
        np_assert_allclose!(
            [value.re, value.im],
            [1.651_579_371_189_262_2e-7, -1.612_740_700_807_058e-8],
            rtol = 1e-14
        );
    }

    let y = DENORM_MIN;
    let value = xsf::hyp2f1(-3.0, 0.125, 0.25, c64(2.0, y));
    assert_eq!(value, c64(0.0, -y));

    // The interior probabilities are subnormal, while both endpoint probabilities are normal.
    let value = xsf::hyp2f1(-32.0, y, 16.0 * y, c64(2.0, y));
    assert_eq!(value.re, 1.0);
}

#[test]
fn test_hyp2f1_c64_uniform_beta_binomial_polynomial() {
    // The 1000 equally weighted powers of -i sum to zero.
    for (a, b) in [(-999.0, 1.0), (1.0, -999.0)] {
        assert_eq!(xsf::hyp2f1(a, b, 2.0, c64(1.0, 1.0)), c64(0.0, 0.0));
    }
    let value = xsf::hyp2f1(-1000.0, 1.0, 2.0, c64(1.0, 1.0));
    assert_eq!(value, c64(1.0 / 1001.0, 0.0));

    // At x=-1 the derivative is -sum(k*2^(k-1), k=1..10)/11.
    let y = DENORM_MIN;
    let value = xsf::hyp2f1(-10.0, 1.0, 2.0, c64(-1.0, y));
    assert_eq!(value.im, -f64::from_bits(838));
}

#[test]
fn test_hyp2f1_c64_nonterminating_beta_endpoints() {
    // Neighbors of the degree-300 beta-binomial polynomial, checked with Arb.
    let degree = (-300.0_f64).to_bits();
    for ([a, b, c], [x, y], [re, im]) in [
        (
            [5.0, f64::from_bits(degree + 1), 10.0],
            [0.5, 0.875],
            [-4.379_935_823_703_022e-8, -3.095_636_837_235_489e-8],
        ),
        (
            [5.0, f64::from_bits(degree - 1), 10.0],
            [0.5, 0.875],
            [-4.379_935_823_702_598e-8, -3.095_636_837_236_046e-8],
        ),
        (
            [4.96, -300.01, 10.26],
            [0.5, 0.875],
            [-8.234_050_476_240_338e-9, -9.331_698_463_661_137e-9],
        ),
        (
            [5.0, -300.01, 10.0],
            [0.3, 0.35],
            [-1.266_497_033_807_951_7e-7, 2.263_173_931_583_558_2e-7],
        ),
        (
            [5.0, -300.01, 10.0],
            [0.1, 0.1],
            [-6.794_542_122_765_471e-5, 3.268_748_884_411_728_4e-5],
        ),
    ] {
        for sign in [1.0, -1.0] {
            let z = c64(x, sign * y);
            let value = xsf::hyp2f1(a, b, c, z);
            np_assert_allclose!([value.re, value.im], [re, sign * im], rtol = 1e-13);
            assert_eq!(value, xsf::hyp2f1(b, a, c, z));
        }
    }
}

#[test]
fn test_hyp2f1_c64_continued_beta_endpoints() {
    // Upper parameters outside (0, c), Arb at the exact binary64 inputs, also at conj(z).
    for ([a, b, c], z, expected) in [
        (
            [-1.111, -186.000_000_01, 6.684],
            c64(0.3024, -0.9576),
            c64(8.200_280_364_814_066, -39.759_998_653_624_79),
        ),
        (
            [9.376, -220.01, 4.227],
            c64(0.641, 0.8724),
            c64(345.211_447_704_520_9, -468.688_508_593_692_35),
        ),
        (
            [-8.919, -290.1, -14.603],
            c64(0.5083, -0.8412),
            c64(4_977_639_494_430.354, 394_482_200_962.892_1),
        ),
    ] {
        for (z, expected) in [(z, expected), (z.conj(), expected.conj())] {
            for actual in [xsf::hyp2f1(a, b, c, z), xsf::hyp2f1(b, a, c, z)] {
                np_assert_allclose!(
                    [actual.re, actual.im],
                    [expected.re, expected.im],
                    rtol = 1e-13
                );
            }
        }
    }
}

#[test]
fn test_hyp2f1_c64_beta_endpoints_tiny_imaginary_part() {
    // The series retains a small imaginary part that is poorly resolved by the contour.
    let value = xsf::hyp2f1(1.0, -48.125, 2.0, c64(0.01, 1e-24));
    np_assert_allclose!(
        [value.re, value.im],
        [0.793_180_702_499_857_2, -1.766_655_712_908_216_4e-23],
        rtol = 1e-13
    );
}

#[test]
fn test_hyp2f1_c64_real_beta_binomial_projection() {
    // This endpoint is exactly halfway between two normal floats.
    let midpoint = xsf::hyp2f1(-30.0, 1.5, 3.0, c64(1.0, 0.0));
    assert_eq!(midpoint.re, 0.012_615_460_793_392_621);

    // At z = 1 the polynomial equals the ratio (c-b)_n / (c)_n.
    let endpoint = xsf::hyp2f1(-1000.0, 49.99, 50.0, c64(1.0, 0.0));
    np_assert_allclose!([endpoint.re], [1.963_950_023_010_045e-90], rtol = 1e-14);
    assert_eq!(endpoint.im.to_bits(), 0);

    // Adjacent lower parameters straddle 1.5 times the smallest subnormal.
    let c = 1_188.019_841_343_689_2_f64;
    let next_c = f64::from_bits(c.to_bits() + 1);
    for (c, units) in [(c, 2), (next_c, 1)] {
        let value = xsf::hyp2f1(-300.0, c - 1.0, c, c64(1.0, 0.0));
        assert_eq!(value.re.to_bits(), units);
        assert_eq!(value.im.to_bits(), 0);
    }
}

#[test]
fn test_hyp2f1_c64_ferrers_half_integer_order() {
    // Hypergeometric factors of the half-integer order Ferrers modes at degree -1/2.
    // Arb references at the exact binary64 inputs.
    for (n, x, y, re, im) in [
        (64, 0.5, 0.0, 1.996_097_839_737_260_4, 0.0),
        (63, 0.5, 1e-9, 3.219_506_205_712_21e-14, -2.535_044_256_466_293_4e-7),
        (63, 0.5, 2.0_f64.powi(-515), 2.798_296_521_731_244e-306, -2.363_402_979_174_737_7e-153),
        (256, 0.1, 0.001, 0.999_902_195_854_822_2, -9.776_095_906_174_828e-7),
        (32, 1.25, 0.0, 0.0, 5.216_678_312_147_949e22),
    ] {
        for sign in [-1.0, 1.0] {
            let actual = xsf::hyp2f1(0.5, 0.5, 0.5 - f64::from(n), c64(x, sign * y));
            np_assert_allclose!([actual.re, actual.im], [re, sign * im], rtol = 2e-13, atol = 0.0);
        }
    }
    assert_eq!(xsf::hyp2f1(0.5, 0.5, -62.5, c64(0.5, 0.0)), c64(0.0, 0.0));
    // The large-order differences are individually smaller than an ulp of the sum.
    assert_eq!(
        xsf::hyp2f1(0.5, 0.5, -299_999.5, c64(-2.0_f64.powi(-20), 0.0)).re,
        1.000_000_000_000_794_7
    );
    // The linear term rounds to exactly 64 subnormal units.
    for sign in [-1.0, 1.0] {
        let z = c64(0.0, sign * f64::from_bits(1 << 20));
        assert_eq!(xsf::hyp2f1(0.5, 0.5, -4095.5, z).im, -sign * f64::from_bits(64));
    }
}

#[test]
fn test_hyp2f1_c64_negative_lower_parameter() {
    // Pearson et al., Table A.2, case 19, and nearby complex arguments (Arb).
    let expected_re = -3.852_027_081_523_918_5e32;
    for y in [0.0, -0.0] {
        let actual = xsf::hyp2f1(10.0, 5.0, -300.5, c64(0.5, y));
        np_assert_allclose!([actual.re], [expected_re], rtol = 4e-14);
        assert_eq!(actual.im.to_bits(), 0);
    }
    let y = DENORM_MIN;
    for sign in [1.0, -1.0] {
        let actual = xsf::hyp2f1(10.0, 5.0, -300.5, c64(0.5, sign * y));
        np_assert_allclose!(
            [actual.re, actual.im],
            [expected_re, -sign * 2.348_714_413_894_915_6e-288],
            rtol = 4e-14,
            atol = 0.0
        );
    }
    for (z, expected) in [
        (c64(0.5, 0.125), c64(-2.349_497_872_024_921e32, -9.112_118_651_455_808e31)),
        (
            c64(0.8, 0.380_580_885_090_989_36),
            c64(-1.208_133_589_671_797_6e128, 9.531_361_570_497_277e127),
        ),
        (
            c64(0.8, 0.380_580_885_090_989_4),
            c64(-1.208_133_589_671_798_1e128, 9.531_361_570_496_575e127),
        ),
    ] {
        let actual = xsf::hyp2f1(10.0, 5.0, -300.5, z);
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = 4e-14);
    }
}

#[test]
fn test_hyp2f1_c64_negative_lower_parameter_legacy_accurate() {
    // Arb at the exact binary64 inputs. The c-down value is less accurate.
    let actual = xsf::hyp2f1(
        19.714_877_719_930_925,
        1.135_255_778_495_805_7,
        -20.875,
        c64(0.01, 0.250_827_529_860_137_1),
    );
    np_assert_allclose!(
        [actual.re, actual.im],
        [-0.708_328_569_167_549_8, -0.060_416_415_146_359_134],
        rtol = 5e-13,
        atol = 0.0
    );
}

#[test]
fn test_hyp2f1_c64_negative_lower_parameter_c_down() {
    // Arb at the exact binary64 inputs. The connection cancels completely, or its gamma ratios
    // overflow.
    let actual = xsf::hyp2f1(0.125, 250.0, -2.75, c64(0.8125, DENORM_MIN));
    np_assert_allclose!(
        [actual.re, actual.im],
        [-4.010_669_060_041_866e189, -2.679_470_116_244_821e-131],
        rtol = 5e-14,
        atol = 0.0
    );
    let actual = xsf::hyp2f1(10.0, 400.0, -600.5, c64(0.351_562_5, 0.0));
    np_assert_allclose!([actual.re], [-1.718_794_060_874_630_4e231], rtol = 1e-14);
}

#[test]
fn test_hyp2f1_c64_negative_lower_parameter_connection() {
    // Arb at the exact binary64 inputs. The power z^(1-c) amplifies rounding for large -c.
    for (a, b, c, z, expected, rtol) in [
        (
            48.182_200_865_203_15,
            36.346_344_282_079_976,
            -462.5,
            c64(0.367_892_340_752_568_55, 0.0),
            c64(-5.961_484_280_097_055_4e32, 0.0),
            1e-13,
        ),
        (
            0.218_038_636_317_081_54,
            17.265_619_719_555_445,
            -499.003_703_751_683_4,
            c64(0.212_418_840_135_297_27, 0.737_066_445_331_071_8),
            c64(0.998_315_176_639_864_1, -0.005_506_388_833_780_157),
            1e-15,
        ),
        (
            37.017_066_355_571_25,
            18.820_456_630_672_53,
            -455.5,
            c64(0.066_319_771_490_485_67, 0.696_745_444_703_667_1),
            c64(0.426_933_370_197_611_1, -0.749_970_908_904_107_3),
            1e-15,
        ),
        (
            35.814_904_418_297_68,
            5.870_706_585_831_542,
            -13_409.206_097_240_702,
            c64(0.321_249_329_126_809, 0.030_792_978_506_227_82),
            c64(0.994_977_840_225_617, -0.000_479_923_719_001_871_2),
            1e-15,
        ),
        // The gamma ratios overflow or underflow.
        (1e-310, 1.0, -600.5, c64(0.82, 0.0), c64(-4.071_648_340_424_102_6e86, 0.0), 2e-13),
        (10.0, 400.0, -3000.5, c64(0.351_562_5, 0.0), c64(0.632_731_795_534_758_4, 0.0), 1e-15),
        (
            10.0,
            250.0,
            -3200.25,
            c64(0.25, 0.125),
            c64(0.820_048_634_846_336_2, -0.078_710_402_110_502_63),
            1e-15,
        ),
        (
            1_625.795_720_232_304_8,
            1_464.348_873_264_744_5,
            -1_804.840_375_426_593_3,
            c64(0.043_341_339_953_302_65, 0.000_380_527_842_563_327_6),
            c64(2.103_216_361_252_056_2e-24, -1.030_477_013_920_713_7e-24),
            2e-13,
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        np_assert_allclose!([actual.re, actual.im], [expected.re, expected.im], rtol = rtol);
    }
}

#[test]
fn test_hyp2f1_c64_negative_lower_parameter_componentwise() {
    // Arb at the exact binary64 inputs, componentwise.
    for im in [1.0, -1.0] {
        let z = c64(0.078_125, im * 2.0_f64.powi(-540));
        let actual = xsf::hyp2f1(2.0_f64.powi(-46), 16.0, -300.25, z);
        np_assert_allclose!(
            [actual.re, actual.im],
            [0.999_999_999_999_999_9, im * -2.094_772_308_230_493_7e-178],
            rtol = 1e-14,
            atol = 0.0
        );
    }
    let actual = xsf::hyp2f1(4.0, 1.125, -300.5, c64(0.078_125, 0.0));
    np_assert_allclose!(
        [actual.re, actual.im],
        [0.998_831_693_479_051_2, 0.0],
        rtol = f64::EPSILON,
        atol = 0.0
    );
}

#[test]
fn test_hyp2f1_c64_negative_lower_parameter_slow_tail() {
    // Arb at the exact binary64 inputs, componentwise. The 1 - z series of the connection
    // converges slowly.
    let (a, b, c) = (73.680_102_833_806_94, 117.940_314_118_962_05, -40.141_963_360_724_26);
    let (x, y) = (0.008_982_543_883_583_44, 0.006_264_459_310_953_478);
    let (re, im) = (0.048_352_207_278_304_735, 0.142_312_799_167_885_72);
    for (a, b, c, z, expected, rtol) in [
        (
            73.741_988_378_154_32,
            77.312_817_555_681_2,
            -27.195_599_170_424_025,
            c64(0.010_675_248_312_002_528, 0.0),
            c64(0.160_856_938_961_155_2, 0.0),
            [1e-14, 0.0],
        ),
        (
            158.287_121_924_868_96,
            144.693_188_697_127_54,
            -39.322_275_314_486_54,
            c64(0.004_380_685_242_286_605, 0.0),
            c64(0.087_574_706_153_352_65, 0.0),
            [5e-14, 0.0],
        ),
        (a, b, c, c64(x, -y), c64(re, im), [2e-13, 1e-14]),
        (a, b, c, c64(x, y), c64(re, -im), [2e-13, 1e-14]),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        np_assert_allclose!([actual.re], [expected.re], rtol = rtol[0], atol = 0.0);
        np_assert_allclose!([actual.im], [expected.im], rtol = rtol[1], atol = 0.0);
    }
}

#[test]
fn test_hyp2f1_c64_maclaurin_tiny_a_growth() {
    // Arb at the exact binary64 inputs, componentwise. An early term is negligible, but later
    // ones grow back.
    for (a, b, c, x, expected, rtol) in [
        (1e-20, 700.0, 240.125, 0.5, 1.000_007_275_456_677, 2e-16),
        (1e-20, 700.0, 240.125, 0.45, 1.000_000_000_000_087_5, 2e-16),
        (-1e-20, 700.0, 240.125, 0.5, 0.999_992_724_543_323_1, 2e-16),
        (6e-29, 700.0, 240.125, 0.5, 1.000_000_000_000_043_7, 2e-16),
        (1e-300, 1600.0, 10.125, 0.5, 8.067_451_349_912_805e154, 1e-14),
        (-1.000_000_000_000_003_8, 0.0534, -1.63, 0.8914, 1.029_202_920_245_400_3, 2e-16),
    ] {
        let actual = xsf::hyp2f1(a, b, c, c64(x, 0.0));
        np_assert_allclose!([actual.re], [expected], rtol = rtol, atol = 0.0);
        assert_eq!(actual.im.to_bits(), 0.0_f64.to_bits());
    }
    // For a complex z, the terms cancel strongly.
    for (a, b, c, z, expected, rtol) in [
        (
            1e-20,
            700.0,
            240.125,
            c64(0.5, 0.1),
            c64(1.000_000_000_003_764_8, -5.632_690_310_646_517e-12),
            [2e-16, 1e-8],
        ),
        (1e-16, 100.0, 50.125, c64(0.5, 0.5), c64(1.0, 1.578_452_919_583_920_6e-16), [0.0, 1e-12]),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        np_assert_allclose!([actual.re], [expected.re], rtol = rtol[0], atol = 0.0);
        np_assert_allclose!([actual.im], [expected.im], rtol = rtol[1], atol = 0.0);
    }
    // Arb: 1.119560506102335e168 and 1.2112660936866382e213, but the later terms take too many
    // steps to sum.
    for x in [0.84, 0.875] {
        assert!(xsf::hyp2f1(1e-20, 700.0, 240.125, c64(x, 0.0)).is_nan());
    }
    // Here the later terms cancel to below their rounding errors.
    let (a, b, c) = (6.077_191_929_165_988e-17, 1_109.068_702_298_363_9, 50.284_871_626_278_594);
    let actual = xsf::hyp2f1(a, b, c, c64(0.332_873_860_275_033_1, -0.666_992_080_834_871_5));
    let expected = c64(0.999_999_999_999_999_8, -1.202_633_650_545_082_3e-16);
    assert!((actual - expected).norm() <= 2e-15, "{actual}");
}

#[test]
fn test_hyp2f1_c64_maclaurin_tail_below_head_rounding() {
    // Arb at the exact binary64 inputs, componentwise. The later terms after a small one are
    // below the rounding of the earlier ones.
    let (a, b, c) = (-1.000_001_899_260_690_3, -1.001_935_782_943_842_2, -16.999_997_374_621_41);
    let actual = xsf::hyp2f1(a, b, c, c64(0.122_588_934_300_364_26, 0.301_389_961_311_846_87));
    let expected = c64(0.992_774_912_227_954_1, -0.017_763_176_883_007_648);
    np_assert_allclose!([actual.re], [expected.re], rtol = 2e-16, atol = 0.0);
    np_assert_allclose!([actual.im], [expected.im], rtol = 2e-13, atol = 0.0);
}
