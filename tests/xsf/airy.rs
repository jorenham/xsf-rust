use core::f64::consts::PI;

use num_complex::{Complex64, c64};

fn as_tuple<T: Clone>(x: &[T; 4]) -> (T, T, T, T) {
    (x[0].clone(), x[1].clone(), x[2].clone(), x[3].clone())
}

xsref::tests! {
    test_airy_f64: "airy", "d-d_d_d_d", |x| as_tuple(&xsf::airy(x[0]));
}

#[test]
fn test_airy_f64_large() {
    // Arb: airyai(x), airyai(x, 1), airybi(x), airybi(x, 1)
    #[rustfmt::skip]
    let cases = [
        (20.0, [1.691_672_868_670_540_4e-27, -7.586_391_625_748_354e-27, 2.103_765_049_651_104e25, 9.381_839_336_133_965e25]),
        (50.0, [4.584_941_724_074_828_5e-104, -3.244_331_819_828_799e-103, 4.909_099_699_444_219_5e101, 3.468_798_779_545_976_5e102]),
        (100.0, [2.634_482_152_088_184_6e-291, -2.635_140_361_604_409_7e-290, 6.041_223_996_670_201e288, 6.039_712_745_310_603e289]),
    ];
    for (x, expected) in cases {
        xsref::np_assert_allclose!(&xsf::airy(x), &expected, rtol = 1e-15);
    }
}

#[test]
fn test_airy_f64_overflow() {
    // mpmath: airyai(x), airyai(x, 1), airybi(x), airybi(x, 1)
    xsref::np_assert_allclose!(
        &xsf::airy(103.5),
        &[
            1.214_559_328_776_030_7e-306,
            -1.235_924_520_706_669e-305,
            1.288_045_734_902_866_3e304,
            1.310_081_374_026_67e305,
        ],
        rtol = 1e-13
    );
    let inf = f64::INFINITY;
    let expected = [
        7.448_752_158_292_226e-309,
        -7.598_056_033_156_867e-308,
        2.095_173_527_033_602e306,
        2.136_162_195_043_275_3e307,
    ];
    xsref::np_assert_allclose!(&xsf::airy(104.0), &expected, rtol = 1e-13);
    let [ai, aip, bi, _] = xsf::airy(106.5);
    let expected = [5.3808e-320, -5.5542e-319];
    xsref::np_assert_allclose!([ai, aip], expected, rtol = 1e-4);
    assert_eq!(bi, inf);
    assert_eq!(xsf::airy(107.7), [0.0, 0.0, inf, inf]);
    let [_, _, bi, bip] = xsf::airy(104.4);
    xsref::np_assert_allclose!([bi], [1.241_898_624_243_905_2e308], rtol = 1e-13);
    assert_eq!(bip, inf);
    assert_eq!(xsf::airy(1e9), [0.0, 0.0, inf, inf]);
}

#[test]
fn test_airy_f64_ai_below_asymptotic() {
    // Arb: airyai(x) and airyai(x, 1), where their Maclaurin series cancel
    for (x, expected) in [
        (0.75, [0.179_336_305_478_645_23, -0.193_175_208_104_376_46]),
        (1.0, [0.135_292_416_312_881_42, -0.159_147_441_296_793_2]),
        (1.5, [0.071_749_497_008_105_41, -0.097_382_012_842_301_32]),
        (2.0, [0.034_924_130_423_274_38, -0.053_090_384_433_653_63]),
    ] {
        let [ai, aip, _, _] = xsf::airy(x);
        xsref::np_assert_allclose!([ai, aip], expected, rtol = 1e-15);
    }
}

#[test]
fn test_airy_f64_large_negative() {
    // mpmath: airyai(x), airyai(x, 1), airybi(x), airybi(x, 1)
    let expected = [
        -0.013_152_978_737_498_165,
        -9.130_078_851_932_488,
        0.028_871_844_285_058_44,
        -4.159_337_010_446_404,
    ];
    xsref::np_assert_allclose!(&xsf::airy(-1e5), &expected, rtol = 1e-14);
    let expected = [
        0.005_418_514_944_210_602,
        -26.701_707_811_976_564,
        0.008_443_821_410_258_535,
        17.134_848_759_576_796,
    ];
    xsref::np_assert_allclose!(&xsf::airy(-1e7), &expected, rtol = 1e-14);
    assert!(xsf::airy(-1e17).iter().all(|y| y.is_nan()));
}

#[test]
fn test_airy_c64_near_overflow() {
    // Arb: airyai(z), airyai(z, 1), airybi(z), and airybi(z, 1) for z = -51.8 + 89.7j
    let expected = [
        c64(1.360_456_903_014_940_5e304, -6.086_330_281_445_445e303),
        c64(-1.228_420_803_265_360_1e305, -8.892_381_505_257_792e304),
        c64(6.086_330_281_445_445e303, 1.360_456_903_014_940_5e304),
        c64(8.892_381_505_257_792e304, -1.228_420_803_265_360_1e305),
    ];
    xsref::np_assert_allclose!(&xsf::airy(c64(-51.8, 89.7)), &expected, rtol = 2e-15);
    // Arb: airyai(z) and airybi(z) for z = -52.2 + 90.2j, where AMOS takes Ai' and Bi' to overflow
    let [ai, _, bi, _] = xsf::airy(c64(-52.2, 90.2));
    let expected = [
        c64(7.975_579_981_739_534e306, 5.040_682_424_706_254e306),
        c64(-5.040_682_424_706_254e306, 7.975_579_981_739_534e306),
    ];
    xsref::np_assert_allclose!(&[ai, bi], &expected, rtol = 2e-15);
    // Arb: airyai(z) and airyai(z, 1) for z = 101.5 exp(2 pi i / 3)
    let [ai, aip, _, _] = xsf::airy(Complex64::from_polar(101.5, 2.0 * PI / 3.0));
    let expected = [
        c64(9.011_428_166_917_66e294, -5.202_750_477_955_84e294),
        c64(-9.076_541_392_795_495e295, -5.240_343_616_438_283e295),
    ];
    xsref::np_assert_allclose!(&[ai, aip], &expected, rtol = 1e-12);
}

#[test]
fn test_airy_c64_large_re() {
    // Arb: acb(z).airy()
    let cases = [
        (
            c64(100.0, 0.0),
            [
                c64(2.634_482_152_088_184_6e-291, 0.0),
                c64(-2.635_140_361_604_409_7e-290, 0.0),
                c64(6.041_223_996_670_201e288, 0.0),
                c64(6.039_712_745_310_603e289, 0.0),
            ],
        ),
        (
            c64(50.0, -10.0),
            [
                c64(-2.947_625_110_756_879e-103, 1.523_716_010_365_468e-102),
                c64(1.025_318_671_408_916_6e-102, -1.104_194_741_809_931_5e-101),
                c64(-1.324_972_145_947_649_7e99, -1.430_009_952_035_945_3e100),
                c64(-1.948_467_003_455_968e100, -1.006_139_393_537_108_6e101),
            ],
        ),
        (
            c64(103.5, 0.5),
            [
                c64(4.482_483_837_040_778e-307, 1.136_862_260_721_522_6e-306),
                c64(-4.533_424_502_327_162e-306, -1.157_965_043_710_214e-305),
                c64(4.666_864_451_708_122_5e303, -1.192_055_426_868_804e304),
                c64(4.776_025_158_307_266e304, -1.211_305_065_370_481_7e305),
            ],
        ),
    ];
    for (z, expected) in cases {
        xsref::np_assert_allclose!(&xsf::airy(z), &expected, rtol = 1e-15);
    }
}

#[test]
fn test_airy_c64_large_phase() {
    // Arb: acb(z).airy() and the scaled ones (e^zeta for Ai and Ai', e^-|Re zeta| for Bi and Bi'),
    // for |ph z| > pi/3, and their conjugates
    let cases = [
        (
            c64(-95.8, 14.4),
            [
                c64(-1.641_283_069_088_02e60, -2.984_218_827_895_379_3e59),
                c64(-1.729_171_707_827_661_6e60, 1.632_703_823_389_417_4e61),
                c64(2.984_218_827_895_379_3e59, -1.641_283_069_088_02e60),
                c64(-1.632_703_823_389_417_4e61, -1.729_171_707_827_661_6e60),
            ],
            [
                c64(0.065_902_575_103_500_7, -0.061_174_096_705_638_61),
                c64(-0.648_587_509_119_586_3, -0.602_103_131_196_073_6),
                c64(0.016_085_546_432_400_23, -0.088_468_495_573_251_79),
                c64(-0.880_060_567_810_631_9, -0.093_205_871_955_009_94),
            ],
        ),
        (
            c64(30.0, 80.0),
            [
                c64(-2.390_496_929_241_727_3e54, -8.108_129_472_483_718e54),
                c64(-2.450_298_168_293_716_2e55, 7.418_817_286_333_621e55),
                c64(8.108_129_472_483_718e54, -2.390_496_929_241_727_3e54),
                c64(-7.418_817_286_333_621e55, -2.450_298_168_293_716_2e55),
            ],
            [
                c64(0.088_564_846_161_249_72, -0.027_675_892_243_092_602),
                c64(-0.818_588_242_770_607_3, -0.255_757_090_837_983_37),
                c64(0.089_000_875_527_667_57, -0.026_239_876_949_515_72),
                c64(-0.814_344_709_349_130_4, -0.268_962_999_446_502_54),
            ],
        ),
        (
            c64(-600.0, 0.5),
            [
                c64(-1_357.798_706_186_429_9, 11_803.262_113_559_258),
                c64(289_133.008_231_432_06, 33_143.593_967_670_23),
                c64(-11_803.262_114_102_576, -1_357.798_706_124_158_1),
                c64(-33_143.593_969_201_34, 289_133.008_218_124_2),
            ],
            [
                c64(0.040_311_560_072_121_576, -0.040_295_338_283_865_154),
                c64(-0.987_424_882_168_184, -0.987_033_128_763_382_6),
                c64(-0.056_624_255_356_529_88, -0.006_513_821_341_515_287),
                c64(-0.159_001_071_924_250_14, 1.387_069_195_877_931_8),
            ],
        ),
    ];
    for (z, expected, scaled) in cases {
        for (z, expected, scaled) in [
            (z, expected, scaled),
            (z.conj(), expected.map(|w| w.conj()), scaled.map(|w| w.conj())),
        ] {
            xsref::np_assert_allclose!(&xsf::airy(z), &expected, rtol = 2e-15);
            xsref::np_assert_allclose!(&xsf::airy_scaled(z), &scaled, rtol = 2e-15);
        }
    }
}

#[test]
fn test_airy_c64_negative_real() {
    // Arb: acb(z).airy(); Bi and Bi' are not scaled, as Re zeta = 0
    let expected = [
        c64(0.100_809_025_019_288_57, 0.0),
        c64(-1.143_627_547_203_056_6, 0.0),
        c64(0.043_226_419_093_977_925, 0.0),
        c64(2.667_171_540_341_03, 0.0),
    ];
    let z = c64(-700.0, 0.0);
    xsref::np_assert_allclose!(&xsf::airy(z), &expected, rtol = 2e-15);
    xsref::np_assert_allclose!(&xsf::airy_scaled(z)[2..], &expected[2..], rtol = 2e-15);
}

#[test]
fn test_airy_c64_stokes_line() {
    // Arb: acb(z).airy(), and Bi e^-|Re zeta|; at ph z = pi/3, where the e^-zeta part of Bi is as
    // large as its e^zeta part, and its conjugate
    let z = c64(50.000_000_000_000_014, 86.602_540_378_443_86);
    let expected = [
        c64(0.054_699_249_999_691_67, -0.070_467_992_455_229_96),
        c64(-0.825_965_298_485_187_2, 0.336_981_030_968_697_7),
        c64(0.235_677_865_045_882_34, 0.122_054_143_239_856_78),
        c64(0.756_651_800_328_036_9, 0.583_668_266_825_105_8),
    ];
    let scaled = [
        c64(0.235_677_865_045_850_56, 0.122_054_143_239_840_32),
        c64(0.756_651_800_327_934_9, 0.583_668_266_825_027),
    ];
    for (z, expected, scaled) in
        [(z, expected, scaled), (z.conj(), expected.map(|w| w.conj()), scaled.map(|w| w.conj()))]
    {
        xsref::np_assert_allclose!(&xsf::airy(z), &expected, rtol = 1e-15);
        xsref::np_assert_allclose!(&xsf::airy_scaled(z)[2..], &scaled, rtol = 1e-15);
    }
}

#[test]
fn test_airy_c64_near_real_axis() {
    // Arb: Im Bi and Im Bi' (also scaled) by component, as their size follows Im z
    for sign in [1.0, -1.0] {
        let z = c64(1.5, sign * 1e-20);
        let [.., bi, bip] = xsf::airy(z);
        let [.., ebi, ebip] = xsf::airy_scaled(z);
        let actual = [bi.im, bip.im, ebi.im, ebip.im];
        let expected = [
            1.886_212_254_848_165_5e-20,
            2.818_412_255_621_842_4e-20,
            5.542_307_563_918_05e-21,
            8.281_415_584_286_761e-21,
        ];
        xsref::np_assert_allclose!(&actual, &expected.map(|w| sign * w), rtol = 2e-15);
    }
}

#[test]
fn test_airy_c64_scaled_zeta_re_below_zero() {
    // Arb: Bi e^-|Re zeta| and Bi' e^-|Re zeta|, at ph z just above pi/3, where Re zeta = -5.5e-15
    // but its rounded value is 0
    let expected = [
        c64(0.054_250_139_932_712_535, 0.139_399_474_750_160_3),
        c64(-1.733_613_607_948_874_9, 1.348_329_937_932_414_5),
    ];
    for (z, expected) in [
        (c64(30.0, 51.961_524_227_066_32), expected),
        (c64(30.0, -51.961_524_227_066_32), expected.map(|w| w.conj())),
    ] {
        xsref::np_assert_allclose!(&xsf::airy_scaled(z)[2..], &expected, rtol = 1e-15);
    }
}

#[test]
fn test_airy_c64_scaled_zeta_re_above_zero() {
    // Arb: Bi e^-|Re zeta| and Bi' e^-|Re zeta|, at ph z just below pi/3, where Re zeta = 9.7e-14
    // but its rounded value is below 0
    let expected = [
        c64(-0.005_327_918_203_229_043, 0.070_112_237_908_936_32),
        c64(-3.082_427_605_065_668_3, 2.005_886_794_538_590_8),
    ];
    for (z, expected) in [
        (c64(186.5, 323.027_475_611_595_6), expected),
        (c64(186.5, -323.027_475_611_595_6), expected.map(|w| w.conj())),
    ] {
        xsref::np_assert_allclose!(&xsf::airy_scaled(z)[2..], &expected, rtol = 1e-15);
    }
}

#[test]
fn test_airy_c64_subnormal() {
    // Arb: acb(z).airy(); Ai'(z) = (3 + 3i) 2^-1074, with e^-zeta below it
    let z = c64(500_000.645_361_992_06, 866_025.031_184_198_1);
    let aip = xsf::airy(z)[1];
    assert_eq!((aip.re, aip.im), (1.5e-323, 1.5e-323));
}

#[test]
fn test_airy_c64_signed_zero_im() {
    // the same bits for real `z` with +0 and -0 as imaginary part, as AMOS
    for x in [2.5, 100.0] {
        let bits = |im: f64| xsf::airy(c64(x, im)).map(|w| (w.re.to_bits(), w.im.to_bits()));
        assert_eq!(bits(0.0), bits(-0.0));
    }
}

#[test]
fn test_airye_c64_large_re() {
    // Arb: acb(z).airy()[2:] e^-|Re zeta|
    let [_, _, bi, bip] = xsf::airy_scaled(c64(50.0, -10.0));
    let expected = [
        c64(-0.019_484_259_483_357_722, -0.210_288_835_538_684_27),
        c64(-0.286_530_073_905_391_7, -1.479_569_293_593_267_2),
    ];
    xsref::np_assert_allclose!(&[bi, bip], &expected, rtol = 1e-15);
}

#[test]
fn test_airye_f64_cephes() {
    // Arb: airyai(x) e^zeta, airyai(x, 1) e^zeta, airybi(x) e^-zeta and airybi(x, 1) e^-zeta
    // for zeta = 2/3 x^(3/2)
    #[rustfmt::skip]
    let cases = [
        // AMOS cancels for Ai and Ai'
        (1.5, [0.244_184_897_671_408_43, -0.331_419_979_686_364, 0.552_094_372_285_784_1, 0.554_230_756_391_805]),
        (2.0, [0.230_164_918_652_511_6, -0.349_888_282_580_087_5, 0.500_437_254_304_095, 0.622_217_997_315_437_6]),
        (5.0, [0.187_002_118_935_943_43, -0.427_035_544_351_945_2, 0.381_108_531_088_877_4, 0.831_878_259_124_801_4]),
        // AMOS fails for x > 1.05e6
        (2e6, [0.007_501_317_867_926_573, -10.608_465_465_431_054, 0.015_002_635_736_958_195, 21.216_930_928_674_223]),
    ];
    for (x, expected) in cases {
        xsref::np_assert_allclose!(&xsf::airy_scaled(x), &expected, rtol = 2e-15);
    }
}

xsref::tests! {
    test_airy_c64: "airy", "cd-cd_cd_cd_cd", |x| as_tuple(&xsf::airy(c64(x[0], x[1])));
    test_airye_f64: "airye", "d-d_d_d_d", |x| as_tuple(&xsf::airy_scaled(x[0]));
    test_airye_c64: "airye", "cd-cd_cd_cd_cd", |x| as_tuple(&xsf::airy_scaled(c64(x[0], x[1])));
    test_itairy: "itairy", "d-d_d_d_d", |x| as_tuple(&xsf::airy_integrals(x[0]));
}

#[test]
fn test_itairy_overflow() {
    // Wolfram: Integrate[AiryBi[t], {t, 0, x}]
    let bpt = [104.3, 104.5, 104.663].map(|x| xsf::airy_integrals(x)[1]);
    let expected =
        [4.382_426_630_792_426e306, 3.377_394_984_522_598_2e307, 1.786_520_403_853_323_6e308];
    xsref::np_assert_allclose!(bpt, expected, rtol = 1e-13);
    assert_eq!(xsf::airy_integrals(104.67)[1], f64::INFINITY);
}

#[test]
fn test_itairy_inf() {
    // the limits as x -> inf
    let [a, b, c, d] = [1.0 / 3.0, f64::INFINITY, 2.0 / 3.0, 0.0];
    for x in [1e205, 1e300, f64::INFINITY] {
        assert_eq!(xsf::airy_integrals(x), [a, b, c, d]);
        assert_eq!(xsf::airy_integrals(-x), [-c, -d, -a, -b]);
    }
    assert!(xsf::airy_integrals(f64::NAN).iter().all(|y| y.is_nan()));
}

#[test]
fn test_itairy_2_20() {
    // mpmath: `_airy_int` in tools/accuracy/specs.py
    let cases = [
        (
            3.0,
            [
                0.329_920_376_007_021_8,
                10.062_003_031_306_201,
                0.801_462_842_671_323_4,
                -0.195_442_489_513_706_02,
            ],
        ),
        (
            9.0,
            [
                0.333_333_332_530_663_66,
                7_377_921.170_522_002,
                0.558_819_748_944_718_8,
                -0.010_330_373_551_680_188,
            ],
        ),
        (
            15.0,
            [
                0.333_333_333_333_333_3,
                4_967_308_973_718_509.0,
                0.683_580_628_361_935_9,
                0.072_026_214_312_908_51,
            ],
        ),
    ];
    for (x, expected) in cases {
        xsref::np_assert_allclose!(xsf::airy_integrals(x), expected, rtol = 2e-15);
    }
}

#[test]
fn test_ai_zeros_0() {
    assert!(xsf::airy_ai_zeros::<0>().iter().all(|a| a.is_empty()));
}

#[test]
fn test_ai_zeros() {
    // >>> from scipy import special
    // >>> a, ap, ai, aip = special.ai_zeros(3)
    let [a, ap, ai, aip] = xsf::airy_ai_zeros::<3>();
    // >>> a
    // array([-2.33810741, -4.08794944, -5.52055983])
    xsref::np_assert_allclose!(&a, &[-2.338_107_41, -4.087_949_44, -5.520_559_83]);
    // >>> ap
    // array([-1.01879297, -3.24819758, -4.82009921])
    xsref::np_assert_allclose!(&ap, &[-1.018_792_97, -3.248_197_58, -4.820_099_21]);
    // >>> ai
    // array([ 0.53565666, -0.41901548,  0.38040647])
    xsref::np_assert_allclose!(&ai, &[0.535_656_66, -0.419_015_48, 0.380_406_47]);
    // >>> aip
    // array([ 0.70121082, -0.80311137,  0.86520403])
    xsref::np_assert_allclose!(&aip, &[0.701_210_82, -0.803_111_37, 0.865_204_03]);
}

#[test]
fn test_ai_zeros_5() {
    // the 5th zeros of Ai and Ai', and Ai and Ai' there (Arb)
    let [a, ap, ai, aip] = xsf::airy_ai_zeros::<5>();
    xsref::np_assert_allclose!(
        [a[4], ap[4], ai[4], aip[4]],
        [
            -7.944_133_587_120_853,
            -7.372_177_255_047_77,
            0.342_301_244_411_623_8,
            0.947_335_709_441_567_8,
        ],
        rtol = 2e-15
    );
}

#[test]
fn test_bi_zeros_0() {
    assert!(xsf::airy_bi_zeros::<0>().iter().all(|a| a.is_empty()));
}

#[test]
fn test_bi_zeros() {
    // >>> from scipy import special
    // >>> b, bp, bi, bip = special.bi_zeros(3)
    let [b, bp, bi, bip] = xsf::airy_bi_zeros::<3>();
    // >>> b
    // array([-1.17371322, -3.2710933 , -4.83073784])
    xsref::np_assert_allclose!(&b, &[-1.173_713_22, -3.271_093_3, -4.830_737_84]);
    // >>> bp
    // array([-2.29443968, -4.07315509, -5.51239573])
    xsref::np_assert_allclose!(&bp, &[-2.294_439_68, -4.073_155_09, -5.512_395_73]);
    // >>> bi
    // array([-0.45494438,  0.39652284, -0.36796916])
    xsref::np_assert_allclose!(&bi, &[-0.454_944_38, 0.396_522_84, -0.367_969_16]);
    // >>> bip
    // array([ 0.60195789, -0.76031014,  0.83699101])
    xsref::np_assert_allclose!(&bip, &[0.601_957_89, -0.760_310_14, 0.836_991_01]);
}

#[test]
fn test_bi_zeros_5() {
    // the 5th zeros of Bi and Bi', and Bi and Bi' there (Arb)
    let [b, bp, bi, bip] = xsf::airy_bi_zeros::<5>();
    xsref::np_assert_allclose!(
        [b[4], bp[4], bi[4], bip[4]],
        [
            -7.376_762_079_367_764,
            -7.940_178_689_168_579,
            -0.336_026_240_133_660_04,
            0.929_983_638_568_026_6,
        ],
        rtol = 2e-15
    );
}
