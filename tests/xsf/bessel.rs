#![allow(clippy::needless_range_loop)]

use num_complex::c64;
use xsref::np_assert_allclose;

xsref::tests! {
    // Bessel J
    test_bessel_j0: "cyl_bessel_j0", "d-d", |x| xsf::bessel_j0(x[0]);
    test_bessel_j1: "cyl_bessel_j1", "d-d", |x| xsf::bessel_j1(x[0]);
    test_bessel_j_f64: "cyl_bessel_j", "d_d-d", |x| xsf::bessel_j(x[0], x[1]);
    test_bessel_j_c64: "cyl_bessel_j", "d_cd-cd", |x| xsf::bessel_j(x[0], c64(x[1], x[2]));
    test_bessel_je_f64: "cyl_bessel_je", "d_d-d", |x| xsf::bessel_je(x[0], x[1]);
    test_bessel_je_c64: "cyl_bessel_je", "d_cd-cd", |x| xsf::bessel_je(x[0], c64(x[1], x[2]));

    // Bessel Y
    test_bessel_y0: "cyl_bessel_y0", "d-d", |x| xsf::bessel_y0(x[0]);
    test_bessel_y1: "cyl_bessel_y1", "d-d", |x| xsf::bessel_y1(x[0]);
    test_bessel_y_f64: "cyl_bessel_y", "d_d-d", |x| xsf::bessel_y(x[0], x[1]);
    test_bessel_y_c64: "cyl_bessel_y", "d_cd-cd", |x| xsf::bessel_y(x[0], c64(x[1], x[2]));
    test_bessel_ye_f64: "cyl_bessel_ye", "d_d-d", |x| xsf::bessel_ye(x[0], x[1]);
    test_bessel_ye_c64: "cyl_bessel_ye", "d_cd-cd", |x| xsf::bessel_ye(x[0], c64(x[1], x[2]));
}

#[test]
fn test_bessel_y_negative_half_integer_tiny_x() {
    // mpmath: bessely(-2.5, 1e-200) = 5.3e-502 and bessely(-50.5, 1e-10) = 2.9e-586
    for (v, x) in [(-2.5, 1e-200), (-50.5, 1e-10)] {
        assert_eq!(xsf::bessel_y(v, x), 0.0);
        assert_eq!(xsf::bessel_ye(v, x), 0.0);
        assert_eq!(xsf::bessel_y(v, c64(x, 0.0)), c64(0.0, 0.0));
        assert_eq!(xsf::bessel_ye(v, c64(x, 0.0)), c64(0.0, 0.0));
    }
}

xsref::tests! {
    // Bessel I
    test_bessel_i0: "cyl_bessel_i0", "d-d", |x| xsf::bessel_i0(x[0]);
    test_bessel_i0e: "cyl_bessel_i0e", "d-d", |x| xsf::bessel_i0e(x[0]);
    test_bessel_i1: "cyl_bessel_i1", "d-d", |x| xsf::bessel_i1(x[0]);
    test_bessel_i1e: "cyl_bessel_i1e", "d-d", |x| xsf::bessel_i1e(x[0]);
    test_bessel_i_f64: "cyl_bessel_i", "d_d-d", |x| xsf::bessel_i(x[0], x[1]);
    test_bessel_i_c64: "cyl_bessel_i", "d_cd-cd", |x| xsf::bessel_i(x[0], c64(x[1], x[2]));
    test_bessel_ie_f64: "cyl_bessel_ie", "d_d-d", |x| xsf::bessel_ie(x[0], x[1]);
    test_bessel_ie_c64: "cyl_bessel_ie", "d_cd-cd", |x| xsf::bessel_ie(x[0], c64(x[1], x[2]));
}

#[test]
fn test_bessel_i_near_negative_integer() {
    // mpmath: besseli(-2.0000000000000004, 1e-150) = 1.7763568394005238e285
    let (v, x) = (-2.000_000_000_000_000_4, 1e-150);
    let actual = [xsf::bessel_ie(v, x), xsf::bessel_i(v, c64(x, 0.0)).re];
    let expected = [1.776_356_839_400_523_8e285; 2];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_bessel_i_small_x_k_overflow() {
    // mpmath: besseli(46.408, 9.696e-6) = 8.8291863067264228e-306,
    // besseli(-40.00000109418565, 4.4076818963921526e-7) = 4.1828779773220576e306, and
    // besseli(-48.000001, 9.813429518583089e-6) = 1.7976931348622338e308
    let actual = [
        xsf::bessel_i(46.408, 9.696e-6),
        xsf::bessel_i(-40.000_001_094_185_65, 4.407_681_896_392_152_6e-7),
        xsf::bessel_i(-48.000_001, 9.813_429_518_583_089e-6),
    ];
    let expected = [
        8.829_186_306_726_423e-306,
        4.182_877_977_322_057_6e306,
        1.797_693_134_862_233_8e308,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);

    // mpmath: besseli(-3.000053769783958, 0.5818348207903663) = 1.8e-19 (cancellation, not overflow)
    let i = xsf::bessel_i(-3.000_053_769_783_958, 0.581_834_820_790_366_3);
    assert!(i.abs() < 3e-19, "{i}");
}

#[test]
fn test_bessel_i_large_v_vanishing_term() {
    // Arb: I_v(x) where u_1(t) = 0 (x = |v| sqrt(2/3)) or u_2(t) = 0 (x = 1.921 |v| or
    // x = 0.1157 |v|, where K_v dominates I_{-v}) in the uniform asymptotic expansion
    let actual = [
        xsf::bessel_i(100.0, 81.649_658_092_772_6),
        xsf::bessel_i(-100.5, 82.057_906_383_236_47),
        xsf::bessel_i(1000.0, 816.496_580_927_726_1),
        xsf::bessel_i(100.0, 192.102_108_828_307_14),
        xsf::bessel_i(-70.5, 8.155_385_051_326_476),
    ];
    let expected = [
        6_392_289_838.263_86,
        7_258_961_449.952_213,
        4.441_674_522_431_231e110,
        6.238_986_295_989_078e70,
        3.311_216_933_797_801_6e55,
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}

#[test]
fn test_bessel_amos_underflow() {
    // mpmath: besselj(7.3156, 1.19e-41) = 2.6872935555094642e-306,
    // besselj(-43, 4.418464059606903e-6) = -1.0500997742421043e-296,
    // besselj(5, -2.9158810273573306e-60) = -5.4893102343051316e-302,
    // besseli(46.408, 9.696e-6) exp(-9.696e-6) = 8.829100699351018e-306, and
    // besseli(155.31716506376495, 1.2644828833671058) exp(-1.2644828833671058) = 1.4129266246231237e-306
    let actual = [
        xsf::bessel_j(7.3156, 1.19e-41),
        xsf::bessel_je(7.3156, 1.19e-41),
        xsf::bessel_j(-43.0, 4.418_464_059_606_903e-6),
        xsf::bessel_j(5.0, -2.915_881_027_357_330_6e-60),
        xsf::bessel_ie(46.408, 9.696e-6),
        xsf::bessel_ie(155.317_165_063_764_95, 1.264_482_883_367_105_8),
    ];
    let expected = [
        2.687_293_555_509_464_2e-306,
        2.687_293_555_509_464_2e-306,
        -1.050_099_774_242_104_3e-296,
        -5.489_310_234_305_132e-302,
        8.829_100_699_351_018e-306,
        1.412_926_624_623_123_7e-306,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-12);

    // an AMOS 0 from cancellation (not underflow), where cephes is NaN
    assert!(xsf::bessel_j(-782.299_556_217_053_6, 794.671_828_363_666_9).is_finite());
    assert!(xsf::bessel_je(-1_744.559_130_346_428, 1_754.183_596_660_477_4).is_finite());
}

#[test]
fn test_bessel_scaled_overflow() {
    // the signs of mpmath's besselj, bessely, and besseli, e.g. besselj(-11.25, 1e-100) = -3.6e1134
    let x = 1e-100;
    let z = c64(x, 0.0);
    let inf = f64::INFINITY;
    assert_eq!(xsf::bessel_ye(10.0, x), -inf);
    assert_eq!(xsf::bessel_ye(10.0, z), c64(-inf, 0.0));
    // Arb: besselj(-1000.25, 300) = 9.2e397
    assert_eq!(xsf::bessel_j(-1000.25, c64(300.0, 0.0)), c64(inf, 0.0));
    for (v, s) in [(-10.25, 1.0), (-11.25, -1.0)] {
        assert_eq!(xsf::bessel_j(v, z), c64(s * inf, 0.0));
        assert_eq!(xsf::bessel_je(v, x), s * inf);
        assert_eq!(xsf::bessel_ye(v, x), -s * inf);
        assert_eq!(xsf::bessel_ie(v, x), s * inf);
        assert_eq!(xsf::bessel_je(v, z), c64(s * inf, 0.0));
        assert_eq!(xsf::bessel_ye(v, z), c64(-s * inf, 0.0));
        assert_eq!(xsf::bessel_ie(v, z), c64(s * inf, 0.0));
    }
}

#[test]
fn test_bessel_jyi_near_overflow() {
    // mpmath: bessely(v, x), besselj(v, x), and besseli(v, x) * exp(-x)
    let actual = [
        xsf::bessel_y(79.762_372_861_572_26, 8.365_870_974_604_77e-3),
        xsf::bessel_ye(-61.536_079_057_536_22, 4.604_376_665_618_018_6e-4),
        xsf::bessel_je(-4.021_221_847_884_007, 5.996_317_823_411_197_6e-77),
        xsf::bessel_ie(-3.819_271_259_116_135_3, 5.989_105_622_767_614e-81),
        xsf::bessel_ie(-2.000_000_000_000_000_4, 4.214_684_851_089_752e-162),
    ];
    let expected = [
        -5.240_453_974_629_458_6e305,
        -1.960_111_897_579_981e305,
        6.803_954_046_133_499e306,
        -2.862_698_759_482_083_5e307,
        1e308,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-12);
    // mpmath: 1.9e308 and 1.99e-37 (exp(-2x) underflows, so no retry)
    for (v, x) in [
        (-2.000_000_000_000_000_4, 3.057_653_761_164_066_6e-162),
        (-840.25, 400.0),
    ] {
        assert_eq!(xsf::bessel_ie(v, x), f64::INFINITY);
    }
    // Arb: 1.4501132090238621e-21
    xsref::np_assert_allclose!(
        &[xsf::bessel_ie(-804.25, 373.0)],
        &[1.450_113_209_023_862_1e-21],
        rtol = 1e-10
    );
}

#[test]
fn test_bessel_huge_order() {
    // Debye expansions (DLMF 10.19.6) at 60 digits, and J(-x) = e^(i pi v) J(x) (DLMF 10.11.1)
    let actual = [
        xsf::bessel_j(2_147_483_647.0, c64(2_149_631_130.647, 0.0)),
        xsf::bessel_j(8_804_682_956.5, c64(-8_813_487_639.456_5, 0.0)),
        xsf::hankel_1(8_804_682_957.25, c64(-8_813_487_640.207_249, 0.0)),
        xsf::bessel_y(36_099_200_123.25, c64(36_135_299_323.373_245, 0.0)),
    ];
    let expected = [
        c64(-1.632_029_240_809_612_2e-5, 0.0),
        c64(0.0, 3.866_168_657_892_145e-5),
        c64(1.953_806_271_557_623_4e-5, -3.513_747_043_055_447_5e-5),
        c64(-8.323_868_312_936_047e-6, 0.0),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-6);
    let j = xsf::bessel_j(2_147_483_648.0, 2_149_631_131.648);
    np_assert_allclose!(&[j], &[-1.632_266_889_229_664_4e-5], rtol = 1e-6);
}

xsref::tests! {
    // Bessel K
    test_bessel_k0: "cyl_bessel_k0", "d-d", |x| xsf::bessel_k0(x[0]);
    test_bessel_k0e: "cyl_bessel_k0e", "d-d", |x| xsf::bessel_k0e(x[0]);
    test_bessel_k1: "cyl_bessel_k1", "d-d", |x| xsf::bessel_k1(x[0]);
    test_bessel_k1e: "cyl_bessel_k1e", "d-d", |x| xsf::bessel_k1e(x[0]);
    test_bessel_k_f64: "cyl_bessel_k", "d_d-d", |x| xsf::bessel_k(x[0], x[1]);
}

#[test]
fn test_bessel_k_near_overflow() {
    // mpmath: besselk(v, x) and besselk(v, x) * exp(x)
    let actual = [
        xsf::bessel_k(1.5, 1e-205),
        xsf::bessel_k(27.344_993_256_619_134, 1.443_805_244_850_114e-10),
        xsf::bessel_ke(1.280_346_742_634_130_8, 2.382_283_514_915_581_2e-239),
        xsf::bessel_ke(46.418_797_513_730_52, 8.621_169_724_258_798e-6),
    ];
    let expected = [
        3.963_327_297_606_011e307,
        1.307_413_220_213_715e304,
        3.623_615_504_892_954e305,
        3.391_349_164_908_478_2e305,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
    let actual = [xsf::bessel_k(115.0, 0.178), xsf::bessel_ke(115.0, 0.178)];
    let expected = [8.404_870_096_416_689e306, 1.004_235_161_288_784_1e307];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-12);
    for (v, x) in [(1.5, 1e-210), (1.0, 1e-310), (-1.0, 1e-310), (0.99, 5e-324)] {
        assert_eq!(xsf::bessel_k(v, x), f64::INFINITY);
        assert_eq!(xsf::bessel_ke(v, x), f64::INFINITY);
    }
}

#[test]
fn test_bessel_amos_tiny_z() {
    // Arb: besselk(0, 1e-306), besselj(-0.25, 5e-324), bessely(0.75, -1e-310 + 1e-310j), and
    // besselj(1e-20, 5e-324) = 1 - 7.4e-18
    let actual = [
        xsf::bessel_k(0.0, c64(1e-306, 0.0)),
        xsf::bessel_j(-0.25, c64(5e-324, 0.0)),
        xsf::bessel_y(0.75, c64(-1e-310, 1e-310)),
        xsf::bessel_j(1e-20, c64(5e-324, 0.0)),
    ];
    let expected = [
        c64(704.706_969_971_836_4, 0.0),
        c64(6.509_198_852_597_573e80, 0.0),
        c64(3.120_727_850_251_321e231, 1.568_895_836_574_547e232),
        c64(1.0, 0.0),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-13);
}

#[test]
fn test_bessel_amos_k_near_overflow() {
    // Arb: besselk(1.001, 1e-305), bessely(2.5, 2.7e-122), and besselk(0.9535, 5e-324 + 5e-324j)
    let actual = [
        xsf::bessel_k(1.001, c64(1e-305, 0.0)),
        xsf::bessel_y(2.5, c64(2.7e-122, 0.0)),
        xsf::bessel_k(0.9535, c64(5e-324, 5e-324)),
    ];
    let expected = [
        c64(2.018_602_028_941_825e305, 0.0),
        c64(-1.998_259_589_560_272_5e304, 0.0),
        c64(9.821_768_931_473_591e307, -9.129_347_072_678_403e307),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-13);
}

#[test]
fn test_bessel_amos_i_near_overflow() {
    // Arb: besseli(0, 705), besseli(50.5, 711), besselj(2.5, 100 + 704j), bessely(0, -700j), and
    // besselk(0.5, -702 + 10j)
    let actual = [
        xsf::bessel_i(0.0, c64(705.0, 0.0)),
        xsf::bessel_i(50.5, c64(711.0, 0.0)),
        xsf::bessel_j(2.5, c64(100.0, 704.0)),
        xsf::bessel_y(0.0, c64(0.0, -700.0)),
        xsf::bessel_k(0.5, c64(-702.0, 10.0)),
    ];
    let expected = [
        c64(2.262_050_552_655_472_5e304, 0.0),
        c64(1.511_240_519_293_431e306, 0.0),
        c64(-1.513_317_749_321_898_5e303, -8.110_033_101_279_623e303),
        c64(-2.972_872_008_947_041e-306, -1.529_593_347_671_873_7e302),
        c64(1.907_231_153_635_185e303, 2.988_026_562_775_276_6e303),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-13);
    // Arb: besseli(0, 720) = 7.3e310
    assert_eq!(xsf::bessel_i(0.0, c64(720.0, 0.0)), c64(f64::INFINITY, 0.0));
}

#[test]
fn test_bessel_amos_subnormal_order() {
    // Arb: besselj(-5e-324, 1e-200) = 1 + 4.6e-321, bessely(5e-324, 0.5 + 0.5j),
    // hankel1(1e-320, -1 + 1j), and besseli(-5e-324, 2)
    let actual = [
        xsf::bessel_j(-5e-324, c64(1e-200, 0.0)),
        xsf::bessel_y(5e-324, c64(0.5, 0.5)),
        xsf::hankel_1(1e-320, c64(-1.0, 1.0)),
        xsf::bessel_i(-5e-324, c64(2.0, 0.0)),
    ];
    let expected = [
        c64(1.0, 0.0),
        c64(-0.227_087_358_059_979_5, 0.614_350_253_501_374_8),
        c64(-0.227_449_894_802_294_76, -0.051_055_458_673_089_62),
        c64(2.279_585_302_336_067_3, 0.0),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-15);
    // Arb: besselk(-1e-310, 1)
    let k = xsf::bessel_k(-1e-310, 1.0);
    np_assert_allclose!(&[k], &[0.421_024_438_240_708_33], rtol = 1e-15);
}

#[test]
fn test_bessel_amos_i_near_max() {
    // Arb: besseli(0, 713.9), besseli(1.5, 712 - 30j), besselj(0, 712j), bessely(0, -712j), and
    // besselk(0, -712 + 0.5j)
    let actual = [
        xsf::bessel_i(0.0, c64(713.9, 0.0)),
        xsf::bessel_i(1.5, c64(712.0, -30.0)),
        xsf::bessel_j(0.0, c64(0.0, 712.0)),
        xsf::bessel_y(0.0, c64(0.0, -712.0)),
        xsf::bessel_k(0.0, c64(-712.0, 0.5)),
    ];
    let expected = [
        c64(1.648_155_186_695_137_8e308, 0.0),
        c64(3.288_042_717_302_672_6e306, 2.441_382_590_454_382e307),
        c64(2.468_411_057_762_752_4e307, 0.0),
        c64(-1.811_143_065_969_3e-311, -2.468_411_057_762_752_4e307),
        c64(-3.715_430_308_432_895e307, -6.806_731_004_735_646e307),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-13);
    // Arb: besseli(0, 714.5) = 3.0e308
    assert_eq!(xsf::bessel_i(0.0, c64(714.5, 0.0)), c64(f64::INFINITY, 0.0));
}

#[test]
fn test_bessel_amos_large_order_near_overflow() {
    // Arb: besselk(100, 0.0627), besselk(86.25, -706), hankel1(86.25, -706j), bessely(100.25,
    // -706j), and besselj(-300.75, -760j)
    let actual = [
        xsf::bessel_k(100.0, c64(0.0627, 0.0)),
        xsf::bessel_k(86.25, c64(-706.0, 0.0)),
        xsf::hankel_1(86.25, c64(0.0, -706.0)),
        xsf::bessel_y(100.25, c64(0.0, -706.0)),
        xsf::bessel_j(-300.75, c64(0.0, -760.0)),
    ];
    let expected = [
        c64(1.109_723_622_584_384_7e306, 0.0),
        c64(0.0, -9.972_378_534_887_147e302),
        c64(-5.865_353_936_519_697e302, 2.429_509_148_624_880_3e302),
        c64(-1.919_476_381_292_648_6e301, -4.634_025_912_371_543e301),
        c64(1.863_353_773_515_49e302, 4.498_533_951_520_18e302),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-12);
}

#[test]
fn test_bessel_amos_huge_z_airy_failure() {
    // AMOS's Airy functions fail for |z| > ~7e8, so these are NaN instead of 0.
    // Arb: hankel2(86.5, 1e15 + 700j) = -2.2e296 + 1.3e296j and bessely(86.5, 1e10 + 10j) = 0.077 +
    // 0.043j
    assert!(xsf::hankel_2(86.5, c64(1e15, 700.0)).re.is_nan());
    assert!(xsf::bessel_y(86.5, c64(1e10, 10.0)).re.is_nan());
    // Arb: bessely(-86.5, 1e15 + 10j) = sin(86.5 pi) besselj(86.5, 1e15 + 10j)
    let y = xsf::bessel_y(-86.5, c64(1e15, 10.0));
    np_assert_allclose!(
        &[y],
        &[c64(-2.384_955_233_412_051_6e-4, 1.426_054_857_376_210_6e-4)],
        rtol = 1e-13
    );
}

xsref::tests! {
    test_bessel_k_c64: "cyl_bessel_k", "d_cd-cd", |x| xsf::bessel_k(x[0], c64(x[1], x[2]));
    test_bessel_ke_f64: "cyl_bessel_ke", "d_d-d", |x| xsf::bessel_ke(x[0], x[1]);
    test_bessel_ke_c64: "cyl_bessel_ke", "d_cd-cd", |x| xsf::bessel_ke(x[0], c64(x[1], x[2]));

    // Hankel
    test_hankel_1_c64: "cyl_hankel_1", "d_cd-cd", |x| xsf::hankel_1(x[0], c64(x[1], x[2]));
    test_hankel_1e_c64: "cyl_hankel_1e", "d_cd-cd", |x| xsf::hankel_1e(x[0], c64(x[1], x[2]));
    test_hankel_2_c64: "cyl_hankel_2", "d_cd-cd", |x| xsf::hankel_2(x[0], c64(x[1], x[2]));
    test_hankel_2e_c64: "cyl_hankel_2e", "d_cd-cd", |x| xsf::hankel_2e(x[0], c64(x[1], x[2]));
}

// Bessel zeros

/// Based on `scipy.special.tests.test_basic.TestBessel.test_jnyn_zeros`
#[test]
fn test_bessel_zeros() {
    let jnz = xsf::bessel_zeros::<5>(1);
    xsref::np_assert_allclose!(
        jnz[0].as_ref(),
        &[3.83171, 7.01559, 10.17347, 13.32369, 16.47063],
        atol = 1.5e-5
    );
    xsref::np_assert_allclose!(
        jnz[1].as_ref(),
        &[1.84118, 5.33144, 8.53632, 11.70600, 14.86359],
        atol = 1.5e-5
    );
    xsref::np_assert_allclose!(
        jnz[2].as_ref(),
        &[2.19714, 5.42968, 8.59601, 11.74915, 14.89744],
        atol = 1.5e-5
    );
    xsref::np_assert_allclose!(
        jnz[3].as_ref(),
        &[3.68302, 6.94150, 10.12340, 13.28576, 16.44006],
        atol = 1.5e-5
    );
}

xsref::tests! {
    // Bessel integrals
    test_it1j0y0: "it1j0y0", "d-d_d", |x| xsf::it1j0y0(x[0]);
}

#[test]
fn test_it2j0y0() {
    // This xsref table contains incorrect values for the 4th integral,
    // see https://github.com/scipy/xsref/issues/10.
    // So we instead use values that were verified using a Wolfram notebook.

    // xsref::test("it2j0y0", "d-d_d", |x| xsf::it2j0y0(x[0]));

    let xs = [-1.0, 10.0, -10.0, 1.0, 0.2];
    let expect_0 = [
        0.121_165_246_995_068_71,
        2.177_866_420_093_336,
        2.177_866_420_093_336,
        0.121_165_246_995_068_71,
        0.004_993_754_627_460_185,
    ];
    let expect_1 = [
        f64::NAN,
        -0.022_987_933_564_659_65,
        f64::NAN,
        0.395_272_901_699_293_2,
        -0.434_230_670_112_316_34,
    ];
    np_assert_allclose!(xs.map(|x| xsf::it2j0y0(x).0), expect_0);
    np_assert_allclose!(xs.map(|x| xsf::it2j0y0(x).1), expect_1);
}

#[test]
fn test_itj0y0_mid_x() {
    // the Neumann series for 0.5 <= x < 40, unlike in xsf (scipy/xsf#446) (mpmath)
    for (actual, expected, rtol) in [
        (xsf::it1j0y0(18.1).1, -0.000_306_924_251_694_864_96, 1e-11),
        (xsf::it1j0y0(19.63).1, -0.179_595_167_019_751_2, 1e-14),
        (xsf::it1j0y0(19.92).0, 1.044_817_073_262_731_5, 1e-15),
        (xsf::it2j0y0(20.0).0, 2.877_310_607_840_809_7, 1e-15),
        (xsf::it2j0y0(21.28).1, 9.723_256_856_027_431e-6, 1e-10),
    ] {
        assert!(
            (actual / expected - 1.0).abs() < rtol,
            "{actual:e} != {expected:e}"
        );
    }
}

xsref::tests! {
    test_it1i0k0: "it1i0k0", "d-d_d", |x| xsf::it1i0k0(x[0]);
    test_it2i0k0: "it2i0k0", "d-d_d", |x| xsf::it2i0k0(x[0]);
}

#[test]
fn test_iti0k0_large_x() {
    // unlike in xsf: quadrature and more asymptotic terms (mpmath)
    for (actual, expected, rtol) in [
        (xsf::it1i0k0(20.0).0, 44_758_598.913_855_75, 1e-12),
        (xsf::it1i0k0(30.0).0, 795_538_858_184.723_6, 1e-12),
        (xsf::it1i0k0(40.0).0, 1.508_971_508_271_920_2e16, 1e-14),
        (xsf::it2i0k0(40.0).0, 387_326_525_512_364.4, 1e-14),
        (
            xsf::it2i0k0(12.014_948_790_755_168).1,
            1.616_101_465_331_337_3e-7,
            1e-14,
        ),
        (
            xsf::it2i0k0(13.525_431_954_096_453).1,
            3.023_995_513_484_374e-8,
            1e-14,
        ),
        (xsf::it2i0k0(1.0625).1, 0.184_084_183_639_630_33, 1e-14),
        (xsf::it2i0k0(30.0).1, 6.780_268_362_546_424e-16, 1e-14),
        (xsf::it2i0k0(39.5).1, 3.399_304_068_768_309e-20, 1e-14),
        (xsf::it2i0k0(40.0).1, 2.024_204_614_539_913_3e-20, 1e-14),
        (xsf::it2i0k0(100.0).1, 4.588_521_432_774_676e-47, 1e-14),
    ] {
        assert!(
            (actual / expected - 1.0).abs() < rtol,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_itk0_mid_x() {
    // unlike in xsf: pi/2 minus the trapezoidal rule on int_x^inf K0 (mpmath)
    for (x, expected) in [
        (1.0, 1.242_509_848_623_778_2),
        (4.0, 1.560_649_570_605_174_2),
        (11.4, 1.570_792_372_273_874_1),
        (12.0, 1.570_794_206_646_576_8),
        (13.0, 1.570_795_574_915_093_5),
    ] {
        let actual = xsf::it1i0k0(x).1;
        assert!(
            (actual / expected - 1.0).abs() < 1e-15,
            "{actual:e} != {expected:e}"
        );
    }
}

xsref::tests! {
    test_besselpoly: "besselpoly", "d_d_d-d", |x| xsf::besselpoly(x[0], x[1], x[2]);
}

#[test]
fn test_besselpoly_large_a() {
    // mpmath: `r_besselpoly` in tools/accuracy/specs.py
    for (a, lambda, v, expected) in [
        (7.0, 0.0, 0.0, 0.080_048_416_665_887_81),
        (9.5, 6.5, 1.5, -0.004_626_673_452_494_62),
        (10.0, 2.0, -3.0, -0.008_181_711_690_050_824),
        (30.0, 1.0, 0.5, 0.001_769_720_528_262_000_3),
        (200.0, 0.5, 5.0, 0.000_376_502_079_621_611_3),
    ] {
        let actual = xsf::besselpoly(a, lambda, v);
        assert!(
            (actual / expected - 1.0).abs() < 1e-14,
            "{actual:e} != {expected:e}"
        );
    }
}

// Riccati-Bessel

/// Based on `scipy.special.tests.test_basic.TestRiccati.test_riccati_jn`
#[test]
fn test_riccati_j() {
    // N, x = 2, 0.2
    const N: usize = 2;
    let x: f64 = 0.2;
    // S = np.empty((N, N))
    let mut s = [[0.0; N], [0.0; N]];
    // for n in range(N):
    for n in 0..N {
        // j = special.spherical_jn(n, x)
        let j = xsf::sph_bessel_j(n.try_into().unwrap(), x);
        // jp = special.spherical_jn(n, x, derivative=True)
        let jp = xsf::sph_bessel_j_prime(n.try_into().unwrap(), x);
        // S[0,n] = x*j
        s[0][n] = x * j;
        // S[1,n] = x*jp + j
        s[1][n] = x * jp + j;
    }

    // assert_allclose(S, special.riccati_jn(n, x), atol=1.5e-8, rtol=0)
    let (jn, jnp) = xsf::riccati_j::<N>(x);
    xsref::np_assert_allclose!(&s[0], &jn, atol = 1.5e-8);
    xsref::np_assert_allclose!(&s[1], &jnp, atol = 1.5e-8);
}

/// Based on `scipy.special.tests.test_basic.TestRiccati.test_riccati_yn`
#[test]
fn test_riccati_y() {
    // N, x = 2, 0.2
    const N: usize = 2;
    let x: f64 = 0.2;
    // S = np.empty((N, N))
    let mut c = [[0.0; N], [0.0; N]];
    // for n in range(N):
    for n in 0..N {
        // j = special.spherical_jn(n, x)
        let j = xsf::sph_bessel_y(n.try_into().unwrap(), x);
        // jp = special.spherical_jn(n, x, derivative=True)
        let jp = xsf::sph_bessel_y_prime(n.try_into().unwrap(), x);
        // S[0,n] = x*j
        c[0][n] = x * j;
        // S[1,n] = x*jp + j
        c[1][n] = x * jp + j;
    }

    // assert_allclose(S, special.riccati_jn(n, x), atol=1.5e-8, rtol=0)
    let (yn, ynp) = xsf::riccati_y::<N>(x);
    xsref::np_assert_allclose!(&c[0], &yn, atol = 1.5e-8);
    xsref::np_assert_allclose!(&c[1], &ynp, atol = 1.5e-8);
}

#[test]
fn test_riccati_n1() {
    // computed with N = 2, since `rctj` and `rcty` write at least 2 terms
    let (j, jp) = xsf::riccati_j::<2>(0.2);
    assert_eq!(xsf::riccati_j::<1>(0.2), ([j[0]], [jp[0]]));
    let (y, yp) = xsf::riccati_y::<2>(0.2);
    assert_eq!(xsf::riccati_y::<1>(0.2), ([y[0]], [yp[0]]));
}

// Based on `scipy.special.tests.test_basic.TestLambda.test_lmbda`
#[test]
fn test_jahnke_emden_lambda() {
    // lam = special.lmbda(1,.1)
    let lam = xsf::jahnke_emden_lambda(1, 0.1);
    // lamr = (
    //     array([special.jn(0,.1), 2*special.jn(1,.1)/.1]),
    //     array([special.jvp(0,.1), -2*special.jv(1,.1)/.01 + 2*special.jvp(1,.1)/.1])
    // )
    let lamr = (
        [xsf::bessel_j(0.0, 0.1), 2e1 * xsf::bessel_j(1.0, 0.1)],
        [
            xsf::bessel_j_prime(0.0, 0.1, 1),
            -2e2 * xsf::bessel_j(1.0, 0.1) + 2e1 * xsf::bessel_j_prime(1.0, 0.1, 1),
        ],
    );
    // assert_allclose(lam, lamr, atol=1.5e-8, rtol=0)
    xsref::np_assert_allclose!(&lam.0, &lamr.0, atol = 1.5e-8);
    xsref::np_assert_allclose!(&lam.1, &lamr.1, atol = 1.5e-8);
}

#[test]
fn test_bessel_zeros_n0() {
    assert_eq!(xsf::bessel_zeros::<0>(0), [[]; 4]);
}

#[test]
fn test_bessel_zeros_large_v() {
    // no `int` overflow of `k * k` in `jynbh`, or in `msta1` and `msta2`
    let _ = xsf::bessel_zeros::<1>(40_000);
    assert!(xsf::bessel_zeros::<1>(u32::MAX)[0][0].is_nan());
}

#[test]
#[should_panic(expected = "v is too large")]
fn test_jahnke_emden_lambda_huge_v() {
    let _ = xsf::jahnke_emden_lambda(2_147_483_648.5, 1.0);
}

#[test]
fn test_jahnke_emden_lambda_nan() {
    // instead of out-of-bounds indexing in `lamn` and `lamv`
    for v in [2.0, 2.25] {
        for x in [f64::NAN, f64::INFINITY, -1e10] {
            let (vl, dl) = xsf::jahnke_emden_lambda(v, x);
            assert_eq!((vl.len(), dl.len()), (3, 3));
            assert!(vl.iter().chain(&dl).all(|y| y.is_nan()), "v = {v}, x = {x}");
        }
    }
}
