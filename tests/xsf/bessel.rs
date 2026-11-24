#![allow(clippy::needless_range_loop)]

use num_complex::{Complex64, c64};
use xsref::{np_assert_allclose, np_assert_equal};

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
fn test_bessel_miller_normalization() {
    // 2 < |z| < 21.8, where AMOS uses the Miller algorithm; Arb references
    let actual = [xsf::bessel_j(-0.3, c64(14.0, 3.0)), xsf::bessel_i(0.3, c64(18.2, 0.5))];
    let expected = [
        c64(0.732_258_715_423_858, -1.985_838_727_686_285_7),
        c64(6_659_440.736_169_231, 3_519_109.449_188_194_2),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-15);
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
    let expected =
        [8.829_186_306_726_423e-306, 4.182_877_977_322_057_6e306, 1.797_693_134_862_233_8e308];
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
    np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_bessel_i_large_v_cancellation() {
    // Arb: the terms of v eta in the uniform asymptotic expansion are large and cancel
    let actual = xsf::bessel_i(2277.875, 1606.5);
    np_assert_allclose!(&[actual], &[2.522_847_039_973_749_8e72], rtol = 1e-15);
}

#[test]
fn test_bessel_i_huge_v_laplace_limit() {
    // x / v near the Laplace limit, where eta = 0, but v eta still over- or underflows
    assert_eq!(xsf::bessel_i(1e100, 6.627_434_193_491_815e99), 0.0);
    assert_eq!(xsf::bessel_i(1e154, 6.627_434_193_491_816e153), f64::INFINITY);
    assert_eq!(xsf::bessel_i(1e155, 6.627_434_193_491_815e154), 0.0);
    // 2 pi v overflows
    assert_eq!(xsf::bessel_i(1e308, 6.627_434_193_491_816e307), f64::INFINITY);
}

#[test]
fn test_bessel_i_large_v_near_overflow() {
    // mpmath: besseli(v, x), where e^(v eta) or e^(-v eta) of the uniform expansion overflows
    let actual = [
        xsf::bessel_i(151.25, 725.75),
        xsf::bessel_i(-110.325, 0.128),
        xsf::bessel_i(60.0, 716.0),
        xsf::bessel_i(-60.5, 716.0),
        xsf::bessel_i(4000.0, 3_066.106_224),
        xsf::bessel_i(2081.0, -1_810.018_939_996_028_6),
        xsf::bessel_i(-1_770.246_227_066_590_8, 1_611.000_032_729_132),
    ];
    let expected = [
        3.426_801_707_970_078_3e306,
        9.211_668_616_148_808e307,
        1.087_560_322_184_507_3e308,
        1.042_772_683_333_504_7e308,
        1.797_692_243_068_206_4e308,
        -8.596_940_648_747_34e306,
        2.391_964_268_782_783e307,
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-15);
    // mpmath: 1.1e310
    assert_eq!(xsf::bessel_i(100.0, 725.0), f64::INFINITY);
}

#[test]
fn test_bessel_i_huge_x() {
    // I_v(x) ~ e^x / sqrt(2 pi x) overflows, also where 2 pi x or (x / v)^2 does
    let inf = f64::INFINITY;
    let x = 2f64.powi(1023);
    for v in [-49.5, -1e-30, 0.0, 0.5, 3.0, 50.0] {
        assert_eq!(xsf::bessel_i(v, x), inf);
    }
    for (v, x) in [(100.0, 1e200), (-100.5, 1e200), (1e5, 1e300)] {
        assert_eq!(xsf::bessel_i(v, x), inf);
    }
    // odd orders
    assert_eq!(xsf::bessel_i(3.0, -x), -inf);
    assert_eq!(xsf::bessel_i(51.0, -1e200), -inf);
}

#[test]
fn test_bessel_amos_underflow() {
    // mpmath: besselj(7.3156, 1.19e-41) = 2.6872935555094642e-306,
    // besselj(-43, 4.4e-6) = -8.770585860904835e-297,
    // besselj(5, -2.9e-60) = -5.341445052083332e-302,
    // besseli(46.408, 9.696e-6) exp(-9.696e-6) = 8.829100699351018e-306, and
    // besseli(155.25, 1.25) exp(-1.25) = 3.4698644197170413e-307
    let actual = [
        xsf::bessel_j(7.3156, 1.19e-41),
        xsf::bessel_je(7.3156, 1.19e-41),
        xsf::bessel_j(-43.0, 4.4e-6),
        xsf::bessel_j(5.0, -2.9e-60),
        xsf::bessel_ie(46.408, 9.696e-6),
        xsf::bessel_ie(155.25, 1.25),
    ];
    let expected = [
        2.687_293_555_509_464_2e-306,
        2.687_293_555_509_464_2e-306,
        -8.770_585_860_904_835e-297,
        -5.341_445_052_083_332e-302,
        8.829_100_699_351_018e-306,
        3.469_864_419_717_041_3e-307,
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
fn test_bessel_jy_negative_order_overflow() {
    // mpmath: bessely(10.25, 1e-100) = -2.5e1033, besselj(-10.25, 1e-100) = 1.8e1033,
    // bessely(-200.5, 1) = 3.9e-437, besselj(-1611.5, 200) = -3.9e1246, and
    // besselj(-10.001, 1e-30) = 3.993385638016181e305
    let inf = f64::INFINITY;
    assert_eq!(xsf::bessel_ye(10.25, 1e-100), -inf);
    assert_eq!(xsf::bessel_je(-10.25, 1e-100), inf);
    assert_eq!(xsf::bessel_y(-200.5, 1.0), 0.0);
    assert_eq!(xsf::bessel_j(-1611.5, 200.0), -inf);
    let j = xsf::bessel_j(-10.001, 1e-30);
    np_assert_allclose!(&[j], &[3.993_385_638_016_181e305], rtol = 1e-13);
}

#[test]
fn test_bessel_j_power_series_overflow_sign() {
    // unlike in xsf, the sign is kept if the power series of cephes `jv` overflows
    let inf = f64::INFINITY;
    for (v, x, expected) in [
        (-123.573_889_191_806_12, 0.028_605_192_190_944_706, -inf),
        (-124.683_117_666_899_78, 0.032_989_216_721_374_115, inf),
        (-111.294_155_042_714_54, 0.060_942_186_932_616_55, -inf),
    ] {
        assert_eq!(xsf::bessel_j(v, x), expected);
    }
}

#[test]
fn test_bessel_ie_negative_order_large_imag() {
    // mpmath (dps = 50): besseli(-7.25, 1e9j), besseli(-0.5, 1e9j), and besseli(-7.25, 1e8j)
    let actual = [
        xsf::bessel_ie(-7.25, c64(0.0, 1e9)),
        xsf::bessel_ie(-0.5, c64(0.0, 1e9)),
        xsf::bessel_ie(-7.25, c64(0.0, 1e8)),
    ];
    let expected = [
        c64(1.773_234_206_304_038_3e-6, 4.280_966_070_123_1e-6),
        c64(1.494_894_726_886_609e-5, -1.494_894_726_886_609e-5),
        c64(3.052_711_541_179_276e-5, 7.369_897_604_727_882e-5),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-12);
}

#[test]
fn test_bessel_ye_large_order_large_imag() {
    // Arb: |z| < v, where e^-2|Im z| underflows (was half the value)
    let z = [
        c64(-56.628_806_048_134_02, 642.064_090_486_218_4),
        c64(0.874_908_836_744_878, 654.491_289_262_904_2),
        c64(-76.629_816_771_216_21, 610.825_715_575_559_7),
    ];
    let expected = [
        c64(2.074_005_539_824_065e-266, -2.802_308_574_285_694_5e-266),
        c64(1.171_146_111_371_311_2e-281, -4.330_334_894_877_472e-280),
        c64(4.744_524_902_874_874e-229, 1.499_209_524_082_626e-228),
    ];
    np_assert_allclose!(&z.map(|z| xsf::bessel_ye(1000.0, z)), &expected, rtol = 3e-13);
    let conj = |z: [Complex64; 3]| z.map(|z| z.conj());
    np_assert_allclose!(&conj(z).map(|z| xsf::bessel_ye(1000.0, z)), &conj(expected), rtol = 3e-13);
}

#[test]
fn test_bessel_k_huge_order_overflow() {
    // K_v(x) and e^x K_v(x) overflow (DLMF 10.41.2, 10.41.4); for the last, ln(e^x K_v(x)) = 3.5e8
    let inf = f64::INFINITY;
    let v30 = 2.0_f64.powi(30);
    for (v, x) in [(v30, 1.0), (-1e12, 1.0), (1e300, 1e200)] {
        assert_eq!(xsf::bessel_k(v, x), inf);
        assert_eq!(xsf::bessel_ke(v, x), inf);
        assert_eq!(xsf::bessel_k(v, c64(x, 0.0)), c64(inf, 0.0));
    }
    assert_eq!(xsf::bessel_ke(v30, 1.5 * v30), inf);
}

#[test]
fn test_bessel_jyi_near_overflow() {
    // mpmath: bessely(v, x), besselj(v, x), and besseli(v, x) * exp(-x)
    let actual = [
        xsf::bessel_y(79.762_372_861_572_26, 8.365_870_974_604_77e-3),
        xsf::bessel_ye(-61.5625, 4.6e-4),
        xsf::bessel_je(-4.02, 6e-77),
        xsf::bessel_ie(-3.82, 6e-81),
        xsf::bessel_ie(-2.000_000_000_000_000_4, 4.214_684_851_089_752e-162),
    ];
    let expected = [
        -5.240_453_974_629_458_6e305,
        -4.986_124_250_198_734e305,
        5.149_951_225_514_283e306,
        -3.245_335_795_067_041e307,
        1e308,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-12);
    // mpmath: 1.9e308 and 1.99e-37 (exp(-2x) underflows, so no retry)
    for (v, x) in [(-2.000_000_000_000_000_4, 3.057_653_761_164_066_6e-162), (-840.25, 400.0)] {
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
    // also for x < 0, as (-1)^v J_v(-x)
    for x in [2_149_631_131.648, -2_149_631_131.648] {
        let j = xsf::bessel_j(2_147_483_648.0, x);
        np_assert_allclose!(&[j], &[-1.632_266_889_229_664_4e-5], rtol = 1e-6);
    }
}

#[test]
fn test_bessel_jy_huge_order_limits() {
    // |J_v(x)| <= e^-t and Y_v(x) ~ -e^t / sqrt(pi v w / 2) for x < v, with t >= v w^3 / 3 and
    // w^2 = 1 - (x / v)^2 (DLMF 10.14.5, 10.19.3); Arb: besselj(1e100, 5) = 1.8e-9.9e101
    let inf = f64::INFINITY;
    // even integers
    for (v, x) in [
        (2_251_799_813_685_250.0, 5.0),
        (1e18, 9.9e17),
        (3.162_277_660_168_379e50, 1e8),
        (1e100, 5.0),
    ] {
        for (v, x) in [(v, x), (v, -x), (-v, x), (-v, -x)] {
            assert_eq!(xsf::bessel_j(v, x), 0.0);
            assert_eq!(xsf::bessel_je(v, x), 0.0);
            assert_eq!(xsf::bessel_j(v, c64(x, 0.0)), c64(0.0, 0.0));
        }
        for v in [v, -v] {
            assert_eq!(xsf::bessel_y(v, x), -inf);
            assert_eq!(xsf::bessel_ye(v, x), -inf);
            assert_eq!(xsf::bessel_y(v, c64(x, 0.0)), c64(-inf, 0.0));
        }
    }
    // the odd integer 2^52 + 1: (-1)^v for -v and -x, also as the sign of the 0
    let odd = 4_503_599_627_370_497.0;
    assert_eq!(xsf::bessel_y(-odd, 5.0), inf);
    for (j, negative) in [
        (xsf::bessel_j(odd, 5.0), false),
        (xsf::bessel_j(odd, -5.0), true),
        (xsf::bessel_j(-odd, 5.0), true),
        (xsf::bessel_j(-odd, -5.0), false),
        (xsf::bessel_j(odd, -0.0), true),
        (xsf::bessel_j(-odd, -0.0), false),
        (xsf::bessel_je(odd, -0.0), true),
        (xsf::bessel_je(-odd, -0.0), false),
        (xsf::bessel_j(odd, c64(-5.0, 0.0)).re, true),
        (xsf::bessel_je(odd, c64(-5.0, 0.0)).re, true),
        (xsf::bessel_j(-odd, c64(-5.0, 0.0)).re, false),
    ] {
        assert!(j == 0.0 && j.is_sign_negative() == negative, "{j}");
    }
    // the half-integer 2^51 + 1/2: J_-v = -Y_v and Y_-v = J_v (DLMF 10.2.3)
    let half = 2_251_799_813_685_248.5;
    assert_eq!(xsf::bessel_j(-half, 5.0), inf);
    assert_eq!(xsf::bessel_y(-half, 5.0), 0.0);
    // Arb: besselj(3e9, -5) = 2.7e-25934660298
    assert_eq!(xsf::bessel_j(3e9, -5.0), 0.0);
}

#[test]
fn test_bessel_y_large_integer_order_huge_x() {
    // Arb: bessely(n, x) for x above the AMOS limit 2^51
    let actual = [
        xsf::bessel_y(2_147_483_647.0, 1e16),
        xsf::bessel_y(-2_147_483_648.0, 1e16),
        xsf::bessel_y(2e8, 1e20),
        xsf::bessel_y(1_000_000_007.0, 2_251_799_813_685_250.0),
        xsf::bessel_y(3.0, 1e300),
    ];
    let expected = [
        7.247_404_648_227_662e-9,
        -3.337_228_298_107_891_7e-9,
        -7.950_547_863_231_89e-11,
        -5.381_323_928_723_016_5e-9,
        -7.860_673_062_724_093e-151,
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-13);
}

#[test]
fn test_bessel_jy_large_x() {
    // Arb: besselj(v, x) and bessely(-v, x) for |v| <= 85.9 < x, where AMOS besj lost up to 1.2e-11
    let (v, x) = (49.75, 676.5);
    let (v2, x2) = (80.375, 3112.25);
    let actual = [
        xsf::bessel_j(v, x),
        xsf::bessel_y(-v, x),
        xsf::bessel_j(85.375, 3305.5),
        xsf::bessel_j(-v2, x2),
        xsf::bessel_y(-v2, x2),
    ];
    let expected = [
        -0.024_514_334_729_217_442,
        0.030_423_395_818_865_59,
        0.003_737_638_601_131_790_8,
        -0.013_938_274_494_678_701,
        0.003_216_465_737_266_603_6,
    ];
    np_assert_allclose!(&actual, &expected, rtol = 2e-14);
}

#[test]
fn test_bessel_j_oscillatory() {
    // Arb: besselj(v, x) for |x| > |v| >= 10, where AMOS lost ~eps |x| through its phase (up to
    // 5e-7 near a zero), and J_v(x) = 0 or the wrong sign for x >= 1e17
    let vx = [
        (2611.25, 2_502_847.0),
        (-2_652.249_251_597_461_4, 598_667.409_649_828),
        (-2213.75, 2_202_672.5),
        (861.625, 1036.375),
        (101.0, -1000.0),
        (-101.0, -1000.0),
        (-102.0, 1000.0),
        (10.0, 80.0),
        (1e7, 1e9),
        (100.0, 1e20),
        (100.5, 1e308),
    ];
    let expected = [
        3.642_473_563_404_470_3e-6,
        -8.122_052_923_768_204e-7,
        -0.000_146_941_885_726_263_02,
        -0.003_007_408_152_349_567_5,
        0.021_152_705_375_840_766,
        -0.021_152_705_375_840_766,
        -0.015_948_981_493_722_39,
        0.024_043_850_978_184_764,
        4.785_055_993_128_89e-6,
        6.698_009_040_703_428e-12,
        3.617_580_596_934_688_3e-155,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
}

#[test]
fn test_bessel_j_turning_point() {
    // Arb: besselj(v, x) for v >= 10 below and around x = v, where AMOS loses ~eps of its exponent
    // or phase
    let vx = [
        (1000.0, 950.0),
        (300.0, 270.0),
        (60.5, 3.25),
        (50.0, 20.0),
        (1000.0, 1005.0),
        (5000.0, 5030.0),
        (3000.5, 3010.0),
        (301.0, -270.0),
    ];
    let expected = [
        4.639_221_366_465_356e-7,
        2.933_415_684_728_316_7e-6,
        8.433_834_271_887_743e-71,
        4.451_039_284_700_681e-16,
        0.062_916_645_435_363_97,
        0.006_749_289_286_719_537,
        0.045_959_772_400_047_88,
        -1.816_491_558_189_553_6e-6,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 1e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
    // subnormal, rounded
    np_assert_allclose!([xsf::bessel_j(903.75, 303.25)], [5.79e-322], rtol = 5e-3);
    // DLMF 10.19.8: Ai(0) (2 / v)^(1/3)
    let expected = [4.473_073_183_964_723e-101];
    np_assert_allclose!([xsf::bessel_j(1e300, 1e300)], expected, rtol = 1e-15);
    // Arb, by DLMF 10.20.4 up to A_1 and B_1: an ulp below x = v, where eta ~ v w^3 / 3 cancels
    // in its closed form
    let expected = [1.891_197_090_318_688e-8];
    np_assert_allclose!([xsf::bessel_j(1e22, 9.999_999_999_999_998e21)], expected, rtol = 1e-15);
}

#[test]
fn test_bessel_y_oscillatory() {
    // Arb: bessely(v, x) for x > |v| >= 10, where AMOS lost ~eps x through its phase, and
    // Y_v(x) = 0 for v = 1e7 or x = 1e308
    let vx = [
        (2611.25, 2_502_847.0),
        (-2652.25, 598_667.5),
        (861.625, 1036.375),
        (-101.0, 1000.0),
        (1e7, 1e9),
        (-100.5, 1e308),
    ];
    let expected = [
        -0.000_504_326_398_870_437_1,
        -0.001_026_972_783_890_225_6,
        -0.033_111_163_526_816_485,
        0.013_872_810_696_637_188,
        -2.477_407_606_524_094_2e-5,
        3.617_580_596_934_688_3e-155,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_ye(v, x)), actual);
}

#[test]
fn test_bessel_y_turning_point() {
    // Arb: bessely(v, x) for |v| >= 10 below and around x = |v|, where AMOS loses ~eps of its
    // exponent or phase
    let vx = [
        (1000.0, 950.0),
        (300.0, 270.0),
        (60.5, 3.25),
        (50.0, 20.0),
        (1000.0, 1010.0),
        (3000.5, 3010.0),
        (-100.5, 95.0),
        (-2000.5, 2020.0),
    ];
    let expected = [
        -2_198.697_783_894_989,
        -830.477_035_670_576_9,
        -6.247_372_440_247_314e67,
        -15_606_426_801_663.736,
        0.006_232_063_932_500_365,
        -0.017_754_079_052_189_897,
        0.019_243_195_235_249_747,
        0.026_065_963_491_770_55,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, rtol = 1e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_ye(v, x)), actual);
    // DLMF 10.19.8: -Bi(0) (2 / v)^(1/3)
    let expected = [-7.747_590_020_600_788e-101];
    np_assert_allclose!([xsf::bessel_y(1e300, 1e300)], expected, rtol = 1e-15);
    // Arb, by DLMF 10.20.5 up to A_1 and B_1: an ulp below x = v
    let expected = [-3.918_787_372_858_491e-8];
    np_assert_allclose!([xsf::bessel_y(1e22, 9.999_999_999_999_998e21)], expected, rtol = 1e-15);
}

#[test]
fn test_bessel_j_turning_point_negative_order() {
    // Arb: besselj(v, x) for v <= -10 below and around x = |v|
    let vx = [(-100.5, 95.0), (-50.5, 20.0), (-2000.5, 2020.0), (-301.0, 270.0)];
    let expected = [
        0.523_152_470_400_604_5,
        34_050_157_434_984.086,
        -0.039_467_973_608_718_98,
        -1.816_491_558_189_553_6e-6,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 1e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
}

#[test]
fn test_bessel_jy_below_turning_point() {
    // Arb: besselj(v, x) and bessely(v, x) for x / |v| in [0.80, 0.87] with eta in [4, 6), where
    // AMOS loses up to ~430 eps
    let vx = [
        (77.264_841_492_886_82, 63.616_263_749_279_72),
        (83.689_103_548_745_25, 71.430_835_724_982_23),
        (60.0, 48.0),
        (110.0, 95.0),
        (-70.5, 58.5),
    ];
    let expected = [
        0.000_163_471_623_958_521_48,
        0.000_516_894_790_852_201_5,
        0.000_245_751_851_522_892_3,
        0.000_200_535_695_968_452_25,
        20.631_660_341_613_38,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
    let expected = [
        -44.498_833_410_646_4,
        -14.169_708_874_182_264,
        -36.063_813_253_052_786,
        -28.692_287_874_628_192,
        0.000_393_271_029_089_817_3,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_ye(v, x)), actual);
}

#[test]
fn test_bessel_jy_below_turning_point_recurrence() {
    // Arb: besselj(v, x) and bessely(v, x) for x / |v| in [0.82, 0.87] with eta < 4, where AMOS
    // loses up to ~280 eps
    let vx = [(83.375, 72.5), (70.0, 60.0), (50.5, 43.0), (40.0, 33.0), (-62.5, 54.0)];
    let expected = [
        0.001_187_826_402_904_425_9,
        0.001_436_386_052_336_167_7,
        0.004_062_204_729_530_515,
        0.004_058_987_686_471_918,
        3.430_481_146_054_990_7,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 3e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
    let expected = [
        -6.542_065_712_945_027,
        -6.180_085_145_468_044,
        -2.988_764_285_772_928,
        -3.501_917_079_225_303,
        0.002_973_733_550_822_862_3,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, rtol = 3e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_ye(v, x)), actual);
}

#[test]
fn test_bessel_jy_below_turning_point_recurrence_edges() {
    // Arb: besselj(v, x) and bessely(v, x); the lowest order of the recurrence is 10 (the least of
    // the expansion), 9 (so AMOS), and 17, and at 70 it is one lower, as order 69 rounds to just
    // outside the power series
    let vx = [(13.0, 10.0), (13.0, 9.5), (20.0, 16.0), (70.0, 60.148_847_869_265_126)];
    let expected = [
        0.028_972_083_926_776_766,
        0.018_156_064_603_798_993,
        0.017_328_746_227_591_996,
        0.001_575_026_341_020_990_7,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 3e-15);
    let expected = [
        -1.363_454_319_803_315,
        -2.018_456_307_835_525_5,
        -1.567_173_990_858_811_6,
        -5.676_985_765_349_779,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, rtol = 3e-15);
}

#[test]
fn test_bessel_jy_above_turning_point() {
    // Arb: besselj(v, x) and bessely(v, x) for x / |v| in [1.11, 1.36] with eta in [3, 6), where
    // AMOS loses up to ~200 eps of their amplitude
    let vx = [
        (81.046_006_211_873_38, 95.862_273_724_479_32),
        (100.0, 115.0),
        (150.0, 168.5),
        (-70.5, 82.0),
        (60.0, 70.0),
    ];
    let expected = [
        0.002_786_135_086_496_126_6,
        -0.038_196_043_224_737_4,
        0.027_066_073_979_626_733,
        0.017_885_916_555_457_456,
        -0.124_230_136_973_084_74,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
    let expected = [
        -0.111_355_875_330_217_51,
        -0.098_607_469_724_625_39,
        -0.086_859_880_177_404_84,
        -0.121_752_801_011_926,
        0.046_233_208_026_780_416,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_ye(v, x)), actual);
}

#[test]
fn test_bessel_jy_above_turning_point_recurrence() {
    // Arb: besselj(v, x) and bessely(v, x) for x / |v| in [1.11, 1.36] with eta < 3, where AMOS
    // loses up to ~80 eps of their amplitude
    let vx = [(20.0, 24.0), (30.5, 35.0), (33.0, 38.0), (58.0, 65.0), (-41.5, 47.5)];
    let expected = [
        0.161_912_651_664_495_28,
        0.145_025_859_216_876_58,
        0.113_311_343_523_371_29,
        0.030_263_581_939_060_775,
        0.154_045_968_240_138_16,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, rtol = 1e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_je(v, x)), actual);
    let expected = [
        0.144_201_886_060_214_56,
        0.123_597_316_078_740_12,
        0.142_700_612_314_663_6,
        0.143_263_581_238_127_85,
        -0.058_755_646_317_514_87,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    // 11 ulp in the last one with the Windows libm (2 ulp on Linux)
    np_assert_allclose!(actual, expected, rtol = 2e-15);
    np_assert_equal!(vx.map(|(v, x)| xsf::bessel_ye(v, x)), actual);
}

#[test]
fn test_bessel_jy_above_turning_point_recurrence_edges() {
    // Arb: besselj(v, x) and bessely(v, x); the lowest order of the recurrence is 10 (the least of
    // the expansion), 9 (so AMOS), and both reflection terms of a negative order are nonzero
    let vx = [(15.0, 18.0), (16.0, 18.4), (16.0, 17.92), (-30.3, 36.0)];
    let expected = [
        0.235_592_357_774_215_22,
        0.257_888_878_317_665_5,
        0.260_626_196_074_620_45,
        -0.117_522_568_584_565_39,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_j(v, x));
    np_assert_allclose!(actual, expected, atol = 1e-15);
    let expected = [
        0.080_480_924_919_380_9,
        -0.006_919_837_599_588_943_5,
        -0.071_686_960_171_829_57,
        0.136_178_378_274_856_35,
    ];
    let actual = vx.map(|(v, x)| xsf::bessel_y(v, x));
    np_assert_allclose!(actual, expected, atol = 1e-15);
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
        xsf::bessel_k(27.375, 1.44e-10),
        xsf::bessel_ke(1.280_346_742_634_130_8, 2.382_283_514_915_581_2e-239),
        xsf::bessel_ke(46.4, 8.625e-6),
    ];
    let expected = [
        3.963_327_297_606_011e307,
        3.125_975_622_874_344_7e304,
        3.623_615_504_892_954e305,
        2.450_893_202_543_499_5e305,
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
fn test_bessel_k_at_amos_overflow() {
    for (v, x, expected) in [
        (1.5, 1e-205, [3.963_327_297_606_011e307, 3.963_327_297_606_011e307]),
        (-70.0, 0.002, [8.556_122_497_405_277e307, 8.573_251_866_058_952e307]),
    ] {
        let actual = [xsf::bessel_k(v, x), xsf::bessel_ke(v, x)];
        xsref::np_assert_allclose!(actual, expected, rtol = 5e-14);
    }
    assert_eq!(xsf::bessel_k(1.5, 1e-206), f64::INFINITY);
    assert_eq!(xsf::bessel_ke(1.5, 1e-206), f64::INFINITY);
}

#[test]
fn test_bessel_k_near_underflow() {
    // mpmath: besselk(v, x)
    let actual = [
        xsf::bessel_k(0.0, 700.6),
        xsf::bessel_k(12.49, 700.6),
        xsf::bessel_k(-20.0, 700.6),
        xsf::bessel_k(0.0, c64(700.6, 1.0)).re,
        xsf::bessel_k(0.0, c64(700.6, 1.0)).im,
    ];
    let expected = [
        2.561_730_385_750_867_3e-306,
        2.863_182_898_540_892_5e-306,
        3.407_332_893_410_498_4e-306,
        1.382_569_916_659_872e-306,
        -2.156_607_594_040_158_8e-306,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
    // subnormal, rounded
    let actual = [xsf::bessel_k(1.0, 708.0), xsf::bessel_k(0.0, 720.0)];
    let expected = [1.558_762_642_090_316e-309, 9.490_549_83e-315];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-9);
    // mpmath: 9.6e-325, 2.4e-325, and below e^-800
    for (v, x) in [(0.0, 743.0), (30.25, 745.0), (100.0, 800.0), (1e9, 1.1e9)] {
        assert_eq!(xsf::bessel_k(v, x), 0.0);
    }
}

#[test]
fn test_bessel_k_normal_near_underflow() {
    // Regression cases from xsf where AMOS returned zero although K is normal, with mpmath references.
    for (v, z, expected) in [
        (0.0, c64(699.0, 0.0), c64(1.270_284_188_032_741_8e-305, 0.0)),
        (50.0, c64(700.0, 0.0), c64(2.779_335_877_012_058_7e-305, 0.0)),
        (0.5, c64(702.0, -50.0), c64(6.138_086_236_570_686e-307, -1.436_738_677_334_688_3e-307)),
        (0.0, c64(700.0, 1e5), c64(-2.869_453_464_726_56e-307, 2.652_560_945_968_814_4e-307)),
    ] {
        np_assert_allclose!([xsf::bessel_k(v, z)], [expected], rtol = 5e-14);
        if z.im == 0.0 {
            np_assert_allclose!([xsf::bessel_k(v, z.re)], [expected.re], rtol = 5e-14);
        }
    }
}

#[test]
fn test_bessel_k_large_order_near_underflow() {
    // Arb: besselk(v, x), where e^x K_v(x) overflows
    let actual = [
        xsf::bessel_k(2000.0, 1749.66),
        xsf::bessel_k(-1e4, 7022.24),
        xsf::bessel_k(1e5, 66659.8),
        xsf::bessel_k(3000.0, c64(2400.63, 0.0)).re,
    ];
    let expected = [
        4.932_831_061_297_247e-306,
        6.630_171_739_413_170_6e-307,
        1.383_450_601_290_715_8e-305,
        1.796_796_534_493_318e-306,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-11);
    assert_eq!(xsf::bessel_k(3000.0, c64(2400.63, 0.0)).im, 0.0);
    // subnormal, rounded
    let actual = [xsf::bessel_k(5000.0, 3731.38)];
    xsref::np_assert_allclose!(&actual, &[9.233_667e-318], rtol = 1e-9);
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
fn test_bessel_amos_i_leading_exponential_overflow() {
    // the leading exponentials exceed `f64::MAX`, unlike the values (Arb)
    let actual = [
        xsf::bessel_i(54.099_467_124_448_5, c64(711.854_673_618_246_2, 0.0)),
        xsf::bessel_i(1_490.491_370_199_187_7, c64(-1_431.443_889_835_916_7, 0.0)),
        xsf::bessel_i(-3_924.611_809_854_627_6, c64(-3_014.674_352_211_202, 0.0)),
    ];
    let expected = [
        c64(2.731_246_300_865_973e306, 0.0),
        c64(4.494_375_122_683_454e304, 1.657_342_432_873_988e306),
        c64(-2.999_368_964_669_362e306, -8.184_750_647_701_384e306),
    ];
    // AMOS loses some accuracy at these large orders, also below the overflow threshold
    np_assert_allclose!(&actual, &expected, rtol = 1e-12);

    // Arb: besselk(100.25, -712.2989872200264 + 1500j), whose pi I term exceeds `f64::MAX`
    // by its leading exponential only
    let k = xsf::bessel_k(100.25, c64(-712.298_987_220_026_4, 1500.0));
    let expected = c64(-8.786_131_234_745_587e306, -1.649_448_687_188_380_2e307);
    np_assert_allclose!(&[k], &[expected], rtol = 1e-12);
}

#[test]
fn test_bessel_amos_i_overflow_by_series() {
    // the real part overflows, unlike the estimates without the series (Arb, halved; a half-ulp
    // change of Im z changes the imaginary part by 1.67e-9)
    let i = xsf::bessel_i(100.25, c64(715.712_422_528_251_9, 1_499.513_771_509_891_7));
    assert!(i.re == f64::INFINITY, "{i}");
    assert!(near_part(i.im, -6.116_621_289_670_503e303, 1.67e-9), "{i}");
    // the exponent has a rounding error of about 1 (DLMF 10.41.3: -2.28e308 (1 + i))
    let i = xsf::bessel_i(2_060_697_274_498_541.2, c64(-1_365_713_557_945_105.0, 0.0));
    assert!(i.re.is_nan() && i.im.is_nan(), "{i}");
}

#[test]
fn test_bessel_amos_i_overflow_by_part() {
    // per part, within the change from a half-ulp change of Re z or Im z (Arb, halved)
    let i = xsf::bessel_i(385.75, c64(716.044_300_663_178_2, 32_760.175_578_041_82));
    assert!(near_part(i.re, 4.098_878_142_449_450_5e307, 4.0e-12), "{i}");
    assert!(i.im == f64::INFINITY, "{i}");
    // the reflection to v < 0 keeps the finite part (the K term is about 1e-312)
    let i = xsf::bessel_i(-100.25, c64(715.712_422_528_251_9, 1_499.513_771_509_891_7));
    assert!(i.re == f64::INFINITY, "{i}");
    assert!(near_part(i.im, -6.116_621_289_670_503e303, 1.67e-9), "{i}");
    // above `f64::MAX` by 6.8e-14, but a half-ulp change of Im z changes the real part by 1.75e-11
    let i = xsf::bessel_i(385.75, c64(718.220_743_426_067_4, 32_760.5));
    assert!(near_part(i.re, 8.988_465_674_312_192e307, 1.75e-11), "{i}");
    assert!(i.im == f64::INFINITY, "{i}");
}

#[test]
fn test_bessel_amos_i_overflow_by_part_lost() {
    // as in xsf, where a part of `e^-|Re z| I_v(z)` underflows to 0 or a subnormal, or cancels in
    // the rotation for `Re z < 0`, so that it lost its digits
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    for (v, z, re, im) in [
        (0.0, c64(720.0, 5e-324), inf, nan),
        (0.0, c64(720.0, 1e-320), inf, inf),
        (0.0, c64(1500.0, 1e-310), inf, inf),
        (0.0, c64(1500.0, 3.16e-322), inf, inf),
        (0.5, c64(-1500.0, 1e-322), inf, inf),
        (1500.25, c64(-1_440.844_949_792_818_9, 0.0), nan, nan),
        (1500.5, c64(-1_441.2, -1e-25), nan, nan),
        (171.5, c64(-734.340_839_421_059_9, -1e-25), -inf, inf),
        (56.670_551_070_453_26, c64(-720.438_283_905_823_6, 0.534_519_658_542_622_3), nan, inf),
        (56.670_551_070_453_31, c64(-720.438_283_905_823_5, 0.534_519_658_542_622_4), -inf, inf),
    ] {
        let i = xsf::bessel_i(v, z);
        // NaN, or an inf of the same sign
        let same = |x: f64, y: f64| {
            if y.is_nan() {
                x.is_nan()
            } else {
                x.is_infinite() && x.is_sign_negative() == y.is_sign_negative()
            }
        };
        assert!(same(i.re, re) && same(i.im, im), "I_{v}({z}) = {i}");
    }
}

/// Whether `x` is within `rtol` of the exact `2 h` (halved to stay finite), or overflows where
/// `rtol` can reach `f64::MAX`.
fn near_part(x: f64, h: f64, rtol: f64) -> bool {
    if x.is_infinite() {
        x.is_sign_negative() == h.is_sign_negative() && h.abs() * (1.0 + rtol) >= f64::MAX / 2.0
    } else {
        (x / 2.0 - h).abs() <= rtol * h.abs()
    }
}

/// Whether each part of `y` is within AMOS's error of that of the exact value `2 half` (halved to
/// stay finite), or overflows where that error can reach `f64::MAX`. The error is at most
/// `4 max(v, |z|) eps |2 half|`, an empirical margin (Arb: at most 1.96 times
/// `max(v, |z|) eps` in the larger part, at 684 values of orders up to 19520.5).
fn is_near(y: Complex64, half: Complex64, v: f64, z: Complex64) -> bool {
    let tol = 4.0 * v.max(z.norm()) * f64::EPSILON * half.norm();
    let near = |x: f64, h: f64| {
        if x.is_infinite() {
            x.is_sign_negative() == h.is_sign_negative() && h.abs() + tol >= f64::MAX / 2.0
        } else {
            (x / 2.0 - h).abs() <= tol
        }
    };
    near(y.re, half.re) && near(y.im, half.im)
}

#[test]
fn test_bessel_amos_i_overflow_by_rounding() {
    // within the error of `f64::MAX` (Arb, halved: above it by 8.9e-12, 2.4e-12, 3.4e-12, and
    // 1.8e-12)
    let cases = [
        (32_767.75, c64(22_114.649_592_571_113, 0.0), c64(8.988_465_674_391_864e307, 0.0)),
        (32_767.75, c64(22_114.649_592_571_11, 0.0), c64(8.988_465_674_333_41e307, 0.0)),
        (
            32_766.999_999_9,
            c64(-22_114.152_585_700_04, 0.0),
            c64(-8.988_465_674_341_9e307, 2.823_831_274_836_541e301),
        ),
        (
            385.75,
            c64(716.044_300_663_177_6, 32_760.175_578_041_82),
            c64(4.098_878_142_447_121e307, 8.988_465_674_327_425e307),
        ),
    ];
    for (v, z, half) in cases {
        let i = xsf::bessel_i(v, z);
        assert!(is_near(i, half, v, z), "I({v}, {z}) = {i}");
    }
    // further below `f64::MAX`, the parts are finite (Arb)
    let actual = [
        xsf::bessel_i(32_767.75, c64(22_114.649_592_571_088, 0.0)),
        xsf::bessel_i(32_766.999_999_9, c64(-22_114.152_585_700_012, 0.0)),
    ];
    let expected = [
        c64(1.797_693_134_796_537_6e308, 0.0),
        c64(-1.797_693_134_774_854e308, 5.647_662_549_379_259e301),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-10);
}

#[test]
fn test_bessel_amos_overflow_by_rotation() {
    // one part overflows once rotated to z (Arb)
    let i = xsf::bessel_i(4.25, c64(-714.1, -49.5));
    let j = xsf::bessel_j(4.25, c64(-49.5, 714.1));
    assert_eq!((i.im, j.re), (f64::NEG_INFINITY, f64::INFINITY));
    let i2 = xsf::bessel_i(69.25, c64(-717.43, -0.75));
    let j2 = xsf::bessel_j(69.25, c64(-1.17, 717.39));
    assert_eq!((i2.re, j2.re), (f64::NEG_INFINITY, f64::NEG_INFINITY));
    // within the error of `f64::MAX` once rotated to z (Arb, halved: above it by 4.9e-12, with
    // parts below 0.913 f64::MAX before)
    let (v, z) = (32_767.75, c64(-22_114.530_640_552_63, 95.5));
    let i3 = xsf::bessel_i(v, z);
    let half = c64(-2.614_861_393_023_947e307, -8.988_465_674_355_52e307);
    assert!(is_near(i3, half, v, z), "{i3}");
    np_assert_allclose!(
        &[i.re, j.im, i2.im, j2.im],
        &[
            -2.744_727_787_941_251e306,
            -7.849_902_985_860_507e307,
            6.434_107_269_133_3e306,
            6.608_834_223_121_121e305,
        ],
        rtol = 1e-12
    );
}

#[test]
fn test_bessel_amos_k_overflow_by_continuation() {
    // within the error of `f64::MAX` by the continuation to Re z < 0 (Arb, halved: below it by
    // 2.2e-12), where the real part has the error of AMOS's phase, as further below `f64::MAX`
    let (v, z) = (673.851_921_288_337_3, c64(-716.234_998_021_642, 8_699.109_351_368_823));
    let k = xsf::bessel_k(v, z);
    let half = c64(8.863_697_362_548_706e306, 8.988_465_674_291_945e307);
    assert!(is_near(k, half, v, z), "{k}");
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
fn test_bessel_amos_huge_z() {
    // up to `|z| = 0.5 / eps`, but the phase is only accurate to about `eps |z|` (Arb)
    let y = xsf::bessel_y(86.5, c64(1e10, 10.0));
    let desired = c64(0.076_723_556_342_248_09, 0.042_838_532_018_034_78);
    np_assert_allclose!(&[y], &[desired], rtol = 1e-5);
    let h2 = xsf::hankel_2(86.5, c64(1e15, 700.0));
    np_assert_allclose!([h2.norm()], [2.559_041_882_177_221e296], rtol = 1e-12);
    // Arb: bessely(-86.5, 1e15 + 10j) = sin(86.5 pi) besselj(86.5, 1e15 + 10j)
    let y = xsf::bessel_y(-86.5, c64(1e15, 10.0));
    np_assert_allclose!(
        &[y],
        &[c64(-2.384_955_233_412_051_6e-4, 1.426_054_857_376_210_6e-4)],
        rtol = 1e-13
    );
}

#[test]
fn test_bessel_amos_continuation_phase_near_axes() {
    // scipy/xsf#516; Arb: besselk(0, -1 + 1e9j) exp(-1 + 1e9j), hankel1(0, 1e15 - 1j)
    // exp(-1j (1e15 - 1j)), hankel2(0.25, -1e15) exp(-1e15j), and bessely(0.25, -1e15)
    let actual = [
        xsf::bessel_ke(0.0, c64(-1.0, 1e9)),
        xsf::hankel_1e(0.0, c64(1e15, -1.0)),
        xsf::hankel_2e(0.25, c64(-1e15, 0.0)),
        xsf::bessel_ye(0.25, c64(-1e15, 0.0)),
    ];
    let expected = [
        c64(2.802_495_607_148_028_4e-5, -2.802_495_609_249_900_3e-5),
        c64(1.784_124_116_152_772e-8, -1.784_124_116_152_770_6e-8),
        c64(-1.219_246_593_649_027_8e-8, -1.322_834_043_959_198_8e-8),
        c64(1.431_895_141_248_916e-8, 6.967_430_772_425_719e-9),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-13);
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
        assert!((actual / expected - 1.0).abs() < rtol, "{actual:e} != {expected:e}");
    }
}

xsref::tests! {
    test_it1i0k0: "it1i0k0", "d-d_d", |x| xsf::it1i0k0(x[0]);
    test_it2i0k0: "it2i0k0", "d-d_d", |x| xsf::it2i0k0(x[0]);
}

#[test]
fn test_iti0k0_series_near_40() {
    // the power series up to x = 40 (mpmath: `x * hyp1f2(1/2, 1, 3/2, x^2 / 4)` and
    // `x^2 / 8 * hyp2f3(1, 1, 2, 2, 2, x^2 / 4)`)
    for (actual, expected) in [
        (xsf::it1i0k0(35.5).0, 178_319_243_668_161.0),
        (xsf::it1i0k0(39.25).0, 7_197_960_272_196_551.0),
        (xsf::it1i0k0(39.75).0, 11_790_021_865_325_739.0),
        (xsf::it2i0k0(35.5).0, 5_175_762_522_389.534),
        (xsf::it2i0k0(39.25).0, 188_389_991_318_682.9),
        (xsf::it2i0k0(39.75).0, 304_585_988_297_001.7),
    ] {
        assert!((actual / expected - 1.0).abs() < 1e-14, "{actual:e} != {expected:e}");
    }
}

#[test]
fn test_iti0k0_large_x() {
    // unlike in xsf: quadrature and more asymptotic terms (mpmath)
    for (actual, expected, rtol) in [
        (xsf::it1i0k0(20.0).0, 44_758_598.913_855_75, 1e-12),
        (xsf::it1i0k0(30.0).0, 795_538_858_184.723_6, 1e-12),
        (xsf::it1i0k0(40.0).0, 1.508_971_508_271_920_2e16, 1e-14),
        (xsf::it2i0k0(40.0).0, 387_326_525_512_364.4, 1e-14),
        (xsf::it2i0k0(12.014_948_790_755_168).1, 1.616_101_465_331_337_3e-7, 1e-14),
        (xsf::it2i0k0(13.5).1, 3.110_087_578_320_996e-8, 1e-14),
        (xsf::it2i0k0(1.0625).1, 0.184_084_183_639_630_33, 1e-14),
        (xsf::it2i0k0(30.0).1, 6.780_268_362_546_424e-16, 1e-14),
        (xsf::it2i0k0(39.5).1, 3.399_304_068_768_309e-20, 1e-14),
        (xsf::it2i0k0(40.0).1, 2.024_204_614_539_913_3e-20, 1e-14),
        (xsf::it2i0k0(100.0).1, 4.588_521_432_774_676e-47, 1e-14),
    ] {
        assert!((actual / expected - 1.0).abs() < rtol, "{actual:e} != {expected:e}");
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
        assert!((actual / expected - 1.0).abs() < 1e-15, "{actual:e} != {expected:e}");
    }
}

#[test]
fn test_bessel_j_negative_order_near_overflow() {
    // AMOS overflows here, and unlike in xsf, cephes `jv` doesn't compute `(x/2)^v / Gamma(v + 1)`
    // as `exp(v ln(x/2) - lgam(v + 1))`, which loses ~|v ln(x/2)| eps (Arb)
    let actual = [
        xsf::bessel_j(-106.25, 0.097),
        xsf::bessel_j(-148.25, 0.9),
        xsf::bessel_j(-150.25, 0.97),
        xsf::bessel_j(-144.75, 0.78),
    ];
    let expected = [
        3.405_745_826_291_524e307,
        3.498_014_073_972_006e307,
        4.953_096_151_443_244e307,
        5.631_648_726_887_811e307,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
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
        assert!((actual / expected - 1.0).abs() < 1e-14, "{actual:e} != {expected:e}");
    }
}

#[test]
fn test_besselpoly_large_order() {
    // mpmath: `a**v / gamma(v + 1) / (l + v + 1) * hyp1f2((l + v + 1) / 2, (l + v + 3) / 2, v + 1,
    // -a**2)`
    for (a, lambda, v, expected) in [
        (92.0, 4.0, 390.5, 3.275_700_191_214_943e-90),
        (64.0, 0.0, 262.5, 1.460_944_853_450_890_9e-58),
        (41.0, -1.0, 398.5, 2.303_483_041_345_722_4e-227),
        (42.23, 0.0, 484.0, 3.238_013_225_871_64e-309),
        (43.75, 0.0, 490.0, 3.625_046_885_527_423_5e-308),
    ] {
        let actual = xsf::besselpoly(a, lambda, v);
        assert!((actual / expected - 1.0).abs() < 1e-12, "{actual:e} != {expected:e}");
    }
}

#[test]
fn test_besselpoly_small_a() {
    // unlike in xsf, `a^v` isn't `exp(v ln(a))`, which loses ~|v ln(a)| eps (Arb)
    let actual = [
        xsf::besselpoly(0.001, 0.0, 60.5),
        xsf::besselpoly(1e-5, 2.0, 25.5),
        xsf::besselpoly(0.001, 0.0, -40.5),
        xsf::besselpoly(0.001, 3.0, -50.25),
    ];
    let expected = [
        7.928_208_231_877_119e-267,
        1.409_647_783_394_057_2e-155,
        -3.277_260_630_252_027e166,
        -4.418_299_195_223_619e211,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_besselpoly_gamma_overflow() {
    // `gamma(v + 1)` overflows for v > 171, outside the Neumann series (v >= a^2 - 1); mpmath:
    // `a**v / gamma(v + 1) / (l + v + 1) * hyp1f2((l + v + 1) / 2, (l + v + 3) / 2, v + 1, -a**2)`
    for (a, lambda, v, expected) in [
        (12.3149, 8.9376, 192.1873, 8.809_530_109_564_947e-151),
        (8.0, 1.0, 175.5, 8.200_286_437_680_841e-164),
        (15.0, 0.0, 300.0, 3.468_324_551_656_663_7e-265),
        (5.0, -0.5, 210.25, 9.512_379_993_975_827e-255),
    ] {
        let actual = xsf::besselpoly(a, lambda, v);
        assert!((actual / expected - 1.0).abs() < 1e-12, "{actual:e} != {expected:e}");
    }
    // `a^v` overflows, but `lambda` is infinite
    assert_eq!(xsf::besselpoly(1e-29, f64::INFINITY, -10.5), 0.0);
}

#[test]
fn test_besselpoly_odd_nu_large() {
    // `static_cast<int>(nu)` is UB in xsf for `|nu| >= 2^31`; this `nu` is odd
    let r = xsf::besselpoly(0.5, 0.0, -2_147_483_649.0);
    assert!(r == 0.0 && r.is_sign_negative(), "{r}");
}

#[test]
fn test_besselpoly_negative_a_and_poles() {
    // mpmath: `quad(lambda t: t**l * besselj(v, 2 * a * t), [0, 0.5, 1])`
    for (a, lambda, v, expected) in [
        (-1.0, -1.0, -1.0, 0.849_045_485_440_153_2),
        (-2.5, 0.5, 2.0, 0.205_682_344_823_084_3),
        (-3.0, 1.0, 3.0, -0.145_444_613_310_487_95),
        (-1.5, 0.0, -2.0, 0.236_483_111_652_664_03),
        (-0.5, 2.0, -1.0, 0.114_903_484_931_900_48),
    ] {
        let actual = xsf::besselpoly(a, lambda, v);
        assert!((actual / expected - 1.0).abs() < 1e-14, "{actual:e} != {expected:e}");
    }
    assert!(xsf::besselpoly(-1.0, 0.0, 0.5).is_nan());

    // lambda + |v| + 1 is a non-positive even integer: the integral diverges
    for (a, lambda, v, expected) in [
        (10.0, -10.0, 1.0, f64::INFINITY),
        (1.0, -10.0, 1.0, f64::INFINITY),
        (1.0, -10.0, -1.0, f64::NEG_INFINITY),
        (10.0, -10.0, -1.0, f64::NEG_INFINITY),
        (-1.0, -10.0, 1.0, f64::NEG_INFINITY),
        (2.0, -1.0, 0.0, f64::INFINITY),
    ] {
        assert_eq!(xsf::besselpoly(a, lambda, v), expected, "{a} {lambda} {v}");
    }
    // lambda + 1 rounds to an even integer, but is not a pole; reference value from mpmath
    np_assert_allclose!(
        [xsf::besselpoly(1.0, -18_014_398_509_481_984.0, 0.0)],
        [-1.242_843_490_019_327_6e-17],
        rtol = 1e-15
    );
    assert_eq!(xsf::besselpoly(1.0, -4_503_599_627_370_499.0, 0.0), f64::INFINITY);
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
    // computed with N = 2, since `rctj` writes at least 2 terms
    let (j, jp) = xsf::riccati_j::<2>(0.2);
    assert_eq!(xsf::riccati_j::<1>(0.2), ([j[0]], [jp[0]]));
    let (y, yp) = xsf::riccati_y::<2>(0.2);
    assert_eq!(xsf::riccati_y::<1>(0.2), ([y[0]], [yp[0]]));
}

#[test]
fn test_riccati_j_all_orders() {
    // Arb; NaN or 0 in xsf above `msta1(x, 200)` (here 177) and for `|x| < 1e-100`, and `rj[1]`
    // cancels for `N = 2`
    let (rj, dj) = xsf::riccati_j::<201>(10.0);
    np_assert_allclose!(
        [rj[200], dj[200]],
        [4.359_433_049_621_081e-236, 8.751_636_342_177_128e-235],
        rtol = 1e-14
    );
    let (rj, dj) = xsf::riccati_j::<3>(1e-150);
    np_assert_allclose!(
        [rj[0], rj[1], rj[2], dj[0], dj[1], dj[2]],
        [1e-150, 3.333_333_333_333_333e-301, 0.0, 1.0, 6.666_666_666_666_667e-151, 2e-301],
        rtol = 1e-15
    );
    let (rj, dj) = xsf::riccati_j::<2>(1e-5);
    np_assert_allclose!(
        [rj[1], dj[1]],
        [3.333_333_333_300_000_5e-11, 6.666_666_666_533_334e-6],
        rtol = 1e-15
    );
}

/// Ported from xsf's `test_rctj.cpp` (rtol 1e-13 instead of 1e-12)
#[test]
fn test_riccati_j_large_x() {
    fn check<const N: usize>(x: f64, rj_n: f64, dj_n: f64) {
        let (rj, dj) = xsf::riccati_j::<N>(x);
        np_assert_allclose!([rj[N - 1], dj[N - 1]], [rj_n, dj_n], rtol = 1e-13);
    }
    check::<3>(5000.0, 0.987_873_519_167_095_8, -0.155_261_148_876_167_32);
    check::<6>(10000.0, 0.951_695_947_040_993_8, -0.307_042_254_603_583_53);
    check::<31>(20000.0, -0.600_733_011_725_298, -0.799_449_513_798_901_8);
    check::<5501>(5000.0, 2.236_532_392_121_206_7e-65, 1.026_728_990_451_556e-65);
}

#[test]
fn test_riccati_y_near_max() {
    fn check<const N: usize>(x: f64, k: usize, expected: [f64; 2]) {
        let (ry, dy) = xsf::riccati_y::<N>(x);
        np_assert_allclose!([ry[k], dy[k]], expected, rtol = 1e-12);
    }
    check::<5>(1e-70, 2, [-3e140, 6e210]);
    check::<151>(1.0, 149, [-1.257_390_673_139_500_4e304, 1.873_469_766_108_334e306]);
    check::<1849>(1000.0, 1848, [-8.119_236_242_111_399e307, 1.261_687_052_821_619_4e308]);
    check::<11768>(9976.0, 11767, [-1.742_111_939_068_904_2e308, 1.089_604_757_174_725_4e308]);

    let (ry, dy) = xsf::riccati_y::<201>(0.5);
    assert_eq!((ry[200], dy[200]), (f64::NEG_INFINITY, f64::INFINITY));
}

#[test]
fn test_riccati_y_order_zero() {
    for x in [1e-70_f64, 1.0] {
        let (ry, dy) = xsf::riccati_y::<1>(x);
        assert_eq!((ry[0], dy[0]), (-x.cos(), x.sin()));
    }
}

#[test]
fn test_riccati_y_overflow() {
    // Arb
    let (ry, dy) = xsf::riccati_y::<3>(1e-61);
    np_assert_allclose!(
        [ry[0], ry[1], ry[2], dy[0], dy[1], dy[2]],
        [-1.0, -1e61, -3e122, 1e-61, 1e122, 6e183],
        rtol = 1e-15
    );
    let (ry, dy) = xsf::riccati_y::<152>(1.0);
    np_assert_allclose!(
        [ry[148], dy[148]],
        [-4.233_686_952_143_891_5e301, 6.265_713_172_700_605e303],
        rtol = 1e-14
    );
    assert_eq!((ry[151], dy[151]), (f64::NEG_INFINITY, f64::INFINITY));
    // `k ry[k - 1]` or `k ry[k]` overflows before the division by `x`
    let (ry, _) = xsf::riccati_y::<522>(100.0);
    np_assert_allclose!([ry[521]], [-1.940_212_970_375_218_5e306], rtol = 1e-14);
    let (_, dy) = xsf::riccati_y::<245>(10.0);
    np_assert_allclose!([dy[244]], [3.878_341_575_822_852_3e307], rtol = 1e-14);
    // the reordered products overflow too, but `ry[k - 1]` cancels them
    let (_, dy) = xsf::riccati_y::<245>(9.937_577_591_883_544);
    np_assert_allclose!([dy[244]], [1.796_242_455_728_229_5e308], rtol = 1e-14);
    let (_, dy) = xsf::riccati_y::<522>(99.439_790_658_471_24);
    np_assert_allclose!([dy[521]], [1.764_702_935_338_309_7e308], rtol = 1e-14);
    let (ry, dy) = xsf::riccati_y::<3>(0.0);
    assert_eq!(ry, [-1.0, f64::NEG_INFINITY, f64::NEG_INFINITY]);
    assert_eq!(dy, [0.0, f64::INFINITY, f64::INFINITY]);
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

/// Based on `scipy.special.tests.test_basic.TestBessel.test_jnjnp_zeros`
#[test]
fn test_bessel_j_jprime_zeros() {
    // as in SciPy 1.18.1
    let (_, n, m, t) = xsf::bessel_j_jprime_zeros::<8>();
    assert_eq!(n, [0, 1, 0, 2, 1, 0, 3, 2]);
    assert_eq!(m, [0, 1, 1, 1, 1, 1, 1, 1]);
    assert_eq!(t, [true, true, false, true, false, true, true, false]);

    let (z, n, _, t) = xsf::bessel_j_jprime_zeros::<1200>();
    for l in 0..1200 {
        let y = xsf::bessel_j_prime(f64::from(n[l]), z[l], u32::from(t[l]));
        assert!(y.abs() < 1e-13, "l = {l}: {y:e}");
    }
    // mpmath.besseljzero(32, 9, derivative=1)
    np_assert_allclose!(&[z[1199]], &[68.563_138_567_071_27], rtol = 1e-15);

    // the 210th zero, which SciPy 1.18.1 misses: mpmath.besseljzero(26, 1, derivative=1)
    let (z, n, m, t) = xsf::bessel_j_jprime_zeros::<210>();
    assert_eq!((n[209], m[209], t[209]), (26, 1, true));
    np_assert_allclose!(&[z[209]], &[28.418_074_834_560_365], rtol = 1e-15);

    // in SciPy 1.18.1, `jdzo` writes past the end of `zo` for these
    let _ = xsf::bessel_j_jprime_zeros::<569>();
    let _ = xsf::bessel_j_jprime_zeros::<572>();
    assert_eq!(xsf::bessel_j_jprime_zeros::<0>(), ([], [], [], []));
}

#[test]
#[should_panic(expected = "N must be at most 1200")]
fn test_bessel_j_jprime_zeros_n1201() {
    let _ = xsf::bessel_j_jprime_zeros::<1201>();
}

/// Checks `f(z_k) = 0` and `g(z_k) = v_k` for the real and complex zeros `z_k` and values `v_k`
fn check_y_zeros(
    zeros: impl Fn(bool) -> [[Complex64; 20]; 2],
    f: impl Fn(Complex64) -> Complex64,
    g: impl Fn(Complex64) -> Complex64,
) {
    for complex in [false, true] {
        let [z, v] = zeros(complex);
        np_assert_allclose!(&z.map(&f), &[c64(0.0, 0.0); 20], atol = 5e-12);
        np_assert_allclose!(&z.map(&g), &v, rtol = 1e-11);
    }
}

/// Based on `scipy.special.tests.test_basic.TestBessel.test_y0_zeros`
#[test]
fn test_bessel_y_zeros() {
    let y0 = |z| xsf::bessel_y(0.0, z);
    let y1 = |z| xsf::bessel_y(1.0, z);
    let y1p = |z| xsf::bessel_y_prime(1.0, z, 1);
    check_y_zeros(xsf::bessel_y0_zeros, y0, y1);
    check_y_zeros(xsf::bessel_y1_zeros, y1, y0);
    check_y_zeros(xsf::bessel_y1_prime_zeros, y1p, y1);

    assert_eq!(xsf::bessel_y0_zeros::<0>(false), [[]; 2]);
}

#[test]
fn test_bessel_y_zeros_mpmath() {
    // mpmath.findroot(lambda z: mpmath.bessely(0, z), z0), and bessely(1, z) there
    let [z, v] = xsf::bessel_y0_zeros::<5>(false);
    let (zr, vr) = (13.361_097_473_872_763, -0.218_358_296_597_671_34);
    // in SciPy 1.18.1: 13.361097493, with relative error 1.4e-9
    np_assert_allclose!(&[z[4].re], &[zr], rtol = 2e-15);
    np_assert_allclose!(&[v[4].re], &[vr], rtol = 2e-15);

    // complex zeros near |z| = 12, as above
    let [z, v] = xsf::bessel_y0_zeros::<4>(true);
    let zr = c64(-11.791_512_030_381_03, 0.548_819_118_378_619_2);
    let vr = c64(-0.009_373_681_680_503_708, 0.402_304_542_938_848_1);
    np_assert_allclose!(&[z[3]], &[zr], rtol = 2e-15);
    np_assert_allclose!(&[v[3]], &[vr], rtol = 2e-15);
    // mpmath.findroot(lambda z: mpmath.bessely(1, z, 1), z0), and bessely(1, z) there
    let [z, v] = xsf::bessel_y1_prime_zeros::<5>(true);
    let zr = c64(-11.706_175_218_894_689, 0.552_859_060_733_942_4);
    let vr = c64(-0.009_544_376_760_384_838, 0.403_753_339_564_625_7);
    np_assert_allclose!(&[z[4]], &[zr], rtol = 2e-15);
    np_assert_allclose!(&[v[4]], &[vr], rtol = 2e-15);

    // the 200th zero, which is NaN in SciPy 1.18.1
    let [z, v] = xsf::bessel_y1_zeros::<200>(true);
    let zr = c64(-625.961_737_150_745_7, 0.549_306_670_047_127_9);
    let vr = c64(2.423_609_448_416_835e-5, -0.055_236_539_410_424_08);
    np_assert_allclose!(&[z[199]], &[zr], rtol = 1e-14);
    np_assert_allclose!(&[v[199]], &[vr], rtol = 1e-14);
}

#[test]
#[should_panic(expected = "v is too large")]
fn test_jahnke_emden_lambda_huge_v() {
    let _ = xsf::jahnke_emden_lambda(2_147_483_648.5, 1.0);
}

#[test]
fn test_jahnke_emden_lambda_tiny_x() {
    // Lambda_v(x) = 1 - x^2 / (4 (v + 1)) + O(x^4), so Lambda_v'(x) = -x / (2 (v + 1)) + O(x^3)
    for v in [3.0, 3.25] {
        for x in [0.0, 1e-101, 1e-50] {
            let (vl, dl) = xsf::jahnke_emden_lambda(v, x);
            let v0 = v - 3.0;
            let dlr: Vec<f64> = (0..4)
                .map(|k| -x / (2.0 * (v0 + f64::from(k) + 1.0)))
                .collect();
            assert_eq!(vl, [1.0; 4], "v = {v}, x = {x}");
            np_assert_allclose!(&dl, &dlr, rtol = 1e-15);
        }
    }
}

#[test]
fn test_jahnke_emden_lambda_high_orders() {
    // mpmath (v = 500) and Arb (v = 500.5); xsf stops at order `msta1(60, 200) = 328`
    for (v, expected) in [
        (
            500.0,
            [
                0.064_658_649_612_803_2,
                0.105_325_332_091_226_67,
                0.165_358_932_985_256_56,
                -0.005_927_303_686_869_632,
                -0.007_924_190_401_755_885,
                -0.009_937_421_238_798_044,
            ],
        ),
        (
            500.5,
            [
                0.064_929_345_091_574_51,
                0.105_622_626_797_108_82,
                0.165_656_964_513_909_68,
                -0.005_942_962_254_042_083,
                -0.007_936_549_366_129_714,
                -0.009_945_334_579_170_939,
            ],
        ),
    ] {
        let (vl, dl) = xsf::jahnke_emden_lambda(v, 60.0);
        let actual = [vl[329], vl[400], vl[500], dl[329], dl[400], dl[500]];
        np_assert_allclose!(actual, expected, rtol = 1e-13);
    }
}

#[test]
fn test_jahnke_emden_lambda_near_zero() {
    // Arb; the highest orders are close to a zero of `J_v(x)`, which the recurrence rounds to 0
    for (v, x, expected) in [
        (370.5, 383.972_410_969_705_27, -7.648_634_385_955_928e-56),
        (10.0, 18.433_463_666_966_58, 1.398_808_575_581_345_8e-4),
        (4.0, 7.588_342_434_503_804, -0.031_079_002_653_488_355),
    ] {
        let (_, dl) = xsf::jahnke_emden_lambda(v, x);
        np_assert_allclose!([dl[dl.len() - 1]], [expected], rtol = 1e-14);
    }
}

#[test]
fn test_jahnke_emden_lambda_subnormal() {
    // Arb; the derivatives are larger than the subnormal `Lambda_{v+1}(x)` they are formed from
    for (v, x, expected) in [(42.5, 1e9, -4e-323), (169.0, 1e4, -1.5e-323)] {
        let (_, dl) = xsf::jahnke_emden_lambda(v, x);
        assert_eq!(dl[dl.len() - 1], expected, "v = {v}, x = {x}");
    }
}

#[test]
fn test_jahnke_emden_lambda_negative_x() {
    // mpmath: Lambda_v(-60) and Lambda_v'(-60) for v = 2 and 1.5
    for (v, expected) in [
        (2.0, [2.067_224_078_837_053_6e-4, -8.977_047_004_812_257e-5]),
        (1.5, [7.894_440_028_306_55e-4, 2.934_810_510_600_466_6e-4]),
    ] {
        let (vl, dl) = xsf::jahnke_emden_lambda(v, -60.0);
        let actual = [vl[vl.len() - 1], dl[dl.len() - 1]];
        np_assert_allclose!(actual, expected, rtol = 1e-13);
    }
    // Lambda_v is even, so Lambda_v' is odd
    for v in [2.0, 1.5, 500.0] {
        for x in [2.0, 20.0, 60.0] {
            let ((vl, dl), (vln, dln)) =
                (xsf::jahnke_emden_lambda(v, x), xsf::jahnke_emden_lambda(v, -x));
            assert_eq!(vln, vl, "v = {v}, x = {x}");
            assert_eq!(dln, dl.iter().map(|d| -d).collect::<Vec<_>>(), "v = {v}, x = {x}");
        }
    }
}

#[test]
fn test_jahnke_emden_lambda_large_x() {
    for (v, x, k, expected) in [
        (2.0, 5000.0, 2, [2.126_507_932_532_124_8e-9, -2.919_271_034_712_797e-9]),
        (2.0, 9215.6, 2, [6.780_524_556_275_264e-10, -3.916_157_000_109_374e-10]),
        (4600.5, 5000.0, 2, [1.185_448_223_000_514_9e-10, -1.870_246_475_852_011e-11]),
        (1480.0, 800.0, 856, [1.798_396_960_028_322_8e-95, -1.232_184_002_870_922_6e-95]),
        (1480.5, 800.0, 856, [2.179_471_070_998_197_6e-95, -1.490_900_613_839_570_3e-95]),
    ] {
        let (vl, dl) = xsf::jahnke_emden_lambda(v, x);
        np_assert_allclose!([vl[k], dl[k]], expected, rtol = 1e-13);
    }
}

#[test]
fn test_jahnke_emden_lambda_large_v_and_x() {
    // Arb; in xsf, the scale of `J_v(x)` underflows around order `x / 2`, and `J_v(x)` itself
    // for `v > x`, which made these orders of magnitude off or 0
    for (v, x, expected) in [
        (2542.59, 1764.81, [2.365_382_472_695_595_7e-143, -9.539_856_558_966_074e-144]),
        (4485.0, 3156.13, [5.456_508_530_235_721e-260, -2.244_000_327_982_604_5e-260]),
        (4421.5, 3132.18, [4.873_607_538_273_819e-260, -2.023_149_554_161_378_5e-260]),
        // above `msta1(x, 200)`, where the power series cancels
        (5016.0, 3363.0, [7.165_465_027_825_305e-262, -2.757_026_551_415_938_5e-262]),
    ] {
        let (vl, dl) = xsf::jahnke_emden_lambda(v, x);
        np_assert_allclose!([vl[vl.len() - 1], dl[dl.len() - 1]], expected, rtol = 1e-13);
    }
}

#[test]
fn test_jahnke_emden_lambda_series_cancellation() {
    // Arb; xsf sums the power series up to `x = 12`, which loses up to about 1e-12 there
    for (v, expected) in [
        (1.0, [-0.038_484_579_775_113_284, 0.010_677_986_802_471_08]),
        (1.5, [-0.017_753_292_661_966_142, -0.008_619_580_458_152_027]),
    ] {
        let (vl, dl) = xsf::jahnke_emden_lambda(v, 11.9);
        np_assert_allclose!([vl[1], dl[1]], expected, rtol = 1e-14);
    }
}

#[test]
fn test_jahnke_emden_lambda_large_x_seeds() {
    // Arb; xsf starts the recurrences from the asymptotic expansions of `J_v0(x)` and
    // `J_{v0+1}(x)`, which lose up to about 4e-8 here
    for (v, x, expected) in [
        (2.25, 12.3, [-0.003_104_530_913_754_87, -0.008_465_513_461_267_032]),
        (13.25, 674_667.78, [-1.015_759_064_595_633_6e-67, -6.505_120_550_041_504e-67]),
        // `Lambda_41(x)` is subnormal
        (40.25, 1e9, [3.152_575_638_733_010_5e-307, 2.160_038_828_520_150_6e-307]),
    ] {
        let (vl, dl) = xsf::jahnke_emden_lambda(v, x);
        np_assert_allclose!([vl[vl.len() - 1], dl[dl.len() - 1]], expected, rtol = 1e-14);
    }
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

#[test]
fn test_amos_large_z() {
    // based on scipy/xsf#509: hankel_1 returned 0, the others NaN (mpmath)
    let h1 = xsf::hankel_1(86.0, c64(1e9, 0.0));
    np_assert_allclose!(
        [h1.re, h1.im],
        [-2.468_749_115_424_341e-5, 5.210_331_359_590_961e-6],
        rtol = 1e-6
    );
    let actual = [xsf::bessel_ie(0.0, 2e9), xsf::bessel_ke(10.0, 2e9)];
    np_assert_allclose!(
        actual,
        [8.920_620_581_321_394e-6, 2.802_495_678_086_199_5e-5],
        rtol = 5e-15
    );
    let [ai, _, bi, _] = xsf::airy_scaled(2e6);
    np_assert_allclose!(
        [ai, bi],
        [0.007_501_317_867_926_573, 0.015_002_635_736_958_194],
        rtol = 5e-15
    );
}

#[test]
fn test_amos_uniform_expansions_large_z() {
    // Regression test for scipy/xsf#502: these returned 0.
    let y = xsf::bessel_y(86.0, 2e15);
    np_assert_allclose!([y], [5.142_851_446_965_517e-9], rtol = 5e-15);
    for value in [
        xsf::hankel_1(86.0, c64(2e15, 0.0)),
        xsf::bessel_ke(86.0, c64(0.0, 2e15)),
        xsf::bessel_j(1e8, c64(2e15, 0.0)),
    ] {
        assert!(value.re.is_nan(), "{value}");
    }

    // Only the airy argument of the analytic continuation in unk2 is past the cutoff
    // (up to rounding).
    let ke = xsf::bessel_ke(86.0, c64(-400.0, 1_501_199_875_790_300.8));
    assert_ne!(ke, c64(0.0, 0.0));
}

#[test]
fn test_amos_orders_above_2_30() {
    // based on scipy/xsf#513: the first was NaN (J_v(v) and Y_v(v) from DLMF 10.19.8)
    let v30 = 2.0_f64.powi(30);
    assert_eq!(xsf::bessel_ie(v30, 1.0), 0.0);
    assert_eq!(xsf::bessel_k(v30, v30 - 1.0), 0.0);
    assert_eq!(xsf::bessel_j(2.0 * v30, -1.0), 0.0);
    let (v, j, y) = (1e15 + 1.0, 4.473_073_183_964_722e-6, -7.747_590_020_600_784e-6);
    let actual = [xsf::bessel_j(v, c64(-v, 0.0)), xsf::hankel_1(v, c64(-v, 0.0))];
    np_assert_allclose!(&actual, &[c64(-j, 0.0), c64(j, -y)], rtol = 1e-13);

    // NaN before; I_v(ix) = e^(i pi v / 2) J_v(x) and K_v(-ix) = i pi / 2 e^(i pi v / 2) H1_v(x)
    // (DLMF 10.27.6, 10.27.8), with J and Y from DLMF 10.19.6 (60 digits). The phase is only
    // accurate to about `eps v`.
    let (v, x) = (2.0 * v30 + 0.5, 4.0 * v30 + 2.0);
    let actual = [xsf::bessel_i(v, c64(0.0, x)), xsf::bessel_k(v, c64(0.0, -x))];
    let i = -9.107_019_807_551_605e-6;
    let desired = [c64(i, i), c64(1.685_736_916_683_522e-5, -1.175_317_735_666_563_8e-5)];
    np_assert_allclose!(&actual, &desired, rtol = 1e-5);
}
