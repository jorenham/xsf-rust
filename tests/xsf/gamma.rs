use num_complex::c64;

xsref::tests! {
    test_gamma_f64: "gamma", "d-d", |x| xsf::gamma(x[0]);
    test_gamma_c64: "gamma", "cd-cd", |x| xsf::gamma(c64(x[0], x[1]));
}

#[test]
fn test_gamma_c64_real_axis() {
    let w = xsf::gamma(c64(-2.0, 0.0));
    assert!(w.re.is_nan() && w.im.is_nan());
    // overflow: inf + 0i instead of inf + NaN i
    let w = xsf::gamma(c64(200.0, 0.0));
    assert_eq!(w.re, f64::INFINITY);
    assert_eq!(w.im.to_bits(), 0.0_f64.to_bits());
}

xsref::tests! {
    test_gammainc: "gammainc", "d_d-d", |x| xsf::gammainc(x[0], x[1]);
    test_gammaincc: "gammaincc", "d_d-d", |x| xsf::gammaincc(x[0], x[1]);
}

#[test]
fn test_gammainc_large_a() {
    // unlike in xsf (scipy/xsf#339): accurate near x ~ a for large `a`, where the power series
    // needs more terms than its iteration limit (mpmath)
    for (a, x, expected) in [
        (300_000.0, 297_500.0, 2.364_252_205_564_727e-6),
        (1e6, 995_000.0, 2.749_580_359_270_071e-7),
        (1e8, 99_950_000.0, 2.854_642_139_958_626e-7),
        (1e10, 9_999_500_000.0, 2.865_326_545_108_891e-7),
        (1e12, 999_970_000_000.0, 4.862_750_805_536_816e-198),
    ] {
        let actual = xsf::gammainc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-13,
            "{actual:e} != {expected:e}"
        );
    }
    let q = xsf::gammaincc(1e8, 99_950_000.0);
    assert!((q / 0.999_999_714_535_786 - 1.0).abs() < 1e-15, "{q}");
    let x = xsf::gammaincinv(1e10, 2.865_326_545_108_891e-7);
    assert!((x / 9_999_500_000.0 - 1.0).abs() < 1e-15, "{x}");
}

#[test]
fn test_gammainc_large_a_tails() {
    // unlike in xsf (scipy/scipy#26535): accurate for |x - a| > 0.4 a, where these were ~2e-11 off
    // (mpmath, Arb)
    for (a, x, expected) in [
        (
            4_539.614_661_308_867,
            2_638.637_810_126_278_5,
            1.053_583_186_114_541_3e-246,
        ),
        (
            5_965.501_074_765_158,
            3_523.337_672_155_038_5,
            2.834_471_063_272_162_8e-306,
        ),
        (
            4_532.911_009_098_821,
            2_712.174_570_596_431_7,
            6.239_827_729_196_435e-223,
        ),
        (
            4_600.020_859_664_592_5,
            2_684.110_636_552_072_3,
            1.001_757_841_288_208e-246,
        ),
    ] {
        let actual = xsf::gammainc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 3e-13,
            "{actual:e} != {expected:e}"
        );
    }
    let q = xsf::gammaincc(4_757.294_281_037_745, 7_144.513_767_033_158);
    let expected = 3.186_365_208_226_036_4e-199; // was 1.4e-11 off
    assert!((q / expected - 1.0).abs() < 3e-13, "{q:e} != {expected:e}");
}

#[test]
fn test_gammainc_huge_a() {
    // Arb: min(P, Q) < 10^(-6e306)
    assert_eq!(xsf::gammainc(f64::MAX, 3.0), 0.0);
    assert_eq!(xsf::gammaincc(4.4e306, 3.7e180), 1.0);
    assert_eq!(xsf::gammaincc(2.532_737_276_080_075_8e305, f64::MAX), 0.0); // MAX / ln(MAX)
    assert_eq!(xsf::chdtr(1e308, 3.0), 0.0);
    assert_eq!(xsf::pdtr(1e307, 1e306), 1.0);
}

#[test]
fn test_gammainc_tiny_a() {
    // unlike in xsf, these are accurate for tiny `a` (Arb)
    assert_eq!(xsf::gammainc(1e-310, 0.5), 1.0); // was 0
    assert_eq!(xsf::gammainc(1e-300, 1.0), 1.0); // was 1.0000000000000238
    for (a, x, expected) in [
        (3e-309, 1e-4, 2.589_967_411_372_412e-308),
        (1e-300, 1.0, 2.193_839_343_955_203e-301),
        (1e-200, 2.0, 4.890_051_070_806_111_5e-202),
    ] {
        let actual = xsf::gammaincc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-14,
            "{actual:e} != {expected:e}"
        );
    }
    let x = xsf::gammainccinv(1e-200, 1e-210);
    assert!((x / 19.984_172_765_540_563 - 1.0).abs() < 1e-15, "{x}");
}

#[test]
fn test_gammainc_tiny_a_near_1() {
    // unlike in xsf, these are not above 1 (mpmath: 1 - 5.2e-17)
    assert_eq!(xsf::gammainc(2e-16, 0.9), 1.0); // was 1 + 4.4e-16
    assert_eq!(xsf::chdtr(4e-16, 1.8), 1.0);
    assert_eq!(xsf::gdtr(2.0, 2e-16, 0.45), 1.0);
    // P = 1 - O(a) (mpmath)
    for (a, x, expected) in [
        (3e-14, 0.5, 0.999_999_999_999_983_2),   // was 1.9e-16 off
        (1e-9, 1e-100, 0.999_999_770_318_732_8), // was 1.9e-16 off
        (1e-6, 0.3, 0.999_999_094_323_186_3),    // was 4.1e-16 off
    ] {
        let actual = xsf::gammainc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-16,
            "{actual} != {expected}"
        );
    }
}

#[test]
fn test_gammainc_small_a() {
    // unlike in xsf, these are accurate for a < 1 (mpmath)
    for (a, x, expected) in [
        (1e-10, 1e-50, 0.999_999_988_544_796_2), // was 4e-15 off
        (1e-12, 0.3, 0.999_999_999_999_094_3),   // was 1.9e-15 off
    ] {
        let actual = xsf::gammainc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 4e-16,
            "{actual} != {expected}"
        );
    }
    for (a, x, expected) in [
        (1e-10, 1.1, 1.859_909_045_564_036_6e-11), // was 1.7e-14 off
        (1e-6, 1.0, 2.193_841_588_705_092_2e-7),   // was 4.7e-15 off
    ] {
        let actual = xsf::gammaincc(a, x);
        assert!(
            (actual / expected - 1.0).abs() < 2e-15,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn test_gammainc_far_from_a() {
    // unlike in xsf, `x^a e^-x / Gamma(a)` isn't `exp(a ln(x) - x - lgam(a))` for |x - a| > 0.4 a,
    // which loses ~eps times its largest term (Arb)
    let actual = [
        xsf::gammainc(120.0, 10.0),
        xsf::gammainc(30.0, 1e-5),
        xsf::gammaincc(100.0, 700.0),
        xsf::gammaincc(160.0, 300.0),
        xsf::gammaincc(0.5, 600.0),
        xsf::gammaincc(0.99, 709.0), // e^-x is subnormal
    ];
    let expected = [
        7.397_669_455_489_574e-84,
        3.769_951_145_241_386e-183,
        5.684_208_283_879_768e-179,
        2.686_206_949_999_588e-19,
        6.099_568_814_808_433_4e-263,
        1.132_809_914_017_541_5e-308,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

xsref::tests! {
    test_gammaincinv: "gammaincinv", "d_d-d", |x| xsf::gammaincinv(x[0], x[1]);
    test_gammainccinv: "gammainccinv", "d_d-d", |x| xsf::gammainccinv(x[0], x[1]);
    test_gammaln: "gammaln", "d-d", |x| xsf::gammaln(x[0]);
    test_gammasgn: "gammasgn", "d-d", |x| xsf::gammasgn(x[0]);
}

#[test]
fn test_gammasgn_large_negative() {
    // unlike in xsf, also for `x < -2^31`: Gamma(x) < 0 if floor(x) is odd, else Gamma(x) > 0
    for (x, s) in [(-2_147_483_648.5, -1.0), (-2_147_483_649.5, 1.0)] {
        assert_eq!(xsf::gammasgn(x), s);
        assert_eq!(xsf::gamma(x).signum(), s);
    }
}

#[test]
fn test_gammaln_near_zeros() {
    // mpmath
    let x = [
        0.999_999_999_999_999_9,
        0.998_622_470_062_837_7,
        1.999_999_999_999_969_1,
        1.15,
        1.85,
        1.999_999_999_9,
        1.000_000_001,
        0.999_999_999,
    ];
    let expected = [
        6.408_381_213_480_008e-17,
        0.000_796_693_611_134_219_1,
        -1.304_890_031_110_462_7e-14,
        -0.069_306_208_671_046_85,
        -0.055_923_813_019_657_244,
        -4.227_843_700_475_531_7e-11,
        -5.772_157_118_381_039e-10,
        5.772_156_493_992_26e-10,
    ];
    xsref::np_assert_allclose!(&x.map(xsf::gammaln), &expected, rtol = 1e-15);
    assert_eq!(x.map(xsf::loggamma), x.map(xsf::gammaln));
    // the zeros are +0
    assert_eq!([1.0, 2.0].map(|x| xsf::gammaln(x).to_bits()), [0, 0]);
    assert_eq!([1.0, 2.0].map(|x| xsf::loggamma(x).to_bits()), [0, 0]);
}
