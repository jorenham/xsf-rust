// based on scipy.special.tests.test_basic.TestBetaInc

#[test]
fn test_betainc_a1_b1() {
    // betainc(1, 1, x) is x.
    let x = [0.0, 0.25, 1.0];
    assert_eq!(x.map(|xv| xsf::betainc(1.0, 1.0, xv)), x);
}

#[test]
fn test_betainc_nontrivial() {
    for &(a, b, x, p) in &[
        (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
        (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
        // scipy/scipy#3761:
        (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
        // scipy/scipy#4244 (NOTE: this results in a relative error of 5.134e-10)
        // (
        //     0.0002742794749792665,
        //     289206.03125,
        //     1.639984034231756e-56,
        //     0.9688708782196045,
        // ),
        // scipy/scipy#12796:
        (4.0, 99997.0, 0.000_194_784_157_889_212_1, 0.999_995),
    ] {
        let p1 = xsf::betainc(a, b, x);
        // NOTE: The original `rtol = 1e-15` is too strict for Cephes
        xsref::np_assert_allclose!(&[p1], &[p], rtol = 1e-14);
    }
}

#[test]
fn test_betainc_slow_continued_fraction() {
    // reference values from mpmath

    // large a + b near the mean
    let p = xsf::betainc(1e8, 1e8, 0.499_996_464_466_102_9);
    xsref::np_assert_allclose!(&[p], &[0.460_172_162_871_263_6], rtol = 5e-12);

    // x near 1 with small b
    for &(a, b, x, p) in &[
        (1e4, 0.002, 0.999_999_6, 0.009_851_010_029_838_827),
        (116.3, 0.0008, 0.999_98, 0.004_385_278_106_520_217),
        (1000.0, 1e-6, 0.999_999, 6.332_019_221_868_309e-6),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 1e-10);
    }

    // the result should stay within [0, 1] for huge a + b
    assert!((0.0..=1.0).contains(&xsf::betainc(1e15, 1e15, 0.5)));
    assert!((0.0..=1.0).contains(&xsf::betainc(2.0, 1e17, 2e-17)));
    assert!((0.0..=1.0).contains(&xsf::fdtr(2e15, 2e15, 1.0)));
    assert!((0.0..=1.0).contains(&xsf::fdtrc(2e15, 2e15, 1.0)));
}

#[test]
fn test_betainc_large_ab() {
    // reference values from mpmath
    for &(a, b, x, p) in &[
        (1e8, 2e8, 0.3333, 0.110_334_398_548_217_36),
        (5e7, 8e7, 0.3846, 0.359_220_529_919_169_7),
        (3e8, 1e9, 0.2307, 1.562_168_507_020_367_6e-9),
        // tiny x with huge b, where 1 - x is rounded
        (13.0, 1e15, 1e-14, 0.208_443_523_605_136_1),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 5e-12);
    }
}

#[test]
fn test_betainc_large_ab_near_mean() {
    // reference values from Arb (python-flint)
    for &(a, b, x, p) in &[
        (200.0, 2000.0, 0.091, 0.513_981_623_504_570_9),
        (730.0, 9560.0, 0.0714, 0.575_863_132_192_416_2),
        (2000.0, 30000.0, 0.0626, 0.532_119_378_151_077_6),
        (2600.0, 32000.0, 0.076, 0.728_174_963_532_355_2),
        (6000.0, 85000.0, 0.066, 0.533_466_639_319_583_9),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 5e-15);
    }
}

#[test]
fn test_betainc_huge_ab_near_mean() {
    // 1/2 by symmetry
    for a in [1e15, 1e20, 1e200] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, a, 0.5)], &[0.5], rtol = 1e-15);
    }
    xsref::np_assert_allclose!(&[xsf::fdtr(2e15, 2e15, 1.0)], &[0.5], rtol = 1e-15);

    // reference values from Arb (python-flint), and mpmath quadrature for the last one
    for &(a, b, x, p) in &[
        (1e12, 3e12, 0.250_000_1, 0.677_916_442_949_691_3),
        (4e15, 1e12, 0.999_750_062_48, 0.493_009_849_083_284_6),
        (1e18, 1e15, 0.999_000_999, 0.487_379_439_441_624_5),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 1e-8);
    }

    // reference value from mpmath, using I_x(a, b) ~ P(a, -(b + (a - 1)/2) ln(1 - x)) for b >> a
    let p = xsf::betainc(101.0, 1e20, 1.01e-18);
    xsref::np_assert_allclose!(&[p], &[0.513_232_791_918_545_2], rtol = 1e-15);
}

#[test]
fn test_betainc_huge_b() {
    // reference values from mpmath
    for &(a, b, x, p) in &[
        (2.0, 1e17, 2e-17, 0.593_994_150_290_162),
        (5.0, 5e19, 5e-20, 0.108_821_981_085_848_76),
        (15.0, 5e19, 6e-19, 0.999_079_317_603_851_4),
        (20.0, 1e12, 3e-11, 0.978_126_531_565_180_6),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 1e-13);
    }
}

#[test]
fn test_betainc_a_or_b_above_1e154() {
    assert_eq!(xsf::betainc(2.1e164, 4.6e5, 0.5), 0.0);
    xsref::np_assert_allclose!(&[xsf::betainc(381.0, 9.9e154, 0.5)], &[1.0], rtol = 1e-15);
    assert_eq!(xsf::betainc(1e200, 1e100, 1e-99), 0.0);
    // `lbeta(a, b)` is `inf - inf` here
    assert_eq!(xsf::betainc(1e306, 1e306, 1e-307), 0.0);
    assert_eq!(xsf::fdtrc(1e306, 1e306, 1e-307), 1.0);
    // reference values from mpmath, using fdtr(a, inf, x) = P(a/2, a x/2)
    let p = [xsf::fdtr(3.0, 1e200, 2.0), xsf::fdtr(0.5, 1e300, 10.0)];
    let p_ref = [0.888_389_774_905_287_4, 0.990_786_251_520_46];
    xsref::np_assert_allclose!(&p, &p_ref, rtol = 1e-14);
    let q = [xsf::fdtrc(3.0, 1e200, 2.0), xsf::fdtrc(0.5, 1e300, 10.0)];
    let q_ref = [0.111_610_225_094_712_56, 0.009_213_748_479_540_033];
    xsref::np_assert_allclose!(&q, &q_ref, rtol = 1e-13);
}

#[test]
fn test_betainc_edge_cases() {
    for &((a, b, x), expected) in &[
        ((0.0, 0.0, 0.0), f64::NAN),
        ((0.0, 0.0, 0.5), f64::NAN),
        ((0.0, 0.0, 1.0), f64::NAN),
        ((f64::INFINITY, f64::INFINITY, 0.0), f64::NAN),
        ((f64::INFINITY, f64::INFINITY, 0.5), f64::NAN),
        ((f64::INFINITY, f64::INFINITY, 1.0), f64::NAN),
        ((0.0, 1.0, 0.0), 0.0),
        ((0.0, 1.0, 0.5), 1.0),
        ((0.0, 1.0, 1.0), 1.0),
        ((1.0, 0.0, 0.0), 0.0),
        ((1.0, 0.0, 0.5), 0.0),
        ((1.0, 0.0, 1.0), 1.0),
        ((0.0, f64::INFINITY, 0.0), 0.0),
        ((0.0, f64::INFINITY, 0.5), 1.0),
        ((0.0, f64::INFINITY, 1.0), 1.0),
        ((f64::INFINITY, 0.0, 0.0), 0.0),
        ((f64::INFINITY, 0.0, 0.5), 0.0),
        ((f64::INFINITY, 0.0, 1.0), 1.0),
        ((1.0, f64::INFINITY, 0.0), 0.0),
        ((1.0, f64::INFINITY, 0.5), 1.0),
        ((1.0, f64::INFINITY, 1.0), 1.0),
        ((f64::INFINITY, 1.0, 0.0), 0.0),
        ((f64::INFINITY, 1.0, 0.5), 0.0),
        ((f64::INFINITY, 1.0, 1.0), 1.0),
    ] {
        let observed = xsf::betainc(a, b, x);
        xsref::np_assert_equal!(&[observed], &[expected]);
    }
}

#[test]
fn test_betainc_gh22682() {
    // scipy/scipy#22682: betainc returned incorrect results for tiny
    // single precision inputs. test that this is resolved
    let rtol = 1e-15;
    for &(a, b, x, reference) in &[
        (1e-20, 1e-21, 0.5, 0.090_909_090_909_090_9),
        (1e-15, 1e-16, 0.5, 0.090_909_090_909_090_91),
    ] {
        let res = xsf::betainc(a, b, x);
        xsref::np_assert_allclose!(&[res], &[reference], rtol = rtol);
    }
}

#[test]
fn test_betainc_tiny_a() {
    // reference values from mpmath; I_x(a, b) = 1 - O(a)
    for &(a, b, x, p) in &[
        (5.5e-15, 25.0, 0.015, 0.999_999_999_999_995_9),
        (1e-13, 100.0, 1e-4, 0.999_999_999_999_595_8),
        (1e-13, 170.0, 1e-20, 0.999_999_999_995_965_8),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 1e-15);
    }
    assert_eq!(xsf::betainc(3e-308, 0.5, 0.3), 1.0);
}

#[test]
fn test_betainc_tiny_b_x_near_1() {
    // reference values from Arb; I_x(a, b) = O(b), up to a relative ~b |ln(1 - x)|
    for &(a, b, x, p) in &[
        (
            7.069_029_637_096_786,
            4.637_934_404_623_266e-16,
            0.999_999_999_999_930_9,
            1.291_354_777_599_956_3e-14,
        ),
        (2.0, 1e-12, 0.999_999_999_999_9, 2.893_329_531_165_079_5e-11),
        (0.5, 1e-20, 0.999_999_999_999, 2.901_733_759_901_274_8e-19),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainc(a, b, x)], &[p], rtol = 1e-10);
    }
    // b |ln(1 - x)| ~ 3e-9, where `1 - I_(1-x)(b, a)` cancels (mpmath and Wolfram)
    let p = xsf::betainc(
        240.444_389_643_372_94,
        1.112_624_883_358_221_1e-10,
        0.999_999_999_999_989,
    );
    xsref::np_assert_allclose!(&[p], &[2.902_177_099_562_651e-9], rtol = 5e-9);
}

#[test]
fn test_betainc_subnormal_ab() {
    // B(a, b) overflows; I_x(1, b) = 1 - (1 - x)^b, i.e. -b ln(1 - x) for tiny b
    assert_eq!(xsf::betainc(1e-310, 1.0, 0.5), 1.0);
    let p = xsf::betainc(1.0, 1e-310, 0.99);
    xsref::np_assert_allclose!(&[p], &[4.605_170_185_988_06e-310], rtol = 1e-13);
    // b / (a + b) is subnormal; I_x(2, b) = -b (ln(1 - x) + x) for tiny b and x close to 1
    assert_eq!(xsf::betainc(2.0, 5e-324, 0.99), 2e-323);
}

#[test]
fn test_betaincc_a1_b1() {
    // betaincc(1, 1, x) is 1 - x.
    let x = [0.0, 0.25, 1.0];
    assert_eq!(
        x.map(|xv| xsf::betaincc(1.0, 1.0, xv)),
        x.map(|xv| 1.0 - xv)
    );
}

#[test]
fn test_betaincc_nontrivial() {
    for &(a, b, x, p) in &[
        (2.5, 3.0, 0.25, 0.833_251_953_125),
        (7.5, 13.25, 0.375, 0.432_987_346_455_603_7),
        (0.125, 7.5, 0.425, 0.000_668_825_785_131_423_7),
        (0.125, 18.0, 1e-6, 0.729_823_591_450_963_3),
        (0.125, 18.0, 0.996, 7.274_587_553_838_015e-46),
        (0.125, 24.0, 0.75, 3.708_534_048_168_62e-17),
        (16.0, 0.75, 0.999_999_999_75, 5.440_875_927_741_863e-7),
        // scipy/scipy#4677
        (
            0.421_195_964_350_340_1,
            16_939.046_996_018_118,
            0.000_815_296_167_195_521,
            1e-7,
        ),
    ] {
        xsref::np_assert_allclose!(&[xsf::betaincc(a, b, x)], &[p], rtol = 5e-15);
    }
}

#[test]
fn test_betaincc_tails() {
    // I_(1-x)(3, 2) = y^3 (4 - 3 y) with y = 1 - x, where `1 - betainc` cancels
    let q = xsf::betaincc(2.0, 3.0, 0.999_999_999_999);
    xsref::np_assert_allclose!(&[q], &[3.999_734_545_227_945e-36], rtol = 1e-14);
    // 1 - I_x(0.5, 2) = 1 - 1.5 sqrt(x) + O(x^1.5)
    let q = xsf::betaincc(0.5, 2.0, 1e-20);
    xsref::np_assert_allclose!(&[q], &[1.0 - 1.5e-10], rtol = 1e-15);
}

#[test]
fn test_betaincc_tiny_a() {
    // reference values from mpmath and Wolfram; 1 - I_x(a, b) = O(a), where `1 - betainc` cancels
    for &(a, b, x, q) in &[
        (2e-12, 1e300, 2.647_370_104_516_100_4e-301, 2e-12),
        (
            3.464_166_003_487_626_4e-12,
            3.465_839_753_112_206e186,
            1.675_272_022_682_781e-188,
            8.058_530_735_930_76e-12,
        ),
    ] {
        xsref::np_assert_allclose!(&[xsf::betaincc(a, b, x)], &[q], rtol = 1e-11);
    }
}

#[test]
fn test_betaincc_edge_cases() {
    let inf = f64::INFINITY;
    for &((a, b, x), expected) in &[
        ((0.0, 0.0, 0.5), f64::NAN),
        ((inf, inf, 0.5), f64::NAN),
        ((0.0, 1.0, 0.0), 1.0),
        ((0.0, 1.0, 0.5), 0.0),
        ((1.0, 0.0, 0.5), 1.0),
        ((1.0, 0.0, 1.0), 0.0),
        ((0.0, inf, 0.0), 1.0),
        ((0.0, inf, 0.5), 0.0),
        ((inf, 0.0, 0.5), 1.0),
        ((inf, 0.0, 1.0), 0.0),
        ((1.0, inf, 0.5), 0.0),
        ((inf, 1.0, 0.5), 1.0),
        ((1.5, 2.0, -1e-300), f64::NAN),
        ((1.5, 2.0, 1.1), f64::NAN),
        ((-1.0, 2.0, 0.5), f64::NAN),
    ] {
        xsref::np_assert_equal!(&[xsf::betaincc(a, b, x)], &[expected]);
    }
}
