//! based on scipy.special.tests.test_basic.TestBetaInc

#[test]
fn test_betaincinv_a1_b1() {
    // betaincinv(1, 1, x) is x.
    let x = [0.0, 0.25, 1.0];
    assert_eq!(x.map(|xv| xsf::betaincinv(1.0, 1.0, xv)), x);
}

#[test]
fn test_betaincinv_nontrivial() {
    for &(a, b, x, p) in &[
        (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
        (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
        // scipy/scipy#3761:
        (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
        // scipy/scipy#4244 (NOTE: this results in a relative error of 4.0e-12)
        // (
        //     0.0002742794749792665,
        //     289206.03125,
        //     1.639984034231756e-56,
        //     0.9688708782196045,
        // ),
        // scipy/scipy#12796:
        (4.0, 99997.0, 0.000_194_784_157_889_212_1, 0.999_995),
    ] {
        let x1 = xsf::betaincinv(a, b, p);
        xsref::np_assert_allclose!(&[x1], &[x], rtol = 5e-13);
    }
}

#[test]
fn test_betaincinv_tiny_y() {
    // Test with extremely small y values.  This test includes
    // a regression test for an issue in the boost code;
    // see https://github.com/boostorg/math/issues/961
    //
    // The reference values were computed with mpmath. For example,
    //
    //   from mpmath import mp
    //   mp.dps = 1000
    //   a = 14.208308325339239
    //   p = 7.703145458496392e-307
    //   x = mp.findroot(lambda t: mp.betainc(a, a, 0, t,
    //                                        regularized=True) - p,
    //                   x0=8.566e-23)
    //   print(float(x))

    for &(a, b, y, ref_) in &[
        (
            14.208_308_325_339_239,
            14.208_308_325_339_239,
            7.703_145_458_496_392e-307,
            8.566_004_561_846_704e-023,
        ),
        (14.0, 14.5, 1e-280, 2.934_391_500_664_242_4e-21),
        (3.5, 15.0, 4e-95, 1.329_075_142_928_922_7e-28),
        (10.0, 1.25, 2e-234, 3.982_659_092_143_654e-24),
        // NOTE: This results in a relative error of 3.806e-11
        // (4.0, 99997.0, 5e-88, 3.309800566862242e-27),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = 1e-14);
    }
}

#[test]
fn test_betaincinv_small_y_near_1() {
    // small y with a root in (0.75, 1 - y); reference values from mpmath
    for &(a, b, y, ref_) in &[
        (509.0, 0.041, 1.3e-48, 0.817_229_092_544_038_5),
        (5_818.0, 0.006, 7.65e-17, 0.995_090_852_864_442_7),
        (100.0, 0.5, 1e-12, 0.775_037_837_698_101_4),
        (10.0, 0.0026, 0.001_05, 0.932_881_143_861_390_6),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = 1e-14);
    }
}

#[test]
fn test_betaincinv_small_a_y_near_1() {
    // a <= 1 and y near 1, with the root far below where I_x(a, b) rounds to 1; reference values
    // from Arb (for a = 1, the root is 1 - (1 - y)^(1/b)), with a condition number of ~1e6
    for &(a, b, y, ref_) in &[
        (0.15, 2.7e6, 0.999_999_7, 4.103_212_725_816_707e-6),
        (0.5, 1e8, 0.999_999_9, 1.418_699_271_051_264e-7),
        (0.25, 1e10, 0.999_999_9, 1.286_131_058_351_140_3e-9),
        (1.0, 1e6, 0.999_999_9, 1.611_796_575_567_885_2e-5),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = 1e-9);
    }
}

#[test]
fn test_betaincinv_y_near_1() {
    // 1 - y = 1e-12 or 1e-13, of which I_x(a, b) near the root keeps only a few digits; reference
    // values from Arb
    for &(a, b, y, ref_) in &[
        (0.001, 100.0, 0.999_999_999_999_9, 0.181_972_188_908_984_5),
        (0.01, 1e3, 0.999_999_999_999, 0.019_829_317_529_540_336),
        (0.9, 1e3, 0.999_999_999_999_9, 0.029_094_942_107_102_114),
        (3.0, 1e3, 0.999_999_999_999, 0.033_446_364_134_901_65),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = 1e-14);
    }
}

#[test]
fn test_betaincinv_huge_a_x_near_1() {
    // huge a, b <= 1 and small y, with the root within ~1e-7 of 1; reference values from Arb
    for &(a, b, y, ref_) in &[
        (4e8, 0.4, 1e-12, 0.999_999_937_792_828_8),
        (4e9, 0.1, 1e-100, 0.999_999_944_216_303_9),
        (1e10, 0.5, 1e-20, 0.999_999_995_641_913_3),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = 1e-15);
    }
}

#[test]
fn test_betaincinv_tiny_a_or_b() {
    // reference values from Arb; the first root is below eps (1 - x rounds to 1), with a condition
    // number of ~2e15
    for &(a, b, y, ref_, rtol) in &[
        (
            4.637_934_404_623_266e-16,
            7.069_029_637_096_786,
            0.999_999_999_999_981_8,
            7.617_553_767_174_632e-19,
            5e-2,
        ),
        (7.07, 1e-10, 1e-10, 0.960_553_238_494_944_4, 1e-14),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = rtol);
    }
}

#[test]
fn test_betaincinv_y_max() {
    // y = 1 - 2^-53, the largest f64 below 1; reference values from Arb, or for b = 1e300 from the
    // limit Q(a, b x) = 1 - y. For a root near 0, the rounding of 1 - x alone moves it by up to
    // ~eps / x, and for a root below eps, the rounding of I_x(a, b) by up to ~1e-2.
    let y = 1.0 - f64::EPSILON / 2.0;
    for &(a, b, ref_, rtol) in &[
        (0.01, 1e3, 0.028_382_173_017_582_47, 1e-14),
        (0.001, 1e10, 2.651_928_480_581_199_2e-9, 1e-7),
        (100.0, 1e5, 0.002_051_315_441_176_488_3, 1e-13),
        (0.5, 1e300, 3.438_162_610_583_420_4e-299, 1e-2),
        (0.1, 3.0, 0.999_985_766_388_452_3, 1e-15),
        (1.0, 10.0, 0.974_617_112_613_867_6, 1e-15),
    ] {
        let x = xsf::betaincinv(a, b, y);
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = rtol);
    }
}

#[test]
fn test_betaincinv_a1_huge_b() {
    // for a = 1, the root is 1 - (1 - y)^(1/b)
    let x = xsf::betaincinv(1.0, 5e15, 0.781_25);
    xsref::np_assert_allclose!(&[x], &[3.039_651_507_488_826e-16], rtol = 1e-14);
}

#[test]
fn test_betaincinv_huge_a_b1() {
    // the root is y^(1/a), which rounds to 1; lgam(a + b) overflows
    let x = xsf::betaincinv(1e306, 1.0, 0.9);
    xsref::np_assert_allclose!(&[x], &[1.0], rtol = 1e-15);
}

#[test]
fn test_betaincinv_gh21426() {
    // Test for scipy/scipy#21426: betaincinv must not return NaN
    let a = 5.0;
    let x = 0.5;
    let result = xsf::betaincinv(a, a, x);
    xsref::np_assert_allclose!(&[result], &[0.5], rtol = 10.0 * f64::EPSILON);
}

#[test]
fn test_betaincinv_tiny_ab() {
    // I_x(a, a) = 1/2 at x = 1/2, and is 1/2 up to ~745 a for all 0 < x < 1
    assert_eq!(xsf::betaincinv(1e-300, 1e-300, 0.5), 0.5);
    assert_eq!(xsf::betaincinv(1e-308, 1e-308, 0.5), 0.5);
}

#[test]
fn test_betaincinv_zero_ab() {
    for y in [0.3, 0.7] {
        assert!(xsf::betaincinv(0.0, 0.5, y).is_nan());
        assert!(xsf::betaincinv(0.5, 0.0, y).is_nan());
    }
}

#[test]
fn test_betainccinv_a1_b1() {
    // betainccinv(1, 1, x) is 1 - x.
    let x = [0.0, 0.25, 1.0];
    assert_eq!(
        x.map(|xv| xsf::betainccinv(1.0, 1.0, xv)),
        x.map(|xv| 1.0 - xv)
    );
}

#[test]
fn test_betainccinv_nontrivial() {
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
        xsref::np_assert_allclose!(&[xsf::betainccinv(a, b, p)], &[x], rtol = 8e-15);
    }
}

#[test]
fn test_betainccinv_tails() {
    // roots far below eps, from 1 - I_x(1, 2) = (1 - x)^2 and I_x(1/2, 1/2) = 2 asin(sqrt(x)) / pi
    let y = 1.0 - f64::EPSILON / 2.0;
    assert_eq!(xsf::betainccinv(1.0, 2.0, y), f64::EPSILON / 4.0);
    xsref::np_assert_allclose!(
        &[xsf::betainccinv(0.5, 0.5, y)],
        &[3.041_306_664_850_247e-32],
        rtol = 1e-15
    );
    // roots near 0 from the complement of a root near 1; for a = 1 the root is
    // 1 - y^(1/b), and for a = 2.5 the reference value is from mpmath
    for &(a, b, y, ref_) in &[
        (1.0, 1e5, 1e-100, 0.002_299_936_177_446_682_6),
        (1.0, 1e5, 1e-10, 0.000_230_232_001_843_413_66),
        (2.5, 1e5, 0.3, 3.032_146_280_274_775_4e-5),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainccinv(a, b, y)], &[ref_], rtol = 1e-14);
    }
}

#[test]
fn test_betainccinv_far_from_incbi() {
    // the roots (1 - y)^(1/a) for b = 1 and 1 - y^(1/b) for a = 1, where the root from `incbi` is
    // off by orders of magnitude, or 0
    for &(a, b, y, ref_) in &[
        (1e-18, 1.0, 1e-16, 3.720_075_976_020_851_6e-44),
        (1.0, 1e306, 0.1, 2.302_585_092_994_045_5e-306),
        // `incbi` is NaN, as `lgam(a + b)` overflows; the root rounds to 1
        (1e306, 1.0, 0.9, 1.0),
    ] {
        xsref::np_assert_allclose!(&[xsf::betainccinv(a, b, y)], &[ref_], rtol = 5e-14);
    }
}

#[test]
fn test_betainccinv_extreme() {
    // the root (1 - y)^(1/a) = 1e-1000 underflows
    assert_eq!(xsf::betainccinv(0.001, 1.0, 0.9), 0.0);
    // 1 - I_x(a, 1) = -a ln(x) rounds to y over a range of x around the root 1/e, where any x is
    // accepted
    let x = xsf::betainccinv(5e-324, 1.0, 5e-324);
    assert_eq!(xsf::betaincc(5e-324, 1.0, x), 5e-324);
}

#[test]
fn test_betainccinv_edge_cases() {
    let inf = f64::INFINITY;
    for &((a, b, y), expected) in &[
        ((2.0, 3.0, 0.0), 1.0),
        ((2.0, 3.0, -0.0), 1.0),
        ((2.0, 3.0, 1.0), 0.0),
        ((0.0, 3.0, 0.0), f64::NAN),
        ((2.0, 0.0, 1.0), f64::NAN),
        ((inf, 2.0, 0.0), 1.0),
        ((2.0, inf, 1.0), 0.0),
        ((1.0, inf, 0.5), 0.0),
        ((inf, 1.0, 0.5), 1.0),
        ((2.0, inf, 0.5), f64::NAN),
        ((inf, inf, 0.5), f64::NAN),
        ((2.0, 3.0, -1e-300), f64::NAN),
        ((2.0, 3.0, 1.1), f64::NAN),
        ((2.0, 3.0, f64::NAN), f64::NAN),
    ] {
        xsref::np_assert_equal!(&[xsf::betainccinv(a, b, y)], &[expected]);
    }
}
