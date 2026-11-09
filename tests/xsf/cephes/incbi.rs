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
        xsref::np_assert_allclose!(&[x], &[ref_], rtol = 1e-10);
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
