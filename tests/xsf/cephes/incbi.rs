//! based on scipy.special.tests.test_basic.TestBetaInc

#[test]
#[allow(clippy::float_cmp)]
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
        // scipy/scipy#4244 (NOTE: this results in a relative error of 1.872e-6)
        // (
        //     0.0002742794749792665,
        //     289206.03125,
        //     1.639984034231756e-56,
        //     0.9688708782196045,
        // ),
        // scipy/scipy#12796 (NOTE: this results in a relative error of 1.787e2)
        // (4.0, 99997.0, 0.0001947841578892121, 0.999995),
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
fn test_betaincinv_gh21426() {
    // Test for scipy/scipy#21426: betaincinv must not return NaN
    let a = 5.0;
    let x = 0.5;
    let result = xsf::betaincinv(a, a, x);
    xsref::np_assert_allclose!(&[result], &[0.5], rtol = 10.0 * f64::EPSILON);
}
