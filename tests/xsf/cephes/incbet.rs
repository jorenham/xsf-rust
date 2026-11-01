// based on scipy.special.tests.test_basic.TestBetaInc

#[test]
#[allow(clippy::float_cmp)]
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
