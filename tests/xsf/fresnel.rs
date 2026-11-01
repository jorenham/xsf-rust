use num_complex::c64;

#[test]
fn test_fresnel_f64() {
    xsref::test("fresnel", "d-d_d", |x| xsf::fresnel(x[0]));
}

#[test]
#[allow(clippy::float_cmp)]
fn test_fresnel_f64_huge() {
    // x * x overflows; S(x) and C(x) round to +-1/2 (mpmath)
    for x in [1e155, -1e200, f64::MAX] {
        let half = 0.5_f64.copysign(x);
        assert_eq!(xsf::fresnel(x), (half, half));
    }
}

#[test]
fn test_fresnel_c64() {
    xsref::test("fresnel", "cd-cd_cd", |x| xsf::fresnel(c64(x[0], x[1])));
}

#[test]
fn test_modified_fresnel_plus_c64() {
    xsref::test("modified_fresnel_plus", "d-cd_cd", |x| {
        xsf::modified_fresnel_plus(x[0])
    });
}

#[test]
fn test_modified_fresnel_minus_c64() {
    xsref::test("modified_fresnel_minus", "d-cd_cd", |x| {
        xsf::modified_fresnel_minus(x[0])
    });
}

#[test]
fn test_modified_fresnel_negative_x() {
    // `ffk` in xsf computes the imaginary part of F±(x) from its real part if x < 0;
    // the reference values are from mpmath
    let (fp, _) = xsf::modified_fresnel_plus(-1.0);
    let (fm, _) = xsf::modified_fresnel_minus(-1.0);
    xsref::np_assert_allclose!(&[fp.im], &[0.936_925_370_381_131_2], rtol = 1e-15);
    xsref::np_assert_allclose!(&[fm.im], &[-0.936_925_370_381_131_2], rtol = 1e-15);
}

/// Ported from `scipy.special.tests.test_basic.TestFresnel.test_fresnel_zeros`
#[test]
fn test_fresnel_zeros() {
    // szo, czo = special.fresnel_zeros(5)
    let (szo, czo) = xsf::fresnel_zeros(5);
    xsref::np_assert_allclose!(
        &szo,
        &[
            c64(2.0093, 0.2885),
            c64(2.8335, 0.2443),
            c64(3.4675, 0.2185),
            c64(4.0026, 0.2009),
            c64(4.4742, 0.1877),
        ],
        atol = 1.5e-3
    );
    xsref::np_assert_allclose!(
        &czo,
        &[
            c64(1.7437, 0.3057),
            c64(2.6515, 0.2529),
            c64(3.3204, 0.2240),
            c64(3.8757, 0.2047),
            c64(4.3611, 0.1907),
        ],
        atol = 1.5e-3
    );

    // vals1 = special.fresnel(szo)[0]
    let vals1 = szo.iter().map(|&z| xsf::fresnel(z).0);
    // vals2 = special.fresnel(czo)[1]
    let vals2 = czo.iter().map(|&z| xsf::fresnel(z).1);
    // assert_allclose(vals1, 0, atol=1.5e-14, rtol=0)
    for val1 in vals1 {
        assert!(val1.norm() < 1.5e-14);
    }
    // assert_allclose(vals2, 0, atol=1.5e-14, rtol=0)
    for val2 in vals2 {
        assert!(val2.norm() < 1.5e-14);
    }
}
