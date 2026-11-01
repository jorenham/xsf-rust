use num_complex::c64;
use xsref::np_assert_allclose;

// Bessel J

#[test]
fn test_bessel_j0() {
    xsref::test("cyl_bessel_j0", "d-d", |x| xsf::bessel_j0(x[0]));
}

#[test]
fn test_bessel_j1() {
    xsref::test("cyl_bessel_j1", "d-d", |x| xsf::bessel_j1(x[0]));
}

#[test]
fn test_bessel_j_f64() {
    xsref::test("cyl_bessel_j", "d_d-d", |x| xsf::bessel_j(x[0], x[1]));
}

#[test]
fn test_bessel_j_c64() {
    xsref::test("cyl_bessel_j", "d_cd-cd", |x| {
        xsf::bessel_j(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_bessel_je_f64() {
    xsref::test("cyl_bessel_je", "d_d-d", |x| xsf::bessel_je(x[0], x[1]));
}

#[test]
fn test_bessel_je_c64() {
    xsref::test("cyl_bessel_je", "d_cd-cd", |x| {
        xsf::bessel_je(x[0], c64(x[1], x[2]))
    });
}

// Bessel Y

#[test]
fn test_bessel_y0() {
    xsref::test("cyl_bessel_y0", "d-d", |x| xsf::bessel_y0(x[0]));
}

#[test]
fn test_bessel_y1() {
    xsref::test("cyl_bessel_y1", "d-d", |x| xsf::bessel_y1(x[0]));
}

#[test]
fn test_bessel_y_f64() {
    xsref::test("cyl_bessel_y", "d_d-d", |x| xsf::bessel_y(x[0], x[1]));
}

#[test]
fn test_bessel_y_c64() {
    xsref::test("cyl_bessel_y", "d_cd-cd", |x| {
        xsf::bessel_y(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_bessel_ye_f64() {
    xsref::test("cyl_bessel_ye", "d_d-d", |x| xsf::bessel_ye(x[0], x[1]));
}

#[test]
fn test_bessel_ye_c64() {
    xsref::test("cyl_bessel_ye", "d_cd-cd", |x| {
        xsf::bessel_ye(x[0], c64(x[1], x[2]))
    });
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

// Bessel I

#[test]
fn test_bessel_i0() {
    xsref::test("cyl_bessel_i0", "d-d", |x| xsf::bessel_i0(x[0]));
}

#[test]
fn test_bessel_i0e() {
    xsref::test("cyl_bessel_i0e", "d-d", |x| xsf::bessel_i0e(x[0]));
}

#[test]
fn test_bessel_i1() {
    xsref::test("cyl_bessel_i1", "d-d", |x| xsf::bessel_i1(x[0]));
}

#[test]
fn test_bessel_i1e() {
    xsref::test("cyl_bessel_i1e", "d-d", |x| xsf::bessel_i1e(x[0]));
}

#[test]
fn test_bessel_i_f64() {
    xsref::test("cyl_bessel_i", "d_d-d", |x| xsf::bessel_i(x[0], x[1]));
}

#[test]
fn test_bessel_i_c64() {
    xsref::test("cyl_bessel_i", "d_cd-cd", |x| {
        xsf::bessel_i(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_bessel_ie_f64() {
    xsref::test("cyl_bessel_ie", "d_d-d", |x| xsf::bessel_ie(x[0], x[1]));
}

#[test]
fn test_bessel_ie_c64() {
    xsref::test("cyl_bessel_ie", "d_cd-cd", |x| {
        xsf::bessel_ie(x[0], c64(x[1], x[2]))
    });
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
fn test_bessel_scaled_overflow() {
    // the signs of mpmath's besselj, bessely, and besseli, e.g. besselj(-11.25, 1e-100) = -3.6e1134
    let x = 1e-100;
    let z = c64(x, 0.0);
    let inf = f64::INFINITY;
    assert_eq!(xsf::bessel_ye(10.0, x), -inf);
    assert_eq!(xsf::bessel_ye(10.0, z), c64(-inf, 0.0));
    for (v, s) in [(-10.25, 1.0), (-11.25, -1.0)] {
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
    // mpmath: 1.9e308, 1.45e-21, and 1.99e-37 (exp(-2x) underflows, so no retry)
    for (v, x) in [
        (-2.000_000_000_000_000_4, 3.057_653_761_164_066_6e-162),
        (-804.25, 373.0),
        (-840.25, 400.0),
    ] {
        assert_eq!(xsf::bessel_ie(v, x), f64::INFINITY);
    }
}

// Bessel K

#[test]
fn test_bessel_k0() {
    xsref::test("cyl_bessel_k0", "d-d", |x| xsf::bessel_k0(x[0]));
}

#[test]
fn test_bessel_k0e() {
    xsref::test("cyl_bessel_k0e", "d-d", |x| xsf::bessel_k0e(x[0]));
}

#[test]
fn test_bessel_k1() {
    xsref::test("cyl_bessel_k1", "d-d", |x| xsf::bessel_k1(x[0]));
}

#[test]
fn test_bessel_k1e() {
    xsref::test("cyl_bessel_k1e", "d-d", |x| xsf::bessel_k1e(x[0]));
}

#[test]
fn test_bessel_k_f64() {
    xsref::test("cyl_bessel_k", "d_d-d", |x| xsf::bessel_k(x[0], x[1]));
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
    for (v, x) in [(1.5, 1e-210), (1.0, 1e-310), (-1.0, 1e-310), (0.5, 1e-306)] {
        assert_eq!(xsf::bessel_k(v, x), f64::INFINITY);
        assert_eq!(xsf::bessel_ke(v, x), f64::INFINITY);
    }
}

#[test]
fn test_bessel_k_c64() {
    xsref::test("cyl_bessel_k", "d_cd-cd", |x| {
        xsf::bessel_k(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_bessel_ke_f64() {
    xsref::test("cyl_bessel_ke", "d_d-d", |x| xsf::bessel_ke(x[0], x[1]));
}

#[test]
fn test_bessel_ke_c64() {
    xsref::test("cyl_bessel_ke", "d_cd-cd", |x| {
        xsf::bessel_ke(x[0], c64(x[1], x[2]))
    });
}

// Hankel

#[test]
fn test_hankel_1_c64() {
    xsref::test("cyl_hankel_1", "d_cd-cd", |x| {
        xsf::hankel_1(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_hankel_1e_c64() {
    xsref::test("cyl_hankel_1e", "d_cd-cd", |x| {
        xsf::hankel_1e(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_hankel_2_c64() {
    xsref::test("cyl_hankel_2", "d_cd-cd", |x| {
        xsf::hankel_2(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_hankel_2e_c64() {
    xsref::test("cyl_hankel_2e", "d_cd-cd", |x| {
        xsf::hankel_2e(x[0], c64(x[1], x[2]))
    });
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

// Bessel integrals

#[test]
fn test_it1j0y0() {
    xsref::test("it1j0y0", "d-d_d", |x| xsf::it1j0y0(x[0]));
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
fn test_it1i0k0() {
    xsref::test("it1i0k0", "d-d_d", |x| xsf::it1i0k0(x[0]));
}

#[test]
fn test_it2i0k0() {
    xsref::test("it2i0k0", "d-d_d", |x| xsf::it2i0k0(x[0]));
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
fn test_besselpoly() {
    xsref::test("besselpoly", "d_d_d-d", |x| {
        xsf::besselpoly(x[0], x[1], x[2])
    });
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
    #[allow(clippy::needless_range_loop)]
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
    #[allow(clippy::needless_range_loop)]
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
