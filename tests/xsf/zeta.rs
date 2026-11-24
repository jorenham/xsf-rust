use num_complex::c64;

xsref::tests! {
    test_riemann_zeta_f64: "riemann_zeta", "d-d", |x| xsf::riemann_zeta(x[0]);
    test_riemann_zeta_c64: "riemann_zeta", "cd-cd", |x| xsf::riemann_zeta(c64(x[0], x[1]));
    test_zeta_f64: "zeta", "d_d-d", |x| xsf::zeta(x[0], x[1]);
    test_zetac_f64: "zetac", "d-d", |x| xsf::zetac(x[0]);
}

#[test]
fn test_riemann_zeta_c64_extreme() {
    // (k + 1)^z overflows in the Borwein sum, and 1 / inf = 0 (C99 Annex G)
    let w = xsf::riemann_zeta(c64(1000.0, 0.5));
    assert_eq!((w.re, w.im), (1.0, -3.170_083_028_960_87e-302));
    // the reflection factor overflows, and zeta(1 - z) is finite
    assert_eq!(xsf::riemann_zeta(c64(-1000.0, 0.5)).re, f64::INFINITY);
    // too many Euler-Maclaurin terms (undefined behavior in xsf)
    for im in [7.4e19, 1e300, f64::INFINITY] {
        let w = xsf::riemann_zeta(c64(3.0, im));
        assert!(w.re.is_nan() && w.im.is_nan());
    }
}

#[test]
fn test_riemann_zeta_c64_moderate_imag() {
    // An eta denominator zero and heights where the Borwein sum is too short (Arb references).
    for (z, expected) in [
        (c64(1.0, 27.194_160_850_963_16), c64(1.901_727_442_929_379_8, 0.482_487_690_810_793_9)),
        (c64(0.5, -40.0), c64(0.793_044_952_561_928_6, 1.041_274_614_651_065)),
        (c64(-0.5, 45.0), c64(7.121_401_947_237_887, 9.435_007_642_389_886)),
        (c64(3.0, 47.0), c64(1.045_054_169_248_156, -0.166_468_454_920_766_37)),
    ] {
        let rel = (xsf::riemann_zeta(z) - expected).norm() / expected.norm();
        assert!(rel < 1e-13, "riemann_zeta({z}): rel = {rel:e}");
    }
    // |zeta| is small here, so the condition number |z zeta'(z) / zeta(z)| is ~1e4
    let expected = c64(-0.006_058_168_095_158_761, -0.005_369_758_200_077_69);
    let rel = (xsf::riemann_zeta(c64(0.5, 48.0)) - expected).norm() / expected.norm();
    assert!(rel < 1e-11, "rel = {rel:e}");
}

#[test]
fn test_riemann_zeta_c64_borwein_height() {
    let s = c64(0.427_718_477_200_037_64, 21.783_747_851_715_468);
    // Reference from Arb at the exact binary64 input.
    let expected = c64(0.476_703_263_844_432_9, 0.689_814_874_460_312_3);
    for (z, f) in [(s, expected), (s.conj(), expected.conj())] {
        let rel = (xsf::riemann_zeta(z) - f).norm() / f.norm();
        assert!(rel < 5e-15, "{z}: relative error {rel:e}");
    }
}

#[test]
fn test_riemann_zeta_c64_near_eta_factor_zeros() {
    // 1 - 2^(1 - z) vanishes at 1 + 2 pi i k / ln(2) (reference values from Arb)
    for (z, expected) in [
        (c64(1.000_001, 1e-6), c64(500_000.577_215_737_7, -500_000.000_041_060_5)),
        (c64(1.5, 1e-100), c64(2.612_375_348_685_488, -3.932_239_737_431_101_5e-100)),
        (c64(1.0, 1e-307), c64(0.577_215_664_901_532_9, -1.000_000_000_000_000_1e307)),
        (c64(0.0, 1e-307), c64(-0.5, -9.189_385_332_046_727e-308)),
        (c64(1e-300, 1e-300), c64(-0.5, -9.189_385_332_046_728e-301)),
        (
            c64(-1.606_574_395_346_453_3e-16, 2.407_500_379_740_984e-254),
            c64(-0.499_999_999_999_999_83, -2.212_344_867_648_871_6e-254),
        ),
        (c64(1.0, 9.064_720_283_654_388), c64(1.346_579_542_836_317_2, 0.109_883_136_796_269_5)),
        (c64(0.0, -9.064_720_283_654_388), c64(1.598_734_526_808_69, -0.278_338_669_639_082_5)),
        (
            c64(1.000_000_000_000_000_2, 18.129_440_567_308_777),
            c64(1.824_744_877_856_936_2, -0.216_035_216_066_775_78),
        ),
    ] {
        // each part separately, as the norm hides e.g. Re zeta(1 + 1e-307i)
        let w = xsf::riemann_zeta(z);
        for (got, want) in [(w.re, expected.re), (w.im, expected.im)] {
            let rel = (got - want).abs() / want.abs();
            assert!(rel < 2e-14, "riemann_zeta({z}) = {w}: rel = {rel:e}");
        }
    }
}

#[test]
fn test_riemann_zeta_c64_small() {
    // the reflection formula evaluates zeta(1 - z) next to its pole (reference values from Arb)
    for (z, expected) in [
        (c64(1e-6, -2e-6), c64(-0.500_000_918_935_523_7, 1.837_881_079_120_255_7e-6)),
        (c64(1e-5, 1e-5), c64(-0.500_009_189_385_33, -9.189_585_969_693_889e-6)),
        (c64(-1e-3, 3e-3), c64(-0.499_073_062_089_067_7, -0.002_750_778_612_088_271)),
        (c64(0.1, -0.2), c64(-0.550_613_484_192_186_5, 0.219_207_672_013_450_7)),
        (c64(-0.3, 1e-20), c64(-0.293_813_068_129_721_24, -5.089_726_696_918_256_5e-21)),
    ] {
        let w = xsf::riemann_zeta(z);
        for (got, want) in [(w.re, expected.re), (w.im, expected.im)] {
            let rel = (got - want).abs() / want.abs();
            assert!(rel < 1e-15, "riemann_zeta({z}) = {w}: rel = {rel:e}");
        }
    }
}

#[test]
fn test_zeta_large_q() {
    // asymptotic expansion for q > 1e8 (reference value from mpmath)
    let (x, q) = (36.0, 1.1e8);
    let expected = 1.016_688_811_413_570_4e-283;
    let rel = (xsf::zeta(x, q) - expected).abs() / expected.abs();
    assert!(rel < 1e-15, "rel = {rel:e}");

    assert_eq!(xsf::zeta(f64::INFINITY, q), 0.0);
    assert_eq!(xsf::zeta(f64::INFINITY, f64::INFINITY), 0.0);
}

#[test]
fn test_zeta_near_1_tiny_q() {
    // the tail sum ~ 1 / (x - 1) isn't negligible (reference values from mpmath)
    for (x, q, expected) in [
        (1.0 + f64::EPSILON, 1e-15, 5.503_599_627_370_504e15),
        (1.000_000_000_001, 1e-17, 1.000_009_999_150_220_6e17),
    ] {
        let rel = (xsf::zeta(x, q) - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "zeta({x}, {q}): rel = {rel:e}");
    }
}

#[test]
fn test_zeta_negative_q() {
    // reference values from mpmath, by summing the terms with q + k < 0 directly
    for (x, q, expected) in [
        (4.0, -431.636_697_197_915_6, 63.992_140_081_094_475),
        (2.0, -1e10 - 0.5, 9.869_604_400_989_358),
        (31.0, -7.5, 1.593_135_539_215_446_3e-29),
        (31.0, -0.500_000_000_000_000_1, 3.304_063_430_367_368e-5),
    ] {
        let rel = (xsf::zeta(x, q) - expected).abs() / expected.abs();
        assert!(rel < 1e-15, "zeta({x}, {q}): rel = {rel:e}");
    }
    // all terms underflow
    assert_eq!(xsf::zeta(1e300, 1.5), 0.0);
    assert_eq!(xsf::zeta(f64::INFINITY, 1.5), 0.0);
    assert_eq!(xsf::zeta(f64::INFINITY, -0.5), f64::INFINITY);
}

#[test]
fn test_zetac_large_x() {
    // zetac(x) = 2^-x to double precision, which underflows to 0 for x >= 1075
    for x in [127.0, 130.0, 1000.0, 1074.0, 1076.0, 1e300] {
        assert_eq!(xsf::zetac(x), (-x).exp2());
    }
}

#[test]
fn test_riemann_zeta_c64_reflection_overflow() {
    let re = -294.866_054_613_497_45;
    // Arb encloses both components beyond the finite binary64 range.
    for (im, expected) in [
        (-16.147_667_799_398_363, c64(f64::INFINITY, f64::NEG_INFINITY)),
        (13.992_648_136_816_397, c64(f64::NEG_INFINITY, f64::INFINITY)),
    ] {
        assert_eq!(xsf::riemann_zeta(c64(re, im)), expected);
        assert_eq!(xsf::riemann_zeta(c64(re, -im)), expected.conj());
    }
}
