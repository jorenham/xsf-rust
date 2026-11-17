xsref::tests! {
    test_binom_f64: "binom", "d_d-d", |x| xsf::binom(x[0], x[1]);
}

#[test]
fn test_binom_large_k() {
    // scipy/xsf#352: the n (n + 1) / (2 k) term, the parity of k >= 2^31, sinpi, and the split
    // |k|^n (which over- or underflows here)
    for (n, k, expected) in [
        (2.5, 300_000_000.7, 1.329_598_785_965_806_2e-30),
        (0.5, 10_000_000_001.25, 1.994_711_401_707_956_7e-16),
        (20.5, 2_100_000_000.25, -9.316_864_069_123_419e-183),
        (-10.5, 500_000_000_000.25, 8.617_145_704_224_517e104),
        (31.687_271_077_096_36, 10_283_293_173.796_093, 5.092_788_374_279_372e-294),
        (-26.5, 1e12, 1.270_443_843_261_605_3e280),
    ] {
        let actual = xsf::binom(n, k);
        assert!((actual / expected - 1.0).abs() < 5e-14, "{actual:e} != {expected:e}");
    }
    assert_eq!(xsf::binom(1.0, 2_147_483_647.0), 0.0);
    assert_eq!(xsf::binom(169.0, 20_000_000_000.0), 0.0);
}

#[test]
fn test_binom_large_n() {
    // unlike in xsf, `n >= 1e10 k` without `exp(-lbeta)`, which errs by up to ~1e3 eps here (Arb)
    let nk = [
        (2_467_198_972_239.182, 19.424_570_770_297_443),
        (2_145_255_138_526.0, 23.0),
        (6.853_882_129_587_282e248, 0.836_181_315_884_165_7),
        (27_132.662_274_677_346, 2.579_539_995_968_121_7e-10),
        (1e13, 0.059),
    ];
    let actual = nk.map(|(n, k)| xsf::binom(n, k));
    let expected = [
        1.197_922_460_130_211_1e223,
        1.627_501_668_305_473e261,
        1.253_563_604_787_780_2e208,
        1.000_000_002_782_221_6,
        6.033_665_105_512_786,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
    // Arb: 1.797_693_134_862_315_99e308, which rounds to inf
    assert_eq!(xsf::binom(1.776_229_249_145_602e280, 1.1), f64::INFINITY);
}

#[test]
fn test_binom_reflection() {
    // unlike in xsf (scipy/xsf#352): the reflection formula for |k| > 2e4 (the last case splits
    // `poch`)
    for (n, k, expected) in [
        (31.69, 1_000_000_000.796, -5.746_138_060_441_207e-261),
        (-31.750_000_000_000_004, 20_000.25, -8.409_833_121_198_072e84),
        (-0.000_21, 20_001.0, -1.052_260_894_432_576_9e-8),
        (-0.000_21, 20_000.999_9, -5.511_843_042_336_944e-9),
        (2.5, -1_000_000.3, 8.558_184_094_049_863e-22),
        (-7.25, -1_000_000_000_000.5, -1.224_023_543_149_038e72),
        (-38.699_131_748_729_36, 2_169_303_697.740_144_3, 6.492_179_234_150_414e307),
    ] {
        let actual = xsf::binom(n, k);
        assert!((actual / expected - 1.0).abs() < 5e-14, "{actual:e} != {expected:e}");
    }
    assert_eq!(xsf::binom(2.5, f64::NEG_INFINITY), 0.0);
}

#[test]
fn test_binom_negative_n() {
    // unlike in xsf (scipy/xsf#352), without the rounding error of 1 + n - k (mpmath, Wolfram)
    // (the last one is near the switch to `lgam` in `beta`)
    let nk = [
        (-75_088_638.198_055_71, 0.805_624_372_712_259_8),
        (-21_271_018.999_750_603, 0.001_101_100_669_294_745),
        (-171.999_999_999_999_97, 0.3),
    ];
    let actual = nk.map(|(n, k)| xsf::binom(n, k));
    let expected = [-47_030.442_475_758_64, -3.481_296_578_119_637, -47_264_014_799_068.09];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
    // B(1 + k, -n) under- or overflows, and n - k or k is an integer in the first three
    assert_eq!(xsf::binom(-10_000.5, 200.5), 0.0);
    assert_eq!(xsf::binom(-173.5, -2.0), 0.0);
    assert_eq!(xsf::binom(-1_000.5, -2.0), 0.0);
    let actual = [xsf::binom(-9_999.999_999_999_998, -135.3)];
    xsref::np_assert_allclose!(&actual, &[5.038_684_715_231_528e-302], rtol = 5e-13);
}

#[test]
fn test_binom_integer_zero() {
    // unlike in xsf, 1 / Gamma(1 + n - k) = 0 also for |k| >= 2^53; the binom table covers
    // 1 / Gamma(1 + k) = 0
    assert_eq!(xsf::binom(171.0, 1e16), 0.0);
}

#[test]
fn test_binom_reflection_negative_a() {
    // unlike in xsf (scipy/xsf#352), without the rounding error of a = 1 + n - k < 0, also near
    // (or onto) a pole of Gamma(a) (Arb)
    let nk = [(-0.333_333_333_333_333_3, 6.666_666_666_666_667), (300.5, 700.75)];
    let actual = nk.map(|(n, k)| xsf::binom(n, k));
    let expected = [-1.043_378_523_381_252_6e-16, 2.591_771_256_021_541_7e-210];
    xsref::np_assert_allclose!(&actual[..1], &expected[..1], rtol = 5e-15);
    xsref::np_assert_allclose!(&actual[1..], &expected[1..], rtol = 2e-13);
}
