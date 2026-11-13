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
        (
            31.687_271_077_096_36,
            10_283_293_173.796_093,
            5.092_788_374_279_372e-294,
        ),
        (-26.5, 1e12, 1.270_443_843_261_605_3e280),
    ] {
        let actual = xsf::binom(n, k);
        assert!(
            (actual / expected - 1.0).abs() < 5e-14,
            "{actual:e} != {expected:e}"
        );
    }
    assert_eq!(xsf::binom(1.0, 2_147_483_647.0), 0.0);
    assert_eq!(xsf::binom(169.0, 20_000_000_000.0), 0.0);
}

#[test]
fn test_binom_reflection() {
    // unlike in xsf (scipy/xsf#352): the reflection formula for |k| > 2e4 (the last case splits
    // `poch`)
    for (n, k, expected) in [
        (31.69, 1_000_000_000.796, -5.746_138_060_441_207e-261),
        (
            -31.750_000_000_000_004,
            20_000.25,
            -8.409_833_121_198_072e84,
        ),
        (-0.000_21, 20_001.0, -1.052_260_894_432_576_9e-8),
        (-0.000_21, 20_000.999_9, -5.511_843_042_336_944e-9),
        (2.5, -1_000_000.3, 8.558_184_094_049_863e-22),
        (-7.25, -1_000_000_000_000.5, -1.224_023_543_149_038e72),
        (
            -38.699_131_748_729_36,
            2_169_303_697.740_144_3,
            6.492_179_234_150_414e307,
        ),
    ] {
        let actual = xsf::binom(n, k);
        assert!(
            (actual / expected - 1.0).abs() < 5e-14,
            "{actual:e} != {expected:e}"
        );
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
    let expected = [
        -47_030.442_475_758_64,
        -3.481_296_578_119_637,
        -47_264_014_799_068.09,
    ];
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
    // B(1 + k, -n) under- or overflows, and n - k or k is an integer in the first three
    assert_eq!(xsf::binom(-10_000.5, 200.5), 0.0);
    assert_eq!(xsf::binom(-173.5, -2.0), 0.0);
    assert_eq!(xsf::binom(-1_000.5, -2.0), 0.0);
    let actual = [xsf::binom(-9_999.999_999_999_998, -135.3)];
    xsref::np_assert_allclose!(&actual, &[5.038_684_715_231_528e-302], rtol = 5e-13);
}
