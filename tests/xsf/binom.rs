#[test]
fn test_binom_f64() {
    xsref::test("binom", "d_d-d", |x| xsf::binom(x[0], x[1]));
}

#[test]
fn test_binom_large_k() {
    // unlike in xsf (scipy/xsf#352): the n (n + 1) / (2 k) term, the parity of k >= 2^31,
    // sinpi, and the split |k|^n (which over- or underflows here)
    for (n, k, expected) in [
        (2.5, 300_000_000.7, 1.329_598_785_965_806_2e-30),
        (0.5, 10_000_000_001.25, 1.994_711_401_707_956_7e-16),
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
