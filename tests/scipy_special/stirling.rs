const STIRLING2_TABLE: [[u64; 11]; 11] = [
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 3, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 7, 6, 1, 0, 0, 0, 0, 0, 0],
    [0, 1, 15, 25, 10, 1, 0, 0, 0, 0, 0],
    [0, 1, 31, 90, 65, 15, 1, 0, 0, 0, 0],
    [0, 1, 63, 301, 350, 140, 21, 1, 0, 0, 0],
    [0, 1, 127, 966, 1_701, 1_050, 266, 28, 1, 0, 0],
    [0, 1, 255, 3_025, 7_770, 6_951, 2_646, 462, 36, 1, 0],
    [0, 1, 511, 9_330, 34_105, 42_525, 22_827, 5_880, 750, 45, 1],
];

#[test]
fn test_stirling2() {
    for (n_idx, row) in STIRLING2_TABLE.iter().enumerate() {
        let n = n_idx.try_into().unwrap();
        for (k_idx, &s0) in row.iter().enumerate() {
            let k = k_idx.try_into().unwrap();
            let s: u64 = xsf::stirling2(n, k);
            assert_eq!(s, s0, "S({n_idx}, {k_idx})");
        }
    }
}

#[test]
fn test_stirling2_k2_large_n() {
    assert_eq!(xsf::stirling2::<u64>(33, 2), (1 << 32) - 1);
    assert_eq!(xsf::stirling2::<u128>(129, 2), u128::MAX);
}

#[test]
fn test_stirling2_f64_small() {
    // `S(n, k)` fits in `u128` for `n <= 43`
    for n in 0..=43 {
        for k in 0..=n + 1 {
            let s0: u128 = xsf::stirling2(n, k);
            #[allow(clippy::cast_precision_loss)]
            let s0 = s0 as f64;
            let s: f64 = xsf::stirling2(n, k);
            assert!((s - s0).abs() <= 4.0 * f64::EPSILON * s0, "S({n}, {k})");
        }
    }
    // exact values from `mpmath.stirling2`
    for (k, s0) in [(10, 2.615_471_651_586_288e43), (25, 7.453_802_153_273_2e42)] {
        let s: f64 = xsf::stirling2(50, k);
        assert!((s - s0).abs() <= 4.0 * f64::EPSILON * s0, "S(50, {k})");
    }
}

#[test]
fn test_stirling2_f64_dp_large_n() {
    // `50 < n <= 219`, exact values
    for (n, k, s0) in [
        (51, 7, 2.491_134_287_812_361_2e39),
        (52, 7, 1.744_467_199_136_095_8e40),
        (150, 12, 1.573_398_590_522_490_7e153),
        (219, 54, 4.190_119_213_575_601e307),
        (219, 110, 2.071_041_889_974_642e260),
    ] {
        let s: f64 = xsf::stirling2(n, k);
        assert!((s - s0).abs() <= 2e-15 * s0, "S({n}, {k}) = {s} != {s0}");
    }
}

#[test]
fn test_stirling2_f64_temme() {
    assert_eq!(xsf::stirling2::<f64>(60, 0), 0.0);
    assert_eq!(xsf::stirling2::<f64>(60, 61), 0.0);
    assert_eq!(xsf::stirling2::<f64>(100, 1), 1.0);
    assert_eq!(xsf::stirling2::<f64>(100, 100), 1.0);
    // exact values from `mpmath.stirling2`
    for (n, k, s0) in [
        (51, 2, 1_125_899_906_842_623.0),
        (51, 25, 2.297_427_976_320_779e44),
        (100, 3, 8.589_625_345_533_522e46),
        (100, 50, 4.309_832_370_093_663_6e101),
        (100, 99, 4_950.0),
        (1000, 990, 2.361_810_141_115_116_3e50),
    ] {
        let s: f64 = xsf::stirling2(n, k);
        assert!((s - s0).abs() <= 2e-5 * s0, "S({n}, {k}) = {s} != {s0}");
    }
}
