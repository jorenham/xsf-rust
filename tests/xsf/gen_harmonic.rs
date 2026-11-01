use xsf::{gen_harmonic, normalized_gen_harmonic};

#[test]
fn test_gen_harmonic() {
    // mpmath: sum(mpf(i)**-a for i in range(1, n + 1))
    let cases: [(u64, f64, f64); 8] = [
        (2, 2.0, 1.25),
        (10, 0.5, 5.020_997_899_292_666_5),
        (100, -1.5, 40_501.224_515_318_936),
        (1000, 1.0, 7.485_470_860_550_345),
        (1000, 2.0, 1.643_934_566_681_56),
        (10, 1.01, 2.902_261_421_419_068_7),
        (1_000_000_000_000_000, 3.0, 1.202_056_903_159_594_2),
        (7, -0.25, 9.606_504_613_320_602),
    ];
    xsref::np_assert_allclose!(
        cases.map(|(n, a, _)| gen_harmonic(n, a)),
        cases.map(|(.., h)| h),
        rtol = 1e-15
    );

    assert!(gen_harmonic(0, 2.0).is_nan());
    assert_eq!(gen_harmonic(1, f64::NAN), 1.0);
    assert!(gen_harmonic(2, f64::NAN).is_nan());
    assert_eq!(gen_harmonic(5, f64::INFINITY), 1.0);
    assert_eq!(gen_harmonic(5, f64::NEG_INFINITY), f64::INFINITY);
    assert_eq!(gen_harmonic(5, 0.0), 5.0);
}

#[test]
fn test_normalized_gen_harmonic() {
    // mpmath: H(j..k, a) / H(1..n, a)
    let cases: [(u64, u64, u64, f64, f64); 7] = [
        (1, 1, 2, 2.0, 0.8),
        (2, 5, 10, 0.5, 0.444_467_552_195_255_47),
        (3, 50, 100, -1.5, 0.178_880_381_185_630_27),
        (10, 100, 1000, 1.0, 0.315_064_918_106_964_85),
        (1, 10, 1000, 2.0, 0.942_718_623_098_786_7),
        (5, 8, 10, 1.01, 0.214_675_893_759_013_25),
        (10, 20, 1_000_000_000_000, 3.0, 0.003_607_030_809_312_731),
    ];
    xsref::np_assert_allclose!(
        cases.map(|(j, k, n, a, _)| normalized_gen_harmonic(j, k, n, a)),
        cases.map(|(.., r)| r),
        rtol = 1e-15
    );

    assert!(normalized_gen_harmonic(0, 1, 1, 2.0).is_nan());
    assert!(normalized_gen_harmonic(2, 1, 3, 2.0).is_nan());
    assert!(normalized_gen_harmonic(1, 3, 2, 2.0).is_nan());
    assert_eq!(normalized_gen_harmonic(1, 1, 1, f64::NAN), 1.0);
    assert_eq!(normalized_gen_harmonic(1, 2, 3, f64::INFINITY), 1.0);
    assert_eq!(normalized_gen_harmonic(2, 2, 3, f64::INFINITY), 0.0);
    assert_eq!(normalized_gen_harmonic(1, 1, 3, f64::NEG_INFINITY), 0.0);
    assert!(normalized_gen_harmonic(1, 2, 3, f64::NEG_INFINITY).is_nan());
    assert_eq!(normalized_gen_harmonic(2, 3, 4, 0.0), 0.5);
}
