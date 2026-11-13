use num_traits::ToPrimitive;

xsref::tests! {
    // The xsref tables also contain negative `n`, which in scipy returns NaN.
    test_expn_u32: "expn", "p_d-d", |x| {
        if x[0] < 0.0 {
            f64::NAN
        } else {
            xsf::expn(x[0].to_u32().unwrap(), x[1])
        }
    };

    // not sure why this table exists; but we might as well use it 🤷🏻
    test_expn_f64: "expn", "d_d-d", |x| {
        if x[0] < 0.0 {
            f64::NAN
        } else {
            xsf::expn(x[0].to_u32().unwrap(), x[1])
        }
    };
}
#[test]
fn test_expn_large_n() {
    // n > i32::MAX, which doesn't fit in the C++ `int n`
    for (n, x, expected) in [
        (2_147_483_648, 1.0, 1.713_072_141_499_455_7e-10),
        (u32::MAX, 1.0, 8.565_360_709_491_557e-11),
    ] {
        let actual = xsf::expn(n, x);
        assert!(
            (actual / expected - 1.0).abs() < 1e-14,
            "{actual:e} != {expected:e}"
        );
    }
    assert_eq!(xsf::expn(u32::MAX, 0.0), 1.0 / f64::from(u32::MAX - 1));
}

#[test]
fn test_expn_subnormal() {
    // MAXLOG < x, where E_n(x) is subnormal (reference values from mpmath and Arb, rounded)
    let x = 714.942_898_659_757_7;
    for (n, expected) in [
        (0, 4.466_556_987e-314),
        (1, 4.460_326_959e-314),
        (50, 4.174_959_247e-314),
        (51, 4.169_514_667_5e-314),
    ] {
        let actual = xsf::expn(n, x);
        assert!(
            (actual - expected).abs() <= 5e-324,
            "expn({n}, {x}) = {actual:e} != {expected:e}"
        );
    }
    // just above half the smallest subnormal, in each of the three paths
    for (n, x) in [
        (0, 738.528_114_325_140_7),
        (10, 738.514_767_187_870_7),
        (52, 738.46),
    ] {
        assert_eq!(xsf::expn(n, x), f64::from_bits(1), "expn({n}, {x})");
    }
    // just below it (by 2.5e-14), where -lambda p (rounded) instead of -x gave 5e-324
    assert_eq!(xsf::expn(69, 738.439_456_627_518_1), 0.0);
}
