use num_traits::ToPrimitive;

// The xsref tables also contain negative `n`, which in scipy returns NaN.

#[test]
fn test_expn_u32() {
    xsref::test("expn", "p_d-d", |x| {
        if x[0] < 0.0 {
            f64::NAN
        } else {
            xsf::expn(x[0].to_u32().unwrap(), x[1])
        }
    });
}

#[test]
fn test_expn_f64() {
    // not sure why this table exists; but we might as well use it 🤷🏻
    xsref::test("expn", "d_d-d", |x| {
        if x[0] < 0.0 {
            f64::NAN
        } else {
            xsf::expn(x[0].to_u32().unwrap(), x[1])
        }
    });
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
