#[test]
fn test_perm_pos() {
    assert_eq!(xsf::perm(10, 0), 1);
    assert_eq!(xsf::perm(10, 1), 10);
    assert_eq!(xsf::perm(10, 2), 90);
    assert_eq!(xsf::perm(10, 3), 720);
    assert_eq!(xsf::perm(10, 4), 5_040);
    assert_eq!(xsf::perm(10, 5), 30_240);
}

#[test]
fn test_perm_neg() {
    assert_eq!(xsf::perm(-10, 0), 1);
    assert_eq!(xsf::perm(-10, 1), -10);
    assert_eq!(xsf::perm(-10, 2), 110);
    assert_eq!(xsf::perm(-10, 3), -1_320);
    assert_eq!(xsf::perm(-10, 4), 17_160);
    assert_eq!(xsf::perm(-10, 5), -240_240);
    // no `n - k` beyond the last factor
    assert_eq!(xsf::perm(i8::MIN, 0), 1);
    assert_eq!(xsf::perm(i8::MIN + 1, 1), -127);
}

#[test]
fn test_perm_k_above_n() {
    assert_eq!(xsf::perm(3_u32, 4), 0);
    assert_eq!(xsf::perm(3_u32, 5), 0);
    assert_eq!(xsf::perm(0_u8, 1), 0);
    assert_eq!(xsf::perm(3, 5), 0);
    assert_eq!(xsf::perm(2.0_f64, 4).to_bits(), 0.0_f64.to_bits());
    assert_eq!(xsf::perm(2.5, 4), -0.9375);
}
