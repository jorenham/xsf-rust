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
}
