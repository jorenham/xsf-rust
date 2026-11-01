#[test]
fn test_cbrt() {
    xsref::test("cbrt", "d-d", |x| xsf::cbrt(x[0]));
}
