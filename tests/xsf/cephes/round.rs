#[test]
fn test_round() {
    xsref::test("round", "d-d", |x| xsf::round(x[0]));
}
