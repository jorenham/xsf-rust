#[test]
fn test_lanczos_sum_expg_scaled() {
    xsref::test("lanczos_sum_expg_scaled", "d-d", |x| {
        xsf::lanczos_sum_expg_scaled(x[0])
    });
}
