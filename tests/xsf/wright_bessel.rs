#[test]
fn test_wright_bessel_f64() {
    xsref::test("wright_bessel", "d_d_d-d", |x| {
        xsf::wright_bessel(x[0], x[1], x[2])
    });
}

#[test]
fn test_log_wright_bessel_f64() {
    xsref::test("log_wright_bessel", "d_d_d-d", |x| {
        xsf::log_wright_bessel(x[0], x[1], x[2])
    });
}
