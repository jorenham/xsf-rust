use num_complex::c64;

#[test]
fn test_sici_f64() {
    xsref::test("sici", "d-d_d", |x| xsf::sici(x[0]));
}

#[test]
fn test_sici_c64() {
    xsref::test("sici", "cd-cd_cd", |x| xsf::sici(c64(x[0], x[1])));
}

#[test]
fn test_shichi_f64() {
    xsref::test("shichi", "d-d_d", |x| xsf::shichi(x[0]));
}

#[test]
fn test_shichi_c64() {
    xsref::test("shichi", "cd-cd_cd", |x| xsf::shichi(c64(x[0], x[1])));
}
