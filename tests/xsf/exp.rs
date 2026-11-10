use num_complex::c64;

#[test]
fn test_expm1_f64() {
    xsref::test("expm1", "d-d", |x| xsf::expm1(x[0]));
}

#[test]
fn test_expm1_c64() {
    xsref::test("expm1", "cd-cd", |x| xsf::expm1(c64(x[0], x[1])));
}

#[test]
fn test_exp2_f64() {
    xsref::test("exp2", "d-d", |x| xsf::exp2(x[0]));
}

#[test]
fn test_exp10_f64() {
    xsref::test("exp10", "d-d", |x| xsf::exp10(x[0]));
}
