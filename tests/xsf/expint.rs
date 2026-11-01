use num_complex::c64;

#[test]
fn test_exp1_f64() {
    xsref::test("exp1", "d-d", |x| xsf::exp1(x[0]));
}

#[test]
fn test_exp1_c64() {
    xsref::test("exp1", "cd-cd", |x| xsf::exp1(c64(x[0], x[1])));
}

#[test]
fn test_expi_f64() {
    xsref::test("expi", "d-d", |x| xsf::expi(x[0]));
}

#[test]
fn test_expi_c64() {
    xsref::test("expi", "cd-cd", |x| xsf::expi(c64(x[0], x[1])));
}

#[test]
fn test_scaled_exp1_f64() {
    xsref::test("scaled_exp1", "d-d", |x| xsf::scaled_exp1(x[0]));
}
