use num_complex::c64;

xsref::tests! {
    test_exp1_f64: "exp1", "d-d", |x| xsf::exp1(x[0]);
    test_exp1_c64: "exp1", "cd-cd", |x| xsf::exp1(c64(x[0], x[1]));
    test_expi_f64: "expi", "d-d", |x| xsf::expi(x[0]);
    test_expi_c64: "expi", "cd-cd", |x| xsf::expi(c64(x[0], x[1]));
    test_scaled_exp1_f64: "scaled_exp1", "d-d", |x| xsf::scaled_exp1(x[0]);
}
