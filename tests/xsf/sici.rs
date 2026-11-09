use num_complex::c64;

xsref::tests! {
    test_sici_f64: "sici", "d-d_d", |x| xsf::sici(x[0]);
    test_sici_c64: "sici", "cd-cd_cd", |x| xsf::sici(c64(x[0], x[1]));
    test_shichi_f64: "shichi", "d-d_d", |x| xsf::shichi(x[0]);
    test_shichi_c64: "shichi", "cd-cd_cd", |x| xsf::shichi(c64(x[0], x[1]));
}
