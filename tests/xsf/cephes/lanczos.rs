xsref::tests! {
    test_lanczos_sum_expg_scaled: "lanczos_sum_expg_scaled", "d-d", |x| {
        xsf::lanczos_sum_expg_scaled(x[0])
    };
}
