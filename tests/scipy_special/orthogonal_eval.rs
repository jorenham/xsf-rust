use xsref::np_assert_allclose;

#[test]
fn test_eval_real_degree_large() {
    // the hypergeometric series cancels catastrophically for large non-integer degrees (mpmath)
    let actual = [
        xsf::eval_legendre(100.5, 0.3),
        xsf::eval_chebyshev_t(-90.7, -0.88),
        xsf::eval_gegenbauer(40.5, 2.5, -0.3),
        xsf::eval_jacobi(49.5, -0.16, -0.98, -0.97),
    ];
    let desired = [
        0.011_650_211_398_786_529,
        0.276_268_747_018_169_5,
        21.990_470_696_631_68,
        -0.010_204_813_296_612_693,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // but not if the starting values change `c - a - b = 0`, or the recurrence is unstable (mpmath)
    let actual = [
        xsf::eval_legendre(7.344_542_190_030_112, -0.996_500_368_048_799_2),
        xsf::eval_jacobi(10.3, 0.0, 30.0, -0.9),
        xsf::eval_jacobi(10.000_001, 0.0, 30.0, -0.5),
        xsf::eval_jacobi(10.000_000_01, -0.99, 40.0, -0.3),
    ];
    let desired = [
        -0.208_532_505_668_321_68,
        -1.648_343_464_373_981_2e28,
        13_912_033.739_818_575,
        34_691_102.898_867_61,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-11);
}
