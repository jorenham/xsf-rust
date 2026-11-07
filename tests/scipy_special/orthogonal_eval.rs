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
    np_assert_allclose!(actual, desired, rtol = 2e-12);

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

    // or near x = -1, where `hyp2f1` was inaccurate (mpmath)
    let actual = [
        xsf::eval_jacobi(15.8, -0.5, 8.7, -0.96),
        xsf::eval_gegenbauer(40.7, 9.3, -0.99),
    ];
    let desired = [37_787_979.694_627_58, -1.002_223_611_030_761_3e15];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // and where `hyp2f1` was inaccurate for c - a - b <= -1 (mpmath)
    let actual = [
        xsf::eval_jacobi(26.39, 1.9, 3.0, -0.99),
        xsf::eval_jacobi(32.06, 3.2, 3.0, -0.88),
        xsf::eval_jacobi(17.32, 0.1, 3.0, -0.99),
    ];
    let desired = [
        64.703_859_369_457_89,
        11.922_922_622_177_607,
        1_132.803_574_409_967_7,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_eval_genlaguerre_alpha_near_minus_one() {
    // `binom(n + alpha, n)` lost the digits of `alpha + 1` (mpmath)
    let actual = [
        xsf::eval_genlaguerre(30_i32, -1.0 + 1e-14, 10.0),
        xsf::eval_genlaguerre(20_i32, -1.0 + 1e-12, 3.0),
        xsf::eval_genlaguerre(7_i32, -1.0 + 1e-15, 2.5),
        xsf::eval_genlaguerre(25.0, -1.0 + 1e-13, 4.0),
        xsf::eval_genlaguerre(40.0, -1.0 + 1e-6, 25.0),
    ];
    let desired = [
        -10.969_942_483_703_025,
        -0.312_622_172_273_621_26,
        -0.554_044_208_829_365,
        -0.131_747_458_557_938_53,
        19_671.412_951_251_726,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // but not for n > 2^53, where `n - 1` rounds (mpmath: 0.999999999999998)
    let actual = xsf::eval_genlaguerre(2.0_f64.powi(54), -2.0_f64.powi(-54), 0.0);
    np_assert_allclose!([actual], [0.999_999_999_999_997_9], rtol = 1e-14);
}
