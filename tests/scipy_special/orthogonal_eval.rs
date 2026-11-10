use num_complex::c64;
use num_traits::ToPrimitive;
use xsf::{
    eval_chebyshev_t, eval_chebyshev_u, eval_gegenbauer, eval_genlaguerre, eval_hermite_h,
    eval_hermite_he, eval_jacobi, eval_laguerre, eval_legendre,
};
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
fn test_eval_jacobi_alpha_near_minus_one() {
    // `binom(n + alpha, n)` lost the digits of `alpha + 1`, and `hyp2f1` treated it as 0 (Arb)
    let actual = [
        xsf::eval_jacobi(7_i32, -1.0 + 1e-15, 2.0, -0.6),
        xsf::eval_jacobi(7.5, -1.0 + 1e-15, 2.0, -0.6),
        xsf::eval_jacobi(7.5, -1.0 + 1e-15, 2.0, 1.0),
        xsf::eval_jacobi(46.5, -1.0 + 1e-14, -0.4, -0.57),
    ];
    let desired = [
        1.370_879_999_999_999_2,
        1.337_487_196_225_991_4,
        1.332_267_629_550_191_2e-16,
        0.040_543_716_001_115_56,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    let actual = xsf::eval_jacobi(7.5, -1.0 + 1e-15, 2.0, c64(-0.6, 0.25));
    let desired = c64(6.901_490_182_319_64, -3.299_333_682_536_697);
    np_assert_allclose!([actual], [desired], rtol = 1e-12);

    // and at alpha = -1 with an integer degree, for any beta (Arb)
    let actual = [
        xsf::eval_jacobi(5_i32, -1.0, -0.9, -0.99),
        xsf::eval_jacobi(30_i32, -1.0, 0.5, 0.3),
        xsf::eval_jacobi(64_i32, -1.0, 9.5, -0.75),
        xsf::eval_sh_jacobi(10_i32, -1.0, 2.5, 0.3),
    ];
    let desired = [
        -5.342_672_184_935_159e-4,
        0.017_321_739_003_478_2,
        -2_235.655_409_379_344_7,
        -0.655_967_500_183_011_9,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // near alpha = -8 the recurrence cancels, but the result stays small (Arb: 3.7e-9)
    assert!(xsf::eval_jacobi(20_i32, -7.999_999_999_999_998, 1.0, 0.9585).abs() < 1e-8);
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
        xsf::eval_genlaguerre(2_i32, -0.999_999_999_999_99, 0.5),
        xsf::eval_genlaguerre(3.0, -0.999_999_999_999_999_8, 2.0),
    ];
    let desired = [
        -10.969_942_483_703_025,
        -0.312_622_172_273_621_26,
        -0.554_044_208_829_365,
        -0.131_747_458_557_938_53,
        19_671.412_951_251_726,
        -0.375,
        0.666_666_666_666_666_5,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // alpha >= -0.5 keeps `binom(n + alpha, n)`, also for n > 2^53 (mpmath: 0.999999999999998)
    let actual = xsf::eval_genlaguerre(2.0_f64.powi(54), -2.0_f64.powi(-54), 0.0);
    np_assert_allclose!([actual], [0.999_999_999_999_997_9], rtol = 1e-14);
}

#[test]
fn test_eval_jacobi_eq_legendre() {
    // manual test cases
    let xs = [-5.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 5.0];
    let ns = [0, 1, 2, 3, 4, 8, 15];
    for n in ns {
        let p0 = xs.map(|x| xsf::legendre_p(n, x));

        let pi = xs.map(|x| eval_jacobi(n, 0.0, 0.0, x));
        np_assert_allclose!(pi, p0, rtol = 1e-13);

        let pf = xs.map(|x| eval_jacobi(f64::from(n), 0.0, 0.0, x));
        np_assert_allclose!(pf, p0, rtol = 1e-12, atol = f64::EPSILON);

        let pc = xs.map(|x| eval_jacobi(f64::from(n), 0.0, 0.0, c64(x, 0.0)).re);
        // the xsf hyp2f1 implementation for complex numbers isn't very accurate
        np_assert_allclose!(pc, p0, rtol = 1e-7, atol = 1e-8);
    }
}

#[test]
fn test_eval_jacobi_scipy() {
    // test data obtained from `scipy.special.eval_jacobi`
    let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
    let ns = [0, 1, 2, 3];

    let expect_1_0 = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-7.0, -1.0, 0.35, 0.5, 0.65, 2.0, 8.0],
        [57.0, 1.0, -0.575, -0.5, -0.375, 3.0, 67.0],
        [-491.0, -1.0, -0.173_125, -0.375, -0.539_375, 4.0, 584.0],
    ];
    let expect_0_1 = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-8.0, -2.0, -0.65, -0.5, -0.35, 1.0, 7.0],
        [67.0, 3.0, -0.375, -0.5, -0.575, 1.0, 57.0],
        [-584.0, -4.0, 0.539_375, 0.375, 0.173_125, 1.0, 491.0],
    ];
    let expect_1_1 = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-10.0, -2.0, -0.2, 0.0, 0.2, 2.0, 10.0],
        [93.0, 3.0, -0.7125, -0.75, -0.7125, 3.0, 93.0],
        [-860.0, -4.0, 0.293, 0.0, -0.293, 4.0, 860.0],
    ];
    let expect_h_h = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-7.5, -1.5, -0.15, 0.0, 0.15, 1.5, 7.5],
        [61.875, 1.875, -0.6, -0.625, -0.6, 1.875, 61.875],
        [
            -535.9375, -2.1875, 0.214_375, 0.0, -0.214_375, 2.1875, 535.9375,
        ],
    ];
    let expect_h_nh = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-4.5, -0.5, 0.4, 0.5, 0.6, 1.5, 5.5],
        [33.375, 0.375, -0.435, -0.375, -0.285, 1.875, 40.875],
        [
            -275.3125, -0.3125, -0.1775, -0.3125, -0.4225, 2.1875, 337.1875,
        ],
    ];
    let expect_nh_nh = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-2.5, -0.5, -0.05, 0.0, 0.05, 0.5, 2.5],
        [18.375, 0.375, -0.3675, -0.375, -0.3675, 0.375, 18.375],
        [-151.5625, -0.3125, 0.0925, 0.0, -0.0925, 0.3125, 151.5625],
    ];
    let expect = [
        (1.0, 0.0, expect_1_0),
        (0.0, 1.0, expect_0_1),
        (1.0, 1.0, expect_1_1),
        (0.5, 0.5, expect_h_h),
        (0.5, -0.5, expect_h_nh),
        (-0.5, -0.5, expect_nh_nh),
    ];

    for (alpha, beta, expected) in expect {
        for (i, &n) in ns.iter().enumerate() {
            let p_i32 = xs.map(|x| eval_jacobi(n, alpha, beta, x));
            np_assert_allclose!(p_i32, expected[i], rtol = 1e-14, atol = f64::EPSILON);

            let p_f64 = xs.map(|x| eval_jacobi(f64::from(n), alpha, beta, x));
            np_assert_allclose!(p_f64, expected[i], rtol = 1e-14, atol = f64::EPSILON);
        }
    }
}

#[test]
fn test_eval_legendre_eq_legendre_f64() {
    // manual test cases
    let xs = [-5.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 5.0];
    let ns = [0, 1, 2, 3, 4, 8, 15];
    for n in ns {
        let p0 = xs.map(|x| xsf::legendre_p(n, x));

        let pi = xs.map(|x| eval_legendre(n, x));
        np_assert_allclose!(pi, p0, rtol = 1e-13);

        let pf = xs.map(|x| eval_legendre(f64::from(n), x));
        np_assert_allclose!(pf, p0, rtol = 1e-12, atol = f64::EPSILON);
    }
}

#[test]
fn test_eval_gegenbauer_eq_self() {
    // manual test cases
    let xs = [
        0.0,
        1.0,
        -1.0,
        1e-6,
        -1e-6,
        2.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ];
    let zs = [c64(0.5, 0.0), c64(-0.5, 0.0), c64(1.0, 0.0), c64(-1.0, 0.0)];
    let ns = [0, 1, 2, 3, 4, 11];
    let alphas = [0.0, 0.5, 1.0, 2.0];
    for &n in &ns {
        for &alpha in &alphas {
            let c_i32 = xs.map(|x| eval_gegenbauer(n, alpha, x));
            let c_f64 = xs.map(|x| eval_gegenbauer(f64::from(n), alpha, x));
            np_assert_allclose!(c_i32, c_f64, rtol = 1.5e-14, atol = 1.5e-14);

            let c_f64 = zs.map(|z| eval_gegenbauer(f64::from(n), alpha, z.re));
            let c_c64 = zs.map(|z| eval_gegenbauer(f64::from(n), alpha, z).re);
            np_assert_allclose!(c_c64, c_f64, rtol = 1e-8, atol = 1e-12);
        }
    }
}

#[test]
fn test_eval_genlaguerre_scipy() {
    // test data obtained from `scipy.special.eval_genlaguerre`
    let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
    let ns = [0, 1, 2, 3];

    let expect_0 = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [6.0, 2.0, 1.1, 1.0, 0.9, 0.0, -4.0],
        [23.5, 3.5, 1.205, 1.0, 0.805, -0.5, 3.5],
        [
            74.333_333_333_333_33,
            5.666_666_666_666_667,
            1.315_166_666_666_666_7,
            1.0,
            0.714_833_333_333_333_3,
            -0.666_666_666_666_666_7,
            2.666_666_666_666_666_7,
        ],
    ];
    let expect_1 = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [7.0, 3.0, 2.1, 2.0, 1.9, 1.0, -3.0],
        [30.5, 6.5, 3.305, 3.0, 2.705, 0.5, 0.5],
        [
            104.833_333_333_333_33,
            12.166_666_666_666_667,
            4.620_166_666_666_667,
            4.0,
            3.419_833_333_333_333_3,
            -0.166_666_666_666_666_7,
            3.166_666_666_666_666_7,
        ],
    ];
    let expect_h = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [6.5, 2.5, 1.6, 1.5, 1.4, 0.5, -3.5],
        [26.875, 4.875, 2.13, 1.875, 1.63, -0.125, 1.875],
        [
            88.645_833_333_333_33,
            8.479_166_666_666_667,
            2.642_666_666_666_666_7,
            2.1875,
            1.767_333_333_333_333_3,
            -0.604_166_666_666_666_7,
            3.229_166_666_666_666_7,
        ],
    ];
    let expect_nh = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [5.5, 1.5, 0.6, 0.5, 0.4, -0.5, -4.5],
        [20.375, 2.375, 0.53, 0.375, 0.23, -0.625, 5.375],
        [
            61.770_833_333_333_333,
            3.604_166_666_666_666_7,
            0.512_666_666_666_666_7,
            0.3125,
            0.137_333_333_333_333_33,
            -0.479_166_666_666_666_7,
            1.354_166_666_666_666_7,
        ],
    ];
    let expect = [
        (0.0, expect_0),
        (1.0, expect_1),
        (0.5, expect_h),
        (-0.5, expect_nh),
    ];

    for (alpha, expected) in expect {
        for (i, &n) in ns.iter().enumerate() {
            let p_i32 = xs.map(|x| eval_genlaguerre(n, alpha, x));
            np_assert_allclose!(p_i32, expected[i], rtol = 1e-15, atol = f64::EPSILON);

            let p_f64 = xs.map(|x| eval_genlaguerre(f64::from(n), alpha, x));
            np_assert_allclose!(p_f64, expected[i], rtol = 1e-14, atol = f64::EPSILON);
        }
    }
}

#[test]
fn test_eval_laguerre_scipy() {
    // test data obtained from `scipy.special.eval_laguerre`
    let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
    let ns = [0, 1, 2, 3, 4, 5, 6, 7, 8];

    let expect = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [6.0, 2.0, 1.1, 1.0, 0.9, 0.0, -4.0],
        [23.5, 3.5, 1.205, 1.0, 0.805, -0.5, 3.5],
        [
            74.333_333_333_333_33,
            5.666_666_666_666_667,
            1.315_166_666_666_666_7,
            1.0,
            0.714_833_333_333_333_3,
            -0.666_666_666_666_666_7,
            2.666_666_666_666_666_7,
        ],
        [
            205.375,
            8.708_333_333_333_333,
            1.430_670_833_333_333_3,
            1.0,
            0.629_337_5,
            -0.625,
            -1.291_666_666_666_666_7,
        ],
        [
            515.583_333_333_333_3,
            12.883_333_333_333_333,
            1.551_687_583_333_333_3,
            1.0,
            0.548_354_083_333_333_3,
            -0.466_666_666_666_666_7,
            -3.166_666_666_666_666_7,
        ],
        [
            1_203.743_055_555_555_6,
            18.509_722_222_222_222,
            1.678_396_334_722_222_2,
            1.0,
            0.471_728_668_055_555_6,
            -0.256_944_444_444_444_4,
            -2.090_277_777_777_777_8,
        ],
        [
            2_653.410_714_285_714_7,
            25.976_587_301_587_3,
            1.810_980_926_408_730_7,
            1.0,
            0.399_310_759_702_381_06,
            -0.040_476_190_476_190_44,
            0.325_396_825_396_825_6,
        ],
        [
            5_580.251_612_103_175,
            35.757_167_658_730_154,
            1.949_629_705_714_534_2,
            1.0,
            0.330_953_705_397_073_47,
            0.153_993_055_555_555_56,
            2.235_739_087_301_587,
        ],
    ];

    for (i, &n) in ns.iter().enumerate() {
        let p_i32 = xs.map(|x| eval_laguerre(n, x));
        np_assert_allclose!(p_i32, expect[i], rtol = 1e-15, atol = f64::EPSILON);

        let p_f64 = xs.map(|x| eval_laguerre(f64::from(n), x));
        np_assert_allclose!(p_f64, expect[i], rtol = 2e-14, atol = f64::EPSILON);
    }
}

fn test_orthogonal_poly<N, F, G>(eval_int: F, eval_f64: G, expected: &[&dyn Fn(f64) -> f64])
where
    N: num_traits::FromPrimitive,
    F: Fn(N, f64) -> f64,
    G: Fn(f64, f64) -> f64,
{
    let xs = [-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0];
    for (n, &expect_fn) in expected.iter().enumerate() {
        let y_expect = xs.map(expect_fn);
        let y_actual_int = xs.map(|x| eval_int(N::from_usize(n).unwrap(), x));
        let n_f64 = n.to_f64().unwrap();
        let y_actual_f64 = xs.map(|x| eval_f64(n_f64, x));
        np_assert_allclose!(y_actual_int, y_expect, rtol = 1e-15, atol = f64::EPSILON);
        np_assert_allclose!(y_actual_f64, y_expect, rtol = 1e-15, atol = f64::EPSILON);
    }
}

#[test]
fn test_eval_hermite_he() {
    test_orthogonal_poly::<u32, _, _>(
        eval_hermite_he,
        eval_hermite_he,
        &[&|_| 1.0, &|x| x, &|x| x * x - 1.0, &|x| x * x * x - 3.0 * x],
    );
}

#[test]
fn test_eval_hermite_h() {
    test_orthogonal_poly::<u32, _, _>(
        eval_hermite_h,
        eval_hermite_h,
        &[&|_| 1.0, &|x| 2.0 * x, &|x| 4.0 * x * x - 2.0, &|x| {
            8.0 * x * x * x - 12.0 * x
        }],
    );
}

#[test]
fn test_eval_hermite_real_degree_neg_x() {
    // mpmath; these were NaN (the tolerance is limited by the accuracy of `hyp1f1`)
    let h = [
        eval_hermite_h(2.529_868_162_602_888_5, -5.382_735_923_248_713),
        eval_hermite_h(3.3, -0.7),
        eval_hermite_h(-15.3, -8.0),
        eval_hermite_h(40.5, -3.0),
    ];
    let h_expected = [
        -22_618_949_487.496_384,
        5.292_284_105_801_645,
        9.165_934_848_962_248e29,
        -6.302_852_583_554_227e31,
    ];
    np_assert_allclose!(h, h_expected, rtol = 1e-10);

    let he = [
        eval_hermite_he(2.5, -1.5),
        eval_hermite_he(17.25, -6.0),
        eval_hermite_he(-0.5, -9.5),
    ];
    let he_expected = [
        2.234_925_546_810_356,
        110_811_784_506.491_07,
        1.824_056_114_304_247_2e19,
    ];
    np_assert_allclose!(he, he_expected, rtol = 1e-10);
}

#[test]
fn test_eval_hermite_real_degree_pos_x() {
    // mpmath; with `hypu` these were NaN or had relative errors of up to 7e-7
    let h = [
        eval_hermite_h(41.745_832_381_668_095, 7.361_880_522_197_993),
        eval_hermite_h(41.400_084_684_563_09, 5.996_470_470_350_591),
        eval_hermite_he(38.573_797_839_021_19, 8.697_439_878_678_9),
        eval_hermite_he(24.5, 3.0),
        eval_hermite_h(46.801_133_288_584_765, 3.655_511_435_592_377),
        eval_hermite_he(44.794_714_303_748_03, 5.925_757_848_223_792_5),
    ];
    let h_expected = [
        1.160_090_641_113_988e43,
        -1.200_809_248_211_909_1e38,
        3.492_391_969_655_658_3e27,
        1_012_043_943_914.947,
        4.759_278_617_763_297e38,
        1.611_367_905_027_905_8e31,
    ];
    np_assert_allclose!(h, h_expected, rtol = 1e-12);
}

#[test]
fn test_eval_chebyshev_t() {
    test_orthogonal_poly::<i32, _, _>(
        eval_chebyshev_t,
        eval_chebyshev_t,
        &[&|_| 1.0, &|x| x, &|x| 2.0 * x * x - 1.0, &|x| {
            4.0 * x * x * x - 3.0 * x
        }],
    );
}

#[test]
fn test_eval_chebyshev_u() {
    test_orthogonal_poly::<i32, _, _>(
        eval_chebyshev_u,
        eval_chebyshev_u,
        &[&|_| 1.0, &|x| 2.0 * x, &|x| 4.0 * x * x - 1.0, &|x| {
            8.0 * x * x * x - 4.0 * x
        }],
    );
}

#[test]
fn test_eval_gegenbauer_negative_alpha() {
    // mpmath: gegenbauer(n, alpha, x)
    let c = xsf::eval_gegenbauer(3_i32, -0.3, 0.5);
    np_assert_allclose!([c], [0.1505], rtol = 1e-14);
    let c = xsf::eval_gegenbauer(4_i32, -0.25, 2.0);
    np_assert_allclose!([c], [-7.093_75], rtol = 1e-14);
    // C_1^alpha(x) = 2 alpha x, also where binom(2 alpha, 1) = 0
    assert_eq!(xsf::eval_gegenbauer(1_i32, -0.5, -0.3), 0.3);
}

#[test]
fn test_eval_gegenbauer_large_alpha() {
    // mpmath: gegenbauer(n, alpha, x)
    let c = [
        eval_gegenbauer(0.5, 100.0, 1.0),
        eval_gegenbauer(2.5, 90.0, -0.7),
        eval_gegenbauer(0.5, 100.0, c64(1.0, 0.0)).re,
        eval_gegenbauer(171.5, 0.01, 1.0),
        eval_gegenbauer(170.624_376_956_302_7, 0.01, 1.0),
    ];
    let expected = [
        15.947_720_785_517_17,
        -3.247_107_273_325_162e25,
        15.947_720_785_517_17,
        1.307_068_781_116_506e-4,
        1.313_641_618_259_32e-4,
    ];
    np_assert_allclose!(c, expected, rtol = 1e-12);
}

#[test]
fn test_eval_jacobi_alpha_m1() {
    // mpmath: jacobi(n, -1 + 1e-35, beta, x)
    np_assert_allclose!([eval_jacobi(3_i32, -1.0, 1.0, 0.3)], [0.1925], rtol = 1e-14);
    np_assert_allclose!([eval_jacobi(3.0, -1.0, 1.0, 0.3)], [0.1925], rtol = 1e-14);
    let p = eval_jacobi(2.5, -1.0, -1.0, 0.3);
    np_assert_allclose!([p], [-0.218_362_414_580_039_7], rtol = 1e-14);
}

#[test]
fn test_shifted_and_scaled() {
    // reference values from scipy.special
    let r = [
        xsf::eval_gengegenbauer(3_i32, 2.5, 1.5, 0.3),
        xsf::eval_gengegenbauer(2.5, 2.5, 1.5, 0.3),
        xsf::eval_sh_legendre(4_i32, 0.3),
        xsf::eval_sh_legendre(2.5, 0.7),
        xsf::eval_chebyshev_c(5_i32, 1.3),
        xsf::eval_chebyshev_s(5_i32, 1.3),
        xsf::eval_chebyshev_c(2.5, 1.3),
        xsf::eval_chebyshev_s(2.5, 1.3),
        xsf::eval_sh_chebyshev_t(4_i32, 0.3),
        xsf::eval_sh_chebyshev_u(4_i32, 0.3),
        xsf::eval_sh_chebyshev_t(1.5, 0.3),
    ];
    let expected = [
        0.013_587_412_587_412_59,
        0.010_779_963_000_470_787,
        -0.113,
        -0.426_559_547_478_622_85,
        -0.772_07,
        -1.175_07,
        -1.108_120_029_599_681_6,
        0.157_988_300_295_633_32,
        -0.0752,
        -0.5104,
        -0.985_900_603_509_299_1,
    ];
    np_assert_allclose!(r, expected, rtol = 1e-13);

    let c = [
        xsf::eval_sh_legendre(3.0, c64(0.2, 0.1)),
        xsf::eval_sh_chebyshev_u(1.5, c64(0.3, 0.2)),
    ];
    let expected = [
        c64(0.54, 0.22),
        c64(-1.550_356_818_307_985, 0.051_799_489_871_846_57),
    ];
    np_assert_allclose!(c, expected, rtol = 1e-13);
}

#[test]
fn test_eval_sh_jacobi() {
    // P_n^(alpha, beta)(2x - 1), see scipy/scipy#18988
    let xs = [0.25, 0.5, 0.75];
    let p = xs.map(|x| xsf::eval_sh_jacobi(4_i32, 1.0, 2.0, x));
    np_assert_allclose!(p, [-1.210_937_5, 0.625, -0.585_937_5], rtol = 1e-14);
    let p = xsf::eval_sh_jacobi(2.5, 1.0, 2.0, 0.3);
    np_assert_allclose!([p], [1.303_920_745_846_762_5], rtol = 1e-14);
}

#[test]
fn test_nan_and_domain() {
    let nan = f64::NAN;
    assert!(eval_legendre(0_i32, nan).is_nan());
    assert!(eval_legendre(0.0, nan).is_nan());
    assert!(xsf::eval_chebyshev_t(0_i32, nan).is_nan());
    assert!(xsf::eval_chebyshev_u(-1_i32, nan).is_nan());
    assert!(eval_laguerre(0_i32, nan).is_nan());
    assert!(xsf::eval_genlaguerre(2_i32, nan, 1.0).is_nan());
    // alpha <= -1
    assert!(xsf::eval_genlaguerre(2_i32, -1.5, 0.3).is_nan());
    assert!(xsf::eval_genlaguerre(2.0, -1.5, 0.3).is_nan());
}
