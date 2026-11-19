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

    // but not if the starting values change `c - a - b = 0`, or the recurrence is unstable (Arb,
    // mpmath)
    let actual = [
        xsf::eval_legendre(7.345, -0.9965),
        xsf::eval_jacobi(10.3, 0.0, 30.0, -0.9),
        xsf::eval_jacobi(10.000_001, 0.0, 30.0, -0.5),
        xsf::eval_jacobi(10.000_000_01, -0.99, 40.0, -0.3),
    ];
    let desired = [
        -0.207_297_681_643_331_4,
        -1.648_343_464_373_981_2e28,
        13_912_033.739_818_575,
        34_691_102.898_867_61,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-11);

    // or near x = -1, where `hyp2f1` was inaccurate (mpmath)
    let actual = [xsf::eval_jacobi(15.8, -0.5, 8.7, -0.96), xsf::eval_gegenbauer(40.7, 9.3, -0.99)];
    let desired = [37_787_979.694_627_58, -1.002_223_611_030_761_3e15];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // and where `hyp2f1` was inaccurate for c - a - b <= -1 (mpmath)
    let actual = [
        xsf::eval_jacobi(26.39, 1.9, 3.0, -0.99),
        xsf::eval_jacobi(32.06, 3.2, 3.0, -0.88),
        xsf::eval_jacobi(17.32, 0.1, 3.0, -0.99),
    ];
    let desired = [64.703_859_369_457_89, 11.922_922_622_177_607, 1_132.803_574_409_967_7];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_eval_integer_degree_reflection() {
    // the integer degrees use the symmetry for x < 0 (Arb)
    let actual = [
        xsf::eval_legendre(742_i32, -1.0),
        xsf::eval_jacobi(3_i32, 1.0, -0.999_999_999_999_999_8, -0.999_999_999_999_99),
    ];
    np_assert_allclose!(actual, [1.0, 1.990_999_957_494_398_3e-14], rtol = 1e-14);

    // only for alpha, beta > -1: the recurrence cancels near their negative integers (Arb)
    let actual = xsf::eval_jacobi(3_i32, 1.0, -1.999_999_999_999_999_8, -0.5);
    np_assert_allclose!([actual], [-0.218_749_999_999_999_97], rtol = 1e-14);

    // and unless the reflected recurrence is nonfinite (Arb)
    let actual = xsf::eval_jacobi(2_i32, 1e100, 0.5, -0.5);
    np_assert_allclose!([actual], [3.125e198], rtol = 1e-14);
}

#[test]
fn test_eval_integer_degree_f64() {
    // integer degrees use the integer-degree recurrence, where `hyp2f1` was NaN or cancelled (Arb)
    let x = 1.0 - f64::EPSILON / 2.0; // the largest f64 below 1
    let actual = [
        xsf::eval_legendre(1843.0, -0.93),
        xsf::eval_legendre(10_000.0, x),
        xsf::eval_jacobi(3569.0, 2.6, 7.6, -0.86),
        xsf::eval_jacobi(4.0, 1.0, -0.999_999_999_999_999_1, -1.0),
        xsf::eval_jacobi(3.0, 1.0, -0.999_999_999_999_999_8, -0.999_999_999_999_99),
        xsf::eval_jacobi(3.0, 1.0, 1e30, 1.0),
        xsf::eval_jacobi(3.0, 1.0, 1e15, x),
    ];
    let desired = [
        0.010_688_711_712_761_17,
        0.999_999_994_448_329_8,
        303.059_497_866_840_5,
        2.220_446_049_250_316_5e-16,
        1.990_999_957_494_398_3e-14,
        4.0,
        3.673_067_558_944_249_6,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);

    // also for complex z (Arb)
    let actual = xsf::eval_jacobi(296.0, 7.0, 6.9, c64(0.23, 0.11));
    let desired = c64(1_666_010_717_672_526.0, -1_439_182_896_291_838.2);
    np_assert_allclose!([actual], [desired], rtol = 1e-12);

    // only for alpha, beta > -1: the recurrence cancels near their negative integers (Arb)
    let actual = [
        xsf::eval_jacobi(2.0, 1.0, -2.999_999_999_999_999_6, 0.25),
        xsf::eval_jacobi(3.0, -1.0, -1.999_999_999_999_999_8, -0.5),
    ];
    let desired = [2.015_624_999_999_999_6, -0.046_874_999_999_999_965];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // the series is used where the reflected recurrence is nonfinite (Arb)
    let actual = xsf::eval_jacobi(2.0, 1e100, 0.5, -0.999_999_999_999_973_4);
    np_assert_allclose!([actual], [8.874_685_183_736_382e171], rtol = 1e-14);

    // complex z with Re z < 0 is reflected too (Arb)
    let actual = xsf::eval_jacobi(1.0, 1e100, 0.0, c64(-0.999_999_999_999_999_1, 1e-20));
    let desired = c64(4.440_892_098_500_626e84, 5e79);
    np_assert_allclose!([actual], [desired], rtol = 1e-14);
}

#[test]
fn test_eval_integer_degree_complex() {
    // complex z with an integer degree uses the recurrence, where the series cancels (mpmath, Arb)
    let actual = [
        xsf::eval_chebyshev_t(97.0, c64(-0.160_374_780_437_770_65, -0.156_540_295_830_070_12)),
        xsf::eval_chebyshev_u(-60.0, c64(0.5, -0.75)),
        xsf::eval_gegenbauer(45.0, 3.75, c64(-0.625, 0.375)),
        xsf::eval_laguerre(50.0, c64(44.296_328_520_227_21, -0.751_454_817_007_129_3)),
        xsf::eval_genlaguerre(40.0, 2.5, c64(30.0, 10.0)),
    ];
    let desired = [
        c64(-615_681.779_834_917_3, 2_153_058.832_413_583_4),
        c64(-1.436_822_892_208_527e18, -4.290_823_921_740_318_5e18),
        c64(214_374_579_713.687_27, -1_922_101_990_974.212_2),
        c64(-194_943_212.130_009_5, 306_564_745.405_629_2),
        c64(-6_132_582_197.759_089, 5_746_717_837.856_203),
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // also for alpha near -1/2 and -1, where `alpha + 1/2` and `alpha + 1` are small (Arb, exact)
    let e = f64::EPSILON / 4.0;
    let actual = [
        xsf::eval_gegenbauer(2.0, -0.499_999_999_999_999_83, c64(0.75, 0.125)),
        xsf::eval_genlaguerre(1.0, -1.0 + 2.0 * e, c64(e, e)),
    ];
    let desired = [c64(0.226_562_499_999_999_83, -0.093_75), c64(e, -e)];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // and for large degrees, where the binomial prefactor rounded to 0 (Arb)
    let actual = xsf::eval_gegenbauer(1000.0, -0.499_999_999_999_999_9, c64(0.75, 0.125));
    let desired = c64(-7.606_673_349_914_896e74, -2.979_424_348_047_212e74);
    np_assert_allclose!([actual], [desired], rtol = 1e-14);

    // and where `alpha + 1/2 - 1` rounds to -1, the prefactor overflows for tiny alpha, or the
    // product overflows before it is scaled down (mpmath)
    let actual = [
        xsf::eval_gegenbauer(171.0, -0.499_999_999_999_999_94, c64(0.0, 0.01)),
        xsf::eval_gegenbauer(1000.0, 1e-310, c64(0.75, 0.125)),
        xsf::eval_gegenbauer(200.0, 0.2, c64(0.0, 17.43)),
    ];
    let desired = [
        c64(0.0, 9.532_255_825_720_248e-4),
        c64(8.465_527_655_970_144e-234, 4.776_235_568_205_383_5e-235),
        c64(1.081_191_070_283_941_4e306, 0.0),
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-14);

    // a fractional degree keeps the rounded prefactor, which hides the cancelled series (mpmath:
    // 2.67e-312 + 2.57e-312i and 2.67e-22 + 2.57e-22i)
    for alpha in [1e-310, 1e-20] {
        let actual = xsf::eval_gegenbauer(200.5, alpha, c64(0.0, 0.01));
        assert!(actual.norm() <= 1e-21, "{actual}");
    }

    // but not with a wrong finite value where `alpha + 1/2` is a negative integer, the
    // overflowed parameters look symmetric, or the degree is negative (mpmath, exact)
    let cases = [
        (xsf::eval_gegenbauer(172.0, -86.5, c64(0.0, 0.01)), 10.419_197_521_792_457),
        (xsf::eval_jacobi(1.0, 1.3e308, 6e307, c64(0.0, 0.0)), 3.5e307),
        (xsf::eval_gegenbauer(-2.0, 2f64.powi(-60), c64(1e100, 0.0)), 0.0),
    ];
    for (actual, desired) in cases {
        let ok = !actual.is_finite() || (actual - desired).norm() <= 1e-13 * desired;
        assert!(ok, "actual: {actual}, desired: {desired}");
    }
    // and the series where the recurrence overflows (mpmath: 12.881371905539124)
    let actual = xsf::eval_gegenbauer(171.0, 1e-310, c64(34.0, 0.0));
    assert!(actual.re.is_nan() && actual.im.is_nan(), "{actual}");
    // and nearly equal parameters aren't symmetric: like the real x, not `-(a + 2) / 4` (exact:
    // 4.109481173084667e208)
    let a = 2f64.powi(400);
    let b = f64::from_bits(a.to_bits() + 1);
    let actual = xsf::eval_jacobi(2.0, a, b, c64(0.0, 0.0));
    assert_eq!(actual, c64(xsf::eval_jacobi(2.0, a, b, 0.0), 0.0));

    // the odd degrees are zero at the origin, and imaginary on the imaginary axis (exact)
    let zero = c64(0.0, 0.0);
    assert_eq!(xsf::eval_legendre(9.0, zero), zero);
    assert_eq!(xsf::eval_legendre(9.0, c64(0.0, 1e-25)), c64(0.0, 2.460_937_5e-25));
    assert_eq!(xsf::eval_legendre(51.0, zero), zero);
    // and the even degrees real
    assert_eq!(xsf::eval_legendre(2.0, c64(-0.0, 0.01)), c64(-0.500_15, 0.0));
    assert_eq!(xsf::eval_chebyshev_u(3.0, zero), zero);
    assert_eq!(xsf::eval_gegenbauer(7.0, 3.75, zero), zero);
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

    // and with beta also near -1, where `alpha + beta + 2` cancelled (Arb)
    let actual = [
        xsf::eval_jacobi(2_i32, -0.999_999_999_999_98, -0.999_999_999_999_999_7, -0.5),
        xsf::eval_jacobi(16_i32, -0.999_999_999_999_996, -0.999_999_999_999_993, 0.58),
    ];
    let desired = [-0.187_500_000_000_003_08, -0.087_684_044_487_768_09];
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

    // also for 0 < n < 1 (Arb)
    let actual = [
        xsf::eval_genlaguerre(0.1, -0.999_999_999_999_999_8, 0.334_97),
        xsf::eval_jacobi(0.2, -0.999_999_999_999_999_7, -0.95, -0.5),
    ];
    let desired = [-0.391_281_618_882_654_9, 0.625_966_941_316_431];
    np_assert_allclose!(actual, desired, rtol = 1e-13);

    // and for subnormal n, where `1 / (n B(alpha + 1, n))` overflows (Arb: 1)
    let actual =
        [xsf::eval_genlaguerre(1e-310, -0.75, 0.0), xsf::eval_jacobi(1e-310, -0.75, 0.0, 1.0)];
    np_assert_allclose!(actual, [1.0, 1.0], rtol = 1e-15);

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
        [-535.9375, -2.1875, 0.214_375, 0.0, -0.214_375, 2.1875, 535.9375],
    ];
    let expect_h_nh = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [-4.5, -0.5, 0.4, 0.5, 0.6, 1.5, 5.5],
        [33.375, 0.375, -0.435, -0.375, -0.285, 1.875, 40.875],
        [-275.3125, -0.3125, -0.1775, -0.3125, -0.4225, 2.1875, 337.1875],
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
    let xs = [0.0, 1.0, -1.0, 1e-6, -1e-6, 2.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
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
    let expect = [(0.0, expect_0), (1.0, expect_1), (0.5, expect_h), (-0.5, expect_nh)];

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
        &[&|_| 1.0, &|x| 2.0 * x, &|x| 4.0 * x * x - 2.0, &|x| 8.0 * x * x * x - 12.0 * x],
    );
}

#[test]
fn test_eval_hermite_real_degree_neg_x() {
    // the first from Arb and mpmath, the rest from mpmath; these were NaN (the tolerance is
    // limited by the accuracy of `hyp1f1`)
    let h = [
        eval_hermite_h(2.625, -5.375),
        eval_hermite_h(3.3, -0.7),
        eval_hermite_h(-15.3, -8.0),
        eval_hermite_h(40.5, -3.0),
    ];
    let h_expected = [
        -18_568_893_765.171_894,
        5.292_284_105_801_645,
        9.165_934_848_962_248e29,
        -6.302_852_583_554_227e31,
    ];
    np_assert_allclose!(h, h_expected, rtol = 1e-10);

    let he =
        [eval_hermite_he(2.5, -1.5), eval_hermite_he(17.25, -6.0), eval_hermite_he(-0.5, -9.5)];
    let he_expected = [2.234_925_546_810_356, 110_811_784_506.491_07, 1.824_056_114_304_247_2e19];
    np_assert_allclose!(he, he_expected, rtol = 1e-10);
}

#[test]
fn test_eval_hermite_negative_integer_degree_neg_x() {
    // mpmath; H_n has no parity for n < 0, e.g. H_{-1}(x) = sqrt(pi)/2 e^(x^2) erfc(x)
    let h =
        [eval_hermite_h(-1.0, -1.999), eval_hermite_h(-3.0, -0.5), eval_hermite_h(-50.0, -10.0)];
    let h_expected = [96.160_046_802_329_64, 0.773_837_912_638_887_6, 1.007_433_655_655_245e32];
    np_assert_allclose!(h, h_expected, rtol = 1e-14);

    let he =
        [eval_hermite_he(-2.0, -2.5), eval_hermite_he(-10.0, -5.0), eval_hermite_he(-5.0, -0.1)];
    let he_expected = [142.740_828_644_349_21, 11_323_616.674_348_935, 0.194_271_792_571_189_13];
    np_assert_allclose!(he, he_expected, rtol = 1e-14);
}

#[test]
fn test_eval_hermite_negative_degree_large_neg_x() {
    // Arb; H was inf or NaN here, since the Kummer M terms overflow (or 1 / Gamma underflows)
    // before their sum comes back into range; the last two use logarithms (hence the tolerance)
    let h = [
        eval_hermite_h(-300.0, -26.0),
        eval_hermite_h(-300.0, -27.0),
        eval_hermite_h(-452.0, -31.5),
        eval_hermite_h(-400.5, -28.25),
    ];
    let h_expected = [
        7.292_328_481_376_246e116,
        1.080_723_992_709_821_7e144,
        3.202_471_036_263_673e122,
        5.967_529_707_942_022e76,
    ];
    np_assert_allclose!(h, h_expected, rtol = 1e-12);

    let he = [eval_hermite_he(-300.0, -20.0), eval_hermite_he(-452.5, -25.5)];
    let he_expected = [1.353_906_919_218_835e-105, 5.300_736_091_746_186e-184];
    np_assert_allclose!(he, he_expected, rtol = 1e-12);

    // n is too small for 1 - (1 - n) / 2 to keep it
    let tiny = [eval_hermite_h(-1e-20, -10.0), eval_hermite_he(-1e-20, -10.0)];
    np_assert_allclose!(tiny, [4.788_753_000_993_948e22, 14.130_199_579_030_07], rtol = 1e-14);
}

#[test]
fn test_eval_hermite_neg_x_overflow() {
    // mpmath; x^2 overflows, and H_n(x) -> inf / Gamma(-n) for x -> -inf
    let h = [
        eval_hermite_h(-1.5, -1e155),
        eval_hermite_h(0.5, -1e200),
        eval_hermite_h(-300.0, f64::NEG_INFINITY),
        eval_hermite_he(1.5, -2e154),
        eval_hermite_he(170.5, f64::NEG_INFINITY),
    ];
    let h_expected =
        [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    assert_eq!(h, h_expected);
}

#[test]
fn test_eval_hermite_real_degree_pos_x() {
    // Arb, mpmath; with `hypu` these were NaN or had relative errors of up to 7e-7
    let h = [
        eval_hermite_h(41.75, 7.375),
        eval_hermite_h(41.375, 6.0),
        eval_hermite_he(38.5625, 8.6875),
        eval_hermite_he(24.5, 3.0),
        eval_hermite_h(46.8125, 3.65625),
        eval_hermite_he(44.8125, 5.9375),
    ];
    let h_expected = [
        1.268_756_886_333_502_9e43,
        -9.263_448_559_927_809e37,
        1.533_866_777_770_011_2e29,
        1_012_043_943_914.947,
        4.828_573_730_463_694e38,
        1.760_567_875_576_497e31,
    ];
    np_assert_allclose!(h, h_expected, rtol = 1e-12);
}

#[test]
fn test_eval_chebyshev_t() {
    test_orthogonal_poly::<i32, _, _>(
        eval_chebyshev_t,
        eval_chebyshev_t,
        &[&|_| 1.0, &|x| x, &|x| 2.0 * x * x - 1.0, &|x| 4.0 * x * x * x - 3.0 * x],
    );
}

#[test]
fn test_eval_chebyshev_u() {
    test_orthogonal_poly::<i32, _, _>(
        eval_chebyshev_u,
        eval_chebyshev_u,
        &[&|_| 1.0, &|x| 2.0 * x, &|x| 4.0 * x * x - 1.0, &|x| 8.0 * x * x * x - 4.0 * x],
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
fn test_eval_jacobi_alpha_lt_minus_one_int_n_plus_beta() {
    // for a non-integer degree with an integer `n + beta`, alpha < -1 and 0 < beta < 1, these were
    // -2.77e8, -3.43e-3, 1.25e6 and 59.7 (Arb, Wolfram)
    let actual = [
        eval_jacobi(40.5, -1.1, 0.5, -0.6),
        eval_jacobi(35.5, -3.3, 0.5, -0.5),
        eval_jacobi(43.3, -3.2, 0.7, -0.6),
        eval_jacobi(30.5, -1.01, 0.5, -0.6),
    ];
    let desired = [
        -0.186_065_602_211_981_8,
        0.007_921_956_643_473_752,
        -0.078_826_187_314_334_7,
        0.210_905_458_121_855_35,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-12);
}

#[test]
fn test_eval_jacobi_alpha_lt_minus_one() {
    // for a non-integer degree, where the recurrence in `a` of `hyp2f1` loses most digits, these
    // were 1.32815e-4 and NaN twice (Arb, Wolfram)
    let actual = [
        eval_jacobi(
            24.445_335_888_413_27,
            -1.161_162_694_355_080_6,
            2.285_495_589_620_469,
            -0.295_343_040_857_555_7,
        ),
        eval_jacobi(
            47.743_147_644_039_865,
            -1.792_148_792_742_589_4,
            3.642_566_215_580_538,
            -0.752_167_001_210_858_3,
        ),
        eval_jacobi(
            27.237_831_793_384_256,
            -2.588_791_137_572_238,
            1.204_770_221_612_834_6,
            -0.541_727_570_598_916_4,
        ),
    ];
    let desired =
        [1.326_943_780_835_686_8e-4, -0.003_454_508_238_536_124_4, -2.685_625_233_703_272_5e-4];
    np_assert_allclose!(actual, desired, rtol = 1e-11);
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

    let c =
        [xsf::eval_sh_legendre(3.0, c64(0.2, 0.1)), xsf::eval_sh_chebyshev_u(1.5, c64(0.3, 0.2))];
    let expected = [c64(0.54, 0.22), c64(-1.550_356_818_307_985, 0.051_799_489_871_846_57)];
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
