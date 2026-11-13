//! Partly translated from `scipy/special/tests/test_orthogonal.py` (scipy/scipy@7d279ba)

use core::f64::consts::PI;
use xsref::np_assert_allclose;

type Rule = (Vec<f64>, Vec<f64>);

/// Like SciPy's `verify_gauss_quad`: the weights sum to `mu0`, and the rule integrates the
/// polynomials `p(k, x)` of degree `0 < k < 2n` that are orthogonal to 1 to (almost) 0, except
/// `p(n, x)`, which vanishes at the nodes
fn assert_gauss((x, w): Rule, p: impl Fn(i32, f64) -> f64, mu0: f64, rtol: f64) {
    np_assert_allclose!([w.iter().sum::<f64>()], [mu0], rtol = rtol);
    let n = i32::try_from(x.len()).unwrap();
    for k in (1..2 * n).filter(|&k| k != n) {
        let terms = x.iter().zip(&w).map(|(&xi, wi)| wi * p(k, xi));
        let (sum, abs) = terms.fold((0.0, 0.0), |(s, a), t| (s + t, a + t.abs()));
        assert!(
            sum.abs() <= rtol * abs,
            "k = {k}: {sum:e} (sum of |terms| = {abs:e})"
        );
    }
}

#[test]
fn test_roots_legendre() {
    // A&S Table 25.4
    let (x, w) = xsf::roots_legendre(3);
    let r = (3.0_f64 / 5.0).sqrt();
    np_assert_allclose!(x, [-r, 0.0, r], rtol = 1e-15);
    np_assert_allclose!(w, [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0], rtol = 1e-15);

    for n in [5, 50, 150] {
        assert_gauss(xsf::roots_legendre(n), xsf::eval_legendre, 2.0, 1e-13);
    }
    let p = |k, x| xsf::eval_sh_legendre(k, x);
    assert_gauss(xsf::roots_sh_legendre(25), p, 1.0, 1e-13);
    // (x + 1) / 2 would cancel near 0 (Arb); 5e-15 covers eigenvalues and f/df off by 8 ulps
    let (x, _) = xsf::roots_sh_legendre(2000);
    np_assert_allclose!([x[0]], [3.612_684_148_443_298e-7], rtol = 5e-15);
}

#[test]
fn test_roots_jacobi() {
    // the example in the SciPy docs
    let (x, w) = xsf::roots_jacobi(3, 0.5, 0.5);
    np_assert_allclose!(x, [-0.5_f64.sqrt(), 0.0, 0.5_f64.sqrt()], atol = 1e-15);
    np_assert_allclose!(w, [PI / 8.0, PI / 4.0, PI / 8.0], rtol = 1e-15);

    for (a, b) in [(-0.5, -0.75), (0.9, 2.0), (18.24, 27.3), (1.0, 658.0)] {
        let p = |k, x| xsf::eval_jacobi(k, a, b, x);
        let mu0 = (a + b + 1.0).exp2() * xsf::beta(a + 1.0, b + 1.0);
        assert_gauss(xsf::roots_jacobi(25, a, b), p, mu0, 1e-12);
    }

    // P_n^(0, 0) = P_n, and P_n^(a, a) is proportional to C_n^(a + 1/2)
    assert_eq!(xsf::roots_jacobi(6, 0.0, 0.0), xsf::roots_legendre(6));
    assert_eq!(
        xsf::roots_jacobi(6, 4.0, 4.0),
        xsf::roots_gegenbauer(6, 4.5)
    );
}

#[test]
fn test_roots_sh_jacobi() {
    let (a, b) = (0.5, -0.5);
    let p = |k, x| xsf::eval_sh_jacobi(k, a, b, x);
    assert_gauss(xsf::roots_sh_jacobi(25, a, b), p, PI / 2.0, 1e-13);

    // G_n^(p, q) is proportional to P_n^*(p - q, q - 1)
    let (p, q) = (47.1, 0.2);
    let expected = xsf::roots_sh_jacobi(25, p - q, q - 1.0);
    assert_eq!(xsf::roots_gengegenbauer(25, p, q), expected);
}

/// Parameters near the ends of their domains (Arb)
#[test]
fn test_roots_jacobi_extreme() {
    // the recurrence coefficients would cancel
    let (x, w) = xsf::roots_jacobi(3, -0.999_999_999_999, -0.999_999_999_9);
    let desired = [
        -0.999_999_999_966_666_7,
        4.125_001_266_302_403e-11,
        0.999_999_999_999_666_7,
    ];
    np_assert_allclose!(x, desired, atol = 1e-15);
    let desired = [
        4_999_999_585.981_552,
        1.333_333_333_292_010_5,
        500_011_061_139.089_5,
    ];
    np_assert_allclose!(w, desired, rtol = 1e-14);

    // the largest node rounds to 1
    let (_, w) = xsf::roots_jacobi(2, -0.999_999_999_999_999_9, 0.0);
    np_assert_allclose!(w, [1.125, 9_007_199_254_740_992.0], rtol = 1e-15);

    // 2^(alpha + beta + 1) would overflow
    let (x, w) = xsf::roots_sh_jacobi(1, 1024.0, 0.0);
    np_assert_allclose!(x, [1.0 / 1026.0], rtol = 1e-13);
    np_assert_allclose!(w, [1.0 / 1025.0], rtol = 1e-15);
    // the largest weight overflows, but mu0 = 2^1035 / 1035 would overflow for all
    let (_, w) = xsf::roots_jacobi(2, 1034.0, 0.0);
    assert_eq!(w[0], f64::INFINITY);
    np_assert_allclose!([w[1]], [5.227_583_763_348_631e307], rtol = 1e-15);
    // mu0 = 2^(1000 - 2e-8) B(1001, 2e-8) overflows, since B(1001, 2e-8) ~ 5e7
    let (_, w) = xsf::roots_jacobi(2, 1000.0, -0.999_999_98);
    np_assert_allclose!([w[1]], [2.686_802_069_732_627e300], rtol = 1e-13);
    let (_, w) = xsf::roots_sh_jacobi(3, 1500.0, 0.0);
    let desired = [
        4.732_855_623_184_923e-4,
        1.859_542_816_943_685e-4,
        6.982_674_308_258_493e-6,
    ];
    np_assert_allclose!(w, desired, rtol = 1e-14);

    // ln(mu0) = 19001 ln(2) + betaln(10001, 9001) would cancel
    let (_, w) = xsf::roots_jacobi(1, 10_000.0, 9_000.0);
    np_assert_allclose!(w, [4_933_761_225.304_621], rtol = 1e-14);

    // a root within 2e-14 of -1 needs more Newton steps; 2e-12 covers eigenvalues and f/df off
    // by 8 ulps
    let (_, w) = xsf::roots_jacobi(2000, -0.999_999_98, -0.999_999_97);
    np_assert_allclose!([w[0]], [16_666_659.794_355_836], rtol = 2e-12);

    // C_n^(alpha + 1/2) would cancel
    let (x, w) = xsf::roots_jacobi(3, -0.999_999_999_999_999_9, -0.999_999_999_999_999_9);
    assert_eq!(x, [-1.0, 0.0, 1.0]);
    let desired = [4_503_599_627_370_496.0, 4.0 / 3.0, 4_503_599_627_370_496.0];
    np_assert_allclose!(w, desired, rtol = 1e-15);

    // the smallest node of `roots_jacobi` rounds to -1
    let (x, w) = xsf::roots_sh_jacobi(1, 10.0, -0.999_999_999_999_999_9);
    np_assert_allclose!(x, [1.009_293_658_750_142_3e-17], rtol = 1e-15);
    np_assert_allclose!(w, [9_007_199_254_740_989.0], rtol = 1e-15);
    let (x, _) = xsf::roots_sh_jacobi(4, 10.0, -0.999_999_999_999_999_9);
    np_assert_allclose!([x[0]], [1.982_541_115_402_065e-18], rtol = 1e-15);
}

#[test]
fn test_roots_gegenbauer() {
    for alpha in [-0.25, 0.1, 10.0, 170.5, 512.5] {
        let p = |k, x| xsf::eval_gegenbauer(k, alpha, x);
        let mu0 = xsf::beta(0.5, alpha + 0.5);
        assert_gauss(xsf::roots_gegenbauer(25, alpha), p, mu0, 1e-12);
    }
    // the limit alpha -> 0
    assert_eq!(xsf::roots_gegenbauer(5, 0.0), xsf::roots_chebyshev_t(5));
    let ((x, w), (y, v)) = (xsf::roots_gegenbauer(2, 1e-16), xsf::roots_chebyshev_t(2));
    np_assert_allclose!(x, y, rtol = 1e-15);
    np_assert_allclose!(w, v, rtol = 1e-15);
    let (x, w) = xsf::roots_gegenbauer(1, 1e-20);
    assert_eq!(x, [0.0]);
    np_assert_allclose!(w, [PI], rtol = 1e-15);
    // `C_n` would underflow
    assert_eq!(xsf::roots_gegenbauer(3, 1e-310), xsf::roots_chebyshev_t(3));

    // a known limitation, see `roots_gegenbauer` (the true weights are about 2^53, 4/3, 2^53)
    let (x, w) = xsf::roots_gegenbauer(3, -0.499_999_999_999_999_94);
    np_assert_allclose!(x, [-1.0, 0.0, 1.0], atol = 1e-15);
    assert!(w[0].is_nan() && w[2].is_nan());
}

#[test]
fn test_roots_chebyshev_t() {
    assert_gauss(xsf::roots_chebyshev_t(25), xsf::eval_chebyshev_t, PI, 1e-14);

    // test_chebyt_symmetry
    let (x, _) = xsf::roots_chebyshev_t(21);
    let neg: Vec<f64> = x[..10].iter().rev().map(|x| -x).collect();
    assert_eq!(x[11..], neg);
    assert_eq!(x[10], 0.0);
}

#[test]
#[should_panic(expected = "n must be a positive integer")]
fn test_roots_legendre_n_zero() {
    let _ = xsf::roots_legendre(0);
}

#[test]
#[should_panic(expected = "alpha and beta must be greater than -1")]
fn test_roots_jacobi_beta_invalid() {
    let _ = xsf::roots_jacobi(3, 1.0, -2.0);
}

#[test]
#[should_panic(expected = "(p - q) must be greater than -1, and q must be greater than 0")]
fn test_roots_gengegenbauer_q_invalid() {
    let _ = xsf::roots_gengegenbauer(3, 2.0, -1.0);
}

#[test]
#[should_panic(expected = "alpha and beta must be greater than -1")]
fn test_roots_gengegenbauer_q_tiny() {
    // q - 1 rounds to -1
    let _ = xsf::roots_gengegenbauer(1, 0.0, 1e-20);
}

#[test]
#[should_panic(expected = "alpha must be greater than -0.5")]
fn test_roots_gegenbauer_alpha_invalid() {
    let _ = xsf::roots_gegenbauer(3, -0.75);
}

#[test]
fn test_roots_laguerre() {
    for alpha in [-0.5, 0.1, 10.0, 50.0] {
        let p = |k, x| xsf::eval_genlaguerre(k, alpha, x);
        let mu0 = xsf::gamma(alpha + 1.0);
        assert_gauss(xsf::roots_genlaguerre(25, alpha), p, mu0, 1e-12);
    }
    assert_eq!(xsf::roots_laguerre(5), xsf::roots_genlaguerre(5, 0.0));
    assert_eq!(xsf::roots_genlaguerre(1, 2.0), (vec![3.0], vec![2.0]));

    // the derivative would cancel (Arb)
    let (_, w) = xsf::roots_genlaguerre(2, -0.999_999_999_999_999_9);
    np_assert_allclose!(w, [9_007_199_254_740_991.0, 0.25], rtol = 1e-15);
    // mu0 / sum(w) would overflow (Arb)
    let (_, w) = xsf::roots_genlaguerre(2, 170.0);
    np_assert_allclose!(
        w,
        [3.905_394_184_064_997e306, 3.352_021_431_243_001_6e306],
        rtol = 1e-14
    );
    // mu0 = gamma(alpha + 1) would overflow (Arb)
    let (_, w) = xsf::roots_genlaguerre(2, 170.7);
    np_assert_allclose!(
        w,
        [1.427_030_187_170_566_4e308, 1.225_209_102_672_075_3e308],
        rtol = 1e-14
    );
    // mu0 / 2^alpha would also overflow (Arb)
    let (_, w) = xsf::roots_genlaguerre(200, 220.0);
    np_assert_allclose!([w[199]], [1.603_580_611_056_378_5e171], rtol = 1e-13);
}

#[test]
fn test_roots_hermite() {
    let sqrt_pi = PI.sqrt();
    let p = |k, x| xsf::eval_hermite_h(k as u32, x);
    assert_gauss(xsf::roots_hermite_h(25), p, sqrt_pi, 1e-12);
    let p = |k, x| xsf::eval_hermite_he(k as u32, x);
    assert_gauss(xsf::roots_hermite_he(25), p, (2.0 * PI).sqrt(), 1e-12);
}

/// test_roots_hermite_asy
#[test]
fn test_roots_hermite_asy() {
    for n in [151, 300, 301, 2000] {
        let (x, w) = xsf::roots_hermite_h(n);
        np_assert_allclose!([w.iter().sum::<f64>()], [PI.sqrt()], rtol = 1e-15);

        // the normalized Hermite functions `psi_k(x) = H_k(x) e^(-x^2 / 2) / sqrt(2^k k! sqrt(pi))`
        // of degree `n` vanish at the nodes
        for xi in x {
            let (mut p0, mut p1) = (0.0, PI.powf(-0.25) * (-xi * xi / 2.0).exp());
            for k in 1..=n {
                let k = k as f64;
                (p0, p1) = (p1, (2.0 / k).sqrt() * xi * p1 - ((k - 1.0) / k).sqrt() * p0);
            }
            assert!(p1.abs() <= 1e-12, "n = {n}, x = {xi}: {p1:e}");
        }
    }

    // He_n(x) = 2^(-n/2) H_n(x / sqrt(2))
    let (x, w) = xsf::roots_hermite_he(200);
    let (y, v) = xsf::roots_hermite_h(200);
    let sqrt_2 = 2.0_f64.sqrt();
    np_assert_allclose!(
        x,
        y.iter().map(|y| y * sqrt_2).collect::<Vec<_>>(),
        rtol = 1e-15
    );
    np_assert_allclose!(
        w,
        v.iter().map(|v| v * sqrt_2).collect::<Vec<_>>(),
        rtol = 1e-15
    );
}

#[test]
#[should_panic(expected = "alpha must be greater than -1")]
fn test_roots_genlaguerre_alpha_invalid() {
    let _ = xsf::roots_genlaguerre(3, -1.1);
}
