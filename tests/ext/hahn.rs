use xsf::hahn_q;
use xsref::np_assert_allclose;

/// `(a)_k`
fn poch(a: f64, k: u32) -> f64 {
    (0..k).map(|j| a + f64::from(j)).product()
}

#[test]
fn test_hahn_q_exact() {
    // exact rational 3F2 sums
    let actual = [
        hahn_q(2, 7, 0.5, 1.25, 2.0),
        hahn_q(3, 7, 0.5, 1.25, 2.5),
        hahn_q(4, 6, -0.5, -0.5, 9.0),
        hahn_q(3, 6, -0.5, -0.5, -1.5),
        hahn_q(40, 50, 3.0, 0.25, -2.0),
        hahn_q(500, 1000, 0.5, 2.0, 300.5),
        // the recurrence in n is unstable here
        hahn_q(900, 1000, 0.5, 2.0, 50.0),
        hahn_q(900, 1000, 0.5, 2.0, 950.0),
        // and the difference equation in x from x = 0
        hahn_q(54, 54, 100.0, -0.875, 27.0),
    ];
    let desired = [
        -0.462_698_412_698_412_7,
        -0.539_445_153_061_224_5,
        337.0,
        15.0,
        334_235_228_600_388.3,
        3.922_564_427_998_434e53,
        1.860_391_435_583_974_5e79,
        1.228_530_893_036_386_2e81,
        -3.584_514_762_843_715_7e-13,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hahn_q_edges() {
    // exact rational 3F2 sums; the recurrence in n is unstable for x near 0 (and N)
    let actual = [
        hahn_q(100, 100, 0.0, 0.0, 1e-10),
        hahn_q(100, 100, 0.0, 0.0, 1e-100),
        hahn_q(100, 100, 0.0, 0.0, 100.0 - 1e-10),
        hahn_q(100, 100, 0.5, 2.0, 1e-300),
        // no premature underflow for subnormal x
        hahn_q(1000, 1000, 0.0, 0.0, -f64::from_bits(1)),
        hahn_q(1000, 1000, 0.3, 1.7, -1e-310),
    ];
    let desired = [
        -1.820_210_918_965_873_8e47,
        1.0,
        -1.820_243_403_011_934_2e47,
        1.0,
        2.024_856_663_013_668e274,
        1.847_782_799_881_995_5e287,
    ];
    np_assert_allclose!(actual, desired, rtol = 1e-13);
}

#[test]
fn test_hahn_q_domain() {
    assert_eq!(hahn_q(0, 3, 0.5, 0.5, 1.5), 1.0);
    assert_eq!(hahn_q(5, 5, 1.0, 2.0, 0.0), 1.0);
    for q in [
        hahn_q(4, 3, 0.5, 0.5, 1.0),
        hahn_q(1, 3, -1.0, 0.5, 1.0),
        hahn_q(1, 3, 0.5, f64::NAN, 1.0),
        hahn_q(1, 3, 0.5, 0.5, f64::NAN),
    ] {
        assert!(q.is_nan());
    }
}

#[test]
fn test_hahn_q_orthogonality() {
    // DLMF Table 18.19.1
    let (alpha, beta, big_n) = (0.5, 1.5, 6);
    let w = |x: u32| {
        poch(alpha + 1.0, x) * poch(beta + 1.0, big_n - x) / (poch(1.0, x) * poch(1.0, big_n - x))
    };
    for n in 0..=big_n {
        // with (-1)^n / (-N)_n = (N - n)! / N!
        let fact = |k: u32| poch(1.0, k);
        let h = poch(f64::from(n) + alpha + beta + 1.0, big_n + 1)
            * poch(beta + 1.0, n)
            * fact(n)
            * fact(big_n - n)
            / ((2.0 * f64::from(n) + alpha + beta + 1.0)
                * poch(alpha + 1.0, n)
                * fact(big_n).powi(2));
        let dots: Vec<f64> = (0..=big_n)
            .map(|m| {
                (0..=big_n)
                    .map(|x| {
                        let x_f = f64::from(x);
                        w(x) * hahn_q(n, big_n, alpha, beta, x_f)
                            * hahn_q(m, big_n, alpha, beta, x_f)
                    })
                    .sum::<f64>()
                    / h
            })
            .collect();
        let delta: Vec<f64> = (0..=big_n).map(|m| f64::from(u8::from(m == n))).collect();
        np_assert_allclose!(&dots, &delta, atol = 1e-13);
    }
}

#[test]
fn test_hahn_q_jacobi_limit() {
    // DLMF 18.21.5
    let (alpha, beta, x, big_n) = (0.5, 1.5, 0.3, 1_000_000);
    for n in 1..5 {
        let q = hahn_q(n, big_n, alpha, beta, x * f64::from(big_n));
        let nf = f64::from(n);
        let p = xsf::eval_jacobi(nf, alpha, beta, 1.0 - 2.0 * x)
            / xsf::eval_jacobi(nf, alpha, beta, 1.0);
        np_assert_allclose!([q], [p], rtol = 1e-4);
    }
}

#[test]
fn test_hahn_q_tl_moments() {
    // The weights of the sample TL-moments (Elamir & Seheult, 2003) are proportional to
    // binom(i - 1, s) binom(n - i, t) Q_{r-1}(i - s - 1; s, t, n - s - t - 1) (Hosking &
    // Balakrishnan, 2015).
    for (r, s, t, n) in [(2, 0, 0, 8), (3, 1, 1, 9), (4, 1, 2, 12)] {
        let binom = |a: u32, b: u32| xsf::binom(f64::from(a), f64::from(b));
        let ratios: Vec<f64> = (s + 1..=n - t)
            .map(|i| {
                let w: f64 = (0..r)
                    .map(|k| {
                        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
                        sign * binom(r - 1, k) * binom(i - 1, r + s - 1 - k) * binom(n - i, t + k)
                    })
                    .sum();
                let q = hahn_q(
                    r - 1,
                    n - s - t - 1,
                    s.into(),
                    t.into(),
                    f64::from(i - s - 1),
                );
                w / (binom(i - 1, s) * binom(n - i, t) * q)
            })
            .collect();
        np_assert_allclose!(&ratios, &vec![ratios[0]; ratios.len()], rtol = 1e-14);
    }
}
