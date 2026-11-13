//! Translated from `scipy/special/tests/test_ellip_harm.py`

use core::f64::consts::PI;
use xsref::np_assert_allclose;

/// The closed forms of `gamma_n^p` for `n <= 3`
fn ellip_normal_known(h2: f64, k2: f64, n: u32, p: u32) -> f64 {
    let (hk, d) = (h2 * k2, k2 - h2);
    let (s2, s3, s4) = (
        h2 * h2 + k2 * k2,
        h2.powi(3) + k2.powi(3),
        h2.powi(4) + k2.powi(4),
    );
    let g2 = |sign: f64| {
        let r = (h2 * h2 + k2 * k2 - hk).sqrt();
        16.0 * PI / 405.0
            * (2.0 * s4 - 4.0 * hk * s2
                + 6.0 * hk * hk
                + sign * r * (2.0 * s3 - 3.0 * hk * (h2 + k2)))
    };
    let g3 = |sign: f64| {
        let r = (4.0 * s2 - 7.0 * hk).sqrt();
        16.0 * PI / 13125.0
            * hk
            * (16.0 * s4 - 36.0 * hk * s2
                + 46.0 * hk * hk
                + sign * r * (8.0 * s3 - 11.0 * hk * (h2 + k2)))
    };
    let g3h = |sign: f64| {
        let r = (h2 * h2 + 4.0 * k2 * k2 - hk).sqrt();
        let p = 6.0 * h2.powi(3) + 8.0 * k2.powi(3) - 9.0 * h2 * h2 * k2 - 13.0 * h2 * k2 * k2;
        let q =
            6.0 * h2.powi(4) + 16.0 * k2.powi(4) - 12.0 * h2.powi(3) * k2 - 28.0 * h2 * k2.powi(3)
                + 34.0 * hk * hk;
        16.0 * PI / 13125.0 * h2 * d * (q + sign * r * p)
    };
    let g3k = |sign: f64| {
        let r = (4.0 * h2 * h2 + k2 * k2 - hk).sqrt();
        let p = 8.0 * h2.powi(3) + 6.0 * k2.powi(3) - 13.0 * h2 * h2 * k2 - 9.0 * h2 * k2 * k2;
        let q =
            16.0 * h2.powi(4) + 6.0 * k2.powi(4) - 28.0 * h2.powi(3) * k2 - 12.0 * h2 * k2.powi(3)
                + 34.0 * hk * hk;
        16.0 * PI / 13125.0 * k2 * d * (q + sign * r * p)
    };
    match (n, p) {
        (0, 1) => 4.0 * PI,
        (1, 1) => 4.0 * PI * hk / 3.0,
        (1, 2) => 4.0 * PI * h2 * d / 3.0,
        (1, 3) => 4.0 * PI * k2 * d / 3.0,
        (2, 1) => g2(1.0),
        (2, 2) => g2(-1.0),
        (2, 3) => 4.0 * PI * h2 * hk * d / 15.0,
        (2, 4) => 4.0 * PI * hk * k2 * d / 15.0,
        (2, 5) => 4.0 * PI * hk * d * d / 15.0,
        (3, 1) => g3(1.0),
        (3, 2) => g3(-1.0),
        (3, 3) => g3h(1.0),
        (3, 4) => g3h(-1.0),
        (3, 5) => g3k(1.0),
        (3, 6) => g3k(-1.0),
        (3, 7) => 4.0 * PI * hk * hk * d * d / 105.0,
        _ => unreachable!(),
    }
}

#[test]
fn test_ellip_normal() {
    for (h2, k2) in [(5.0, 8.0), (0.3, 30.0), (1e3, 1.1e3)] {
        for n in 0..=3 {
            for p in 1..=2 * n + 1 {
                np_assert_allclose!(
                    [xsf::ellip_normal(h2, k2, n, p)],
                    [ellip_normal_known(h2, k2, n, p)],
                    rtol = 1e-12
                );
            }
        }
    }
}

#[test]
fn test_ellip_harm_2() {
    // sum_p F_1^p / (3 E_1^p) = 1 / (s sqrt((s^2 - h^2)(s^2 - k^2)))
    let (h2, k2, s) = (5.0, 8.0, 10.0);
    let i1: f64 = (1..=3)
        .map(|p| {
            xsf::ellip_harm_2(h2, k2, 1, p, s) / (3.0 * xsf::ellip_harm(h2, k2, 1, p, s, 1.0, 1.0))
        })
        .sum();
    np_assert_allclose!(
        [i1],
        [1.0 / (s * ((s * s - h2) * (s * s - k2)).sqrt())],
        rtol = 1e-14
    );

    // values produced by code from arXiv:1204.0267
    let expected = [
        0.001_080_568_533_82,
        0.001_058_205_138_09,
        0.001_060_583_847_43,
        0.001_067_744_923_06,
        0.001_079_763_564_54,
    ];
    for (p, y) in (1..=5).zip(expected) {
        np_assert_allclose!(
            [xsf::ellip_harm_2(5.0, 8.0, 2, p, 10.0)],
            [y],
            atol = 1.5e-14,
            rtol = 0.0
        );
    }
}

#[test]
fn test_ellip_harm_2_s_near_k() {
    // at and near the endpoint singularity of the integral, for an exact `k = 2`; mpmath
    np_assert_allclose!(
        [
            xsf::ellip_harm_2(1.0, 4.0, 1, 2, 2.0),
            xsf::ellip_harm_2(1.0, 4.0, 2, 1, 2.0),
            xsf::ellip_harm_2(1.0, 4.0, 0, 1, 2.000_000_000_000_01),
            xsf::ellip_harm_2(1.0, 4.0, 1, 3, 2.000_000_000_002),
        ],
        [
            0.703_730_314_407_309_7,
            0.924_100_945_358_312_6,
            0.842_875_119_056_628_4,
            0.866_023_328_389_255,
        ],
        rtol = 1e-14
    );
}

#[test]
fn test_ellip_harm_2_s_near_rounded_k() {
    // `k = sqrt(k2)` is rounded, and `s - k ~ 1e-12`; mpmath
    np_assert_allclose!(
        [xsf::ellip_harm_2(
            1.0,
            1.000_000_02,
            0,
            1,
            1.000_000_010_001
        )],
        [10.240_061_121_404_25],
        rtol = 1e-14
    );
}

#[test]
fn test_ellip_harm_2_large_s() {
    // `E_4^1(s)` overflows, and the result (~1e-450) underflows
    assert_eq!(xsf::ellip_harm_2(1.0, 4.0, 4, 1, 1e90), 0.0);
    // `F_1^1(s) ~ 1 / s^2`, where the integral `~ 1 / (3 s^3)` alone underflows
    np_assert_allclose!(
        [xsf::ellip_harm_2(1.0, 4.0, 1, 1, 1e110)],
        [1e-220],
        rtol = 1e-14
    );
    // `E_2^5(s)^2 ~ s^8` overflows; mpmath
    np_assert_allclose!(
        [xsf::ellip_harm_2(1.0, 4.0, 2, 5, 1e90)],
        [1e-270],
        rtol = 1e-14
    );
}

#[test]
fn test_ellip_normal_h_near_k() {
    // `t^2 - h^2` and `t^2 - k^2` cancel for `t in [h, k]`
    let (h2, k2) = (1.0, 1.000_000_001);
    let expected = 4.0 * PI * h2 * k2 * (k2 - h2) * (k2 - h2) / 15.0;
    np_assert_allclose!([xsf::ellip_normal(h2, k2, 2, 5)], [expected], rtol = 1e-14);
    // `1 - t^2 / h2` cancels; mpmath
    np_assert_allclose!(
        [xsf::ellip_normal(1.0, 1.000_000_02, 2, 1)],
        [3.351_032_264_526_018e-16],
        rtol = 1e-14
    );
}

#[test]
fn test_ellip_harm_2_invalid() {
    assert!(xsf::ellip_harm_2(5.0, 8.0, 2, 6, 10.0).is_nan());
    assert!(xsf::ellip_normal(5.0, 8.0, 2, 0).is_nan());
}

/// The potential `1 / |x1 - x2|` as a series of ellipsoidal harmonics
#[test]
fn test_ellip_potential() {
    let cartesian = |l: f64, m: f64, n: f64, h2: f64, k2: f64| {
        let x = (l * l * m * m * n * n / (h2 * k2)).sqrt();
        let y = ((l * l - h2) * (m * m - h2) * (h2 - n * n) / (h2 * (k2 - h2))).sqrt();
        let z = ((l * l - k2) * (k2 - m * m) * (k2 - n * n) / (k2 * (k2 - h2))).sqrt();
        [x, y, z]
    };
    for [l1, m1, n1, l2, m2, n2, h2, k2] in [
        [
            120.0,
            19.0_f64.sqrt(),
            2.0,
            41.0,
            17.0_f64.sqrt(),
            2.0,
            15.0,
            25.0,
        ],
        [
            120.0,
            16.0_f64.sqrt(),
            3.2,
            21.0,
            11.0_f64.sqrt(),
            2.9,
            11.0,
            20.0,
        ],
    ] {
        let e = |n, p, s| xsf::ellip_harm(h2, k2, n, p, s, 1.0, 1.0);
        let mut sum = 0.0_f64;
        for n in 0..20 {
            let xsum: f64 = (1..=2 * n + 1)
                .map(|p| {
                    let inner = e(n, p, l2) * e(n, p, m2) * e(n, p, n2);
                    let outer = xsf::ellip_harm_2(h2, k2, n, p, l1) * e(n, p, m1) * e(n, p, n1);
                    4.0 * PI * inner * outer
                        / (xsf::ellip_normal(h2, k2, n, p) * f64::from(2 * n + 1))
                })
                .sum();
            if xsum.abs() < 1e-9 * sum.abs() {
                break;
            }
            sum += xsum;
        }
        let (x1, x2) = (cartesian(l1, m1, n1, h2, k2), cartesian(l2, m2, n2, h2, k2));
        let dist = x1
            .iter()
            .zip(x2)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        np_assert_allclose!([sum], [1.0 / dist], rtol = 1e-8);
    }
}
