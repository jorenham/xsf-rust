//! Translated from `scipy/special/tests/test_ellip_harm.py`

use xsref::np_assert_allclose;

/// The closed forms of `E_n^p(s)` for `n <= 3`
fn ellip_harm_known(h2: f64, k2: f64, n: u32, p: u32, s: f64) -> f64 {
    let s2 = s * s;
    let (sh, sk) = ((s2 - h2).abs().sqrt(), (s2 - k2).abs().sqrt());
    let root2 = ((h2 + k2) * (h2 + k2) - 3.0 * h2 * k2).abs().sqrt();
    let root3 = (4.0 * (h2 + k2) * (h2 + k2) - 15.0 * h2 * k2).sqrt();
    let root3h = ((h2 + 2.0 * k2) * (h2 + 2.0 * k2) - 5.0 * h2 * k2)
        .abs()
        .sqrt();
    let root3k = ((2.0 * h2 + k2) * (2.0 * h2 + k2) - 5.0 * h2 * k2)
        .abs()
        .sqrt();
    match (n, p) {
        (0, 1) => 1.0,
        (1, 1) => s,
        (1, 2) => sh,
        (1, 3) => sk,
        (2, 1) => s2 - 1.0 / 3.0 * ((h2 + k2) + root2),
        (2, 2) => s2 - 1.0 / 3.0 * ((h2 + k2) - root2),
        (2, 3) => s * sh,
        (2, 4) => s * sk,
        (2, 5) => sh * sk,
        (3, 1) => s * s2 - (s / 5.0) * (2.0 * (h2 + k2) + root3),
        (3, 2) => s * s2 - (s / 5.0) * (2.0 * (h2 + k2) - root3),
        (3, 3) => sh * (s2 - 1.0 / 5.0 * ((h2 + 2.0 * k2) + root3h)),
        (3, 4) => sh * (s2 - 1.0 / 5.0 * ((h2 + 2.0 * k2) - root3h)),
        (3, 5) => sk * (s2 - 1.0 / 5.0 * ((2.0 * h2 + k2) + root3k)),
        (3, 6) => sk * (s2 - 1.0 / 5.0 * ((2.0 * h2 + k2) - root3k)),
        (3, 7) => s * sh * sk,
        _ => unreachable!(),
    }
}

#[test]
fn test_ellip_harm() {
    for (h2, k2, s) in [
        (0.5, 2.0, 0.2),
        (5.0, 8.0, 2.5),
        (3.0, 100.0, 12.0),
        (1e-3, 4e3, 0.7),
    ] {
        for n in 0..=3 {
            for p in 1..=2 * n + 1 {
                np_assert_allclose!(
                    [xsf::ellip_harm(h2, k2, n, p, s, 1.0, 1.0)],
                    [ellip_harm_known(h2, k2, n, p, s)],
                    rtol = 1e-12
                );
            }
        }
    }
}

#[test]
fn test_ellip_harm_signs() {
    // n = 2, p = 5 is of type N, with the sign `signm * signn`
    let e = xsf::ellip_harm(5.0, 8.0, 2, 5, 2.5, 1.0, 1.0);
    assert_eq!(xsf::ellip_harm(5.0, 8.0, 2, 5, 2.5, -1.0, -1.0), e);
    assert_eq!(xsf::ellip_harm(5.0, 8.0, 2, 5, 2.5, -1.0, 1.0), -e);
    assert!(xsf::ellip_harm(5.0, 8.0, 2, 5, 2.5, 0.5, 1.0).is_nan());
}

#[test]
fn test_ellip_harm_invalid_p() {
    assert!(xsf::ellip_harm(0.5, 2.0, 4, 10, 0.2, 1.0, 1.0).is_nan());
    assert!(xsf::ellip_harm(0.5, 2.0, 4, 0, 0.2, 1.0, 1.0).is_nan());
}

#[test]
fn test_ellip_harm_graded() {
    // the eigenvector of the coefficients is graded, with `|v_0 / v_1| ~ 1e20`
    np_assert_allclose!(
        [xsf::ellip_harm(1e-40, 1.0, 2, 1, 2.0, 1.0, 1.0)],
        [10.0 / 3.0],
        rtol = 1e-15
    );
}

#[test]
fn test_ellip_harm_large_s() {
    // `(s^2 - h2) (s^2 - k2)` overflows
    np_assert_allclose!(
        [xsf::ellip_harm(1.0, 4.0, 2, 5, 1e90, 1.0, 1.0)],
        [9.999_999_999_999_999e179],
        rtol = 1e-15
    );
}

#[test]
fn test_ellip_harm_s_near_k() {
    // `s^2 - k2` cancels
    np_assert_allclose!(
        [xsf::ellip_harm(
            1.0,
            2.0,
            1,
            3,
            1.414_213_562_374_509_3,
            1.0,
            1.0
        )],
        [2.000_060_244_115_435e-6],
        rtol = 1e-14
    );
}

#[test]
fn test_ellip_harm_extreme() {
    // `(-h2)^10` underflows, and `s^2` overflows
    np_assert_allclose!(
        [
            xsf::ellip_harm(1e-40, 2e-40, 20, 1, 1.0, 1.0, 1.0),
            xsf::ellip_harm(1e20, 4e20, 1, 2, 1e160, 1.0, 1.0),
        ],
        [1.0, 1e160],
        rtol = 1e-15
    );
}
