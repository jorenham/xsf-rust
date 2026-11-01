#![allow(clippy::float_cmp)]

// there are no xsref tables, so we use smoke tests instead

use core::f64::consts;
use num_complex::c64;
const I: num_complex::Complex<f64> = num_complex::Complex::new(0.0, 1.0);
const LN_3: f64 = 1.098_612_288_668_109_8_f64;

// legendre_p

#[test]
fn test_legendre_p_f64() {
    assert_eq!(
        [
            xsf::legendre_p(0, 0.0),
            xsf::legendre_p(1, 0.0),
            xsf::legendre_p(2, 0.0),
        ],
        [1.0, 0.0, -0.5],
    );
}

#[test]
fn test_legendre_p_c64() {
    assert_eq!(
        [
            xsf::legendre_p(0, I),
            xsf::legendre_p(1, I),
            xsf::legendre_p(2, I),
        ],
        [c64(1.0, 0.0), I, c64(-2.0, 0.0)],
    );
}

// legendre_p_all

#[test]
fn test_legendre_p_all_f64() {
    assert_eq!(xsf::legendre_p_all(2, 0.0), vec![1.0, 0.0, -0.5]);
}

#[test]
fn test_legendre_p_all_c64() {
    assert_eq!(
        xsf::legendre_p_all(2, I),
        vec![c64(1.0, 0.0), I, c64(-2.0, 0.0)],
    );
}

// sph_legendre_p

#[test]
fn test_sph_legendre_p_f64() {
    assert_eq!(
        [
            xsf::sph_legendre_p(0, 0, 0.0),
            xsf::sph_legendre_p(1, 0, 0.0),
            xsf::sph_legendre_p(1, 1, 0.0),
        ],
        [0.282_094_791_773_878_14, 0.488_602_511_902_919_87, 0.0],
    );
}

#[test]
fn test_sph_legendre_p_c64() {
    xsref::np_assert_allclose!(
        &[
            xsf::sph_legendre_p(0, 0, I),
            xsf::sph_legendre_p(1, 0, I),
            xsf::sph_legendre_p(1, 1, I),
        ],
        &[
            c64(0.282_094_791_773_878_2, 0.0),
            c64(0.753_953_074_239_480_4, 0.0),
            c64(-0.406_025_136_855_663_4, 0.0),
        ],
        atol = f64::EPSILON
    );
}

// sph_legendre_p_all

#[test]
fn test_sph_legendre_p_all_f64() {
    assert_eq!(
        xsf::sph_legendre_p_all(1, 1, 0.0),
        vec![
            vec![0.282_094_791_773_878_14, 0.0, 0.0],
            vec![0.488_602_511_902_919_87, 0.0, 0.0],
        ]
    );
}

#[test]
fn test_sph_legendre_p_all_c64() {
    let pnm = xsf::sph_legendre_p_all(1, 1, I);
    assert_eq!(pnm.len(), 2);

    xsref::np_assert_allclose!(
        &pnm[0],
        &[
            c64(0.282_094_791_773_878_2, 0.0),
            c64(0.0, 0.0),
            c64(0.0, 0.0)
        ],
        atol = f64::EPSILON
    );
    xsref::np_assert_allclose!(
        &pnm[1],
        &[
            c64(0.753_953_074_239_480_4, 0.0),
            c64(-0.406_025_136_855_663_4, 0.0),
            c64(0.406_025_136_855_663_4, 0.0),
        ],
        atol = f64::EPSILON
    );
}

// assoc_legendre_p

#[test]
fn test_assoc_legendre_p_f64() {
    assert_eq!(
        [
            xsf::assoc_legendre_p(0, 0, 0.0),
            xsf::assoc_legendre_p(1, 0, 0.0),
            xsf::assoc_legendre_p(0, 1, 0.0),
            xsf::assoc_legendre_p(1, 1, 0.0),
        ],
        [1.0, 0.0, 0.0, -1.0],
    );
}

#[test]
fn test_assoc_legendre_p_c64() {
    assert_eq!(
        [
            xsf::assoc_legendre_p(0, 0, I),
            xsf::assoc_legendre_p(1, 0, I),
            xsf::assoc_legendre_p(0, 1, I),
            xsf::assoc_legendre_p(1, 1, I),
        ],
        [c64(1.0, 0.0), I, c64(0.0, 0.0), c64(-consts::SQRT_2, 0.0)],
    );
}

// assoc_legendre_p_all

#[test]
fn test_assoc_legendre_p_all_f64() {
    assert_eq!(
        xsf::assoc_legendre_p_all(1, 1, 0.0),
        vec![vec![1.0, 0.0, 0.0], vec![0.0, -1.0, 0.5]],
    );
}

#[test]
fn test_assoc_legendre_p_all_c64() {
    assert_eq!(
        xsf::assoc_legendre_p_all(1, 1, I),
        vec![
            vec![c64(1.0, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)],
            vec![
                c64(0.0, 1.0),
                c64(-consts::SQRT_2, 0.0),
                c64(consts::FRAC_1_SQRT_2, 0.0),
            ],
        ]
    );
}

// assoc_legendre_p_norm

#[test]
fn test_assoc_legendre_p_norm_f64() {
    xsref::np_assert_allclose!(
        &[
            xsf::assoc_legendre_p_norm(0, 0, 0.0),
            xsf::assoc_legendre_p_norm(1, 0, 0.0),
            xsf::assoc_legendre_p_norm(0, 1, 0.0),
            xsf::assoc_legendre_p_norm(1, 1, 0.0),
        ],
        &[
            consts::FRAC_1_SQRT_2,
            0.0,
            0.0,
            -0.866_025_403_784_438_6_f64, // -sqrt(3) / 2
        ],
        atol = f64::EPSILON
    );
}

#[test]
fn test_assoc_legendre_p_norm_c64() {
    xsref::np_assert_allclose!(
        &[
            xsf::assoc_legendre_p_norm(0, 0, I),
            xsf::assoc_legendre_p_norm(1, 0, I),
            xsf::assoc_legendre_p_norm(0, 1, I),
            xsf::assoc_legendre_p_norm(1, 1, I),
        ],
        &[
            c64(consts::FRAC_1_SQRT_2, 0.0),
            c64(0.0, 1.224_744_871_391_589),
            c64(0.0, 0.0),
            c64(-1.224_744_871_391_589_2, 0.0),
        ],
        atol = f64::EPSILON
    );
}

// legendre_q_all

#[test]
fn test_legendre_q_all_f64() {
    let (qn, qd) = xsf::legendre_q_all(4, 0.5);

    assert_eq!(qn.len(), 5);
    assert_eq!(qd.len(), 5);

    assert_eq!(qn[..2], [LN_3 * 0.5, LN_3 * 0.25 - 1.0]);
    assert_eq!(qd[..2], [4.0 / 3.0, 2.0 / 3.0 + 0.5 * LN_3]);
}

#[test]
fn test_legendre_q_all_c64() {
    let (qn, qd) = xsf::legendre_q_all(4, c64(0.0, 1.0));

    assert_eq!(qn.len(), 5);
    assert_eq!(qd.len(), 5);

    assert_eq!(
        qn[..2],
        [
            c64(0.0, consts::FRAC_PI_4),
            c64(-1.0 - consts::FRAC_PI_4, 0.0),
        ]
    );
}

#[test]
fn test_legendre_q_all_n0() {
    // computed with n = 1, since `lqn` writes at least 2 terms
    let (qn, qd) = xsf::legendre_q_all(1, 1.5);
    assert_eq!(
        xsf::legendre_q_all(0, 1.5),
        (qn[..1].to_vec(), qd[..1].to_vec())
    );

    let z = c64(1.5, 0.0);
    let (qn, qd) = xsf::legendre_q_all(1, z);
    assert_eq!(
        xsf::legendre_q_all(0, z),
        (qn[..1].to_vec(), qd[..1].to_vec())
    );
}

// assoc_legendre_q_all

#[test]
fn assoc_test_legendre_q_all_f64() {
    assert_eq!(
        xsf::assoc_legendre_q_all(1, 2, 0.0),
        (
            vec![
                vec![0.0, -1.0, 0.0], // order 0,1,2 for degree 0
                vec![-1.0, 0.0, 2.0], // order 0,1,2 for degree 1
            ],
            vec![
                vec![1.0, 0.0, 2.0],  // order 0,1,2 for degree 0
                vec![0.0, -2.0, 0.0], // order 0,1,2 for degree 1
            ],
        ),
    );
}

#[test]
fn assoc_test_legendre_q_all_c64() {
    assert_eq!(
        xsf::assoc_legendre_q_all(1, 2, c64(0.0, 0.0)),
        (
            vec![
                vec![c64(0.0, 0.0), c64(-1.0, 0.0), c64(0.0, 0.0)],
                vec![c64(-1.0, 0.0), c64(0.0, 0.0), c64(2.0, 0.0)],
            ],
            vec![
                vec![c64(1.0, 0.0), c64(0.0, 0.0), c64(2.0, 0.0)],
                vec![c64(0.0, 0.0), c64(-2.0, 0.0), c64(0.0, 0.0)],
            ],
        ),
    );
}

#[test]
fn test_assoc_legendre_q_all_n0_m0() {
    // computed with n = 1 and m = 1, since `lqmn` writes at least 2 x 2 terms
    macro_rules! check {
        ($z:expr) => {
            let (qm, qd) = xsf::assoc_legendre_q_all(1, 1, $z);
            for (n, m) in [(0, 0), (0, 1), (1, 0)] {
                let sub = |q: &[Vec<_>]| q[..=n].iter().map(|r| r[..=m].to_vec()).collect();
                assert_eq!(xsf::assoc_legendre_q_all(n, m, $z), (sub(&qm), sub(&qd)));
            }
        };
    }
    check!(1.5);
    check!(0.5);
    check!(c64(1.5, 0.0));
    check!(c64(0.3, 0.4));
}

#[test]
fn test_assoc_legendre_q_all_int_zero() {
    // `static_cast<T>(ls * (j + i - 1) * (j - i + 2))` with `ls = -1, i = 2, j = 0` is `+0.0`
    let q = xsf::assoc_legendre_q_all(2, 2, c64(1.00005, 0.0)).0[0][2];
    assert_eq!(q.im.to_bits(), 0.0_f64.to_bits());
}

// deviations from xsf

#[test]
fn test_large_n_normalized() {
    // NaN in xsf because of `int` overflow, for `n >= 646`
    let p = xsf::sph_legendre_p(700, 0, 0.5);
    let expected = -0.337_154_106_114_745_7; // mpmath
    assert!(((p - expected) / expected).abs() < 1e-13, "{p}");

    let p = xsf::assoc_legendre_p_norm(700, 0, 0.3);
    let expected = 0.802_058_211_710_411_2; // mpmath
    assert!(((p - expected) / expected).abs() < 1e-13, "{p}");
}

#[test]
fn test_assoc_legendre_q_all_near_minus_one() {
    // garbage in xsf, since `(int)log(x - 1)` is used for `-1.1 <= x < -1`
    for x in [-1.0001, -1.05, -1.1] {
        let (qm, qd) = xsf::assoc_legendre_q_all(3, 2, x);
        let (cqm, cqd) = xsf::assoc_legendre_q_all(3, 2, c64(x, 0.0));
        for (q, cq) in qm
            .iter()
            .flatten()
            .chain(qd.iter().flatten())
            .zip(cqm.iter().flatten().chain(cqd.iter().flatten()))
        {
            assert_eq!(*q, cq.re);
        }
    }
    // mpmath: legenq(2, 1, -1.05, type=3) = 2.211389887900321167... (with the opposite sign)
    let q = xsf::assoc_legendre_q_all(2, 1, -1.05).0[2][1];
    assert!((q + 2.211_389_887_900_321_2).abs() < 1e-15, "{q}");
}
