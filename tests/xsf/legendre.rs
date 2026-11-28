// there are no xsref tables, so we use smoke tests instead

use core::f64::consts;
use num_complex::c64;

use crate::DENORM_MIN;
const I: num_complex::Complex<f64> = num_complex::Complex::new(0.0, 1.0);
const LN_3: f64 = 1.098_612_288_668_109_8_f64;

// Same sign and within one least subnormal, for results Arb rounds just off a midpoint of the
// leading term
fn assert_near_subnormal(actual: f64, expected: f64) {
    assert_eq!(actual.is_sign_negative(), expected.is_sign_negative(), "{actual:e}");
    assert!((actual - expected).abs() <= DENORM_MIN, "{actual:e} != {expected:e}");
}

// legendre_p

#[test]
fn test_legendre_p_f64() {
    assert_eq!(
        [xsf::legendre_p(0, 0.0), xsf::legendre_p(1, 0.0), xsf::legendre_p(2, 0.0),],
        [1.0, 0.0, -0.5],
    );
}

#[test]
fn test_legendre_p_c64() {
    assert_eq!(
        [xsf::legendre_p(0, I), xsf::legendre_p(1, I), xsf::legendre_p(2, I),],
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
    assert_eq!(xsf::legendre_p_all(2, I), vec![c64(1.0, 0.0), I, c64(-2.0, 0.0)],);
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
fn test_sph_legendre_p_scaled_diagonal() {
    // Arb applied to the spherical normalization of P_n^m(cos(theta)).
    let expected = -5.828_985_299_854_025e-252;
    let scalar = xsf::sph_legendre_p(600, 475, 0.2);
    let all = xsf::sph_legendre_p_all(600, 475, 0.2);
    xsref::np_assert_allclose!(
        [scalar, all[600][475], -all[600][476]],
        [expected; 3],
        rtol = 1e-12
    );
}

#[test]
fn test_sph_legendre_p_c64() {
    xsref::np_assert_allclose!(
        &[xsf::sph_legendre_p(0, 0, I), xsf::sph_legendre_p(1, 0, I), xsf::sph_legendre_p(1, 1, I)],
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
        vec![vec![0.282_094_791_773_878_14, 0.0, 0.0], vec![0.488_602_511_902_919_87, 0.0, 0.0],]
    );
}

#[test]
fn test_sph_legendre_p_all_c64() {
    let pnm = xsf::sph_legendre_p_all(1, 1, I);
    assert_eq!(pnm.len(), 2);

    xsref::np_assert_allclose!(
        &pnm[0],
        &[c64(0.282_094_791_773_878_2, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)],
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

#[test]
fn test_assoc_legendre_p_negative_degree() {
    // P_n^m = P_{-n-1}^m (DLMF 14.9.5); mpmath `legenp(n, m, z, type=2)`
    xsref::np_assert_allclose!(
        [
            xsf::assoc_legendre_p(-2, 0, 0.5),
            xsf::assoc_legendre_p(-3, 1, 0.5),
            xsf::assoc_legendre_p(-3, -1, 0.5),
            xsf::assoc_legendre_p(-4, 2, -0.3),
        ],
        [0.5, -1.299_038_105_676_658, 0.216_506_350_946_109_65, -4.095_000_000_000_001],
        rtol = f64::EPSILON
    );
    xsref::np_assert_allclose!(
        [xsf::assoc_legendre_p(-3, 1, c64(0.5, 0.2))],
        [c64(-1.410_682_430_103_991, -0.370_033_672_257_135_8)],
        rtol = f64::EPSILON
    );
    xsref::np_assert_allclose!(
        [xsf::legendre_p(-1, 0.5), xsf::legendre_p(-3, 0.5)],
        [1.0, -0.125],
        rtol = f64::EPSILON
    );
}

#[test]
fn test_assoc_legendre_p_zero_and_infinity() {
    // P_n^m(0) = 0 for odd n + m (DLMF 14.5.1), also where the diagonal (-1)^m (2m-1)!! overflows;
    // P_n(z) ~ z^n as z -> +-inf
    let inf = f64::INFINITY;
    assert_eq!(
        [
            xsf::assoc_legendre_p(201, 200, 0.0),
            xsf::assoc_legendre_p(153, 152, 0.0),
            xsf::assoc_legendre_p(200, 200, 0.0),
            xsf::assoc_legendre_p(201, 201, 0.0),
            xsf::assoc_legendre_p(3, 0, -inf),
            xsf::legendre_p(3, inf),
            xsf::legendre_p(3, -inf),
            xsf::legendre_p(4, -inf),
        ],
        [0.0, 0.0, inf, -inf, -inf, inf, -inf, inf],
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
            vec![c64(0.0, 1.0), c64(-consts::SQRT_2, 0.0), c64(consts::FRAC_1_SQRT_2, 0.0),],
        ]
    );
}

// assoc_legendre_p_norm

#[test]
fn test_assoc_legendre_p_norm_small_diagonal() {
    // Arb references at the exact binary64 arguments.
    for (n, m, x, expected) in [
        (2100, 790, 0.925, 3.345_564_164_239_196),
        (2000, 750, 0.925, 2.961_386_133_824_163_4),
        (600, 475, 0.98, -3.169_217_138_093_830_4e-251),
    ] {
        let value = xsf::assoc_legendre_p_norm(n, m, x);
        xsref::np_assert_allclose!([value], [expected], rtol = 4e-13);
    }
    let table = xsf::assoc_legendre_p_norm_all(600, 475, 0.98);
    let value = xsf::assoc_legendre_p_norm(600, 475, 0.98);
    assert_eq!(table[600][475], value);
    assert_eq!(table[600][476], -value);
}

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

    xsref::np_assert_allclose!(&qn[..2], &[LN_3 * 0.5, LN_3 * 0.25 - 1.0], rtol = f64::EPSILON);
    xsref::np_assert_allclose!(&qd[..2], &[4.0 / 3.0, 2.0 / 3.0 + 0.5 * LN_3], rtol = f64::EPSILON);
}

#[test]
fn test_legendre_q_all_c64() {
    let (qn, qd) = xsf::legendre_q_all(4, c64(0.0, 1.0));

    assert_eq!(qn.len(), 5);
    assert_eq!(qd.len(), 5);

    assert_eq!(qn[..2], [c64(0.0, consts::FRAC_PI_4), c64(-1.0 - consts::FRAC_PI_4, 0.0),]);
}

#[test]
fn test_legendre_q_all_c64_pole_at_one() {
    for sign in [-1.0, 1.0] {
        let z = c64(1.0, sign * 0.25);
        let (q, d) = xsf::legendre_q_all(0, z);
        // Q_0(z) = log(1 - 8i) / 2 on the upper side and Q_0'(z) = 1 / (1 - z^2).
        let expected = c64(65.0_f64.ln() / 4.0, -sign * 8.0_f64.atan() / 2.0);
        let expected_derivative = c64(16.0 / 65.0, sign * (128.0 / 65.0));
        xsref::np_assert_allclose!([q[0], d[0]], [expected, expected_derivative], rtol = 1e-15);

        let (q, _) = xsf::legendre_q_all(0, c64(1.0, sign * DENORM_MIN));
        let expected = c64(1075.0 * consts::LN_2 / 2.0, sign * consts::FRAC_PI_4);
        xsref::np_assert_allclose!([q[0]], [expected], rtol = f64::EPSILON);

        let (q, d) = xsf::legendre_q_all(0, c64(1.0, sign * 0.0));
        assert_eq!([q[0], d[0]], [c64(1e300, 0.0); 2]);
    }
}

#[test]
fn test_legendre_q_c64_degree_zero_derivative_near_poles() {
    let overflow_ordinate = 1_u64 << 49;
    let first_finite = f64::from_bits(0x7fef_ffff_ffff_fff0);
    for x in [-1.0, 1.0] {
        for side in [-1.0, 1.0] {
            for (y, magnitude) in [
                (DENORM_MIN, f64::INFINITY),
                (f64::from_bits(overflow_ordinate - 1), f64::INFINITY),
                (f64::from_bits(overflow_ordinate), f64::INFINITY),
                (f64::from_bits(overflow_ordinate + 1), first_finite),
                (1e-10, 5e9),
            ] {
                let z = c64(x, side * y);
                let (_, qd) = xsf::legendre_q_all(0, z);
                let (_, qmd) = xsf::assoc_legendre_q_all(0, 0, z);
                let expected = c64(0.25, x * side * magnitude);
                assert_eq!(qd[0], expected, "{z:?}");
                assert_eq!(qmd[0][0], expected, "{z:?}");
            }
        }
    }
}

#[test]
fn test_legendre_q_all_backward_seeds() {
    // mpmath legenq with type=3, reflected by Q_k(-x) = (-1)^(k+1) Q_k(x).
    let (negative, _) = xsf::legendre_q_all(100, -1.02101);
    let (positive, _) = xsf::legendre_q_all(200, 1.02101);
    xsref::np_assert_allclose!(
        [negative[0], positive[0], negative[100]],
        [-2.283_177_075_289_996_7, 2.283_177_075_289_996_7, -3.205_778_187_395_51e-10],
        rtol = 1e-13
    );
    assert_eq!(xsf::legendre_q_all(3000, 1.2).0[3000], 0.0);
}

#[test]
fn test_legendre_q_all_c64_q0_subnormal_near_pole() {
    // At y = minsub, the real magnitude rounds to log(2 / y) / 2.
    let re = 1075.0 * consts::LN_2 / 2.0;
    for lip in [-1.0, 1.0] {
        let y = lip * DENORM_MIN;
        let im = lip * consts::FRAC_PI_4;
        let (q, _) = xsf::legendre_q_all(0, c64(-1.0, y));
        xsref::np_assert_allclose!([q[0]], [c64(-re, im)], rtol = f64::EPSILON);
        for sign in [-1.0, 1.0] {
            let (q, _) = xsf::assoc_legendre_q_all(0, 0, c64(sign, y));
            xsref::np_assert_allclose!([q[0][0]], [c64(sign * re, im)], rtol = f64::EPSILON);
        }
    }
}

#[test]
fn test_legendre_q_all_f64_reflection() {
    // Q_n(-x) = (-1)^(n+1) Q_n(x) and Q_n'(-x) = (-1)^n Q_n'(x), see DLMF 14.24.2
    for x in [1.5, 951.41] {
        let (qn, qd) = xsf::legendre_q_all(30, x);
        let (qn_neg, qd_neg) = xsf::legendre_q_all(30, -x);
        for n in 0..=30 {
            let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
            assert_eq!(qn_neg[n], sign * qn[n], "n = {n}, x = -{x}");
            assert_eq!(qd_neg[n], -sign * qd[n], "n = {n}, x = -{x}");
        }
    }
}

#[test]
fn test_legendre_q_all_n0() {
    // computed with n = 1, since `lqn` writes at least 2 terms
    let (qn, qd) = xsf::legendre_q_all(1, 1.5);
    assert_eq!(xsf::legendre_q_all(0, 1.5), (qn[..1].to_vec(), qd[..1].to_vec()));

    let z = c64(1.5, 0.0);
    let (qn, qd) = xsf::legendre_q_all(1, z);
    assert_eq!(xsf::legendre_q_all(0, z), (qn[..1].to_vec(), qd[..1].to_vec()));
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
fn test_assoc_legendre_q_all_q0_tiny_and_large_x() {
    // mpmath: legenq(0, 0, x, type=2), and type=3 for |x| > 1
    for (x, q0) in [(1e-10, 1e-10), (-1e-6, -1.000_000_000_000_333_3e-6)] {
        let (qn, _) = xsf::legendre_q_all(1, x);
        xsref::np_assert_allclose!([qn[0]], [q0], rtol = 5e-16);
    }
    let (qm, _) = xsf::assoc_legendre_q_all(1, 1, 1e6);
    xsref::np_assert_allclose!([qm[0][0]], [1.000_000_000_000_333_3e-6], rtol = 5e-16);
}

#[test]
fn test_legendre_q_all_extreme_complex() {
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            // Arb at exact binary64 arguments.
            let z = c64(sx * 1e160, sy * 1e159);
            let q1 = c64(f64::from_bits(0x28f), -sx * sy * f64::from_bits(0x84));
            let dq0 = c64(-f64::from_bits(0x7ac), sx * sy * f64::from_bits(0x18d));
            let (q, dq) = xsf::legendre_q_all(1, z);
            assert_eq!((q[1], dq[0]), (q1, dq0));
            let (q, dq) = xsf::assoc_legendre_q_all(1, 0, z);
            assert_eq!((q[1][0], dq[0][0]), (q1, dq0));

            let z = c64(sx * 1e160, sy * 0.0);
            let (q, dq) = xsf::legendre_q_all(1, z);
            assert_eq!(q[1].re.to_bits(), 0x2a3);
            assert_eq!(q[1].im.to_bits(), (-sx * sy * 0.0).to_bits());
            assert_eq!(dq[0].re.to_bits(), 0x8000_0000_0000_07e8);
            assert_eq!(dq[0].im.to_bits(), (sx * sy * 0.0).to_bits());

            let z = c64(sx * 2.0_f64.powi(511), sy * 0.0);
            let (q, dq) = xsf::assoc_legendre_q_all(1, 0, z);
            assert_eq!(q[1][0].re.to_bits(), 0x0005_5555_5555_5555);
            assert_eq!(q[1][0].im.to_bits(), (-sx * sy * 0.0).to_bits());
            assert_eq!(dq[0][0].re.to_bits(), (-f64::MIN_POSITIVE).to_bits());
            assert_eq!(dq[0][0].im.to_bits(), (sx * sy * 0.0).to_bits());

            let z = c64(sx * 0.0, sy * 1e160);
            let (q, dq) = xsf::legendre_q_all(1, z);
            let (qm, dqm) = xsf::assoc_legendre_q_all(1, 0, z);
            for (q, dq) in [(q[0], dq[0]), (qm[0][0], dqm[0][0])] {
                assert_eq!(q.re.to_bits(), (sy * 0.0).to_bits());
                assert_eq!(q.im, -sy * 1e-160);
                assert_eq!(dq.re, 1e-320);
                assert_eq!(dq.im.to_bits(), (sx * sy * 0.0).to_bits());
            }
            for q in [q[1], qm[1][0]] {
                assert_eq!(q.re, -f64::from_bits(0x2a3));
                assert_eq!(q.im.to_bits(), (-sx * sy * 0.0).to_bits());
            }
        }
    }

    // The small imaginary-coordinate correction rounds the derivative down to max subnormal.
    let z = c64(2.0_f64.powi(511), 5.048_054_719_547_867e145);
    let (_, dq) = xsf::legendre_q_all(0, z);
    assert_eq!(dq[0].re.to_bits(), 0x800f_ffff_ffff_ffff);

    let x = 2.0_f64.powi(537);
    let (q, dq) = xsf::legendre_q_all(1, c64(x, x));
    // Here the leading imaginary derivative is the half-minsub midpoint.
    assert_eq!(dq[0].re.to_bits(), 0);
    assert_eq!(dq[0].im.to_bits(), 0);
    assert_eq!(q[1].re.to_bits(), (-0.0_f64).to_bits());
    assert_eq!(q[1].im.to_bits(), (-0.0_f64).to_bits());

    // Arb: the real part cancels to about 2^-51 of the imaginary one.
    let x = 1.25 * 2.0_f64.powi(479);
    // Two ulps above x.
    let (_, dq) = xsf::legendre_q_all(0, c64(x, f64::from_bits(x.to_bits() + 2)));
    assert_eq!(dq[0].re, 4.666_318_092_516_091e-305);
}

#[test]
fn test_legendre_q_all_extreme_degree_two() {
    let check = |z: num_complex::Complex64, expected: [f64; 4], order_one: [f64; 4]| {
        let (q, d) = xsf::legendre_q_all(2, z);
        assert_eq!(
            [q[2].re, q[2].im, d[2].re, d[2].im].map(f64::to_bits),
            expected.map(f64::to_bits)
        );
        let (q, d) = xsf::assoc_legendre_q_all(2, 73, z);
        for (order, (&q, &d)) in q[2].iter().zip(&d[2]).enumerate() {
            let expected = if order % 2 == 0 { expected } else { order_one };
            assert_eq!([q.re, q.im, d.re, d.im].map(f64::to_bits), expected.map(f64::to_bits));
        }
    };
    let scale = 2.0_f64.powi(479);
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            // Arb at the exact binary64 arguments. The last pair has x^2 - 2xy - y^2 = -1
            // before scaling, so the order-zero derivative's real projection is positive.
            for (x, y, q_re, q_im, d_re, d_im) in [
                (1.0, 0.1, 0.0, -0.0, -0.0, 0.0),
                (1.0, 1.0, -0.0, -0.0, 0.0, -0.0),
                (4_217_293_152_016_490.0, 1_746_860_020_068_409.0, 0.0, -0.0, 0.0, 0.0),
            ] {
                check(
                    c64(sx * x * scale, sy * y * scale),
                    [sx * q_re, sy * q_im, d_re, sx * sy * d_im],
                    [-q_re, -sx * sy * q_im, -sx * d_re, -sy * d_im],
                );
            }
            check(
                c64(sx * 0.0, sy * scale),
                [0.0, sy * 0.0, -0.0, 0.0],
                [0.0, -sx * sy * 0.0, sx * 0.0, 0.0],
            );
        }
    }
}

#[test]
fn test_assoc_legendre_q_all_extreme_order_one() {
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            // Arb at the exact binary64 arguments.
            let (q, dq) = xsf::assoc_legendre_q_all(0, 1, c64(sx * 1e160, sy * 1e159));
            let value = c64(-9.900_990_099_009_902e-161, sx * sy * 9.900_990_099_009_9e-162);
            let slope = c64(sx * f64::from_bits(0x7ac), -sy * f64::from_bits(0x18d));
            assert_eq!((q[0][1], dq[0][1]), (value, slope));

            let (q, dq) = xsf::assoc_legendre_q_all(0, 1, c64(sx * 0.0, sy * 1e160));
            assert_eq!(q[0][1].re.to_bits(), 0);
            assert_eq!(q[0][1].im, sx * sy * 1e-160);
            assert_eq!(dq[0][1].re, -sx * f64::from_bits(0x7e8));
            assert_eq!(dq[0][1].im.to_bits(), 0);
        }
    }

    // Round tiny projections, including ties where the first tail correction matters.
    for (x, y, bits) in [
        (1.207_381_461_363_012_6e208, 2.0_f64.powi(308), 1),
        (3.375, 3.0 * 2.0_f64.powi(536), (-DENORM_MIN).to_bits()),
        (1.0 - f64::EPSILON / 2.0, 2.0_f64.powi(511), 0x800f_ffff_ffff_ffff),
    ] {
        let (q, _) = xsf::assoc_legendre_q_all(0, 1, c64(x, y));
        let projection = if x < y { q[0][1].re } else { q[0][1].im };
        assert_eq!(projection.to_bits(), bits);
    }
    let (q, _) = xsf::assoc_legendre_q_all(0, 1, c64(3.0 * 2.0_f64.powi(536), 1.125));
    assert_near_subnormal(q[0][1].im, DENORM_MIN);
}

#[test]
fn test_legendre_q_extreme_derivative_subnormal_midpoint() {
    // Arb at the exact binary64 inputs: the imaginary magnitude is just below
    // the leading value of 3.5 times the smallest subnormal.
    let a = 5991.0_f64;
    let major = a * 2.0_f64.powi(467);
    let minor = 7.0 * a.powi(3) * 2.0_f64.powi(325);
    for (x, y) in [(major, minor), (minor, major)] {
        for (sx, sy) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)] {
            let z = c64(sx * x, sy * y);
            let (_, derivative) = xsf::legendre_q_all(0, z);
            let (_, associated_derivative) = xsf::assoc_legendre_q_all(0, 1, z);
            let magnitude = f64::from_bits(3);
            assert_near_subnormal(derivative[0].im, sx * sy * magnitude);
            assert_near_subnormal(associated_derivative[0][1].im, -sy * magnitude);
        }
    }
}

#[test]
fn test_assoc_legendre_q_all_extreme_order_two() {
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            // Exact rational order-two values and derivatives at binary64 inputs.
            let (q, dq) = xsf::assoc_legendre_q_all(1, 2, c64(sx * 1e160, sy * 1e159));
            let value = c64(sx * 1.980_198_019_801_980_3e-160, -sy * 1.980_198_019_801_98e-161);
            let slope = c64(-f64::from_bits(0xf59), sx * sy * f64::from_bits(0x31a));
            assert_eq!((q[0][2], dq[0][2]), (value, slope));
            assert_eq!(q[1][2], -slope);
            assert_eq!(dq[1][2].re.to_bits(), (-sx * 0.0).to_bits());
            assert_eq!(dq[1][2].im.to_bits(), (sy * 0.0).to_bits());

            let (q, dq) = xsf::assoc_legendre_q_all(1, 2, c64(sx * 0.0, sy * 1e160));
            assert_eq!(q[0][2].re.to_bits(), (sx * 0.0).to_bits());
            assert_eq!(q[0][2].im, -sy * 2e-160);
            assert_eq!(dq[0][2].re.to_bits(), 0xfd0);
            assert_eq!(dq[0][2].im.to_bits(), 0);
            assert_eq!(q[1][2].re, -f64::from_bits(0xfd0));
            assert_eq!(q[1][2].im.to_bits(), 0);
            assert_eq!(dq[1][2].re.to_bits(), 0);
            assert_eq!(dq[1][2].im.to_bits(), (-sy * 0.0).to_bits());
        }
    }

    // The omitted terms resolve the half-minsubnormal derivative in opposite directions.
    let x = 2.0_f64.powi(538);
    let (_, dq) = xsf::assoc_legendre_q_all(0, 2, c64(x, 0.5));
    assert_eq!(dq[0][2].re.to_bits(), (-DENORM_MIN).to_bits());
    let (q, dq) = xsf::assoc_legendre_q_all(0, 2, c64(x, 1.0));
    assert_near_subnormal(dq[0][2].re, -0.0);
    assert_eq!(q[0][2].im.to_bits(), (-DENORM_MIN).to_bits());

    // The real Q_1^2 midpoint changes direction at 1/sqrt(3).
    let cases = [(0.577_350_269_189_625_7, DENORM_MIN), (0.577_350_269_189_625_8, 0.0)];
    for (y, expected) in cases {
        let (q, _) = xsf::assoc_legendre_q_all(1, 2, c64(x, y));
        assert_near_subnormal(q[1][2].re, expected);
    }

    // The leading imaginary derivative at this input is 3.5 minsubnormals.
    let a = 5991.0_f64;
    let (_, dq) = xsf::assoc_legendre_q_all(
        0,
        2,
        c64(a * 2.0_f64.powi(467), 7.0 * a.powi(3) * 2.0_f64.powi(324)),
    );
    assert_near_subnormal(dq[0][2].im, f64::from_bits(3));
}

#[test]
fn test_assoc_legendre_q_all_extreme_degree_one_order_one() {
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            // Arb at the exact binary64 arguments and at both imaginary-axis lips.
            let (q, dq) = xsf::assoc_legendre_q_all(1, 1, c64(sx * 1e160, sy * 1e159));
            assert_eq!(q[1][1], c64(-sx * f64::from_bits(0x51e), sy * f64::from_bits(0x109)));
            assert_eq!(dq[1][1].re.to_bits(), 0);
            assert_eq!(dq[1][1].im.to_bits(), (-sx * sy * 0.0).to_bits());

            let (q, dq) = xsf::assoc_legendre_q_all(1, 1, c64(sx * 0.0, sy * 1e160));
            assert_eq!(q[1][1].re, sx * f64::from_bits(0x545));
            assert_eq!(q[1][1].im.to_bits(), (sy * 0.0).to_bits());
            assert_eq!(dq[1][1].re.to_bits(), 0);
            assert_eq!(dq[1][1].im.to_bits(), (sx * sy * 0.0).to_bits());
        }
    }

    // The axis approximation is a 3.5-minsubnormal imaginary projection.
    let a = 5991.0_f64;
    let z = c64(a * 2.0_f64.powi(467), 21.0 * a.powi(3) * 2.0_f64.powi(324));
    let (q, _) = xsf::assoc_legendre_q_all(1, 1, z);
    assert_near_subnormal(q[1][1].im, f64::from_bits(3));
}

// Arb at the exact binary64 inputs, with `Q_n^m(z) = sqrt(z^2 - 1)^m d^m Q_n(z) / dz^m`.

#[test]
fn test_assoc_legendre_q_all_extreme_imaginary_axis() {
    let (q, _) = xsf::assoc_legendre_q_all(1, 74, c64(0.0, 1e155));
    xsref::np_assert_allclose!(
        [q[0][3].im, q[0][40].im, q[1][3].re, q[1][40].re, q[1][74].re],
        [
            6e-155,
            -8.159_152_832_478_977e-108,
            8e-310,
            -1.115_084_220_438_793_6e-261,
            -8.269_713_603_798_466e-202,
        ],
        rtol = 1e-15
    );

    // The real part is about 1e-355 relative to the imaginary one.
    let (q, _) = xsf::assoc_legendre_q_all(0, 170, c64(1e-200, 1e155));
    xsref::np_assert_allclose!(
        [q[0][121].re, q[0][170].re],
        [-8.094_298_525_273_46e-310, 7.257_415_615_307_999e-204],
        rtol = 1e-15
    );
}

#[test]
fn test_assoc_legendre_q_all_extreme_dispatch_boundary_and_diagonal() {
    // The next order need not underflow at the dispatch boundary (6.69 least subnormals).
    let x = 2.0_f64.powi(479);
    let (q, _) = xsf::assoc_legendre_q_all(2, 74, c64(x, 1.0));
    assert_eq!(q[2][74].re, f64::from_bits(7));

    // On the diagonal the leading real parts vanish, and the next terms in 1/z^2 give them.
    let (q, dq) = xsf::assoc_legendre_q_all(1, 200, c64(x, x));
    assert_eq!(q[1][200].re, -8.903_420_595_702_894e-198);
    assert_eq!(dq[0][180].re, 1.370_823_570_613_042_2e-244);
}

#[test]
fn test_assoc_legendre_q_all_extreme_double_double_resolution() {
    // Re z^-4 is 7e-32 of its modulus, near the double-double resolution.
    let x = 2.0_f64.powi(479);
    let z = c64(4_217_293_152_016_490.0 * x, 1_746_860_020_068_409.0 * x);
    let (q, dq) = xsf::assoc_legendre_q_all(3, 333, z);
    xsref::np_assert_allclose!(
        [dq[2][200].re, dq[2][333].re, q[3][200].re],
        [1.686_547_975_285_891e-292, -6.090_345_402_686_098e30, -1.630_329_709_443_028e-291],
        rtol = 1e-15
    );
}

#[test]
fn test_assoc_legendre_q_all_extreme_pell_point() {
    // x^2 - 3y^2 = 1 before scaling, so Re z^-3 cancels to about 2^-105.
    let x = 2.0_f64.powi(479);
    let z = c64(5_170_128_475_599_457.0 * x, 2_984_975_067_132_296.0 * x);
    let (q, dq) = xsf::assoc_legendre_q_all(2, 200, z);
    xsref::np_assert_allclose!(
        [q[2][200].re, q[2][200].im, dq[1][200].re, dq[1][200].im],
        [
            6.410_887_451_325_874e-134,
            -2.638_328_859_165_686_3e-102,
            -3.173_706_659_072_214_6e-135,
            1.306_103_395_626_577_4e-103,
        ],
        rtol = 1e-15
    );

    // The same point at scale 2^427.
    let x = 2.0_f64.powi(427);
    let z = c64(5_170_128_475_599_457.0 * x, 2_984_975_067_132_296.0 * x);
    let (_, dq) = xsf::assoc_legendre_q_all(1, 100, z);
    assert_eq!(dq[1][100].re, -1.723_798_258_841_282_8e-305);
}

#[test]
fn test_assoc_legendre_q_all_extreme_cancelled_leading_part() {
    // Re z^-28 cancels beyond double-double precision. Arb: +3.5e-4090, which underflows to +0.
    let z = c64(4.271_237_208_723_019e144, 7.257_123_905_997_274e143);
    let (q, _) = xsf::assoc_legendre_q_all(27, 1, z);
    assert_eq!(q[27][1].re.to_bits(), 0);
}

#[test]
fn test_assoc_legendre_q_all_extreme_z_squared_overflow() {
    let (q, dq) = xsf::assoc_legendre_q_all(1, 100, c64(-6e199, 8e199));
    xsref::np_assert_allclose!(
        [q[0][5].re, q[0][5].im, q[1][100].re, q[1][100].im, dq[0][100].re, dq[0][100].im],
        [
            -7.2e-199,
            -9.6e-199,
            -8.797_551_242_515_803e-242,
            3.016_303_283_148_274_5e-241,
            2.613_134_032_430_436_5e-243,
            -8.959_316_682_618_637e-243,
        ],
        rtol = 1e-15
    );
}

#[test]
fn test_assoc_legendre_q_all_large_x_and_z() {
    // where the backward recurrences are rescaled; mpmath: legenq(n, m, x, type=3), and Arb for
    // complex z
    for (x, n, m, q) in [
        (1e7, 0, 0, 1.000_000_000_000_003_4e-7),
        (1e8, 0, 0, 1e-8),
        (907.57, 30, 26, 8.065_789_716_549_577e-60),
    ] {
        let (qm, _) = xsf::assoc_legendre_q_all(n, m, x);
        xsref::np_assert_allclose!([qm[n][m]], [q], rtol = 5e-15);
    }
    let (qm, _) = xsf::assoc_legendre_q_all(3, 2, c64(1e8, 3e7));
    xsref::np_assert_allclose!(
        [qm[3][2]],
        [c64(3.789_872_471_995_64e-33, -8.841_146_634_093_655e-33)],
        rtol = 5e-15
    );
}

#[test]
fn test_assoc_legendre_q_all_exterior_forward_seed() {
    // Arb at the exact binary64 inputs, with the principal sqrt(z² - 1).
    for re in [-0.6_f64, 0.6] {
        let z = c64(re, 0.800_000_001);
        let (q, dq) = xsf::assoc_legendre_q_all(3, 1, z);
        let expected = c64(0.025_564_854_500_179_927_f64.copysign(re), -0.132_252_482_592_822_88);
        let expected_dq = c64(0.364_058_770_999_769_7, 0.185_370_137_953_699_2_f64.copysign(re));
        xsref::np_assert_allclose!([q[3][1], dq[3][1]], [expected, expected_dq], rtol = 2e-14);
    }
}

#[test]
fn test_assoc_legendre_q_all_unit_circle_derivative() {
    for z in [I, -I] {
        let (q, dq) = xsf::assoc_legendre_q_all(0, 1, z);
        assert_eq!(dq[0][0], c64(0.5, 0.0));
        // Both branch conventions give the same magnitude for Q_0^1.
        xsref::np_assert_allclose!([q[0][1].norm()], [consts::FRAC_1_SQRT_2], rtol = 1e-15);
    }

    // The norm rounds to one on both sides of the exact unit circle.
    for (x, y) in [(0.6, 0.8), (0.28, 0.96)] {
        for re in [-x, x] {
            for im in [-y, y] {
                let z = c64(re, im);
                let (_, dq) = xsf::assoc_legendre_q_all(0, 0, z);
                xsref::np_assert_allclose!([dq[0][0]], [1.0 / (1.0 - z * z)], rtol = 1e-15);
            }
        }
    }
}

#[test]
fn test_assoc_legendre_q_all_q0_tiny_and_large_z() {
    // Arb: atanh(z), and acoth(z) for |z| > 1
    for (z, q0) in [
        (c64(1e-9, 2e-9), c64(1e-9, 2e-9)),
        (c64(-3e-7, 1e-7), c64(-3.000_000_000_000_059_7e-7, 1.000_000_000_000_086_6e-7)),
        (c64(1e8, 3e7), c64(9.174_311_926_605_505e-9, -2.752_293_577_981_651_6e-9)),
        (c64(0.0, 1e-10), c64(0.0, 1e-10)),
    ] {
        let (qm, _) = xsf::assoc_legendre_q_all(0, 0, z);
        xsref::np_assert_allclose!([qm[0][0]], [q0], rtol = 1e-15);
    }
}

#[test]
fn test_legendre_q_all_c64_inf() {
    // `Q_0(z) -> 0` as `|z| -> inf`
    let inf = f64::INFINITY;
    for z in [c64(inf, 0.0), c64(-inf, 1.0), c64(0.5, inf), c64(inf, -inf)] {
        assert_eq!(xsf::legendre_q_all(1, z).0[0], c64(0.0, 0.0), "{z}");
        assert_eq!(xsf::assoc_legendre_q_all(1, 1, z).0[0][0], c64(0.0, 0.0), "{z}");
    }
    // NaN wins over infinity
    let z = c64(inf, f64::NAN);
    assert!(xsf::legendre_q_all(1, z).0[0].is_nan());
    assert!(xsf::assoc_legendre_q_all(1, 1, z).0[0][0].is_nan());
}

#[test]
fn test_complex_legendre_q_exterior_near_unit_circle() {
    let z = c64(0.154_255_747_404_549_34, -0.988_032_379_898_993_6);
    let (q, d) = xsf::legendre_q_all(30, z);
    let z = c64(0.702_729_738_621_847_4, 0.711_585_757_122_401_1);
    let (aq, ad) = xsf::assoc_legendre_q_all(13, 7, z);
    let (inner, _) = xsf::assoc_legendre_q_all(32, 32, c64(0.0, 1.00005));
    let (outer, _) = xsf::assoc_legendre_q_all(32, 32, c64(0.0, 1.02));
    // Arb values, derivatives, and the oblate radial ratio Q_n^m(i u/E) / Q_n^m(i b/E).
    xsref::np_assert_allclose!(
        [q[30], d[30], aq[13][7], ad[13][7], outer[32][32] / inner[32][32]],
        [
            c64(1.096_037_748_414_413_7e-13, 4.462_716_575_483_255e-13),
            c64(9.583_271_348_237_42e-12, -3.148_034_361_847_163_4e-12),
            c64(1_425.676_174_740_597_4, 2_272.285_117_041_535_6),
            c64(-28_530.805_247_220_02, 7_356.536_755_705_318),
            c64(0.712_916_713_273_161_5, 0.0),
        ],
        rtol = 1e-13
    );
}

#[test]
fn test_assoc_legendre_q_all_interior_diagonals() {
    // Arb type 2 at the exact binary64 argument.
    let z = c64(0.0, 0.95);
    let (q, d) = xsf::assoc_legendre_q_all(41, 40, z);
    xsref::np_assert_allclose!(
        [q[40][40], d[40][40], q[41][40], d[41][40]],
        [
            c64(0.0, 4.836_108_252_798_383e64),
            c64(9.659_506_628_466_28e65, 0.0),
            c64(-3.721_385_300_530_679_5e66, 0.0),
            c64(0.0, 7.824_715_119_076_829e67),
        ],
        rtol = 2e-14
    );

    let (prefix_q, prefix_d) = xsf::assoc_legendre_q_all(40, 40, z);
    assert_eq!((prefix_q[40][40], prefix_d[40][40]), (q[40][40], d[40][40]));
    let (q_conj, d_conj) = xsf::assoc_legendre_q_all(41, 40, z.conj());
    for n in [40, 41] {
        assert_eq!((q_conj[n][40], d_conj[n][40]), (q[n][40].conj(), d[n][40].conj()));
    }

    // The unscaled order recurrence overflows before this finite diagonal.
    let (q, d) = xsf::assoc_legendre_q_all(133, 133, z);
    xsref::np_assert_allclose!(
        [q[133][133], d[133][133]],
        [c64(0.0, -4.673_913_265_621_371e283), c64(-3.104_068_021_609_778e285, 0.0)],
        rtol = 3e-14
    );

    // Arb, with the real part of z subnormal.
    let x = 1e-321;
    let (_, d) = xsf::assoc_legendre_q_all(2, 2, c64(x, 0.1));
    assert_eq!(d[2][2].im, -3e-323);
}

#[test]
fn test_legendre_q_all_c64_off_axis() {
    // NaN in xsf: `km` is not positive for `z = -0.389 - 1.0005i`, and the backward recurrence
    // overflows for `z = 1e8 + 3e7i`; Arb references
    let (q, d) = xsf::legendre_q_all(2, c64(-0.389, -1.0005));
    xsref::np_assert_allclose!(
        [q[0], q[1], q[2], d[0]],
        [
            c64(-0.189_288_021_445_159_6, 0.747_410_144_189_767_6),
            c64(-0.178_583_110_395_970_47, -0.101_359_880_633_937_43),
            c64(0.046_731_414_777_246_995, -0.046_552_978_818_228_626),
            c64(0.459_296_280_510_478_77, 0.193_282_793_484_476_33),
        ],
        rtol = 1e-13
    );
    let (q, d) = xsf::legendre_q_all(2, c64(1e8, 3e7));
    xsref::np_assert_allclose!(
        [q[0], q[1], q[2], d[0]],
        [
            c64(9.174_311_926_605_505e-9, -2.752_293_577_981_651_6e-9),
            c64(2.553_095_979_575_232e-17, -1.683_359_986_533_12e-17),
            c64(7.515_919_205_927_693e-26, -8.988_215_707_910_789e-26),
            c64(-7.659_287_938_725_697e-17, 5.050_079_959_599_361e-17),
        ],
        rtol = 1e-14
    );
    // finite before, but 4.7e-14 off in `Q_0`, as the recurrence overflowed; mpmath references
    let (q, d) = xsf::legendre_q_all(7, c64(1e6, 3e5));
    xsref::np_assert_allclose!(
        [q[0], q[7], d[0], d[7]],
        [
            c64(9.174_311_926_607_384e-7, -2.752_293_577_983_898_4e-7),
            c64(-1.214_582_083_191_448e-51, -1.275_706_366_917_626_5e-51),
            c64(-7.659_287_938_729_013e-13, 5.050_079_959_607_096e-13),
            c64(1.172_325_866_618_442e-56, 6.688_673_335_487_226e-57),
        ],
        rtol = 1e-14
    );
}

#[test]
fn test_legendre_q_all_c64_tiny_im() {
    // the imaginary parts of the seeds were subnormal, so these are compared by component; Arb
    // references
    let parts = |q: num_complex::Complex64, d: num_complex::Complex64| [q.re, q.im, d.re, d.im];

    let (cqm, cqd) = xsf::assoc_legendre_q_all(1, 140, c64(10000.0, 1e-304));
    xsref::np_assert_allclose!(
        parts(cqm[1][140], cqd[1][140]),
        [
            6.327_269_914_478_720_4e234,
            -1.265_478_793_247_977_9e-73,
            -1.265_478_793_247_978e231,
            3.796_535_621_848_62e-77,
        ],
        rtol = 1e-12
    );
    let (cqm, cqd) = xsf::assoc_legendre_q_all(60, 120, c64(10000.0, 1e-100));
    xsref::np_assert_allclose!(
        parts(cqm[60][120], cqd[60][120]),
        [
            2.381_042_336_285_591_4e-16,
            -1.452_435_860_331_049_8e-118,
            -1.452_435_860_331_049_7e-18,
            9.005_102_555_792_6e-121,
        ],
        rtol = 1e-13
    );
    let (cqm, cqd) = xsf::assoc_legendre_q_all(40, 40, c64(3.0, 1e-280));
    xsref::np_assert_allclose!(
        parts(cqm[40][40], cqd[40][40]),
        [
            3.206_432_995_060_525e39,
            -4.916_818_632_369_847e-240,
            -4.916_818_632_369_847e40,
            7.743_557_286_066_945e-239,
        ],
        rtol = 1e-14
    );
    let (cqn, cqd) = xsf::legendre_q_all(1, c64(10.0, 1e-280));
    xsref::np_assert_allclose!(
        parts(cqn[1], cqd[1]),
        [
            0.003_353_477_310_755_806,
            -6.747_532_790_254_294e-284,
            -0.000_674_753_279_025_429_4,
            2.040_608_101_214_161_7e-284,
        ],
        rtol = 1e-14
    );
}

#[test]
fn test_legendre_q_all_c64_tiny_im_real_parts() {
    // the real parts are those at `z` itself, which keep `Re z = 1e-40` and their values near the
    // zeros of `Re Q_1` and `Re Q_2'`; Arb references (Ferrers), by part
    let z = c64(1e-40, 1e-100);
    let expected = [
        [1e-40, 1e-100, 1.0, 2e-140],
        [-1.0, 2e-140, 2e-40, 2e-100],
        [-2e-40, -2e-100, -2.0, 8e-140],
    ];
    let (cqn, cqd) = xsf::legendre_q_all(2, z);
    let (cqm, cqmd) = xsf::assoc_legendre_q_all(2, 0, z);
    for (j, e) in expected.iter().enumerate() {
        let (q, d) = (cqn[j], cqd[j]);
        xsref::np_assert_allclose!([q.re, q.im, d.re, d.im], e, rtol = 1e-15);
        assert_eq!((cqm[j][0], cqmd[j][0]), (q, d));
    }

    // these cancel in the recurrence (as do `Im Q_2` and `Im Q_2'` near 0.6232), so all parts
    // are only accurate relative to the largest real or imaginary part
    for (x, re, im) in [
        (
            0.833_556_559_600_964_8,
            [
                1.199_678_640_257_734,
                2.131_319_626_955_773_2e-16,
                -0.599_839_320_128_866_8,
                3.276_717_531_228_073_6,
                3.931_008_032_372_374,
                3.276_717_531_228_074,
            ],
            [
                3.276_717_531_228_073_7e-100,
                3.931_008_032_372_374e-100,
                3.276_717_531_228_074e-100,
                1.789_958_980_540_111_4e-99,
                2.147_375_555_891_480_2e-99,
                2.969_261_390_251_824e-99,
            ],
        ),
        (
            0.623_174_709_782_724_5,
            [
                0.730_178_774_530_167_4,
                -0.544_971_054_092_657_4,
                -0.874_507_654_976_349_6,
                1.634_913_162_277_973,
                1.749_015_309_952_699_5,
                4.725_715_661_779_876e-16,
            ],
            [
                1.634_913_162_277_973_1e-100,
                1.749_015_309_952_699_6e-100,
                4.725_715_661_779_876e-116,
                3.331_418_523_943_972e-100,
                5.345_882_096_379_523e-100,
                8.578_464_453_802_071e-100,
            ],
        ),
    ] {
        let (cqn, cqd) = xsf::legendre_q_all(2, c64(x, 1e-100));
        let parts = |f: fn(&num_complex::Complex64) -> f64| {
            cqn.iter().chain(&cqd).map(f).collect::<Vec<_>>()
        };
        xsref::np_assert_allclose!(parts(|q| q.re), re, rtol = 1e-15, atol = 4.0 * f64::EPSILON);
        xsref::np_assert_allclose!(parts(|q| q.im), im, rtol = 1e-15, atol = 1e-99 * f64::EPSILON);
    }
}

#[test]
fn test_legendre_q_all_c64_tiny_re() {
    // the quotient in `Q_0` loses the digits of `Re z`, so these are compared by component; Arb
    // references (mpmath's `legenq` type 2 for `|z| < 1`, else 3)
    let parts = |q: &[num_complex::Complex64]| -> Vec<f64> {
        q.iter().flat_map(|c| [c.re, c.im]).collect()
    };
    let z = c64(1e-10, 1.5);
    let (q, d) = xsf::legendre_q_all(2, z);
    let (qm, _) = xsf::assoc_legendre_q_all(2, 0, z);
    let expected = [
        3.076_923_076_923_077e-11,
        -0.588_002_603_547_567_5,
        -0.117_996_094_678_648_67,
        -1.264_641_420_091_060_2e-11,
        -4.629_597_634_363_833e-12,
        0.028_510_088_746_824_26,
    ];
    xsref::np_assert_allclose!(parts(&q), expected, rtol = 1e-14);
    xsref::np_assert_allclose!(parts(&[qm[0][0], qm[1][0], qm[2][0]]), expected, rtol = 1e-14);
    xsref::np_assert_allclose!(
        parts(&d),
        [
            0.307_692_307_692_307_7,
            2.840_236_686_390_532_7e-11,
            1.893_491_124_260_355e-11,
            -0.126_464_142_009_106_02,
            -0.046_295_976_343_638_33,
            -9.536_875_738_826_48e-12,
        ],
        rtol = 1e-14
    );

    let (q, _) = xsf::legendre_q_all(3, c64(1e-12, -0.3));
    xsref::np_assert_allclose!(
        parts(&[q[0], q[3]]),
        [
            9.174_311_926_605_504e-13,
            -0.291_456_794_477_867_1,
            1.042_495_557_808_963,
            2.608_689_170_191_196e-12,
        ],
        rtol = 1e-14
    );
    let (q, _) = xsf::legendre_q_all(3, c64(-3e-9, 50.0));
    xsref::np_assert_allclose!(
        parts(&[q[0], q[2], q[3]]),
        [
            -1.199_520_191_923_230_7e-12,
            -0.019_997_333_973_150_53,
            1.918_903_368_919_533_4e-16,
            1.066_301_074_244_355e-6,
            9.138_795_201_736_145e-9,
            -2.192_823_601_526_560_7e-18,
        ],
        rtol = 1e-14
    );
}

#[test]
fn test_legendre_q_all_c64_tiny_re_backward_rescaling() {
    // Compare the small parts separately against Arb on both sides of the imaginary axis.
    let values = [
        9.999_000_099_990_002e-255,
        -6.665_866_752_372_065e-257,
        -3.999_428_638_087_966e-259,
        2.285_333_381_812_588_3e-261,
    ];
    let derivatives = [
        1.999_600_059_992_001e-256,
        1.999_600_059_992_001e-258,
        -1.599_657_196_183_204_4e-260,
        -1.142_590_519_820_686_9e-262,
    ];
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            let (q, dq) = xsf::legendre_q_all(3, c64(sx * 1e-250, sy * 100.0));
            xsref::np_assert_allclose!(
                [q[0].re, q[1].im, q[2].re, q[3].im],
                [sx * values[0], sx * sy * values[1], sx * values[2], sx * sy * values[3]],
                rtol = 1e-14
            );
            xsref::np_assert_allclose!(
                [dq[0].im, dq[1].re, dq[2].im, dq[3].re],
                [
                    sx * sy * derivatives[0],
                    sx * derivatives[1],
                    sx * sy * derivatives[2],
                    sx * derivatives[3],
                ],
                rtol = 1e-14
            );
        }
    }
}

#[test]
fn test_legendre_q_all_c64_tiny_re_subnormal() {
    // The first-order terms preserve subnormal parts at higher degrees (Arb).
    let (q, dq) = xsf::legendre_q_all(95, c64(1e-280, 1.125));
    xsref::np_assert_allclose!(
        [q[94].re, q[95].im, dq[94].im, dq[95].re],
        [-1.3544e-319, 5.1763e-320, -8.603_877e-318, -3.32274e-318],
        rtol = 0.0,
        atol = DENORM_MIN
    );

    // At this degree, rounding the small parts of Q before differentiating loses digits.
    let (q, dq) = xsf::legendre_q_all(750, c64(5.9e-15, 1.125));
    assert_eq!(q[739].im.to_bits(), 1);
    assert_eq!(dq[739].re.to_bits(), (-(f64::from_bits(284))).to_bits());
}

#[test]
fn test_legendre_q_all_c64_tiny_relative_re() {
    let z = c64(2.0_f64.powi(-100), 2.0_f64.powi(38));
    let (q, _) = xsf::legendre_q_all(3, z);
    xsref::np_assert_allclose!([q[0].re], [1.044_048_714_879_764e-53], rtol = 1e-15);

    let (_, dq) = xsf::legendre_q_all(1, c64(-1e-250, 1e155));
    assert_eq!(dq[0].im.to_bits(), (-0.0_f64).to_bits());
}

#[test]
fn test_assoc_legendre_q_all_c64_tiny_im_axis() {
    // on and near the imaginary axis; Arb references (Ferrers), with
    // `Q_20^20(z) = 2^20 20! z + O(z^3)`
    for (n, m, z, q, d) in [
        (20, 20, c64(0.0, 1e-320), 2.551_054_255_388_657e-296, 2.551_082_656_125_828_5e24),
        (20, 20, c64(1e-300, 1e-320), 2.551_054_255_388_657e-296, 2.551_082_656_125_828_5e24),
        (20, 20, c64(0.0, 1e-200), 2.551_082_656_125_828_5e-176, 2.551_082_656_125_828_5e24),
    ] {
        let (cqm, cqd) = xsf::assoc_legendre_q_all(n, m, z);
        xsref::np_assert_allclose!([cqm[n][m].im, cqd[n][m].re], [q, d], rtol = 1e-15);
    }
    let (cqm, cqd) = xsf::assoc_legendre_q_all(3, 2, c64(0.0, 1e-320));
    assert_eq!((cqm[3][2], cqd[3][2]), (c64(-8.0, 0.0), c64(0.0, 6.399_93e-319)));
}

#[test]
fn test_assoc_legendre_q_all_c64_tiny_im_lift_overflow() {
    // the lifted run overflows here, so this imaginary part is kept; Arb reference (Ferrers)
    let d = xsf::assoc_legendre_q_all(0, 33, c64(0.999_999_999_999_999_9, 1e-320)).1[0][33].im;
    xsref::np_assert_allclose!([d], [-5.088_204_835_535_055e17], rtol = 1e-14);
    // NaN without the lift, but about -5.6e650 (Arb)
    let q = xsf::assoc_legendre_q_all(0, 100, c64(1.000_000_000_000_000_2, 1e-320)).0[0][100].im;
    assert_eq!(q, f64::NEG_INFINITY);
}

#[test]
fn test_assoc_legendre_q_all_derivative_underflow_sign() {
    let x = 1.560_874_275_157_996e144_f64;
    let (_, positive) = xsf::assoc_legendre_q_all(11, 1, x);
    let (_, negative) = xsf::assoc_legendre_q_all(11, 1, -x);
    for n in [4, 5, 6, 10, 11] {
        let expected = if n % 2 == 0 { -0.0_f64 } else { 0.0 };
        assert_eq!(negative[n][0].to_bits(), expected.to_bits());
        assert_eq!(negative[n][1].to_bits(), expected.to_bits());
        assert_eq!(positive[n][0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(positive[n][1].to_bits(), 0.0_f64.to_bits());
    }
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
            // not bit-identical, since only the real `Q_0` uses `ln_1p`
            assert!((q - cq.re).abs() <= 1e-15 * q.abs(), "{q} != {}", cq.re);
        }
    }
    // mpmath: legenq(2, 1, -1.05, type=3) = 2.211389887900321167... (with the opposite sign)
    let q = xsf::assoc_legendre_q_all(2, 1, -1.05).0[2][1];
    assert!((q + 2.211_389_887_900_321_2).abs() < 1e-15, "{q}");
}

#[test]
fn test_assoc_legendre_q_all_extreme_real_x() {
    for x in [1.0e160_f64, -1.0e160] {
        let (q, dq) = xsf::assoc_legendre_q_all(1, 1, x);
        assert_eq!(q[0][0], x.recip());
        let inv = x.abs().recip();
        assert_eq!(q[0][1], -inv);
        assert_eq!(q[1][0], inv / (3.0 * x.abs()));
        assert_eq!(dq[0][0], -inv / x.abs());
        assert_eq!(dq[0][1], inv / x);
    }
    // At x = 2^538 the leading derivative magnitude 2/x² is 2^-1075,
    // halfway to the least subnormal. The exact derivative has greater magnitude.
    let x = 2.0_f64.powi(538);
    let (_, dq) = xsf::assoc_legendre_q_all(0, 2, x);
    assert_eq!(dq[0][2], -DENORM_MIN);
}

#[test]
fn test_assoc_legendre_q_all_extreme_complex_axis() {
    let square_overflow = 2.0_f64.powi(512);
    let below = f64::from_bits(square_overflow.to_bits() - 1);
    let above = f64::from_bits(square_overflow.to_bits() + 1);
    for sign in [-1.0, 1.0] {
        for lip in [0.0_f64, -0.0] {
            let (q, _) = xsf::assoc_legendre_q_all(1, 0, c64(sign * below, lip));
            // Q_1 is finite just before the square overflows.
            let expected = (1.0 / below) / (3.0 * below);
            assert!((q[1][0].re - expected).abs() <= DENORM_MIN);
            for x in [square_overflow, above, 1e160] {
                let (q, d) = xsf::assoc_legendre_q_all(2, 0, c64(sign * x, lip));
                let inverse = 1.0 / x;
                assert_eq!(q[1][0].re, inverse / (3.0 * x));
                assert_eq!(q[1][0].im.to_bits(), (-sign * lip).to_bits());
                assert_eq!(d[0][0].re, -inverse / x);
                assert_eq!(d[0][0].im.to_bits(), (sign * lip).to_bits());
                // The leading Q_2 and Q_1' terms scale as x^-3.
                assert_eq!(q[2][0].re.to_bits(), (sign * 0.0).to_bits());
                assert_eq!(d[1][0].re.to_bits(), (-sign * 0.0).to_bits());
                assert_eq!(d[1][0].im.to_bits(), lip.to_bits());
            }
        }
    }
}

#[test]
fn test_assoc_legendre_q_all_large_x() {
    // NaN in xsf, since the backward recurrences in `lqmn` overflow; Arb references
    let (qm, qd) = xsf::assoc_legendre_q_all(30, 26, 907.57);
    xsref::np_assert_allclose!(
        [qm[30][26], qd[30][26]],
        [8.065_789_716_549_577e-60, -2.755_046_337_773_331_8e-61],
        rtol = 1e-14
    );
    let (cqm, cqd) = xsf::assoc_legendre_q_all(30, 26, c64(600.0, 400.0));
    xsref::np_assert_allclose!(
        [cqm[30][26], cqd[30][26]],
        [
            c64(8.184_368_733_381_842e-57, 5.860_850_039_088_338e-57),
            c64(-4.325_076_007_856_238e-58, -1.447_173_501_838_919e-59),
        ],
        rtol = 1e-13
    );
    // the orders 0 and 1 underflow, but the order recurrence raises them back into range
    let (qm, qd) = xsf::assoc_legendre_q_all(31, 26, 1e10);
    xsref::np_assert_allclose!(
        [qm[31][26], qd[31][26]],
        [3.609_593_573_829_175e-288, -1.155_069_943_625_335_8e-296],
        rtol = 1e-14
    );
    let (cqm, cqd) = xsf::assoc_legendre_q_all(110, 80, c64(1000.0, 500.0));
    xsref::np_assert_allclose!(
        [cqm[110][80], cqd[110][80]],
        [
            c64(6.002_372_596_505_569e-200, -1.541_295_560_633_932_2e-199),
            c64(1.513_253_937_967_822_2e-201, 1.635_176_324_583_510_4e-200),
        ],
        rtol = 1e-13
    );
    // the lifted orders grow past the overflow threshold before they are scaled back
    let (qm, qd) = xsf::assoc_legendre_q_all(400, 250, 10.0);
    xsref::np_assert_allclose!(
        [qm[400][250], qd[400][250]],
        [8.776_894_642_465_389e156, -3.544_215_960_415_300_4e158],
        rtol = 1e-13
    );
    let (cqm, cqd) = xsf::assoc_legendre_q_all(400, 250, c64(10.0, 1.0));
    xsref::np_assert_allclose!(
        [cqm[400][250], cqd[400][250]],
        [
            c64(-9.462_528_891_795_308e155, -6.445_726_315_936_386e155),
            c64(4.043_140_341_411_489e157, 2.192_738_476_456_626_5e157),
        ],
        rtol = 1e-13
    );
    // the products overflow before they cancel, after the lift is undone; mpmath references
    let (qm, qd) = xsf::assoc_legendre_q_all(111, 315, 4960.0);
    xsref::np_assert_allclose!(
        [qm[111][315], qd[111][315]],
        [-1.591_573_603_547_800_6e308, 3.593_876_527_546_9e306],
        rtol = 1e-13
    );
    let (cqm, cqd) = xsf::assoc_legendre_q_all(111, 315, c64(4960.0, 1.0));
    xsref::np_assert_allclose!(
        [cqm[111][315], cqd[111][315]],
        [
            c64(-1.591_164_238_540_185e308, 3.593_562_895_647_709_5e306),
            c64(3.592_935_648_855_753_2e306, -8.186_936_570_392_784e304),
        ],
        rtol = 1e-13
    );
}
