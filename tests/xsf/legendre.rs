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

    xsref::np_assert_allclose!(
        &qn[..2],
        &[LN_3 * 0.5, LN_3 * 0.25 - 1.0],
        rtol = f64::EPSILON
    );
    xsref::np_assert_allclose!(
        &qd[..2],
        &[4.0 / 3.0, 2.0 / 3.0 + 0.5 * LN_3],
        rtol = f64::EPSILON
    );
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
fn test_assoc_legendre_q_all_q0_tiny_and_large_z() {
    // Arb: atanh(z), and acoth(z) for |z| > 1
    for (z, q0) in [
        (c64(1e-9, 2e-9), c64(1e-9, 2e-9)),
        (
            c64(-3e-7, 1e-7),
            c64(-3.000_000_000_000_059_7e-7, 1.000_000_000_000_086_6e-7),
        ),
        (
            c64(1e8, 3e7),
            c64(9.174_311_926_605_505e-9, -2.752_293_577_981_651_6e-9),
        ),
        (c64(0.0, 1e-10), c64(0.0, 1e-10)),
    ] {
        let (qm, _) = xsf::assoc_legendre_q_all(0, 0, z);
        xsref::np_assert_allclose!([qm[0][0]], [q0], rtol = 1e-15);
    }
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
