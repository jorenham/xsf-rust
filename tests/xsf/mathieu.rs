use core::f64::consts::FRAC_1_SQRT_2;

xsref::tests! {
    test_mathieu_a: "cem_cva", "d_d-d", |x| xsref::with_int(x[0], |m| xsf::mathieu_a(m, x[1]));
    test_mathieu_b: "sem_cva", "d_d-d", |x| xsref::with_int(x[0], |m| xsf::mathieu_b(m, x[1]));
    test_mathieu_cem: "cem", "d_d_d-d_d", |x| {
        xsref::with_int(x[0], |m| xsf::mathieu_cem(m, x[1], x[2]))
    };
    test_mathieu_sem: "sem", "d_d_d-d_d", |x| {
        xsref::with_int(x[0], |m| xsf::mathieu_sem(m, x[1], x[2]))
    };
}

#[test]
fn test_mathieu_sem_even_m() {
    // only accurate to ~2e-9 in xsf (scipy/xsf#343); the expected values are from mpmath
    let (y, yp) = xsf::mathieu_sem(50, -1e-5, 143.239_448_782_705_8);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.616_040_534_405_098_7, 39.385_722_677_066_29],
        rtol = 1e-13
    );

    // zeros in xsf, because of an overflow in `fcoef`
    let (y, yp) = xsf::mathieu_sem(96, 0.01, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.866_016_937_460_664, -48.001_456_348_112_84],
        rtol = 1e-13
    );
}

#[test]
fn test_mathieu_small_q_low_modes() {
    // References from the normalized Fourier eigenvector.
    let at_zero = xsf::mathieu_cem(2, 1e-3, 0.0).0;
    let derivative = xsf::mathieu_cem(2, 1e-3, 22.5).1;
    let cv = xsf::mathieu_a(2, 1e-3);
    xsref::np_assert_allclose!(
        [at_zero, derivative, cv],
        [1.000_166_603_258_439_3, -1.413_880_146_799_069_4, 4.000_000_416_666_611],
        rtol = 2e-15
    );
}

#[test]
fn test_mathieu_cem_sem_large_m() {
    // zeros in xsf, because of an overflow in `fcoef` (scipy/xsf#344); the expected values
    // are from mpmath
    let (y, yp) = xsf::mathieu_cem(92, 0.01, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.939_699_500_563_974_8, 31.464_228_932_692_56],
        rtol = 1e-13
    );
    let (y, yp) = xsf::mathieu_cem(93, 0.01, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.866_035_069_310_786_4, 46.498_493_909_628_46],
        rtol = 1e-13
    );
    let (y, yp) = xsf::mathieu_sem(93, 0.01, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.499_984_345_080_641_2, -80.541_173_731_920_72],
        rtol = 1e-13
    );

    // NaN in xsf, because the coefficients in `fcoef` themselves overflow; the expected values
    // are from mpmath
    let (y, yp) = xsf::mathieu_cem(200, 1.0, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.939_995_772_283_597_3, 68.242_508_215_617_67],
        rtol = 1e-13
    );
    let (y, yp) = xsf::mathieu_cem(161, 1e-3, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.984_807_586_411_290_2, -27.957_524_503_128_47],
        rtol = 1e-13
    );
    let (y, yp) = xsf::mathieu_sem(100, 1e-6, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.984_807_753_355_469_1, 17.364_817_597_481_87],
        rtol = 1e-13
    );
    let (y, yp) = xsf::mathieu_sem(351, 300.0, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-0.990_457_823_187_051_1, -51.090_578_748_552_91],
        rtol = 1e-13
    );
}

#[test]
fn test_mathieu_cem_sem_km_gt_251() {
    // NaN in xsf, which only has room for 251 coefficients; the expected values are from the
    // tridiagonal eigenproblem (DLMF 28.4.5) in mpmath
    let (y, yp) = xsf::mathieu_cem(1000, 1.0, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [0.173_479_844_385_168_09, 984.836_971_453_448_2],
        rtol = 1e-12
    );
    let (y, yp) = xsf::mathieu_sem(10_000, 100.0, 37.3);
    xsref::np_assert_allclose!(
        [y, yp],
        [0.639_087_540_712_494, 7_691.339_718_009_938],
        rtol = 1e-11
    );
    // also needs the backward sum of squares to leave room for the forward one
    let (y, yp) = xsf::mathieu_cem(131, 2100.0, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [0.400_830_572_645_642_4, -113.482_637_316_450_79],
        rtol = 1e-12
    );
}

#[test]
fn test_mathieu_cem_large_q() {
    // zeros in xsf, because the backward sum of squares in `fcoef` overflows; the expected values
    // are from the tridiagonal eigenproblem (DLMF 28.4.5) in mpmath
    let (y, yp) = xsf::mathieu_cem(170, 1500.0, 10.0);
    xsref::np_assert_allclose!(
        [y, yp],
        [-1.015_989_436_317_424_7, -21.795_261_688_798_792],
        rtol = 1e-13
    );
}

#[test]
fn test_mathieu_a_b_huge_m() {
    // 0 resp. NaN in xsf, because `m * m` overflows `int`; the second expected value is from
    // the series DLMF 28.6.14
    let m = 1 << 30;
    assert_eq!(xsf::mathieu_a(m, 0.0), 2f64.powi(60));
    assert_eq!(xsf::mathieu_b(m, 0.0), 2f64.powi(60));
    xsref::np_assert_allclose!(
        [xsf::mathieu_a(46341, 1e6)],
        [2_147_488_513.830_157_2],
        rtol = 1e-15
    );

    // `m` doesn't fit in `int`
    assert!(xsf::mathieu_a(u32::MAX, 0.0).is_nan());
    assert!(xsf::mathieu_b(u32::MAX, 0.0).is_nan());
}

xsref::tests! {
    test_mathieu_modcem1: "mcm1", "d_d_d-d_d", |x| {
        xsref::with_int(x[0], |m| xsf::mathieu_modcem1(m, x[1], x[2]))
    };
    test_mathieu_modsem1: "msm1", "d_d_d-d_d", |x| {
        xsref::with_int(x[0], |m| xsf::mathieu_modsem1(m, x[1], x[2]))
    };
    test_mathieu_modcem2: "mcm2", "d_d_d-d_d", |x| {
        xsref::with_int(x[0], |m| xsf::mathieu_modcem2(m, x[1], x[2]))
    };
    test_mathieu_modsem2: "msm2", "d_d_d-d_d", |x| {
        xsref::with_int(x[0], |m| xsf::mathieu_modsem2(m, x[1], x[2]))
    };
}

#[test]
fn test_mathieu_a_nan() {
    assert!(xsf::mathieu_a(13, f64::NAN).is_nan());
}

#[test]
fn test_mathieu_mod_first_kind_near_origin() {
    for sign in [-1.0, 0.0, 1.0] {
        let x = sign * 1.0e-9;
        let (mc, dmc) = xsf::mathieu_modcem1(9, 1.0, x);
        let (ms, dms) = xsf::mathieu_modsem1(9, 1.0, x);
        xsref::np_assert_allclose!(
            &[mc, dmc, ms, dms],
            &[
                1.083_083_254_725_196_7e-8,
                sign * 8.557_034_674_740_288e-16,
                sign * 9.625_470_732_010_198e-17,
                9.625_470_732_010_198e-8,
            ],
            rtol = 4e-15
        );
    }
}

#[test]
fn test_mathieu_modcem1_near_linearized_coefficients() {
    let (_, dy) = xsf::mathieu_modcem1(1, 1.001_442_352_383_129_1e-7, f64::from_bits(1580));
    assert_eq!(dy.to_bits(), 0);
    let (_, dy) = xsf::mathieu_modcem1(1, f64::from_bits(1.0e-7_f64.to_bits() + 1), 0.01);
    xsref::np_assert_allclose!(&[dy], &[3.162_329_930_191_721e-6], rtol = 2e-14);
}

#[test]
fn test_mathieu_modcem1_large_x() {
    // `(f^2 + (f'/u)^2) * u -> 2/pi` with `u = sqrt(q) * exp(x)`, i.e. `|Mc1| ~ sqrt(2 / (pi u))`
    let x = 22.0;
    let u = f64::exp(x);
    let (f, d) = xsf::mathieu_modcem1(0, 1.0, x);
    let env = (f * f + (d / u) * (d / u)) * u;
    assert!((env * core::f64::consts::FRAC_PI_2 - 1.0).abs() < 1e-6, "{env}");
}

#[test]
fn test_mathieu_coef() {
    // q = 0: ce_0 = 1 / sqrt(2) and se_1 = sin(z) (DLMF 28.4.13 and 28.4.15)
    for (coefs, first) in
        [(xsf::mathieu_even_coef(0, 0.0), FRAC_1_SQRT_2), (xsf::mathieu_odd_coef(1, 0.0), 1.0)]
    {
        let coefs = coefs.unwrap();
        xsref::np_assert_allclose!(&[coefs[0]], &[first], rtol = 1e-15);
        assert!(coefs[1..].iter().all(|&c| c == 0.0), "{coefs:?}");
    }

    // q = 1: the normalizations of DLMF 28.4.9 and 28.4.11
    let a = xsf::mathieu_even_coef(2, 1.0).unwrap();
    let b = xsf::mathieu_odd_coef(1, 1.0).unwrap();
    let norm_a = a[0] * a[0] + a.iter().map(|c| c * c).sum::<f64>();
    let norm_b = b.iter().map(|c| c * c).sum::<f64>();
    xsref::np_assert_allclose!(&[norm_a, norm_b], &[1.0, 1.0], rtol = 1e-15);
}

#[test]
fn test_mathieu_even_coef_match_underflow() {
    let coefficients = xsf::mathieu_even_coef(310, 892.1875).unwrap();
    let (radial, _) = xsf::mathieu_modcem1(310, 892.1875, 0.0);
    xsref::np_assert_allclose!(
        [coefficients[0], radial],
        [9.231_601_377_464_515e-276, 9.274_457_443_718_162e-276],
        rtol = 1e-13
    );
}

#[test]
fn test_mathieu_coef_huge_q() {
    // `km` is NaN, inf, or too large for `usize`
    for q in [f64::NAN, f64::INFINITY, 1e15, 1e300] {
        assert!(xsf::mathieu_even_coef(0, q).is_none(), "q = {q}");
        assert!(xsf::mathieu_odd_coef(1, q).is_none(), "q = {q}");
    }
}

#[test]
#[should_panic(expected = "q must be non-negative")]
fn test_mathieu_coef_negative_q() {
    let _ = xsf::mathieu_even_coef(1, -1.0);
}

#[test]
fn test_mathieu_coef_large_m() {
    // NaN in xsf, because the coefficients in `fcoef` overflow; the expected `A_200` and
    // `B_201` are from mpmath, and `A_0` and `B_1` (~1e-435) underflow to 0
    for (fc, expected) in [
        (xsf::mathieu_even_coef(200, 1.0), 0.999_998_437_383_417_9),
        (xsf::mathieu_odd_coef(201, 1.0), 0.999_998_452_894_309_6),
    ] {
        let fc = fc.unwrap();
        xsref::np_assert_allclose!([fc[100]], [expected], rtol = 1e-13);
        assert_eq!(fc[0], 0.0);
    }
}

#[test]
#[should_panic(expected = "m must be positive")]
fn test_mathieu_odd_coef_m0() {
    let _ = xsf::mathieu_odd_coef(0, 0.0);
}

#[test]
fn test_mathieu_cem_sem_km251() {
    // `|q| > 2102`, for which `mtu0` in xsf would read out of bounds, or `m > 65536`
    for (y, yp) in [
        xsf::mathieu_cem(0, 2103.0, 10.0),
        xsf::mathieu_cem(0, -2103.0, 10.0),
        xsf::mathieu_sem(1, 2103.0, 10.0),
        xsf::mathieu_sem(1, -2103.0, 10.0),
        xsf::mathieu_cem(65_537, 1.0, 10.0),
        xsf::mathieu_sem(65_537, 1.0, 10.0),
    ] {
        assert!(y.is_nan() && yp.is_nan());
    }
    assert_eq!(xsf::mathieu_sem(0, 2103.0, 10.0), (0.0, 0.0));
}

#[test]
fn test_mathieu_huge_q() {
    // `km` would overflow `int`
    for q in [1e10, f64::INFINITY, f64::NAN] {
        for (y, yp) in [
            xsf::mathieu_cem(0, q, 10.0),
            xsf::mathieu_sem(1, q, 10.0),
            xsf::mathieu_modcem1(0, q, 1.0),
            xsf::mathieu_modsem1(1, q, 1.0),
            xsf::mathieu_modcem2(0, q, 1.0),
            xsf::mathieu_modsem2(1, q, 1.0),
        ] {
            assert!(y.is_nan() && yp.is_nan(), "q = {q}");
        }
    }
}

#[test]
fn test_mathieu_ce_se_radians() {
    // the same functions as `mathieu_cem` and `mathieu_sem`, of an angle in radians
    for (m, q, deg) in [(0, 1.0, 30.0), (3, 5.0, 100.0), (4, -2.5, 45.0), (5, -10.0, 200.0)] {
        let rad = f64::to_radians(deg);
        for (actual, desired) in [
            (xsf::mathieu_ce(m, q, rad), xsf::mathieu_cem(m, q, deg)),
            (xsf::mathieu_se(m, q, rad), xsf::mathieu_sem(m, q, deg)),
        ] {
            xsref::np_assert_allclose!(
                <[f64; 2]>::from(actual),
                <[f64; 2]>::from(desired),
                rtol = 1e-13
            );
        }
    }
    // NaN, as in SciPy, unlike `mathieu_sem`
    for (q, x) in [(f64::NAN, 1.0), (f64::INFINITY, 1.0), (1.0, f64::NAN)] {
        let (f, d) = xsf::mathieu_se(0, q, x);
        assert!(f.is_nan() && d.is_nan());
    }
}
