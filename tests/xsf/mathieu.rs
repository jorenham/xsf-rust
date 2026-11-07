use num_traits::ToPrimitive;

#[test]
fn test_mathieu_a() {
    xsref::test("cem_cva", "d_d-d", |x| {
        xsf::mathieu_a(x[0].to_u32().unwrap(), x[1])
    });
}

#[test]
fn test_mathieu_b() {
    xsref::test("sem_cva", "d_d-d", |x| {
        xsf::mathieu_b(x[0].to_u32().unwrap(), x[1])
    });
}

#[test]
fn test_mathieu_cem() {
    xsref::test("cem", "d_d_d-d_d", |x| {
        xsf::mathieu_cem(x[0].to_u32().unwrap(), x[1], x[2])
    });
}

#[test]
fn test_mathieu_sem() {
    xsref::test("sem", "d_d_d-d_d", |x| {
        xsf::mathieu_sem(x[0].to_u32().unwrap(), x[1], x[2])
    });
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

#[test]
fn test_mathieu_modcem1() {
    xsref::test("mcm1", "d_d_d-d_d", |x| {
        xsf::mathieu_modcem1(x[0].to_u32().unwrap(), x[1], x[2])
    });
}

#[test]
fn test_mathieu_modsem1() {
    xsref::test("msm1", "d_d_d-d_d", |x| {
        xsf::mathieu_modsem1(x[0].to_u32().unwrap(), x[1], x[2])
    });
}

#[test]
fn test_mathieu_modcem2() {
    xsref::test("mcm2", "d_d_d-d_d", |x| {
        xsf::mathieu_modcem2(x[0].to_u32().unwrap(), x[1], x[2])
    });
}

#[test]
fn test_mathieu_modsem2() {
    xsref::test("msm2", "d_d_d-d_d", |x| {
        xsf::mathieu_modsem2(x[0].to_u32().unwrap(), x[1], x[2])
    });
}

#[test]
fn test_mathieu_a_nan() {
    assert!(xsf::mathieu_a(13, f64::NAN).is_nan());
}

#[test]
fn test_mathieu_modcem1_large_x() {
    // `(f^2 + (f'/u)^2) * u -> 2/pi` with `u = sqrt(q) * exp(x)`, i.e. `|Mc1| ~ sqrt(2 / (pi u))`
    let x = 22.0;
    let u = f64::exp(x);
    let (f, d) = xsf::mathieu_modcem1(0, 1.0, x);
    let env = (f * f + (d / u) * (d / u)) * u;
    assert!(
        (env * core::f64::consts::FRAC_PI_2 - 1.0).abs() < 1e-6,
        "{env}"
    );
}

#[test]
fn test_mathieu_even_coef() {
    // Smoke test: verify function returns Ok and non-empty coefficients
    let result = xsf::mathieu_even_coef(2, 1.0);
    assert!(result.is_ok());
    let coefs = result.unwrap();
    assert_ne!(coefs, [] as [f64; 0]);
    assert!(coefs.iter().all(|&c| c.is_finite()));

    // Test with q = 0
    let result = xsf::mathieu_even_coef(0, 0.0);
    assert!(result.is_ok());

    // Test error case: negative q
    let result = xsf::mathieu_even_coef(1, -1.0);
    assert!(result.is_err());
}

#[test]
fn test_mathieu_odd_coef() {
    // Smoke test: verify function returns Ok and non-empty coefficients
    let result = xsf::mathieu_odd_coef(1, 1.0);
    assert!(result.is_ok());
    let coefs = result.unwrap();
    assert_ne!(coefs, [] as [f64; 0]);
    assert!(coefs.iter().all(|&c| c.is_finite()));

    // Test with q = 0
    let result = xsf::mathieu_odd_coef(1, 0.0);
    assert!(result.is_ok());

    // Test error case: negative q
    let result = xsf::mathieu_odd_coef(1, -1.0);
    assert!(result.is_err());
}

#[test]
fn test_mathieu_coef_huge_q() {
    // `km` is NaN, inf, or too large for `usize`
    for q in [f64::NAN, f64::INFINITY, 1e15, 1e300] {
        assert!(xsf::mathieu_even_coef(0, q).is_err(), "q = {q}");
        assert!(xsf::mathieu_odd_coef(1, q).is_err(), "q = {q}");
    }
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
fn test_mathieu_odd_coef_m0() {
    assert!(xsf::mathieu_odd_coef(0, 0.0).is_err());
}

#[test]
fn test_mathieu_cem_sem_km251() {
    // `mtu0` would read out of bounds
    for (y, yp) in [
        xsf::mathieu_cem(0, 2103.0, 10.0),
        xsf::mathieu_cem(0, -2103.0, 10.0),
        xsf::mathieu_sem(1, 2100.0, 10.0),
        xsf::mathieu_sem(1, -2100.0, 10.0),
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
