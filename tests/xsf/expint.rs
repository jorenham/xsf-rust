use core::f64::consts::PI;

use num_complex::c64;
use xsref::np_assert_allclose;

xsref::tests! {
    test_exp1_f64: "exp1", "d-d", |x| xsf::exp1(x[0]);
    test_exp1_c64: "exp1", "cd-cd", |x| xsf::exp1(c64(x[0], x[1]));
    test_expi_f64: "expi", "d-d", |x| xsf::expi(x[0]);
    test_expi_c64: "expi", "cd-cd", |x| xsf::expi(c64(x[0], x[1]));
    test_scaled_exp1_f64: "scaled_exp1", "d-d", |x| xsf::scaled_exp1(x[0]);
}

#[test]
fn test_exp1_c64_negative_real_cut() {
    // the sign of `Im z` picks the side of the cut, also for `|z| >= 40` (continued fraction)
    for x in [-39.0, -40.0, -100.0] {
        for y in [0.0, -0.0] {
            assert_eq!(xsf::exp1(c64(x, y)).im, -core::f64::consts::PI.copysign(y));
        }
    }
}

#[test]
fn test_exp1_c64_near_negative_real_cut() {
    // and just off it (mpmath)
    for (x, y, im) in [
        (-40.0, 1e-30, 3.141_592_653_589_787_3),
        (-40.0, 1e-10, -588_460.025_499_896_4),
        (-45.0, 1e-17, -4.621_578_692_518_006),
        (-100.0, 1e-50, 3.141_592_650_901_676_3),
        (-39.999_999_999_999_99, 1e-6, -5_884_631_667.782_933),
        (-40.125, 0.24, -1_580_841_922_636_904.8),
        (-40.0, 0.99, -4_961_908_935_273_592.0),
    ] {
        for (y, im) in [(-y, im), (y, -im)] {
            np_assert_allclose!([xsf::exp1(c64(x, y)).im], [im], rtol = 1e-15);
        }
    }
}

#[test]
fn test_exp1_c64_near_negative_real_cut_exp_overflow() {
    // `exp(710)` overflows, but E1 does not, on both sides of `|Im z| = 1` (Arb)
    let below_one = f64::from_bits(1.0_f64.to_bits() - 1);
    for (x, y, re, im) in [
        (-710.0, 0.75, -2.307_759_350_779_169_3e305, 2.145_345_040_598_677_4e305),
        (-709.8, below_one, -1.397_299_739_654_388e305, 2.169_427_256_952_295e305),
        (-709.8, 1.0, -1.397_299_739_654_387_9e305, 2.169_427_256_952_295_3e305),
    ] {
        for (y, im) in [(y, im), (-y, -im)] {
            let res = xsf::exp1(c64(x, y));
            np_assert_allclose!([res.re, res.im], [re, im], rtol = 1e-15);
        }
    }
    for (y, im) in [(0.0, -PI), (-0.0, PI)] {
        assert_eq!(xsf::exp1(c64(-710.0, y)).im, im);
    }
    // `y / x` is subnormal (Arb)
    for (x, y, im) in [
        (-750.0, 1.48e-321, 100.780_067_417_471_7),
        (-716.35, 1.48e-321, -3.141_592_653_589_528_4),
        (-710.0, 5e-310, -3.141_435_330_014_711_5),
        (-1417.0, 3.152_929_657_504_704e-305, 5.528_728_815_880_115e307),
    ] {
        for (y, im) in [(y, im), (-y, -im)] {
            np_assert_allclose!([xsf::exp1(c64(x, y)).im], [im], rtol = 1e-15);
        }
    }
}
