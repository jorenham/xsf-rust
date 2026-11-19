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
