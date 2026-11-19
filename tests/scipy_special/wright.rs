//! based on <https://github.com/scipy/scipy/blob/main/scipy/special/tests/test_wrightomega.py>

use core::f64::consts::PI;
use num_complex::c64;
use xsref::np_assert_allclose;

#[test]
fn test_wrightomega_nan() {
    let pts = [
        c64(f64::NAN, 0.0),
        c64(0.0, f64::NAN),
        c64(f64::NAN, f64::NAN),
        c64(f64::NAN, 1.0),
        c64(1.0, f64::NAN),
    ];
    for p in pts {
        let res = xsf::wrightomega(p);
        assert!(res.re.is_nan());
        assert!(res.im.is_nan());
    }
}

#[test]
fn test_wrightomega_inf_branch() {
    let pts = [
        c64(f64::NEG_INFINITY, PI / 4.0),
        c64(f64::NEG_INFINITY, -PI / 4.0),
        c64(f64::NEG_INFINITY, 3.0 * PI / 4.0),
        c64(f64::NEG_INFINITY, -3.0 * PI / 4.0),
    ];
    let expected_results = [c64(0.0, 0.0), c64(0.0, -0.0), c64(-0.0, 0.0), c64(-0.0, -0.0)];
    for i in 0..4 {
        let (p, expected) = (pts[i], expected_results[i]);
        let res = xsf::wrightomega(p);
        assert_eq!(res.re.signum(), expected.re.signum());
        assert_eq!(res.im.signum(), expected.im.signum());
    }
}

#[test]
fn test_wrightomega_inf() {
    let pts = [
        c64(f64::INFINITY, 10.0),
        c64(f64::NEG_INFINITY, 10.0),
        c64(10.0, f64::INFINITY),
        c64(10.0, f64::NEG_INFINITY),
    ];
    for p in pts {
        assert_eq!(xsf::wrightomega(p), p);
    }
}

#[test]
fn test_wrightomega_huge() {
    // reference values from scipy.special.wrightomega
    for (im, re) in [(1e51, -117.431_839_742_696_33), (1e308, -709.196_208_642_166)] {
        for im in [im, -im] {
            let res = xsf::wrightomega(c64(0.0, im));
            np_assert_allclose!([res.re], [re], rtol = 1e-15);
            assert_eq!(res.im, im);
        }
    }
    // |z| overflows, but omega(z) rounds to z
    let z = c64(1.3e308, 1.3e308);
    assert_eq!(xsf::wrightomega(z), z);
}

#[test]
fn test_wrightomega_pi_inside_strip() {
    // `y = +-PI` is inside the strip `|y| < pi`, so that `omega(conj(z)) = conj(omega(z))` (mpmath)
    for (x, re, im) in [
        (-60.0, -8.756_510_762_696_52e-27, 1.072_363_287_723_564_2e-42),
        (-3.0, -0.052_469_097_457_714_87, 6.781_426_556_467_746e-18),
        (-2.0, -0.158_594_339_563_039_37, 2.308_304_536_576_337_7e-17),
        (-1.5, -0.301_709_562_684_336, 5.291_317_630_438_808e-17),
    ] {
        for (y, im) in [(PI, im), (-PI, -im)] {
            let res = xsf::wrightomega(c64(x, y));
            np_assert_allclose!([res.re, res.im], [re, im], rtol = 1e-15);
        }
    }
    // underflowing to the signed zeros of `exp(z)`
    let neg_zero = (-0.0_f64).to_bits();
    for x in [-800.0, f64::NEG_INFINITY] {
        let res = xsf::wrightomega(c64(x, -PI));
        assert_eq!([res.re.to_bits(), res.im.to_bits()], [neg_zero; 2]);
    }

    // the wing beyond the cut, `omega(z) = z - pi i - ln(pi i - z)` (mpmath)
    let res = xsf::wrightomega(c64(-1e300, f64::from_bits(PI.to_bits() + 1)));
    np_assert_allclose!([res.im], [3.216_245_299_353_273e-16], rtol = 1e-15);
}

#[test]
fn test_wrightomega_singular() {
    let pts = [c64(-1.0, PI), c64(-1.0, -PI)];
    for p in pts {
        let res = xsf::wrightomega(p);
        assert_eq!(res, c64(-1.0, 0.0));
        assert!(res.im.is_sign_positive());
    }
}

#[test]
fn test_wrightomega_real_infinities() {
    assert_eq!(xsf::wrightomega(f64::NEG_INFINITY), 0.0);
    assert_eq!(xsf::wrightomega(f64::INFINITY), f64::INFINITY);
}

#[test]
fn test_wrightomega_real_nan() {
    assert!(xsf::wrightomega(f64::NAN).is_nan());
}

#[test]
fn test_wrightomega_real_series_crossover() {
    // the neighbours of the crossover at 1e20 (mpmath)
    let x = 1e20_f64.to_bits();
    let actual = [x - 1, x + 1].map(|b| xsf::wrightomega(f64::from_bits(b)));
    let desired = [99_999_999_999_999_983_569.948, 100_000_000_000_000_016_337.948];
    np_assert_allclose!(actual, desired, atol = 0.0, rtol = 2.0 * f64::EPSILON);
}

#[test]
fn test_wrightomega_exp_approximation_crossover() {
    // the neighbours of the crossover at -50, towards and away from 0 (mpmath)
    let x = (-50.0_f64).to_bits();
    let actual = [x - 1, x + 1].map(|b| xsf::wrightomega(f64::from_bits(b)));
    let desired = [1.928_749_847_963_931_5e-22, 1.928_749_847_963_904_2e-22];
    np_assert_allclose!(actual, desired, atol = 0.0, rtol = 2.0 * f64::EPSILON);
}

#[test]
fn test_wrightomega_real_versus_complex() {
    // x = np.linspace(-500, 500, 1001)
    let res_f: Vec<f64> = (-500..=500)
        .map(|x| xsf::wrightomega(f64::from(x)))
        .collect();
    let res_c: Vec<f64> = (-500..=500)
        .map(|x| xsf::wrightomega(c64(f64::from(x), 0.0)).re)
        .collect();
    np_assert_allclose!(res_c, res_f, atol = 0.0, rtol = 1e-14);
}
