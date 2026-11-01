//! based on <https://github.com/scipy/scipy/blob/main/scipy/special/tests/test_wrightomega.py>

#![allow(clippy::float_cmp)]

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
    let expected_results = [
        c64(0.0, 0.0),
        c64(0.0, -0.0),
        c64(-0.0, 0.0),
        c64(-0.0, -0.0),
    ];
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
    for (im, re) in [
        (1e51, -117.431_839_742_696_33),
        (1e308, -709.196_208_642_166),
    ] {
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
