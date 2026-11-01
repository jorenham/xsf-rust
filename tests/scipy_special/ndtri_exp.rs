//! Translated from `scipy.special.tests.TestNdtriExp` at
//! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_ndtri_exp.py>

use xsref::np_assert_allclose;

const UNIFORM_RANDOM_POINTS: [f64; 20] = [
    0.191_519_45,
    0.622_108_77,
    0.437_727_74,
    0.785_358_58,
    0.779_975_81,
    0.272_592_61,
    0.276_464_26,
    0.801_872_18,
    0.958_139_35,
    0.875_932_63,
    0.357_817_27,
    0.500_995_13,
    0.683_462_94,
    0.712_702_03,
    0.370_250_75,
    0.561_196_19,
    0.503_083_17,
    0.013_768_45,
    0.772_826_62,
    0.882_641_19,
];

fn log_ndtr_ndtri_exp(y: f64) -> f64 {
    xsf::log_ndtr(xsf::ndtri_exp(y))
}

#[test]
fn test_very_small_arg() {
    for scale in [-1e1, -1e2, -1e10, -1e20, -f64::MAX] {
        let points = UNIFORM_RANDOM_POINTS.map(|p| scale * (0.5 * p + 0.5));
        np_assert_allclose!(points.map(log_ndtr_ndtri_exp), points, rtol = 1e-14);
    }
}

#[test]
#[allow(clippy::float_cmp)]
fn test_asymptotes() {
    assert_eq!(xsf::ndtri_exp(f64::NEG_INFINITY), f64::NEG_INFINITY);
    assert_eq!(xsf::ndtri_exp(0.0), f64::INFINITY);
}

#[test]
fn test_outside_domain() {
    assert!(xsf::ndtri_exp(1.0).is_nan());
}
