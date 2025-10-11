use crate::xsf::fp_error_metrics::{extended_absolute_error, extended_relative_error};

/// Similar to `numpy.testing.allclose`
#[inline(always)]
pub(crate) fn np_assert_allclose(actual: &[f64], expected: &[f64], rtol: f64, atol: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &e) in actual.iter().zip(expected.iter()) {
        if e.is_nan() {
            assert!(a.is_nan(), "expected NaN but got {a}");
        } else {
            let (err, tol) = if atol == 0.0 {
                (extended_relative_error(a, e), rtol)
            } else if rtol == 0.0 {
                (extended_absolute_error(a, e), atol)
            } else {
                (extended_absolute_error(a, e), atol + rtol * e.abs())
            };
            assert!(
                err <= tol,
                "actual: {a}, desired: {e}, error: {err:.3e}, tol: {tol:.3e}"
            );
        }
    }
}

/// Similar to `numpy.testing.assert_equal`
#[inline(always)]
pub(crate) fn np_assert_equal(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &e) in actual.iter().zip(expected.iter()) {
        if e.is_nan() {
            assert!(a.is_nan(), "expected NaN but got {a}");
        } else {
            assert_eq!(a, e, "desired {e} but got {a}");
        }
    }
}
