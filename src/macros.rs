/// Similar to `numpy.testing.assert_allclose`
#[macro_export]
macro_rules! np_assert_allclose {
    ($actual:expr, $expected:expr $(, rtol = $rtol:expr)? $(, atol = $atol:expr)?) => {{
        let actual = $actual;
        let expected = $expected;
        let rtol = $crate::np_assert_allclose!(@default_rtol $($rtol)?);
        let atol = $crate::np_assert_allclose!(@default_atol $($atol)?);

        assert_eq!(actual.len(), expected.len());
        for (&a, &e) in actual.iter().zip(expected.iter()) {
            let a: f64 = a;
            let e: f64 = e;
            if e.is_nan() {
                assert!(a.is_nan(), "expected NaN but got {}", a);
            } else {
                let (err, tol) = if atol == 0.0 {
                    ($crate::xsf::extended_relative_error(a, e), rtol)
                } else if rtol == 0.0 {
                    ($crate::xsf::extended_absolute_error(a, e), atol)
                } else {
                    ($crate::xsf::extended_absolute_error(a, e), atol + rtol * e.abs())
                };
                assert!(
                    err <= tol,
                    "actual: {}, desired: {}, error: {:.3e}, tol: {:.3e}", a, e, err, tol
                );
            }
        }
    }};
    (@default_rtol) => { 1e-7 };
    (@default_rtol $rtol:expr) => { $rtol };
    (@default_atol) => { 0.0 };
    (@default_atol $atol:expr) => { $atol };
}

/// Similar to `numpy.testing.assert_equal`
#[macro_export]
macro_rules! np_assert_equal {
    ($actual:expr, $expected:expr) => {{
        let actual = $actual;
        let expected = $expected;
        assert_eq!(actual.len(), expected.len());
        for (&a, &e) in actual.iter().zip(expected.iter()) {
            if e.is_nan() {
                assert!(a.is_nan(), "expected NaN but got {}", a);
            } else {
                assert_eq!(a, e, "desired {} but got {}", e, a);
            }
        }
    }};
}
