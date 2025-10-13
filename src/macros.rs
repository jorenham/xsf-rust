/// Similar to `numpy.testing.assert_allclose`
#[macro_export]
macro_rules! np_assert_allclose {
    // Both rtol and atol provided
    ($actual:expr, $expected:expr, rtol = $rtol:expr, atol = $atol:expr) => {{
        $crate::np_assert_allclose!(@impl $actual, $expected, $rtol, $atol)
    }};
    ($actual:expr, $expected:expr, atol = $atol:expr, rtol = $rtol:expr) => {{
        $crate::np_assert_allclose!(@impl $actual, $expected, $rtol, $atol)
    }};
    // Only rtol provided - atol defaults to 0.0
    ($actual:expr, $expected:expr, rtol = $rtol:expr) => {{
        $crate::np_assert_allclose!(@impl $actual, $expected, $rtol, 0.0)
    }};
    // Only atol provided - rtol defaults to 0.0
    ($actual:expr, $expected:expr, atol = $atol:expr) => {{
        $crate::np_assert_allclose!(@impl $actual, $expected, 0.0, $atol)
    }};
    // Neither provided - rtol = 1e-7, atol = 0.0
    ($actual:expr, $expected:expr) => {{
        $crate::np_assert_allclose!(@impl $actual, $expected, 1e-7, 0.0)
    }};
    // Internal implementation
    (@impl $actual:expr, $expected:expr, $rtol:expr, $atol:expr) => {{
        let actual = $actual;
        let expected = $expected;
        let rtol: f64 = $rtol;
        let atol: f64 = $atol;

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
