/// Wright's generalized Bessel function
///
/// Corresponds to [`scipy.special.wright_bessel`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wright_bessel.html
///
/// # See also
/// - [`log_wright_bessel`]: natural logarithm of this function
#[must_use]
#[inline]
pub fn wright_bessel(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::wright_bessel(a, b, x) }
}

/// Natural logarithm of [`wright_bessel`]
///
/// Corresponds to [`scipy.special.log_wright_bessel`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_wright_bessel.html
#[must_use]
#[inline]
pub fn log_wright_bessel(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::log_wright_bessel(a, b, x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_wright_bessel_f64() {
        xsref::test("wright_bessel", "d_d_d-d", |x| {
            crate::wright_bessel(x[0], x[1], x[2])
        });
    }

    #[test]
    fn test_log_wright_bessel_f64() {
        xsref::test("log_wright_bessel", "d_d_d-d", |x| {
            crate::log_wright_bessel(x[0], x[1], x[2])
        });
    }
}
