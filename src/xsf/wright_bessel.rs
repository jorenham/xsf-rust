/// Wright's generalized Bessel function
///
/// See [`log_wright_bessel`] for the natural logarithm of this function.
pub fn wright_bessel(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::wright_bessel(a, b, x) }
}

/// Natural logarithm of Wright's generalized Bessel function, see [`wright_bessel`]
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
