/// Wright's generalized Bessel function
pub fn wright_bessel(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::wright_bessel(a, b, x) }
}

/// Logarithm of [`wright_bessel`]
pub fn log_wright_bessel(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::log_wright_bessel(a, b, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_wright_bessel_f64() {
        testing::test::<f64, _>("wright_bessel", "d_d_d-d", |x: &[f64]| {
            wright_bessel(x[0], x[1], x[2])
        });
    }

    #[test]
    fn test_log_wright_bessel_f64() {
        testing::test::<f64, _>("log_wright_bessel", "d_d_d-d", |x: &[f64]| {
            log_wright_bessel(x[0], x[1], x[2])
        });
    }
}
