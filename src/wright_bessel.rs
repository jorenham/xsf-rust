use crate::bindings::xsf_impl;

xsf_impl!(wright_bessel, (a: f64, b: f64, x: f64), "Wright's generalized Bessel function");
xsf_impl!(log_wright_bessel, (a: f64, b: f64, x: f64), "Logarithm of `wright_bessel`");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_wright_bessel_f64() {
        xsref::test::<f64, _>("wright_bessel", "d_d_d-d", |x: &[f64]| {
            wright_bessel(x[0], x[1], x[2])
        });
    }

    #[test]
    fn test_log_wright_bessel_f64() {
        xsref::test::<f64, _>("log_wright_bessel", "d_d_d-d", |x: &[f64]| {
            log_wright_bessel(x[0], x[1], x[2])
        });
    }
}
