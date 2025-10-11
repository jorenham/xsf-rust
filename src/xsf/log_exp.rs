/// Expit function, `1/(1 + exp(-x))`
pub fn expit(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::expit(x) }
}

/// Relative error exponential, `(exp(x) - 1)/x`
pub fn exprel(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::exprel(x) }
}

/// Logit function, `log(x / (1 - x))`
pub fn logit(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::logit(x) }
}

/// Log of the expit function, `log(expit(x))`
pub fn log_expit(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::log_expit(x) }
}

/// Compute `log(1 - exp(x))`
pub fn log1mexp(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::log1mexp(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_expit() {
        testing::test::<f64, _>("expit", "d-d", |x: &[f64]| expit(x[0]));
    }

    #[test]
    fn test_exprel() {
        testing::test::<f64, _>("exprel", "d-d", |x: &[f64]| exprel(x[0]));
    }

    #[test]
    fn test_logit() {
        testing::test::<f64, _>("logit", "d-d", |x: &[f64]| logit(x[0]));
    }

    #[test]
    fn test_log_expit() {
        testing::test::<f64, _>("log_expit", "d-d", |x: &[f64]| log_expit(x[0]));
    }
}
