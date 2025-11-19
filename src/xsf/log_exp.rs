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
    #[test]
    fn test_expit() {
        xsref::test("expit", "d-d", |x| crate::expit(x[0]));
    }

    #[test]
    fn test_exprel() {
        xsref::test("exprel", "d-d", |x| crate::exprel(x[0]));
    }

    #[test]
    fn test_logit() {
        xsref::test("logit", "d-d", |x| crate::logit(x[0]));
    }

    #[test]
    fn test_log_expit() {
        xsref::test("log_expit", "d-d", |x| crate::log_expit(x[0]));
    }
}
