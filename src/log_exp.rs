use crate::bindings::xsf_impl;

xsf_impl!(expit, (x: f64), "Expit function, `1/(1 + exp(-x))`");
xsf_impl!(exprel, (x: f64), "Relative error exponential, `(exp(x) - 1)/x`");

xsf_impl!(logit, (x: f64), "Logit function, `log(x / (1 - x))`");
xsf_impl!(log_expit, (x: f64), "Log of the expit function, `log(expit(x))`");
xsf_impl!(log1mexp, (x: f64), "Compute `log(1 - exp(x))`");

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
