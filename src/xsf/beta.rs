/// Beta function
pub fn beta(a: f64, b: f64) -> f64 {
    unsafe { crate::ffi::xsf::beta(a, b) }
}

/// Logarithm of the absolute value of [`beta`]
#[doc(alias = "lbeta")]
#[doc(alias = "logbeta")]
#[doc(alias = "beta_ln")]
pub fn betaln(a: f64, b: f64) -> f64 {
    unsafe { crate::ffi::xsf::betaln(a, b) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_beta() {
        xsref::test::<f64, _>("beta", "d_d-d", |x: &[f64]| beta(x[0], x[1]));
    }

    #[test]
    fn test_betaln() {
        xsref::test::<f64, _>("betaln", "d_d-d", |x: &[f64]| betaln(x[0], x[1]));
    }
}
