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
    #[test]
    fn test_beta() {
        crate::xsref::test("beta", "d_d-d", |x| crate::beta(x[0], x[1]));
    }

    #[test]
    fn test_betaln() {
        crate::xsref::test("betaln", "d_d-d", |x| crate::betaln(x[0], x[1]));
    }
}
