/// Beta function
#[must_use]
#[inline]
pub fn beta(a: f64, b: f64) -> f64 {
    unsafe { crate::ffi::xsf::beta(a, b) }
}

/// Logarithm of the absolute value of [`beta`]
#[doc(alias = "lbeta", alias = "logbeta", alias = "beta_ln")]
#[must_use]
#[inline]
pub fn betaln(a: f64, b: f64) -> f64 {
    unsafe { crate::ffi::xsf::betaln(a, b) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_beta() {
        xsref::test("beta", "d_d-d", |x| crate::beta(x[0], x[1]));
    }

    #[test]
    fn test_betaln() {
        xsref::test("betaln", "d_d-d", |x| crate::betaln(x[0], x[1]));
    }
}
