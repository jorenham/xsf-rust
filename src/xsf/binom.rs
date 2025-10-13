/// Binomial coefficient
#[doc(alias = "comb")]
pub fn binom(n: f64, k: f64) -> f64 {
    unsafe { crate::ffi::xsf::binom(n, k) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_binom() {
        crate::xsref::test("binom", "d_d-d", |x| crate::binom(x[0], x[1]));
    }
}
