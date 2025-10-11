/// Binomial coefficient
pub fn binom(n: f64, k: f64) -> f64 {
    unsafe { crate::ffi::xsf::binom(n, k) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_binom() {
        testing::test::<f64, _>("binom", "d_d-d", |x: &[f64]| binom(x[0], x[1]));
    }
}
