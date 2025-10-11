use crate::ffi::xsf_impl;

xsf_impl!(binom, (n: f64, k: f64), "Binomial coefficient");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_binom() {
        testing::test::<f64, _>("binom", "d_d-d", |x: &[f64]| binom(x[0], x[1]));
    }
}
