use crate::bindings::xsf_impl;

xsf_impl!(binom, (n: f64, k: f64), "Binomial coefficient");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_binom() {
        xsref::test::<f64, _>("binom", "d_d-d", |x: &[f64]| binom(x[0], x[1]));
    }
}
