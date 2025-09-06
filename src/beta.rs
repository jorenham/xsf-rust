use crate::bindings::xsf_impl;

xsf_impl!(beta, (a: f64, b: f64), "Beta function");
xsf_impl!(betaln, (a: f64, b: f64), "Natural log of `|beta|`");

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
