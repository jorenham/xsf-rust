use crate::bindings::xsf_impl;

xsf_impl!(gdtrib, (a: f64, p: f64, x: f64), "Inverse of `p = gdtr(a, b, x)` with respect to `b`");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_gdtrib() {
        xsref::test::<f64, _>("gdtrib", "d_d_d-d", |x: &[f64]| gdtrib(x[0], x[1], x[2]));
    }
}
