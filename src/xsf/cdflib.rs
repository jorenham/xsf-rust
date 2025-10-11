/// Inverse of `p = gdtr(a, b, x)` with respect to `b`
pub fn gdtrib(a: f64, p: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gdtrib(a, p, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_gdtrib() {
        testing::test::<f64, _>("gdtrib", "d_d_d-d", |x: &[f64]| gdtrib(x[0], x[1], x[2]));
    }
}
