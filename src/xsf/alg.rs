/// Cube root
pub fn cbrt(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cbrt(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_cbrt() {
        xsref::test::<f64, _>("cbrt", "d-d", |x: &[f64]| cbrt(x[0]));
    }
}
