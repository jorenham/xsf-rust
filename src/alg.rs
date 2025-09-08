use crate::bindings;

/// Cube root
pub fn cbrt(x: f64) -> f64 {
    unsafe { bindings::cbrt(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_cbrt() {
        testing::test::<f64, _>("cbrt", "d-d", |x: &[f64]| cbrt(x[0]));
    }
}
