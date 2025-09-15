use crate::bindings;

/// Round to nearest or even integer-valued float
///
/// Returns the nearest integer to x as a f64 precision floating point result.
/// If x ends in 0.5 exactly, the nearest even integer is chosen.
pub fn round(x: f64) -> f64 {
    unsafe { bindings::round(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_round() {
        testing::test::<f64, _>("round", "d-d", |x: &[f64]| round(x[0]));
    }
}
