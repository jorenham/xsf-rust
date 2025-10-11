use crate::ffi;

#[doc(hidden)]
pub fn lanczos_sum_expg_scaled(x: f64) -> f64 {
    unsafe { ffi::lanczos_sum_expg_scaled(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_lanczos_sum_expg_scaled() {
        testing::test::<f64, _>("lanczos_sum_expg_scaled", "d-d", |x: &[f64]| {
            lanczos_sum_expg_scaled(x[0])
        });
    }
}
