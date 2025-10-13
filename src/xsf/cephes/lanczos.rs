#[doc(hidden)]
pub fn lanczos_sum_expg_scaled(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::lanczos_sum_expg_scaled(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_lanczos_sum_expg_scaled() {
        crate::xsref::test("lanczos_sum_expg_scaled", "d-d", |x| {
            crate::lanczos_sum_expg_scaled(x[0])
        });
    }
}
