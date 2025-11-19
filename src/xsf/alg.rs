/// Cube root
pub fn cbrt(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cbrt(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_cbrt() {
        xsref::test("cbrt", "d-d", |x| crate::cbrt(x[0]));
    }
}
