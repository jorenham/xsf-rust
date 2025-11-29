/// Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind
#[must_use]
#[inline]
pub fn iv_ratio(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::iv_ratio(v, x) }
}

/// Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind
#[must_use]
#[inline]
pub fn iv_ratio_c(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::iv_ratio_c(v, x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_iv_ratio() {
        xsref::test("iv_ratio", "d_d-d", |x| crate::iv_ratio(x[0], x[1]));
    }

    #[test]
    fn test_iv_ratio_c() {
        xsref::test("iv_ratio_c", "d_d-d", |x| crate::iv_ratio_c(x[0], x[1]));
    }
}
