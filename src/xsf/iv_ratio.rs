/// Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind
pub fn iv_ratio(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::iv_ratio(v, x) }
}

/// Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind
pub fn iv_ratio_c(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::iv_ratio_c(v, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_iv_ratio() {
        xsref::test::<f64, _>("iv_ratio", "d_d-d", |x: &[f64]| iv_ratio(x[0], x[1]));
    }

    #[test]
    fn test_iv_ratio_c() {
        xsref::test::<f64, _>("iv_ratio_c", "d_d-d", |x: &[f64]| iv_ratio_c(x[0], x[1]));
    }
}
