/// Integral of the Struve function of order 0
#[doc(alias = "itstruve_h0")]
#[doc(alias = "it1struve0")]
pub fn itstruve0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::itstruve0(x) }
}

/// Integral related to the Struve function of order 0
#[doc(alias = "it2struve_l0")]
pub fn it2struve0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::it2struve0(x) }
}

/// Integral of the modified Struve function of order 0
#[doc(alias = "itstruve_l0")]
pub fn itmodstruve0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::itmodstruve0(x) }
}

/// Struve `H` function
#[doc(alias = "struve")]
pub fn struve_h(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::struve_h(v, x) }
}

/// Struve `L` function
#[doc(alias = "modstruve")]
pub fn struve_l(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::struve_l(v, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_itstruve0_f64() {
        xsref::test::<f64, _>("itstruve0", "d-d", |x: &[f64]| itstruve0(x[0]));
    }

    #[test]
    fn test_it2struve0_f64() {
        xsref::test::<f64, _>("it2struve0", "d-d", |x: &[f64]| it2struve0(x[0]));
    }

    #[test]
    fn test_itmodstruve0_f64() {
        xsref::test::<f64, _>("itmodstruve0", "d-d", |x: &[f64]| itmodstruve0(x[0]));
    }

    #[test]
    fn test_struve_h_f64() {
        xsref::test::<f64, _>("struve_h", "d_d-d", |x: &[f64]| struve_h(x[0], x[1]));
    }

    #[test]
    fn test_struve_l_f64() {
        xsref::test::<f64, _>("struve_l", "d_d-d", |x: &[f64]| struve_l(x[0], x[1]));
    }
}
