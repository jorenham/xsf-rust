use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExpArg: sealed::Sealed {
    fn expm1(self) -> Self;
}

impl ExpArg for f64 {
    #[inline(always)]
    fn expm1(self) -> Self {
        unsafe { crate::ffi::xsf::expm1(self) }
    }
}

impl ExpArg for Complex<f64> {
    #[inline(always)]
    fn expm1(self) -> Self {
        unsafe { crate::ffi::xsf::expm1_1(self.into()) }.into()
    }
}

/// libc `exp` function
pub(crate) fn exp(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::exp(x) }
}

/// `exp(x) - 1` for real or complex input
pub fn expm1<T: ExpArg>(z: T) -> T {
    z.expm1()
}

/// `2^x`
pub fn exp2(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::exp2(x) }
}

/// `10^x`
pub fn exp10(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::exp10(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::{Complex, c64};

    #[test]
    fn test_expm1_f64() {
        xsref::test::<f64, _>("expm1", "d-d", |x: &[f64]| expm1(x[0]));
    }

    #[test]
    fn test_expm1_c64() {
        xsref::test::<Complex<f64>, _>("expm1", "cd-cd", |x: &[f64]| expm1(c64(x[0], x[1])));
    }

    #[test]
    fn test_exp2_f64() {
        xsref::test::<f64, _>("exp2", "d-d", |x: &[f64]| exp2(x[0]));
    }

    #[test]
    fn test_exp10_f64() {
        xsref::test::<f64, _>("exp10", "d-d", |x: &[f64]| exp10(x[0]));
    }
}
