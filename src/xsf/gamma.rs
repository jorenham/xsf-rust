use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait GammaArg: sealed::Sealed {
    fn xsf_gamma(self) -> Self;
}

impl GammaArg for f64 {
    #[inline(always)]
    fn xsf_gamma(self) -> f64 {
        unsafe { crate::ffi::xsf::gamma(self) }
    }
}

impl GammaArg for Complex<f64> {
    #[inline(always)]
    fn xsf_gamma(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::gamma_1(self.into()) }.into()
    }
}

/// Gamma function for real or complex input
pub fn gamma<T: GammaArg>(z: T) -> T {
    z.xsf_gamma()
}

/// Regularized lower incomplete gamma function
#[doc(alias = "inc_gamma")]
#[doc(alias = "gamma_lr")]
pub fn gammainc(a: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammainc(a, x) }
}

/// Regularized upper incomplete gamma function
#[doc(alias = "gamma_ur")]
pub fn gammaincc(a: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaincc(a, x) }
}

/// Inverse of [`gammainc`]
pub fn gammaincinv(a: f64, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaincinv(a, p) }
}

/// Inverse of [`gammaincc`]
pub fn gammainccinv(a: f64, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammainccinv(a, p) }
}

/// Logarithm of the absolute value of the gamma function
pub fn gammaln(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaln(x) }
}

/// Sign of the Gamma function
pub fn gammasgn(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammasgn(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::{Complex, c64};

    #[test]
    fn test_gamma_f64() {
        xsref::test::<f64, _>("gamma", "d-d", |x: &[f64]| gamma(x[0]));
    }

    #[test]
    fn test_gamma_c64() {
        xsref::test::<Complex<f64>, _>("gamma", "cd-cd", |x: &[f64]| gamma(c64(x[0], x[1])));
    }

    #[test]
    fn test_gammainc() {
        xsref::test::<f64, _>("gammainc", "d_d-d", |x: &[f64]| gammainc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincc() {
        xsref::test::<f64, _>("gammaincc", "d_d-d", |x: &[f64]| gammaincc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincinv() {
        xsref::test::<f64, _>("gammaincinv", "d_d-d", |x: &[f64]| gammaincinv(x[0], x[1]));
    }

    #[test]
    fn test_gammainccinv() {
        xsref::test::<f64, _>("gammainccinv", "d_d-d", |x: &[f64]| {
            gammainccinv(x[0], x[1])
        });
    }

    #[test]
    fn test_gammaln() {
        xsref::test::<f64, _>("gammaln", "d-d", |x: &[f64]| gammaln(x[0]));
    }

    #[test]
    fn test_gammasgn() {
        xsref::test::<f64, _>("gammasgn", "d-d", |x: &[f64]| gammasgn(x[0]));
    }
}
