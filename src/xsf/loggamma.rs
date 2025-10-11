use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait LogGammaArg: sealed::Sealed {
    fn xsf_loggamma(self) -> Self;
    fn xsf_rgamma(self) -> Self;
}

impl LogGammaArg for f64 {
    #[inline(always)]
    fn xsf_loggamma(self) -> f64 {
        unsafe { crate::ffi::xsf::loggamma(self) }
    }
    #[inline(always)]
    fn xsf_rgamma(self) -> f64 {
        unsafe { crate::ffi::xsf::rgamma(self) }
    }
}

impl LogGammaArg for Complex<f64> {
    #[inline(always)]
    fn xsf_loggamma(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::loggamma_1(self.into()) }.into()
    }
    #[inline(always)]
    fn xsf_rgamma(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::rgamma_1(self.into()) }.into()
    }
}

/// Principal branch of the logarithm of `gamma(z)`
pub fn loggamma<T: LogGammaArg>(z: T) -> T {
    z.xsf_loggamma()
}

/// Reciprocal Gamma function `1 / gamma(z)`
pub fn rgamma<T: LogGammaArg>(z: T) -> T {
    z.xsf_rgamma()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::{Complex, c64};

    // loggamma

    #[test]
    fn test_loggamma_f64() {
        xsref::test::<f64, _>("loggamma", "d-d", |x: &[f64]| loggamma(x[0]));
    }

    #[test]
    fn test_loggamma_c64() {
        xsref::test::<Complex<f64>, _>("loggamma", "cd-cd", |x: &[f64]| loggamma(c64(x[0], x[1])));
    }

    // rgamma

    #[test]
    fn test_rgamma_f64() {
        xsref::test::<f64, _>("rgamma", "d-d", |x: &[f64]| rgamma(x[0]));
    }

    #[test]
    fn test_rgamma_c64() {
        xsref::test::<Complex<f64>, _>("rgamma", "cd-cd", |x: &[f64]| rgamma(c64(x[0], x[1])));
    }
}
