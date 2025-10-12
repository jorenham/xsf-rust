use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait DigammaArg: sealed::Sealed {
    fn digamma(self) -> Self;
}

impl DigammaArg for f64 {
    #[inline(always)]
    fn digamma(self) -> f64 {
        unsafe { crate::ffi::xsf::digamma(self) }
    }
}

impl DigammaArg for Complex<f64> {
    #[inline(always)]
    fn digamma(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::digamma_1(self.into()) }.into()
    }
}

/// Digamma function for real or complex input
#[doc(alias = "psi")]
pub fn digamma<T: DigammaArg>(x: T) -> T {
    x.digamma()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_digamma_f64() {
        xsref::test::<f64, _>("digamma", "d-d", |x: &[f64]| digamma(x[0]));
    }

    #[test]
    fn test_digamma_c64() {
        xsref::test::<num_complex::Complex<f64>, _>("digamma", "cd-cd", |x: &[f64]| {
            digamma(num_complex::c64(x[0], x[1]))
        });
    }
}
