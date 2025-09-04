use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait LogGammaArg: sealed::Sealed {
    type Output;
    fn xsf_loggamma(self) -> Self::Output;
    fn xsf_rgamma(self) -> Self::Output;
}

impl LogGammaArg for f64 {
    type Output = f64;
    fn xsf_loggamma(self) -> f64 {
        unsafe { bindings::loggamma(self) }
    }
    fn xsf_rgamma(self) -> f64 {
        unsafe { bindings::rgamma(self) }
    }
}

impl LogGammaArg for Complex<f64> {
    type Output = Complex<f64>;
    fn xsf_loggamma(self) -> Complex<f64> {
        unsafe { bindings::loggamma_1(self.into()) }.into()
    }
    fn xsf_rgamma(self) -> Complex<f64> {
        unsafe { bindings::rgamma_1(self.into()) }.into()
    }
}

/// Principal branch of the logarithm of `gamma(z)`
pub fn loggamma<T: LogGammaArg>(z: T) -> T::Output {
    z.xsf_loggamma()
}

/// Reciprocal Gamma function `1 / gamma(z)`
pub fn rgamma<T: LogGammaArg>(z: T) -> T::Output {
    z.xsf_rgamma()
}
