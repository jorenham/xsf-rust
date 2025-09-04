use crate::bindings;
use crate::bindings::xsf_impl;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait GammaArg: sealed::Sealed {
    type Output;
    fn xsf_gamma(self) -> Self::Output;
}

impl GammaArg for f64 {
    type Output = f64;
    fn xsf_gamma(self) -> f64 {
        unsafe { bindings::gamma_(self) }
    }
}

impl GammaArg for Complex<f64> {
    type Output = Complex<f64>;
    fn xsf_gamma(self) -> Complex<f64> {
        unsafe { bindings::gamma__1(self.into()) }.into()
    }
}

/// Gamma function
pub fn gamma<T: GammaArg>(z: T) -> T::Output {
    z.xsf_gamma()
}

xsf_impl!(gammainc, (a: f64, x: f64), "Incomplete Gamma integral");
xsf_impl!(gammaincc, (a: f64, x: f64), "Complemented incomplete Gamma integral");
xsf_impl!(gammaincinv, (a: f64, p: f64), "Inverse of `gammainc`");
xsf_impl!(gammainccinv, (a: f64, p: f64), "Inverse of `gammaincc`");
xsf_impl!(gammaln, (x: f64), "Natural logarithm of Gamma function");
xsf_impl!(gammasgn, (x: f64), "Sign of the Gamma function");
