use std::os::raw::c_int;

use num_complex::Complex;

use crate::bindings;
use crate::bindings::xsf_impl;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait StatsArg: sealed::Sealed {
    type Output;
    fn ndtr(self) -> Self::Output;
    fn log_ndtr(self) -> Self::Output;
}

impl StatsArg for f64 {
    type Output = f64;
    fn ndtr(self) -> f64 {
        unsafe { bindings::ndtr(self) }
    }
    fn log_ndtr(self) -> Self::Output {
        unsafe { bindings::log_ndtr(self) }
    }
}

impl StatsArg for Complex<f64> {
    type Output = Complex<f64>;
    fn ndtr(self) -> Complex<f64> {
        unsafe { bindings::ndtr_1(self.into()) }.into()
    }
    fn log_ndtr(self) -> Self::Output {
        unsafe { bindings::log_ndtr_1(self.into()) }.into()
    }
}

/// Normal distribution function `F(z)` for real or complex `z`
pub fn ndtr<T: StatsArg>(z: T) -> T::Output {
    z.ndtr()
}
/// Log of `ndtr(z)`
pub fn log_ndtr<T: StatsArg>(z: T) -> T::Output {
    z.log_ndtr()
}
xsf_impl!(ndtri, (x: f64), "Inverse of `ndtr`");

xsf_impl!(kolmogorov, (x: f64), "Kolmogorov survival function");
xsf_impl!(kolmogi, (x: f64), "Inverse of `kolmogorov`");
xsf_impl!(kolmogc, (x: f64), "Kolmogorov distribution function");
xsf_impl!(kolmogci, (x: f64), "Inverse of `kolmogc`");
xsf_impl!(kolmogp, (x: f64), "Derivative of `kolmogorov`");

xsf_impl!(smirnov, (n: c_int, x: f64), "Kolmogorov-Smirnov survival function");
xsf_impl!(smirnovc, (n: c_int, x: f64), "Kolmogorov-Smirnov distribution function");
xsf_impl!(smirnovi, (n: c_int, x: f64), "Inverse of `smirnov`");
xsf_impl!(smirnovci, (n: c_int, x: f64), "Inverse of `smirnovc`");
xsf_impl!(smirnovp, (n: c_int, x: f64), "Derivative of `smirnov`");

xsf_impl!(tukeylambdacdf, (x: f64, lmbda: f64), "Tukey-Lambda distribution function");

xsf_impl!(owens_t, (h: f64, a: f64), "Owen's T function");

xsf_impl!(chdtr, (df: f64, x: f64), "Chi-squared distribution function");
xsf_impl!(chdtrc, (df: f64, x: f64), "Chi-squared survival function");
xsf_impl!(chdtri, (df: f64, y: f64), "Chi-squared quantile function");

xsf_impl!(fdtr, (a: f64, b: f64, x: f64), "F distribution function");
xsf_impl!(fdtrc, (a: f64, b: f64, x: f64), "F survival function");
xsf_impl!(fdtri, (a: f64, b: f64, y: f64), "F quantile function");

xsf_impl!(gdtr, (a: f64, b: f64, x: f64), "Gamma distribution function");
xsf_impl!(gdtrc, (a: f64, b: f64, x: f64), "Gamma survival function");

xsf_impl!(pdtr, (k: f64, m: f64), "Poisson distribution function");
xsf_impl!(pdtri, (k: c_int, y: f64), "Poisson quantile function");
xsf_impl!(pdtrc, (k: f64, m: f64), "Poisson survival function");

xsf_impl!(bdtr, (k: f64, n: c_int, p: f64), "Binomial distribution function");
xsf_impl!(bdtrc, (k: f64, n: c_int, p: f64), "Binomial survival function");
xsf_impl!(bdtri, (k: f64, n: c_int, y: f64), "Binomial quantile function");

xsf_impl!(nbdtr, (k: c_int, n: c_int, p: f64), "Negative binomial distribution function");
xsf_impl!(nbdtrc, (k: c_int, n: c_int, p: f64), "Negative binomial survival function");
xsf_impl!(nbdtri, (k: c_int, n: c_int, p: f64), "Negative binomial quantile function");
