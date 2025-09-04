use crate::bindings;
use crate::bindings::xsf_impl;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait LogArg: sealed::Sealed {
    type Output;
    fn xsf_log1p(self) -> Self::Output;
    fn xsf_xlogy<T: Into<Self::Output>>(self, x: T) -> Self::Output;
    fn xsf_xlog1py<T: Into<Self::Output>>(self, x: T) -> Self::Output;
}

impl LogArg for f64 {
    type Output = f64;
    fn xsf_log1p(self) -> f64 {
        unsafe { bindings::log1p_(self) }
    }
    fn xsf_xlogy<T: Into<f64>>(self, x: T) -> f64 {
        unsafe { bindings::xlogy(x.into(), self) }
    }
    fn xsf_xlog1py<T: Into<f64>>(self, x: T) -> f64 {
        unsafe { bindings::xlog1py(x.into(), self) }
    }
}

impl LogArg for Complex<f64> {
    type Output = Complex<f64>;
    fn xsf_log1p(self) -> Complex<f64> {
        unsafe { bindings::log1p__1(self.into()) }.into()
    }
    fn xsf_xlogy<T: Into<Complex<f64>>>(self, x: T) -> Complex<f64> {
        unsafe { bindings::xlogy_1(x.into().into(), self.into()) }.into()
    }
    fn xsf_xlog1py<T: Into<Complex<f64>>>(self, x: T) -> Complex<f64> {
        unsafe { bindings::xlog1py_1(x.into().into(), self.into()) }.into()
    }
}

/// `log(z + 1)` for real or complex input
pub fn log1p<T: LogArg>(z: T) -> T::Output {
    z.xsf_log1p()
}

xsf_impl!(log1pmx, (x: f64), "Compute `log(1 + x) - x` for real input");

/// Compute `x * log(y)` for real or complex input
pub fn xlogy<X, Y>(x: X, y: Y) -> Y::Output
where
    X: Into<Y::Output>,
    Y: LogArg,
{
    y.xsf_xlogy(x)
}
/// Compute `x * log(1 + y)` for real or complex input
pub fn xlog1py<X, Y>(x: X, y: Y) -> Y::Output
where
    X: Into<Y::Output>,
    Y: LogArg,
{
    y.xsf_xlog1py(x)
}
