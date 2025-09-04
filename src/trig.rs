use num_complex::Complex;

use crate::bindings;
use crate::bindings::xsf_impl;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait TrigArg: sealed::Sealed {
    type Output;
    fn sinpi(self) -> Self::Output;
    fn cospi(self) -> Self::Output;
}

impl TrigArg for f64 {
    type Output = f64;
    fn sinpi(self) -> f64 {
        unsafe { bindings::sinpi(self) }
    }
    fn cospi(self) -> Self::Output {
        unsafe { bindings::cospi(self) }
    }
}

impl TrigArg for Complex<f64> {
    type Output = Complex<f64>;
    fn sinpi(self) -> Complex<f64> {
        unsafe { bindings::sinpi_1(self.into()) }.into()
    }
    fn cospi(self) -> Self::Output {
        unsafe { bindings::cospi_1(self.into()) }.into()
    }
}

/// Compute `sin(pi*z)` for real or complex `z`
pub fn sinpi<T: TrigArg>(z: T) -> T::Output {
    z.sinpi()
}
/// Compute `cos(pi*z)` for real or complex `z`
pub fn cospi<T: TrigArg>(z: T) -> T::Output {
    z.cospi()
}

xsf_impl!(sindg, (x: f64), "Circular sine of angle in degrees");
xsf_impl!(cosdg, (x: f64), "Circular cosine of angle in degrees");
xsf_impl!(tandg, (x: f64), "Circular tangent of argument in degrees");
xsf_impl!(cotdg, (x: f64), "Circular cotangent of argument in degrees");

xsf_impl!(cosm1, (x: f64), "Compute `cos(x) - 1`");

xsf_impl!(radian, (d: f64, m: f64, s: f64), "Degrees, minutes, seconds to radians");
