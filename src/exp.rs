use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait ExpArg: sealed::Sealed {
    type Output;
    fn expm1(self) -> Self::Output;
}

impl ExpArg for f64 {
    type Output = f64;

    fn expm1(self) -> f64 {
        unsafe { bindings::expm1_(self) }
    }
}

impl ExpArg for Complex<f64> {
    type Output = Complex<f64>;

    fn expm1(self) -> Complex<f64> {
        unsafe { bindings::expm1__1(self.into()) }.into()
    }
}

/// `exp(x) - 1`
pub fn expm1<T: ExpArg>(x: T) -> T::Output {
    x.expm1()
}

/// `2^x`
pub fn exp2(x: f64) -> f64 {
    unsafe { bindings::exp2_(x) }
}

/// `10^x`
pub fn exp10(x: f64) -> f64 {
    unsafe { bindings::exp10_(x) }
}
