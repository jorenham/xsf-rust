use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait Hyp2F1Arg: sealed::Sealed {
    type Output;
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self::Output;
}

impl Hyp2F1Arg for f64 {
    type Output = f64;
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> f64 {
        unsafe { bindings::hyp2f1(self, a, b, c) }
    }
}

impl Hyp2F1Arg for Complex<f64> {
    type Output = Complex<f64>;
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Complex<f64> {
        unsafe { bindings::hyp2f1_1(a, b, c, self.into()) }.into()
    }
}

/// Gauss hypergeometric function `2F1(a, b; c; z)`
pub fn hyp2f1<T: Hyp2F1Arg>(a: f64, b: f64, c: f64, z: T) -> T::Output {
    z.hyp2f1(a, b, c)
}
