use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait DigammaArg: sealed::Sealed {
    type Output;
    fn digamma(self) -> Self::Output;
}

impl DigammaArg for f64 {
    type Output = f64;

    fn digamma(self) -> f64 {
        unsafe { bindings::digamma(self) }
    }
}

impl DigammaArg for Complex<f64> {
    type Output = Complex<f64>;

    fn digamma(self) -> Complex<f64> {
        unsafe { bindings::digamma_1(self.into()) }.into()
    }
}

/// Digamma function
pub fn digamma<T: DigammaArg>(x: T) -> T::Output {
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
