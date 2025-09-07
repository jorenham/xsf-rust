use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait ExpIntArg: sealed::Sealed {
    type Output;
    fn expi(self) -> Self::Output;
    fn exp1(self) -> Self::Output;
}

impl ExpIntArg for f64 {
    type Output = f64;

    fn expi(self) -> f64 {
        unsafe { bindings::expi(self) }
    }

    fn exp1(self) -> f64 {
        unsafe { bindings::exp1(self) }
    }
}

impl ExpIntArg for Complex<f64> {
    type Output = Complex<f64>;

    fn expi(self) -> Complex<f64> {
        unsafe { bindings::expi_1(self.into()) }.into()
    }

    fn exp1(self) -> Complex<f64> {
        unsafe { bindings::exp1_1(self.into()) }.into()
    }
}

/// Exponential integral `E_i(x)`
pub fn expi<T: ExpIntArg>(x: T) -> T::Output {
    x.expi()
}

/// Exponential integral `E_i(x)`
pub fn exp1<T: ExpIntArg>(x: T) -> T::Output {
    x.exp1()
}

/// Scaled version of the exponential integral `E_1(x)`
pub fn scaled_exp1(x: f64) -> f64 {
    unsafe { bindings::scaled_exp1(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    // expi

    #[test]
    fn test_expi_f64() {
        testing::test::<f64, _>("expi", "d-d", |x: &[f64]| expi(x[0]));
    }

    #[test]
    fn test_expi_c64() {
        testing::test::<Complex<f64>, _>("expi", "cd-cd", |x: &[f64]| expi(c64(x[0], x[1])));
    }

    // exp1

    #[test]
    fn test_exp1_f64() {
        testing::test::<f64, _>("exp1", "d-d", |x: &[f64]| exp1(x[0]));
    }

    #[test]
    fn test_exp1_c64() {
        testing::test::<Complex<f64>, _>("exp1", "cd-cd", |x: &[f64]| exp1(c64(x[0], x[1])));
    }

    // scaled_exp1

    #[test]
    fn test_scaled_exp1_f64() {
        testing::test::<f64, _>("scaled_exp1", "d-d", |x: &[f64]| scaled_exp1(x[0]));
    }
}
