use crate::bindings;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExpIntArg: sealed::Sealed {
    fn expi(self) -> Self;
    fn exp1(self) -> Self;
}

impl ExpIntArg for f64 {
    #[inline(always)]
    fn expi(self) -> f64 {
        unsafe { bindings::expi(self) }
    }
    #[inline(always)]
    fn exp1(self) -> f64 {
        unsafe { bindings::exp1(self) }
    }
}

impl ExpIntArg for Complex<f64> {
    #[inline(always)]
    fn expi(self) -> Complex<f64> {
        unsafe { bindings::expi_1(self.into()) }.into()
    }
    #[inline(always)]
    fn exp1(self) -> Complex<f64> {
        unsafe { bindings::exp1_1(self.into()) }.into()
    }
}

/// Exponential integral `E_i(x)` for real or complex input
pub fn expi<T: ExpIntArg>(x: T) -> T {
    x.expi()
}

/// Exponential integral `E_i(x)` for real or complex input
pub fn exp1<T: ExpIntArg>(x: T) -> T {
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
