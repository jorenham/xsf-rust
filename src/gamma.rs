use crate::bindings;
use crate::bindings::xsf_impl;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait GammaArg: sealed::Sealed {
    fn xsf_gamma(self) -> Self;
}

impl GammaArg for f64 {
    #[inline(always)]
    fn xsf_gamma(self) -> f64 {
        unsafe { bindings::gamma(self) }
    }
}

impl GammaArg for Complex<f64> {
    #[inline(always)]
    fn xsf_gamma(self) -> Complex<f64> {
        unsafe { bindings::gamma_1(self.into()) }.into()
    }
}

/// Gamma function for real or complex input
pub fn gamma<T: GammaArg>(z: T) -> T {
    z.xsf_gamma()
}

xsf_impl!(gammainc, (a: f64, x: f64), "Incomplete Gamma integral");
xsf_impl!(gammaincc, (a: f64, x: f64), "Complemented incomplete Gamma integral");
xsf_impl!(gammaincinv, (a: f64, p: f64), "Inverse of `gammainc`");
xsf_impl!(gammainccinv, (a: f64, p: f64), "Inverse of `gammaincc`");
xsf_impl!(gammaln, (x: f64), "Natural logarithm of Gamma function");
xsf_impl!(gammasgn, (x: f64), "Sign of the Gamma function");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    #[test]
    fn test_gamma_f64() {
        testing::test::<f64, _>("gamma", "d-d", |x: &[f64]| gamma(x[0]));
    }

    #[test]
    fn test_gamma_c64() {
        testing::test::<Complex<f64>, _>("gamma", "cd-cd", |x: &[f64]| gamma(c64(x[0], x[1])));
    }

    #[test]
    fn test_gammainc() {
        testing::test::<f64, _>("gammainc", "d_d-d", |x: &[f64]| gammainc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincc() {
        testing::test::<f64, _>("gammaincc", "d_d-d", |x: &[f64]| gammaincc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincinv() {
        testing::test::<f64, _>("gammaincinv", "d_d-d", |x: &[f64]| gammaincinv(x[0], x[1]));
    }

    #[test]
    fn test_gammainccinv() {
        testing::test::<f64, _>("gammainccinv", "d_d-d", |x: &[f64]| {
            gammainccinv(x[0], x[1])
        });
    }

    #[test]
    fn test_gammaln() {
        testing::test::<f64, _>("gammaln", "d-d", |x: &[f64]| gammaln(x[0]));
    }

    #[test]
    fn test_gammasgn() {
        testing::test::<f64, _>("gammasgn", "d-d", |x: &[f64]| gammasgn(x[0]));
    }
}
