use num_complex::Complex;

use crate::bindings;
use crate::bindings::xsf_impl;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait ZetaArg: sealed::Sealed {
    type Output;
    fn riemann_zeta(self) -> Self::Output;
    fn zeta(self, q: f64) -> Self::Output;
}

impl ZetaArg for f64 {
    type Output = f64;
    fn riemann_zeta(self) -> f64 {
        unsafe { bindings::riemann_zeta(self) }
    }
    fn zeta(self, q: f64) -> f64 {
        unsafe { bindings::zeta(self, q) }
    }
}

impl ZetaArg for Complex<f64> {
    type Output = Complex<f64>;
    fn riemann_zeta(self) -> Complex<f64> {
        unsafe { bindings::riemann_zeta_1(self.into()) }.into()
    }
    fn zeta(self, q: f64) -> Complex<f64> {
        unsafe { bindings::zeta_1(self.into(), q) }.into()
    }
}

/// Riemann zeta function for real or complex input
pub fn riemann_zeta<T: ZetaArg>(z: T) -> T::Output {
    z.riemann_zeta()
}

/// Riemann zeta function of two arguments
pub fn zeta<T: ZetaArg>(z: T, q: f64) -> T::Output {
    z.zeta(q)
}

xsf_impl!(zetac, (x: f64), "Riemann zeta function, minus one");
