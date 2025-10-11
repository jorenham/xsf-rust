use num_complex::Complex;

use crate::ffi;
use crate::ffi::xsf_impl;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ZetaArg: sealed::Sealed {
    fn riemann_zeta(self) -> Self;
    fn zeta(self, q: f64) -> Self;
}

impl ZetaArg for f64 {
    #[inline(always)]
    fn riemann_zeta(self) -> f64 {
        unsafe { ffi::riemann_zeta(self) }
    }
    #[inline(always)]
    fn zeta(self, q: f64) -> f64 {
        unsafe { ffi::zeta(self, q) }
    }
}

impl ZetaArg for Complex<f64> {
    #[inline(always)]
    fn riemann_zeta(self) -> Complex<f64> {
        unsafe { ffi::riemann_zeta_1(self.into()) }.into()
    }
    #[inline(always)]
    fn zeta(self, q: f64) -> Complex<f64> {
        unsafe { ffi::zeta_1(self.into(), q) }.into()
    }
}

/// Riemann zeta function for real or complex input
pub fn riemann_zeta<T: ZetaArg>(z: T) -> T {
    z.riemann_zeta()
}

/// Riemann zeta function of two arguments for real or complex `z`
pub fn zeta<T: ZetaArg>(z: T, q: f64) -> T {
    z.zeta(q)
}

xsf_impl!(zetac, (x: f64), "Riemann zeta function, minus one");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    #[test]
    fn test_riemann_zeta_f64() {
        testing::test::<f64, _>("riemann_zeta", "d-d", |x: &[f64]| riemann_zeta(x[0]));
    }

    #[test]
    fn test_riemann_zeta_c64() {
        testing::test::<Complex<f64>, _>("riemann_zeta", "cd-cd", |x: &[f64]| {
            riemann_zeta(c64(x[0], x[1]))
        });
    }

    #[test]
    fn test_zeta_f64() {
        testing::test::<f64, _>("zeta", "d_d-d", |x: &[f64]| zeta(x[0], x[1]));
    }

    #[test]
    fn test_zetac_f64() {
        testing::test::<f64, _>("zetac", "d-d", |x: &[f64]| zetac(x[0]));
    }
}
