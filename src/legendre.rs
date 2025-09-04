use crate::bindings;
use std::os::raw::c_int;

use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait LegendreArg: sealed::Sealed {
    type Output;
    fn legendre_p(self, n: c_int) -> Self::Output;
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self::Output;
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> Self::Output;
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> Self::Output;
}

impl LegendreArg for f64 {
    type Output = f64;
    fn legendre_p(self, n: c_int) -> f64 {
        unsafe { bindings::legendre_p(n, self) }
    }
    fn sph_legendre_p(self, n: c_int, m: c_int) -> f64 {
        unsafe { bindings::sph_legendre_p(n, m, self) }
    }
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> f64 {
        unsafe { bindings::assoc_legendre_p_0(n, m, self, branch_cut) }
    }
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> f64 {
        unsafe { bindings::assoc_legendre_p_1(n, m, self, branch_cut) }
    }
}

impl LegendreArg for Complex<f64> {
    type Output = Complex<f64>;
    fn legendre_p(self, n: c_int) -> Complex<f64> {
        unsafe { bindings::legendre_p_1(n, self.into()) }.into()
    }
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Complex<f64> {
        unsafe { bindings::sph_legendre_p_1(n, m, self.into()) }.into()
    }
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> Complex<f64> {
        unsafe { bindings::assoc_legendre_p_0_1(n, m, self.into(), branch_cut) }.into()
    }
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> Complex<f64> {
        unsafe { bindings::assoc_legendre_p_1_1(n, m, self.into(), branch_cut) }.into()
    }
}

/// Legendre polynomial of degree n
pub fn legendre_p<T: LegendreArg>(n: i32, z: T) -> T::Output {
    z.legendre_p(n as c_int)
}

/// Spherical Legendre polynomial of degree n and order m
pub fn sph_legendre_p<T: LegendreArg>(n: i32, m: i32, z: T) -> T::Output {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// Associated Legendre polynomial of the first kind
pub fn assoc_legendre_p<T: LegendreArg>(n: i32, m: i32, z: T, branch_cut: u8) -> T::Output {
    z.assoc_legendre_p(n as c_int, m as c_int, branch_cut as c_int)
}

/// Normalized associated Legendre polynomial of the first kind
pub fn assoc_legendre_p_norm<T: LegendreArg>(n: i32, m: i32, z: T, branch_cut: u8) -> T::Output {
    z.assoc_legendre_p_norm(n as c_int, m as c_int, branch_cut as c_int)
}
