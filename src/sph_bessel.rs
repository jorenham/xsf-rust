use std::os::raw::c_long;

use crate::bindings;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait SphBesselArg: sealed::Sealed {
    type Output;
    fn sph_bessel_j(self, n: c_long) -> Self::Output;
    fn sph_bessel_y(self, n: c_long) -> Self::Output;
    fn sph_bessel_i(self, n: c_long) -> Self::Output;
    fn sph_bessel_k(self, n: c_long) -> Self::Output;
    fn sph_bessel_j_jac(self, n: c_long) -> Self::Output;
    fn sph_bessel_y_jac(self, n: c_long) -> Self::Output;
    fn sph_bessel_i_jac(self, n: c_long) -> Self::Output;
    fn sph_bessel_k_jac(self, n: c_long) -> Self::Output;
}

impl SphBesselArg for f64 {
    type Output = f64;
    fn sph_bessel_j(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_j(n, self) }
    }
    fn sph_bessel_y(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_y(n, self) }
    }
    fn sph_bessel_i(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_i(n, self) }
    }
    fn sph_bessel_k(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_k(n, self) }
    }
    fn sph_bessel_j_jac(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_j_jac(n, self) }
    }
    fn sph_bessel_y_jac(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_y_jac(n, self) }
    }
    fn sph_bessel_i_jac(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_i_jac(n, self) }
    }
    fn sph_bessel_k_jac(self, n: c_long) -> f64 {
        unsafe { bindings::sph_bessel_k_jac(n, self) }
    }
}

impl SphBesselArg for Complex<f64> {
    type Output = Complex<f64>;
    fn sph_bessel_j(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_j_1(n, self.into()) }.into()
    }
    fn sph_bessel_y(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_y_1(n, self.into()) }.into()
    }
    fn sph_bessel_i(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_i_1(n, self.into()) }.into()
    }
    fn sph_bessel_k(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_k_1(n, self.into()) }.into()
    }
    fn sph_bessel_j_jac(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_j_jac_1(n, self.into()) }.into()
    }
    fn sph_bessel_y_jac(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_y_jac_1(n, self.into()) }.into()
    }
    fn sph_bessel_i_jac(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_i_jac_1(n, self.into()) }.into()
    }
    fn sph_bessel_k_jac(self, n: c_long) -> Self::Output {
        unsafe { bindings::sph_bessel_k_jac_1(n, self.into()) }.into()
    }
}

/// Spherical Bessel function of the first kind
pub fn sph_bessel_j<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_j(n as c_long)
}

/// Spherical Bessel function of the second kind
pub fn sph_bessel_y<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_y(n as c_long)
}

/// Modified spherical Bessel function of the first kind
pub fn sph_bessel_i<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_i(n as c_long)
}

/// Modified spherical Bessel function of the second kind
pub fn sph_bessel_k<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_k(n as c_long)
}

/// Derivative of `sph_bessel_j` w.r.t `z`
pub fn sph_bessel_j_jac<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_j_jac(n as c_long)
}

/// Derivative of `sph_bessel_y` w.r.t `z`
pub fn sph_bessel_y_jac<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_y_jac(n as c_long)
}

/// Derivative of `sph_bessel_i` w.r.t `z`
pub fn sph_bessel_i_jac<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_i_jac(n as c_long)
}

/// Derivative of `sph_bessel_k` w.r.t `z`
pub fn sph_bessel_k_jac<T: SphBesselArg>(n: i64, z: T) -> T::Output {
    z.sph_bessel_k_jac(n as c_long)
}
