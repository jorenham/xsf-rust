use core::ffi::c_long;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait SphBesselArg: sealed::Sealed {
    fn sph_bessel_j(self, n: c_long) -> Self;
    fn sph_bessel_y(self, n: c_long) -> Self;
    fn sph_bessel_i(self, n: c_long) -> Self;
    fn sph_bessel_k(self, n: c_long) -> Self;
    fn sph_bessel_j_prime(self, n: c_long) -> Self;
    fn sph_bessel_y_prime(self, n: c_long) -> Self;
    fn sph_bessel_i_prime(self, n: c_long) -> Self;
    fn sph_bessel_k_prime(self, n: c_long) -> Self;
}

impl SphBesselArg for f64 {
    #[inline(always)]
    fn sph_bessel_j(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_j(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_y(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_y(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_i(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_i(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_k(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_k(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_j_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_j_jac(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_y_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_y_jac(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_i_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_i_jac(n, self) }
    }
    #[inline(always)]
    fn sph_bessel_k_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_k_jac(n, self) }
    }
}

impl SphBesselArg for Complex<f64> {
    #[inline(always)]
    fn sph_bessel_j(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_j_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_y(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_y_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_i(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_i_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_k(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_k_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_j_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_j_jac_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_y_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_y_jac_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_i_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_i_jac_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_bessel_k_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_k_jac_1(n, self.into()) }.into()
    }
}

/// Spherical Bessel function of the first kind
#[doc(alias = "spherical_jn")]
pub fn sph_bessel_j<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_j(n as c_long)
}

/// Spherical Bessel function of the second kind
#[doc(alias = "spherical_yn")]
pub fn sph_bessel_y<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_y(n as c_long)
}

/// Modified spherical Bessel function of the first kind
#[doc(alias = "spherical_in")]
pub fn sph_bessel_i<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_i(n as c_long)
}

/// Modified spherical Bessel function of the second kind
#[doc(alias = "spherical_kn")]
pub fn sph_bessel_k<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_k(n as c_long)
}

/// Derivative of [`sph_bessel_j`]
#[doc(alias = "spherical_jnp")]
pub fn sph_bessel_j_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_j_prime(n as c_long)
}

/// Derivative of [`sph_bessel_y`]
#[doc(alias = "spherical_ynp")]
pub fn sph_bessel_y_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_y_prime(n as c_long)
}

/// Derivative of [`sph_bessel_i`]
#[doc(alias = "spherical_inp")]
pub fn sph_bessel_i_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_i_prime(n as c_long)
}

/// Derivative of [`sph_bessel_k`]
#[doc(alias = "spherical_knp")]
pub fn sph_bessel_k_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_k_prime(n as c_long)
}

#[cfg(test)]
mod tests {
    // TODO(@jorenham): jorenham/xsf-rust#91
}
