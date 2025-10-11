use crate::ffi;
use crate::ffi::xsf_impl;
use alloc::vec::Vec;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait BesselArg: sealed::Sealed {
    fn cyl_bessel_j(self, v: f64) -> Self;
    fn cyl_bessel_je(self, v: f64) -> Self;
    fn cyl_bessel_y(self, v: f64) -> Self;
    fn cyl_bessel_ye(self, v: f64) -> Self;
    fn cyl_bessel_i(self, v: f64) -> Self;
    fn cyl_bessel_ie(self, v: f64) -> Self;
    fn cyl_bessel_k(self, v: f64) -> Self;
    fn cyl_bessel_ke(self, v: f64) -> Self;
}

impl BesselArg for f64 {
    #[inline(always)]
    fn cyl_bessel_j(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_j(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_je(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_je(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_y(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_y(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_ye(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_ye(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_i(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_i(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_ie(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_ie(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_k(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_k(v, self) }
    }
    #[inline(always)]
    fn cyl_bessel_ke(self, v: f64) -> f64 {
        unsafe { ffi::cyl_bessel_ke(v, self) }
    }
}

impl BesselArg for Complex<f64> {
    #[inline(always)]
    fn cyl_bessel_j(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_j_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_je(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_je_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_y(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_y_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_ye(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_ye_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_i(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_i_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_ie(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_ie_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_k(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_k_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_bessel_ke(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_bessel_ke_1(v, self.into()) }.into()
    }
}

pub trait HankelArg: sealed::Sealed {
    fn cyl_hankel_1(self, v: f64) -> Complex<f64>;
    fn cyl_hankel_1e(self, v: f64) -> Complex<f64>;
    fn cyl_hankel_2(self, v: f64) -> Complex<f64>;
    fn cyl_hankel_2e(self, v: f64) -> Complex<f64>;
}

impl HankelArg for f64 {
    #[inline(always)]
    fn cyl_hankel_1(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_1(v, Complex::new(self, 0.0).into()) }.into()
    }
    #[inline(always)]
    fn cyl_hankel_1e(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_1e(v, Complex::new(self, 0.0).into()) }.into()
    }
    #[inline(always)]
    fn cyl_hankel_2(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_2(v, Complex::new(self, 0.0).into()) }.into()
    }
    #[inline(always)]
    fn cyl_hankel_2e(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_2e(v, Complex::new(self, 0.0).into()) }.into()
    }
}

impl HankelArg for Complex<f64> {
    #[inline(always)]
    fn cyl_hankel_1(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_1(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_hankel_1e(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_1e(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_hankel_2(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_2(v, self.into()) }.into()
    }
    #[inline(always)]
    fn cyl_hankel_2e(self, v: f64) -> Complex<f64> {
        unsafe { ffi::cyl_hankel_2e(v, self.into()) }.into()
    }
}

xsf_impl!(cyl_bessel_j0, (x: f64), "Bessel function, 1st kind, order 0");
xsf_impl!(cyl_bessel_j1, (x: f64), "Bessel function, 1st kind, order 1");
xsf_impl!(cyl_bessel_y0, (x: f64), "Bessel function, 2nd kind, order 0");
xsf_impl!(cyl_bessel_y1, (x: f64), "Bessel function, 2nd kind, order 1");
xsf_impl!(cyl_bessel_i0, (x: f64), "Modified Bessel function, 1st kind, order 0");
xsf_impl!(
    cyl_bessel_i0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 0"
);
xsf_impl!(cyl_bessel_i1, (x: f64), "Modified Bessel function, 1st kind, order 1");
xsf_impl!(
    cyl_bessel_i1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 1"
);
xsf_impl!(cyl_bessel_k0, (x: f64), "Modified Bessel function, 2nd kind, order 0");
xsf_impl!(
    cyl_bessel_k0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 0"
);
xsf_impl!(cyl_bessel_k1, (x: f64), "Modified Bessel function, 2nd kind, order 1");
xsf_impl!(
    cyl_bessel_k1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 1"
);

/// Bessel function, 1st kind
pub fn cyl_bessel_j<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_j(v)
}

/// Exponentially scaled Bessel function, 1st kind
pub fn cyl_bessel_je<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_je(v)
}

/// Bessel function, 2nd kind
pub fn cyl_bessel_y<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_y(v)
}

/// Exponentially scaled Bessel function, 2nd kind
pub fn cyl_bessel_ye<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_ye(v)
}

/// Modified Bessel function, 1st kind
pub fn cyl_bessel_i<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_i(v)
}

/// Exponentially scaled modified Bessel function, 1st kind
pub fn cyl_bessel_ie<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_ie(v)
}

/// Modified Bessel function, 2nd kind
pub fn cyl_bessel_k<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_k(v)
}

/// Exponentially scaled modified Bessel function, 2nd kind
pub fn cyl_bessel_ke<T: BesselArg>(v: f64, x: T) -> T {
    x.cyl_bessel_ke(v)
}

/// Hankel function, 1st kind
pub fn cyl_hankel_1<T: HankelArg>(v: f64, z: T) -> Complex<f64> {
    z.cyl_hankel_1(v)
}

/// Exponentially scaled Hankel function, 1st kind
pub fn cyl_hankel_1e<T: HankelArg>(v: f64, z: T) -> Complex<f64> {
    z.cyl_hankel_1e(v)
}

/// Hankel function, 2nd kind
pub fn cyl_hankel_2<T: HankelArg>(v: f64, z: T) -> Complex<f64> {
    z.cyl_hankel_2(v)
}

/// Exponentially scaled Hankel function, 2nd kind
pub fn cyl_hankel_2e<T: HankelArg>(v: f64, z: T) -> Complex<f64> {
    z.cyl_hankel_2e(v)
}

xsf_impl!(
    besselpoly,
    (a: f64, lambda: f64, nu: f64),
    "Weighted integral of the Bessel function of the first kind"
);

/// Integrals of Bessel functions of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ J₀(t) dt
/// - ∫₀ˣ Y₀(t) dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `j0int`: Integral for J₀
/// - `y0int`: Integral for Y₀
pub fn it1j0y0(x: f64) -> (f64, f64) {
    let mut j0int = f64::NAN;
    let mut y0int = f64::NAN;
    unsafe {
        ffi::it1j0y0(x, &mut j0int, &mut y0int);
    }
    (j0int, y0int)
}

/// Integrals related to Bessel functions of the first kind of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ (1 - J₀(t))/t dt
/// - ∫ₓ∞ Y₀(t)/t dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `j0int`: Integral for J₀
/// - `y0int`: Integral for Y₀
pub fn it2j0y0(x: f64) -> (f64, f64) {
    let mut j0int = f64::NAN;
    let mut y0int = f64::NAN;
    unsafe {
        ffi::it2j0y0(x, &mut j0int, &mut y0int);
    }
    (j0int, y0int)
}

/// Integrals of modified Bessel functions of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ I₀(t) dt
/// - ∫₀ˣ K₀(t) dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `i0int`: The integral for I₀
/// - `k0int`: The integral for K₀
pub fn it1i0k0(x: f64) -> (f64, f64) {
    let mut i0int = f64::NAN;
    let mut k0int = f64::NAN;
    unsafe {
        ffi::it1i0k0(x, &mut i0int, &mut k0int);
    }
    (i0int, k0int)
}

/// Integrals related to modified Bessel functions of order 0 up to `n`
///
/// Computes the integrals:
///
/// - ∫₀ˣ (I₀(t) - 1)/t dt
/// - ∫ₓ∞ K₀(t)/t dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `i0int`: The integral for I₀
/// - `k0int`: The integral for K₀
pub fn it2i0k0(x: f64) -> (f64, f64) {
    let mut i0int = f64::NAN;
    let mut k0int = f64::NAN;
    unsafe {
        ffi::it2i0k0(x, &mut i0int, &mut k0int);
    }
    (i0int, k0int)
}

/// Compute Ricatti-Bessel function of the 1st kind and their derivatives for up to `n`
///
/// The Ricatti-Bessel function of the first kind is defined as `x j_n(x)`, where `j_n` is the
/// spherical Bessel function of the first kind of order `n`.
///
/// This function computes the value and first derivative of the
/// Ricatti-Bessel function for all orders up to and including `n`.
///
/// # Arguments
/// - `n` - Maximum order of function to compute
/// - `x` - Argument at which to evaluate
///
/// # Returns
/// - `jn`: Value of *j0(x), ..., jn(x)*
/// - `jnp`:  First derivative *j0'(x), ..., jn'(x)*
/// - `nmax`: Highest order computed
pub fn rctj(n: usize, x: f64) -> (Vec<f64>, Vec<f64>, i32) {
    let mut rj = alloc::vec![0.0; n + 1];
    let mut dj = alloc::vec![0.0; n + 1];
    let nm = unsafe { ffi::rctj(n, x, rj.as_mut_ptr(), dj.as_mut_ptr()) } as i32;
    (rj, dj, nm)
}

/// Compute Ricatti-Bessel function of the 2nd kind and their derivatives for up to `n`
///
/// The Ricatti-Bessel function of the second kind is defined here as `+x y_n(x)`, where `y_n` is
/// the spherical Bessel function of the second kind of order `n`. *Note that this is in contrast
/// to a common convention that includes a minus sign in the definition.*
///
/// This function computes the value and first derivative of the function for all orders up to and
/// including `n`.
///
/// # Arguments
/// - `n` - Maximum order of function to compute
/// - `x` - Argument at which to evaluate
///
/// # Returns
/// - `yn`: Value of y0(x), ..., yn(x)
/// - `ynp`:  First derivative y0'(x), ..., yn'(x)
/// - `nmax`: Highest order computed
pub fn rcty(n: usize, x: f64) -> (Vec<f64>, Vec<f64>, i32) {
    let mut ry = alloc::vec![0.0; n + 1];
    let mut dy = alloc::vec![0.0; n + 1];
    let nm = unsafe { ffi::rcty(n, x, ry.as_mut_ptr(), dy.as_mut_ptr()) } as i32;
    (ry, dy, nm as i32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    // bessel j

    #[test]
    fn test_cyl_bessel_j0() {
        testing::test::<f64, _>("cyl_bessel_j0", "d-d", |x: &[f64]| cyl_bessel_j0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_j1() {
        testing::test::<f64, _>("cyl_bessel_j1", "d-d", |x: &[f64]| cyl_bessel_j1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_j_f64() {
        testing::test::<f64, _>("cyl_bessel_j", "d_d-d", |x: &[f64]| {
            cyl_bessel_j(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_j_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_j", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_j(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_je_f64() {
        testing::test::<f64, _>("cyl_bessel_je", "d_d-d", |x: &[f64]| {
            cyl_bessel_je(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_je_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_je", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_je(x[0], c64(x[1], x[2]))
        });
    }

    // bessel y

    #[test]
    fn test_cyl_bessel_y0() {
        testing::test::<f64, _>("cyl_bessel_y0", "d-d", |x: &[f64]| cyl_bessel_y0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_y1() {
        testing::test::<f64, _>("cyl_bessel_y1", "d-d", |x: &[f64]| cyl_bessel_y1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_y_f64() {
        testing::test::<f64, _>("cyl_bessel_y", "d_d-d", |x: &[f64]| {
            cyl_bessel_y(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_y_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_y", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_y(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_ye_f64() {
        testing::test::<f64, _>("cyl_bessel_ye", "d_d-d", |x: &[f64]| {
            cyl_bessel_ye(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_ye_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_ye", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_ye(x[0], c64(x[1], x[2]))
        });
    }

    // bessel i

    #[test]
    fn test_cyl_bessel_i0() {
        testing::test::<f64, _>("cyl_bessel_i0", "d-d", |x: &[f64]| cyl_bessel_i0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i0e() {
        testing::test::<f64, _>("cyl_bessel_i0e", "d-d", |x: &[f64]| cyl_bessel_i0e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i1() {
        testing::test::<f64, _>("cyl_bessel_i1", "d-d", |x: &[f64]| cyl_bessel_i1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i1e() {
        testing::test::<f64, _>("cyl_bessel_i1e", "d-d", |x: &[f64]| cyl_bessel_i1e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i_f64() {
        testing::test::<f64, _>("cyl_bessel_i", "d_d-d", |x: &[f64]| {
            cyl_bessel_i(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_i_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_i", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_i(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_ie_f64() {
        testing::test::<f64, _>("cyl_bessel_ie", "d_d-d", |x: &[f64]| {
            cyl_bessel_ie(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_ie_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_ie", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_ie(x[0], c64(x[1], x[2]))
        });
    }

    // bessel k

    #[test]
    fn test_cyl_bessel_k0() {
        testing::test::<f64, _>("cyl_bessel_k0", "d-d", |x: &[f64]| cyl_bessel_k0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k0e() {
        testing::test::<f64, _>("cyl_bessel_k0e", "d-d", |x: &[f64]| cyl_bessel_k0e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k1() {
        testing::test::<f64, _>("cyl_bessel_k1", "d-d", |x: &[f64]| cyl_bessel_k1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k1e() {
        testing::test::<f64, _>("cyl_bessel_k1e", "d-d", |x: &[f64]| cyl_bessel_k1e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k_f64() {
        testing::test::<f64, _>("cyl_bessel_k", "d_d-d", |x: &[f64]| {
            cyl_bessel_k(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_k_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_k", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_k(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_ke_f64() {
        testing::test::<f64, _>("cyl_bessel_ke", "d_d-d", |x: &[f64]| {
            cyl_bessel_ke(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_ke_c64() {
        testing::test::<Complex<f64>, _>("cyl_bessel_ke", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_ke(x[0], c64(x[1], x[2]))
        });
    }

    // hankel

    #[test]
    fn test_cyl_hankel_1_c64() {
        testing::test::<Complex<f64>, _>("cyl_hankel_1", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_1(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_hankel_1e_c64() {
        testing::test::<Complex<f64>, _>("cyl_hankel_1e", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_1e(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_hankel_2_c64() {
        testing::test::<Complex<f64>, _>("cyl_hankel_2", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_2(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_hankel_2e_c64() {
        testing::test::<Complex<f64>, _>("cyl_hankel_2e", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_2e(x[0], c64(x[1], x[2]))
        });
    }

    // besselpoly

    #[test]
    fn test_besselpoly_f64() {
        testing::test::<f64, _>("besselpoly", "d_d_d-d", |x: &[f64]| {
            besselpoly(x[0], x[1], x[2])
        });
    }

    // bessel integrals

    #[test]
    fn test_it1j0y0_f64() {
        testing::test::<(f64, f64), _>("it1j0y0", "d-d_d", |x: &[f64]| it1j0y0(x[0]));
    }

    #[test]
    fn test_it2j0y0_f64() {
        testing::test::<(f64, f64), _>("it2j0y0", "d-d_d", |x: &[f64]| it2j0y0(x[0]));
    }

    #[test]
    fn test_it1i0k0_f64() {
        testing::test::<(f64, f64), _>("it1i0k0", "d-d_d", |x: &[f64]| it1i0k0(x[0]));
    }

    #[test]
    fn test_it2i0k0_f64() {
        testing::test::<(f64, f64), _>("it2i0k0", "d-d_d", |x: &[f64]| it2i0k0(x[0]));
    }
}
