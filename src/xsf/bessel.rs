use core::ffi::c_int;
use num_complex::Complex;

mod sealed {
    use core::ops::{Add, Mul};

    // the additional bounds are needed for the arithmetic in `_bessel_diff_formula`
    pub trait Sealed: Copy + From<f64> + Add<Output = Self> + Mul<Output = Self> {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait BesselArg: sealed::Sealed {
    fn bessel_j(self, v: f64) -> Self;
    fn bessel_je(self, v: f64) -> Self;
    fn bessel_y(self, v: f64) -> Self;
    fn bessel_ye(self, v: f64) -> Self;
    fn bessel_i(self, v: f64) -> Self;
    fn bessel_ie(self, v: f64) -> Self;
    fn bessel_k(self, v: f64) -> Self;
    fn bessel_ke(self, v: f64) -> Self;
    fn hankel_1(self, v: f64) -> Complex<f64>;
    fn hankel_1e(self, v: f64) -> Complex<f64>;
    fn hankel_2(self, v: f64) -> Complex<f64>;
    fn hankel_2e(self, v: f64) -> Complex<f64>;
}

impl BesselArg for f64 {
    #[inline(always)]
    fn bessel_j(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_j(v, self) }
    }
    #[inline(always)]
    fn bessel_je(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_je(v, self) }
    }
    #[inline(always)]
    fn bessel_y(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_y(v, self) }
    }
    #[inline(always)]
    fn bessel_ye(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ye(v, self) }
    }
    #[inline(always)]
    fn bessel_i(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_i(v, self) }
    }
    #[inline(always)]
    fn bessel_ie(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ie(v, self) }
    }
    #[inline(always)]
    fn bessel_k(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_k(v, self) }
    }
    #[inline(always)]
    fn bessel_ke(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ke(v, self) }
    }
    #[inline(always)]
    fn hankel_1(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1(v, Complex::from(self)) }
    }
    #[inline(always)]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1e(v, Complex::from(self)) }
    }
    #[inline(always)]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2(v, Complex::from(self)) }
    }
    #[inline(always)]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2e(v, Complex::from(self)) }
    }
}

impl BesselArg for Complex<f64> {
    #[inline(always)]
    fn bessel_j(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_j_1(v, self) }
    }
    #[inline(always)]
    fn bessel_je(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_je_1(v, self) }
    }
    #[inline(always)]
    fn bessel_y(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_y_1(v, self) }
    }
    #[inline(always)]
    fn bessel_ye(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ye_1(v, self) }
    }
    #[inline(always)]
    fn bessel_i(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_i_1(v, self) }
    }
    #[inline(always)]
    fn bessel_ie(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ie_1(v, self) }
    }
    #[inline(always)]
    fn bessel_k(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_k_1(v, self) }
    }
    #[inline(always)]
    fn bessel_ke(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ke_1(v, self) }
    }
    #[inline(always)]
    fn hankel_1(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1(v, self) }
    }
    #[inline(always)]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1e(v, self) }
    }
    #[inline(always)]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2(v, self) }
    }
    #[inline(always)]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2e(v, self) }
    }
}

// Bessel J

/// Bessel function of the first kind of order 0
#[doc(alias = "j0")]
#[doc(alias = "cyl_bessel_j0")]
pub fn bessel_j0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_j0(x) }
}

/// Bessel function of the first kind of order 1
#[doc(alias = "j1")]
#[doc(alias = "cyl_bessel_j1")]
pub fn bessel_j1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_j1(x) }
}

/// Bessel function of the first kind
#[doc(alias = "jv")]
#[doc(alias = "cyl_bessel_j")]
pub fn bessel_j<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_j(v)
}

/// Exponentially scaled Bessel function of the first kind
#[doc(alias = "jve")]
#[doc(alias = "cyl_bessel_je")]
pub fn bessel_je<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_je(v)
}

// Bessel Y (Neumann)

/// Bessel function, second kind of order 0
#[doc(alias = "y0")]
#[doc(alias = "cyl_neumann_0")]
#[doc(alias = "cyl_bessel_y0")]
pub fn bessel_y0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_y0(x) }
}

/// Bessel function, second kind of order 1
#[doc(alias = "y1")]
#[doc(alias = "cyl_neumann_1")]
#[doc(alias = "cyl_bessel_y1")]
pub fn bessel_y1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_y1(x) }
}

/// Bessel function of the second kind
#[doc(alias = "yv")]
#[doc(alias = "cyl_neumann")]
#[doc(alias = "cyl_bessel_y")]
pub fn bessel_y<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_y(v)
}

/// Exponentially scaled Bessel function of the second kind
#[doc(alias = "yve")]
#[doc(alias = "cyl_neumann_e")]
#[doc(alias = "cyl_bessel_ye")]
pub fn bessel_ye<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ye(v)
}

// Bessel I

/// Modified Bessel function of the first kind of order 0
#[doc(alias = "i0")]
#[doc(alias = "cyl_bessel_i0")]
pub fn bessel_i0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i0(x) }
}

/// Exponentially scaled modified Bessel function of the first kind of order 0
#[doc(alias = "i0e")]
#[doc(alias = "cyl_bessel_i0e")]
pub fn bessel_i0e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i0e(x) }
}

/// Modified Bessel function of the first kind of order 1
#[doc(alias = "i1")]
#[doc(alias = "cyl_bessel_i1")]
pub fn bessel_i1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i1(x) }
}

/// Exponentially scaled modified Bessel function of the first kind of order 1
#[doc(alias = "i1e")]
#[doc(alias = "cyl_bessel_i1e")]
pub fn bessel_i1e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i1e(x) }
}

/// Modified Bessel function of the first kind
#[doc(alias = "iv")]
#[doc(alias = "cyl_bessel_i")]
pub fn bessel_i<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_i(v)
}

/// Exponentially scaled modified Bessel function of the first kind
#[doc(alias = "ive")]
#[doc(alias = "cyl_bessel_ie")]
pub fn bessel_ie<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ie(v)
}

// Bessel K

/// Modified Bessel function, second kind of order 0
#[doc(alias = "k0")]
#[doc(alias = "cyl_bessel_k0")]
pub fn bessel_k0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k0(x) }
}

/// Exponentially scaled modified Bessel function, second kind of order 0
#[doc(alias = "k0e")]
#[doc(alias = "cyl_bessel_k0e")]
pub fn bessel_k0e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k0e(x) }
}

/// Modified Bessel function, second kind of order 1
#[doc(alias = "k1")]
#[doc(alias = "cyl_bessel_k1")]
pub fn bessel_k1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k1(x) }
}

/// Exponentially scaled modified Bessel function, second kind of order 1
#[doc(alias = "k1e")]
#[doc(alias = "cyl_bessel_k1e")]
pub fn bessel_k1e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k1e(x) }
}

/// Modified Bessel function of the second kind
#[doc(alias = "kv")]
#[doc(alias = "cyl_bessel_k")]
pub fn bessel_k<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_k(v)
}

/// Exponentially scaled modified Bessel function of the second kind
#[doc(alias = "kve")]
#[doc(alias = "cyl_bessel_ke")]
pub fn bessel_ke<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ke(v)
}

// Zeros

/// Compute `NT` zeros of Bessel functions *J<sub>n</sub>(x)*, *J'<sub>n</sub>(x)*,
/// *Y<sub>n</sub>(x)*, and *Y'<sub>n</sub>(x)*
///
/// Corresponds to [`scipy.special.jnyn_zeros`][jnyn_zeros].
///
/// [jnyn_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnyn_zeros.html
#[doc(alias = "jnyn_zeros")]
pub fn bessel_zeros<const NT: usize>(n: u32) -> [[f64; NT]; 4] {
    assert!(NT <= 1200, "NT must be at most 1200");

    let mut jn = [f64::NAN; NT];
    let mut jnp = [f64::NAN; NT];
    let mut yn = [f64::NAN; NT];
    let mut ynp = [f64::NAN; NT];
    unsafe {
        crate::ffi::xsf::jyzo(
            n as c_int,
            NT as c_int,
            jn.as_mut_ptr(),
            jnp.as_mut_ptr(),
            yn.as_mut_ptr(),
            ynp.as_mut_ptr(),
        );
    }
    [jn, jnp, yn, ynp]
}

// Hankel 1

/// Hankel function of the first kind
#[doc(alias = "h1v")]
#[doc(alias = "hankel1")]
#[doc(alias = "cyl_hankel_1")]
pub fn hankel_1<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_1(v)
}

/// Exponentially scaled Hankel function of the first kind
#[doc(alias = "h1ve")]
#[doc(alias = "hankel1e")]
#[doc(alias = "cyl_hankel_1e")]
pub fn hankel_1e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_1e(v)
}

// Hankel 2

/// Hankel function of the second kind
#[doc(alias = "h2v")]
#[doc(alias = "hankel2")]
#[doc(alias = "cyl_hankel_2")]
pub fn hankel_2<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_2(v)
}

/// Exponentially scaled Hankel function of the second kind
#[doc(alias = "h2ve")]
#[doc(alias = "hankel2e")]
#[doc(alias = "cyl_hankel_2e")]
pub fn hankel_2e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_2e(v)
}

// Integrals

/// Weighted integral of the Bessel function of the first kind
///
/// Computes
///
/// *∫<sub>0</sub><sup>1</sup> x<sup>λ</sup>J<sub>ν</sub>(2ax) dx*
///
/// where *J<sub>ν</sub>* is a Bessel function and *λ*=`lambda`, *ν*=`nu`.
pub fn besselpoly(a: f64, lambda: f64, nu: f64) -> f64 {
    unsafe { crate::ffi::xsf::besselpoly(a, lambda, nu) }
}

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
#[doc(alias = "itj0y0")]
pub fn it1j0y0(x: f64) -> (f64, f64) {
    let (mut j0int, mut y0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it1j0y0(x, &mut j0int, &mut y0int);
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
    let (mut j0int, mut y0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it2j0y0(x, &mut j0int, &mut y0int);
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
#[doc(alias = "iti0k0")]
pub fn it1i0k0(x: f64) -> (f64, f64) {
    let (mut i0int, mut k0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it1i0k0(x, &mut i0int, &mut k0int);
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
    let (mut i0int, mut k0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it2i0k0(x, &mut i0int, &mut k0int);
    }
    (i0int, k0int)
}

// Riccati-Bessel

/// Compute Riccati-Bessel function of the first kind and derivatives for the first `NT` orders
///
/// The Riccati-Bessel function of the first kind is defined as `x j_n(x)`, where `j_n` is the
/// spherical Bessel function of the first kind of order `n`.
///
/// This function computes the value and first derivative of the function for all orders up to `NT`.
///
/// Corresponds to [`scipy.special.riccati_jn`][scipy-rctj].
///
/// # Returns
/// - [ *j<sub>0</sub>(x), ..., j<sub>N-1</sub>(x)* ]
/// - [ *j'<sub>0</sub>(x), ..., j'<sub>N-1</sub>(x)* ]
///
/// # See also
/// - [`riccati_y`]
/// - [`bessel_j`]
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
/// [scipy-rctj]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_jn.html
#[doc(alias = "rctj")]
#[doc(alias = "riccati_jn")]
pub fn riccati_j<const NT: usize>(x: f64) -> ([f64; NT], [f64; NT]) {
    let (mut rj, mut dj) = ([f64::NAN; NT], [f64::NAN; NT]);
    if NT == 0 {
        return (rj, dj);
    }

    let _ = unsafe { crate::ffi::xsf::rctj(NT, x, rj.as_mut_ptr(), dj.as_mut_ptr()) };
    (rj, dj)
}

/// Compute Riccati-Bessel function of the second kind and derivatives for the first `NT` orders
///
/// The Riccati-Bessel function of the second kind is defined here as `x y_n(x)`, where `y_n` is
/// the spherical Bessel function of the second kind of order `n`. *Note that this is in contrast
/// to a common convention that includes a minus sign in the definition.*
///
/// This function computes the value and first derivative of the function for all orders up to `NT`.
///
/// Corresponds to [`scipy.special.riccati_yn`][scipy-rcty].
///
/// # Returns
/// - [ *y<sub>0</sub>(x), ..., y<sub>N-1</sub>(x)* ]
/// - [ *y'<sub>0</sub>(x), ..., y'<sub>N-1</sub>(x)* ]
///
/// # See also
/// - [`riccati_j`]
/// - [`bessel_y`](crate::bessel_y)
/// - [`sph_bessel_y`](crate::sph_bessel_y)
///
/// [scipy-rcty]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_yn.html
#[doc(alias = "rcty")]
#[doc(alias = "riccati_yn")]
pub fn riccati_y<const NT: usize>(x: f64) -> ([f64; NT], [f64; NT]) {
    let (mut ry, mut dy) = ([f64::NAN; NT], [f64::NAN; NT]);
    if NT == 0 {
        return (ry, dy);
    }

    let _ = unsafe { crate::ffi::xsf::rcty(NT, x, ry.as_mut_ptr(), dy.as_mut_ptr()) };
    (ry, dy)
}

/// Evaluate the Jahnke-Emden Lambda function Λ<sub>v</sub>(x) and its derivatives
///
/// The Jahnke-Emden Lambda function is defined as
///
/// Λ<sub>v</sub>(x) = Γ(v+1) (2/x)<sup>v</sup> J<sub>v</sub>(x)
///
/// where Γ is the [Gamma function](crate::gamma) and J<sub>v</sub> is the
/// [Bessel function of the first kind](crate::bessel_j).
///
/// Corresponds to [`scipy.special.lmbda`][scipy-lmbda] in SciPy, and calls the FFI functions
/// `xsf::specfun::lamn` and `xsf::specfun::lamv`.
///
/// [scipy-lmbda]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lmbda.html
///
/// # Arguments
/// - `v`: Order of the Lambda function as a non-negative real number
/// - `x`: Value at which to evaluate the function and derivatives
///
/// # Returns
/// - `vl`: Values of Λ<sub>v<sub>i</sub></sub>(x) for v<sub>i</sub> = v-⌊v⌋, v-⌊v⌋+1, ..., v
/// - `dl`: Derivatives Λ'<sub>v<sub>i</sub></sub>(x) for v<sub>i</sub> = v-⌊v⌋, v-⌊v⌋+1, ..., v
///
/// Both vectors have length ⌊v⌋+1, i.e. `(v as usize) + 1`.
///
/// # See also
/// - [`bessel_j`]
/// - [`gamma`](crate::gamma)
/// - [`sph_bessel_j`](crate::sph_bessel_j)
#[doc(alias = "lmbda")]
pub fn jahnke_emden_lambda<V: Into<f64>>(v: V, x: f64) -> (Vec<f64>, Vec<f64>) {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L2008-L2055
    let v: f64 = v.into();
    assert!(v >= 0.0, "v must be non-negative");

    let n = v as usize;
    let v0 = v - n as f64;
    let n1 = if n < 1 { 1 } else { n };
    let v1 = v0 + n1 as f64;

    let size = v1 as usize + 1;
    let (mut vl, mut dl) = (vec![f64::NAN; size], vec![f64::NAN; size]);
    unsafe {
        if v != n as f64 {
            crate::ffi::xsf::lamv(v1, x, vl.as_mut_ptr(), dl.as_mut_ptr());
        } else {
            crate::ffi::xsf::lamn(v1 as c_int, x, vl.as_mut_ptr(), dl.as_mut_ptr());
        }
    };

    vl.resize(n + 1, f64::NAN);
    dl.resize(n + 1, f64::NAN);
    (vl, dl)
}

// Tests

#[cfg(test)]
mod tests {
    use num_complex::c64;

    // Bessel J

    #[test]
    fn test_bessel_j0() {
        crate::xsref::test("cyl_bessel_j0", "d-d", |x| crate::bessel_j0(x[0]));
    }

    #[test]
    fn test_bessel_j1() {
        crate::xsref::test("cyl_bessel_j1", "d-d", |x| crate::bessel_j1(x[0]));
    }

    #[test]
    fn test_bessel_j_f64() {
        crate::xsref::test("cyl_bessel_j", "d_d-d", |x| crate::bessel_j(x[0], x[1]));
    }

    #[test]
    fn test_bessel_j_c64() {
        crate::xsref::test("cyl_bessel_j", "d_cd-cd", |x| {
            crate::bessel_j(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_je_f64() {
        crate::xsref::test("cyl_bessel_je", "d_d-d", |x| crate::bessel_je(x[0], x[1]));
    }

    #[test]
    fn test_bessel_je_c64() {
        crate::xsref::test("cyl_bessel_je", "d_cd-cd", |x| {
            crate::bessel_je(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel Y

    #[test]
    fn test_bessel_y0() {
        crate::xsref::test("cyl_bessel_y0", "d-d", |x| crate::bessel_y0(x[0]));
    }

    #[test]
    fn test_bessel_y1() {
        crate::xsref::test("cyl_bessel_y1", "d-d", |x| crate::bessel_y1(x[0]));
    }

    #[test]
    fn test_bessel_y_f64() {
        crate::xsref::test("cyl_bessel_y", "d_d-d", |x| crate::bessel_y(x[0], x[1]));
    }

    #[test]
    fn test_bessel_y_c64() {
        crate::xsref::test("cyl_bessel_y", "d_cd-cd", |x| {
            crate::bessel_y(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ye_f64() {
        crate::xsref::test("cyl_bessel_ye", "d_d-d", |x| crate::bessel_ye(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ye_c64() {
        crate::xsref::test("cyl_bessel_ye", "d_cd-cd", |x| {
            crate::bessel_ye(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel I

    #[test]
    fn test_bessel_i0() {
        crate::xsref::test("cyl_bessel_i0", "d-d", |x| crate::bessel_i0(x[0]));
    }

    #[test]
    fn test_bessel_i0e() {
        crate::xsref::test("cyl_bessel_i0e", "d-d", |x| crate::bessel_i0e(x[0]));
    }

    #[test]
    fn test_bessel_i1() {
        crate::xsref::test("cyl_bessel_i1", "d-d", |x| crate::bessel_i1(x[0]));
    }

    #[test]
    fn test_bessel_i1e() {
        crate::xsref::test("cyl_bessel_i1e", "d-d", |x| crate::bessel_i1e(x[0]));
    }

    #[test]
    fn test_bessel_i_f64() {
        crate::xsref::test("cyl_bessel_i", "d_d-d", |x| crate::bessel_i(x[0], x[1]));
    }

    #[test]
    fn test_bessel_i_c64() {
        crate::xsref::test("cyl_bessel_i", "d_cd-cd", |x| {
            crate::bessel_i(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ie_f64() {
        crate::xsref::test("cyl_bessel_ie", "d_d-d", |x| crate::bessel_ie(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ie_c64() {
        crate::xsref::test("cyl_bessel_ie", "d_cd-cd", |x| {
            crate::bessel_ie(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel K

    #[test]
    fn test_bessel_k0() {
        crate::xsref::test("cyl_bessel_k0", "d-d", |x| crate::bessel_k0(x[0]));
    }

    #[test]
    fn test_bessel_k0e() {
        crate::xsref::test("cyl_bessel_k0e", "d-d", |x| crate::bessel_k0e(x[0]));
    }

    #[test]
    fn test_bessel_k1() {
        crate::xsref::test("cyl_bessel_k1", "d-d", |x| crate::bessel_k1(x[0]));
    }

    #[test]
    fn test_bessel_k1e() {
        crate::xsref::test("cyl_bessel_k1e", "d-d", |x| crate::bessel_k1e(x[0]));
    }

    #[test]
    fn test_bessel_k_f64() {
        crate::xsref::test("cyl_bessel_k", "d_d-d", |x| crate::bessel_k(x[0], x[1]));
    }

    #[test]
    fn test_bessel_k_c64() {
        crate::xsref::test("cyl_bessel_k", "d_cd-cd", |x| {
            crate::bessel_k(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ke_f64() {
        crate::xsref::test("cyl_bessel_ke", "d_d-d", |x| crate::bessel_ke(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ke_c64() {
        crate::xsref::test("cyl_bessel_ke", "d_cd-cd", |x| {
            crate::bessel_ke(x[0], c64(x[1], x[2]))
        });
    }

    // Hankel

    #[test]
    fn test_hankel_1_c64() {
        crate::xsref::test("cyl_hankel_1", "d_cd-cd", |x| {
            crate::hankel_1(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_1e_c64() {
        crate::xsref::test("cyl_hankel_1e", "d_cd-cd", |x| {
            crate::hankel_1e(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_2_c64() {
        crate::xsref::test("cyl_hankel_2", "d_cd-cd", |x| {
            crate::hankel_2(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_2e_c64() {
        crate::xsref::test("cyl_hankel_2e", "d_cd-cd", |x| {
            crate::hankel_2e(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel zeros

    /// Based on `scipy.special.tests.test_basic.TestBessel.test_jnyn_zeros`
    #[test]
    fn test_bessel_zeros() {
        // jnz = special.jnyn_zeros(1, 5)
        let jnz = crate::bessel_zeros::<5>(1);
        // assert_allclose(jnz, (array([3.83171,
        //                              7.01559,
        //                              10.17347,
        //                              13.32369,
        //                              16.47063]),
        //                       array([1.84118,
        //                              5.33144,
        //                              8.53632,
        //                              11.70600,
        //                              14.86359]),
        //                       array([2.19714,
        //                              5.42968,
        //                              8.59601,
        //                              11.74915,
        //                              14.89744]),
        //                       array([3.68302,
        //                              6.94150,
        //                              10.12340,
        //                              13.28576,
        //                              16.44006])),
        //                 atol=1.5e-5, rtol=0)
        crate::np_assert_allclose!(
            jnz[0].as_ref(),
            &[3.83171, 7.01559, 10.17347, 13.32369, 16.47063],
            atol = 1.5e-5
        );
        crate::np_assert_allclose!(
            jnz[1].as_ref(),
            &[1.84118, 5.33144, 8.53632, 11.70600, 14.86359],
            atol = 1.5e-5
        );
        crate::np_assert_allclose!(
            jnz[2].as_ref(),
            &[2.19714, 5.42968, 8.59601, 11.74915, 14.89744],
            atol = 1.5e-5
        );
        crate::np_assert_allclose!(
            jnz[3].as_ref(),
            &[3.68302, 6.94150, 10.12340, 13.28576, 16.44006],
            atol = 1.5e-5
        );
    }

    // Bessel integrals

    #[test]
    fn test_it1j0y0() {
        crate::xsref::test("it1j0y0", "d-d_d", |x| crate::it1j0y0(x[0]));
    }

    #[test]
    fn test_it2j0y0() {
        crate::xsref::test("it2j0y0", "d-d_d", |x| crate::it2j0y0(x[0]));
    }

    #[test]
    fn test_it1i0k0() {
        crate::xsref::test("it1i0k0", "d-d_d", |x| crate::it1i0k0(x[0]));
    }

    #[test]
    fn test_it2i0k0() {
        crate::xsref::test("it2i0k0", "d-d_d", |x| crate::it2i0k0(x[0]));
    }

    #[test]
    fn test_besselpoly() {
        crate::xsref::test("besselpoly", "d_d_d-d", |x| {
            crate::besselpoly(x[0], x[1], x[2])
        });
    }

    // Riccati-Bessel

    /// Based on `scipy.special.tests.test_basic.TestRiccati.test_riccati_jn`
    #[test]
    fn test_riccati_j() {
        // N, x = 2, 0.2
        const N: usize = 2;
        let x: f64 = 0.2;
        // S = np.empty((N, N))
        let mut s = [[0.0; N], [0.0; N]];
        // for n in range(N):
        #[allow(clippy::needless_range_loop)]
        for n in 0..N {
            // j = special.spherical_jn(n, x)
            let j = crate::sph_bessel_j(n as i64, x);
            // jp = special.spherical_jn(n, x, derivative=True)
            let jp = crate::sph_bessel_j_prime(n as i64, x);
            // S[0,n] = x*j
            s[0][n] = x * j;
            // S[1,n] = x*jp + j
            s[1][n] = x * jp + j;
        }

        // assert_allclose(S, special.riccati_jn(n, x), atol=1.5e-8, rtol=0)
        let (jn, jnp) = crate::riccati_j::<N>(x);
        crate::np_assert_allclose!(&s[0], &jn, atol = 1.5e-8);
        crate::np_assert_allclose!(&s[1], &jnp, atol = 1.5e-8);
    }

    /// Based on `scipy.special.tests.test_basic.TestRiccati.test_riccati_yn`
    #[test]
    fn test_riccati_y() {
        // N, x = 2, 0.2
        const N: usize = 2;
        let x: f64 = 0.2;
        // S = np.empty((N, N))
        let mut c = [[0.0; N], [0.0; N]];
        // for n in range(N):
        #[allow(clippy::needless_range_loop)]
        for n in 0..N {
            // j = special.spherical_jn(n, x)
            let j = crate::sph_bessel_y(n as i64, x);
            // jp = special.spherical_jn(n, x, derivative=True)
            let jp = crate::sph_bessel_y_prime(n as i64, x);
            // S[0,n] = x*j
            c[0][n] = x * j;
            // S[1,n] = x*jp + j
            c[1][n] = x * jp + j;
        }

        // assert_allclose(S, special.riccati_jn(n, x), atol=1.5e-8, rtol=0)
        let (yn, ynp) = crate::riccati_y::<N>(x);
        crate::np_assert_allclose!(&c[0], &yn, atol = 1.5e-8);
        crate::np_assert_allclose!(&c[1], &ynp, atol = 1.5e-8);
    }

    // Based on `scipy.special.tests.test_basic.TestLambda.test_lmbda`
    #[test]
    fn test_jahnke_emden_lambda() {
        // lam = special.lmbda(1,.1)
        let lam = crate::jahnke_emden_lambda(1, 0.1);
        // lamr = (
        //     array([special.jn(0,.1), 2*special.jn(1,.1)/.1]),
        //     array([special.jvp(0,.1), -2*special.jv(1,.1)/.01 + 2*special.jvp(1,.1)/.1])
        // )
        let lamr = (
            [crate::bessel_j(0.0, 0.1), 2e1 * crate::bessel_j(1.0, 0.1)],
            [
                crate::bessel_j_prime(0.0, 0.1, 1),
                -2e2 * crate::bessel_j(1.0, 0.1) + 2e1 * crate::bessel_j_prime(1.0, 0.1, 1),
            ],
        );
        // assert_allclose(lam, lamr, atol=1.5e-8, rtol=0)
        crate::np_assert_allclose!(&lam.0, &lamr.0, atol = 1.5e-8);
        crate::np_assert_allclose!(&lam.1, &lamr.1, atol = 1.5e-8);
    }
}
