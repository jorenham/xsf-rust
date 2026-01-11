use num_complex::Complex;
use num_traits::ToPrimitive;

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
    #[inline]
    fn bessel_j(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_j(v, self) }
    }
    #[inline]
    fn bessel_je(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_je(v, self) }
    }
    #[inline]
    fn bessel_y(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_y(v, self) }
    }
    #[inline]
    fn bessel_ye(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ye(v, self) }
    }
    #[inline]
    fn bessel_i(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_i(v, self) }
    }
    #[inline]
    fn bessel_ie(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ie(v, self) }
    }
    #[inline]
    fn bessel_k(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_k(v, self) }
    }
    #[inline]
    fn bessel_ke(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ke(v, self) }
    }
    #[inline]
    fn hankel_1(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1(v, self.into()) }
    }
    #[inline]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1e(v, self.into()) }
    }
    #[inline]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2(v, self.into()) }
    }
    #[inline]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2e(v, self.into()) }
    }
}

impl BesselArg for Complex<f64> {
    #[inline]
    fn bessel_j(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_j_1(v, self) }
    }
    #[inline]
    fn bessel_je(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_je_1(v, self) }
    }
    #[inline]
    fn bessel_y(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_y_1(v, self) }
    }
    #[inline]
    fn bessel_ye(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ye_1(v, self) }
    }
    #[inline]
    fn bessel_i(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_i_1(v, self) }
    }
    #[inline]
    fn bessel_ie(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ie_1(v, self) }
    }
    #[inline]
    fn bessel_k(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_k_1(v, self) }
    }
    #[inline]
    fn bessel_ke(self, v: f64) -> Self {
        unsafe { crate::ffi::xsf::cyl_bessel_ke_1(v, self) }
    }
    #[inline]
    fn hankel_1(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1(v, self) }
    }
    #[inline]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_1e(v, self) }
    }
    #[inline]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2(v, self) }
    }
    #[inline]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        unsafe { crate::ffi::xsf::cyl_hankel_2e(v, self) }
    }
}

// Bessel J

/// Bessel function of the first kind, $J_v(z)$
///
/// Corresponds to [`scipy.special.jv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html
///
/// # See also
/// - [`bessel_j0`]: faster version of $J_0$ for real argument
/// - [`bessel_j1`]: faster version of $J_1$ for real argument
/// - [`bessel_je`]: $J_v$ without the leading exponential factor
/// - [`bessel_j_prime`](crate::bessel_j_prime): derivative $J_v\'(z)$
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "jv", alias = "cyl_bessel_j")]
#[must_use]
#[inline]
pub fn bessel_j<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_j(v)
}

/// Bessel function of the first kind of order 0, $J_0(x)$
///
/// Corresponds to [`scipy.special.j0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j0.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$ of real order and complex argument
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "j0", alias = "cyl_bessel_j0")]
#[must_use]
#[inline]
pub fn bessel_j0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_j0(x) }
}

/// Bessel function of the first kind of order 1, $J_1(x)$
///
/// Corresponds to [`scipy.special.j1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j1.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$ of real order and complex argument
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "j1", alias = "cyl_bessel_j1")]
#[must_use]
#[inline]
pub fn bessel_j1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_j1(x) }
}

/// Exponentially scaled Bessel function of the first kind
///
/// Corresponds to [`scipy.special.jve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jve.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$
#[doc(alias = "jve", alias = "cyl_bessel_je")]
#[must_use]
#[inline]
pub fn bessel_je<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_je(v)
}

// Bessel Y (Neumann)

/// Bessel function of the second kind of order 0, $Y_0(x)$
///
/// Corresponds to [`scipy.special.y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y0.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$ of real order and complex argument
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "y0", alias = "cyl_bessel_y0", alias = "cyl_neumann_0")]
#[must_use]
#[inline]
pub fn bessel_y0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_y0(x) }
}

/// Bessel function of the second kind of order 1, $Y_1(x)$
///
/// Corresponds to [`scipy.special.y1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$ of real order and complex argument
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "y1", alias = "cyl_bessel_y1", alias = "cyl_neumann_1")]
#[must_use]
#[inline]
pub fn bessel_y1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_y1(x) }
}

/// Bessel function of the second kind, $Y_v(z)$
///
/// Corresponds to [`scipy.special.yv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yv.html
///
/// # See also
/// - [`bessel_y0`]: faster version of $Y_0$ for real argument
/// - [`bessel_y1`]: faster version of $Y_1$ for real argument
/// - [`bessel_ye`]: $Y_v$ without the leading exponential factor
/// - [`bessel_y_prime`](crate::bessel_y_prime): derivative $Y_v\'(z)$
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "yv", alias = "cyl_neumann", alias = "cyl_bessel_y")]
#[must_use]
#[inline]
pub fn bessel_y<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_y(v)
}

/// Exponentially scaled Bessel function of the second kind
///
/// Corresponds to [`scipy.special.yve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yve.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$
#[doc(alias = "yve", alias = "cyl_neumann_e", alias = "cyl_bessel_ye")]
#[must_use]
#[inline]
pub fn bessel_ye<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ye(v)
}

// Bessel I

/// Modified Bessel function of the first kind of order 0, $I_0(x)$
///
/// Corresponds to [`scipy.special.i0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$ of real order and complex argument
/// - [`bessel_i0e`]: exponentially scaled $I_0$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "i0", alias = "cyl_bessel_i0")]
#[must_use]
#[inline]
pub fn bessel_i0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i0(x) }
}

/// Exponentially scaled modified Bessel function of the first kind of order 0, $e^{-|x|} I_0(x)$
///
/// Corresponds to [`scipy.special.i0e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0e.html
///
/// # See also
/// - [`bessel_i0`]: modified Bessel function $I_0(x)$
/// - [`bessel_ie`]: exponentially scaled $I_v(z)$
#[doc(alias = "i0e", alias = "cyl_bessel_i0e")]
#[must_use]
#[inline]
pub fn bessel_i0e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i0e(x) }
}

/// Modified Bessel function of the first kind of order 1, $I_1(x)$
///
/// Corresponds to [`scipy.special.i1`][i1].
///
/// [i1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$ of real order and complex argument
/// - [`bessel_i1e`]: exponentially scaled $I_1$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "i1", alias = "cyl_bessel_i1")]
#[must_use]
#[inline]
pub fn bessel_i1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i1(x) }
}

/// Exponentially scaled modified Bessel function of the first kind of order 1, $e^{-|x|} I_1(x)$
///
/// Corresponds to [`scipy.special.i1e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1e.html
///
/// # See also
/// - [`bessel_i1`]: modified Bessel function $I_1(x)$
/// - [`bessel_ie`]: exponentially scaled $I_v(z)$
#[doc(alias = "i1e", alias = "cyl_bessel_i1e")]
#[must_use]
#[inline]
pub fn bessel_i1e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_i1e(x) }
}

/// Modified Bessel function of the first kind, $I_v(z)$
///
/// Corresponds to [`scipy.special.iv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iv.html
///
/// # See also
/// - [`bessel_i0`]: faster version of $I_0$ for real argument
/// - [`bessel_i1`]: faster version of $I_1$ for real argument
/// - [`bessel_ie`]: $I_v$ without the leading exponential factor
/// - [`bessel_i_prime`](crate::bessel_i_prime): derivative $I_v\'(z)$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "iv", alias = "cyl_bessel_i")]
#[must_use]
#[inline]
pub fn bessel_i<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_i(v)
}

/// Exponentially scaled modified Bessel function of the first kind
///
/// Corresponds to [`scipy.special.ive`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ive.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$
#[doc(alias = "ive", alias = "cyl_bessel_ie")]
#[must_use]
#[inline]
pub fn bessel_ie<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ie(v)
}

// Bessel K

/// Modified Bessel function of the second kind of order 0, $K_0(x)$
///
/// Corresponds to [`scipy.special.k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$ of real order and complex argument
/// - [`bessel_k0e`]: exponentially scaled $K_0$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "k0", alias = "cyl_bessel_k0")]
#[must_use]
#[inline]
pub fn bessel_k0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k0(x) }
}

/// Exponentially scaled modified Bessel function of the second kind of order 0, $e^x K_0(x)$
///
/// Corresponds to [`scipy.special.k0e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0e.html
///
/// # See also
/// - [`bessel_k0`]: modified Bessel function $K_0(x)$
/// - [`bessel_ke`]: exponentially scaled $K_v(z)$
#[doc(alias = "k0e", alias = "cyl_bessel_k0e")]
#[must_use]
#[inline]
pub fn bessel_k0e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k0e(x) }
}

/// Modified Bessel function of the second kind of order 1, $K_1(x)$
///
/// Corresponds to [`scipy.special.k1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$ of real order and complex argument
/// - [`bessel_k1e`]: exponentially scaled $K_1$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "k1", alias = "cyl_bessel_k1")]
#[must_use]
#[inline]
pub fn bessel_k1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k1(x) }
}

/// Exponentially scaled modified Bessel function of the second kind of order 1, $e^x K_1(x)$
///
/// Corresponds to [`scipy.special.k1e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1e.html
///
/// # See also
/// - [`bessel_k1`]: modified Bessel function $K_1(x)$
/// - [`bessel_ke`]: exponentially scaled $K_v(z)$
#[doc(alias = "k1e", alias = "cyl_bessel_k1e")]
#[must_use]
#[inline]
pub fn bessel_k1e(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cyl_bessel_k1e(x) }
}

/// Modified Bessel function of the second kind, $K_v(z)$
///
/// Corresponds to [`scipy.special.kv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kv.html
///
/// # See also
/// - [`bessel_k0`]: faster version of $K_0$ for real argument
/// - [`bessel_k1`]: faster version of $K_1$ for real argument
/// - [`bessel_ke`]: $K_v$ without the leading exponential factor
/// - [`bessel_k_prime`](crate::bessel_k_prime): derivative $K_v\'(z)$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "kv", alias = "cyl_bessel_k")]
#[must_use]
#[inline]
pub fn bessel_k<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_k(v)
}

/// Exponentially scaled modified Bessel function of the second kind
///
/// Corresponds to [`scipy.special.kve`][kve].
///
/// [kve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kve.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$
#[doc(alias = "kve", alias = "cyl_bessel_ke")]
#[must_use]
#[inline]
pub fn bessel_ke<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ke(v)
}

// Zeros

/// Compute $N$ zeros of Bessel functions $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$
///
/// Corresponds to [`scipy.special.jnyn_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnyn_zeros.html
///
/// # Panics
/// - Panics if `N` is greater than 1200
#[doc(alias = "jnyn_zeros")]
#[must_use]
#[inline]
pub fn bessel_zeros<const N: usize>(v: u32) -> [[f64; N]; 4] {
    assert!(N <= 1200, "N must be at most 1200");

    let (mut j, mut jp) = ([f64::NAN; N], [f64::NAN; N]);
    let (mut y, mut yp) = ([f64::NAN; N], [f64::NAN; N]);
    unsafe {
        crate::ffi::xsf::jyzo(
            v.try_into().unwrap(),
            N.try_into().unwrap(),
            j.as_mut_ptr(),
            jp.as_mut_ptr(),
            y.as_mut_ptr(),
            yp.as_mut_ptr(),
        );
    }
    [j, jp, y, yp]
}

// Hankel 1

/// Hankel function of the first kind, $H_v^{(1)}(z)$
///
/// Corresponds to [`scipy.special.hankel1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1.html
///
/// # See also
/// - [`hankel_1e`]: $H_v^{(1)}$ without the leading exponential factor
/// - [`hankel_1_prime`](crate::hankel_1_prime): $n$th derivative, ${d^n\over dz^n}H_v^{(1)}(z)$
/// - [`hankel_2`]: Hankel function of the second kind $H_v^{(2)}(z)$
#[doc(alias = "h1v", alias = "hankel1", alias = "cyl_hankel_1")]
#[must_use]
#[inline]
pub fn hankel_1<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_1(v)
}

/// Exponentially scaled Hankel function of the first kind
///
/// Corresponds to [`scipy.special.hankel1e`][h1ve].
///
/// [h1ve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1e.html
///
/// # See also
/// - [`hankel_1`]: Hankel function $H_v^{(1)}(z)$
#[doc(alias = "h1ve", alias = "hankel1e", alias = "cyl_hankel_1e")]
#[must_use]
#[inline]
pub fn hankel_1e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_1e(v)
}

// Hankel 2

/// Hankel function of the second kind, $H_v^{(2)}(z)$
///
/// Corresponds to [`scipy.special.hankel2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2.html
///
/// # See also
/// - [`hankel_2e`]: $H_v^{(2)}$ without the leading exponential factor
/// - [`hankel_2_prime`](crate::hankel_2_prime): $n$th derivative, ${d^n\over dz^n}H_v^{(2)}(z)$
/// - [`hankel_1`]: Hankel function of the first kind $H_v^{(1)}(z)$
#[doc(alias = "h2v", alias = "hankel2", alias = "cyl_hankel_2")]
#[must_use]
#[inline]
pub fn hankel_2<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_2(v)
}

/// Exponentially scaled Hankel function of the second kind
///
/// Corresponds to [`scipy.special.hankel2e`][h2ve].
///
/// [h2ve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2e.html
///
/// # See also
/// - [`hankel_2`]: Hankel function $H_v^{(2)}(z)$
#[doc(alias = "h2ve", alias = "hankel2e", alias = "cyl_hankel_2e")]
#[must_use]
#[inline]
pub fn hankel_2e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_2e(v)
}

// Integrals

/// Weighted integral of the Bessel function of the first kind,
/// $\int_0^1 x^\lambda \mathop{J}_v(2ax) \dd x$
///
/// Corresponds to [`scipy.special.besselpoly`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.besselpoly.html
#[must_use]
#[inline]
pub fn besselpoly(a: f64, lambda: f64, v: f64) -> f64 {
    unsafe { crate::ffi::xsf::besselpoly(a, lambda, v) }
}

/// Integrals of Bessel functions of order 0
///
/// Corresponds to [`scipy.special.itj0y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itj0y0.html
///
/// # Returns
/// - $\int_0^x J_0(t) \dd t$
/// - $\int_0^x Y_0(t) \dd t$
///
/// # See also
/// - [`it2j0y0`]
/// - [`bessel_j0`]
/// - [`bessel_y0`]
#[doc(alias = "itj0y0")]
#[must_use]
#[inline]
pub fn it1j0y0(x: f64) -> (f64, f64) {
    let (mut j0int, mut y0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it1j0y0(x, &raw mut j0int, &raw mut y0int);
    }
    (j0int, y0int)
}

/// Integrals related to Bessel functions of the first kind of order 0
///
/// Corresponds to [`scipy.special.it2j0y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2j0y0.html
///
/// # Returns
/// - $\int_0^x {1 - J_0(t) \over t} \dd t$
/// - $\int_x^\infty {Y_0(t) \over t} \dd t$
///
/// # See also
/// - [`it1j0y0`]
/// - [`bessel_j0`]
/// - [`bessel_y0`]
#[must_use]
#[inline]
pub fn it2j0y0(x: f64) -> (f64, f64) {
    let (mut j0int, mut y0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it2j0y0(x, &raw mut j0int, &raw mut y0int);
    }
    (j0int, y0int)
}

/// Integrals of modified Bessel functions of order 0
///
/// Corresponds to [`scipy.special.iti0k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iti0k0.html
///
/// # Returns
/// - $\int_0^x I_0(t) \dd t$
/// - $\int_0^x K_0(t) \dd t$
///
/// # See also
/// - [`it2i0k0`]
/// - [`bessel_i0`]
/// - [`bessel_k0`]
#[doc(alias = "iti0k0")]
#[must_use]
#[inline]
pub fn it1i0k0(x: f64) -> (f64, f64) {
    let (mut i0int, mut k0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it1i0k0(x, &raw mut i0int, &raw mut k0int);
    }
    (i0int, k0int)
}

/// Integrals related to modified Bessel functions of order 0.
///
/// Corresponds to [`scipy.special.it2i0k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2i0k0.html
///
/// # Returns
/// - $\int_0^x {I_0(t) - 1 \over t} \dd t$
/// - $\int_x^\infty {K_0(t) \over t} \dd t$
///
/// # See also
/// - [`it1i0k0`]
/// - [`bessel_i0`]
/// - [`bessel_k0`]
#[must_use]
#[inline]
pub fn it2i0k0(x: f64) -> (f64, f64) {
    let (mut i0int, mut k0int) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::it2i0k0(x, &raw mut i0int, &raw mut k0int);
    }
    (i0int, k0int)
}

// Riccati-Bessel

/// Compute Riccati-Bessel function of the first kind and derivatives for the first $N$ orders
///
/// The Riccati-Bessel function of the first kind is defined as $x \mathop{j}_n(x)$, where $j_n$ is
/// the spherical Bessel function of the first kind of order $n$.
///
/// This function computes the value and first derivative of the function for all orders up to $N$.
///
/// Corresponds to [`scipy.special.riccati_jn`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_jn.html
///
/// # Returns
/// - $\[j_0(x), \ldots, j_{N-1}(x)\]$
/// - $\[j_0\'(x), \ldots, j_{N-1}\'(x)\]$
///
/// # See also
/// - [`riccati_y`]
/// - [`bessel_j`]
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
#[doc(alias = "rctj", alias = "riccati_jn")]
#[must_use]
#[inline]
pub fn riccati_j<const N: usize>(x: f64) -> ([f64; N], [f64; N]) {
    let (mut rj, mut dj) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 {
        return (rj, dj);
    }

    let _ = unsafe { crate::ffi::xsf::rctj(N, x, rj.as_mut_ptr(), dj.as_mut_ptr()) };
    (rj, dj)
}

/// Compute Riccati-Bessel function of the second kind and derivatives for the first $N$ orders
///
/// The Riccati-Bessel function of the second kind is defined here as $x \mathop{y}_n(x)$, where
/// $y_n$ is the spherical Bessel function of the second kind of order $n$.
/// *Note that this is in contrast to a common convention that includes a minus sign in the
/// definition.*
///
/// This function computes the value and first derivative of the function for all orders up to $N$.
///
/// Corresponds to [`scipy.special.riccati_yn`][scipy].
///
/// # Returns
/// - $\[y_0(x), \ldots, y_{N-1}(x)\]$
/// - $\[y_0\'(x), \ldots, y_{N-1}\'(x)\]$
///
/// # See also
/// - [`riccati_j`]
/// - [`bessel_y`](crate::bessel_y)
/// - [`sph_bessel_y`](crate::sph_bessel_y)
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_yn.html
#[doc(alias = "rcty", alias = "riccati_yn")]
#[must_use]
#[inline]
pub fn riccati_y<const N: usize>(x: f64) -> ([f64; N], [f64; N]) {
    let (mut ry, mut dy) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 {
        return (ry, dy);
    }

    let _ = unsafe { crate::ffi::xsf::rcty(N, x, ry.as_mut_ptr(), dy.as_mut_ptr()) };
    (ry, dy)
}

/// Evaluate the Jahnke-Emden Lambda function $\Lambda_v(x)$ and its derivatives
///
/// The Jahnke-Emden Lambda function is defined as
///
/// $$\Lambda_v(x) = \Gamma(v+1) \left({2 \over x}\right)^v J_v(x)$$
///
/// where $\Gamma$ is the [Gamma function](crate::gamma) and $J_v$ is the
/// [Bessel function of the first kind](crate::bessel_j).
///
/// Corresponds to [`scipy.special.lmbda`][scipy] in SciPy, and calls the FFI functions
/// `xsf::specfun::lamn` and `xsf::specfun::lamv`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lmbda.html
///
/// # Arguments
/// - `v`: Order of the Lambda function as a non-negative real number
/// - `x`: Value at which to evaluate the function and derivatives
///
/// # Returns
/// - `vl`: Values of $\Lambda_{v_i}(x)$ for
///   $v_i = v-\lfloor v\rfloor, v-\lfloor v\rfloor+1, \ldots, v$
/// - `dl`: Derivatives $\Lambda_{v_i}\'(x)$ for
///   $v_i = v-\lfloor v\rfloor, v-\lfloor v\rfloor+1, \ldots, v$
///
/// Both vectors have length $\lfloor v\rfloor+1$, i.e. `(v as usize) + 1`.
///
/// # See also
/// - [`bessel_j`]
/// - [`gamma`](crate::gamma)
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
/// # Panics
/// - Panics if `v` is negative
/// - Panics if `v` is too large to convert to `c_int`
#[doc(alias = "lmbda")]
#[must_use]
pub fn jahnke_emden_lambda<V: Into<f64>>(v: V, x: f64) -> (Vec<f64>, Vec<f64>) {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L2008-L2055
    let v: f64 = v.into();
    assert!(v >= 0.0, "v must be non-negative");

    let n = v.to_usize().unwrap();
    let v0 = v - n.to_f64().unwrap();
    let n1 = if n < 1 { 1 } else { n };
    let v1 = v0 + n1.to_f64().unwrap();

    let size = v1.to_usize().unwrap() + 1;
    let (mut vl, mut dl) = (vec![f64::NAN; size], vec![f64::NAN; size]);
    unsafe {
        #[allow(clippy::float_cmp)]
        if v == n.to_f64().unwrap() {
            crate::ffi::xsf::lamn(
                v1.to_usize().unwrap().try_into().unwrap(),
                x,
                vl.as_mut_ptr(),
                dl.as_mut_ptr(),
            );
        } else {
            crate::ffi::xsf::lamv(v1, x, vl.as_mut_ptr(), dl.as_mut_ptr());
        }
    };

    vl.resize(n + 1, f64::NAN);
    dl.resize(n + 1, f64::NAN);
    (vl, dl)
}

// Tests

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;
    use num_complex::c64;

    // Bessel J

    #[test]
    fn test_bessel_j0() {
        xsref::test("cyl_bessel_j0", "d-d", |x| crate::bessel_j0(x[0]));
    }

    #[test]
    fn test_bessel_j1() {
        xsref::test("cyl_bessel_j1", "d-d", |x| crate::bessel_j1(x[0]));
    }

    #[test]
    fn test_bessel_j_f64() {
        xsref::test("cyl_bessel_j", "d_d-d", |x| crate::bessel_j(x[0], x[1]));
    }

    #[test]
    fn test_bessel_j_c64() {
        xsref::test("cyl_bessel_j", "d_cd-cd", |x| {
            crate::bessel_j(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_je_f64() {
        xsref::test("cyl_bessel_je", "d_d-d", |x| crate::bessel_je(x[0], x[1]));
    }

    #[test]
    fn test_bessel_je_c64() {
        xsref::test("cyl_bessel_je", "d_cd-cd", |x| {
            crate::bessel_je(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel Y

    #[test]
    fn test_bessel_y0() {
        xsref::test("cyl_bessel_y0", "d-d", |x| crate::bessel_y0(x[0]));
    }

    #[test]
    fn test_bessel_y1() {
        xsref::test("cyl_bessel_y1", "d-d", |x| crate::bessel_y1(x[0]));
    }

    #[test]
    fn test_bessel_y_f64() {
        xsref::test("cyl_bessel_y", "d_d-d", |x| crate::bessel_y(x[0], x[1]));
    }

    #[test]
    fn test_bessel_y_c64() {
        xsref::test("cyl_bessel_y", "d_cd-cd", |x| {
            crate::bessel_y(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ye_f64() {
        xsref::test("cyl_bessel_ye", "d_d-d", |x| crate::bessel_ye(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ye_c64() {
        xsref::test("cyl_bessel_ye", "d_cd-cd", |x| {
            crate::bessel_ye(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel I

    #[test]
    fn test_bessel_i0() {
        xsref::test("cyl_bessel_i0", "d-d", |x| crate::bessel_i0(x[0]));
    }

    #[test]
    fn test_bessel_i0e() {
        xsref::test("cyl_bessel_i0e", "d-d", |x| crate::bessel_i0e(x[0]));
    }

    #[test]
    fn test_bessel_i1() {
        xsref::test("cyl_bessel_i1", "d-d", |x| crate::bessel_i1(x[0]));
    }

    #[test]
    fn test_bessel_i1e() {
        xsref::test("cyl_bessel_i1e", "d-d", |x| crate::bessel_i1e(x[0]));
    }

    #[test]
    fn test_bessel_i_f64() {
        xsref::test("cyl_bessel_i", "d_d-d", |x| crate::bessel_i(x[0], x[1]));
    }

    #[test]
    fn test_bessel_i_c64() {
        xsref::test("cyl_bessel_i", "d_cd-cd", |x| {
            crate::bessel_i(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ie_f64() {
        xsref::test("cyl_bessel_ie", "d_d-d", |x| crate::bessel_ie(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ie_c64() {
        xsref::test("cyl_bessel_ie", "d_cd-cd", |x| {
            crate::bessel_ie(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel K

    #[test]
    fn test_bessel_k0() {
        xsref::test("cyl_bessel_k0", "d-d", |x| crate::bessel_k0(x[0]));
    }

    #[test]
    fn test_bessel_k0e() {
        xsref::test("cyl_bessel_k0e", "d-d", |x| crate::bessel_k0e(x[0]));
    }

    #[test]
    fn test_bessel_k1() {
        xsref::test("cyl_bessel_k1", "d-d", |x| crate::bessel_k1(x[0]));
    }

    #[test]
    fn test_bessel_k1e() {
        xsref::test("cyl_bessel_k1e", "d-d", |x| crate::bessel_k1e(x[0]));
    }

    #[test]
    fn test_bessel_k_f64() {
        xsref::test("cyl_bessel_k", "d_d-d", |x| crate::bessel_k(x[0], x[1]));
    }

    #[test]
    fn test_bessel_k_c64() {
        xsref::test("cyl_bessel_k", "d_cd-cd", |x| {
            crate::bessel_k(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ke_f64() {
        xsref::test("cyl_bessel_ke", "d_d-d", |x| crate::bessel_ke(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ke_c64() {
        xsref::test("cyl_bessel_ke", "d_cd-cd", |x| {
            crate::bessel_ke(x[0], c64(x[1], x[2]))
        });
    }

    // Hankel

    #[test]
    fn test_hankel_1_c64() {
        xsref::test("cyl_hankel_1", "d_cd-cd", |x| {
            crate::hankel_1(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_1e_c64() {
        xsref::test("cyl_hankel_1e", "d_cd-cd", |x| {
            crate::hankel_1e(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_2_c64() {
        xsref::test("cyl_hankel_2", "d_cd-cd", |x| {
            crate::hankel_2(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_2e_c64() {
        xsref::test("cyl_hankel_2e", "d_cd-cd", |x| {
            crate::hankel_2e(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel zeros

    /// Based on `scipy.special.tests.test_basic.TestBessel.test_jnyn_zeros`
    #[test]
    fn test_bessel_zeros() {
        let jnz = crate::bessel_zeros::<5>(1);
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
        xsref::test("it1j0y0", "d-d_d", |x| crate::it1j0y0(x[0]));
    }

    #[test]
    fn test_it2j0y0() {
        // This xsref table contains incorrect values for the 4th integral,
        // see https://github.com/scipy/xsref/issues/10.
        // So we instead use values that were verified using a Wolfram notebook.

        // xsref::test("it2j0y0", "d-d_d", |x| crate::it2j0y0(x[0]));

        let xs = [-1.0, 10.0, -10.0, 1.0, 0.2];
        let expect_0 = [
            0.121_165_246_995_068_71,
            2.177_866_420_093_336,
            2.177_866_420_093_336,
            0.121_165_246_995_068_71,
            0.004_993_754_627_460_185,
        ];
        let expect_1 = [
            f64::NAN,
            -0.022_987_933_564_659_65,
            f64::NAN,
            0.395_272_901_699_293_2,
            -0.434_230_670_112_316_34,
        ];
        np_assert_allclose!(xs.map(|x| crate::it2j0y0(x).0), expect_0);
        np_assert_allclose!(xs.map(|x| crate::it2j0y0(x).1), expect_1);
    }

    #[test]
    fn test_it1i0k0() {
        xsref::test("it1i0k0", "d-d_d", |x| crate::it1i0k0(x[0]));
    }

    #[test]
    fn test_it2i0k0() {
        xsref::test("it2i0k0", "d-d_d", |x| crate::it2i0k0(x[0]));
    }

    #[test]
    fn test_besselpoly() {
        xsref::test("besselpoly", "d_d_d-d", |x| {
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
            let j = crate::sph_bessel_j(n.try_into().unwrap(), x);
            // jp = special.spherical_jn(n, x, derivative=True)
            let jp = crate::sph_bessel_j_prime(n.try_into().unwrap(), x);
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
            let j = crate::sph_bessel_y(n.try_into().unwrap(), x);
            // jp = special.spherical_jn(n, x, derivative=True)
            let jp = crate::sph_bessel_y_prime(n.try_into().unwrap(), x);
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
