//! Derivatives of Bessel functions
use crate::xsf::bessel::BesselArg;

/// Translated from <https://github.com/scipy/scipy/blob/9531cc5/scipy/special/_basic.py#L803-L814>
#[inline]
fn bessel_diff_formula<T, L>(v: f64, n: u32, bessel_fn: L, phase: f64) -> T
where
    T: BesselArg,
    L: Fn(f64) -> T,
{
    if n == 0 {
        return bessel_fn(v);
    }

    // p = 1.0
    let mut p = 1.0;
    // s = L(v-n, z)
    let mut s = bessel_fn(v - f64::from(n));
    // for i in range(1, n+1):
    for i in 1..=n {
        let (i, n) = (f64::from(i), f64::from(n));
        // p = phase * (p * (n-i+1)) / i   # = choose(k, i)
        p = phase * (p * (n - i + 1.0)) / i;
        // s += p*L(v-n + i*2, z)
        s = s + bessel_fn(v - (n - i * 2.0)) * p.into();
    }
    // return s / (2.**n)
    // (we can't use f64::powi or f64::powf because of no_std)
    s * (1.0 / f64::from(1 << n)).into()
}

/// Compute the $n$th derivative of [`bessel_j(v, z)`](crate::bessel_j) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.jvp`][scipy-jvp].
///
/// [scipy-jvp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jvp.html
///
/// # See also
/// - [`bessel_j`](crate::bessel_j): Bessel function $J_v(z)$
/// - [`sph_bessel_j_prime`](crate::sph_bessel_j_prime): spherical Bessel derivative $j_n\'(z)$
#[doc(alias = "jvp")]
#[inline]
pub fn bessel_j_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    bessel_diff_formula(v, n, |v| z.bessel_j(v), -1.0)
}

/// Compute the $n$th derivative of [`bessel_y(v, z)`](crate::bessel_y) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.yvp`][scipy-yvp].
///
/// [scipy-yvp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yvp.html
///
/// # See also
/// - [`bessel_y`](crate::bessel_y): Bessel function $Y_v(z)$
/// - [`sph_bessel_y_prime`](crate::sph_bessel_y_prime): spherical Bessel derivative $y_n\'(z)$
#[doc(alias = "yvp")]
#[inline]
pub fn bessel_y_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    bessel_diff_formula(v, n, |v| z.bessel_y(v), -1.0)
}

/// Compute the $n$th derivative of [`bessel_i(v, z)`](crate::bessel_i) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.ivp`][scipy-ivp].
///
/// [scipy-ivp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ivp.html
///
/// # See also
/// - [`bessel_i`](crate::bessel_i): modified Bessel function $I_v(z)$
/// - [`sph_bessel_i_prime`](crate::sph_bessel_i_prime): spherical modified Bessel derivative $i_n\'(z)$
#[doc(alias = "ivp")]
#[inline]
pub fn bessel_i_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    bessel_diff_formula(v, n, |v| z.bessel_i(v), 1.0)
}

/// Compute the $n$th derivative of [`bessel_k(v, z)`](crate::bessel_k) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.kvp`][scipy-kvp].
///
/// [scipy-kvp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kvp.html
///
/// # See also
/// - [`bessel_k`](crate::bessel_k): modified Bessel function $K_v(z)$
/// - [`sph_bessel_k_prime`](crate::sph_bessel_k_prime): spherical modified Bessel derivative $k_n\'(z)$
#[doc(alias = "kvp")]
#[inline]
pub fn bessel_k_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
    bessel_diff_formula(v, n, |v| z.bessel_k(v), 1.0) * sign.into()
}

/// Compute the $n$th derivative of [`hankel_1(v, z)`](crate::hankel_1) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.h1vp`][scipy-h1vp].
///
/// [scipy-h1vp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.h1vp.html
///
/// # See also
/// - [`hankel_1`](crate::hankel_1): Hankel function $H_v^{(1)}(z)$
#[doc(alias = "h1vp")]
#[inline]
pub fn hankel_1_prime<T: BesselArg>(v: f64, z: T, n: u32) -> num_complex::Complex<f64> {
    bessel_diff_formula(v, n, |v| z.hankel_1(v), -1.0)
}

/// Compute the $n$th derivative of [`hankel_2(v, z)`](crate::hankel_2) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.h2vp`][scipy-h2vp].
///
/// [scipy-h2vp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.h2vp.html
///
/// # See also
/// - [`hankel_2`](crate::hankel_2): Hankel function $H_v^{(2)}(z)$
#[doc(alias = "h2vp")]
#[inline]
pub fn hankel_2_prime<T: BesselArg>(v: f64, z: T, n: u32) -> num_complex::Complex<f64> {
    bessel_diff_formula(v, n, |v| z.hankel_2(v), -1.0)
}
