//! Derivatives of Bessel functions
use num_complex::Complex;

use crate::xsf::bessel::BesselArg;

/// Translated from https://github.com/scipy/scipy/blob/9531cc5/scipy/special/_basic.py#L803-L814
#[inline(always)]
fn _bessel_diff_formula<T, L>(v: f64, n: u32, bessel_fn: L, phase: f64) -> T
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
    let mut s = bessel_fn(v - n as f64);
    // for i in range(1, n+1):
    for i in 1..n + 1 {
        let (i, n) = (i as f64, n as f64);
        // p = phase * (p * (n-i+1)) / i   # = choose(k, i)
        p = phase * (p * (n - i + 1.0)) / i;
        // s += p*L(v-n + i*2, z)
        s = s + bessel_fn(v - (n - i * 2.0)) * p.into();
    }
    // return s / (2.**n)
    s * (0.5f64).powi(n as i32).into()
}

/// Compute the *n*<sup>th</sup> derivative of [`bessel_j(v, z)`](crate::bessel_j) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.jvp`][scipy-jvp].
///
/// [scipy-jvp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jvp.html
#[doc(alias = "jvp")]
pub fn bessel_j_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    _bessel_diff_formula(v, n, |v| z.bessel_j(v), -1.0)
}

/// Compute the *n*<sup>th</sup> derivative of [`bessel_y(v, z)`](crate::bessel_y) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.yvp`][scipy-yvp].
///
/// [scipy-yvp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yvp.html
#[doc(alias = "yvp")]
pub fn bessel_y_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    _bessel_diff_formula(v, n, |v| z.bessel_y(v), -1.0)
}

/// Compute the *n*<sup>th</sup> derivative of [`bessel_i(v, z)`](crate::bessel_i) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.ivp`][scipy-ivp].
///
/// [scipy-ivp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ivp.html
#[doc(alias = "ivp")]
pub fn bessel_i_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    _bessel_diff_formula(v, n, |v| z.bessel_i(v), 1.0)
}

/// Compute the *n*<sup>th</sup> derivative of [`bessel_k(v, z)`](crate::bessel_k) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.kvp`][scipy-kvp].
///
/// [scipy-kvp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kvp.html
#[doc(alias = "kvp")]
pub fn bessel_k_prime<T: BesselArg>(v: f64, z: T, n: u32) -> T {
    _bessel_diff_formula(v, n, |v| z.bessel_k(v), 1.0) * (-1.0f64).powi(n as i32).into()
}

/// Compute the *n*<sup>th</sup> derivative of [`hankel_1(v, z)`](crate::hankel_1) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.h1vp`][scipy-h1vp].
///
/// [scipy-h1vp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.h1vp.html
#[doc(alias = "h1vp")]
pub fn hankel_1_prime<T: BesselArg>(v: f64, z: T, n: u32) -> Complex<f64> {
    _bessel_diff_formula::<Complex<f64>, _>(v, n, |v| z.hankel_1(v), -1.0)
}

/// Compute the *n*<sup>th</sup> derivative of [`hankel_2(v, z)`](crate::hankel_2) w.r.t. `z`
///
/// Pure rust implementation of [`scipy.special.h2vp`][scipy-h2vp].
///
/// [scipy-h2vp]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.h2vp.html
#[doc(alias = "h2vp")]
pub fn hankel_2_prime<T: BesselArg>(v: f64, z: T, n: u32) -> Complex<f64> {
    _bessel_diff_formula::<Complex<f64>, _>(v, n, |v| z.hankel_2(v), -1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::np_assert_allclose;

    const ATOL: f64 = 1.5e-10;

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_jvp`
    #[test]
    fn test_bessel_j_prime() {
        // jvprim = special.jvp(2,2)
        let jvprim = bessel_j_prime(2.0, 2.0, 1);
        // jv0 = (special.jv(1,2)-special.jv(3,2))/2
        let jv0 = (crate::bessel_j(1.0, 2.0) - crate::bessel_j(3.0, 2.0)) * 0.5;
        // assert_allclose(jvprim, jv0, atol=1.5e-10, rtol=0)
        np_assert_allclose(&[jvprim], &[jv0], 0.0, ATOL);
    }

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_yvp`
    #[test]
    fn test_bessel_y_prime() {
        // yvpr = (special.yv(1,.2) - special.yv(3,.2))/2.0
        let yvpr = (crate::bessel_y(1.0, 0.2) - crate::bessel_y(3.0, 0.2)) * 0.5;
        // yvp1 = special.yvp(2,.2)
        let yvp1 = bessel_y_prime(2.0, 0.2, 1);
        // assert_allclose(yvp1, yvpr, atol=1.5e-10, rtol=0)
        np_assert_allclose(&[yvp1], &[yvpr], 0.0, ATOL);
    }

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_ivp0`
    #[test]
    fn test_bessel_i_prime_0() {
        // assert_allclose(special.iv(1, 2), special.ivp(0, 2), atol=1.5e-10, rtol=0)
        np_assert_allclose(
            &[crate::bessel_i(1.0, 2.0)],
            &[bessel_i_prime(0.0, 2.0, 1)],
            0.0,
            ATOL,
        );
    }

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_ivp`
    #[test]
    fn test_bessel_i_prime() {
        // y = (special.iv(0,2) + special.iv(2,2))/2
        let y = (crate::bessel_i(0.0, 2.0) + crate::bessel_i(2.0, 2.0)) * 0.5;
        // x = special.ivp(1,2)
        let x = bessel_i_prime(1.0, 2.0, 1);
        // assert_allclose(x, y, atol=1.5e-10, rtol=0)
        np_assert_allclose(&[x], &[y], 0.0, ATOL);
    }

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_kvp_v0n1`
    #[test]
    fn test_bessel_k_prime_v0n1() {
        // z = 2.2
        const Z: f64 = 2.2;
        // assert_allclose(-special.kv(1, z), special.kvp(0, z, n=1), atol=1.5e-10, rtol=0)
        np_assert_allclose(
            &[-crate::bessel_k(1.0, Z)],
            &[bessel_k_prime(0.0, Z, 1)],
            0.0,
            ATOL,
        );
    }

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_kvp_n1`
    #[test]
    fn test_bessel_k_prime_n1() {
        // v = 3.
        const V: f64 = 3.0;
        // z = 2.2
        const Z: f64 = 2.2;
        // xc = -special.kv(v+1,z) + v/z*special.kv(v,z)
        let xc = -crate::bessel_k(V + 1.0, Z) + V / Z * crate::bessel_k(V, Z);
        // x = special.kvp(v,z, n=1)
        let x = bessel_k_prime(V, Z, 1);
        // assert_allclose(xc, x, atol=1.5e-10, rtol=0)
        np_assert_allclose(&[x], &[xc], 0.0, ATOL);
    }

    /// Translated from `scipy.special.tests.test_basic.TestBessel.test_kvp_n2`
    #[test]
    fn test_bessel_k_prime_n2() {
        // v = 3.
        const V: f64 = 3.0;
        // z = 2.2
        const Z: f64 = 2.2;
        // xc = (z**2 + v**2 - v) / z**2 * special.kv(v, z) + special.kv(v+1, z) / z
        const Z2: f64 = Z * Z;
        let xc = (Z2 + V * V - V) / Z2 * crate::bessel_k(V, Z) + crate::bessel_k(V + 1.0, Z) / Z;
        // x = special.kvp(v, z, n=2)
        let x = bessel_k_prime(V, Z, 2);
        // assert_allclose(xc, x, atol=1.5e-10, rtol=0)
        np_assert_allclose(&[x], &[xc], 0.0, ATOL);
    }
}
