use core::ffi::c_long;

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
    #[inline]
    fn sph_bessel_j(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_j(n, self) }
    }
    #[inline]
    fn sph_bessel_y(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_y(n, self) }
    }
    #[inline]
    fn sph_bessel_i(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_i(n, self) }
    }
    #[inline]
    fn sph_bessel_k(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_k(n, self) }
    }
    #[inline]
    fn sph_bessel_j_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_j_jac(n, self) }
    }
    #[inline]
    fn sph_bessel_y_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_y_jac(n, self) }
    }
    #[inline]
    fn sph_bessel_i_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_i_jac(n, self) }
    }
    #[inline]
    fn sph_bessel_k_prime(self, n: c_long) -> f64 {
        unsafe { crate::ffi::xsf::sph_bessel_k_jac(n, self) }
    }
}

impl SphBesselArg for num_complex::Complex<f64> {
    #[inline]
    fn sph_bessel_j(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_j_1(n, self) }
    }
    #[inline]
    fn sph_bessel_y(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_y_1(n, self) }
    }
    #[inline]
    fn sph_bessel_i(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_i_1(n, self) }
    }
    #[inline]
    fn sph_bessel_k(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_k_1(n, self) }
    }
    #[inline]
    fn sph_bessel_j_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_j_jac_1(n, self) }
    }
    #[inline]
    fn sph_bessel_y_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_y_jac_1(n, self) }
    }
    #[inline]
    fn sph_bessel_i_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_i_jac_1(n, self) }
    }
    #[inline]
    fn sph_bessel_k_prime(self, n: c_long) -> Self {
        unsafe { crate::ffi::xsf::sph_bessel_k_jac_1(n, self) }
    }
}

/// Spherical Bessel function of the first kind, $j_n(z)$
///
/// Corresponds to [`scipy.special.spherical_jn(n, z)`][jn].
///
/// [jn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_jn.html
///
/// # See also
/// - [`sph_bessel_j_prime`]: derivative $j_n\'(z)$
/// - [`bessel_j`](crate::bessel_j): Bessel function $J_v(z)$
#[doc(alias = "spherical_jn")]
pub fn sph_bessel_j<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_j(n as c_long)
}

/// Derivative of [`sph_bessel_j`], $j_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_jn(n, z, derivative=True)`][jn].
///
/// [jn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_jn.html
///
/// # See also
/// - [`sph_bessel_j`]: spherical Bessel function $j_n(z)$
/// - [`bessel_j_prime`](crate::bessel_j_prime): derivative of Bessel function $J_v\'(z)$
#[doc(alias = "spherical_jnp")]
pub fn sph_bessel_j_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_j_prime(n as c_long)
}

/// Spherical Bessel function of the second kind, $y_n(z)$
///
/// Corresponds to [`scipy.special.spherical_yn(n, z)`][yn].
///
/// [yn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_yn.html
///
/// # See also
/// - [`sph_bessel_y_prime`]: derivative $y_n\'(z)$
/// - [`bessel_y`](crate::bessel_y): Bessel function $Y_v(z)$
#[doc(alias = "spherical_yn")]
pub fn sph_bessel_y<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_y(n as c_long)
}

/// Derivative of [`sph_bessel_y`], $y_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_yn(n, z, derivative=True)`][yn].
///
/// [yn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_yn.html
///
/// # See also
/// - [`sph_bessel_y`]: spherical Bessel function $y_n(z)$
/// - [`bessel_y_prime`](crate::bessel_y_prime): derivative of Bessel function $Y_v\'(z)$
#[doc(alias = "spherical_ynp")]
pub fn sph_bessel_y_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_y_prime(n as c_long)
}

/// Modified spherical Bessel function of the first kind, $i_n(z)$
///
/// Corresponds to [`scipy.special.spherical_in(n, z)`][in_].
///
/// [in_]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_in.html
///
/// # See also
/// - [`sph_bessel_i_prime`]: derivative $i_n\'(z)$
/// - [`bessel_i`](crate::bessel_i): modified Bessel function $I_v(z)$
#[doc(alias = "spherical_in")]
pub fn sph_bessel_i<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_i(n as c_long)
}

/// Derivative of [`sph_bessel_i`], $i_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_in(n, z, derivative=True)`][in_].
///
/// [in_]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_in.html
///
/// # See also
/// - [`sph_bessel_i`]: spherical modified Bessel function $i_n(z)$
/// - [`bessel_i_prime`](crate::bessel_i_prime): derivative of modified Bessel function $I_v\'(z)$
#[doc(alias = "spherical_inp")]
pub fn sph_bessel_i_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_i_prime(n as c_long)
}

/// Modified spherical Bessel function of the second kind, $k_n(z)$
///
/// Corresponds to [`scipy.special.spherical_kn(n, z)`][kn].
///
/// [kn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_kn.html
///
/// # See also
/// - [`sph_bessel_k_prime`]: derivative $k_n\'(z)$
/// - [`bessel_k`](crate::bessel_k): modified Bessel function $K_v(z)$
#[doc(alias = "spherical_kn")]
pub fn sph_bessel_k<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_k(n as c_long)
}

/// Derivative of [`sph_bessel_k`], $k_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_kn(n, z, derivative=True)`][kn].
///
/// [kn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_kn.html
///
/// # See also
/// - [`sph_bessel_k`]: spherical modified Bessel function $k_n(z)$
/// - [`bessel_k_prime`](crate::bessel_k_prime): derivative of modified Bessel function $K_v\'(z)$
#[doc(alias = "spherical_knp")]
pub fn sph_bessel_k_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_k_prime(n as c_long)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalJn.test_spherical_jn_exact`
    #[test]
    fn test_sph_bessel_j_exact() {
        // x = np.array([0.12, 1.23, 12.34, 123.45, 1234.5])
        let x: [f64; 5] = [0.12, 1.23, 12.34, 123.45, 1234.5];
        // assert_allclose(spherical_jn(2, x),
        //                 (-1/x + 3/x**3)*sin(x) - 3/x**2*cos(x))
        let actual = x.map(|x| crate::sph_bessel_j(2, x));
        let expected = x.map(|x| ((-1.0 + 3.0 / (x * x)) * x.sin() - 3.0 / x * x.cos()) / x);
        crate::np_assert_allclose!(&actual, &expected);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalYn.test_spherical_yn_exact`
    #[test]
    fn test_sph_bessel_y_exact() {
        // x = np.array([0.12, 1.23, 12.34, 123.45, 1234.5])
        let x: [f64; 5] = [0.12, 1.23, 12.34, 123.45, 1234.5];
        // assert_allclose(spherical_yn(2, x),
        //                 (1/x - 3/x**3)*cos(x) - 3/x**2*sin(x))
        let actual = x.map(|x| crate::sph_bessel_y(2, x));
        let expected = x.map(|x| ((1.0 - 3.0 / (x * x)) * x.cos() - 3.0 / x * x.sin()) / x);
        crate::np_assert_allclose!(&actual, &expected);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalIn.test_spherical_in_exact`
    #[test]
    fn test_sph_bessel_i_exact() {
        // x = np.array([0.12, 1.23, 12.34, 123.45])
        let x: [f64; 4] = [0.12, 1.23, 12.34, 123.45];
        // assert_allclose(spherical_in(2, x),
        //                 (1/x + 3/x**3)*sinh(x) - 3/x**2*cosh(x))
        let actual = x.map(|x| crate::sph_bessel_i(2, x));
        let expected = x.map(|x| ((1.0 + 3.0 / (x * x)) * x.sinh() - 3.0 / x * x.cosh()) / x);
        crate::np_assert_allclose!(&actual, &expected);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalKn.test_spherical_kn_exact`
    #[test]
    fn test_sph_bessel_k_exact() {
        // x = np.array([0.12, 1.23, 12.34, 123.45])
        let x: [f64; 4] = [0.12, 1.23, 12.34, 123.45];
        // assert_allclose(spherical_kn(2, x),
        //                 pi/2*exp(-x)*(1/x + 3/x**2 + 3/x**3))
        let actual = x.map(|x| crate::sph_bessel_k(2, x));
        let expected = x.map(|x| FRAC_PI_2 * (-x).exp() * (1.0 + 3.0 / x + 3.0 / (x * x)) / x);
        crate::np_assert_allclose!(&actual, &expected);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalOld.test_sph_jn`
    #[test]
    fn test_sph_bessel_jn_prime() {
        // s1 = np.empty((2,3))
        let mut s1 = [[0.0_f64; 3]; 2];
        // x = 0.2
        let x = 0.2_f64;

        // s1[0][0] = spherical_jn(0, x)
        s1[0][0] = crate::sph_bessel_j(0, x);
        // s1[0][1] = spherical_jn(1, x)
        s1[0][1] = crate::sph_bessel_j(1, x);
        // s1[0][2] = spherical_jn(2, x)
        s1[0][2] = crate::sph_bessel_j(2, x);
        // s1[1][0] = spherical_jn(0, x, derivative=True)
        s1[1][0] = crate::sph_bessel_j_prime(0, x);
        // s1[1][1] = spherical_jn(1, x, derivative=True)
        s1[1][1] = crate::sph_bessel_j_prime(1, x);
        // s1[1][2] = spherical_jn(2, x, derivative=True)
        s1[1][2] = crate::sph_bessel_j_prime(2, x);

        // s10 = -s1[0][1]
        let s10 = -s1[0][1];
        // s11 = s1[0][0]-2.0/0.2*s1[0][1]
        let s11 = s1[0][0] - 2.0 / 0.2 * s1[0][1];
        // s12 = s1[0][1]-3.0/0.2*s1[0][2]
        let s12 = s1[0][1] - 3.0 / 0.2 * s1[0][2];
        // assert_allclose(s1[0], [0.99334665397530607731,
        //                         0.066400380670322230863,
        //                         0.0026590560795273856680],
        //                 atol=1.5e-12, rtol=0)
        crate::np_assert_allclose!(
            &s1[0],
            &[
                0.993_346_653_975_306_1,
                0.066_400_380_670_322_24,
                0.002_659_056_079_527_385_5,
            ],
            atol = 1.5e-12
        );
        // assert_allclose(s1[1], [s10, s11, s12], atol=1.5e-12, rtol=0)
        crate::np_assert_allclose!(&s1[1], &[s10, s11, s12], atol = 1.5e-12);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalOld.test_sph_yn`
    #[test]
    fn test_sph_bessel_yn_prime() {
        // sy1 = spherical_yn(2, 0.2)
        let sy1 = crate::sph_bessel_y(2, 0.2);
        // sy2 = spherical_yn(0, 0.2)
        let sy2 = crate::sph_bessel_y(0, 0.2);

        // # previous values in the system
        // assert_allclose(sy1, -377.52483, atol=1.5e-5, rtol=0)
        crate::np_assert_allclose!(&[sy1], &[-377.52483], atol = 1.5e-5);
        // assert_allclose(sy2, -4.9003329, atol=1.5e-5, rtol=0)
        crate::np_assert_allclose!(&[sy2], &[-4.900_332_9], atol = 1.5e-5);

        // sphpy = (spherical_yn(0, 0.2) - 2*spherical_yn(2, 0.2))/3
        let sphpy = (sy2 - 2.0 * sy1) / 3.0;
        // sy3 = spherical_yn(1, 0.2, derivative=True)
        let sy3 = crate::sph_bessel_y_prime(1, 0.2);

        // # compare correct derivative val. (correct =-system val).
        // assert_allclose(sy3, sphpy, atol=1.5e-4, rtol=0)
        crate::np_assert_allclose!(&[sy3], &[sphpy], atol = 1.5e-4);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalOld.test_sph_in`
    #[test]
    fn test_sph_bessel_in_prime() {
        // i1n = np.empty((2,2))
        let mut i1n = [[0.0_f64; 2]; 2];
        // x = 0.2
        let x = 0.2_f64;

        // i1n[0][0] = spherical_in(0, x)
        i1n[0][0] = crate::sph_bessel_i(0, x);
        // i1n[0][1] = spherical_in(1, x)
        i1n[0][1] = crate::sph_bessel_i(1, x);
        // i1n[1][0] = spherical_in(0, x, derivative=True)
        i1n[1][0] = crate::sph_bessel_i_prime(0, x);
        // i1n[1][1] = spherical_in(1, x, derivative=True)
        i1n[1][1] = crate::sph_bessel_i_prime(1, x);

        // inp0 = (i1n[0][1])
        let inp0 = i1n[0][1];
        // inp1 = (i1n[0][0] - 2.0/0.2 * i1n[0][1])
        let inp1 = i1n[0][0] - 2.0 / 0.2 * i1n[0][1];
        // assert_allclose(i1n[0], np.array([1.0066800127054699381,
        //                                   0.066933714568029540839]),
        //                 atol=1.5e-12, rtol=0.0)
        crate::np_assert_allclose!(
            &i1n[0],
            &[1.006_680_012_705_47, 0.066_933_714_568_029_54],
            atol = 1.5e-12
        );
        // assert_allclose(i1n[1], [inp0, inp1], atol=1.5e-12, rtol=0)
        crate::np_assert_allclose!(&i1n[1], &[inp0, inp1], atol = 1.5e-12);
    }

    /// From `scipy.special.tests.test_spherical_bessel.TestSphericalOld.test_sph_kn`
    #[test]
    fn test_sph_bessel_kn_prime() {
        // kn = np.empty((2,3))
        let mut kn = [[0.0_f64; 3]; 2];
        // x = 0.2
        let x = 0.2_f64;

        // kn[0][0] = spherical_kn(0, x)
        kn[0][0] = crate::sph_bessel_k(0, x);
        // kn[0][1] = spherical_kn(1, x)
        kn[0][1] = crate::sph_bessel_k(1, x);
        // kn[0][2] = spherical_kn(2, x)
        kn[0][2] = crate::sph_bessel_k(2, x);
        // kn[1][0] = spherical_kn(0, x, derivative=True)
        kn[1][0] = crate::sph_bessel_k_prime(0, x);
        // kn[1][1] = spherical_kn(1, x, derivative=True)
        kn[1][1] = crate::sph_bessel_k_prime(1, x);
        // kn[1][2] = spherical_kn(2, x, derivative=True)
        kn[1][2] = crate::sph_bessel_k_prime(2, x);

        // kn0 = -kn[0][1]
        let kn0 = -kn[0][1];
        // kn1 = -kn[0][0]-2.0/0.2*kn[0][1]
        let kn1 = -kn[0][0] - 2.0 / 0.2 * kn[0][1];
        // kn2 = -kn[0][1]-3.0/0.2*kn[0][2]
        let kn2 = -kn[0][1] - 3.0 / 0.2 * kn[0][2];
        // assert_allclose(kn[0], [6.4302962978445670140,
        //                         38.581777787067402086,
        //                         585.15696310385559829],
        //                 atol=1.5e-12, rtol=0)
        crate::np_assert_allclose!(
            &kn[0],
            &[
                6.430_296_297_844_567,
                38.581_777_787_067_4,
                585.156_963_103_855_6
            ],
            atol = 1.5e-12
        );
        // assert_allclose(kn[1], [kn0, kn1, kn2], atol=1.5e-9, rtol=0)
        crate::np_assert_allclose!(&kn[1], &[kn0, kn1, kn2], atol = 1.5e-9);
    }
}
