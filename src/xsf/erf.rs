use num_complex::Complex64;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex64 {}
}

pub trait ErfArg: sealed::Sealed {
    fn erf(self) -> Self;
    fn erfc(self) -> Self;
    fn erfcx(self) -> Self;
    fn erfi(self) -> Self;
    fn dawsn(self) -> Self;
}

impl ErfArg for f64 {
    #[inline]
    fn erf(self) -> Self {
        unsafe { crate::ffi::xsf::erf(self) }
    }
    #[inline]
    fn erfc(self) -> Self {
        unsafe { crate::ffi::xsf::erfc(self) }
    }
    #[inline]
    fn erfcx(self) -> Self {
        unsafe { crate::ffi::xsf::erfcx(self) }
    }
    #[inline]
    fn erfi(self) -> Self {
        unsafe { crate::ffi::xsf::erfi(self) }
    }
    #[inline]
    fn dawsn(self) -> Self {
        unsafe { crate::ffi::xsf::dawsn(self) }
    }
}

impl ErfArg for Complex64 {
    #[inline]
    fn erf(self) -> Self {
        unsafe { crate::ffi::xsf::erf_1(self) }
    }
    #[inline]
    fn erfc(self) -> Self {
        unsafe { crate::ffi::xsf::erfc_1(self) }
    }
    #[inline]
    fn erfcx(self) -> Self {
        unsafe { crate::ffi::xsf::erfcx_1(self) }
    }
    #[inline]
    fn erfi(self) -> Self {
        unsafe { crate::ffi::xsf::erfi_1(self) }
    }
    #[inline]
    fn dawsn(self) -> Self {
        unsafe { crate::ffi::xsf::dawsn_1(self) }
    }
}

/// Error function $\erf(z)$ for real or complex input
///
/// Corresponds to [`scipy.special.erf`][erf].
///
/// [erf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf.html
///
/// # See also
/// - [`erfc`]: Complementary error function
/// - [`erfi`]: Imaginary error function
/// - [`erfinv`](crate::erfinv): Inverse error function
/// - [`erf_zeros`]: Zeros of the error function
#[inline]
pub fn erf<T: ErfArg>(z: T) -> T {
    z.erf()
}

/// Complementary error function $1 - \erf(z)$ for real or complex input
///
/// Corresponds to [`scipy.special.erfc`][erfc].
///
/// [erfc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfc.html
///
/// # See also
/// - [`erf`]: Error function
/// - [`erfcinv`](crate::erfcinv): Inverse of the complementary error function
/// - [`erfcx`]: Scaled complementary error function
#[inline]
pub fn erfc<T: ErfArg>(z: T) -> T {
    z.erfc()
}

/// Scaled complementary error function $e^{z^2} \erfc(z)$ for real or complex input
///
/// Corresponds to [`scipy.special.erfcx`][erfcx].
///
/// [erfcx]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcx.html
///
/// # See also
/// - [`erfc`]: Complementary error function
/// - [`erf`]: Error function
#[inline]
pub fn erfcx<T: ErfArg>(z: T) -> T {
    z.erfcx()
}

/// Imaginary error function `-i erf(i z)` for real or complex input
///
/// Corresponds to [`scipy.special.erfi`][erfi].
///
/// [erfi]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfi.html
///
/// # See also
/// - [`erf`]: Error function
#[inline]
pub fn erfi<T: ErfArg>(z: T) -> T {
    z.erfi()
}

/// Zeros of the error function `erf(z)` in the first quadrant of the complex plane
///
/// Zeros in the other quadrants can be obtained by using the symmetries *erf(-z) = erf(z)* and
/// *erf(conj(z)) = conj(erf(z))*.
///
/// See [`erf`](crate::erf) for the error function itself.
///
/// Corresponds to [`scipy.special.erf_zeros`][erf_zeros].
///
/// [erf_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf_zeros.html
///
/// # Panics
/// - if `N` is too large to fit into a [`c_int`](core::ffi::c_int)
///
/// # See also
/// - [`erf`]: Error functionq
/// - [`erfinv`](crate::erfinv): Inverse error function
#[doc(alias = "cerzo")]
#[must_use]
#[inline]
pub fn erf_zeros<const N: usize>() -> [Complex64; N] {
    let mut zo = [f64::NAN.into(); N];
    unsafe {
        crate::ffi::xsf::cerzo(N.try_into().unwrap(), zo.as_mut_ptr());
    }
    zo
}

/// Dawson function `sqrt(pi)/2 * exp(-z^2) * erfi(z)` for real or complex input
#[doc(alias = "dawson")]
#[must_use]
#[inline]
pub fn dawsn<T: ErfArg>(z: T) -> T {
    z.dawsn()
}

/// Faddeeva function `exp(-z^2) * erfc(-i z)`
#[doc(alias = "faddeeva")]
#[must_use]
#[inline]
pub fn wofz(z: Complex64) -> Complex64 {
    unsafe { crate::ffi::xsf::wofz(z) }
}

/// Voigt profile
#[must_use]
#[inline]
pub fn voigt_profile(x: f64, sigma: f64, gamma: f64) -> f64 {
    unsafe { crate::ffi::xsf::voigt_profile(x, sigma, gamma) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    // erf

    #[test]
    fn test_erf_f64() {
        xsref::test("erf", "d-d", |x| crate::erf(x[0]));
    }

    #[test]
    fn test_erf_c64() {
        xsref::test("erf", "cd-cd", |x| crate::erf(c64(x[0], x[1])));
    }

    // erfc

    #[test]
    fn test_erfc_f64() {
        xsref::test("erfc", "d-d", |x| crate::erfc(x[0]));
    }

    #[test]
    fn test_erfc_c64() {
        xsref::test("erfc", "cd-cd", |x| crate::erfc(c64(x[0], x[1])));
    }

    // erfcx

    #[test]
    fn test_erfcx_f64() {
        xsref::test("erfcx", "d-d", |x| crate::erfcx(x[0]));
    }

    #[test]
    fn test_erfcx_c64() {
        xsref::test("erfcx", "cd-cd", |x| crate::erfcx(c64(x[0], x[1])));
    }

    // erfi

    #[test]
    fn test_erfi_f64() {
        xsref::test("erfi", "d-d", |x| crate::erfi(x[0]));
    }

    #[test]
    fn test_erfi_c64() {
        xsref::test("erfi", "cd-cd", |x| crate::erfi(c64(x[0], x[1])));
    }

    // erf_zeros (from `scipy.speceial.tests.test_basic.TestErf.test_erf_zeros`)

    #[test]
    fn test_erf_zeros() {
        let erz = crate::erf_zeros::<5>();
        let erzr = [
            c64(1.450_616_16, 1.880_943_00),
            c64(2.244_659_28, 2.616_575_14),
            c64(2.839_741_05, 3.175_628_10),
            c64(3.335_460_74, 3.646_174_38),
            c64(3.769_005_57, 4.060_697_23),
        ];
        crate::np_assert_allclose!(&erz, &erzr, atol = 1.5e-4);
    }

    // dawsn

    #[test]
    fn test_dawsn_f64() {
        xsref::test("dawsn", "d-d", |x| crate::dawsn(x[0]));
    }

    #[test]
    fn test_dawsn_c64() {
        xsref::test("dawsn", "cd-cd", |x| crate::dawsn(c64(x[0], x[1])));
    }

    // wofz

    #[test]
    fn test_wofz() {
        xsref::test("wofz", "cd-cd", |x| crate::wofz(c64(x[0], x[1])));
    }

    // voigt_profile

    #[test]
    fn test_voigt_profile() {
        xsref::test("voigt_profile", "d_d_d-d", |x| {
            crate::voigt_profile(x[0], x[1], x[2])
        });
    }
}
