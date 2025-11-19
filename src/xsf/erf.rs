use core::ffi::c_int;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ErfArg: sealed::Sealed {
    fn erf(self) -> Self;
    fn erfc(self) -> Self;
    fn erfcx(self) -> Self;
    fn erfi(self) -> Self;
    fn dawsn(self) -> Self;
}

impl ErfArg for f64 {
    #[inline(always)]
    fn erf(self) -> Self {
        unsafe { crate::ffi::xsf::erf(self) }
    }
    #[inline(always)]
    fn erfc(self) -> Self {
        unsafe { crate::ffi::xsf::erfc(self) }
    }
    #[inline(always)]
    fn erfcx(self) -> Self {
        unsafe { crate::ffi::xsf::erfcx(self) }
    }
    #[inline(always)]
    fn erfi(self) -> Self {
        unsafe { crate::ffi::xsf::erfi(self) }
    }
    #[inline(always)]
    fn dawsn(self) -> Self {
        unsafe { crate::ffi::xsf::dawsn(self) }
    }
}

impl ErfArg for Complex<f64> {
    #[inline(always)]
    fn erf(self) -> Self {
        unsafe { crate::ffi::xsf::erf_1(self) }
    }
    #[inline(always)]
    fn erfc(self) -> Self {
        unsafe { crate::ffi::xsf::erfc_1(self) }
    }
    #[inline(always)]
    fn erfcx(self) -> Self {
        unsafe { crate::ffi::xsf::erfcx_1(self) }
    }
    #[inline(always)]
    fn erfi(self) -> Self {
        unsafe { crate::ffi::xsf::erfi_1(self) }
    }
    #[inline(always)]
    fn dawsn(self) -> Self {
        unsafe { crate::ffi::xsf::dawsn_1(self) }
    }
}

/// Error function `erf(z)` for real or complex input
pub fn erf<T: ErfArg>(z: T) -> T {
    z.erf()
}

/// Complementary error function `1 - erf(z)` for real or complex input
pub fn erfc<T: ErfArg>(z: T) -> T {
    z.erfc()
}

/// Scaled complementary error function `exp(z^2) * erfc(z)` for real or complex input
pub fn erfcx<T: ErfArg>(z: T) -> T {
    z.erfcx()
}

/// Imaginary error function `-i erf(i z)` for real or complex input
pub fn erfi<T: ErfArg>(z: T) -> T {
    z.erfi()
}

/// Zeros of the error function `erf(z)` in the first quadrant of the complex plane
///
/// Zeros in the other quadrants can be obtained by using the symmetries *erf(-z) = erf(z)* and
/// *erf(conj(z)) = conj(erf(z))*.
///
/// See [`erf`](crate::erf) for the error function itself.
#[doc(alias = "cerzo")]
pub fn erf_zeros<const NT: usize>() -> [Complex<f64>; NT] {
    let mut zo = [Complex::new(f64::NAN, 0.0); NT];
    unsafe {
        crate::ffi::xsf::cerzo(NT as c_int, zo.as_mut_ptr());
    }
    zo
}

/// Dawson function `sqrt(pi)/2 * exp(-z^2) * erfi(z)` for real or complex input
#[doc(alias = "dawson")]
pub fn dawsn<T: ErfArg>(z: T) -> T {
    z.dawsn()
}

/// Faddeeva function `exp(-z^2) * erfc(-i z)`
#[doc(alias = "faddeeva")]
pub fn wofz(z: Complex<f64>) -> Complex<f64> {
    unsafe { crate::ffi::xsf::wofz(z) }
}

/// Voigt profile
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
        // erz = special.erf_zeros(5)
        let erz = crate::erf_zeros::<5>();
        // erzr = array([1.45061616+1.88094300j,
        //               2.24465928+2.61657514j,
        //               2.83974105+3.17562810j,
        //               3.33546074+3.64617438j,
        //               3.76900557+4.06069723j])
        let erzr = [
            c64(1.45061616, 1.88094300),
            c64(2.24465928, 2.61657514),
            c64(2.83974105, 3.17562810),
            c64(3.33546074, 3.64617438),
            c64(3.76900557, 4.06069723),
        ];
        // assert_allclose(erz, erzr, atol=1.5e-4, rtol=0)
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
