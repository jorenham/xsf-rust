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
    fn erf(self) -> f64 {
        unsafe { crate::ffi::xsf::erf(self) }
    }
    #[inline(always)]
    fn erfc(self) -> f64 {
        unsafe { crate::ffi::xsf::erfc(self) }
    }
    #[inline(always)]
    fn erfcx(self) -> f64 {
        unsafe { crate::ffi::xsf::erfcx(self) }
    }
    #[inline(always)]
    fn erfi(self) -> f64 {
        unsafe { crate::ffi::xsf::erfi(self) }
    }
    #[inline(always)]
    fn dawsn(self) -> f64 {
        unsafe { crate::ffi::xsf::dawsn(self) }
    }
}

impl ErfArg for Complex<f64> {
    #[inline(always)]
    fn erf(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::erf_1(self.into()) }.into()
    }
    #[inline(always)]
    fn erfc(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::erfc_1(self.into()) }.into()
    }
    #[inline(always)]
    fn erfcx(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::erfcx_1(self.into()) }.into()
    }
    #[inline(always)]
    fn erfi(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::erfi_1(self.into()) }.into()
    }
    #[inline(always)]
    fn dawsn(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::dawsn_1(self.into()) }.into()
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

/// Dawson function `sqrt(pi)/2 * exp(-z^2) * erfi(z)` for real or complex input
pub fn dawsn<T: ErfArg>(z: T) -> T {
    z.dawsn()
}

/// Faddeeva function `exp(-z^2) * erfc(-i z)`
pub fn wofz(z: Complex<f64>) -> Complex<f64> {
    unsafe { crate::ffi::xsf::wofz(z.into()) }.into()
}

/// Voigt profile
pub fn voigt_profile(x: f64, sigma: f64, gamma: f64) -> f64 {
    unsafe { crate::ffi::xsf::voigt_profile(x, sigma, gamma) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    // erf

    #[test]
    fn test_erf_f64() {
        testing::test::<f64, _>("erf", "d-d", |x: &[f64]| erf(x[0]));
    }

    #[test]
    fn test_erf_c64() {
        testing::test::<Complex<f64>, _>("erf", "cd-cd", |x: &[f64]| erf(c64(x[0], x[1])));
    }

    // erfc

    #[test]
    fn test_erfc_f64() {
        testing::test::<f64, _>("erfc", "d-d", |x: &[f64]| erfc(x[0]));
    }

    #[test]
    fn test_erfc_c64() {
        testing::test::<Complex<f64>, _>("erfc", "cd-cd", |x: &[f64]| erfc(c64(x[0], x[1])));
    }

    // erfcx

    #[test]
    fn test_erfcx_f64() {
        testing::test::<f64, _>("erfcx", "d-d", |x: &[f64]| erfcx(x[0]));
    }

    #[test]
    fn test_erfcx_c64() {
        testing::test::<Complex<f64>, _>("erfcx", "cd-cd", |x: &[f64]| erfcx(c64(x[0], x[1])));
    }

    // erfi

    #[test]
    fn test_erfi_f64() {
        testing::test::<f64, _>("erfi", "d-d", |x: &[f64]| erfi(x[0]));
    }

    #[test]
    fn test_erfi_c64() {
        testing::test::<Complex<f64>, _>("erfi", "cd-cd", |x: &[f64]| erfi(c64(x[0], x[1])));
    }

    // dawsn

    #[test]
    fn test_dawsn_f64() {
        testing::test::<f64, _>("dawsn", "d-d", |x: &[f64]| dawsn(x[0]));
    }

    #[test]
    fn test_dawsn_c64() {
        testing::test::<Complex<f64>, _>("dawsn", "cd-cd", |x: &[f64]| dawsn(c64(x[0], x[1])));
    }

    // wofz

    #[test]
    fn test_wofz() {
        testing::test::<Complex<f64>, _>("wofz", "cd-cd", |x: &[f64]| wofz(c64(x[0], x[1])));
    }

    // voigt_profile

    #[test]
    fn test_voigt_profile() {
        testing::test::<f64, _>("voigt_profile", "d_d_d-d", |x: &[f64]| {
            voigt_profile(x[0], x[1], x[2])
        });
    }
}
