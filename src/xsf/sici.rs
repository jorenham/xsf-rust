use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait SiciArg: sealed::Sealed + Sized {
    fn sici(self) -> (Self, Self);
    fn shichi(self) -> (Self, Self);
}

impl SiciArg for f64 {
    #[inline(always)]
    fn sici(self) -> (Self, Self) {
        let (mut si, mut ci) = (f64::NAN, f64::NAN);
        unsafe {
            crate::ffi::xsf::sici(self, &mut si, &mut ci);
        }
        (si, ci)
    }

    #[inline(always)]
    fn shichi(self) -> (Self, Self) {
        let (mut shi, mut chi) = (f64::NAN, f64::NAN);
        unsafe {
            crate::ffi::xsf::shichi(self, &mut shi, &mut chi);
        }
        (shi, chi)
    }
}

impl SiciArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn sici(self) -> (Self, Self) {
        let (mut si, mut ci) = (Complex::new(f64::NAN, 0.0), Complex::new(f64::NAN, 0.0));
        unsafe {
            crate::ffi::xsf::sici_1(self, &mut si, &mut ci);
        }
        (si, ci)
    }

    #[inline(always)]
    fn shichi(self) -> (Self, Self) {
        let (mut shi, mut chi) = (Complex::new(f64::NAN, 0.0), Complex::new(f64::NAN, 0.0));
        unsafe {
            crate::ffi::xsf::shichi_1(self, &mut shi, &mut chi);
        }
        (shi, chi)
    }
}

/// Sine and cosine integrals.
///
/// The sine integral is:
///
/// Si(z) = ∫₀ᶻ (sin t)/t dt
///
/// The cosine integral is:
///
/// Ci(z) = γ + ln(z) + ∫₀ᶻ (cos t - 1)/t dt
///
/// where γ is Euler's constant and ln is the principal branch of the logarithm.
///
/// # Arguments
///
/// - `z` - real- or complex-valued argument
///
/// # Returns
///
/// - `Si(z)` - Sine integral
/// - `Ci(z)` - Cosine integral
pub fn sici<T: SiciArg>(z: T) -> (T, T) {
    z.sici()
}

/// Hyperbolic sine and cosine integrals.
///
/// The hyperbolic sine integral is:
///
/// Shi(z) = ∫₀ᶻ (sinh t)/t dt
///
/// The hyperbolic cosine integral is:
///
/// Chi(z) = γ + ln(z) + ∫₀ᶻ (cosh t - 1)/t dt
///
/// where γ is Euler's constant and ln is the principal branch of the logarithm.
///
/// # Arguments
///
/// - `z` - real- or complex-valued argument
///
/// # Returns
///
/// - `Shi(z)` - Hyperbolic sine integral
/// - `Chi(z)` - Hyperbolic cosine integral
pub fn shichi<T: SiciArg>(z: T) -> (T, T) {
    z.shichi()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_sici_f64() {
        crate::xsref::test("sici", "d-d_d", |x| crate::sici(x[0]));
    }

    #[test]
    fn test_sici_c64() {
        crate::xsref::test("sici", "cd-cd_cd", |x| crate::sici(c64(x[0], x[1])));
    }

    #[test]
    fn test_shichi_f64() {
        crate::xsref::test("shichi", "d-d_d", |x| crate::shichi(x[0]));
    }

    #[test]
    fn test_shichi_c64() {
        crate::xsref::test("shichi", "cd-cd_cd", |x| crate::shichi(c64(x[0], x[1])));
    }
}
