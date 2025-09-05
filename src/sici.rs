use crate::bindings;
use num_complex::Complex;

#[inline(always)]
fn c_c64_nan() -> bindings::root::std::complex<f64> {
    Complex::new(f64::NAN, f64::NAN).into()
}

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait SiciArg: sealed::Sealed {
    type Output;

    fn sici(self) -> (Self::Output, Self::Output);
    fn shichi(self) -> (Self::Output, Self::Output);
}

impl SiciArg for f64 {
    type Output = f64;

    #[inline(always)]
    fn sici(self) -> (Self::Output, Self::Output) {
        let mut si = f64::NAN;
        let mut ci = f64::NAN;

        unsafe {
            bindings::sici(self, &mut si, &mut ci);
        }
        (si, ci)
    }

    #[inline(always)]
    fn shichi(self) -> (Self::Output, Self::Output) {
        let mut shi = f64::NAN;
        let mut chi = f64::NAN;

        unsafe {
            bindings::shichi(self, &mut shi, &mut chi);
        }
        (shi, chi)
    }
}

impl SiciArg for Complex<f64> {
    type Output = Complex<f64>;

    #[inline(always)]
    fn sici(self) -> (Self::Output, Self::Output) {
        let mut si = c_c64_nan();
        let mut ci = c_c64_nan();

        unsafe {
            bindings::sici_1(self.into(), &mut si, &mut ci);
        }
        (si.into(), ci.into())
    }

    #[inline(always)]
    fn shichi(self) -> (Self::Output, Self::Output) {
        let mut shi = c_c64_nan();
        let mut chi = c_c64_nan();

        unsafe {
            bindings::shichi_1(self.into(), &mut shi, &mut chi);
        }
        (shi.into(), chi.into())
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
pub fn sici<T: SiciArg>(z: T) -> (T::Output, T::Output) {
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
pub fn shichi<T: SiciArg>(z: T) -> (T::Output, T::Output) {
    z.shichi()
}
