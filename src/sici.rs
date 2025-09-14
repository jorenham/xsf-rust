use crate::bindings;
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
        let mut si = f64::NAN;
        let mut ci = f64::NAN;

        unsafe {
            bindings::sici(self, &mut si, &mut ci);
        }
        (si, ci)
    }

    #[inline(always)]
    fn shichi(self) -> (Self, Self) {
        let mut shi = f64::NAN;
        let mut chi = f64::NAN;

        unsafe {
            bindings::shichi(self, &mut shi, &mut chi);
        }
        (shi, chi)
    }
}

impl SiciArg for Complex<f64> {
    #[inline(always)]
    fn sici(self) -> (Self, Self) {
        let mut si = f64::NAN.into();
        let mut ci = f64::NAN.into();

        unsafe {
            bindings::sici_1(self.into(), &mut si, &mut ci);
        }
        (si.into(), ci.into())
    }

    #[inline(always)]
    fn shichi(self) -> (Self, Self) {
        let mut shi = f64::NAN.into();
        let mut chi = f64::NAN.into();

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
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    // sici

    #[test]
    fn test_sici_f64() {
        testing::test::<(f64, f64), _>("sici", "d-d_d", |x: &[f64]| sici(x[0]));
    }

    #[test]
    fn test_sici_c64() {
        testing::test::<(Complex<f64>, Complex<f64>), _>("sici", "cd-cd_cd", |x: &[f64]| {
            sici(c64(x[0], x[1]))
        });
    }

    // shichi

    #[test]
    fn test_shichi_f64() {
        testing::test::<(f64, f64), _>("shichi", "d-d_d", |x: &[f64]| shichi(x[0]));
    }

    #[test]
    fn test_shichi_c64() {
        testing::test::<(Complex<f64>, Complex<f64>), _>("shichi", "cd-cd_cd", |x: &[f64]| {
            shichi(c64(x[0], x[1]))
        });
    }
}
