use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait FresnelArg: sealed::Sealed {
    type Output;
    fn fresnel(self) -> (Self::Output, Self::Output);
}

#[inline(always)]
fn c_c64_nan() -> bindings::root::std::complex<f64> {
    Complex::new(f64::NAN, f64::NAN).into()
}

impl FresnelArg for f64 {
    type Output = f64;

    #[inline(always)]
    fn fresnel(self) -> (Self::Output, Self::Output) {
        let mut fs = f64::NAN;
        let mut fc = f64::NAN;
        unsafe {
            bindings::fresnel(self, &mut fs, &mut fc);
        }
        (fs, fc)
    }
}

impl FresnelArg for Complex<f64> {
    type Output = Complex<f64>;

    #[inline(always)]
    fn fresnel(self) -> (Self::Output, Self::Output) {
        let mut fs = c_c64_nan();
        let mut fc = c_c64_nan();

        unsafe {
            bindings::fresnel_1(self.into(), &mut fs, &mut fc);
        }
        (fs.into(), fc.into())
    }
}

/// Fresnel integrals
///
/// The Fresnel integrals are defined as:
///
/// - S(z) = ∫₀ᶻ sin(π t² / 2) dt
/// - C(z) = ∫₀ᶻ cos(π t² / 2) dt
///
/// where the integrals are taken from 0 to z.
///
/// # Arguments
///
/// - `z` - real- or complex-valued argument
///
/// # Returns
///
/// - `S` - S(z)
/// - `C` - C(z)
///
/// # See also
///
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals
pub fn fresnel<T: FresnelArg>(z: T) -> (T::Output, T::Output) {
    z.fresnel()
}

/// Modified Fresnel positive integrals
///
/// Computes the modified Fresnel integrals:
///
/// - F₊(x) = ∫ₓ∞ exp(i t²) dt
/// - K₊(x) = F₊(x) × exp(-i (x² + π/4)) / √π
///
/// # Arguments
///
/// - `x` - real-valued argument
///
/// # Returns
///
/// - `fp` - Integral F₊(x)
/// - `kp` - Integral K₊(x)
///
/// # See also
///
/// - [`fresnel`] - Standard Fresnel integrals
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals
pub fn modified_fresnel_plus(x: f64) -> (Complex<f64>, Complex<f64>) {
    let mut fp = c_c64_nan();
    let mut kp = c_c64_nan();
    unsafe {
        bindings::modified_fresnel_plus(x, &mut fp, &mut kp);
    }
    (fp.into(), kp.into())
}

/// Modified Fresnel negative integrals
///
/// Computes the modified Fresnel integrals:
///
/// - F₋(x) = ∫ₓ∞ exp(-i t²) dt
/// - K₋(x) = F₋(x) × exp(i (x² + π/4)) / √π
///
/// # Arguments
///
/// - `x` - real-valued argument
///
/// # Returns
///
/// - `fm` - Integral F₋(x)
/// - `km` - Integral K₋(x)
///
/// # See also
///
/// - [`fresnel`] - Standard Fresnel integrals
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals
pub fn modified_fresnel_minus(x: f64) -> (Complex<f64>, Complex<f64>) {
    let mut fm = c_c64_nan();
    let mut km = c_c64_nan();
    unsafe {
        bindings::modified_fresnel_minus(x, &mut fm, &mut km);
    }
    (fm.into(), km.into())
}
