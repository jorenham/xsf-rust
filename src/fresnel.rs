use crate::bindings;
use alloc::vec::Vec;
use core::ffi::c_int;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait FresnelArg: sealed::Sealed {
    type Output;
    fn fresnel(self) -> (Self::Output, Self::Output);
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
        let mut fs = bindings::complex_nan();
        let mut fc = bindings::complex_nan();
        unsafe {
            bindings::fresnel_1(self.into(), &mut fs, &mut fc);
        }
        (fs.into(), fc.into())
    }
}

/// Fresnel integrals S(z) and C(z) for real or complex argument
///
/// The Fresnel integrals are defined as:
/// - S(z) = ∫₀ᶻ sin(π t² / 2) dt
/// - C(z) = ∫₀ᶻ cos(π t² / 2) dt
///
/// where the integrals are taken from 0 to z.
///
/// # See also
/// - [`fresnel_zeros`] - Zeros of S(z) and C(z)
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals
pub fn fresnel<T: FresnelArg>(z: T) -> (T::Output, T::Output) {
    z.fresnel()
}

/// Modified Fresnel positive integrals
///
/// Computes the modified Fresnel integrals:
/// - F₊(x) = ∫ₓ∞ exp(i t²) dt
/// - K₊(x) = F₊(x) × exp(-i (x² + π/4)) / √π
///
/// # See also
/// - [`fresnel`] - Standard Fresnel integrals
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals
pub fn modified_fresnel_plus(x: f64) -> (Complex<f64>, Complex<f64>) {
    let mut fp = bindings::complex_nan();
    let mut kp = bindings::complex_nan();
    unsafe {
        bindings::modified_fresnel_plus(x, &mut fp, &mut kp);
    }
    (fp.into(), kp.into())
}

/// Modified Fresnel negative integrals
///
/// Computes the modified Fresnel integrals:
/// - F₋(x) = ∫ₓ∞ exp(-i t²) dt
/// - K₋(x) = F₋(x) × exp(i (x² + π/4)) / √π
///
/// # See also
/// - [`fresnel`] - Standard Fresnel integrals S(z) and C(z)
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals
pub fn modified_fresnel_minus(x: f64) -> (Complex<f64>, Complex<f64>) {
    let mut fm = bindings::complex_nan();
    let mut km = bindings::complex_nan();
    unsafe {
        bindings::modified_fresnel_minus(x, &mut fm, &mut km);
    }
    (fm.into(), km.into())
}

/// Zeros of Fresnel integrals S(z) and C(z)
///
/// Compute `nt` complex zeros of the Fresnel integrals S(z) and C(z).
///
/// # See also
/// - [`fresnel`] - Fresnel integrals S(z) and C(z)
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals F₊(x) and K₊(x)
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals F₋(x) and K₋(x)
pub fn fresnel_zeros(nt: usize) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    assert!(nt <= c_int::MAX as usize);

    let mut szo = bindings::complex_zeros(nt);
    let mut czo = bindings::complex_zeros(nt);
    unsafe {
        bindings::fcszo(2, nt as c_int, szo.as_mut_ptr());
        bindings::fcszo(1, nt as c_int, czo.as_mut_ptr());
    }
    (bindings::cvec_into(szo), bindings::cvec_into(czo))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    // fresnel

    #[test]
    fn test_fresnel_f64() {
        testing::test::<(f64, f64), _>("fresnel", "d-d_d", |x: &[f64]| fresnel(x[0]));
    }

    #[test]
    fn test_fresnel_c64() {
        testing::test::<(Complex<f64>, Complex<f64>), _>("fresnel", "cd-cd_cd", |x: &[f64]| {
            fresnel(c64(x[0], x[1]))
        });
    }

    // modified_fresnel_plus

    #[test]
    fn test_modified_fresnel_plus_c64() {
        testing::test::<(Complex<f64>, Complex<f64>), _>(
            "modified_fresnel_plus",
            "d-cd_cd",
            |x: &[f64]| modified_fresnel_plus(x[0]),
        );
    }

    // modified_fresnel_minus

    #[test]
    fn test_modified_fresnel_minus_c64() {
        testing::test::<(Complex<f64>, Complex<f64>), _>(
            "modified_fresnel_minus",
            "d-cd_cd",
            |x: &[f64]| modified_fresnel_minus(x[0]),
        );
    }

    // fresnel_zeros

    fn assert_allclose(actual: &[Complex<f64>], expected: &[Complex<f64>], atol: f64) {
        assert_eq!(actual.len(), expected.len());
        for (&a, &e) in actual.iter().zip(expected.iter()) {
            assert!((a - e).norm() <= atol);
        }
    }

    /// Ported from `scipy.special.tests.test_basic.TestFresnel.test_fresnel_zeros`
    #[test]
    fn test_fresnel_zeros() {
        // szo, czo = special.fresnel_zeros(5)
        let (szo, czo) = fresnel_zeros(5);
        assert_allclose(
            &szo,
            &[
                c64(2.0093, 0.2885),
                c64(2.8335, 0.2443),
                c64(3.4675, 0.2185),
                c64(4.0026, 0.2009),
                c64(4.4742, 0.1877),
            ],
            1.5e-3,
        );
        assert_allclose(
            &czo,
            &[
                c64(1.7437, 0.3057),
                c64(2.6515, 0.2529),
                c64(3.3204, 0.2240),
                c64(3.8757, 0.2047),
                c64(4.3611, 0.1907),
            ],
            1.5e-3,
        );

        // vals1 = special.fresnel(szo)[0]
        let vals1 = szo.iter().map(|&z| fresnel(z).0);
        // vals2 = special.fresnel(czo)[1]
        let vals2 = czo.iter().map(|&z| fresnel(z).1);
        // assert_allclose(vals1, 0, atol=1.5e-14, rtol=0)
        for val1 in vals1 {
            assert!(val1.norm() < 1.5e-14);
        }
        // assert_allclose(vals2, 0, atol=1.5e-14, rtol=0)
        for val2 in vals2 {
            assert!(val2.norm() < 1.5e-14);
        }
    }
}
