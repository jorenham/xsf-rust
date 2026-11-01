//! The `voigt_profile` is translated into pure Rust from `xsf/erf.h` (xsf v0.2.2), and the other
//! functions use `xsf/faddeeva.rs`.

use core::f64::consts::{FRAC_1_SQRT_2, PI};

use num_complex::{Complex64, c64};

use crate::utils::SQRT_TAU;

pub trait ErfArg: crate::sealed::Sealed {
    fn erf(self) -> Self;
    fn erfc(self) -> Self;
    fn erfcx(self) -> Self;
    fn erfi(self) -> Self;
    fn dawsn(self) -> Self;
}

impl ErfArg for f64 {
    #[inline]
    fn erf(self) -> Self {
        crate::xsf::cephes::erf(self)
    }
    #[inline]
    fn erfc(self) -> Self {
        crate::xsf::cephes::erfc(self)
    }
    #[inline]
    fn erfcx(self) -> Self {
        crate::xsf::faddeeva::erfcx(self)
    }
    #[inline]
    fn erfi(self) -> Self {
        crate::xsf::faddeeva::erfi(self)
    }
    #[inline]
    fn dawsn(self) -> Self {
        crate::xsf::faddeeva::dawson(self)
    }
}

impl ErfArg for Complex64 {
    #[inline]
    fn erf(self) -> Self {
        crate::xsf::faddeeva::erf_c64(self, 0.0)
    }
    #[inline]
    fn erfc(self) -> Self {
        crate::xsf::faddeeva::erfc_c64(self, 0.0)
    }
    #[inline]
    fn erfcx(self) -> Self {
        crate::xsf::faddeeva::erfcx_c64(self, 0.0)
    }
    #[inline]
    fn erfi(self) -> Self {
        crate::xsf::faddeeva::erfi_c64(self, 0.0)
    }
    #[inline]
    fn dawsn(self) -> Self {
        crate::xsf::faddeeva::dawson_c64(self, 0.0)
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
/// Zeros in the other quadrants can be obtained by using the symmetries *erf(-z) = -erf(z)* and
/// *erf(conj(z)) = conj(erf(z))*.
///
/// Corresponds to [`scipy.special.erf_zeros`][erf_zeros].
///
/// [erf_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf_zeros.html
///
/// # Panics
/// - if `N` is too large to fit into a [`c_int`](core::ffi::c_int)
///
/// # See also
/// - [`erf`]: Error function
/// - [`erfinv`](crate::erfinv): Inverse error function
#[doc(alias = "cerzo")]
#[must_use]
#[inline]
pub fn erf_zeros<const N: usize>() -> [Complex64; N] {
    let mut zo = [f64::NAN.into(); N];
    crate::xsf::specfun::cerzo(N.try_into().unwrap(), &mut zo);
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
    crate::xsf::faddeeva::w(z, 0.0)
}

/// Voigt profile
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
pub fn voigt_profile(x: f64, sigma: f64, gamma: f64) -> f64 {
    const INV_SQRT_2: f64 = FRAC_1_SQRT_2;

    if sigma == 0.0 {
        if gamma == 0.0 {
            if x.is_nan() {
                return x;
            }
            if x == 0.0 {
                return f64::INFINITY;
            }
            return 0.0;
        }
        return gamma / PI / (x * x + gamma * gamma);
    }
    if gamma == 0.0 {
        return 1.0 / SQRT_TAU / sigma * (-(x / sigma) * (x / sigma) / 2.0).exp();
    }

    let zreal = x / sigma * INV_SQRT_2;
    let zimag = gamma / sigma * INV_SQRT_2;
    let z = c64(zreal, zimag);
    let w = crate::xsf::faddeeva::w(z, 0.0);
    w.re / sigma / SQRT_TAU
}
