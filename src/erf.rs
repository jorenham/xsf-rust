use crate::bindings;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

macro_rules! define_erf_functions {
    ($(($method:ident, $binding:ident, $doc:expr)),*) => {
        pub trait ErfArg: sealed::Sealed {
            type Output;
            $(
                fn $method(self) -> Self::Output;
            )*
        }

        impl ErfArg for f64 {
            type Output = f64;
            $(
                fn $method(self) -> Self::Output {
                    unsafe { bindings::$binding(self) }
                }
            )*
        }

        impl ErfArg for Complex<f64> {
            type Output = Complex<f64>;
            $(
                fn $method(self) -> Self::Output {
                    paste::paste! {
                        unsafe { bindings::[<$binding _1>](self.into()) }.into()
                    }
                }
            )*
        }

        $(
            #[doc = $doc]
            pub fn $method<T: ErfArg>(z: T) -> T::Output {
                z.$method()
            }
        )*
    };
}

define_erf_functions! {
    (erf, erf_, "Error function `erf(z)`"),
    (erfc, erfc_, "Complementary error function `1 - erf(z)`"),
    (erfcx, erfcx, "Scaled complementary error function `exp(z^2) * erfc(z)`"),
    (erfi, erfi, "Imaginary error function `-i erf(i z)`"),
    (dawsn, dawsn, "Dawson function `sqrt(pi)/2 * exp(-z^2) * erfi(z)`")
}

/// Faddeeva function `exp(-z^2) * erfc(-i z)`
pub fn wofz(z: Complex<f64>) -> Complex<f64> {
    unsafe { bindings::wofz(z.into()) }.into()
}

/// Voigt profile
pub fn voigt_profile(x: f64, sigma: f64, gamma: f64) -> f64 {
    unsafe { bindings::voigt_profile(x, sigma, gamma) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::{Complex, c64};

    // erf

    #[test]
    fn test_erf_f64() {
        xsref::test::<f64, _>("erf", "d-d", |x: &[f64]| erf(x[0]));
    }

    #[test]
    fn test_erf_c64() {
        xsref::test::<Complex<f64>, _>("erf", "cd-cd", |x: &[f64]| erf(c64(x[0], x[1])));
    }

    // erfc

    #[test]
    fn test_erfc_f64() {
        xsref::test::<f64, _>("erfc", "d-d", |x: &[f64]| erfc(x[0]));
    }

    #[test]
    fn test_erfc_c64() {
        xsref::test::<Complex<f64>, _>("erfc", "cd-cd", |x: &[f64]| erfc(c64(x[0], x[1])));
    }

    // erfcx

    #[test]
    fn test_erfcx_f64() {
        xsref::test::<f64, _>("erfcx", "d-d", |x: &[f64]| erfcx(x[0]));
    }

    #[test]
    fn test_erfcx_c64() {
        xsref::test::<Complex<f64>, _>("erfcx", "cd-cd", |x: &[f64]| erfcx(c64(x[0], x[1])));
    }

    // erfi

    #[test]
    fn test_erfi_f64() {
        xsref::test::<f64, _>("erfi", "d-d", |x: &[f64]| erfi(x[0]));
    }

    #[test]
    fn test_erfi_c64() {
        xsref::test::<Complex<f64>, _>("erfi", "cd-cd", |x: &[f64]| erfi(c64(x[0], x[1])));
    }

    // dawsn

    #[test]
    fn test_dawsn_f64() {
        xsref::test::<f64, _>("dawsn", "d-d", |x: &[f64]| dawsn(x[0]));
    }

    #[test]
    fn test_dawsn_c64() {
        xsref::test::<Complex<f64>, _>("dawsn", "cd-cd", |x: &[f64]| dawsn(c64(x[0], x[1])));
    }

    // wofz

    #[test]
    fn test_wofz() {
        xsref::test::<Complex<f64>, _>("wofz", "cd-cd", |x: &[f64]| wofz(c64(x[0], x[1])));
    }

    // voigt_profile

    #[test]
    fn test_voigt_profile() {
        xsref::test::<f64, _>("voigt_profile", "d_d_d-d", |x: &[f64]| {
            voigt_profile(x[0], x[1], x[2])
        });
    }
}
