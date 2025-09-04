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
