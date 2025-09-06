use crate::bindings;
use crate::bindings::xsf_impl;
use num_complex::Complex;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

macro_rules! impl_methods {
    (bessel for f64: { $($fn_name:ident),* }) => {
        $(
            fn $fn_name(self, v: f64) -> f64 {
                unsafe { bindings::$fn_name(v, self) }
            }
        )*
    };

    (bessel for Complex<f64>: { $($fn_name:ident),* }) => {
        $(
            fn $fn_name(self, v: f64) -> Complex<f64> {
                paste::paste! {
                    unsafe { bindings::[<$fn_name _1>](v, self.into()) }.into()
                }
            }
        )*
    };

    (hankel for f64: { $($fn_name:ident),* }) => {
        $(
            fn $fn_name(self, v: f64) -> Complex<f64> {
                unsafe { bindings::$fn_name(v, Complex::new(self, 0.0).into()) }.into()
            }
        )*
    };

    (hankel for Complex<f64>: { $($fn_name:ident),* }) => {
        $(
            fn $fn_name(self, v: f64) -> Complex<f64> {
                unsafe { bindings::$fn_name(v, self.into()) }.into()
            }
        )*
    };
}

// Macro to define BesselArg trait, implement it for both types, and generate public functions
macro_rules! impl_bessel_args {
    (
        bessel: [$(($bessel_fn:ident, $bessel_doc:expr)),* $(,)?],
        hankel: [$(($hankel_fn:ident, $hankel_doc:expr)),* $(,)?] $(,)?
    ) => {
        /// Trait for types that can be used with Bessel functions
        pub trait BesselArg: sealed::Sealed {
            type Output;
            $(fn $bessel_fn(self, v: f64) -> Self::Output;)*
            $(fn $hankel_fn(self, v: f64) -> Complex<f64>;)*
        }

        impl BesselArg for f64 {
            type Output = f64;
            impl_methods!(bessel for f64: { $($bessel_fn),* });
            impl_methods!(hankel for f64: { $($hankel_fn),* });
        }

        impl BesselArg for Complex<f64> {
            type Output = Complex<f64>;
            impl_methods!(bessel for Complex<f64>: { $($bessel_fn),* });
            impl_methods!(hankel for Complex<f64>: { $($hankel_fn),* });
        }

        // Generate public Bessel functions
        $(
            #[doc = $bessel_doc]
            pub fn $bessel_fn<T: BesselArg>(v: f64, x: T) -> T::Output {
                x.$bessel_fn(v)
            }
        )*

        // Generate public Hankel functions
        $(
            #[doc = $hankel_doc]
            pub fn $hankel_fn<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
                z.$hankel_fn(v)
            }
        )*
    };
}

xsf_impl!(cyl_bessel_j0, (x: f64), "Bessel function, 1st kind, order 0");
xsf_impl!(cyl_bessel_j1, (x: f64), "Bessel function, 1st kind, order 1");
xsf_impl!(cyl_bessel_y0, (x: f64), "Bessel function, 2nd kind, order 0");
xsf_impl!(cyl_bessel_y1, (x: f64), "Bessel function, 2nd kind, order 1");
xsf_impl!(cyl_bessel_i0, (x: f64), "Modified Bessel function, 1st kind, order 0");
xsf_impl!(
    cyl_bessel_i0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 0"
);
xsf_impl!(cyl_bessel_i1, (x: f64), "Modified Bessel function, 1st kind, order 1");
xsf_impl!(
    cyl_bessel_i1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 1"
);
xsf_impl!(cyl_bessel_k0, (x: f64), "Modified Bessel function, 2nd kind, order 0");
xsf_impl!(
    cyl_bessel_k0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 0"
);
xsf_impl!(cyl_bessel_k1, (x: f64), "Modified Bessel function, 2nd kind, order 1");
xsf_impl!(
    cyl_bessel_k1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 1"
);

impl_bessel_args! {
    bessel: [
        (cyl_bessel_j, "Bessel function, 1st kind"),
        (cyl_bessel_je, "Exponentially scaled Bessel function, 1st kind"),
        (cyl_bessel_y, "Bessel function, 2nd kind"),
        (cyl_bessel_ye, "Exponentially scaled Bessel function, 2nd kind"),
        (cyl_bessel_i, "Modified Bessel function, 1st kind"),
        (cyl_bessel_ie, "Exponentially scaled modified Bessel function, 1st kind"),
        (cyl_bessel_k, "Modified Bessel function, 2nd kind"),
        (cyl_bessel_ke, "Exponentially scaled modified Bessel function, 2nd kind"),
    ],
    hankel: [
        (cyl_hankel_1, "Hankel function, 1st kind"),
        (cyl_hankel_1e, "Exponentially scaled Hankel function, 1st kind"),
        (cyl_hankel_2, "Hankel function, 2nd kind"),
        (cyl_hankel_2e, "Exponentially scaled Hankel function, 2nd kind"),
    ],
}

xsf_impl!(
    besselpoly,
    (a: f64, lambda: f64, nu: f64),
    "Weighted integral of the Bessel function of the first kind"
);

/// Integrals of Bessel functions of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ J₀(t) dt
/// - ∫₀ˣ Y₀(t) dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `j0int`: Integral for J₀
/// - `y0int`: Integral for Y₀
pub fn it1j0y0(x: f64) -> (f64, f64) {
    let mut j0int = f64::NAN;
    let mut y0int = f64::NAN;
    unsafe {
        bindings::it1j0y0(x, &mut j0int, &mut y0int);
    }
    (j0int, y0int)
}

/// Integrals related to Bessel functions of the first kind of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ (1 - J₀(t))/t dt
/// - ∫ₓ∞ Y₀(t)/t dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `j0int`: Integral for J₀
/// - `y0int`: Integral for Y₀
pub fn it2j0y0(x: f64) -> (f64, f64) {
    let mut j0int = f64::NAN;
    let mut y0int = f64::NAN;
    unsafe {
        bindings::it2j0y0(x, &mut j0int, &mut y0int);
    }
    (j0int, y0int)
}

/// Integrals of modified Bessel functions of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ I₀(t) dt
/// - ∫₀ˣ K₀(t) dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `i0int`: The integral for I₀
/// - `k0int`: The integral for K₀
pub fn it1i0k0(x: f64) -> (f64, f64) {
    let mut i0int = f64::NAN;
    let mut k0int = f64::NAN;
    unsafe {
        bindings::it1i0k0(x, &mut i0int, &mut k0int);
    }
    (i0int, k0int)
}

/// Integrals related to modified Bessel functions of order 0
///
/// Computes the integrals:
///
/// - ∫₀ˣ (I₀(t) - 1)/t dt
/// - ∫ₓ∞ K₀(t)/t dt
///
/// # Arguments
/// - `x` - The input value
///
/// # Returns
/// - `i0int`: The integral for I₀
/// - `k0int`: The integral for K₀
pub fn it2i0k0(x: f64) -> (f64, f64) {
    let mut i0int = f64::NAN;
    let mut k0int = f64::NAN;
    unsafe {
        bindings::it2i0k0(x, &mut i0int, &mut k0int);
    }
    (i0int, k0int)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::{Complex, c64};

    // bessel j

    #[test]
    fn test_cyl_bessel_j0() {
        xsref::test::<f64, _>("cyl_bessel_j0", "d-d", |x: &[f64]| cyl_bessel_j0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_j1() {
        xsref::test::<f64, _>("cyl_bessel_j1", "d-d", |x: &[f64]| cyl_bessel_j1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_j_f64() {
        xsref::test::<f64, _>("cyl_bessel_j", "d_d-d", |x: &[f64]| {
            cyl_bessel_j(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_j_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_j", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_j(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_je_f64() {
        xsref::test::<f64, _>("cyl_bessel_je", "d_d-d", |x: &[f64]| {
            cyl_bessel_je(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_je_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_je", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_je(x[0], c64(x[1], x[2]))
        });
    }

    // bessel y

    #[test]
    fn test_cyl_bessel_y0() {
        xsref::test::<f64, _>("cyl_bessel_y0", "d-d", |x: &[f64]| cyl_bessel_y0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_y1() {
        xsref::test::<f64, _>("cyl_bessel_y1", "d-d", |x: &[f64]| cyl_bessel_y1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_y_f64() {
        xsref::test::<f64, _>("cyl_bessel_y", "d_d-d", |x: &[f64]| {
            cyl_bessel_y(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_y_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_y", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_y(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_ye_f64() {
        xsref::test::<f64, _>("cyl_bessel_ye", "d_d-d", |x: &[f64]| {
            cyl_bessel_ye(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_ye_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_ye", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_ye(x[0], c64(x[1], x[2]))
        });
    }

    // bessel i

    #[test]
    fn test_cyl_bessel_i0() {
        xsref::test::<f64, _>("cyl_bessel_i0", "d-d", |x: &[f64]| cyl_bessel_i0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i0e() {
        xsref::test::<f64, _>("cyl_bessel_i0e", "d-d", |x: &[f64]| cyl_bessel_i0e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i1() {
        xsref::test::<f64, _>("cyl_bessel_i1", "d-d", |x: &[f64]| cyl_bessel_i1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i1e() {
        xsref::test::<f64, _>("cyl_bessel_i1e", "d-d", |x: &[f64]| cyl_bessel_i1e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_i_f64() {
        xsref::test::<f64, _>("cyl_bessel_i", "d_d-d", |x: &[f64]| {
            cyl_bessel_i(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_i_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_i", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_i(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_ie_f64() {
        xsref::test::<f64, _>("cyl_bessel_ie", "d_d-d", |x: &[f64]| {
            cyl_bessel_ie(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_ie_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_ie", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_ie(x[0], c64(x[1], x[2]))
        });
    }

    // bessel k

    #[test]
    fn test_cyl_bessel_k0() {
        xsref::test::<f64, _>("cyl_bessel_k0", "d-d", |x: &[f64]| cyl_bessel_k0(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k0e() {
        xsref::test::<f64, _>("cyl_bessel_k0e", "d-d", |x: &[f64]| cyl_bessel_k0e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k1() {
        xsref::test::<f64, _>("cyl_bessel_k1", "d-d", |x: &[f64]| cyl_bessel_k1(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k1e() {
        xsref::test::<f64, _>("cyl_bessel_k1e", "d-d", |x: &[f64]| cyl_bessel_k1e(x[0]));
    }

    #[test]
    fn test_cyl_bessel_k_f64() {
        xsref::test::<f64, _>("cyl_bessel_k", "d_d-d", |x: &[f64]| {
            cyl_bessel_k(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_k_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_k", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_k(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_bessel_ke_f64() {
        xsref::test::<f64, _>("cyl_bessel_ke", "d_d-d", |x: &[f64]| {
            cyl_bessel_ke(x[0], x[1])
        });
    }

    #[test]
    fn test_cyl_bessel_ke_c64() {
        xsref::test::<Complex<f64>, _>("cyl_bessel_ke", "d_cd-cd", |x: &[f64]| {
            cyl_bessel_ke(x[0], c64(x[1], x[2]))
        });
    }

    // hankel

    #[test]
    fn test_cyl_hankel_1_c64() {
        xsref::test::<Complex<f64>, _>("cyl_hankel_1", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_1(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_hankel_1e_c64() {
        xsref::test::<Complex<f64>, _>("cyl_hankel_1e", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_1e(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_hankel_2_c64() {
        xsref::test::<Complex<f64>, _>("cyl_hankel_2", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_2(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_cyl_hankel_2e_c64() {
        xsref::test::<Complex<f64>, _>("cyl_hankel_2e", "d_cd-cd", |x: &[f64]| {
            cyl_hankel_2e(x[0], c64(x[1], x[2]))
        });
    }

    // besselpoly

    #[test]
    fn test_besselpoly_f64() {
        xsref::test::<f64, _>("besselpoly", "d_d_d-d", |x: &[f64]| {
            besselpoly(x[0], x[1], x[2])
        });
    }

    // bessel integrals

    #[test]
    fn test_it1j0y0_f64() {
        xsref::test::<(f64, f64), _>("it1j0y0", "d-d_d", |x: &[f64]| it1j0y0(x[0]));
    }

    #[test]
    fn test_it2j0y0_f64() {
        xsref::test::<(f64, f64), _>("it2j0y0", "d-d_d", |x: &[f64]| it2j0y0(x[0]));
    }

    #[test]
    fn test_it1i0k0_f64() {
        xsref::test::<(f64, f64), _>("it1i0k0", "d-d_d", |x: &[f64]| it1i0k0(x[0]));
    }

    #[test]
    fn test_it2i0k0_f64() {
        xsref::test::<(f64, f64), _>("it2i0k0", "d-d_d", |x: &[f64]| it2i0k0(x[0]));
    }
}
