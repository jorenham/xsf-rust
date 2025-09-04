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
