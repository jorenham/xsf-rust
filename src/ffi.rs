#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(unreachable_pub)]

/// bindgen bindings to wrapped xsf C++ functions
pub(crate) mod xsf {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
    pub use root::xsf_wrapper::*;
}

// C++ std::complex type wrapper

pub(crate) use xsf::cdouble;

impl From<f64> for cdouble {
    #[inline(always)]
    fn from(x: f64) -> Self {
        if x.is_nan() {
            Self { re: x, im: x }
        } else {
            Self { re: x, im: 0.0 }
        }
    }
}

impl From<num_complex::Complex<f64>> for cdouble {
    #[inline(always)]
    fn from(z: num_complex::Complex<f64>) -> Self {
        Self { re: z.re, im: z.im }
    }
}

impl From<cdouble> for num_complex::Complex<f64> {
    #[inline(always)]
    fn from(z: cdouble) -> Self {
        Self::new(z.re, z.im)
    }
}
