#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub(crate) use root::xsf_wrapper::*;

// C++ std::complex type wrapper

impl cdouble {
    pub(crate) fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
}

impl Default for cdouble {
    fn default() -> Self {
        Self::new(0.0, 0.0)
    }
}

impl From<f64> for cdouble {
    fn from(x: f64) -> Self {
        Self::new(x, 0.0)
    }
}

impl From<num_complex::Complex<f64>> for cdouble {
    fn from(z: num_complex::Complex<f64>) -> Self {
        Self::new(z.re, z.im)
    }
}

impl From<cdouble> for num_complex::Complex<f64> {
    fn from(z: cdouble) -> Self {
        Self::new(z.re, z.im)
    }
}

pub(crate) const CNAN: cdouble = cdouble {
    re: f64::NAN,
    im: f64::NAN,
};

// macros

macro_rules! xsf_impl {
    ($name:ident, ($($param:ident: $type:ty),*), $docs:expr) => {
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { crate::ffi::$name($($param),*) }
        }
    };
}

pub(crate) use xsf_impl;
