#![allow(non_snake_case)]
// needed for std::complex_value_type
#![allow(non_camel_case_types)]
#![allow(dead_code)]

use num_complex::Complex;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

impl std_complex<f64> {
    pub(crate) fn new(re: f64, im: f64) -> Self {
        unsafe { complex__new(re, im) }
    }
    pub(crate) fn real(self) -> f64 {
        unsafe { complex__re(self) }
    }
    pub(crate) fn imag(self) -> f64 {
        unsafe { complex__im(self) }
    }
}

impl From<std_complex<f64>> for Complex<f64> {
    fn from(z: std_complex<f64>) -> Self {
        Complex::new(z.real(), z.imag())
    }
}

impl From<Complex<f64>> for std_complex<f64> {
    fn from(z: Complex<f64>) -> Self {
        Self::new(z.re, z.im)
    }
}

macro_rules! xsf_impl {
    ($name:ident, ($($param:ident: $type:ty),*), $docs:expr) => {
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { crate::bindings::$name($($param),*) }
        }
    };
}

pub(crate) use xsf_impl;
