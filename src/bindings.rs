#![allow(non_snake_case)]
// needed for std::complex_value_type
#![allow(non_camel_case_types)]
#![allow(dead_code)]

use num_complex::Complex;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub(crate) use root::xsf_wrapper::*;

impl root::std::complex<f64> {
    pub(crate) fn new(re: f64, im: f64) -> Self {
        unsafe { complex__new(re, im) }
    }
}

impl From<Complex<f64>> for root::std::complex<f64> {
    fn from(z: Complex<f64>) -> Self {
        Self::new(z.re, z.im)
    }
}

impl From<root::std::complex<f64>> for Complex<f64> {
    fn from(z: root::std::complex<f64>) -> Self {
        let mut re: f64 = 0.0;
        let mut im: f64 = 0.0;
        unsafe {
            complex__values(z, &mut re, &mut im);
        }
        Self::new(re, im)
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
