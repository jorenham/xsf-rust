#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub(crate) use root::std::complex;
pub(crate) use root::xsf_wrapper::*;

pub(crate) type cdouble = complex<f64>;

impl cdouble {
    pub(crate) fn new(re: f64, im: f64) -> Self {
        unsafe { complex__new(re, im) }
    }
}

impl From<num_complex::Complex<f64>> for cdouble {
    fn from(z: num_complex::Complex<f64>) -> Self {
        Self::new(z.re, z.im)
    }
}

impl From<cdouble> for num_complex::Complex<f64> {
    fn from(z: cdouble) -> Self {
        let mut re: f64 = 0.0;
        let mut im: f64 = 0.0;
        unsafe {
            complex__values(z, &mut re, &mut im);
        }
        Self::new(re, im)
    }
}

#[inline(always)]
pub(crate) fn complex_nan() -> cdouble {
    complex::new(f64::NAN, f64::NAN)
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
