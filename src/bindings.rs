#![allow(non_snake_case)]

use num_complex::Complex;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

macro_rules! xsf_impl {
    ($name:ident, ($($param:ident: $type:ty),*), $docs:expr) => {
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { crate::bindings::$name($($param),*) }
        }
    };
}

pub(crate) use xsf_impl;

impl<T> std_complex<T> {
    pub(crate) fn new(re: T, im: T) -> Self {
        Self {
            _phantom_0: std::marker::PhantomData,
            _M_real: re,
            _M_imag: im,
        }
    }
}

impl<T> From<std_complex<T>> for Complex<T> {
    fn from(val: std_complex<T>) -> Self {
        Complex::new(val._M_real, val._M_imag)
    }
}

impl<T> From<Complex<T>> for std_complex<T> {
    fn from(c: Complex<T>) -> Self {
        Self::new(c.re, c.im)
    }
}
