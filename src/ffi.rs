#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(unreachable_pub)]

/// bindgen bindings to wrapped xsf C++ functions
pub(crate) mod xsf {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
    pub use root::xsf_wrapper::*;
}

// c_complex is type-aliased to num_complex::Complex<f64> in the generated bindings
// Both have the same memory layout (two consecutive f64 fields: re, im)
// This allows us to use Complex<f64> directly in FFI without any conversions!
