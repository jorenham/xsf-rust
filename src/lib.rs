#![cfg_attr(not(test), no_std)]
#![warn(
    missing_debug_implementations,
    // missing_docs, // TODO
    rust_2018_idioms,
    unreachable_pub
)]

extern crate alloc;

#[cfg(test)]
mod testing;

mod ffi;
mod utils;

mod xsf;
pub use xsf::*;

pub use xsf::cephes::betainc;
pub use xsf::cephes::betaincinv;
pub use xsf::cephes::expn;
pub use xsf::cephes::lanczos_sum_expg_scaled;
pub use xsf::cephes::lgam1p;
pub use xsf::cephes::round;
pub use xsf::cephes::spence;
pub use xsf::cephes::{erfcinv, erfinv};
pub use xsf::cephes::{pow_falling, pow_rising};

mod scipy_special;
pub use scipy_special::*;
