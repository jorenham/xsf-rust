//! Pure-rust implementations of some mathematical convenience functions from
//! [NumPy](https://github.com/numpy/numpy)

mod npymath;
pub(crate) use npymath::LogAddExpArg;
pub use npymath::{logaddexp, logaddexp2};
