//! Pure-rust implementations of some mathematical convenience functions from
//! [NumPy](https://github.com/numpy/numpy)

mod function_base;
mod npymath;
pub use function_base::sinc;
pub(crate) use npymath::LogAddExpArg;
pub use npymath::{logaddexp, logaddexp2};
