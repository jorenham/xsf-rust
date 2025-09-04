use crate::bindings;

/// Cube root
pub fn cbrt(x: f64) -> f64 {
    unsafe { bindings::cbrt_(x) }
}
