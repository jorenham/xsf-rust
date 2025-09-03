use std::os::raw::c_long;

use crate::bindings;
use num_complex::Complex;

/// Lambert W function.
///
/// The Lambert W function `W(z)` is defined as the inverse function of `w * exp(w)`. In other
/// words, the value of `W(z)` is such that `z = W(z) * exp(W(z))` for any complex number `z`.
///
/// The Lambert W function is a multivalued function with infinitely many branches. Each branch
/// gives a separate solution of the equation `z = w exp(w)`. Here, the branches are indexed by the
/// integer `k`.
pub fn lambertw(z: Complex<f64>, k: c_long, tol: f64) -> Complex<f64> {
    unsafe { bindings::lambertw(z.into(), k, tol) }.into()
}
