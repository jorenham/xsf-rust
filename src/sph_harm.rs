use std::os::raw::c_int;

use crate::bindings;
use num_complex::Complex;

/// Spherical harmonics
pub fn sph_harm_y(n: i32, m: i32, theta: f64, phi: f64) -> Complex<f64> {
    unsafe { bindings::sph_harm_y(n as c_int, m as c_int, theta, phi) }.into()
}

#[cfg(test)]
mod tests {
    // TODO: no xsref table -> need manual smoketests
}
