use crate::bindings;
use crate::bindings::xsf_impl;
use crate::utils::c_complex64_nan;
use num_complex::Complex;

xsf_impl!(ber, (x: f64), "Kelvin function `ber`");
xsf_impl!(bei, (x: f64), "Kelvin function `bei`");
xsf_impl!(ker, (x: f64), "Kelvin function `ker`");
xsf_impl!(kei, (x: f64), "Kelvin function `kei`");

xsf_impl!(berp, (x: f64), "Derivative of the Kelvin function `ber`");
xsf_impl!(beip, (x: f64), "Derivative of the Kelvin function `bei`");
xsf_impl!(kerp, (x: f64), "Derivative of the Kelvin function `ker`");
xsf_impl!(keip, (x: f64), "Derivative of the Kelvin function `kei`");

/// Kelvin functions as complex numbers
///
/// # Arguments
/// - `x` - Real argument
///
/// # Returns
/// - *Be*: Value of the Kelvin function [`ber`] + *i* [`bei`]
/// - *Ke*: Value of the Kelvin function [`ker`] + *i* [`kei`]
/// - *Be*': Derivative of the Kelvin function [`berp`] + *i* [`beip`]
/// - *Ke*': Derivative of the Kelvin function [`kerp`] + *i* [`keip`]
pub fn kelvin(x: f64) -> (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>) {
    let mut be = c_complex64_nan();
    let mut ke = c_complex64_nan();
    let mut bep = c_complex64_nan();
    let mut kep = c_complex64_nan();
    unsafe {
        bindings::kelvin(x, &mut be, &mut ke, &mut bep, &mut kep);
    }
    (be.into(), ke.into(), bep.into(), kep.into())
}
