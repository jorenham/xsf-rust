use crate::bindings;
use crate::bindings::xsf_impl;

xsf_impl!(cem_cva, (m: f64, q: f64), "Characteristic value of even Mathieu functions");
xsf_impl!(sem_cva, (m: f64, q: f64), "Characteristic value of odd Mathieu functions");

/// Even Mathieu function and its derivative
///
/// Returns the even Mathieu function, *ce_m(x, q)*, of order `m` and parameter `q` evaluated at `x`
/// (given in degrees). Also returns the derivative with respect to `x` of *ce_m(x, q)*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn cem(m: f64, q: f64, x: f64) -> (f64, f64) {
    let mut y = f64::NAN;
    let mut yp = f64::NAN;
    unsafe {
        bindings::cem(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Odd Mathieu function and its derivative
///
/// Returns the odd Mathieu function, *se_m(x, q)*, of order `m` and parameter `q` evaluated at `x`
/// (given in degrees). Also returns the derivative with respect to `x` of *se_m(x, q)*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn sem(m: f64, q: f64, x: f64) -> (f64, f64) {
    let mut y = f64::NAN;
    let mut yp = f64::NAN;
    unsafe {
        bindings::sem(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Even modified Mathieu function of the first kind and its derivative
///
/// Evaluates the even modified Mathieu function of the first kind, `Mc1m(x, q)`, and its
/// derivative at `x` for order `m` and parameter `q`.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn mcm1(m: f64, q: f64, x: f64) -> (f64, f64) {
    let mut y = f64::NAN;
    let mut yp = f64::NAN;
    unsafe {
        bindings::mcm1(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Odd modified Mathieu function of the first kind and its derivative
///
/// Evaluates the odd modified Mathieu function of the first kind, `Ms1m(x, q)`, and its
/// derivative at `x` for order `m` and parameter `q`.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn msm1(m: f64, q: f64, x: f64) -> (f64, f64) {
    let mut y = f64::NAN;
    let mut yp = f64::NAN;
    unsafe {
        bindings::msm1(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Even modified Mathieu function of the second kind and its derivative
///
/// Evaluates the even modified Mathieu function of the second kind, `Mc2m(x, q)`, and its
/// derivative at `x` for order `m` and parameter `q`.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn mcm2(m: f64, q: f64, x: f64) -> (f64, f64) {
    let mut y = f64::NAN;
    let mut yp = f64::NAN;
    unsafe {
        bindings::mcm2(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Odd modified Mathieu function of the second kind and its derivative
///
/// Evaluates the odd modified Mathieu function of the second kind, `Ms2m(x, q)`, and its
/// derivative at `x` for order `m` and parameter `q`.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn msm2(m: f64, q: f64, x: f64) -> (f64, f64) {
    let mut y = f64::NAN;
    let mut yp = f64::NAN;
    unsafe {
        bindings::msm2(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}
