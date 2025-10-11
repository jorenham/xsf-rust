/// Characteristic value of even Mathieu functions
pub fn cem_cva(m: f64, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::cem_cva(m, q) }
}

/// Characteristic value of odd Mathieu functions
pub fn sem_cva(m: f64, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::sem_cva(m, q) }
}

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
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::cem(m, q, x, &mut y, &mut yp);
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
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::sem(m, q, x, &mut y, &mut yp);
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
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::mcm1(m, q, x, &mut y, &mut yp);
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
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm1(m, q, x, &mut y, &mut yp);
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
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::mcm2(m, q, x, &mut y, &mut yp);
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
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm2(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_cem_cva() {
        xsref::test::<f64, _>("cem_cva", "d_d-d", |x: &[f64]| cem_cva(x[0], x[1]));
    }

    #[test]
    fn test_sem_cva() {
        xsref::test::<f64, _>("sem_cva", "d_d-d", |x: &[f64]| sem_cva(x[0], x[1]));
    }

    #[test]
    fn test_cem() {
        xsref::test::<(f64, f64), _>("cem", "d_d_d-d_d", |x: &[f64]| cem(x[0], x[1], x[2]));
    }

    #[test]
    fn test_sem() {
        xsref::test::<(f64, f64), _>("sem", "d_d_d-d_d", |x: &[f64]| sem(x[0], x[1], x[2]));
    }

    #[test]
    fn test_mcm1() {
        xsref::test::<(f64, f64), _>("mcm1", "d_d_d-d_d", |x: &[f64]| mcm1(x[0], x[1], x[2]));
    }

    #[test]
    fn test_msm1() {
        xsref::test::<(f64, f64), _>("msm1", "d_d_d-d_d", |x: &[f64]| msm1(x[0], x[1], x[2]));
    }

    #[test]
    fn test_mcm2() {
        xsref::test::<(f64, f64), _>("mcm2", "d_d_d-d_d", |x: &[f64]| mcm2(x[0], x[1], x[2]));
    }

    #[test]
    fn test_msm2() {
        xsref::test::<(f64, f64), _>("msm2", "d_d_d-d_d", |x: &[f64]| msm2(x[0], x[1], x[2]));
    }
}
