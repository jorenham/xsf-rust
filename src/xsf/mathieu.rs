/// Characteristic value of even Mathieu functions
#[doc(alias = "mathieu_a")]
pub fn cem_cva(m: f64, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::cem_cva(m, q) }
}

/// Characteristic value of odd Mathieu functions
#[doc(alias = "mathieu_b")]
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
#[doc(alias = "mathieu_cem")]
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
#[doc(alias = "mathieu_sem")]
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
#[doc(alias = "modcem1")]
#[doc(alias = "mathieu_modcem1")]
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
#[doc(alias = "modsem1")]
#[doc(alias = "mathieu_modsem1")]
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
#[doc(alias = "modcem2")]
#[doc(alias = "mathieu_modcem2")]
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
#[doc(alias = "modsem2")]
#[doc(alias = "mathieu_modsem2")]
pub fn msm2(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm2(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_cem_cva() {
        crate::xsref::test("cem_cva", "d_d-d", |x| crate::cem_cva(x[0], x[1]));
    }

    #[test]
    fn test_sem_cva() {
        crate::xsref::test("sem_cva", "d_d-d", |x| crate::sem_cva(x[0], x[1]));
    }

    #[test]
    fn test_cem() {
        crate::xsref::test("cem", "d_d_d-d_d", |x| crate::cem(x[0], x[1], x[2]));
    }

    #[test]
    fn test_sem() {
        crate::xsref::test("sem", "d_d_d-d_d", |x| crate::sem(x[0], x[1], x[2]));
    }

    #[test]
    fn test_mcm1() {
        crate::xsref::test("mcm1", "d_d_d-d_d", |x| crate::mcm1(x[0], x[1], x[2]));
    }

    #[test]
    fn test_msm1() {
        crate::xsref::test("msm1", "d_d_d-d_d", |x| crate::msm1(x[0], x[1], x[2]));
    }

    #[test]
    fn test_mcm2() {
        crate::xsref::test("mcm2", "d_d_d-d_d", |x| crate::mcm2(x[0], x[1], x[2]));
    }

    #[test]
    fn test_msm2() {
        crate::xsref::test("msm2", "d_d_d-d_d", |x| crate::msm2(x[0], x[1], x[2]));
    }
}
