/// Characteristic value of even Mathieu functions
#[doc(alias = "cem_cva")]
pub fn mathieu_a(m: f64, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::cem_cva(m, q) }
}

/// Characteristic value of odd Mathieu functions
#[doc(alias = "sem_cva")]
pub fn mathieu_b(m: f64, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::sem_cva(m, q) }
}

/// Even Mathieu function and its derivative
///
/// Returns the even Mathieu function, *ce<sub>m</sub>(x, q)*, of order *m* and parameter *q*
/// evaluated at *x* (given in degrees).
/// Also returns the derivative with respect to *x* of *ce<sub>m</sub>(x, q)*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn mathieu_cem(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::cem(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Odd Mathieu function and its derivative
///
/// Returns the odd Mathieu function, *se<sub>m</sub>(x, q)*, of order *m* and parameter *q*
/// evaluated at *x* (given in degrees).
/// Also returns the derivative with respect to *x* of *se<sub>m</sub>(x, q)*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. `x`
pub fn mathieu_sem(m: f64, q: f64, x: f64) -> (f64, f64) {
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
#[doc(alias = "mcm1")]
pub fn mathieu_modcem1(m: f64, q: f64, x: f64) -> (f64, f64) {
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
#[doc(alias = "msm1")]
pub fn mathieu_modsem1(m: f64, q: f64, x: f64) -> (f64, f64) {
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
#[doc(alias = "mcm2")]
pub fn mathieu_modcem2(m: f64, q: f64, x: f64) -> (f64, f64) {
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
#[doc(alias = "msm2")]
pub fn mathieu_modsem2(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm2(m, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_mathieu_a() {
        crate::xsref::test("cem_cva", "d_d-d", |x| crate::mathieu_a(x[0], x[1]));
    }

    #[test]
    fn test_mathieu_b() {
        crate::xsref::test("sem_cva", "d_d-d", |x| crate::mathieu_b(x[0], x[1]));
    }

    #[test]
    fn test_mathieu_cem() {
        crate::xsref::test("cem", "d_d_d-d_d", |x| crate::mathieu_cem(x[0], x[1], x[2]));
    }

    #[test]
    fn test_mathieu_sem() {
        crate::xsref::test("sem", "d_d_d-d_d", |x| crate::mathieu_sem(x[0], x[1], x[2]));
    }

    #[test]
    fn test_mathieu_modcem1() {
        crate::xsref::test("mcm1", "d_d_d-d_d", |x| {
            crate::mathieu_modcem1(x[0], x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modsem1() {
        crate::xsref::test("msm1", "d_d_d-d_d", |x| {
            crate::mathieu_modsem1(x[0], x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modcem2() {
        crate::xsref::test("mcm2", "d_d_d-d_d", |x| {
            crate::mathieu_modcem2(x[0], x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modsem2() {
        crate::xsref::test("msm2", "d_d_d-d_d", |x| {
            crate::mathieu_modsem2(x[0], x[1], x[2])
        });
    }
}
