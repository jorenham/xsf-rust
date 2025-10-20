/// Characteristic value of even Mathieu functions
///
/// # See also
/// - [`mathieu_b`]: Characteristic value of odd Mathieu functions
/// - [`mathieu_cem`]: Even Mathieu function
#[doc(alias = "cem_cva")]
pub fn mathieu_a(m: u32, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::cem_cva(m as f64, q) }
}

/// Characteristic value of odd Mathieu functions
///
/// # See also
/// - [`mathieu_a`]: Characteristic value of even Mathieu functions
/// - [`mathieu_sem`]: Odd Mathieu function
#[doc(alias = "sem_cva")]
pub fn mathieu_b(m: u32, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::sem_cva(m as f64, q) }
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
/// - `yp`: Value of the derivative w.r.t. *x*
///
/// # See also
/// - [`mathieu_sem`]: Odd Mathieu function
/// - [`mathieu_a`]: Characteristic value of even Mathieu functions
#[doc(alias = "cem")]
pub fn mathieu_cem(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::cem(m as f64, q, x, &mut y, &mut yp);
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
/// - `yp`: Value of the derivative w.r.t. *x*
///
/// # See also
/// - [`mathieu_cem`]: Even Mathieu function
/// - [`mathieu_b`]: Characteristic value of odd Mathieu functions
#[doc(alias = "sem")]
pub fn mathieu_sem(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::sem(m as f64, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Even modified Mathieu function of the first kind and its derivative
///
/// Evaluates the even modified Mathieu function of the first kind, *Mc1<sub>m</sub>(x, q)*,
/// and its derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. *x*
///
/// # See also
/// - [`mathieu_cem`]: Even Mathieu function
/// - [`mathieu_modcem2`]: Even modified Mathieu function of the second kind
/// - [`mathieu_modsem1`]: Odd modified Mathieu function of the first kind
#[doc(alias = "mcm1")]
pub fn mathieu_modcem1(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::mcm1(m as f64, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Odd modified Mathieu function of the first kind and its derivative
///
/// Evaluates the odd modified Mathieu function of the first kind, *Ms1<sub>m</sub>(x, q)*, and its
/// derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. *x*
///
/// # See also
/// - [`mathieu_sem`]: Even Mathieu function
/// - [`mathieu_modsem2`]: Even modified Mathieu function of the second kind
/// - [`mathieu_modcem1`]: Odd modified Mathieu function of the first kind
#[doc(alias = "msm1")]
pub fn mathieu_modsem1(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm1(m as f64, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Even modified Mathieu function of the second kind and its derivative
///
/// Evaluates the even modified Mathieu function of the second kind, *Mc2<sub>m</sub>(x, q)*, and
/// its derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. *x*
///
/// # See also
/// - [`mathieu_cem`]: Even Mathieu function
/// - [`mathieu_modcem1`]: Even modified Mathieu function of the first kind
/// - [`mathieu_modsem2`]: Odd modified Mathieu function of the second kind
#[doc(alias = "mcm2")]
pub fn mathieu_modcem2(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::mcm2(m as f64, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

/// Odd modified Mathieu function of the second kind and its derivative
///
/// Evaluates the odd modified Mathieu function of the second kind, *Ms2<sub>m</sub>(x, q)*, and
/// its derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// # Arguments
/// - `m`: The order of the function
/// - `q`: The parameter of the function.
/// - `x`: Argument of the function, *given in degrees, not radians*
///
/// # Returns
/// - `y`: value of the function
/// - `yp`: Value of the derivative w.r.t. *x*
///
/// # See also
/// - [`mathieu_sem`]: Odd Mathieu function
/// - [`mathieu_modsem1`]: Odd modified Mathieu function of the first kind
/// - [`mathieu_modcem2`]: Even modified Mathieu function of the second kind
#[doc(alias = "msm2")]
pub fn mathieu_modsem2(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm2(m as f64, q, x, &mut y, &mut yp);
    }
    (y, yp)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_mathieu_a() {
        crate::xsref::test("cem_cva", "d_d-d", |x| crate::mathieu_a(x[0] as u32, x[1]));
    }

    #[test]
    fn test_mathieu_b() {
        crate::xsref::test("sem_cva", "d_d-d", |x| crate::mathieu_b(x[0] as u32, x[1]));
    }

    #[test]
    fn test_mathieu_cem() {
        crate::xsref::test("cem", "d_d_d-d_d", |x| {
            crate::mathieu_cem(x[0] as u32, x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_sem() {
        crate::xsref::test("sem", "d_d_d-d_d", |x| {
            crate::mathieu_sem(x[0] as u32, x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modcem1() {
        crate::xsref::test("mcm1", "d_d_d-d_d", |x| {
            crate::mathieu_modcem1(x[0] as u32, x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modsem1() {
        crate::xsref::test("msm1", "d_d_d-d_d", |x| {
            crate::mathieu_modsem1(x[0] as u32, x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modcem2() {
        crate::xsref::test("mcm2", "d_d_d-d_d", |x| {
            crate::mathieu_modcem2(x[0] as u32, x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modsem2() {
        crate::xsref::test("msm2", "d_d_d-d_d", |x| {
            crate::mathieu_modsem2(x[0] as u32, x[1], x[2])
        });
    }
}
