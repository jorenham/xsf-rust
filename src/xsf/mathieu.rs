use core::ffi::c_int;

use num_traits::ToPrimitive;

/// Characteristic value of even Mathieu functions
///
/// # See also
/// - [`mathieu_b`]: Characteristic value of odd Mathieu functions
/// - [`mathieu_cem`]: Even Mathieu function
#[doc(alias = "cem_cva")]
#[must_use]
#[inline]
pub fn mathieu_a(m: u32, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::cem_cva(m.into(), q) }
}

/// Characteristic value of odd Mathieu functions
///
/// # See also
/// - [`mathieu_a`]: Characteristic value of even Mathieu functions
/// - [`mathieu_sem`]: Odd Mathieu function
#[doc(alias = "sem_cva")]
#[must_use]
#[inline]
pub fn mathieu_b(m: u32, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::sem_cva(m.into(), q) }
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
#[must_use]
#[inline]
pub fn mathieu_cem(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::cem(m.into(), q, x, &raw mut y, &raw mut yp);
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
#[must_use]
#[inline]
pub fn mathieu_sem(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::sem(m.into(), q, x, &raw mut y, &raw mut yp);
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
#[must_use]
#[inline]
pub fn mathieu_modcem1(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::mcm1(m.into(), q, x, &raw mut y, &raw mut yp);
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
/// - [`mathieu_sem`]: Odd Mathieu function
/// - [`mathieu_modsem2`]: Odd modified Mathieu function of the second kind
/// - [`mathieu_modcem1`]: Even modified Mathieu function of the first kind
#[doc(alias = "msm1")]
#[must_use]
#[inline]
pub fn mathieu_modsem1(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm1(m.into(), q, x, &raw mut y, &raw mut yp);
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
#[must_use]
#[inline]
pub fn mathieu_modcem2(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::mcm2(m.into(), q, x, &raw mut y, &raw mut yp);
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
#[must_use]
#[inline]
pub fn mathieu_modsem2(m: u32, q: f64, x: f64) -> (f64, f64) {
    let (mut y, mut yp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::msm2(m.into(), q, x, &raw mut y, &raw mut yp);
    }
    (y, yp)
}

#[inline]
fn fcoef<const KD: c_int>(m: u32, q: f64) -> Result<([f64; 251], usize), String> {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L1582-L1639

    // if (q < 0):
    if q < 0.0 {
        // raise ValueError("q >=0")
        return Err("q >=0".to_string());
    }

    let q_sqrt = q.sqrt();
    let qm = if q <= 1.0 {
        7.5 + 56.1 * q_sqrt - 134.7 * q + 90.7 * q_sqrt * q
    } else {
        17.0 + 3.1 * q_sqrt - 0.126 * q + 0.0037 * q_sqrt * q
    };
    // km = int(qm + 0.5*m)
    let km = (qm + 0.5 * f64::from(m)).floor().to_usize().unwrap();
    // if km > 251:
    if km > 251 {
        // warnings.warn("Too many predicted coefficients.", RuntimeWarning, stacklevel=2)
        return Err(format!("Too many predicted coefficients ({km})"));
    }

    let cv = if KD == 1 || KD == 2 {
        mathieu_a(m, q)
    } else if KD == 3 || KD == 4 {
        mathieu_b(m, q)
    } else {
        unreachable!()
    };

    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_specfun.pyx#L158-L171
    let mut fc = [f64::NAN; 251];
    unsafe {
        crate::ffi::xsf::fcoef(KD, m.try_into().unwrap(), q, cv, fc.as_mut_ptr());
    }
    if fc[0].is_nan() {
        assert!(fc.iter().all(|c| c.is_nan()), "expected all NaN");
        Err("Computation of Fourier coefficients failed".to_string())
    } else {
        Ok((fc, km))
    }
}

/// Fourier coefficients for even Mathieu and modified Mathieu functions
///
/// See the SciPy documentation for `scipy.special.mathieu_even_coef` for details:
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_even_coef.html>
///
/// # Errors
/// - `q < 0`
/// - too many predicted coefficients
///
/// # See also
/// - [`mathieu_odd_coef`]: Fourier coefficients for odd Mathieu and modified Mathieu functions
#[inline]
pub fn mathieu_even_coef(m: u32, q: f64) -> Result<Vec<f64>, String> {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L1582-L1639

    // kd = 1
    // m = int(floor(m))
    // if m % 2:
    //     kd = 2
    // fc = _specfun.fcoef(kd, m, q, a)
    let (fc, km) = if m % 2 == 0 {
        fcoef::<1>(m, q)?
    } else {
        fcoef::<2>(m, q)?
    };

    // return fc[:km]
    Ok(fc[..km].to_vec())
}

/// Fourier coefficients for odd Mathieu and modified Mathieu functions
///
/// See the SciPy documentation for `scipy.special.mathieu_odd_coef` for details:
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_odd_coef.html>
///
/// # Errors
/// - `q < 0`
/// - too many predicted coefficients
///
/// # See also
/// - [`mathieu_even_coef`]: Fourier coefficients for even Mathieu and modified Mathieu functions
#[inline]
pub fn mathieu_odd_coef(m: u32, q: f64) -> Result<Vec<f64>, String> {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L1642-L1698

    // kd = 4
    // m = int(floor(m))
    // if m % 2:
    //     kd = 3
    // fc = _specfun.fcoef(kd, m, q, b)
    let (fc, km) = if m % 2 == 0 {
        fcoef::<4>(m, q)?
    } else {
        fcoef::<3>(m, q)?
    };

    // return fc[:km]
    Ok(fc[..km].to_vec())
}

#[cfg(test)]
mod tests {
    use num_traits::ToPrimitive;

    #[test]
    fn test_mathieu_a() {
        xsref::test("cem_cva", "d_d-d", |x| {
            crate::mathieu_a(x[0].to_u32().unwrap(), x[1])
        });
    }

    #[test]
    fn test_mathieu_b() {
        xsref::test("sem_cva", "d_d-d", |x| {
            crate::mathieu_b(x[0].to_u32().unwrap(), x[1])
        });
    }

    #[test]
    fn test_mathieu_cem() {
        xsref::test("cem", "d_d_d-d_d", |x| {
            crate::mathieu_cem(x[0].to_u32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_sem() {
        xsref::test("sem", "d_d_d-d_d", |x| {
            crate::mathieu_sem(x[0].to_u32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modcem1() {
        xsref::test("mcm1", "d_d_d-d_d", |x| {
            crate::mathieu_modcem1(x[0].to_u32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modsem1() {
        xsref::test("msm1", "d_d_d-d_d", |x| {
            crate::mathieu_modsem1(x[0].to_u32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modcem2() {
        xsref::test("mcm2", "d_d_d-d_d", |x| {
            crate::mathieu_modcem2(x[0].to_u32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_modsem2() {
        xsref::test("msm2", "d_d_d-d_d", |x| {
            crate::mathieu_modsem2(x[0].to_u32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_mathieu_even_coef() {
        // Smoke test: verify function returns Ok and non-empty coefficients
        let result = crate::mathieu_even_coef(2, 1.0);
        assert!(result.is_ok());
        let coefs = result.unwrap();
        assert!(!coefs.is_empty());
        assert!(coefs.iter().all(|&c| c.is_finite()));

        // Test with q = 0
        let result = crate::mathieu_even_coef(0, 0.0);
        assert!(result.is_ok());

        // Test error case: negative q
        let result = crate::mathieu_even_coef(1, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_mathieu_odd_coef() {
        // Smoke test: verify function returns Ok and non-empty coefficients
        let result = crate::mathieu_odd_coef(1, 1.0);
        assert!(result.is_ok());
        let coefs = result.unwrap();
        assert!(!coefs.is_empty());
        assert!(coefs.iter().all(|&c| c.is_finite()));

        // Test with q = 0
        let result = crate::mathieu_odd_coef(1, 0.0);
        assert!(result.is_ok());

        // Test error case: negative q
        let result = crate::mathieu_odd_coef(1, -1.0);
        assert!(result.is_err());
    }
}
