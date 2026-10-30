//! The Mathieu functions are translated into pure Rust from `xsf/mathieu.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the `double` → `int` conversions saturate instead of being undefined behavior.
//! This way, the modified Mathieu functions don't return 0 or NaN if `√q·exp(|x|) > 2.39e9`, and
//! the characteristic values are NaN instead of 0 if `q` is NaN and `m > 12`.

use num_traits::ToPrimitive;

use crate::xsf::specfun::{Status, cva2, fcoef as specfun_fcoef, mtu0, mtu12};

/* Mathieu functions */
/* Characteristic values */
#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn cem_cva(m: f64, q: f64) -> f64 {
    let mut kd = 1;

    if (m < 0.0) || (m != m.floor()) {
        return f64::NAN;
    }
    let int_m = m as i32;
    if q < 0.0 {
        /* https://dlmf.nist.gov/28.2#E26 */
        if int_m % 2 == 0 {
            return cem_cva(m, -q);
        }
        return sem_cva(m, -q);
    }

    if int_m % 2 != 0 {
        kd = 2;
    }
    cva2(kd, int_m, q)
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn sem_cva(m: f64, q: f64) -> f64 {
    let mut kd = 4;

    if (m <= 0.0) || (m != m.floor()) {
        return f64::NAN;
    }
    let int_m = m as i32;
    if q < 0.0 {
        /* https://dlmf.nist.gov/28.2#E26 */
        if int_m % 2 == 0 {
            return sem_cva(m, -q);
        }
        return cem_cva(m, -q);
    }
    if int_m % 2 != 0 {
        kd = 3;
    }
    cva2(kd, int_m, q)
}

/* Mathieu functions */
#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn cem(m: f64, q: f64, x: f64) -> (f64, f64) {
    let kf = 1;
    let (mut csf, mut csd) = (0.0, 0.0);
    if (m < 0.0) || (m != m.floor()) {
        return (f64::NAN, f64::NAN);
    }
    let int_m = m as i32;
    if q < 0.0 {
        /* https://dlmf.nist.gov/28.2#E34 */
        let sgn = f64::from(if (int_m / 2) % 2 == 0 { 1 } else { -1 });
        let (f, d) = if int_m % 2 == 0 {
            cem(m, -q, 90.0 - x)
        } else {
            sem(m, -q, 90.0 - x)
        };
        return (sgn * f, -sgn * d);
    }
    let status = mtu0(kf, int_m, q, x, &mut csf, &mut csd);
    if status != Status::Ok {
        return (f64::NAN, f64::NAN);
    }
    (csf, csd)
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn sem(m: f64, q: f64, x: f64) -> (f64, f64) {
    let kf = 2;
    let (mut csf, mut csd) = (0.0, 0.0);
    if (m < 0.0) || (m != m.floor()) {
        return (f64::NAN, f64::NAN);
    }
    let int_m = m as i32;
    if int_m == 0 {
        return (0.0, 0.0);
    }
    if q < 0.0 {
        /* https://dlmf.nist.gov/28.2#E34 */
        let (sgn, (f, d)) = if int_m % 2 == 0 {
            let sgn = if (int_m / 2) % 2 == 0 { -1 } else { 1 };
            (sgn, sem(m, -q, 90.0 - x))
        } else {
            let sgn = if (int_m / 2) % 2 == 0 { 1 } else { -1 };
            (sgn, cem(m, -q, 90.0 - x))
        };
        let sgn = f64::from(sgn);
        return (sgn * f, -sgn * d);
    }
    let status = mtu0(kf, int_m, q, x, &mut csf, &mut csd);
    if status != Status::Ok {
        return (f64::NAN, f64::NAN);
    }
    (csf, csd)
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn mcm1(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (kf, kc) = (1, 1);
    let (mut f1r, mut d1r, mut f2r, mut d2r) = (0.0, 0.0, 0.0, 0.0);

    if (m < 0.0) || (m != m.floor()) || (q < 0.0) {
        return (f64::NAN, f64::NAN);
    }
    let int_m = m as i32;
    let status = mtu12(kf, kc, int_m, q, x, &mut f1r, &mut d1r, &mut f2r, &mut d2r);
    if status != Status::Ok {
        return (f64::NAN, f64::NAN);
    }
    (f1r, d1r)
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn msm1(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (kf, kc) = (2, 1);
    let (mut f1r, mut d1r, mut f2r, mut d2r) = (0.0, 0.0, 0.0, 0.0);

    if (m < 1.0) || (m != m.floor()) || (q < 0.0) {
        return (f64::NAN, f64::NAN);
    }
    let int_m = m as i32;
    let status = mtu12(kf, kc, int_m, q, x, &mut f1r, &mut d1r, &mut f2r, &mut d2r);
    if status != Status::Ok {
        return (f64::NAN, f64::NAN);
    }
    (f1r, d1r)
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn mcm2(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (kf, kc) = (1, 2);
    let (mut f1r, mut d1r, mut f2r, mut d2r) = (0.0, 0.0, 0.0, 0.0);

    if (m < 0.0) || (m != m.floor()) || (q < 0.0) {
        return (f64::NAN, f64::NAN);
    }
    let int_m = m as i32;
    let status = mtu12(kf, kc, int_m, q, x, &mut f1r, &mut d1r, &mut f2r, &mut d2r);
    if status != Status::Ok {
        return (f64::NAN, f64::NAN);
    }
    (f2r, d2r)
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn msm2(m: f64, q: f64, x: f64) -> (f64, f64) {
    let (kf, kc) = (2, 2);
    let (mut f1r, mut d1r, mut f2r, mut d2r) = (0.0, 0.0, 0.0, 0.0);

    if (m < 1.0) || (m != m.floor()) || (q < 0.0) {
        return (f64::NAN, f64::NAN);
    }
    let int_m = m as i32;
    let status = mtu12(kf, kc, int_m, q, x, &mut f1r, &mut d1r, &mut f2r, &mut d2r);
    if status != Status::Ok {
        return (f64::NAN, f64::NAN);
    }
    (f2r, d2r)
}

/// Characteristic value of even Mathieu functions
///
/// Corresponds to [`scipy.special.mathieu_a`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_a.html
///
/// # See also
/// - [`mathieu_b`]: Characteristic value of odd Mathieu functions
/// - [`mathieu_cem`]: Even Mathieu function
#[doc(alias = "cem_cva")]
#[must_use]
#[inline]
pub fn mathieu_a(m: u32, q: f64) -> f64 {
    cem_cva(f64::from(m), q)
}

/// Characteristic value of odd Mathieu functions
///
/// Corresponds to [`scipy.special.mathieu_b`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_b.html
///
/// # See also
/// - [`mathieu_a`]: Characteristic value of even Mathieu functions
/// - [`mathieu_sem`]: Odd Mathieu function
#[doc(alias = "sem_cva")]
#[must_use]
#[inline]
pub fn mathieu_b(m: u32, q: f64) -> f64 {
    sem_cva(f64::from(m), q)
}

/// Even Mathieu function and its derivative
///
/// Returns the even Mathieu function, *ce<sub>m</sub>(x, q)*, of order *m* and parameter *q*
/// evaluated at *x* (given in degrees).
/// Also returns the derivative with respect to *x* of *ce<sub>m</sub>(x, q)*.
///
/// Corresponds to [`scipy.special.mathieu_cem`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_cem.html
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
/// Both are NaN if more than 250 Fourier coefficients are needed, i.e. for large `|q|` or `m`.
///
/// # See also
/// - [`mathieu_sem`]: Odd Mathieu function
/// - [`mathieu_a`]: Characteristic value of even Mathieu functions
#[doc(alias = "cem")]
#[must_use]
#[inline]
pub fn mathieu_cem(m: u32, q: f64, x: f64) -> (f64, f64) {
    if mathieu_km_too_large(m, q) {
        // (unlike in scipy) NaN, as for km > 251, instead of indexing out of bounds
        return (f64::NAN, f64::NAN);
    }
    cem(f64::from(m), q, x)
}

/// Odd Mathieu function and its derivative
///
/// Returns the odd Mathieu function, *se<sub>m</sub>(x, q)*, of order *m* and parameter *q*
/// evaluated at *x* (given in degrees).
/// Also returns the derivative with respect to *x* of *se<sub>m</sub>(x, q)*.
///
/// Corresponds to [`scipy.special.mathieu_sem`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_sem.html
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
/// Both are NaN if more than 250 Fourier coefficients are needed, i.e. for large `|q|` or `m`.
///
/// # See also
/// - [`mathieu_cem`]: Even Mathieu function
/// - [`mathieu_b`]: Characteristic value of odd Mathieu functions
#[doc(alias = "sem")]
#[must_use]
#[inline]
pub fn mathieu_sem(m: u32, q: f64, x: f64) -> (f64, f64) {
    if m > 0 && mathieu_km_too_large(m, q) {
        // (unlike in scipy) NaN, as for km > 251, instead of indexing out of bounds
        return (f64::NAN, f64::NAN);
    }
    sem(f64::from(m), q, x)
}

/// Even modified Mathieu function of the first kind and its derivative
///
/// Evaluates the even modified Mathieu function of the first kind, *Mc1<sub>m</sub>(x, q)*,
/// and its derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// Corresponds to [`scipy.special.mathieu_modcem1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modcem1.html
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
/// Both are NaN if more than 250 Fourier coefficients are needed, i.e. for large `|q|` or `m`.
///
/// # See also
/// - [`mathieu_cem`]: Even Mathieu function
/// - [`mathieu_modcem2`]: Even modified Mathieu function of the second kind
/// - [`mathieu_modsem1`]: Odd modified Mathieu function of the first kind
#[doc(alias = "mcm1")]
#[must_use]
#[inline]
pub fn mathieu_modcem1(m: u32, q: f64, x: f64) -> (f64, f64) {
    if mathieu_km_too_large(m, q) {
        // (unlike in scipy) instead of indexing out of bounds if `km` overflows `int`
        return (f64::NAN, f64::NAN);
    }
    mcm1(f64::from(m), q, x)
}

/// Odd modified Mathieu function of the first kind and its derivative
///
/// Evaluates the odd modified Mathieu function of the first kind, *Ms1<sub>m</sub>(x, q)*, and its
/// derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// Corresponds to [`scipy.special.mathieu_modsem1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modsem1.html
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
/// Both are NaN if more than 250 Fourier coefficients are needed, i.e. for large `|q|` or `m`.
///
/// # See also
/// - [`mathieu_sem`]: Odd Mathieu function
/// - [`mathieu_modsem2`]: Odd modified Mathieu function of the second kind
/// - [`mathieu_modcem1`]: Even modified Mathieu function of the first kind
#[doc(alias = "msm1")]
#[must_use]
#[inline]
pub fn mathieu_modsem1(m: u32, q: f64, x: f64) -> (f64, f64) {
    if mathieu_km_too_large(m, q) {
        // (unlike in scipy) instead of indexing out of bounds if `km` overflows `int`
        return (f64::NAN, f64::NAN);
    }
    msm1(f64::from(m), q, x)
}

/// Even modified Mathieu function of the second kind and its derivative
///
/// Evaluates the even modified Mathieu function of the second kind, *Mc2<sub>m</sub>(x, q)*, and
/// its derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// Corresponds to [`scipy.special.mathieu_modcem2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modcem2.html
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
/// Both are NaN if more than 250 Fourier coefficients are needed, i.e. for large `|q|` or `m`.
///
/// # See also
/// - [`mathieu_cem`]: Even Mathieu function
/// - [`mathieu_modcem1`]: Even modified Mathieu function of the first kind
/// - [`mathieu_modsem2`]: Odd modified Mathieu function of the second kind
#[doc(alias = "mcm2")]
#[must_use]
#[inline]
pub fn mathieu_modcem2(m: u32, q: f64, x: f64) -> (f64, f64) {
    if mathieu_km_too_large(m, q) {
        // (unlike in scipy) instead of indexing out of bounds if `km` overflows `int`
        return (f64::NAN, f64::NAN);
    }
    mcm2(f64::from(m), q, x)
}

/// Odd modified Mathieu function of the second kind and its derivative
///
/// Evaluates the odd modified Mathieu function of the second kind, *Ms2<sub>m</sub>(x, q)*, and
/// its derivative at *x* (given in degrees) for order *m* and parameter *q*.
///
/// Corresponds to [`scipy.special.mathieu_modsem2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modsem2.html
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
/// Both are NaN if more than 250 Fourier coefficients are needed, i.e. for large `|q|` or `m`.
///
/// # See also
/// - [`mathieu_sem`]: Odd Mathieu function
/// - [`mathieu_modsem1`]: Odd modified Mathieu function of the first kind
/// - [`mathieu_modcem2`]: Even modified Mathieu function of the second kind
#[doc(alias = "msm2")]
#[must_use]
#[inline]
pub fn mathieu_modsem2(m: u32, q: f64, x: f64) -> (f64, f64) {
    if mathieu_km_too_large(m, q) {
        // (unlike in scipy) instead of indexing out of bounds if `km` overflows `int`
        return (f64::NAN, f64::NAN);
    }
    msm2(f64::from(m), q, x)
}

/// `qm` in `mathieu_{even,odd}_coef`, and in `fcoef` and `mtu0` of `specfun.h`, where the number
/// of predicted Fourier coefficients is `km = int(qm + 0.5*m)`
fn mathieu_qm(q: f64) -> f64 {
    let q_sqrt = q.sqrt();
    if q <= 1.0 {
        7.5 + 56.1 * q_sqrt - 134.7 * q + 90.7 * q_sqrt * q
    } else {
        17.0 + 3.1 * q_sqrt - 0.126 * q + 0.0037 * q_sqrt * q
    }
}

/// Whether `km >= 251`, in which case `mtu0` (used by `cem` and `sem`) and `mtu12` (used by the
/// modified Mathieu functions) would index out of bounds: `mtu0` only returns NaN if `km > 251`,
/// and both cast `km` to `int`, which overflows for huge or NaN `q`
fn mathieu_km_too_large(m: u32, q: f64) -> bool {
    // `q < 0` is reflected to `-q`, with the same `m`
    let km = mathieu_qm(q.abs()) + 0.5 * f64::from(m);
    km.is_nan() || km >= 251.0
}

#[inline]
fn fcoef<const KD: i32>(m: u32, q: f64) -> Result<([f64; 252], usize), String> {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L1582-L1639

    // if (q < 0):
    if q < 0.0 {
        // raise ValueError("q >=0")
        return Err("q >=0".to_string());
    }

    let qm = mathieu_qm(q);
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
    // (unlike in scipy) 252 instead of 251, since `fcoef` reads `fc[km]`, and `km` can be 251
    let mut fc = [f64::NAN; 252];
    specfun_fcoef(KD, m.try_into().unwrap(), q, cv, &mut fc);
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
/// - `m == 0`
/// - `q < 0`
/// - too many predicted coefficients
///
/// # See also
/// - [`mathieu_even_coef`]: Fourier coefficients for even Mathieu and modified Mathieu functions
#[inline]
pub fn mathieu_odd_coef(m: u32, q: f64) -> Result<Vec<f64>, String> {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L1642-L1698

    // if (m != floor(m)) or (m <= 0):
    if m == 0 {
        // raise ValueError("m must be an integer > 0")
        return Err("m must be an integer > 0".to_string());
    }

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
    fn test_mathieu_a_nan() {
        assert!(crate::mathieu_a(13, f64::NAN).is_nan());
    }

    #[test]
    fn test_mathieu_modcem1_large_x() {
        // `(f² + (f'/u)²)·u → 2/π` with `u = √q·exp(x)`, i.e. `|Mc1| ~ √(2 / (π u))`
        let x = 22.0;
        let u = f64::exp(x);
        let (f, d) = crate::mathieu_modcem1(0, 1.0, x);
        let env = (f * f + (d / u) * (d / u)) * u;
        assert!(
            (env * core::f64::consts::FRAC_PI_2 - 1.0).abs() < 1e-6,
            "{env}"
        );
    }

    #[test]
    fn test_mathieu_even_coef() {
        // Smoke test: verify function returns Ok and non-empty coefficients
        let result = crate::mathieu_even_coef(2, 1.0);
        assert!(result.is_ok());
        let coefs = result.unwrap();
        assert_ne!(coefs, [] as [f64; 0]);
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
        assert_ne!(coefs, [] as [f64; 0]);
        assert!(coefs.iter().all(|&c| c.is_finite()));

        // Test with q = 0
        let result = crate::mathieu_odd_coef(1, 0.0);
        assert!(result.is_ok());

        // Test error case: negative q
        let result = crate::mathieu_odd_coef(1, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_mathieu_odd_coef_m0() {
        assert!(crate::mathieu_odd_coef(0, 0.0).is_err());
    }

    #[test]
    fn test_mathieu_cem_sem_km251() {
        // `mtu0` would read out of bounds
        for (y, yp) in [
            crate::mathieu_cem(0, 2103.0, 10.0),
            crate::mathieu_cem(0, -2103.0, 10.0),
            crate::mathieu_sem(1, 2100.0, 10.0),
            crate::mathieu_sem(1, -2100.0, 10.0),
        ] {
            assert!(y.is_nan() && yp.is_nan());
        }
        assert_eq!(crate::mathieu_sem(0, 2103.0, 10.0), (0.0, 0.0));
    }

    #[test]
    fn test_mathieu_huge_q() {
        // `km` would overflow `int`
        for q in [1e10, f64::INFINITY, f64::NAN] {
            for (y, yp) in [
                crate::mathieu_cem(0, q, 10.0),
                crate::mathieu_sem(1, q, 10.0),
                crate::mathieu_modcem1(0, q, 1.0),
                crate::mathieu_modsem1(1, q, 1.0),
                crate::mathieu_modcem2(0, q, 1.0),
                crate::mathieu_modsem2(1, q, 1.0),
            ] {
                assert!(y.is_nan() && yp.is_nan(), "q = {q}");
            }
        }
    }
}
