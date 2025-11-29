#![allow(clippy::cast_precision_loss)]

/// Characteristic value of prolate spheroidal function
#[doc(alias = "pro_cv")]
#[must_use]
#[inline]
pub fn prolate_segv(m: u64, n: u64, c: f64) -> f64 {
    unsafe { crate::ffi::xsf::prolate_segv(m as f64, n as f64, c) }
}

/// Characteristic value of oblate spheroidal function
#[doc(alias = "obl_cv")]
#[must_use]
#[inline]
pub fn oblate_segv(m: u64, n: u64, c: f64) -> f64 {
    unsafe { crate::ffi::xsf::oblate_segv(m as f64, n as f64, c) }
}

/// Prolate spheroidal angular function of the 1st kind and its derivative
///
/// Computes the prolate spheroidal angular function of the first kind and its derivative (with
/// respect to `x`).
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument (*|x| < 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "pro_ang1")]
#[must_use]
#[inline]
pub fn prolate_aswfa_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_aswfa_nocv(m as f64, n as f64, c, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Oblate spheroidal angular function of the 1st kind and its derivative
///
/// Computes the oblate spheroidal angular function of the first kind and its derivative (with
/// respect to `x`).
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument (*|x| < 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "obl_ang1")]
#[must_use]
#[inline]
pub fn oblate_aswfa_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_aswfa_nocv(m as f64, n as f64, c, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Prolate spheroidal radial function of the 1st kind and its derivative
///
/// Computes the prolate spheroidal radial function of the first kind and its derivative (with
/// respect to `x`).
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "pro_rad1")]
#[must_use]
#[inline]
pub fn prolate_radial1_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial1_nocv(m as f64, n as f64, c, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Oblate spheroidal radial function of the 1st kind and its derivative
///
/// Computes the oblate spheroidal radial function of the first kind and its derivative (with
/// respect to `x`).
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "obl_rad1")]
#[must_use]
#[inline]
pub fn oblate_radial1_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial1_nocv(m as f64, n as f64, c, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Prolate spheroidal radial function of the 2nd kind and its derivative
///
/// Computes the prolate spheroidal radial function of the second kind and its derivative (with
/// respect to `x`).
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "pro_rad2")]
#[must_use]
#[inline]
pub fn prolate_radial2_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial2_nocv(m as f64, n as f64, c, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Oblate spheroidal radial function of the second kind and its derivative
///
/// Computes the oblate spheroidal radial function of the second kind and its derivative (with
/// respect to `x`).
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "obl_rad2")]
#[must_use]
#[inline]
pub fn oblate_radial2_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial2_nocv(m as f64, n as f64, c, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Prolate spheroidal angular function for precomputed characteristic value
///
/// Computes the prolate spheroidal angular function of the first kind and its derivative (with
/// respect to `x`). Requires pre-computed characteristic value.
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument (*|x| < 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "pro_ang1_cv")]
#[must_use]
#[inline]
pub fn prolate_aswfa(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_aswfa(m as f64, n as f64, c, cv, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Oblate spheroidal angular function for precomputed characteristic value
///
/// Computes the oblate spheroidal angular function of the first kind and its derivative (with
/// respect to `x`). Requires pre-computed characteristic value.
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument (*|x| < 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "obl_ang1_cv")]
#[must_use]
#[inline]
pub fn oblate_aswfa(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_aswfa(m as f64, n as f64, c, cv, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Prolate spheroidal radial function of the 1st kind for precomputed characteristic value
///
/// Computes the prolate spheroidal radial function of the first kind and its derivative (with
/// respect to `x`). Requires pre-computed characteristic value.
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "pro_rad1_cv")]
#[must_use]
#[inline]
pub fn prolate_radial1(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial1(m as f64, n as f64, c, cv, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Oblate spheroidal radial function of the 1st kind for precomputed characteristic value
///
/// Computes the oblate spheroidal radial function of the first kind and its derivative (with
/// respect to `x`). Requires pre-computed characteristic value.
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "obl_rad1_cv")]
#[must_use]
#[inline]
pub fn oblate_radial1(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial1(m as f64, n as f64, c, cv, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Prolate spheroidal radial function of the 2nd kind for precomputed characteristic value
///
/// Computes the prolate spheroidal radial function of the second kind and its derivative (with
/// respect to `x`). Requires pre-computed characteristic value.
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "pro_rad2_cv")]
#[must_use]
#[inline]
pub fn prolate_radial2(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial2(m as f64, n as f64, c, cv, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

/// Oblate spheroidal radial function of the 2nd kind for precomputed characteristic value
///
/// Computes the oblate spheroidal radial function of the second kind and its derivative (with
/// respect to `x`). Requires pre-computed characteristic value.
///
/// # Arguments
/// - `m`: Mode parameter (*m >= 0*)
/// - `n`: Mode parameter (*n >= m*)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument (*x > 1*)
///
/// # Returns
/// - `s`: Value of the function
/// - `sp` : Value of the derivative w.r.t. `x`
#[doc(alias = "obl_rad2_cv")]
#[must_use]
#[inline]
pub fn oblate_radial2(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial2(m as f64, n as f64, c, cv, x, &raw mut s, &raw mut sp);
    }
    (s, sp)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

    #[test]
    fn test_prolate_segv() {
        xsref::test("prolate_segv", "d_d_d-d", |x| {
            crate::prolate_segv(x[0] as u64, x[1] as u64, x[2])
        });
    }

    // https://github.com/scipy/xsref/issues/9
    // #[test]
    // fn test_oblate_segv() {
    //     testing::test("oblate_segv", "d_d_d-d", |x| {
    //         crate::oblate_segv(x[0] as u64, x[1] as u64, x[2])
    //     });
    // }

    #[test]
    fn test_prolate_aswfa_nocv() {
        xsref::test("prolate_aswfa_nocv", "d_d_d_d-d_d", |x| {
            crate::prolate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_aswfa_nocv() {
        xsref::test("oblate_aswfa_nocv", "d_d_d_d-d_d", |x| {
            crate::oblate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_radial1_nocv() {
        xsref::test("prolate_radial1_nocv", "d_d_d_d-d_d", |x| {
            crate::prolate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_radial1_nocv() {
        xsref::test("oblate_radial1_nocv", "d_d_d_d-d_d", |x| {
            crate::oblate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_radial2_nocv() {
        xsref::test("prolate_radial2_nocv", "d_d_d_d-d_d", |x| {
            crate::prolate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_radial2_nocv() {
        xsref::test("oblate_radial2_nocv", "d_d_d_d-d_d", |x| {
            crate::oblate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_aswfa() {
        xsref::test("prolate_aswfa", "d_d_d_d_d-d_d", |x| {
            crate::prolate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_aswfa() {
        xsref::test("oblate_aswfa", "d_d_d_d_d-d_d", |x| {
            crate::oblate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial1() {
        xsref::test("prolate_radial1", "d_d_d_d_d-d_d", |x| {
            crate::prolate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_radial1() {
        xsref::test("oblate_radial1", "d_d_d_d_d-d_d", |x| {
            crate::oblate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial2() {
        xsref::test("prolate_radial2", "d_d_d_d_d-d_d", |x| {
            crate::prolate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_radial2() {
        xsref::test("oblate_radial2", "d_d_d_d_d-d_d", |x| {
            crate::oblate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }
}
