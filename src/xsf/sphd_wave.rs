/// Characteristic value of prolate spheroidal function
#[doc(alias = "pro_cv")]
pub fn prolate_segv(m: u64, n: u64, c: f64) -> f64 {
    unsafe { crate::ffi::xsf::prolate_segv(m as f64, n as f64, c) }
}

/// Characteristic value of oblate spheroidal function
#[doc(alias = "obl_cv")]
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
pub fn prolate_aswfa_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_aswfa_nocv(m as f64, n as f64, c, x, &mut s, &mut sp);
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
pub fn oblate_aswfa_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_aswfa_nocv(m as f64, n as f64, c, x, &mut s, &mut sp);
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
pub fn prolate_radial1_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial1_nocv(m as f64, n as f64, c, x, &mut s, &mut sp);
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
pub fn oblate_radial1_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial1_nocv(m as f64, n as f64, c, x, &mut s, &mut sp);
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
pub fn prolate_radial2_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial2_nocv(m as f64, n as f64, c, x, &mut s, &mut sp);
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
pub fn oblate_radial2_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial2_nocv(m as f64, n as f64, c, x, &mut s, &mut sp);
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
pub fn prolate_aswfa(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_aswfa(m as f64, n as f64, c, cv, x, &mut s, &mut sp);
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
pub fn oblate_aswfa(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_aswfa(m as f64, n as f64, c, cv, x, &mut s, &mut sp);
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
pub fn prolate_radial1(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial1(m as f64, n as f64, c, cv, x, &mut s, &mut sp);
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
pub fn oblate_radial1(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial1(m as f64, n as f64, c, cv, x, &mut s, &mut sp);
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
pub fn prolate_radial2(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::prolate_radial2(m as f64, n as f64, c, cv, x, &mut s, &mut sp);
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
pub fn oblate_radial2(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    let (mut s, mut sp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::oblate_radial2(m as f64, n as f64, c, cv, x, &mut s, &mut sp);
    }
    (s, sp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_prolate_segv() {
        xsref::test::<f64, _>("prolate_segv", "d_d_d-d", |x: &[f64]| {
            prolate_segv(x[0] as u64, x[1] as u64, x[2])
        });
    }

    // https://github.com/scipy/xsref/issues/9
    // #[test]
    // fn test_oblate_segv() {
    //     testing::test::<f64, _>("oblate_segv", "d_d_d-d", |x: &[f64]| {
    //         oblate_segv(x[0] as u64, x[1] as u64, x[2])
    //     });
    // }

    #[test]
    fn test_prolate_aswfa_nocv() {
        xsref::test::<(f64, f64), _>("prolate_aswfa_nocv", "d_d_d_d-d_d", |x: &[f64]| {
            prolate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_aswfa_nocv() {
        xsref::test::<(f64, f64), _>("oblate_aswfa_nocv", "d_d_d_d-d_d", |x: &[f64]| {
            oblate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_radial1_nocv() {
        xsref::test::<(f64, f64), _>("prolate_radial1_nocv", "d_d_d_d-d_d", |x: &[f64]| {
            prolate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_radial1_nocv() {
        xsref::test::<(f64, f64), _>("oblate_radial1_nocv", "d_d_d_d-d_d", |x: &[f64]| {
            oblate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_radial2_nocv() {
        xsref::test::<(f64, f64), _>("prolate_radial2_nocv", "d_d_d_d-d_d", |x: &[f64]| {
            prolate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_radial2_nocv() {
        xsref::test::<(f64, f64), _>("oblate_radial2_nocv", "d_d_d_d-d_d", |x: &[f64]| {
            oblate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_aswfa() {
        xsref::test::<(f64, f64), _>("prolate_aswfa", "d_d_d_d_d-d_d", |x: &[f64]| {
            prolate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_aswfa() {
        xsref::test::<(f64, f64), _>("oblate_aswfa", "d_d_d_d_d-d_d", |x: &[f64]| {
            oblate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial1() {
        xsref::test::<(f64, f64), _>("prolate_radial1", "d_d_d_d_d-d_d", |x: &[f64]| {
            prolate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_radial1() {
        xsref::test::<(f64, f64), _>("oblate_radial1", "d_d_d_d_d-d_d", |x: &[f64]| {
            oblate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial2() {
        xsref::test::<(f64, f64), _>("prolate_radial2", "d_d_d_d_d-d_d", |x: &[f64]| {
            prolate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_radial2() {
        xsref::test::<(f64, f64), _>("oblate_radial2", "d_d_d_d_d-d_d", |x: &[f64]| {
            oblate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }
}
