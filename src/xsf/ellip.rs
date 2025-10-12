/// Complete elliptic integral of the first kind
#[doc(alias = "ellip_k")]
pub fn ellipk(m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipk(m) }
}

/// Complete elliptic integral of the first kind around `m = 1`
pub fn ellipkm1(p: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipkm1(p) }
}

/// Incomplete elliptic integral of the first kind
#[doc(alias = "ellip_k_inc")]
pub fn ellipkinc(phi: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipkinc(phi, m) }
}

/// Complete elliptic integral of the second kind
#[doc(alias = "ellip_e")]
pub fn ellipe(m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipe(m) }
}

/// Incomplete elliptic integral of the second kind
#[doc(alias = "ellip_e_inc")]
pub fn ellipeinc(phi: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipeinc(phi, m) }
}

/// Jacobi elliptic functions
///
/// # Arguments
/// - `u` - Real argument
/// - `m` - Parameter between 0 and 1
///
/// # Returns
/// - *sn(u | m)* - sine amplitude
/// - *cn(u | m)* - cosine amplitude
/// - *dn(u | m)* - delta amplitude
/// - phase *φ* s.t. *sn(u | m) = sin(φ)* and *cn(u | m) = cos(φ)*
///
/// # See Also
/// - [`ellipk`] - Complete elliptic integral of the first kind
/// - [`ellipkinc`] - Incomplete elliptic integral of the first kind
pub fn ellipj(u: f64, m: f64) -> (f64, f64, f64, f64) {
    let mut sn = f64::NAN;
    let mut cn = f64::NAN;
    let mut dn = f64::NAN;
    let mut phi = f64::NAN;
    unsafe {
        crate::ffi::xsf::ellipj(u, m, &mut sn, &mut cn, &mut dn, &mut phi);
    }
    (sn, cn, dn, phi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    // ellipk*

    #[test]
    fn test_ellipk() {
        xsref::test::<f64, _>("ellipk", "d-d", |x: &[f64]| ellipk(x[0]));
    }

    #[test]
    fn test_ellipkm1() {
        xsref::test::<f64, _>("ellipkm1", "d-d", |x: &[f64]| ellipkm1(x[0]));
    }

    #[test]
    fn test_ellipkinc() {
        xsref::test::<f64, _>("ellipkinc", "d_d-d", |x: &[f64]| ellipkinc(x[0], x[1]));
    }

    // ellipe*

    #[test]
    fn test_ellipe() {
        xsref::test::<f64, _>("ellipe", "d-d", |x: &[f64]| ellipe(x[0]));
    }

    #[test]
    fn test_ellipeinc() {
        xsref::test::<f64, _>("ellipeinc", "d_d-d", |x: &[f64]| ellipeinc(x[0], x[1]));
    }

    // ellipj

    #[test]
    fn test_ellipj() {
        xsref::test::<(f64, f64, f64, f64), _>("ellipj", "d_d-d_d_d_d", |x: &[f64]| {
            ellipj(x[0], x[1])
        });
    }
}
