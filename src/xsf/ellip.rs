/// Complete elliptic integral of the first kind
///
/// Corresponds to [`scipy.special.ellipk`][ellipk] in SciPy
///
/// [ellipk]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html
///
/// # See also
/// - [`ellipe`](crate::ellipe): Complete elliptic integral of the second kind
/// - [`ellipkinc`](crate::ellipkinc): Incomplete elliptic integral of the first kind
/// - [`ellipkm1`](crate::ellipkm1): Complete elliptic integral of the first kind around $m = 1$
/// - [`ellipj`](crate::ellipj): Jacobi elliptic functions
#[doc(alias = "ellip_k")]
#[must_use]
#[inline]
pub fn ellipk(m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipk(m) }
}

/// Complete elliptic integral of the first kind around $m = 1$
///
/// Corresponds to [`scipy.special.ellipkm1`][ellipkm1] in SciPy
///
/// [ellipkm1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkm1.html
///
/// # See also
/// - [`ellipk`](crate::ellipk): Complete elliptic integral of the first kind
/// - [`ellipkinc`](crate::ellipkinc): Incomplete elliptic integral of the first kind
/// - [`ellipj`](crate::ellipj): Jacobi elliptic functions
#[doc(alias = "ellip_k_m1")]
#[must_use]
#[inline]
pub fn ellipkm1(p: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipkm1(p) }
}

/// Incomplete elliptic integral of the first kind
///
/// Corresponds to [`scipy.special.ellipkinc`][ellipkinc] in SciPy
///
/// [ellipkinc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkinc.html
///
/// # See also
/// - [`ellipk`](crate::ellipk): Complete elliptic integral of the first kind
/// - [`ellipkm1`](crate::ellipkm1): Complete elliptic integral of the first kind around $m = 1$
/// - [`ellipj`](crate::ellipj): Jacobi elliptic functions
#[doc(alias = "ellip_k_inc")]
#[must_use]
#[inline]
pub fn ellipkinc(phi: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipkinc(phi, m) }
}

/// Complete elliptic integral of the second kind
///
/// Corresponds to [`scipy.special.ellip_e`][ellipe] in SciPy
///
/// [ellipe]: https://docs.scipy.org/doc/scipy/reference/generated/scipy
///
/// # See also
/// - [`ellipk`](crate::ellipk): Complete elliptic integral of the first kind
/// - [`ellipeinc`](crate::ellipeinc): Incomplete elliptic integral of the second kind
/// - [`ellipj`](crate::ellipj): Jacobi elliptic functions
#[doc(alias = "ellip_e")]
#[must_use]
#[inline]
pub fn ellipe(m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipe(m) }
}

/// Incomplete elliptic integral of the second kind
///
/// Corresponds to [`scipy.special.ellip_e_inc`][ellip_e_inc] in SciPy
///
/// [ellip_e_inc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special
///
/// # See also
/// - [`ellipe`](crate::ellipe): Complete elliptic integral of the second kind
/// - [`ellipk`](crate::ellipk): Complete elliptic integral of the first kind
/// - [`ellipj`](crate::ellipj): Jacobi elliptic functions
#[doc(alias = "ellip_e_inc")]
#[must_use]
#[inline]
pub fn ellipeinc(phi: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::ellipeinc(phi, m) }
}

/// Jacobi elliptic functions
///
/// Corresponds to [`scipy.special.ellipj`][ellipj] in SciPy
///
/// [ellipj]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipj.html
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
/// - [`ellipe`] - Complete elliptic integral of the second kind
/// - [`ellipkinc`] - Incomplete elliptic integral of the first kind
/// - [`ellipeinc`] - Incomplete elliptic integral of the second kind
#[doc(alias = "ellip_j")]
#[must_use]
#[inline]
pub fn ellipj(u: f64, m: f64) -> (f64, f64, f64, f64) {
    let (mut sn, mut cn, mut dn, mut am) = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::ellipj(u, m, &raw mut sn, &raw mut cn, &raw mut dn, &raw mut am);
    }
    (sn, cn, dn, am)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_ellipk() {
        xsref::test("ellipk", "d-d", |x| crate::ellipk(x[0]));
    }

    #[test]
    fn test_ellipkm1() {
        xsref::test("ellipkm1", "d-d", |x| crate::ellipkm1(x[0]));
    }

    #[test]
    fn test_ellipkinc() {
        xsref::test("ellipkinc", "d_d-d", |x| crate::ellipkinc(x[0], x[1]));
    }

    #[test]
    fn test_ellipe() {
        xsref::test("ellipe", "d-d", |x| crate::ellipe(x[0]));
    }

    #[test]
    fn test_ellipeinc() {
        xsref::test("ellipeinc", "d_d-d", |x| crate::ellipeinc(x[0], x[1]));
    }

    #[test]
    fn test_ellipj() {
        // Workaround for https://github.com/scipy/xsref/issues/11:
        // The xsref table contains incorrect values for the phase φ (am),
        // so we reconstruct it from sn using sn.asin()
        xsref::test("ellipj", "d_d-d_d_d_d", |x| {
            let (sn, cn, dn, am) = crate::ellipj(x[0], x[1]);
            (sn, cn, dn, am.sin().asin())
        });
    }
}
