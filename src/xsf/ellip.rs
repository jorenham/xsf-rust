/// Complete elliptic integral of the first kind
///
/// Corresponds to [`scipy.special.ellipk`][ellipk] in SciPy
///
/// [ellipk]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html
///
/// # See also
/// - [`ellipe`]: Complete elliptic integral of the second kind
/// - [`ellipkinc`]: Incomplete elliptic integral of the first kind
/// - [`ellipkm1`]: Complete elliptic integral of the first kind around $m = 1$
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_k")]
#[must_use]
#[inline]
pub fn ellipk(m: f64) -> f64 {
    crate::xsf::cephes::ellpk(1.0 - m)
}

/// Complete elliptic integral of the first kind around $m = 1$
///
/// Corresponds to [`scipy.special.ellipkm1`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkm1.html
///
/// # See also
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipkinc`]: Incomplete elliptic integral of the first kind
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_k_m1")]
#[must_use]
#[inline]
pub fn ellipkm1(p: f64) -> f64 {
    crate::xsf::cephes::ellpk(p)
}

/// Incomplete elliptic integral of the first kind
///
/// Corresponds to [`scipy.special.ellipkinc`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkinc.html
///
/// # See also
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipkm1`]: Complete elliptic integral of the first kind around $m = 1$
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_k_inc")]
#[must_use]
#[inline]
pub fn ellipkinc(phi: f64, m: f64) -> f64 {
    crate::xsf::cephes::ellik(phi, m)
}

/// Complete elliptic integral of the second kind
///
/// Corresponds to [`scipy.special.ellipe`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html
///
/// # See also
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipeinc`]: Incomplete elliptic integral of the second kind
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_e")]
#[must_use]
#[inline]
pub fn ellipe(m: f64) -> f64 {
    crate::xsf::cephes::ellpe(m)
}

/// Incomplete elliptic integral of the second kind
///
/// Corresponds to [`scipy.special.ellipeinc`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipeinc.html
///
/// # See also
/// - [`ellipe`]: Complete elliptic integral of the second kind
/// - [`ellipk`]: Complete elliptic integral of the first kind
/// - [`ellipj`]: Jacobi elliptic functions
#[doc(alias = "ellip_e_inc")]
#[must_use]
#[inline]
pub fn ellipeinc(phi: f64, m: f64) -> f64 {
    crate::xsf::cephes::ellie(phi, m)
}

/// Jacobi elliptic functions
///
/// Corresponds to [`scipy.special.ellipj`][scipy] in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipj.html
///
/// # Arguments
/// - `u` - Real argument
/// - `m` - Parameter, `m <= 1`
///
/// # Returns
/// - *sn(u | m)* - sine amplitude
/// - *cn(u | m)* - cosine amplitude
/// - *dn(u | m)* - delta amplitude
/// - phase *φ* s.t. *sn(u | m) = sin(φ)* and *cn(u | m) = cos(φ)*
///
/// # See also
/// - [`ellipk`] - Complete elliptic integral of the first kind
/// - [`ellipe`] - Complete elliptic integral of the second kind
/// - [`ellipkinc`] - Incomplete elliptic integral of the first kind
/// - [`ellipeinc`] - Incomplete elliptic integral of the second kind
#[doc(alias = "ellip_j")]
#[must_use]
#[inline]
pub fn ellipj(u: f64, m: f64) -> (f64, f64, f64, f64) {
    crate::xsf::cephes::ellpj(u, m)
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
    fn test_ellipkinc_ellipeinc_huge_phi() {
        // phi - npio2 * pi/2 isn't in [-pi/2, pi/2] (reference values from mpmath)
        let (phi, m) = (-1.788_575_894_471_164_1e19, -1.244_349_504_000_301_7e35);
        for (got, expected) in [
            (crate::ellipkinc(phi, m), -1_348.954_827_564_235_2),
            (crate::ellipeinc(phi, m), -4.016_600_941_581_635_3e36),
        ] {
            let rel = (got - expected).abs() / expected.abs();
            assert!(rel < 1e-15, "rel = {rel:e}");
        }
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
