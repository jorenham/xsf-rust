pub trait ZetaArg: crate::sealed::Sealed {
    fn riemann_zeta(self) -> Self;
    fn zeta(self, q: f64) -> Self;
}

impl ZetaArg for f64 {
    #[inline]
    fn riemann_zeta(self) -> Self {
        unsafe { crate::ffi::xsf::riemann_zeta(self) }
    }

    #[inline]
    fn zeta(self, q: f64) -> Self {
        crate::xsf::cephes::zeta(self, q)
    }
}

impl ZetaArg for num_complex::Complex<f64> {
    #[inline]
    fn riemann_zeta(self) -> Self {
        unsafe { crate::ffi::xsf::riemann_zeta_1(self) }
    }

    #[inline]
    fn zeta(self, q: f64) -> Self {
        if self.im == 0.0 {
            num_complex::Complex::new(self.re.zeta(q), 0.0)
        } else {
            // Complex input for Hurwitz Zeta is not currently supported.
            num_complex::Complex::new(f64::NAN, f64::NAN)
        }
    }
}

/// Riemann zeta function $\zeta(z)$ for real or complex $z$
///
/// Corresponds to [`scipy.special.zeta(z)`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html
///
/// # Definition
///
/// $$ \zeta(z) = \sum_{k=1}^\infty k^{-z} $$
///
/// # See also
/// - [`zeta`]: Hurwitz zeta function $\zeta(z, q)$
/// - [`zetac`]: Riemann zeta function minus one, $\zeta(x) - 1$
#[must_use]
#[inline]
pub fn riemann_zeta<T: ZetaArg>(z: T) -> T {
    z.riemann_zeta()
}

/// Hurwitz zeta function $\zeta(z, q)$ for real or complex $z$
///
/// Corresponds to [`scipy.special.zeta(z, q)`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html
///
/// # Definition
///
/// $$ \zeta(z, q) = \sum_{k=0}^\infty (k + q)^{-z} $$
///
/// # See also
/// - [`riemann_zeta`]: Riemann zeta function $\zeta(z)$
/// - [`zetac`]: Riemann zeta function minus one, $\zeta(x) - 1$
#[doc(alias = "hurwitz_zeta")]
#[must_use]
#[inline]
pub fn zeta<T: ZetaArg>(z: T, q: f64) -> T {
    z.zeta(q)
}

/// Riemann zeta function minus one, $\zeta(x) - 1$
///
/// Corresponds to [`scipy.special.zetac`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zetac.html
///
/// # Definition
///
/// $$ \zeta(x) - 1 = \sum_{k=2}^\infty k^{-x} $$
///
/// # See also
/// - [`riemann_zeta`]: Riemann zeta function $\zeta(z)$
/// - [`zeta`]: Hurwitz zeta function $\zeta(z, q)$
#[doc(alias = "riemann_zetac")]
#[must_use]
#[inline]
pub fn zetac(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::zetac(x) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_riemann_zeta_f64() {
        xsref::test("riemann_zeta", "d-d", |x| crate::riemann_zeta(x[0]));
    }

    #[test]
    fn test_riemann_zeta_c64() {
        xsref::test("riemann_zeta", "cd-cd", |x| {
            crate::riemann_zeta(c64(x[0], x[1]))
        });
    }

    #[test]
    fn test_zeta_f64() {
        xsref::test("zeta", "d_d-d", |x| crate::zeta(x[0], x[1]));
    }

    #[test]
    fn test_zeta_large_q() {
        // asymptotic expansion for q > 1e8 (reference value from mpmath)
        let (x, q) = (36.0, 1.1e8);
        let expected = 1.016_688_811_413_570_4e-283;
        let rel = (crate::zeta(x, q) - expected).abs() / expected;
        assert!(rel < 1e-15, "rel = {rel:e}");

        assert_eq!(crate::zeta(f64::INFINITY, q), 0.0);
        assert_eq!(crate::zeta(f64::INFINITY, f64::INFINITY), 0.0);
    }

    #[test]
    fn test_zetac_f64() {
        xsref::test("zetac", "d-d", |x| crate::zetac(x[0]));
    }
}
