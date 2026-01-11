//! NOTE: The Cephes library only provides a real implementation of the Spence function.
//! For complex argument, `scipy.special.spence` uses the Cython implementation at
//! `scipy/special/_spence.pxd`.

/// Spence's function, also known as the dilogarithm.
///
/// $$ \int_1^x {\log(t) \over 1 - t} \dd t $$
///
/// Corresponds to [`scipy.special.spence`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spence.html
#[doc(alias = "dilogarithm")]
#[must_use]
#[inline]
pub fn spence(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::spence(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_spence() {
        xsref::test("spence", "d-d", |x| crate::spence(x[0]));
    }
}
