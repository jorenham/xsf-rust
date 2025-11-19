/// Inverse of the error function [*erf(x)*](crate::erf)
///
/// In the complex domain, there is no unique complex number w satisfying erf(w)=z.
/// This indicates a true inverse function would be multivalued. When the domain restricts to the
/// real, -1 < x < 1, there is a unique real number satisfying erf(erfinv(x))=x.
///
/// Note that unlike [`scipy.special.erfinv`][scipy], which uses the Boost implementation,
/// this function uses the Cephes implementation, which can be less accurate in certain regions.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfinv.html
#[doc(alias = "erf_inv")]
pub fn erfinv(y: f64) -> f64 {
    unsafe { crate::ffi::xsf::erfinv(y) }
}

/// Inverse of the complementary error function [*erfc(x)*](crate::erfc)
///
/// In the complex domain, there is no unique complex number w satisfying erfc(w)=z.
/// This indicates a true inverse function would be multivalued. When the domain restricts to the
/// real, 0 < x < 2, there is a unique real number satisfying erfc(erfcinv(x))=erfcinv(erfc(x)).
///
/// It is related to inverse of the error function by erfcinv(1 - y) = [erfinv(y)](fn.erfinv.html).
///
/// Note that [`scipy.special.erfcinv`][scipy] also uses the Cephes implementation.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcinv.html
#[doc(alias = "erfc_inv")]
pub fn erfcinv(y: f64) -> f64 {
    unsafe { crate::ffi::xsf::erfcinv(y) }
}

#[cfg(test)]
mod tests {
    /// Based on `scipy.special.tests.test_TestInverseErrorFunction.test_literal_values`
    #[test]
    fn test_erfinv_literal_values() {
        let y = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
        let actual = y.map(crate::erfinv);
        let expected = [
            0.0,
            0.08885599049425769,
            0.1791434546212917,
            0.2724627147267543,
            0.37080715859355795,
            0.4769362762044699,
            0.5951160814499948,
            0.7328690779592167,
            0.9061938024368233,
            1.1630871536766743,
        ];
        crate::np_assert_allclose!(&actual, &expected, atol = 1e-15);
    }

    #[test]
    fn test_erfcinv() {
        xsref::test("erfcinv", "d-d", |x| crate::erfcinv(x[0]));
    }
}
