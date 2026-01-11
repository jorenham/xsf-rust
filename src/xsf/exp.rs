mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExpArg: sealed::Sealed {
    fn expm1(self) -> Self;
}

impl ExpArg for f64 {
    #[inline]
    fn expm1(self) -> Self {
        unsafe { crate::ffi::xsf::expm1(self) }
    }
}

impl ExpArg for num_complex::Complex<f64> {
    #[inline]
    fn expm1(self) -> Self {
        unsafe { crate::ffi::xsf::expm1_1(self) }
    }
}

/// $e^x - 1$ for real or complex input
///
/// Corresponds to [`scipy.special.expm1`][expm1]
///
/// [expm1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expm1.html
///
/// # See Also
/// - [`exp2`] for $2^x$
/// - [`exp10`] for $10^x$
#[doc(alias = "exp_m1")]
#[must_use]
#[inline]
pub fn expm1<T: ExpArg>(z: T) -> T {
    z.expm1()
}

/// $2^x$
///
/// Corresponds to [`scipy.special.exp2`][exp2]
///
/// [exp2]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp2.html
///
/// # See Also
/// - [`expm1`] for $e^x - 1$
/// - [`exp10`] for $10^x$
#[must_use]
#[inline]
pub fn exp2(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::exp2(x) }
}

/// $10^x$
///
/// Corresponds to [`scipy.special.exp10`][exp10]
///
/// [exp10]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp10.html
///
/// # See Also
/// - [`expm1`] for $e^x - 1$
/// - [`exp2`] for $2^x$
#[must_use]
#[inline]
pub fn exp10(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::exp10(x) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_expm1_f64() {
        xsref::test("expm1", "d-d", |x| crate::expm1(x[0]));
    }

    #[test]
    fn test_expm1_c64() {
        xsref::test("expm1", "cd-cd", |x| crate::expm1(c64(x[0], x[1])));
    }

    #[test]
    fn test_exp2_f64() {
        xsref::test("exp2", "d-d", |x| crate::exp2(x[0]));
    }

    #[test]
    fn test_exp10_f64() {
        xsref::test("exp10", "d-d", |x| crate::exp10(x[0]));
    }
}
