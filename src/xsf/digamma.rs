mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait DigammaArg: sealed::Sealed {
    fn digamma(self) -> Self;
}

impl DigammaArg for f64 {
    #[inline(always)]
    fn digamma(self) -> Self {
        unsafe { crate::ffi::xsf::digamma(self) }
    }
}

impl DigammaArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn digamma(self) -> Self {
        unsafe { crate::ffi::xsf::digamma_1(self) }
    }
}

/// Digamma function for real or complex input
///
/// Corresponds to [`scipy.special.digamma`][digamma] in scipy
///
/// See also [`polygamma`](crate::polygamma).
///
/// [digamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.digamma.html
#[doc(alias = "psi")]
#[inline]
pub fn digamma<T: DigammaArg>(x: T) -> T {
    x.digamma()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_digamma_f64() {
        xsref::test("digamma", "d-d", |x| crate::digamma(x[0]));
    }

    #[test]
    fn test_digamma_c64() {
        xsref::test("digamma", "cd-cd", |x| crate::digamma(c64(x[0], x[1])));
    }
}
