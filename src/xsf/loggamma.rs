mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait LogGammaArg: sealed::Sealed {
    fn xsf_loggamma(self) -> Self;
    fn xsf_rgamma(self) -> Self;
}

impl LogGammaArg for f64 {
    #[inline]
    fn xsf_loggamma(self) -> f64 {
        unsafe { crate::ffi::xsf::loggamma(self) }
    }

    #[inline]
    fn xsf_rgamma(self) -> f64 {
        unsafe { crate::ffi::xsf::rgamma(self) }
    }
}

impl LogGammaArg for num_complex::Complex<f64> {
    #[inline]
    fn xsf_loggamma(self) -> Self {
        unsafe { crate::ffi::xsf::loggamma_1(self) }
    }

    #[inline]
    fn xsf_rgamma(self) -> Self {
        unsafe { crate::ffi::xsf::rgamma_1(self) }
    }
}

/// Principal branch of the logarithm of `gamma(z)`
#[doc(alias = "lgamma", alias = "ln_gamma", alias = "log_gamma")]
#[inline]
pub fn loggamma<T: LogGammaArg>(z: T) -> T {
    z.xsf_loggamma()
}

/// Reciprocal Gamma function `1 / gamma(z)`
#[inline]
pub fn rgamma<T: LogGammaArg>(z: T) -> T {
    z.xsf_rgamma()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_loggamma_f64() {
        xsref::test("loggamma", "d-d", |x| crate::loggamma(x[0]));
    }

    #[test]
    fn test_loggamma_c64() {
        xsref::test("loggamma", "cd-cd", |x| crate::loggamma(c64(x[0], x[1])));
    }

    #[test]
    fn test_rgamma_f64() {
        xsref::test("rgamma", "d-d", |x| crate::rgamma(x[0]));
    }

    #[test]
    fn test_rgamma_c64() {
        xsref::test("rgamma", "cd-cd", |x| crate::rgamma(c64(x[0], x[1])));
    }
}
