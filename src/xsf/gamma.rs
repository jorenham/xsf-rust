mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait GammaArg: sealed::Sealed {
    fn xsf_gamma(self) -> Self;
}

impl GammaArg for f64 {
    #[inline(always)]
    fn xsf_gamma(self) -> Self {
        unsafe { crate::ffi::xsf::gamma(self) }
    }
}

impl GammaArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn xsf_gamma(self) -> Self {
        unsafe { crate::ffi::xsf::gamma_1(self) }
    }
}

/// Gamma function for real or complex input
pub fn gamma<T: GammaArg>(z: T) -> T {
    z.xsf_gamma()
}

/// Regularized lower incomplete gamma function
#[doc(alias = "inc_gamma")]
#[doc(alias = "gamma_lr")]
pub fn gammainc(a: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammainc(a, x) }
}

/// Regularized upper incomplete gamma function
#[doc(alias = "gamma_ur")]
pub fn gammaincc(a: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaincc(a, x) }
}

/// Inverse of [`gammainc`]
pub fn gammaincinv(a: f64, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaincinv(a, p) }
}

/// Inverse of [`gammaincc`]
pub fn gammainccinv(a: f64, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammainccinv(a, p) }
}

/// Logarithm of the absolute value of the gamma function
pub fn gammaln(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaln(x) }
}

/// Sign of the Gamma function
pub fn gammasgn(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammasgn(x) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_gamma_f64() {
        xsref::test("gamma", "d-d", |x| crate::gamma(x[0]));
    }

    #[test]
    fn test_gamma_c64() {
        xsref::test("gamma", "cd-cd", |x| crate::gamma(c64(x[0], x[1])));
    }

    #[test]
    fn test_gammainc() {
        xsref::test("gammainc", "d_d-d", |x| crate::gammainc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincc() {
        xsref::test("gammaincc", "d_d-d", |x| crate::gammaincc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincinv() {
        xsref::test("gammaincinv", "d_d-d", |x| crate::gammaincinv(x[0], x[1]));
    }

    #[test]
    fn test_gammainccinv() {
        xsref::test("gammainccinv", "d_d-d", |x| crate::gammainccinv(x[0], x[1]));
    }

    #[test]
    fn test_gammaln() {
        xsref::test("gammaln", "d-d", |x| crate::gammaln(x[0]));
    }

    #[test]
    fn test_gammasgn() {
        xsref::test("gammasgn", "d-d", |x| crate::gammasgn(x[0]));
    }
}
