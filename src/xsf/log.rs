mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait LogArg: sealed::Sealed {
    fn xsf_log1p(self) -> Self;
    fn xsf_xlogy(self, x: Self) -> Self;
    fn xsf_xlog1py(self, x: Self) -> Self;
}

impl LogArg for f64 {
    #[inline]
    fn xsf_log1p(self) -> Self {
        unsafe { crate::ffi::xsf::log1p(self) }
    }

    #[inline]
    fn xsf_xlogy(self, x: Self) -> Self {
        unsafe { crate::ffi::xsf::xlogy(x, self) }
    }

    #[inline]
    fn xsf_xlog1py(self, x: Self) -> Self {
        unsafe { crate::ffi::xsf::xlog1py(x, self) }
    }
}

impl LogArg for num_complex::Complex<f64> {
    #[inline]
    fn xsf_log1p(self) -> Self {
        unsafe { crate::ffi::xsf::log1p_1(self) }
    }

    #[inline]
    fn xsf_xlogy(self, x: Self) -> Self {
        unsafe { crate::ffi::xsf::xlogy_1(x, self) }
    }

    #[inline]
    fn xsf_xlog1py(self, x: Self) -> Self {
        unsafe { crate::ffi::xsf::xlog1py_1(x, self) }
    }
}

/// `log(z + 1)` for real or complex input
#[doc(alias = "ln_1p", alias = "log_1p")]
#[must_use]
#[inline]
pub fn log1p<T: LogArg>(z: T) -> T {
    z.xsf_log1p()
}

/// Compute `log(1 + x) - x` for real input
#[must_use]
#[inline]
pub fn log1pmx(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::log1pmx(x) }
}

/// Compute `x * log(y)` for real or complex input
#[doc(alias = "x_ln_y", alias = "x_log_y")]
#[must_use]
#[inline]
pub fn xlogy<T: LogArg>(x: T, y: T) -> T {
    y.xsf_xlogy(x)
}

/// Compute `x * log(1 + y)` for real or complex input
#[doc(alias = "x_ln_1py", alias = "x_log_1py")]
#[must_use]
#[inline]
pub fn xlog1py<T: LogArg>(x: T, y: T) -> T {
    y.xsf_xlog1py(x)
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_log1p_f64() {
        xsref::test("log1p", "d-d", |x| crate::log1p(x[0]));
    }

    #[test]
    fn test_log1p_c64() {
        xsref::test("log1p", "cd-cd", |x| crate::log1p(c64(x[0], x[1])));
    }

    #[test]
    fn test_log1pmx_f64() {
        xsref::test("log1pmx", "d-d", |x| crate::log1pmx(x[0]));
    }

    #[test]
    fn test_xlogy_f64() {
        xsref::test("xlogy", "d_d-d", |x| crate::xlogy(x[0], x[1]));
    }

    #[test]
    fn test_xlogy_c64() {
        xsref::test("xlogy", "cd_cd-cd", |x| {
            crate::xlogy(c64(x[0], x[1]), c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_xlog1py_f64() {
        xsref::test("xlog1py", "d_d-d", |x| crate::xlog1py(x[0], x[1]));
    }

    #[test]
    fn test_xlog1py_c64() {
        xsref::test("xlog1py", "cd_cd-cd", |x| {
            crate::xlog1py(c64(x[0], x[1]), c64(x[2], x[3]))
        });
    }
}
