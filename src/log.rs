use crate::bindings;
use crate::bindings::xsf_impl;
use num_complex::Complex;

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
    #[inline(always)]
    fn xsf_log1p(self) -> f64 {
        unsafe { bindings::log1p(self) }
    }
    #[inline(always)]
    fn xsf_xlogy(self, x: Self) -> Self {
        unsafe { bindings::xlogy(x, self) }
    }
    #[inline(always)]
    fn xsf_xlog1py(self, x: Self) -> Self {
        unsafe { bindings::xlog1py(x, self) }
    }
}

impl LogArg for Complex<f64> {
    #[inline(always)]
    fn xsf_log1p(self) -> Complex<f64> {
        unsafe { bindings::log1p_1(self.into()) }.into()
    }
    #[inline(always)]
    fn xsf_xlogy(self, x: Self) -> Self {
        unsafe { bindings::xlogy_1(x.into(), self.into()) }.into()
    }
    #[inline(always)]
    fn xsf_xlog1py(self, x: Self) -> Self {
        unsafe { bindings::xlog1py_1(x.into(), self.into()) }.into()
    }
}

/// `log(z + 1)` for real or complex input
pub fn log1p<T: LogArg>(z: T) -> T {
    z.xsf_log1p()
}

xsf_impl!(log1pmx, (x: f64), "Compute `log(1 + x) - x` for real input");

/// Compute `x * log(y)` for real or complex input
pub fn xlogy<T: LogArg>(x: T, y: T) -> T {
    y.xsf_xlogy(x)
}
/// Compute `x * log(1 + y)` for real or complex input
pub fn xlog1py<T: LogArg>(x: T, y: T) -> T {
    y.xsf_xlog1py(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    #[test]
    fn test_log1p_f64() {
        testing::test::<f64, _>("log1p", "d-d", |x: &[f64]| log1p(x[0]));
    }

    #[test]
    fn test_log1p_c64() {
        testing::test::<Complex<f64>, _>("log1p", "cd-cd", |x: &[f64]| log1p(c64(x[0], x[1])));
    }

    #[test]
    fn test_log1pmx_f64() {
        testing::test::<f64, _>("log1pmx", "d-d", |x: &[f64]| log1pmx(x[0]));
    }

    #[test]
    fn test_xlogy_f64() {
        testing::test::<f64, _>("xlogy", "d_d-d", |x: &[f64]| xlogy(x[0], x[1]));
    }

    #[test]
    fn test_xlogy_c64() {
        testing::test::<Complex<f64>, _>("xlogy", "cd_cd-cd", |x: &[f64]| {
            xlogy(c64(x[0], x[1]), c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_xlog1py_f64() {
        testing::test::<f64, _>("xlog1py", "d_d-d", |x: &[f64]| xlog1py(x[0], x[1]));
    }

    #[test]
    fn test_xlog1py_c64() {
        testing::test::<Complex<f64>, _>("xlog1py", "cd_cd-cd", |x: &[f64]| {
            xlog1py(c64(x[0], x[1]), c64(x[2], x[3]))
        });
    }
}
