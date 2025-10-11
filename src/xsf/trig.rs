use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait TrigArg: sealed::Sealed {
    fn sinpi(self) -> Self;
    fn cospi(self) -> Self;
}

impl TrigArg for f64 {
    #[inline(always)]
    fn sinpi(self) -> f64 {
        unsafe { crate::ffi::xsf::sinpi(self) }
    }

    #[inline(always)]
    fn cospi(self) -> Self {
        unsafe { crate::ffi::xsf::cospi(self) }
    }
}

impl TrigArg for Complex<f64> {
    #[inline(always)]
    fn sinpi(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::sinpi_1(self.into()) }.into()
    }

    #[inline(always)]
    fn cospi(self) -> Self {
        unsafe { crate::ffi::xsf::cospi_1(self.into()) }.into()
    }
}

/// Compute `sin(pi*z)` for real or complex `z`
pub fn sinpi<T: TrigArg>(z: T) -> T {
    z.sinpi()
}

/// Compute `cos(pi*z)` for real or complex `z`
pub fn cospi<T: TrigArg>(z: T) -> T {
    z.cospi()
}

/// Circular sine of angle in degrees
pub fn sindg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::sindg(x) }
}

/// Circular cosine of angle in degrees
pub fn cosdg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cosdg(x) }
}

/// Circular tangent of argument in degrees
pub fn tandg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::tandg(x) }
}

/// Circular cotangent of argument in degrees
pub fn cotdg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cotdg(x) }
}

/// Compute `cos(x) - 1`
pub fn cosm1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cosm1(x) }
}

/// Degrees, minutes, seconds to radians
pub fn radian(d: f64, m: f64, s: f64) -> f64 {
    unsafe { crate::ffi::xsf::radian(d, m, s) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    #[test]
    fn test_sinpi_f64() {
        testing::test::<f64, _>("sinpi", "d-d", |x: &[f64]| sinpi(x[0]));
    }

    #[test]
    fn test_sinpi_c64() {
        testing::test::<Complex<f64>, _>("sinpi", "cd-cd", |x: &[f64]| sinpi(c64(x[0], x[1])));
    }

    #[test]
    fn test_cospi_f64() {
        testing::test::<f64, _>("cospi", "d-d", |x: &[f64]| cospi(x[0]));
    }

    #[test]
    fn test_cospi_c64() {
        testing::test::<Complex<f64>, _>("cospi", "cd-cd", |x: &[f64]| cospi(c64(x[0], x[1])));
    }

    #[test]
    fn test_sindg() {
        testing::test::<f64, _>("sindg", "d-d", |x: &[f64]| sindg(x[0]));
    }

    #[test]
    fn test_cosdg() {
        testing::test::<f64, _>("cosdg", "d-d", |x: &[f64]| cosdg(x[0]));
    }

    #[test]
    fn test_tandg() {
        testing::test::<f64, _>("tandg", "d-d", |x: &[f64]| tandg(x[0]));
    }

    #[test]
    fn test_cotdg() {
        testing::test::<f64, _>("cotdg", "d-d", |x: &[f64]| cotdg(x[0]));
    }

    #[test]
    fn test_cosm1() {
        testing::test::<f64, _>("cosm1", "d-d", |x: &[f64]| cosm1(x[0]));
    }

    #[test]
    fn test_radian() {
        testing::test::<f64, _>("radian", "d_d_d-d", |x: &[f64]| radian(x[0], x[1], x[2]));
    }
}
