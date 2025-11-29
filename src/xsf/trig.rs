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
    #[inline]
    fn sinpi(self) -> Self {
        unsafe { crate::ffi::xsf::sinpi(self) }
    }

    #[inline]
    fn cospi(self) -> Self {
        unsafe { crate::ffi::xsf::cospi(self) }
    }
}

impl TrigArg for num_complex::Complex<f64> {
    #[inline]
    fn sinpi(self) -> Self {
        unsafe { crate::ffi::xsf::sinpi_1(self) }
    }

    #[inline]
    fn cospi(self) -> Self {
        unsafe { crate::ffi::xsf::cospi_1(self) }
    }
}

/// Compute `sin(pi*z)` for real or complex `z`
#[must_use]
#[inline]
pub fn sinpi<T: TrigArg>(z: T) -> T {
    z.sinpi()
}

/// Compute `cos(pi*z)` for real or complex `z`
#[must_use]
#[inline]
pub fn cospi<T: TrigArg>(z: T) -> T {
    z.cospi()
}

/// Circular sine of angle in degrees
#[must_use]
#[inline]
pub fn sindg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::sindg(x) }
}

/// Circular cosine of angle in degrees
#[must_use]
#[inline]
pub fn cosdg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cosdg(x) }
}

/// Circular tangent of argument in degrees
#[must_use]
#[inline]
pub fn tandg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::tandg(x) }
}

/// Circular cotangent of argument in degrees
#[must_use]
#[inline]
pub fn cotdg(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cotdg(x) }
}

/// Compute `cos(x) - 1`
#[must_use]
#[inline]
pub fn cosm1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cosm1(x) }
}

/// Degrees, minutes, seconds to radians
#[must_use]
#[inline]
pub fn radian(d: f64, m: f64, s: f64) -> f64 {
    unsafe { crate::ffi::xsf::radian(d, m, s) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_sinpi_f64() {
        xsref::test("sinpi", "d-d", |x| crate::sinpi(x[0]));
    }

    #[test]
    fn test_sinpi_c64() {
        xsref::test("sinpi", "cd-cd", |x| crate::sinpi(c64(x[0], x[1])));
    }

    #[test]
    fn test_cospi_f64() {
        xsref::test("cospi", "d-d", |x| crate::cospi(x[0]));
    }

    #[test]
    fn test_cospi_c64() {
        xsref::test("cospi", "cd-cd", |x| crate::cospi(c64(x[0], x[1])));
    }

    #[test]
    fn test_sindg() {
        xsref::test("sindg", "d-d", |x| crate::sindg(x[0]));
    }

    #[test]
    fn test_cosdg() {
        xsref::test("cosdg", "d-d", |x| crate::cosdg(x[0]));
    }

    #[test]
    fn test_tandg() {
        xsref::test("tandg", "d-d", |x| crate::tandg(x[0]));
    }

    #[test]
    fn test_cotdg() {
        xsref::test("cotdg", "d-d", |x| crate::cotdg(x[0]));
    }

    #[test]
    fn test_cosm1() {
        xsref::test("cosm1", "d-d", |x| crate::cosm1(x[0]));
    }

    #[test]
    fn test_radian() {
        xsref::test("radian", "d_d_d-d", |x| crate::radian(x[0], x[1], x[2]));
    }
}
