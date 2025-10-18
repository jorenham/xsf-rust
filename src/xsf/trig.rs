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
    fn sinpi(self) -> Self {
        unsafe { crate::ffi::xsf::sinpi(self) }
    }

    #[inline(always)]
    fn cospi(self) -> Self {
        unsafe { crate::ffi::xsf::cospi(self) }
    }
}

impl TrigArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn sinpi(self) -> Self {
        unsafe { crate::ffi::xsf::sinpi_1(self) }
    }

    #[inline(always)]
    fn cospi(self) -> Self {
        unsafe { crate::ffi::xsf::cospi_1(self) }
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
    use num_complex::c64;

    #[test]
    fn test_sinpi_f64() {
        crate::xsref::test("sinpi", "d-d", |x| crate::sinpi(x[0]));
    }

    #[test]
    fn test_sinpi_c64() {
        crate::xsref::test("sinpi", "cd-cd", |x| crate::sinpi(c64(x[0], x[1])));
    }

    #[test]
    fn test_cospi_f64() {
        crate::xsref::test("cospi", "d-d", |x| crate::cospi(x[0]));
    }

    #[test]
    fn test_cospi_c64() {
        crate::xsref::test("cospi", "cd-cd", |x| crate::cospi(c64(x[0], x[1])));
    }

    #[test]
    fn test_sindg() {
        crate::xsref::test("sindg", "d-d", |x| crate::sindg(x[0]));
    }

    #[test]
    fn test_cosdg() {
        crate::xsref::test("cosdg", "d-d", |x| crate::cosdg(x[0]));
    }

    #[test]
    fn test_tandg() {
        crate::xsref::test("tandg", "d-d", |x| crate::tandg(x[0]));
    }

    #[test]
    fn test_cotdg() {
        crate::xsref::test("cotdg", "d-d", |x| crate::cotdg(x[0]));
    }

    #[test]
    fn test_cosm1() {
        crate::xsref::test("cosm1", "d-d", |x| crate::cosm1(x[0]));
    }

    #[test]
    fn test_radian() {
        crate::xsref::test("radian", "d_d_d-d", |x| crate::radian(x[0], x[1], x[2]));
    }
}
