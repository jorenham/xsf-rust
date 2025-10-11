use crate::ffi;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExtendedErrorArg: sealed::Sealed {
    fn xsf_extended_absolute_error(self, other: Self) -> f64;
    fn xsf_extended_relative_error(self, other: Self) -> f64;
}

impl ExtendedErrorArg for f64 {
    #[inline(always)]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        unsafe { ffi::extended_absolute_error(self, other) }
    }
    #[inline(always)]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        unsafe { ffi::extended_relative_error(self, other) }
    }
}

impl ExtendedErrorArg for Complex<f64> {
    #[inline(always)]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        unsafe { ffi::extended_absolute_error_1(self.into(), other.into()) }
    }
    #[inline(always)]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        unsafe { ffi::extended_relative_error_1(self.into(), other.into()) }
    }
}

/// Extended absolute error metric between two `f64` or `Complex<f64>` values
pub fn extended_absolute_error<T: ExtendedErrorArg>(actual: T, expected: T) -> f64 {
    actual.xsf_extended_absolute_error(expected)
}

/// Extended relative error metric between two `f64` or `Complex<f64>` values
pub fn extended_relative_error<T: ExtendedErrorArg>(actual: T, expected: T) -> f64 {
    actual.xsf_extended_relative_error(expected)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::c64;

    #[test]
    fn test_extended_absolute_error_f64() {
        assert_eq!(extended_absolute_error(0.0, 0.0), 0.0);
        assert_eq!(extended_absolute_error(1.0, 0.0), 1.0);
        assert_eq!(extended_absolute_error(1.0, 2.0), 1.0);
        assert_eq!(extended_absolute_error(2.0, 1.0), 1.0);
        assert_eq!(extended_absolute_error(3.0, 1.0), 2.0);
    }

    #[test]
    fn test_extended_absolute_error_c64() {
        assert_eq!(extended_absolute_error(c64(1.0, 1.0), c64(1.0, 1.0)), 0.0);
        assert_eq!(extended_absolute_error(c64(0.0, 0.0), c64(3.0, 4.0)), 5.0);
    }

    #[test]
    fn test_extended_relative_error_f64() {
        assert_eq!(extended_relative_error(0.0, 0.0), 0.0);
        assert_eq!(extended_relative_error(1.0, 0.0), f64::INFINITY);
        assert_eq!(extended_relative_error(1.0, 2.0), 0.5);
        assert_eq!(extended_relative_error(2.0, 1.0), 1.0);
        assert_eq!(extended_relative_error(3.0, 1.0), 2.0);
    }

    #[test]
    fn test_extended_relative_error_c64() {
        assert_eq!(extended_relative_error(c64(1.0, 1.0), c64(1.0, 1.0)), 0.0);
        assert_eq!(extended_relative_error(c64(0.0, 0.0), c64(3.0, 4.0)), 1.0);
    }
}
