//! Translated into pure Rust from `xsf/fp_error_metrics.h` (xsf v0.2.2).

use num_complex::Complex;

/// `std::numeric_limits<double>::denorm_min()`
const DENORM_MIN: f64 = f64::from_bits(1);

// max_float * 2**-(mantissa_bits + 1) = ulp(max_float)
const ULP_MAX: f64 = f64::EPSILON / 2.0 * f64::MAX;

#[allow(clippy::float_cmp)]
fn extended_absolute_error_f64(actual: f64, desired: f64) -> f64 {
    if actual == desired || (actual.is_nan() && desired.is_nan()) {
        return 0.0;
    }
    if desired.is_nan() || actual.is_nan() {
        /* If expected nan but got non-NaN or expected non-NaN but got NaN
         * we consider this to be an infinite error. */
        return f64::INFINITY;
    }
    if actual.is_infinite() {
        /* We don't want to penalize early overflow too harshly, so instead
         * compare with the mythical value nextafter(max_float). */
        let sgn = actual.signum();
        return ((sgn * f64::MAX - desired) + sgn * ULP_MAX).abs();
    }
    if desired.is_infinite() {
        let sgn = desired.signum();
        return ((sgn * f64::MAX - actual) + sgn * ULP_MAX).abs();
    }
    (actual - desired).abs()
}

fn extended_absolute_error_c64(actual: Complex<f64>, desired: Complex<f64>) -> f64 {
    extended_absolute_error_f64(actual.re, desired.re)
        .hypot(extended_absolute_error_f64(actual.im, desired.im))
}

fn extended_relative_error_f64(actual: f64, desired: f64) -> f64 {
    let abs_error = extended_absolute_error_f64(actual, desired);
    let mut abs_desired = desired.abs();
    if desired == 0.0 {
        /* If the desired result is 0.0, normalize by smallest subnormal instead
         * of zero. */
        abs_desired = DENORM_MIN;
    } else if desired.is_infinite() {
        abs_desired = f64::MAX;
    } else if desired.is_nan() {
        /* This ensures extended_relative_error(nan, nan) = 0 but
         * extended_relative_error(x0, x1) is infinite if one but not both of
         * x0 and x1 equals NaN */
        abs_desired = 1.0;
    }
    abs_error / abs_desired
}

fn extended_relative_error_c64(actual: Complex<f64>, mut desired: Complex<f64>) -> f64 {
    let abs_error = extended_absolute_error_c64(actual, desired);

    if desired.re == 0.0 {
        desired.re = DENORM_MIN.copysign(desired.re);
    } else if desired.re.is_infinite() {
        desired.re = f64::MAX.copysign(desired.re);
    } else if desired.re.is_nan() {
        /* In this case, the value used for desired doesn't matter. If desired.real() is NaN
         * but actual.real() isn't NaN, then the extended_absolute_error will be inf already
         * anyway. */
        desired.re = 1.0;
    }

    if desired.im == 0.0 {
        desired.im = DENORM_MIN.copysign(desired.im);
    } else if desired.im.is_infinite() {
        desired.im = f64::MAX.copysign(desired.im);
    } else if desired.im.is_nan() {
        /* In this case, the value used for desired doesn't matter. If desired.imag() is NaN
         * but actual.imag() isn't NaN, then the extended_absolute_error will be inf already
         * anyway. */
        desired.im = 1.0;
    }

    if !desired.is_infinite() && desired.norm().is_infinite() {
        /* Rescale to avoid overflow */
        return (abs_error / 2.0) / (desired / 2.0).norm();
    }

    abs_error / desired.norm()
}

pub trait ExtendedErrorArg: crate::sealed::Sealed {
    fn xsf_extended_absolute_error(self, other: Self) -> f64;
    fn xsf_extended_relative_error(self, other: Self) -> f64;
}

impl ExtendedErrorArg for f64 {
    #[inline]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        extended_absolute_error_f64(self, other)
    }

    #[inline]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        extended_relative_error_f64(self, other)
    }
}

impl ExtendedErrorArg for Complex<f64> {
    #[inline]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        extended_absolute_error_c64(self, other)
    }

    #[inline]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        extended_relative_error_c64(self, other)
    }
}

/// Extended absolute error metric between two `f64` or `Complex<f64>` values
#[inline]
pub fn extended_absolute_error<T: ExtendedErrorArg>(actual: T, expected: T) -> f64 {
    actual.xsf_extended_absolute_error(expected)
}

/// Extended relative error metric between two `f64` or `Complex<f64>` values
#[inline]
pub fn extended_relative_error<T: ExtendedErrorArg>(actual: T, expected: T) -> f64 {
    actual.xsf_extended_relative_error(expected)
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_extended_absolute_error_f64() {
        assert_eq!(crate::extended_absolute_error(0.0, 0.0), 0.0);
        assert_eq!(crate::extended_absolute_error(1.0, 0.0), 1.0);
        assert_eq!(crate::extended_absolute_error(1.0, 2.0), 1.0);
        assert_eq!(crate::extended_absolute_error(2.0, 1.0), 1.0);
        assert_eq!(crate::extended_absolute_error(3.0, 1.0), 2.0);
    }

    #[test]
    fn test_extended_absolute_error_c64() {
        assert_eq!(
            crate::extended_absolute_error(c64(1.0, 1.0), c64(1.0, 1.0)),
            0.0
        );
        assert_eq!(
            crate::extended_absolute_error(c64(0.0, 0.0), c64(3.0, 4.0)),
            5.0
        );
    }

    #[test]
    fn test_extended_relative_error_f64() {
        assert_eq!(crate::extended_relative_error(0.0, 0.0), 0.0);
        assert_eq!(crate::extended_relative_error(1.0, 0.0), f64::INFINITY);
        assert_eq!(crate::extended_relative_error(1.0, 2.0), 0.5);
        assert_eq!(crate::extended_relative_error(2.0, 1.0), 1.0);
        assert_eq!(crate::extended_relative_error(3.0, 1.0), 2.0);
    }

    #[test]
    fn test_extended_relative_error_c64() {
        assert_eq!(
            crate::extended_relative_error(c64(1.0, 1.0), c64(1.0, 1.0)),
            0.0
        );
        assert_eq!(
            crate::extended_relative_error(c64(0.0, 0.0), c64(3.0, 4.0)),
            1.0
        );
    }
}
