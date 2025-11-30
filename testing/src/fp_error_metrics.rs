use num_complex::Complex;

const DENORM_MIN: f64 = f64::from_bits(1);
const ULP_AT_MAX: f64 = (f64::EPSILON * f64::MAX) / 2.0;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExtendedErrorArg: Copy + sealed::Sealed {
    fn xsf_extended_absolute_error(self, other: Self) -> f64;
    fn xsf_extended_relative_error(self, other: Self) -> f64;
    fn xsf_magnitude(self) -> f64;
    fn xsf_is_nan(self) -> bool;
}

impl ExtendedErrorArg for f32 {
    #[inline]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        extended_absolute_error_scalar(self as f64, other as f64)
    }

    #[inline]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        extended_relative_error_scalar(self as f64, other as f64)
    }

    #[inline]
    fn xsf_magnitude(self) -> f64 {
        self.abs() as f64
    }

    #[inline]
    fn xsf_is_nan(self) -> bool {
        self.is_nan()
    }
}

impl ExtendedErrorArg for f64 {
    #[inline]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        extended_absolute_error_scalar(self, other)
    }

    #[inline]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        extended_relative_error_scalar(self, other)
    }

    #[inline]
    fn xsf_magnitude(self) -> f64 {
        self.abs()
    }

    #[inline]
    fn xsf_is_nan(self) -> bool {
        self.is_nan()
    }
}

impl ExtendedErrorArg for Complex<f64> {
    #[inline]
    fn xsf_extended_absolute_error(self, other: Self) -> f64 {
        extended_absolute_error_complex(self, other)
    }

    #[inline]
    fn xsf_extended_relative_error(self, other: Self) -> f64 {
        extended_relative_error_complex(self, other)
    }

    #[inline]
    fn xsf_magnitude(self) -> f64 {
        self.l1_norm()
    }

    #[inline]
    fn xsf_is_nan(self) -> bool {
        self.is_nan()
    }
}

#[inline]
fn extended_absolute_error_scalar(actual: f64, desired: f64) -> f64 {
    if actual == desired || (actual.is_nan() && desired.is_nan()) {
        return 0.0;
    }

    if desired.is_nan() || actual.is_nan() {
        return f64::INFINITY;
    }

    if actual.is_infinite() {
        let sgn = actual.signum();
        let max_float = f64::MAX;
        return ((sgn * max_float - desired) + sgn * ULP_AT_MAX).abs();
    }

    if desired.is_infinite() {
        let sgn = desired.signum();
        let max_float = f64::MAX;
        return ((sgn * max_float - actual) + sgn * ULP_AT_MAX).abs();
    }

    (actual - desired).abs()
}

#[inline]
fn extended_relative_error_scalar(actual: f64, desired: f64) -> f64 {
    let abs_error = extended_absolute_error_scalar(actual, desired);
    let mut abs_desired = desired.abs();

    if desired == 0.0 {
        abs_desired = DENORM_MIN;
    } else if desired.is_infinite() {
        abs_desired = f64::MAX;
    } else if desired.is_nan() {
        abs_desired = 1.0;
    }

    abs_error / abs_desired
}

#[inline]
fn extended_absolute_error_complex(actual: Complex<f64>, desired: Complex<f64>) -> f64 {
    let real_err = extended_absolute_error_scalar(actual.re, desired.re);
    let imag_err = extended_absolute_error_scalar(actual.im, desired.im);
    real_err.hypot(imag_err)
}

#[inline]
fn adjust_component(mut value: f64) -> f64 {
    if value == 0.0 {
        value = f64::copysign(DENORM_MIN, value);
    } else if value.is_infinite() {
        value = f64::copysign(f64::MAX, value);
    } else if value.is_nan() {
        value = 1.0;
    }

    value
}

#[inline]
fn extended_relative_error_complex(actual: Complex<f64>, desired: Complex<f64>) -> f64 {
    let abs_error = extended_absolute_error_complex(actual, desired);
    let mut desired = desired;

    desired.re = adjust_component(desired.re);
    desired.im = adjust_component(desired.im);

    let magnitude = desired.norm();

    if !desired.re.is_infinite() && !desired.im.is_infinite() && magnitude.is_infinite() {
        return (abs_error / 2.0) / ((desired / 2.0).norm());
    }

    abs_error / magnitude
}

pub fn extended_absolute_error<T: ExtendedErrorArg>(actual: T, expected: T) -> f64 {
    actual.xsf_extended_absolute_error(expected)
}

pub fn extended_relative_error<T: ExtendedErrorArg>(actual: T, expected: T) -> f64 {
    actual.xsf_extended_relative_error(expected)
}
