use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExpIntArg: sealed::Sealed {
    fn expi(self) -> Self;
    fn exp1(self) -> Self;
}

impl ExpIntArg for f64 {
    #[inline(always)]
    fn expi(self) -> f64 {
        unsafe { crate::ffi::xsf::expi(self) }
    }
    #[inline(always)]
    fn exp1(self) -> f64 {
        unsafe { crate::ffi::xsf::exp1(self) }
    }
}

impl ExpIntArg for Complex<f64> {
    #[inline(always)]
    fn expi(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::expi_1(self.into()) }.into()
    }
    #[inline(always)]
    fn exp1(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::exp1_1(self.into()) }.into()
    }
}

/// Exponential integral E₁ for real or complex input
///
/// For complex z ≠ 0 the exponential integral can be defined as
///
/// E₁(z) = ∫[z,∞) e⁻ᵗ / t dt,
///
/// where the path of the integral does not cross the negative real axis or pass through the origin.
///
/// ## See also:
/// - [`expi`]: exponential integral Ei
/// - [`expn`](fn.expn.html): generalization of E₁
pub fn exp1<T: ExpIntArg>(z: T) -> T {
    z.exp1()
}

/// Exponential integral Ei for real or complex input
///
/// For real x, the exponential integral is defined as
///
/// Ei(x) = ∫(−∞,x] eᵗ / t dt.
///
/// For x > 0 the integral is understood as a Cauchy principal value.
///
/// It is extended to the complex plane by analytic continuation of the function on the interval
/// (0, ∞). The complex variant has a branch cut on the negative real axis.
///
/// ## See also:
/// - [`exp1`]: exponential integral E₁
/// - [`expn`](fn.expn.html): generalization of E₁
pub fn expi<T: ExpIntArg>(z: T) -> T {
    z.expi()
}

/// Scaled version of the exponential integral E₁ for real input
///
/// This function computes F(x), where F is the factor remaining in E₁(x)
/// when e⁻ˣ/x is factored out. That is:
///
/// E₁(x) = x⁻¹ e⁻ˣ F(x)
///
/// or
///
/// F(x) = x eˣ E₁(x)
///
/// The function is defined for real x ≥ 0. For x < 0, NaN is returned.
///
/// F has the properties:
///
/// - F(0) = 0
/// - F(x) is increasing on [0, ∞).
/// - F(x) = 1 in the limit as x → ∞.
pub fn scaled_exp1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::scaled_exp1(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    // exp1

    #[test]
    fn test_exp1_f64() {
        testing::test::<f64, _>("exp1", "d-d", |x: &[f64]| exp1(x[0]));
    }

    #[test]
    fn test_exp1_c64() {
        testing::test::<Complex<f64>, _>("exp1", "cd-cd", |x: &[f64]| exp1(c64(x[0], x[1])));
    }

    // expi

    #[test]
    fn test_expi_f64() {
        testing::test::<f64, _>("expi", "d-d", |x: &[f64]| expi(x[0]));
    }

    #[test]
    fn test_expi_c64() {
        testing::test::<Complex<f64>, _>("expi", "cd-cd", |x: &[f64]| expi(c64(x[0], x[1])));
    }

    // scaled_exp1

    #[test]
    fn test_scaled_exp1_f64() {
        testing::test::<f64, _>("scaled_exp1", "d-d", |x: &[f64]| scaled_exp1(x[0]));
    }
}
