//! <https://github.com/numpy/numpy/blob/4992eaf/numpy/lib/_function_base_impl.py>
use core::f64::consts::PI;

/// Normalized sinc function
///
/// The sinc function is equal to sin(π*x*)/(π*x*) for any argument *x* ≠ 0.
/// sinc(0) takes the limit value 1, making *sinc* not only everywhere continuous but also
/// infinitely differentiable.
///
/// See [NumPy's documentation](https://numpy.org/doc/stable/reference/generated/numpy.sinc.html)
/// for more details.
#[must_use]
#[inline]
pub fn sinc(x: f64) -> f64 {
    if x == 0.0 {
        1.0
    } else {
        // Unlike in NumPy, sin(pi x) is computed with `sinpi`, as the rounding error of y = pi x in
        // sin(y) / y is amplified for large x and near the nonzero integers, e.g.
        // sinc(3.0) was 4e-17 instead of 0.
        crate::sinpi(x) / (PI * x)
    }
}
