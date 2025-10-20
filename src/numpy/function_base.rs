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
pub fn sinc(x: f64) -> f64 {
    if x == 0.0 {
        1.0
    } else {
        let pi_x = PI * x;
        pi_x.sin() / pi_x
    }
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    #[test]
    fn test_sinc() {
        // based on `numpy.lib.tests.test_function_base.TestSinc.test_simple`
        let w = (-1..=1)
            .map(|i| i as f64 * 0.01)
            .map(crate::sinc)
            .collect::<Vec<f64>>();
        np_assert_allclose!(
            w.clone(),
            w.iter().copied().rev().collect::<Vec<_>>(),
            atol = 1e-15
        );
    }
}
