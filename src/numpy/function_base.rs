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
        // Unlike in NumPy, sin(pi x) is computed with `sinpi`, which reduces x modulo 2 first.
        // NumPy computes sin(y) / y with y = pi x, whose rounding error is amplified for large x
        // and near the zeros at the nonzero integers.
        //
        // E.g. sinc(-9781324477.118902) was -1.1875257554868476e-11 instead of
        // -1.1875312717547689e-11.
        crate::sinpi(x) / (PI * x)
    }
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    #[test]
    fn test_sinc() {
        // based on `numpy.lib.tests.test_function_base.TestSinc.test_simple`
        let w = (-1..=1)
            .map(|i| f64::from(i) * 0.01)
            .map(crate::sinc)
            .collect::<Vec<f64>>();
        np_assert_allclose!(
            w.clone(),
            w.iter().copied().rev().collect::<Vec<_>>(),
            atol = 1e-15
        );
    }

    #[test]
    fn test_sinc_large_or_near_int() {
        // mpmath reference values
        let x = [
            -9_781_324_477.118_902,
            1.000_000_000_1,
            2.999_999_75,
            1e15 + 0.5,
        ];
        let y = [
            -1.187_531_271_754_768_9e-11,
            -1.000_000_082_640_371e-10,
            8.333_334_021_540_307e-8,
            3.183_098_861_837_905e-16,
        ];
        np_assert_allclose!(x.map(crate::sinc), y, rtol = 1e-15);
    }
}
