//! Translated from `scipy/special/_logsumexp.py`

const NINF: f64 = f64::NEG_INFINITY;

/// Return the maximum value in a slice of f64, propagating NaNs
#[inline]
fn fmax(xs: &[f64]) -> f64 {
    xs.iter().fold(NINF, |a, &b| if a >= b { a } else { b })
}

/// Compute the log of the sum of exponentials of input elements
///
/// Pure Rust translation of [`scipy.special.logsumexp`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html
///
/// # Definition
///
/// For a sequence of real values $\vec{\bm{x}}$, `logsumexp` is defined as
///
/// $$ \ln \sum_{x \in \vec{\bm{x}}} e^x $$
///
/// The implementation uses a numerically stable equivalent algorithm which factors out the maximum
/// value.
///
/// # Notes
///
/// Unlike the SciPy version, this implementation only supports 1-dimensional slices, and doesn't
/// support complex numbers at the moment.
///
/// # See also
/// - [`logaddexp`](crate::logaddexp)
/// - [`softmax`]
#[must_use]
#[allow(clippy::float_cmp, clippy::missing_panics_doc)]
pub fn logsumexp(xs: &[f64]) -> f64 {
    let x_max = fmax(xs);

    // fast-path for -/+inf and NaN
    if !x_max.is_finite() {
        return x_max;
    }

    // for precision, these terms are separated out of the main sum.
    let m: f64 = xs
        .iter()
        .fold(0.0, |acc, &x| if x == x_max { acc + 1.0 } else { acc });
    let xs = xs.iter().map(|&x| if x == x_max { NINF } else { x });

    // Shift, exponentiate, and sum
    let s = xs.map(|x| (x - x_max).exp()).sum::<f64>();

    // The log functions need positive arguments
    let s = if s < -1.0 { -s - 2.0 } else { s };

    // Take log and undo shift
    s.ln_1p() + m.ln() + x_max
}

/// Compute the softmax function
///
/// Pure Rust translation of [`scipy.special.softmax`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.softmax.html
///
/// # Definition
///
/// The formula for the softmax function $\sigma(\bm{x})$ for a vector
/// $\bm{x} = \begin{bmatrix} x_0 & x_1 & \cdots & x_{n-1} \end{bmatrix}$ is
///
/// $$
/// \sigma(\bm{x})_j = {e^{x_j} \over \sum_k e^{x_k}}
/// $$
///
/// The `softmax` function is the gradient of [`logsumexp`].
///
/// The implementation uses shifting to improve numerical stability.
///
/// # See also
/// - [`log_softmax`]
/// - [`logsumexp`](crate::logsumexp)
#[must_use]
pub fn softmax(xs: &[f64]) -> Vec<f64> {
    let x_max = fmax(xs);
    let exs0 = xs.iter().map(|&x| (x - x_max).exp()).collect::<Vec<f64>>();
    let s = exs0.iter().sum::<f64>().recip();
    exs0.iter().map(|ex| ex * s).collect()
}

/// Compute the logarithm of the [`softmax`] function
///
/// Pure Rust translation of [`scipy.special.log_softmax`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_softmax.html
///
/// # Notes
///
/// `log_softwmax` is more accurate than taking the log of [`softmax`] with inputs that make
/// [`softmax`] saturate.
///
/// # See also
/// - [`softmax`]
/// - [`logsumexp`](crate::logsumexp)
#[must_use]
pub fn log_softmax(xs: &[f64]) -> Vec<f64> {
    let x_max = fmax(xs);
    let x_max = if x_max.is_finite() { x_max } else { 0.0 };

    let xs0 = xs.iter().map(|&x| x - x_max).collect::<Vec<f64>>();
    let exs0 = xs0.iter().map(|&x| x.exp()).collect::<Vec<f64>>();

    let s = exs0.iter().sum::<f64>().ln();
    xs0.iter().map(|&x| x - s).collect()
}

#[cfg(test)]
mod tests {
    //! Translated from `scipy.special.tests.test_logsumexp` at
    //! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_logsumexp.py>

    use crate::np_assert_allclose;
    use core::f64::consts::{E, LN_2};

    /// Corresponds to `TestLogSumExp.test_logsumexp`.
    #[test]
    #[allow(clippy::float_cmp, clippy::cast_precision_loss)]
    fn test_logsumexp() {
        // default rtol of `xp_assert_close`
        let rtol: f64 = f64::EPSILON.sqrt() * 4.0;

        // Test with zero-size array
        assert_eq!(crate::logsumexp(&[]), f64::NEG_INFINITY);

        // Test whether logsumexp() function correctly handles large inputs.
        let xs = (0..200).map(f64::from).collect::<Vec<f64>>();
        let desired = xs.iter().copied().map(f64::exp).sum::<f64>().ln();
        np_assert_allclose!([crate::logsumexp(&xs)], [desired], rtol = rtol);

        // Now test with large numbers
        let xs = [1000.0, 1000.0];
        let desired = 1000.0 + LN_2;
        np_assert_allclose!([crate::logsumexp(&xs)], [desired], rtol = rtol);

        let (x, n) = (10_000.0, 1_000);
        let xs = vec![x; n];
        let desired = x + (n as f64).ln();
        np_assert_allclose!([crate::logsumexp(&xs)], [desired], rtol = rtol);

        let ys = vec![1e-40; 1_000_000];
        let xs = ys.iter().copied().map(f64::ln).collect::<Vec<f64>>();
        np_assert_allclose!(
            [crate::logsumexp(&xs).exp()],
            [ys.iter().sum::<f64>()],
            rtol = rtol
        );

        // Handling special values properly
        assert!(crate::logsumexp(&[f64::NAN]).is_nan());
        assert_eq!(crate::logsumexp(&[f64::INFINITY]), f64::INFINITY);
        assert_eq!(crate::logsumexp(&[f64::NEG_INFINITY]), f64::NEG_INFINITY);
        assert_eq!(
            crate::logsumexp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]),
            f64::NEG_INFINITY
        );
    }

    /// Corresponds to `TestSoftmax.test_softmax_fixtures`.
    #[test]
    fn test_softmax() {
        np_assert_allclose!(
            crate::softmax(&[1000.0, 0.0, 0.0, 0.0]),
            [1.0, 0.0, 0.0, 0.0],
            rtol = 1e-13
        );
        np_assert_allclose!(crate::softmax(&[1.0, 1.0]), [0.5, 0.5], rtol = 1e-13);
        np_assert_allclose!(
            crate::softmax(&[0.0, 1.0]),
            [1.0 / (1.0 + E), E / (1.0 + E)],
            rtol = 1e-13
        );

        // Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
        // converted to float.
        let x = [0.0, 1.0, 2.0, 3.0];
        let expected = [
            0.032_058_603_280_084_99,
            0.087_144_318_742_032_56,
            0.236_882_818_089_910_13,
            0.643_914_259_887_972_2,
        ];
        np_assert_allclose!(crate::softmax(&x), expected, rtol = 1e-13);

        // Translation property.  If all the values are changed by the same amount,
        // the softmax result does not change.
        np_assert_allclose!(
            crate::softmax(&x.map(|x| x + 100.0)),
            expected,
            rtol = 1e-13
        );
    }

    /// Corresponds to `TestLogSoftmax.*`
    #[test]
    fn test_log_softmax() {
        // test_log_softmax_basic
        np_assert_allclose!(
            crate::log_softmax(&[1000.0, 1.0]),
            [0., -999.],
            rtol = 1e-13
        );

        // test_log_softmax_scalar
        np_assert_allclose!(crate::log_softmax(&[1.0]), [0.], rtol = 1e-13);

        // test_log_softmax_translation
        let x = [0.0, 1.0, 2.0, 3.0];
        // Expected value computed using mpmath (with mpmath.mp.dps = 200)
        let expect = [
            -3.440_189_698_561_195_3,
            -2.440_189_698_561_195_3,
            -1.440_189_698_561_195_3,
            -0.440_189_698_561_195_33,
        ];
        np_assert_allclose!(crate::log_softmax(&x), expect, rtol = 1e-13);

        // Translation property.  If all the values are changed by the same amount,
        // the softmax result does not change.
        np_assert_allclose!(
            crate::log_softmax(&x.map(|x| x + 100.0)),
            expect,
            rtol = 1e-13
        );
    }
}
