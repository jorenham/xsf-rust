//! Translated from `scipy/special/_logsumexp.py`

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
/// - [`logaddexp2`](crate::logaddexp2)
#[must_use]
#[allow(clippy::float_cmp, clippy::missing_panics_doc)]
pub fn logsumexp(xs: &[f64]) -> f64 {
    // Find element with maximum real part, since this is what affects the magnitude
    // of the exponential. Possible enhancement: include log of `b` magnitude in `a`.
    let x_max = xs
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| if a >= b { a } else { b }); // >= will propagate NaNs

    // fast-path for -/+inf and NaN
    if !x_max.is_finite() {
        return x_max;
    }

    // for precision, these terms are separated out of the main sum.
    let m: f64 = xs
        .iter()
        .fold(0.0, |acc, &x| if x == x_max { acc + 1.0 } else { acc });
    let xs = xs
        .iter()
        .map(|&x| if x == x_max { f64::NEG_INFINITY } else { x });

    // Shift, exponentiate, and sum
    let s = xs.map(|x| (x - x_max).exp()).sum::<f64>();

    // The log functions need positive arguments
    let s = if s < -1.0 { -s - 2.0 } else { s };

    // Take log and undo shift
    s.ln_1p() + m.ln() + x_max
}

#[cfg(test)]
mod tests {
    //! Translated from `scipy.special.tests.test_logsumexp.TestLogSumExp` at
    //! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_logsumexp.py>

    use crate::np_assert_allclose;
    use core::f64::consts::LN_2;

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
}
