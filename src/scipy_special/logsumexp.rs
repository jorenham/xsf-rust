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
/// - [`logsumexp`]
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
/// `log_softmax` is more accurate than taking the log of [`softmax`] with inputs that make
/// [`softmax`] saturate.
///
/// # See also
/// - [`softmax`]
/// - [`logsumexp`]
#[must_use]
pub fn log_softmax(xs: &[f64]) -> Vec<f64> {
    let x_max = fmax(xs);
    let x_max = if x_max.is_finite() { x_max } else { 0.0 };

    let xs0 = xs.iter().map(|&x| x - x_max).collect::<Vec<f64>>();
    let s = xs0.iter().map(|&x| x.exp()).sum::<f64>().ln();
    xs0.iter().map(|&x| x - s).collect()
}
