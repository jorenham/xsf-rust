//! Translated from <https://github.com/scipy/scipy/blob/38ecfef/scipy/special/_basic.py#L93-L192>

use core::f64::consts::TAU;

/// Dirichlet kernel, also known as the periodic sinc function.
///
/// It is defined as
///
/// *D<sub>n</sub>(x) = sin(n x/2) / (n sin(x/2))*
///
/// for *n > 0*.
///
/// Corresponds to [`scipy.special.diric`][diric] in SciPy.
///
/// [diric]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.diric.html
///
/// **Note:** The argument order in this crate is `(n, x)`, whereas SciPy uses `(x, n)`.
///
/// # Panics
/// - Panics if `n` is zero.
#[must_use]
#[inline]
pub fn diric(n: u32, x: f64) -> f64 {
    // Empirical minval for 64-bit float computations
    // where sin(x/2) < minval, result is fixed at +1 or -1
    // For f64: eps < 1e-15, so minval = 1e-7
    const MINVAL: f64 = 1e-7;

    // Handle invalid n (n <= 0 is prevented by u32, but check n == 0)
    assert!(n != 0, "diric: n must be a positive integer, got n = 0");

    let x_half = x / 2.0;
    let denom = x_half.sin();

    // Check if |sin(x/2)| < minval
    if denom.abs() < MINVAL {
        return (-1.0_f64).powf((x / TAU).round() * (f64::from(n) - 1.0));
    }

    // Normal case: sin(n*x/2) / (n*sin(x/2))
    let n_f64 = f64::from(n);
    (n_f64 * x_half).sin() / (n_f64 * denom)
}
