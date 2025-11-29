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

#[cfg(test)]
mod tests {
    use core::f64::consts::TAU;

    const ATOL: f64 = 1.5e-15;

    /// Translated from `scipy.special.tests.test_basic.TestSystematic.test_diric`
    #[test]
    fn test_diric() {
        // Test behavior near multiples of tau. Regression test for scipy/scipy#4001.

        // n_odd = [1, 5, 25]
        let n_odd = [1, 5, 25];
        // x = np.array(2*np.pi + 1e-9).astype(np.float64)
        let x = TAU + 1e-9;
        // assert_allclose(special.diric(x, n_odd), 1.0, atol=ATOL, rtol=0)
        crate::np_assert_allclose!(n_odd.map(|n| crate::diric(n, x)), &[1.0; 3], atol = ATOL);

        // x = np.array(2*np.pi + 1e-15).astype(np.float64)
        let x = TAU + 1e-15;
        // assert_allclose(special.diric(x, n_odd), 1.0, atol=ATOL, rtol=0)
        crate::np_assert_allclose!(n_odd.map(|n| crate::diric(n, x)), &[1.0; 3], atol = ATOL);

        // n_even = [2, 4, 24]
        let n_even = [2, 4, 24];
        // x = np.array(2*np.pi + 1e-9).astype(np.float64)
        let x = TAU + 1e-9;
        // assert_allclose(special.diric(x, n_even), -1.0, atol=ATOL, rtol=0)
        crate::np_assert_allclose!(n_even.map(|n| crate::diric(n, x)), &[-1.0; 3], atol = ATOL);

        // Test at some values not near a multiple of pi
        // x = np.arange(0.2*np.pi, 1.0*np.pi, 0.2*np.pi)
        let x = [0.1 * TAU, 0.2 * TAU, 0.3 * TAU, 0.4 * TAU];
        // octave_result = [0.872677996249965, 0.539344662916632,
        //                  0.127322003750035, -0.206011329583298]
        let octave_result = [
            0.872_677_996_249_965,
            0.539_344_662_916_632,
            0.127_322_003_750_035,
            -0.206_011_329_583_298,
        ];
        // assert_allclose(special.diric(x, 3), octave_result, atol=ATOL, rtol=0)
        crate::np_assert_allclose!(x.map(|x| crate::diric(3, x)), &octave_result, atol = ATOL);
    }
}
