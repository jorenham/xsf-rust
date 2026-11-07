//! Translated from <https://github.com/scipy/scipy/blob/38ecfef/scipy/special/_basic.py#L93-L192>

use crate::xsf::cephes::dd_real::two_prod;

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
    assert!(n != 0, "diric: n must be a positive integer, got n = 0");

    let x_half = x / 2.0;
    if x_half == 0.0 {
        return 1.0;
    }

    // Unlike in SciPy, sin(n x / 2) is evaluated at the exact product n (x / 2) = hi + lo, instead
    // of fixing the result at +1 or -1 where |sin(x / 2)| < 1e-7.
    //
    // - Near x = 2 pi k the rounding error of n (x / 2) is large relative to sin(n x / 2), e.g. 3e-10
    //   for diric(3, 2 pi + 1e-6), which the cutoff only hides for x even closer to 2 pi k.
    // - The cutoff is wrong for large n, e.g. it gives 1 for diric(2^30, 1e-7) ~ -0.005.
    let n_f64 = f64::from(n);
    let (hi, lo) = two_prod(n_f64, x_half);
    let ((sin_hi, cos_hi), (sin_lo, cos_lo)) = (hi.sin_cos(), lo.sin_cos());
    (sin_hi * cos_lo + cos_hi * sin_lo) / (n_f64 * x_half.sin())
}
