//! `cbrt` uses the implementation from the Rust standard library, which (at least on Linux) is more
//! accurate than the Cephes one in xsf.

/// Cube root of $x$, $\sqrt\[3\]{x}$
///
/// This corresponds to [`scipy.special.cbrt`][cbrt] in SciPy.
///
/// [cbrt]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cbrt.html
#[must_use]
#[inline]
pub fn cbrt(x: f64) -> f64 {
    x.cbrt()
}
