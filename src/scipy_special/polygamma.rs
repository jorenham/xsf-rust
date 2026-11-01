//! Translated from <https://github.com/scipy/scipy/blob/38ecfef/scipy/special/_basic.py#L1539-L1579>

/// Polygamma function ψ<sup>(n)</sup>(x)
///
/// The n-th derivative of the [`digamma`](crate::digamma) function.
///
/// Corresponds to [`scipy.special.polygamma`][polygamma] in scipy
///
/// [polygamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.polygamma.html
///
/// # See also
/// - [`digamma`](crate::digamma): Digamma function (0th derivative)
#[must_use]
#[inline]
pub fn polygamma(n: u32, x: f64) -> f64 {
    if n == 0 {
        crate::digamma(x)
    } else {
        let n1p = f64::from(n) + 1.0;
        let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
        sign * crate::gamma(n1p) * crate::zeta(n1p, x)
    }
}
