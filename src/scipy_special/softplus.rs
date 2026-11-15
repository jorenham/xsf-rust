use crate::LogAddExpArg;

/// $\ln(1 + e^x)$
///
/// Corresponds to [`scipy.special.softplus`][softplus] in SciPy.
///
/// [softplus]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.softplus.html
#[inline]
pub fn softplus<T: LogAddExpArg>(x: T) -> T {
    x.npy_logaddexp(T::ZERO)
}
