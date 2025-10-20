use crate::LogAddExpArg;

/// ln(1 + e<sup>x</sup>)
///
/// Corresponds to [`scipy.special.softplus`][softplus] in SciPy.
///
/// [softplus]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.softplus.html
#[inline]
pub fn softplus<T: LogAddExpArg>(x: T) -> T {
    x.npy_logaddexp(T::ZERO)
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    #[test]
    fn test_softplus_f64() {
        let x = [-100.0, -1.0, 0.0, 1.0, 100.0];
        let y = x.map(|xi: f64| xi.exp().ln_1p());
        np_assert_allclose!(x.map(crate::softplus), y, atol = 1e-15);
    }
}
