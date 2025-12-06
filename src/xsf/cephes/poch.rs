/// Rising factorial $\rpow x m$
///
/// $$\rpow x m = {\Gamma(x+m) \over \Gamma(x)}$$
///
/// Corresponds to [`scipy.special.poch`][poch] in SciPy.
///
/// [poch]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.poch.html
///
/// # See also
/// - [`pow_falling`]: falling factorial $\fpow x m$
/// - [`gamma`](crate::gamma): gamma function $\Gamma(x)$
#[doc(alias = "poch")]
#[must_use]
#[inline]
pub fn pow_rising(x: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::poch(x, m) }
}

/// Falling factorial $\fpow x m$
///
/// $$\fpow x m = {\Gamma(x+1) \over \Gamma(x-m+1)}$$
///
/// Note that there is no `scipy.special` analogue for this function, but it can be expressed in
/// terms of the rising factorial as `pow_rising(x - m + 1, m)`.
///
/// # See also
/// - [`pow_rising`]: rising factorial $\rpow x m$
/// - [`gamma`](crate::gamma): gamma function $\Gamma(x)$
#[must_use]
#[inline]
pub fn pow_falling(x: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::poch(x - m + 1.0, m) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_pow_rising() {
        xsref::test("poch", "d_d-d", |x| crate::pow_rising(x[0], x[1]));
    }
}
