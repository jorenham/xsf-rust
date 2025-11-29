/// Binomial coefficient considered as a function of two real variables
///
/// For real arguments, the binomial coefficient is defined as
///
/// <sub>*n*</sub>C<sub>*k*</sub> = Γ(*n* + 1) / (Γ(*k* + 1) Γ(*n* - *k* + 1))
/// = 1 / ( (*n* + 1) Β(*n* - *k* + 1, *k* + 1) )
///
/// Where Γ is the Gamma function ([`gamma`](crate::gamma)) and Β is the Beta function ([`beta`](crate::beta)).
///
/// Corresponds to [`scipy.special.binom`][binom].
///
/// [binom]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.binom.html
///
/// # See also
/// - [`comb`](crate::comb) for the integer version
#[must_use]
#[inline]
pub fn binom(n: f64, k: f64) -> f64 {
    unsafe { crate::ffi::xsf::binom(n, k) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_binom_f64() {
        xsref::test("binom", "d_d-d", |x| crate::binom(x[0], x[1]));
    }
}
