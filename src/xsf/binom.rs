/// Binomial coefficient considered as a function of two real variables
///
/// For real arguments, the binomial coefficient is defined as
///
/// $$
/// \begin{align*}
/// \binom{n}{k}
/// &= {\Gamma(n+1) \over \Gamma(k+1) \\ \Gamma(n-k+1) } \\\\
/// &= {1 \over (n+1) \\ \Beta(k+1,\\, n-k+1)}
/// \end{align*}
/// $$
///
/// Where $\Gamma$ is the Gamma function ([`gamma`](crate::gamma)) and $\Beta$ the Beta function
/// ([`beta`](crate::beta)).
///
/// Corresponds to [`scipy.special.binom`][binom].
///
/// [binom]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.binom.html
///
/// # See also
/// - [`comb`](crate::comb) integer version "n choose k"
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
