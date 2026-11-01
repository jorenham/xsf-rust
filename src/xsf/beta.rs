/// Beta function $\B(a,b)$
///
/// This function is defined as [^DLMF]:
///
/// $$
/// \B(a,b)
/// = \int_0^1 t^{a-1} (1-t)^{b-1} \dd t
/// = {\Gamma(a)\Gamma(b) \over \Gamma(a+b)} ,
/// $$
///
/// where $\Gamma$ is the Gamma function.
///
/// # See also
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function $\I_x(a,b)$
/// - [`betaln`]: Natural logarithm of the absolute value of the Beta function
/// - [`scipy.special.beta`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.beta.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/5.12.E1>.
#[must_use]
#[inline]
pub fn beta(a: f64, b: f64) -> f64 {
    crate::xsf::cephes::beta(a, b)
}

/// Natural logarithm of the absolute value of [`beta`], $\ln{\abs{\B(a,b)}}$
///
/// # See also
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function $\I_x(a,b)$
/// - [`beta`]: Beta function $\B(a,b)$
/// - [`scipy.special.betaln`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaln.html
#[doc(alias = "lbeta", alias = "logbeta", alias = "beta_ln")]
#[must_use]
#[inline]
pub fn betaln(a: f64, b: f64) -> f64 {
    crate::xsf::cephes::lbeta(a, b)
}
