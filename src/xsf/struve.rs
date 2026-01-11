/// Integral of the Struve $H_0$ function
///
/// $$ \int_0^x H_0(t) \dd t $$
///
/// Corresponds to [`scipy.special.itstruve0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itstruve0.html
///
/// # See also
/// - [`struve_h`]: Struve $H_v$ function
/// - [`it2struve0`]
#[doc(alias = "itstruve_h0", alias = "it1struve0")]
#[must_use]
#[inline]
pub fn itstruve0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::itstruve0(x) }
}

/// Integral related to the Struve $H_0$ function
///
/// $$ \int_0^x {H_0(t) \over t} \dd t $$
///
/// Corresponds to [`scipy.special.it2struve0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2struve0.html
///
/// # See also
/// - [`struve_h`]: Struve $H_v$ function
/// - [`itstruve0`]
#[doc(alias = "it2struve_l0")]
#[must_use]
#[inline]
pub fn it2struve0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::it2struve0(x) }
}

/// Integral of the modified Struve $L_0$ function
///
/// $$ \int_0^x L_0(t) \dd t $$
///
/// Corresponds to [`scipy.special.itmodstruve0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itmodstruve0.html
///
/// # See also
/// - [`struve_l`]: Struve $L_v$ function
#[doc(alias = "itstruve_l0")]
#[must_use]
#[inline]
pub fn itmodstruve0(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::itmodstruve0(x) }
}

/// Struve $H_v$ function
///
/// Corresponds to [`scipy.special.struve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.struve.html
///
/// # Definition
///
/// $$
/// H_v(x) =
/// (x/2)^{v+1}
/// \sum_{n=0}^\infty {(-1)^n (x/2)^{2n} \over \Gamma(n+3/2) \\ \Gamma(n+v+3/2)} \\ ,
/// $$
///
/// where $\Gamma$ is the Gamma function.
///
/// # See also
/// - [`struve_l`]: Struve $L$ function
/// - [`itstruve0`]: Integral of $H_0$
/// - [`it2struve0`]: Integral related to $H_0$
#[doc(alias = "struve")]
#[must_use]
#[inline]
pub fn struve_h(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::struve_h(v, x) }
}

/// Modified Struve $L_v$ function
///
/// Corresponds to [`scipy.special.modstruve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modstruve.html
///
/// # Definition
///
/// $$
/// \begin{align*}
/// L_v(x)
/// &= (z/2)^{v+1} \sum_{n=0}^\infty {(-z/2)^{2n} \over \Gamma(n+3/2) \\ \Gamma(n+v+3/2)} \\\\
/// &= -i e^{-\pi iv/2} H_v(ix) \\ ,
/// \end{align*}
/// $$
///
/// where $\Gamma$ is the Gamma function, and $H_v$ the Struve function of order $v$ ([`struve_h`]).
///
/// # See also
/// - [`struve_h`]: Struve $H$ function
/// - [`itmodstruve0`]: Integral of $L_0$
#[doc(alias = "modstruve")]
#[must_use]
#[inline]
pub fn struve_l(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::struve_l(v, x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_itstruve0_f64() {
        xsref::test("itstruve0", "d-d", |x| crate::itstruve0(x[0]));
    }

    #[test]
    fn test_it2struve0_f64() {
        xsref::test("it2struve0", "d-d", |x| crate::it2struve0(x[0]));
    }

    #[test]
    fn test_itmodstruve0_f64() {
        xsref::test("itmodstruve0", "d-d", |x| crate::itmodstruve0(x[0]));
    }

    #[test]
    fn test_struve_h_f64() {
        xsref::test("struve_h", "d_d-d", |x| crate::struve_h(x[0], x[1]));
    }

    #[test]
    fn test_struve_l_f64() {
        xsref::test("struve_l", "d_d-d", |x| crate::struve_l(x[0], x[1]));
    }
}
