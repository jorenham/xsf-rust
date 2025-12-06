mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait GammaArg: sealed::Sealed {
    fn xsf_gamma(self) -> Self;
}

impl GammaArg for f64 {
    #[inline]
    fn xsf_gamma(self) -> Self {
        unsafe { crate::ffi::xsf::gamma(self) }
    }
}

impl GammaArg for num_complex::Complex<f64> {
    #[inline]
    fn xsf_gamma(self) -> Self {
        unsafe { crate::ffi::xsf::gamma_1(self) }
    }
}

/// Gamma function $\Gamma(z)$ for real or complex argument.
///
/// It is defined as [^DLMF]:
///
/// $$
/// \Gamma(z) = \int_0^\infty t^{z-1} e^{-t} \dd t ,
/// $$
///
/// for $\Re(z) > 0$, and is extended to the complex plane by analytic continuation.
///
/// # Type Parameters
/// - `T`: either [`f64`] or [`Complex<f64>`](num_complex::Complex)
///
/// # See also
/// - [`gammainc`]: Regularized lower incomplete gamma function, $P(a,x)$
/// - [`gammaincc`]: Regularized upper incomplete gamma function, $Q(a,x)$
/// - [`gammaln`]: Natural logarithm of the absolute value of the gamma function
/// - [`rgamma`](crate::rgamma): Reciprocal of the gamma function, $1 \over \Gamma(z)$
/// - [`scipy.special.gamma`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gamma.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/5.2.E1>.
#[must_use]
#[inline]
pub fn gamma<T: GammaArg>(z: T) -> T {
    z.xsf_gamma()
}

/// Regularized lower incomplete gamma function $P(a,x)$
///
/// It is defined as [^DLMF]:
///
/// $$
/// P(a,x)
/// = {1 \over \Gamma(a)} \int_0^x t^{a-1} e^{-t} \dd t
/// = 1 - Q(a,x) ,
/// $$
///
/// for $a > 0$ and $x \ge 0$, with $Q(a,x)$ the regularized upper incomplete gamma function.
///
/// # Notes
/// Unlike [`scipy.special.gammainc`][scipy], this wraps the `gammainc` Cephes routine [^CEPHES]
/// rather than the Boost `gamma_p` routine [^BOOST].
///
/// # See also
/// - [`gammaincc`]: Regularized upper incomplete gamma function, $Q(a,x)$
/// - [`gammaincinv`]: Inverse of the regularized lower incomplete gamma function
/// - [`gammainccinv`]: Inverse of the regularized upper incomplete gamma function
/// - [`gamma`]: Gamma function $\Gamma(z)$
/// - [`scipy.special.gammainc`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.2.E4>.
/// [^CEPHES]: Cephes Math Library Release 2.4. Translated into C++ by SciPy developers in 2024.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(alias = "gamma_lr", alias = "gamma_p")]
#[must_use]
#[inline]
pub fn gammainc(a: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammainc(a, x) }
}

/// Regularized upper incomplete gamma function $Q(a,x)$
///
/// It is defined as [^DLMF]:
///
/// $$
/// Q(a,x)
/// = {1 \over \Gamma(a)} \int_x^\infty t^{a-1} e^{-t} \dd t
/// = 1 - P(a,x) ,
/// $$
///
/// for $a > 0$ and $x \ge 0$, with $P(a,x)$ the regularized lower incomplete gamma function.
///
/// # Notes
/// Unlike [`scipy.special.gammaincc`][scipy], this wraps the `gammaincc` Cephes routine [^CEPHES]
/// rather than the Boost `gamma_q` routine [^BOOST].
///
/// # See also
/// - [`gammainc`]: Regularized lower incomplete gamma function, $P(a,x)$
/// - [`gammaincinv`]: Inverse of the regularized lower incomplete gamma function
/// - [`gammainccinv`]: Inverse of the regularized upper incomplete gamma function
/// - [`gamma`]: Gamma function $\Gamma(z)$
/// - [`scipy.special.gammaincc`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaincc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.2.E4>.
/// [^CEPHES]: Cephes Math Library Release 2.4. Translated into C++ by SciPy developers in 2024.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(alias = "gamma_ur", alias = "gamma_q")]
#[must_use]
#[inline]
pub fn gammaincc(a: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaincc(a, x) }
}

/// Inverse of [`gammainc`]
///
/// Given an input $y \in \[0, 1\]$, returns $x$ such that $y = P(a,x)$.
/// Here, $P$ is the regularized lower incomplete gamma function; see [`gammainc`].
/// This is well-defined because the lower incomplete gamma function is monotonic [^DLMF].
///
/// # See also
/// - [`gammainc`]: Regularized lower incomplete gamma function, $P(a,x)$
/// - [`gammaincc`]: Regularized upper incomplete gamma function, $Q(a,x)$
/// - [`gammainccinv`]: Inverse of the regularized upper incomplete gamma function
/// - [`gamma`]: Gamma function $\Gamma(z)$
/// - [`scipy.special.gammaincinv`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaincinv.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.2.E4>.
#[must_use]
#[inline]
pub fn gammaincinv(a: f64, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaincinv(a, p) }
}

/// Inverse of [`gammaincc`]
///
/// Given an input $p \in \[0, 1\]$, returns $x$ such that $p = Q(a,x)$.
/// Here, $Q$ is the regularized upper incomplete gamma function; see [`gammaincc`].
/// This is well-defined because the upper incomplete gamma function is monotonic [^DLMF].
///
/// # See also
/// - [`gammainc`]: Regularized lower incomplete gamma function, $P(a,x)$
/// - [`gammaincc`]: Regularized upper incomplete gamma function, $Q(a,x)$
/// - [`gammaincinv`]: Inverse of the regularized lower incomplete gamma function
/// - [`gamma`]: Gamma function $\Gamma(z)$
/// - [`scipy.special.gammainccinv`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainccinv.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.2.E4>.
#[must_use]
#[inline]
pub fn gammainccinv(a: f64, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammainccinv(a, p) }
}

/// Logarithm of the absolute value of [`gamma`], $\ln{\abs{\Gamma(x)}}$
///
/// # See also
/// - [`gamma`]: Gamma function $\Gamma(z)$
/// - [`gammasgn`]: $\sgn \Gamma(x)$
/// - [`scipy.special.gammaln`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaln.html
#[doc(alias = "lgamma")]
#[must_use]
#[inline]
pub fn gammaln(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammaln(x) }
}

/// Sign of the Gamma function, $\sgn \Gamma(x)$
///
/// Since the Gamma function is never zero, this is always either `1.0` or `-1.0` (or [`f64::NAN`]).
///
/// # See also
/// - [`gamma`]: Gamma function $\Gamma(z)$
/// - [`gammaln`]: Natural logarithm of the absolute value of the gamma function
/// - [`scipy.special.gammasgn`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammasgn.html
#[must_use]
#[inline]
pub fn gammasgn(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gammasgn(x) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_gamma_f64() {
        xsref::test("gamma", "d-d", |x| crate::gamma(x[0]));
    }

    #[test]
    fn test_gamma_c64() {
        xsref::test("gamma", "cd-cd", |x| crate::gamma(c64(x[0], x[1])));
    }

    #[test]
    fn test_gammainc() {
        xsref::test("gammainc", "d_d-d", |x| crate::gammainc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincc() {
        xsref::test("gammaincc", "d_d-d", |x| crate::gammaincc(x[0], x[1]));
    }

    #[test]
    fn test_gammaincinv() {
        xsref::test("gammaincinv", "d_d-d", |x| crate::gammaincinv(x[0], x[1]));
    }

    #[test]
    fn test_gammainccinv() {
        xsref::test("gammainccinv", "d_d-d", |x| crate::gammainccinv(x[0], x[1]));
    }

    #[test]
    fn test_gammaln() {
        xsref::test("gammaln", "d-d", |x| crate::gammaln(x[0]));
    }

    #[test]
    fn test_gammasgn() {
        xsref::test("gammasgn", "d-d", |x| crate::gammasgn(x[0]));
    }
}
