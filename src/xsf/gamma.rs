use num_complex::Complex;

use crate::utils::is_nonpos_int;
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::loggamma::cloggamma;

pub trait GammaArg: crate::sealed::Sealed {
    fn xsf_gamma(self) -> Self;
}

impl GammaArg for f64 {
    #[inline]
    fn xsf_gamma(self) -> Self {
        crate::xsf::cephes::gamma(self)
    }
}

impl GammaArg for Complex<f64> {
    #[inline]
    fn xsf_gamma(self) -> Self {
        self.via(cgamma)
    }
}

/// `xsf::gamma` for complex `z`, translated into pure Rust from `xsf/gamma.h` (scipy/xsf@56a9b96)
pub(crate) fn cgamma(z: C64) -> C64 {
    // Compute Gamma(z) using loggamma.
    if z.is_real() && is_nonpos_int(z.re) {
        // Poles
        return C64::NAN;
    }
    // Unlike in xsf, evaluate negative real non-poles on the real axis.
    if z.is_real() && z.re < 0.0 {
        let psi = crate::xsf::cephes::psi(z.re);
        // Near a zero of ψ its reflection terms (of order log(1 − x)) nearly cancel.
        // Keep the complex path when their rounding may obscure the derivative sign.
        let sign_margin = 32.0 * f64::EPSILON * (1.0 + (-z.re).ln_1p());
        if psi.abs() > sign_margin {
            let mut re = crate::xsf::cephes::gamma(z.re);
            // The real kernel's reflection denominator can overflow before Gamma underflows.
            if re == 0.0 {
                re = cloggamma(z).re.exp().copysign(re);
            }
            let derivative_negative = re.is_sign_negative() ^ psi.is_sign_negative();
            return c64(re, if derivative_negative { -z.im } else { z.im });
        }
    }
    cloggamma(z).exp()
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
/// On the negative real axis away from poles and stationary points, the imaginary zero follows
/// the sign of `y * Γ′(x)` for input `x + i y` with `y = ±0`.
/// A small imaginary residual can remain near a stationary point.
///
/// # Types
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
#[doc(alias = "gamma_lr", alias = "gamma_p")]
#[must_use]
#[inline]
pub fn gammainc(a: f64, x: f64) -> f64 {
    crate::xsf::cephes::igam(a, x)
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
#[doc(alias = "gamma_ur", alias = "gamma_q")]
#[must_use]
#[inline]
pub fn gammaincc(a: f64, x: f64) -> f64 {
    crate::xsf::cephes::igamc(a, x)
}

/// Inverse of [`gammainc`]
///
/// Given an input $p \in \[0, 1\]$, returns $x$ such that $p = P(a,x)$.
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
    crate::xsf::cephes::igami(a, p)
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
    crate::xsf::cephes::igamci(a, p)
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
    crate::xsf::cephes::lgam(x)
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
    crate::xsf::cephes::gammasgn(x)
}
