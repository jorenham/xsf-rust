//! The complex `sici` and `shichi` are translated into pure Rust from `xsf/sici.h` (xsf v0.2.2),
//! which was translated from Cython into C++ by SciPy developers in 2024.
//!
//! Implementation of sin/cos/sinh/cosh integrals for complex arguments
//!
//! Sources
//! \[1\] Fredrik Johansson and others. mpmath: a Python library for
//!     arbitrary-precision floating-point arithmetic (version 0.19),
//!     December 2013. <http://mpmath.org/>.
//! \[2\] NIST, "Digital Library of Mathematical Functions",
//!     <https://dlmf.nist.gov/>
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use core::f64::consts::{FRAC_PI_2, PI};

use num_complex::{Complex, c64};

use crate::utils::SCIPY_EULER;
use crate::xsf::complex::cmul;
use crate::xsf::log::clog;

fn sici_power_series(sgn: i32, z: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    /* DLMF 6.6.5 and 6.6.6. If sgn = -1 computes si/ci, and if sgn = 1
     * computes shi/chi.
     */
    const TOL: f64 = f64::EPSILON;

    let mut fac = z;
    let mut s = fac;
    let mut c = c64(0.0, 0.0);
    for n in 1..100 {
        let n = f64::from(n);
        fac = cmul(fac, f64::from(sgn) * z / (2.0 * n));
        let term2 = fac / (2.0 * n);
        c += term2;
        fac = cmul(fac, z / (2.0 * n + 1.0));
        let term1 = fac / (2.0 * n + 1.0);
        s += term1;
        if term1.norm() < TOL * s.norm() && term2.norm() < TOL * c.norm() {
            break;
        }
    }
    (s, c)
}

#[allow(clippy::float_cmp)]
fn sici_c64(z: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    /* Compute sin/cos integrals at complex arguments. The algorithm
     * largely follows that of [1].
     */

    const EULER: f64 = SCIPY_EULER;

    if z == c64(f64::INFINITY, 0.0) {
        return (c64(FRAC_PI_2, 0.0), c64(0.0, 0.0));
    }
    if z == c64(f64::NEG_INFINITY, 0.0) {
        return (c64(-FRAC_PI_2, 0.0), c64(0.0, PI));
    }

    if z.norm() < 0.8 {
        // Use the series to avoid cancellation in si
        let (si, mut ci) = sici_power_series(-1, z);

        if z == c64(0.0, 0.0) {
            ci = c64(f64::NEG_INFINITY, f64::NAN);
        } else {
            ci += EULER + clog(z);
        }
        return (si, ci);
    }

    // DLMF 6.5.5/6.5.6 plus DLMF 6.4.4/6.4.6/6.4.7
    let jz = cmul(c64(0.0, 1.0), z);
    let term1 = crate::expi(jz);
    let term2 = crate::expi(-jz);
    let mut si = cmul(c64(0.0, -0.5), term1 - term2);
    let mut ci = 0.5 * (term1 + term2);
    if z.re == 0.0 {
        if z.im > 0.0 {
            ci += c64(0.0, FRAC_PI_2);
        } else if z.im < 0.0 {
            ci -= c64(0.0, FRAC_PI_2);
        }
    } else if z.re > 0.0 {
        si -= FRAC_PI_2;
    } else {
        si += FRAC_PI_2;
        if z.im >= 0.0 {
            ci += c64(0.0, PI);
        } else {
            ci -= c64(0.0, PI);
        }
    }
    (si, ci)
}

#[allow(clippy::float_cmp)]
fn shichi_c64(z: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    /* Compute sinh/cosh integrals at complex arguments. The algorithm
     * largely follows that of [1].
     */
    const EULER: f64 = SCIPY_EULER;

    if z == c64(f64::INFINITY, 0.0) {
        return (c64(f64::INFINITY, 0.0), c64(f64::INFINITY, 0.0));
    }
    if z == c64(f64::NEG_INFINITY, 0.0) {
        return (c64(f64::NEG_INFINITY, 0.0), c64(f64::INFINITY, 0.0));
    }
    if z.norm() < 0.8 {
        // Use the series to avoid cancellation in shi
        let (shi, mut chi) = sici_power_series(1, z);
        if z == c64(0.0, 0.0) {
            chi = c64(f64::NEG_INFINITY, f64::NAN);
        } else {
            chi += EULER + clog(z);
        }
        return (shi, chi);
    }

    let term1 = crate::expi(z);
    let term2 = crate::expi(-z);
    let mut shi = 0.5 * (term1 - term2);
    let mut chi = 0.5 * (term1 + term2);
    if z.im > 0.0 {
        shi -= c64(0.0, FRAC_PI_2);
        chi += c64(0.0, FRAC_PI_2);
    } else if z.im < 0.0 {
        shi += c64(0.0, FRAC_PI_2);
        chi -= c64(0.0, FRAC_PI_2);
    } else if z.re < 0.0 {
        chi += c64(0.0, PI);
    }
    (shi, chi)
}

pub trait SiciArg: crate::sealed::Sealed + Sized {
    fn sici(self) -> (Self, Self);
    fn shichi(self) -> (Self, Self);
}

impl SiciArg for f64 {
    #[inline]
    fn sici(self) -> (Self, Self) {
        crate::xsf::cephes::sici(self)
    }

    #[inline]
    fn shichi(self) -> (Self, Self) {
        crate::xsf::cephes::shichi(self)
    }
}

impl SiciArg for num_complex::Complex64 {
    #[inline]
    fn sici(self) -> (Self, Self) {
        sici_c64(self)
    }

    #[inline]
    fn shichi(self) -> (Self, Self) {
        shichi_c64(self)
    }
}

/// Sine and cosine integrals.
///
/// Corresponds to [`scipy.special.sici`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sici.html
///
/// # Definitions
///
/// $$
/// \begin{align*}
/// \Si(z) &= \int_0^z {\sin t \over t} \dd t \\\\
/// \Ci(z) &= \gamma + \ln(z) + \int_0^z {\cos t - 1 \over t} \dd t
/// \end{align*}
/// $$
///
/// where $\gamma$ is Euler's constant and $\ln$ is the principal branch of the logarithm.
///
/// # Arguments
/// - $z$ - real- or complex-valued argument
///
/// # Returns
/// - $\Si(z)$ - Sine integral
/// - $\Ci(z)$ - Cosine integral
///
/// # See also
/// - [`shichi`]: $\Shi(z)$ and $\Chi(z)$
/// - [`exp1`](crate::exp1): Exponential integral $E_1$
/// - [`expi`](crate::expi): Exponential integral $Ei$
#[inline]
pub fn sici<T: SiciArg>(z: T) -> (T, T) {
    z.sici()
}

/// Hyperbolic sine and cosine integrals.
///
/// Corresponds to [`scipy.special.shichi`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.shichi.html
///
/// # Definitions
///
/// $$
/// \begin{align*}
/// \Shi(z) &= \int_0^z {\sinh t \over t} \dd t \\\\
/// \Chi(z) &= \gamma + \ln(z) + \int_0^z {\cosh t - 1 \over t} \dd t
/// \end{align*}
/// $$
///
/// where $\gamma$ is Euler's constant and $\ln$ is the principal branch of the logarithm.
///
/// # Arguments
/// - $z$ - real- or complex-valued argument
///
/// # Returns
/// - $\Shi(z)$ - Hyperbolic sine integral
/// - $\Chi(z)$ - Hyperbolic cosine integral
///
/// # See also
/// - [`sici`]: $\Si(z)$ and $\Ci(z)$
/// - [`exp1`](crate::exp1): Exponential integral $E_1$
/// - [`expi`](crate::expi): Exponential integral $Ei$
#[inline]
pub fn shichi<T: SiciArg>(z: T) -> (T, T) {
    z.shichi()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_sici_f64() {
        xsref::test("sici", "d-d_d", |x| crate::sici(x[0]));
    }

    #[test]
    fn test_sici_c64() {
        xsref::test("sici", "cd-cd_cd", |x| crate::sici(c64(x[0], x[1])));
    }

    #[test]
    fn test_shichi_f64() {
        xsref::test("shichi", "d-d_d", |x| crate::shichi(x[0]));
    }

    #[test]
    fn test_shichi_c64() {
        xsref::test("shichi", "cd-cd_cd", |x| crate::shichi(c64(x[0], x[1])));
    }
}
