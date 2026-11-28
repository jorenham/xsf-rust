//! The complex `sici` and `shichi` are translated into pure Rust from `xsf/sici.h`
//! (scipy/xsf@7fb9276).
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

use crate::utils::EGAMMA;
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::consts::MAXLOG;
use crate::xsf::expint::{exp1_c64, exp1_half_c64, expi_c64};

fn sici_power_series(sgn: i32, z: C64) -> (C64, C64) {
    /* DLMF 6.6.5 and 6.6.6. If sgn = -1 computes si/ci, and if sgn = 1
     * computes shi/chi.
     */
    const TOL: f64 = f64::EPSILON;

    let mut fac = z;
    let mut s = fac;
    let mut c = C64::ZERO;
    for n in 1..100 {
        let n = f64::from(n);
        fac *= f64::from(sgn) * z / (2.0 * n);
        let term2 = fac / (2.0 * n);
        c += term2;
        fac *= z / (2.0 * n + 1.0);
        let term1 = fac / (2.0 * n + 1.0);
        s += term1;
        if term1.norm() < TOL * s.norm() && term2.norm() < TOL * c.norm() {
            break;
        }
    }
    (s, c)
}

fn sici_c64(z: C64) -> (C64, C64) {
    /* Compute sin/cos integrals at complex arguments. The algorithm
     * largely follows that of [1].
     */

    const EULER: f64 = EGAMMA;

    if z == C64::INFINITY {
        return (c64(FRAC_PI_2, 0.0), C64::ZERO);
    }
    if z == C64::NEG_INFINITY {
        return (c64(-FRAC_PI_2, 0.0), c64(0.0, PI.copysign(z.im)));
    }

    if z.norm() < 0.8 {
        // Use the series to avoid cancellation in si
        let (si, mut ci) = sici_power_series(-1, z);

        if z == C64::ZERO {
            ci = c64(f64::NEG_INFINITY, f64::NAN);
        } else {
            ci += EULER + z.ln();
        }
        return (si, ci);
    }

    let jz = C64::I * z;
    let scale_e1 = z.im.abs() > MAXLOG && z.is_finite();
    if z.re != 0.0 {
        // Unlike in xsf, DLMF 6.5.5/6.5.6 with `E1` itself, not through `Ei(w) = -E1(-w) +- pi i`,
        // whose `pi i` cancel only after absorbing the imaginary part of a small `E1`, e.g. for
        // huge `Re z`. Odd `Si` and DLMF 6.4.4 extend them to `Re z < 0`.
        let (rotation, mut ci) = if scale_e1 {
            let e1p = exp1_half_c64(jz);
            let e1m = exp1_half_c64(-jz);
            let difference = e1m - e1p;
            // A component swap avoids multiplying an overflowing part by zero.
            (c64(-difference.im, difference.re), -(e1p + e1m))
        } else {
            let e1p = exp1_c64(jz);
            let e1m = exp1_c64(-jz);
            (c64(0.0, 0.5) * (e1m - e1p), -0.5 * (e1p + e1m))
        };
        let si = rotation + FRAC_PI_2.copysign(z.re);
        if z.re < 0.0 {
            ci += c64(0.0, PI.copysign(z.im));
        }
        return (si, ci);
    }

    if scale_e1 {
        // On the imaginary axis the decaying term is negligible.
        let integral = -exp1_half_c64(c64(-z.im.abs(), 0.0)).re;
        let si = c64(0.0_f64.copysign(-z.im), integral.copysign(z.im));
        let ci = c64(integral, FRAC_PI_2.copysign(z.im));
        return (si, ci);
    }

    // on the imaginary axis, `E1(+-iz)` is on its branch cut, so through `Ei` (DLMF 6.5.5/6.5.6 plus
    // DLMF 6.4.4/6.4.6/6.4.7)
    let term1 = expi_c64(jz);
    let term2 = expi_c64(-jz);
    let si = c64(0.0, -0.5) * (term1 - term2);
    let mut ci = 0.5 * (term1 + term2);
    if z.im > 0.0 {
        ci += c64(0.0, FRAC_PI_2);
    } else if z.im < 0.0 {
        ci -= c64(0.0, FRAC_PI_2);
    }
    (si, ci)
}

fn shichi_c64(z: C64) -> (C64, C64) {
    /* Compute sinh/cosh integrals at complex arguments. The algorithm
     * largely follows that of [1].
     */
    const EULER: f64 = EGAMMA;

    if z == C64::INFINITY {
        return (C64::INFINITY, C64::INFINITY);
    }
    if z == C64::NEG_INFINITY {
        return (C64::NEG_INFINITY, C64::INFINITY);
    }
    if z.norm() < 0.8 {
        // Use the series to avoid cancellation in shi
        let (shi, mut chi) = sici_power_series(1, z);
        if z == C64::ZERO {
            chi = c64(f64::NEG_INFINITY, f64::NAN);
        } else {
            chi += EULER + z.ln();
        }
        return (shi, chi);
    }

    let term1 = expi_c64(z);
    let term2 = expi_c64(-z);
    let mut shi = 0.5 * (term1 - term2);
    let mut chi = 0.5 * (term1 + term2);
    if z.im > 0.0 {
        shi -= c64(0.0, FRAC_PI_2);
        chi += c64(0.0, FRAC_PI_2);
    } else if z.im < 0.0 {
        shi += c64(0.0, FRAC_PI_2);
        chi -= c64(0.0, FRAC_PI_2);
    } else if z.re < 0.0 {
        chi += c64(0.0, PI.copysign(z.im));
    }
    // Unlike in xsf, retry overflowing components with the dominant half-Ei.
    if z.re.abs() > MAXLOG && z.is_finite() && !(shi.is_finite() && chi.is_finite()) {
        return shichi_overflow(z, shi, chi);
    }
    (shi, chi)
}

#[cold]
fn shichi_overflow(z: C64, shi: C64, chi: C64) -> (C64, C64) {
    let sign = z.re.signum();
    let w = z * sign;
    let h = (0.5 * w.re).exp();
    let inv = 1.0 / w;
    let scale = inv.norm();
    if !h.is_finite() || !scale.is_normal() {
        return (shi, chi);
    }

    // DLMF 6.12.2 gives the Ei expansion. Rotating the E1 contour (DLMF 6.2.1)
    // by pi/4 away from the pole bounds the normalized remainder by 32 * 8! / |w|^8.
    let factorials = [5040.0, 720.0, 120.0, 24.0, 6.0, 2.0, 1.0, 1.0];
    let series = factorials.into_iter().fold(C64::ZERO, |s, c| s * inv + c);
    // The remainder of `q`: the normalized bound, with 8! = 8 * 7!, times `|1/w|`.
    let remainder = 32.0 * (8.0 * factorials[0]) * scale.powi(9);
    let (rotation, q) = (C64::cis(w.im), inv * series);
    let phase = rotation * q;
    // For `Re w > 0` the leading `1/w` term dominates `q`, so its evaluation has no severe
    // cancellation and both its parts are accurate. Given accurate `cis` and `exp`, each
    // projection is then accurate relative to the magnitudes of its two products.
    let re_size = (rotation.re * q.re).abs() + (rotation.im * q.im).abs();
    let im_size = (rotation.re * q.im).abs() + (rotation.im * q.re).abs();
    let half_ei = |t: f64| ((t * h) * 0.5) * h;
    let recover = |old: f64, projection: f64, size: f64| {
        if old.is_infinite() {
            // Include the phase, exponential and scaling roundoff, plus the bounded
            // complementary E1 and cut terms (DLMF 6.5.3-4).
            let error = half_ei(16.0 * f64::EPSILON * size + remainder) + 4.0;
            let value = half_ei(projection);
            // Near a phase zero even a finite projection can be inaccurate.
            if error < 256.0 * f64::EPSILON * value.abs() && (value.abs() + error).is_finite() {
                return value;
            }
        }
        old
    };
    (
        c64(recover(shi.re, sign * phase.re, re_size), recover(shi.im, sign * phase.im, im_size)),
        c64(recover(chi.re, phase.re, re_size), recover(chi.im, phase.im, im_size)),
    )
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
        let (si, ci) = sici_c64(self.into());
        (si.into(), ci.into())
    }

    #[inline]
    fn shichi(self) -> (Self, Self) {
        let (shi, chi) = shichi_c64(self.into());
        (shi.into(), chi.into())
    }
}

/// Sine and cosine integrals.
///
/// Corresponds to [`scipy.special.sici`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sici.html
///
/// # Definition
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
/// # Definition
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
