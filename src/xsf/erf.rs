//! The `voigt_profile` is translated into pure Rust from `xsf/erf.h` (scipy/xsf@09c93c7), and the
//! other functions use `xsf/faddeeva.rs`.

use core::f64::consts::{FRAC_1_SQRT_2, FRAC_2_PI, PI};

use num_complex::Complex64;

use crate::utils::{SQRT_FRAC_2_PI, SQRT_TAU_DOWN, mul_exp, pow2};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, FRAC_PI_2_DD, two_prod};
use crate::xsf::phase::{self, Quadrant};

pub trait ErfArg: crate::sealed::Sealed {
    fn erf(self) -> Self;
    fn erfc(self) -> Self;
    fn erfcx(self) -> Self;
    fn erfi(self) -> Self;
    fn dawsn(self) -> Self;
}

impl ErfArg for f64 {
    #[inline]
    fn erf(self) -> Self {
        crate::xsf::cephes::erf(self)
    }
    #[inline]
    fn erfc(self) -> Self {
        crate::xsf::cephes::erfc(self)
    }
    #[inline]
    fn erfcx(self) -> Self {
        crate::xsf::faddeeva::erfcx(self)
    }
    #[inline]
    fn erfi(self) -> Self {
        crate::xsf::faddeeva::erfi(self)
    }
    #[inline]
    fn dawsn(self) -> Self {
        crate::xsf::faddeeva::dawson(self)
    }
}

impl ErfArg for Complex64 {
    #[inline]
    fn erf(self) -> Self {
        self.via(crate::xsf::faddeeva::erf_c64)
    }
    #[inline]
    fn erfc(self) -> Self {
        self.via(erfc_c64)
    }
    #[inline]
    fn erfcx(self) -> Self {
        self.via(crate::xsf::faddeeva::erfcx_c64)
    }
    #[inline]
    fn erfi(self) -> Self {
        self.via(crate::xsf::faddeeva::erfi_c64)
    }
    #[inline]
    fn dawsn(self) -> Self {
        self.via(crate::xsf::faddeeva::dawson_c64)
    }
}

fn erfc_c64(z: C64) -> C64 {
    gaussian_tail(z, 2.0)
}

// k/2 * erfc(z * sqrt(k/2)) for k = 1 or 2, with the Gaussian factor exp(-k*z^2/2) of the
// unscaled z, so that ndtr shares the overflow and zero-sign recovery of erfc.
#[inline]
pub(super) fn gaussian_tail(z: C64, k: f64) -> C64 {
    let raw_erfc = |z: C64| crate::xsf::faddeeva::erfc_c64(z * (0.5 * k).sqrt()) * (0.5 * k);

    let (x_abs, y_abs) = (z.re.abs(), z.im.abs());
    // Unlike in xsf, try certifying overflow before evaluation in this wedge,
    // where the exact Gaussian factor exceeds MAX.
    if y_abs >= 64.0 && x_abs <= y_abs * 0.5 {
        if let Some(result) = erfc_overflow(z, k) {
            return result;
        }
    }
    let result = raw_erfc(z);
    if result.is_nan() && z.is_finite() && x_abs == y_abs && (k * z.re * z.re).is_infinite() {
        return erfc_overrange_diagonal(z, k, result);
    }
    if result.is_nan() {
        return erfc_overflow(z, k)
            .or_else(|| {
                // Unlike in xsf, recover a finite component independently of its infinite companion.
                if result.re.is_infinite() || result.im.is_infinite() {
                    erfc_rescale_overflow(z, k, result)
                } else {
                    None
                }
            })
            .unwrap_or(result);
    }
    if result.re.is_infinite() && result.im.is_infinite() {
        // Unlike in xsf, include the Faddeeva factor when deciding overflow.
        // The rounded phase has at least one radian between representable values here.
        return if (k * x_abs) * y_abs > pow2(53) {
            erfc_overflow(z, k)
        } else {
            erfc_rescale_overflow(z, k, result)
        }
        .unwrap_or(result);
    }
    // Unlike in xsf, keep the sign of a component that rounds to zero in the Gaussian tail.
    if result.im == 0.0 && z.im != 0.0 {
        return gaussian_tail_zeros(z, k, result);
    }
    result
}

#[cold]
fn erfc_rescale_overflow(z: C64, k: f64, result: C64) -> Option<C64> {
    // The continued fraction takes the scaled argument, the phase and exponent the unscaled one.
    let t = z * (0.5 * k).sqrt();
    let (x, y) = (t.re.abs(), t.im);
    let partial = result.is_nan();
    let min_real = if partial { 8.0 } else { 32.0 };
    // In this region the general Gaussian continued fraction has degree four through nine.
    if !(x >= min_real && x < y.abs() && x + y.abs() <= 512.0) {
        return None;
    }
    let (w, degree) = crate::xsf::faddeeva::w_continued_fraction(-y, x);
    let (x0, y0) = (z.re.abs(), z.im);
    let theta = (Dd::new(-x0) * y0).mul_pwr2(k);
    let (sin, cos) = theta.sin_cos();
    let projected = c64(cos, sin) * w;
    let exp = ((Dd::new(y0.abs()) - x0) * (Dd::new(y0.abs()) + x0)).mul_pwr2(0.5 * k);
    // Keep the upper and lower overflow tests at the same scaled magnitude.
    let scale = mul_exp(f64::MIN_POSITIVE, exp.hi) * (1.0 + exp.lo) * (0.5 * k);
    let limit = f64::MAX * f64::MIN_POSITIVE;

    // Bound the remainder of the Gaussian Cauchy transform (DLMF 7.7.2)
    // by 1/(x*(|z|-4)^8).
    // Pad for the continued fraction and its rotation and scale.
    let distance_squared = (t.norm() - 4.0).powi(2);
    let remainder = 1.0 / (x * distance_squared.powi(4));
    let roundoff = w.l1_norm() * pow2(-44);
    let error = scale * (remainder + roundoff);
    let scaled = projected * scale;
    let finite = |val: f64, error: f64| {
        let a = val.abs();
        (a > error && a + error < limit).then(|| val / f64::MIN_POSITIVE)
    };
    let reflect_real = |re| if z.re < 0.0 { k - re } else { re };
    let certify = |val: f64, error: f64| {
        if val.abs() - error > limit {
            Some(f64::INFINITY.copysign(val))
        } else {
            finite(val, error)
        }
    };
    if partial {
        let companion = |val: f64| {
            certify(val, error).or_else(|| {
                // Unlike in xsf, tighten the overflow test for an infinite companion
                // whose sign is already known. The degree-n remainder uses the squared
                // monic Hermite norm n!/2^n for the normalized weight (DLMF Table 18.3.1)
                // and the squared distance to its zeros per level.
                if !(val.abs() > error && degree > 4.0) {
                    return None;
                }
                let mut remainder = remainder;
                let mut n = 5.0;
                while n <= degree {
                    remainder *= n / (2.0 * distance_squared);
                    n += 1.0;
                }
                certify(val, scale * (remainder + roundoff))
            })
        };
        // Unlike in xsf, recover a finite part even when the NaN part overflows.
        let repair = |missing, overflowed| {
            let value = certify(missing, error)?;
            let other = companion(overflowed);
            (value.is_finite() || other.is_some()).then_some((value, other))
        };
        return if result.re.is_nan() {
            repair(scaled.re, scaled.im)
                .map(|(re, im)| c64(reflect_real(re), im.unwrap_or(result.im)))
        } else {
            repair(scaled.im, scaled.re).map(|(im, re)| c64(re.map_or(result.re, reflect_real), im))
        };
    }
    Some(c64(reflect_real(certify(scaled.re, error)?), certify(scaled.im, error)?))
}

// Recover the zero signs of k/2 * erfc(z * sqrt(k/2)), with k = 1 or 2.
#[cold]
fn gaussian_tail_zeros(z: C64, k: f64, result: C64) -> C64 {
    let reflected = z.re < 0.0;
    let limit = if reflected { k } else { 0.0 };
    if result.re != limit {
        return result;
    }

    // Reflection and conjugation leave the imaginary part unchanged.
    let (x, y) = (z.re.abs(), z.im);
    let ay = y.abs();
    // The Laplace integral bounds the reflected tail's magnitude by
    // exp(-k(x²-y²)/2)/(sqrt(pi)*x*sqrt(k/2)). The monotone lower bound
    // x(x-|y|) on x²-y² leaves room for rounding below half the smallest subnormal.
    if !z.is_finite() || k * ((x - ay) * x) < 1500.0 {
        return result;
    }

    let hi = -x * y;
    let phase = k * hi;
    let signs = if phase.abs() <= 1.0 {
        // The integral along the vertical segment has positive cosine throughout.
        // The leading real part exceeds cos(1) - |y|/x, well above the remainder.
        c64(1.0, -y)
    } else {
        let (rotation, error) = if phase.abs() > pow2(53) {
            // Unlike in xsf, reduce the original product even when it overflows.
            // The reducer doubles its product; x/2 is exact in this wide-phase region.
            let (quadrant, r) = phase::reduce_product_phase((-0.5 * k) * x, y);
            let (sin, cos) = r.sin_cos();
            let (re, im) = quadrant.rotate(sin, cos);
            (c64(re, im), 0.0)
        } else {
            let lo = (-x).mul_add(y, -hi);
            let (quadrant, r) = recentered_phase(phase, k * lo);
            let Some(rotation) = phase_rotation(quadrant, r) else {
                return result;
            };
            rotation
        };
        let leading = c64(1.0, -y / x) * rotation;
        // Integration by parts bounds the normalized remainder by sqrt(2)/(k*x²).
        let error = error + (2.0 / k / x) / x + pow2(-44);
        if leading.re.abs().min(leading.im.abs()) <= error {
            return result;
        }
        leading
    };
    let im = 0.0_f64.copysign(signs.im);
    if reflected {
        c64(k, im)
    } else {
        c64(0.0_f64.copysign(signs.re), im)
    }
}

#[cold]
fn erfc_overrange_diagonal(z: C64, k: f64, result: C64) -> C64 {
    let x = z.re.abs();
    let (quadrant, r) = phase::reduce_diagonal_phase(x, k);
    let (sin, cos) = r.sin_cos();
    let (re, im) = quadrant.rotate(sin.abs().copysign(r), cos.abs());
    // On |r| <= pi/4, Taylor bounds give 7|r|/8 <= |sin r| <= |r|
    // and 11/16 <= cos r <= 1. Odd quadrants exchange the magnitudes.
    let sin_bounds = (0.875 * r.abs(), r.abs());
    let cos_bounds = (0.6875, 1.0);
    let (re_bounds, im_bounds) = if matches!(quadrant, Quadrant::Zero | Quadrant::HalfTurn) {
        (cos_bounds, sin_bounds)
    } else {
        (sin_bounds, cos_bounds)
    };

    // Unlike in xsf, rotate the diagonal tail before dividing to retain subnormals.
    let component = |factor: f64, (lo, hi): (f64, f64), old: f64| {
        const PAD: f64 = 8.0 * f64::EPSILON;
        if !old.is_nan() {
            return old;
        }
        let scale = 0.5 * SQRT_FRAC_2_PI * (0.5 * k).sqrt();
        let value = (scale * factor) / x;
        let lower = ((scale * lo) * (1.0 - PAD)) / x;
        let upper = ((scale * hi) * (1.0 + PAD)) / x;
        if lower.classify() == upper.classify() && lower.classify() == value.classify() {
            value
        } else {
            old
        }
    };
    let re = if z.re < 0.0 && result.re.is_nan() {
        k
    } else {
        component(re, re_bounds, result.re)
    };
    c64(re, component(-z.im.signum() * im, im_bounds, result.im))
}

#[cold]
fn erfc_overflow(z: C64, k: f64) -> Option<C64> {
    // Leave roundoff margin above ln(6*MAX).
    const LOG_OVERFLOW_BOUND: f64 = 713.0;

    let (x, y) = (z.re, z.im);
    let (x_abs, y_abs) = (x.abs(), y.abs());
    if !z.is_finite() || x == 0.0 || x_abs >= y_abs {
        return None;
    }

    // Unlike in xsf, recover signed infinities when the leading term of erfcx
    // certifies both signs. Use the original argument for the phase.
    let (hi, lo) = two_prod(-x, y);
    let phase = k * hi;
    let exponent = if y_abs >= 64.0 && x_abs <= y_abs * 0.5 {
        if phase.abs() < 0.5 {
            // Integrate from 0 to iy, then to z: the small phase keeps the horizontal
            // cosine positive and its sine term below the vertical integral.
            return Some(c64(-f64::INFINITY.copysign(x), -f64::INFINITY.copysign(y)));
        }
        // In this wedge the exponent is at least 1536k; the sign test suffices.
        f64::INFINITY
    } else {
        let exponent = 0.5 * k * (y_abs - x_abs) * (y_abs + x_abs);
        // A positive clearance is below one, and ln(|y|) > 3 whenever
        // the exponent exceeds LOG_OVERFLOW_BOUND.
        if exponent <= LOG_OVERFLOW_BOUND + 3.0 {
            return None;
        }
        exponent
    };
    let (quadrant, r) = if phase.abs() > pow2(53) {
        phase::reduce_product_phase((-0.5 * k) * x, y)
    } else {
        // Keep the rounded quadrant within the exact-integer range of f64.
        recentered_phase(phase, k * lo)
    };
    if r.abs() > 1.0 {
        return None;
    }

    // Fifth-degree sine and fourth-degree cosine bound the component
    // error by |r|^6/600 for |r| <= 1 and |x/y| <= 1.
    let r2 = r * r;
    let sin = r * (1.0 - r2 * (1.0 / 6.0 - r2 / 120.0));
    let cos = 1.0 - r2 * (0.5 - r2 / 24.0);
    let (re, im) = quadrant.rotate(sin, cos);

    let leading = c64(x / y_abs, -y.signum()) * c64(re, im);
    // The normalized erfcx remainder is bounded by 6/(k*y^2).
    // Allow for phase reduction and floating-point evaluation.
    let error = r2 * r2 * r2 / 600.0 + (8.0 / k / y_abs) / y_abs + pow2(-44);
    let clearance = leading.re.abs().min(leading.im.abs()) - error;
    // The exact scale is at least sqrt(k/2) * exp(exponent) / (2*sqrt(pi)*|y|),
    // which exceeds exp(exponent) / (6*|y|).
    // Positive clearance gives true normalized magnitudes above 2^-45, so
    // exponent >= 1536 suffices even at |y| = MAX.
    if clearance > 0.0
        && (exponent >= 1536.0 || exponent - y_abs.ln() + clearance.ln() > LOG_OVERFLOW_BOUND)
    {
        Some(c64(f64::INFINITY.copysign(leading.re), f64::INFINITY.copysign(leading.im)))
    } else {
        None
    }
}

fn phase_rotation(quadrant: Quadrant, r: f64) -> Option<(C64, f64)> {
    if r.abs() > 1.0 {
        return None;
    }
    let r2 = r * r;
    let sin = r * (1.0 - r2 / 6.0);
    let cos = 1.0 - 0.5 * r2;
    let (re, im) = quadrant.rotate(sin, cos);
    Some((c64(re, im), r2 * r2 / 20.0))
}

fn recentered_phase(phase: f64, lo: f64) -> (Quadrant, f64) {
    let quadrant = (phase * FRAC_2_PI).round();
    let r = (-quadrant).mul_add(FRAC_PI_2_DD.hi, phase) + (-quadrant).mul_add(FRAC_PI_2_DD.lo, lo);
    // Unlike in xsf, recenter if the rounded product picked a neighboring quadrant.
    if r.abs() > 1.0 {
        let step = 1.0_f64.copysign(r);
        (Quadrant::from_index((quadrant + step) as i64 as u64), r - step * FRAC_PI_2_DD.hi)
    } else {
        (Quadrant::from_index(quadrant as i64 as u64), r)
    }
}

/// Error function $\erf(z)$ for real or complex input
///
/// This function is defined as [^DLMF]:
///
/// $$
/// \erf(z) = {2 \over \sqrt\pi} \int_0^z e^{-t^2} \dd t .
/// $$
///
/// Corresponds to [`scipy.special.erf`][erf].
///
/// [erf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/7.2.E1>.
///
/// # See also
/// - [`erfc`]: Complementary error function
/// - [`erfi`]: Imaginary error function
/// - [`erfinv`](crate::erfinv): Inverse error function
/// - [`erf_zeros`]: Zeros of the error function
#[inline]
pub fn erf<T: ErfArg>(z: T) -> T {
    z.erf()
}

/// Complementary error function $1 - \erf(z)$ for real or complex input
///
/// Corresponds to [`scipy.special.erfc`][erfc].
///
/// [erfc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfc.html
///
/// # See also
/// - [`erf`]: Error function
/// - [`erfcinv`](crate::erfcinv): Inverse of the complementary error function
/// - [`erfcx`]: Scaled complementary error function
#[inline]
pub fn erfc<T: ErfArg>(z: T) -> T {
    z.erfc()
}

/// Scaled complementary error function $e^{z^2} \erfc(z)$ for real or complex input
///
/// Corresponds to [`scipy.special.erfcx`][erfcx].
///
/// [erfcx]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcx.html
///
/// # See also
/// - [`erfc`]: Complementary error function
/// - [`erf`]: Error function
#[inline]
pub fn erfcx<T: ErfArg>(z: T) -> T {
    z.erfcx()
}

/// Imaginary error function $-i \erf(i z)$ for real or complex input
///
/// Corresponds to [`scipy.special.erfi`][erfi].
///
/// [erfi]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfi.html
///
/// # See also
/// - [`erf`]: Error function
#[inline]
pub fn erfi<T: ErfArg>(z: T) -> T {
    z.erfi()
}

/// Zeros of the error function $\erf(z)$ in the first quadrant of the complex plane
///
/// Zeros in the other quadrants can be obtained by using the symmetries $\erf(-z) = -\erf(z)$
/// and $\erf(\overline z) = \overline{\erf(z)}$.
///
/// Corresponds to [`scipy.special.erf_zeros`][erf_zeros].
///
/// [erf_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf_zeros.html
///
/// # Panics
/// - if `N` is too large to fit into a [`c_int`](core::ffi::c_int)
///
/// # See also
/// - [`erf`]: Error function
/// - [`erfinv`](crate::erfinv): Inverse error function
#[doc(alias = "cerzo")]
#[must_use]
#[inline]
pub fn erf_zeros<const N: usize>() -> [Complex64; N] {
    let mut zo = [C64::NAN_RE; N];
    crate::xsf::specfun::cerzo(N.try_into().unwrap(), &mut zo);
    zo.map(Into::into)
}

/// Dawson function $\sqrt\pi e^{-z^2} \erfi(z) / 2$ for real or complex input
///
/// Corresponds to [`scipy.special.dawsn`][dawsn].
///
/// [dawsn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.dawsn.html
///
/// # See also
/// - [`erfi`]: Imaginary error function
#[doc(alias = "dawson")]
#[must_use]
#[inline]
pub fn dawsn<T: ErfArg>(z: T) -> T {
    z.dawsn()
}

/// Faddeeva function $e^{-z^2} \erfc(-i z)$
///
/// Corresponds to [`scipy.special.wofz`][wofz].
///
/// [wofz]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wofz.html
///
/// # See also
/// - [`erfc`]: Complementary error function
/// - [`erfcx`]: Scaled complementary error function
/// - [`voigt_profile`]: Voigt profile, which is computed from this function
#[doc(alias = "faddeeva")]
#[must_use]
#[inline]
pub fn wofz(z: Complex64) -> Complex64 {
    z.via(crate::xsf::faddeeva::w)
}

/// Voigt profile
///
/// The convolution of the densities of a normal distribution with standard deviation `sigma` and a
/// Cauchy distribution with half-width at half-maximum `gamma`, evaluated at `x`.
///
/// Corresponds to [`scipy.special.voigt_profile`][voigt].
///
/// [voigt]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.voigt_profile.html
///
/// # See also
/// - [`wofz`]: Faddeeva function
#[must_use]
#[inline]
pub fn voigt_profile(x: f64, sigma: f64, gamma: f64) -> f64 {
    const INV_SQRT_2: f64 = FRAC_1_SQRT_2;

    if sigma == 0.0 {
        if gamma == 0.0 {
            if x.is_nan() {
                return x;
            }
            if x == 0.0 {
                return f64::INFINITY;
            }
            return 0.0;
        }
        return gamma / PI / (x * x + gamma * gamma);
    }
    if gamma == 0.0 {
        return 1.0 / SQRT_TAU_DOWN / sigma * (-(x / sigma) * (x / sigma) / 2.0).exp();
    }

    let zreal = x / sigma * INV_SQRT_2;
    let zimag = gamma / sigma * INV_SQRT_2;
    let z = c64(zreal, zimag);
    let w = crate::xsf::faddeeva::w(z);
    w.re / sigma / SQRT_TAU_DOWN
}
