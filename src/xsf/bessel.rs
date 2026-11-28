#![allow(clippy::cast_possible_wrap)]

use core::ops::{Add, Mul};

use core::f64::consts::{FRAC_2_PI, FRAC_PI_2, FRAC_PI_4, LN_2, PI, TAU};

use num_complex::Complex;
use num_traits::ToPrimitive;

use crate::utils::{EGAMMA, FRAC_PI_SQ_6, LN_PI, mul_exp, neg_if_odd, pow2};
use crate::xsf::airy::{airy, airy_scaled};
use crate::xsf::amos::{self, ierr_to_sferr, set_error_and_nan};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::consts::{MAXLOG, MINLOG};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, FRAC_PI_2_DD, LOG2, two_sum};
use crate::xsf::cephes::ikv_asymptotic_uniform_sums;
use crate::xsf::complex::Num as _;
use crate::xsf::expint::scaled_exp1;
use crate::xsf::phase::Quadrant;
use crate::xsf::specfun;
use crate::xsf::specfun::{SPECFUN_BIG, SPECFUN_INF};
use crate::xsf::trig::{cephes_cospi, cephes_sinpi};

// the arithmetic bounds are needed for `bessel_diff_formula` in `scipy_special::bessel_prime`, and
// `Into` for the Hankel functions, which are complex for real arguments too
pub trait BesselArg:
    crate::sealed::Sealed
    + Copy
    + From<f64>
    + Into<Complex<f64>>
    + Add<Output = Self>
    + Mul<Output = Self>
{
    fn bessel_j(self, v: f64) -> Self;
    fn bessel_je(self, v: f64) -> Self;
    fn bessel_y(self, v: f64) -> Self;
    fn bessel_ye(self, v: f64) -> Self;
    fn bessel_i(self, v: f64) -> Self;
    /// `bessel_i` as in xsf, where the overflow of a part makes both inf or NaN
    fn bessel_i_xsf(self, v: f64) -> Self;
    fn bessel_ie(self, v: f64) -> Self;
    fn bessel_k(self, v: f64) -> Self;
    fn bessel_ke(self, v: f64) -> Self;
}

impl BesselArg for f64 {
    #[inline]
    fn bessel_j(self, v: f64) -> Self {
        cyl_bessel_j_f64(v, self)
    }
    #[inline]
    fn bessel_je(self, v: f64) -> Self {
        cyl_bessel_je_f64(v, self)
    }
    #[inline]
    fn bessel_y(self, v: f64) -> Self {
        cyl_bessel_y_f64(v, self)
    }
    #[inline]
    fn bessel_ye(self, v: f64) -> Self {
        cyl_bessel_ye_f64(v, self)
    }
    #[inline]
    fn bessel_i(self, v: f64) -> Self {
        crate::xsf::cephes::iv(v, self)
    }
    #[inline]
    fn bessel_i_xsf(self, v: f64) -> Self {
        self.bessel_i(v)
    }
    #[inline]
    fn bessel_ie(self, v: f64) -> Self {
        cyl_bessel_ie_f64(v, self)
    }
    #[inline]
    fn bessel_k(self, v: f64) -> Self {
        cyl_bessel_k_f64(v, self)
    }
    #[inline]
    fn bessel_ke(self, v: f64) -> Self {
        cyl_bessel_ke_f64(v, self)
    }
}

impl BesselArg for Complex<f64> {
    #[inline]
    fn bessel_j(self, v: f64) -> Self {
        self.via(|z| {
            let j = cyl_bessel_j_c64(v, z);
            // Unlike in xsf, recover a result that AMOS underflows to 0.
            if j == C64::ZERO && z != C64::ZERO {
                bessel_j_series_underflow(v, z).unwrap_or(j)
            } else {
                j
            }
        })
    }
    #[inline]
    fn bessel_je(self, v: f64) -> Self {
        self.via(|z| cyl_bessel_je_c64(v, z))
    }
    #[inline]
    fn bessel_y(self, v: f64) -> Self {
        self.via(|z| {
            let y = cyl_bessel_y_c64(v, z);
            // Unlike in xsf, recover signed overflow at small complex arguments.
            if y.re.is_nan() && y.im.is_nan() {
                return bessel_small_argument_overflow(v, z, 1.0).map_or(y, |k| -k);
            }
            y
        })
    }
    #[inline]
    fn bessel_ye(self, v: f64) -> Self {
        self.via(|z| cyl_bessel_ye_c64(v, z))
    }
    #[inline]
    fn bessel_i(self, v: f64) -> Self {
        self.via(|z| cyl_bessel_i_amos(v, z, true))
    }
    #[inline]
    fn bessel_i_xsf(self, v: f64) -> Self {
        self.via(|z| cyl_bessel_i_c64(v, z))
    }
    #[inline]
    fn bessel_ie(self, v: f64) -> Self {
        self.via(|z| cyl_bessel_ie_c64(v, z))
    }
    #[inline]
    fn bessel_k(self, v: f64) -> Self {
        self.via(|z| {
            let k = cyl_bessel_k_c64(v, z);
            // Unlike in xsf, return signed overflow when both components can be resolved.
            if k.re.is_nan() && k.im.is_nan() {
                return bessel_small_argument_overflow(v.abs(), z, 1.0).unwrap_or(k);
            }
            k
        })
    }
    #[inline]
    fn bessel_ke(self, v: f64) -> Self {
        self.via(|z| cyl_bessel_ke_c64(v, z))
    }
}

// Translated from `xsf/cpu/bessel.h` (scipy/xsf@0f44f5d)

fn rotate(z: C64, v: f64) -> C64 {
    let c = cephes_cospi(v);
    let s = cephes_sinpi(v);

    c64(c * z.re - s * z.im, s * z.re + c * z.im)
}

fn rotate_jy(j: C64, y: C64, v: f64) -> C64 {
    let c = cephes_cospi(v);
    let s = cephes_sinpi(v);

    if c == 0.0 {
        // half-integer v: c * j is 0, but NaN if j overflowed
        return -s * y;
    }
    c * j - s * y
}

fn reflect_jy(jy: &mut C64, v: f64) -> bool {
    /* NB: Y_v may be huge near negative integers -- so handle exact
     *     integers carefully
     */
    if v != v.floor() {
        return false;
    }

    let i = (v - 16384.0 * (v / 16384.0).floor()) as i32;
    *jy = neg_if_odd(i, *jy);
    true
}

fn reflect_i(_ik: &mut C64, v: f64) -> bool {
    v == v.floor() /* I is symmetric for integer v */
}

fn rotate_i(i: C64, k: C64, v: f64) -> C64 {
    // Unlike in xsf (scipy/xsf#383), sinpi(v) instead of the inaccurate sin(pi v).
    let s = cephes_sinpi(v) * FRAC_2_PI;

    i + s * k
}

/// The signature of `amos::besi`, `besj`, `besk` and `besy`.
type AmosBessel = fn(C64, f64, i32, i32, &mut [C64], &mut i32) -> i32;

/// Returns the order `v` result of the AMOS routine `bes` and its `ierr`, with the result set to
/// NaN on the errors of `set_error_and_nan`.
fn amos1(bes: AmosBessel, z: C64, v: f64, kode: i32) -> (C64, i32) {
    let n = 1;
    let mut ierr = 0;
    let mut cy = C64::NAN;
    let nz = bes(z, v, kode, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    (cy, ierr)
}

/// Whether AMOS overflowed (`ierr = 2`) on the real axis `z >= 0`, where the result is `+-inf`.
fn overflowed_real(ierr: i32, z: C64) -> bool {
    ierr == 2 && z.re >= 0.0 && z.is_real()
}

/// Whether `J_v(x) = 0` and `Y_v(x) = -inf` in f64, for real `z = x` and a finite order `v` above
/// the AMOS limit of 2^51.
///
/// With `w = sqrt(1 - (x / v)^2)` and `t = v (atanh(w) - w) >= v w^3 / 3`, `|J_v(x)| <= e^-t`
/// (DLMF 10.14.5) and `Y_v(x) ~ -e^t / sqrt(pi v w / 2)` (DLMF 10.19.3), where
/// `e^t > f64::MAX^2` leaves room for the square root.
fn jy_huge_order(v: f64, z: C64) -> bool {
    let s = z.re.abs() / v;
    let w2 = (1.0 - s) * (1.0 + s);
    let huge = v > 0.5 / f64::EPSILON && v.is_finite();
    huge && z.is_real() && v * w2 * w2.sqrt() / 3.0 > 2.0 * MAXLOG
}

/// `amos1` of `besj`. Unlike in xsf, 0 instead of NaN for a `jy_huge_order`, with the sign of
/// `J_v(-x) = (-1)^v J_v(x)` for an integer `v`.
fn amos_j(z: C64, v: f64, kode: i32) -> (C64, i32) {
    if jy_huge_order(v, z) {
        let sign = if v % 2.0 == 1.0 { z.re } else { 1.0 };
        return (c64(0.0_f64.copysign(sign), 0.0), 0);
    }
    // Unlike in xsf, J_v(x) = Re H1_v(x) (DLMF 10.4.3) for real x > fnul >= v, where besj
    // can recur down from the order fnul (85.9) and lose up to 3e-12 of sqrt(J^2 + Y^2),
    // but besh does not.
    //
    // The scaling e^-|Im z| is 1 here, so this ignores `kode`.
    let (.., fnul) = amos::MACH;
    if z.is_real() && z.re > fnul && v <= fnul {
        let besh1 =
            |z, v, kode, n, cy: &mut [C64], ierr: &mut i32| amos::besh(z, v, kode, 1, n, cy, ierr);
        let (h, ierr) = amos1(besh1, z, v, 1);
        return (c64(h.re, 0.0), ierr);
    }
    amos1(amos::besj, z, v, kode)
}

/// `amos1` of `besy`. Unlike in xsf, an overflow (`ierr = 2`) instead of `ierr = 4` for a
/// `jy_huge_order`.
fn amos_y(z: C64, v: f64, kode: i32) -> (C64, i32) {
    if jy_huge_order(v, z) {
        return (C64::NAN, 2);
    }
    amos1(amos::besy, z, v, kode)
}

fn cyl_bessel_ie_c64(mut v: f64, z: C64) -> C64 {
    let kode = 2;
    let mut sign = 1;

    if v.is_nan() || z.is_nan() {
        return C64::NAN;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let (mut cy, _) = amos1(amos::besi, z, v, kode);

    if sign == -1 && !reflect_i(&mut cy, v) {
        let (mut cy_k, ierr) = amos1(amos::besk, z, v, kode);
        // Unlike in xsf, scale the reflection term before the K recurrence if K overflowed.
        if ierr == 2 {
            if let Some(value) = i_reflection_from_lower_orders(cy, v, z, kode) {
                return value;
            }
        }
        /* adjust scaling to match zbesi */
        cy_k *= C64::cis(z.im).conj();
        if z.re > 0.0 {
            // Unlike in xsf, apply e^-2x as (e^-x)^2 if it underflows, as the result may not.
            let e2 = (-2.0 * z.re).exp();
            cy_k = if e2 >= f64::MIN_POSITIVE {
                cy_k * e2
            } else {
                let e = (-z.re).exp();
                cy_k * e * e
            };
        }
        // Unlike in xsf (scipy/xsf#383), an AMOS overflow of K_v(x) is inf instead of NaN.
        if overflowed_real(ierr, z) {
            cy_k = C64::INFINITY;
        }
        /* v -> -v */
        cy = rotate_i(cy, cy_k, v);
    }

    cy
}

fn cyl_bessel_ie_f64(v: f64, x: f64) -> f64 {
    if v != v.floor() && x < 0.0 {
        return f64::NAN;
    }

    let ie = cyl_bessel_ie_c64(v, c64(x, 0.0)).re;
    // the e^-x K_{-v}(x) term of the lower orders is lost if exp(-2x) underflows
    if ie.is_infinite() && (-2.0 * x).exp() >= f64::MIN_POSITIVE {
        return cyl_bessel_overflow(cyl_bessel_ie_c64, 1.0, v, x, ie);
    }
    // Unlike in xsf (scipy/xsf#403), use cephes for e^-|x| I_v(x) if AMOS underflowed and e^-|x| is
    // normal.
    if amos_underflow(ie, v, x) {
        let scale = (-x.abs()).exp();
        if scale >= f64::MIN_POSITIVE {
            return crate::xsf::cephes::iv(v, x) * scale;
        }
    }
    ie
}

fn cyl_bessel_ke_c64(v: f64, z: C64) -> C64 {
    cyl_bessel_k_kode(v, z, 2)
}

pub(crate) fn cyl_bessel_ke_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return f64::INFINITY;
    }

    let ke = cyl_bessel_ke_c64(v, c64(x, 0.0)).re;
    cyl_bessel_overflow(cyl_bessel_ke_c64, 1.0, v.abs(), x, ke)
}

// AMOS overflows above about 1e304; recur from lower orders (DLMF 10.29.1).
// Unlike in xsf, retry on an infinite wrapper result, try starting more than one order lower,
// and share the recurrence with I, J, and Y (DLMF 10.6.1, 10.29.1).
fn cyl_bessel_overflow(f: fn(f64, C64) -> C64, c: f64, v: f64, x: f64, res: f64) -> f64 {
    if !res.is_infinite() {
        return res;
    }
    let (d, z) = (v.signum(), c64(x, 0.0));
    for m in [1, 2, 4, 8, 16_u8] {
        if f64::from(m) >= v.abs() {
            break;
        }
        let mut f1 = f(v - f64::from(m) * d, z).re;
        if !f1.is_finite() {
            continue;
        }
        let mut f0 = f(v - f64::from(m + 1) * d, z).re;
        for k in (1..=m).rev() {
            (f0, f1) = (f1, 2.0 * (v - f64::from(k) * d) / x * f1 + c * f0);
        }
        return if f1.is_finite() { f1 } else { res };
    }
    res
}

/// Whether a real AMOS J or I result underflowed to 0, which AMOS does early: from about 1e-290
/// for integer v, else 1e-304.
///
/// For negative non-integer v, a 0 is a cancellation in the reflection formula instead.
fn amos_underflow(res: f64, v: f64, x: f64) -> bool {
    res == 0.0 && x != 0.0 && (v >= 0.0 || v == v.floor())
}

/// The smallest component to recover when a rotation of `I_v` or its reflection term
/// has an error relative to the magnitude of the complex result.
const I_ROTATED_PART_MIN: f64 = 1.0 / 16.0;

#[cold]
fn i_reflection_from_lower_orders(i: C64, v: f64, z: C64, kode: i32) -> Option<C64> {
    // Keep the first small-argument correction |z|²/(4(v−1)) at most 1/2.
    if !(i.is_finite() && v >= 2.0 && z.re > 0.0 && !z.is_real() && z.norm_sqr() <= v) {
        return None;
    }
    // Use the order reflection (DLMF 10.27.2). The seeds are unscaled K,
    // so Ie only needs exp(-|Re z|).
    let scale = if kode == 2 { (-z.re).exp() } else { 1.0 };
    let weight = (FRAC_2_PI * cephes_sinpi(v)) * scale;
    let weight = weight.is_normal().then_some(weight)?;

    let evaluate = |offset: i32| {
        let mut ierr = 0;
        let mut seeds = [C64::NAN; 2];
        let nz = amos::besk(z, v - f64::from(offset), 1, 2, &mut seeds, &mut ierr);
        (ierr == 0 && nz == 0).then_some(seeds).ok_or(ierr)
    };
    // Unlike in xsf, try another pair of lower orders if the first pair overflows.
    let (offset, seeds) = match evaluate(2) {
        Ok(seeds) => (2, seeds),
        Err(2) if v >= 4.0 => {
            // The K integral (DLMF 10.32.9) bounds its relative displacement from K_v(Re z)
            // by |Im z| * (v/Re z + 1).
            // If I is also nearly real, their sum cannot pass the final imaginary-part check.
            let margin = I_ROTATED_PART_MIN / 2.0;
            if z.im.abs() * (v + z.re) <= margin * z.re && i.im.abs() <= margin * i.norm() {
                return None;
            }
            (4, evaluate(4).ok()?)
        },
        Err(_) => return None,
    };
    let weighted = seeds.map(|seed| seed * weight);
    let [lower, upper] = (2..offset).rev().try_fold(weighted, |[lower, upper], k| {
        let step = (upper / z) * (2.0 * (v - f64::from(k)));
        let next = lower + step;
        // Bound amplification from cancellation in each additional step.
        (next.norm() >= 0.5 * (lower.norm() + step.norm())).then_some([upper, next])
    })?;
    let step = (upper / z) * (2.0 * (v - 1.0));
    let value = i + lower + step;
    let amplitude = i.norm() + lower.norm() + step.norm();
    let minimum = (I_ROTATED_PART_MIN * amplitude).max(2.0 * f64::MIN_POSITIVE);
    let range = minimum..=0.5 * f64::MAX;
    (range.contains(&value.re.abs()) && range.contains(&value.im.abs())).then_some(value)
}

/// `I_v(z)` for `v >= 0` where AMOS overflowed, off the real axis or for a non-integer `v`.
///
/// xsf multiplies both parts of `e^-|Re z| I_v(z)` by inf, which also loses a part below
/// `f64::MAX`. So this unscales each normal part instead, unless the rotation for `Re z < 0`
/// cancelled it.
fn cyl_bessel_i_overflow(v: f64, z: C64) -> C64 {
    let ie = cyl_bessel_ie_c64(v, z);
    let cy = ie.mul_exp(z.re.abs());
    let min = if z.re < 0.0 {
        I_ROTATED_PART_MIN * ie.norm()
    } else {
        0.0
    };
    // as in xsf for a part that lost its digits to underflow or cancellation
    let part = |p: f64, y: f64| {
        if p.is_normal() && p.abs() >= min {
            y
        } else {
            p * f64::INFINITY
        }
    };
    c64(part(ie.re, cy.re), part(ie.im, cy.im))
}

pub(crate) fn cyl_bessel_i_c64(v: f64, z: C64) -> C64 {
    cyl_bessel_i_amos(v, z, false)
}

/// With `by_part`, unlike in xsf, a part below `f64::MAX` stays finite where the other overflows.
/// Only for the value of `bessel_i`, as a sum of such values (e.g. a derivative) can be inf where
/// xsf is NaN.
fn cyl_bessel_i_amos(mut v: f64, z: C64, by_part: bool) -> C64 {
    let kode = 1;
    let mut sign = 1.0;

    if v.is_nan() || z.is_nan() {
        return C64::NAN;
    }
    if v < 0.0 {
        v = -v;
        sign = -1.0;
    }
    let (mut cy, ierr) = amos1(amos::besi, z, v, kode);
    if ierr == 2 {
        /* overflow */
        if z.is_real() && (z.re >= 0.0 || v == v.floor()) {
            if z.re < 0.0 && v / 2.0 != (v / 2.0).floor() {
                cy.re = f64::NEG_INFINITY;
            } else {
                cy.re = f64::INFINITY;
            }
            cy.im = 0.0;
        } else if by_part {
            cy = cyl_bessel_i_overflow(v, z);
        } else {
            cy = cyl_bessel_ie_c64(v * sign, z);
            cy.re *= f64::INFINITY;
            cy.im *= f64::INFINITY;
        }
    }

    if sign == -1.0 && !reflect_i(&mut cy, v) {
        let (cy_k, ierr_k) = amos1(amos::besk, z, v, kode);
        // Unlike in xsf, scale the reflection term before the K recurrence if K overflowed.
        if ierr_k == 2 {
            if let Some(value) = i_reflection_from_lower_orders(cy, v, z, kode) {
                return value;
            }
            if by_part {
                // Unlike in xsf, recover both reflection signs when the weighted K overflows.
                // At these small arguments the positive-order I is bounded by 2.
                let weight = FRAC_2_PI * cephes_sinpi(v);
                if let Some(value) = bessel_small_argument_overflow(v, z, weight) {
                    return value;
                }
            }
        }
        cy = rotate_i(cy, cy_k, v);
    }

    cy
}

pub(crate) fn cyl_bessel_k_c64(v: f64, z: C64) -> C64 {
    cyl_bessel_k_kode(v, z, 1)
}

/// Unlike in xsf, `K_v` and `e^z K_v` share this implementation (xsf duplicates it), with
/// `kode = 2` for exponential scaling.
fn cyl_bessel_k_kode(mut v: f64, z: C64, kode: i32) -> C64 {
    if v.is_nan() || z.is_nan() {
        return C64::NAN;
    }

    if v < 0.0 {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }

    if cyl_bessel_k_overflows(v, z, kode == 2) {
        return C64::INFINITY;
    }

    let (mut cy, ierr) = amos1(amos::besk, z, v, kode);
    if overflowed_real(ierr, z) {
        /* overflow */
        cy = C64::INFINITY;
    }
    // AMOS flushes K to zero near Re(z) = 700, also where K is normal. Scale K exp(z) by
    // exp(-z) in two steps.
    // Unlike in xsf, check for a zero AMOS result and a finite scaled value in place of the flags.
    // Rotate first, use componentwise mul_exp with its log fallback, and use the uniform expansion
    // for real z if the scaled result overflows.
    if kode == 1 && cy == C64::ZERO {
        let w = cyl_bessel_k_kode(v, z, 2) * C64::cis(z.im).conj();
        if w.is_finite() {
            cy = c64(mul_exp(w.re, -z.re), mul_exp(w.im, -z.re));
        } else if z.is_real() {
            crate::xsf::cephes::ikv_asymptotic_uniform(v, z.re, None, Some(&mut cy.re));
        }
    }

    cy
}

/* AMOS rejects orders above 0.5 / tol (~2e15). For real z = x > 0, check if K_v(x)
 * (times exp(x) if scaled) overflows, with the lower bound
 * K_v(x) >= sqrt(pi / (2 v)) (2 v / (e x))^v exp(-x^2 / (4 (v - 1))),
 * from Jensen's inequality on https://dlmf.nist.gov/10.32.E10 and
 * https://dlmf.nist.gov/5.6.E1.
 */
fn cyl_bessel_k_overflows(v: f64, z: C64, scaled: bool) -> bool {
    let x = z.re;
    if v <= f64::from(i32::MAX) * 0.5 || !z.is_real() || !(x > 0.0) {
        return false;
    }
    let logk =
        v * ((v / x * 2.0).ln() - 1.0) - x / (v - 1.0) * x / 4.0 - (v / FRAC_PI_2).ln() / 2.0;
    logk + if scaled { x } else { 0.0 } > MAXLOG
}

pub(crate) fn cyl_bessel_k_f64(v: f64, z: f64) -> f64 {
    if z < 0.0 {
        return f64::NAN;
    }

    if z == 0.0 {
        return f64::INFINITY;
    }

    // Unlike in xsf, this is a strict bound instead of x > 710 (1 + |v|), where K_v(x) can be
    // subnormal: K_v(x) <= sqrt(pi / (2 x)) e^(v^2 / (2 x) - x) < 2^-1075, by
    // cosh(t) >= 1 + t^2 / 2 in https://dlmf.nist.gov/10.32.E9.
    if z - v / z * v / 2.0 > -MINLOG {
        /* Underflow. See uniform expansion https://dlmf.nist.gov/10.41
         * We are here working around a restriction in AMOS.
         */
        return 0.0;
    }

    let k = cyl_bessel_k_c64(v, c64(z, 0.0)).re;
    cyl_bessel_overflow(cyl_bessel_k_c64, 1.0, v.abs(), z, k)
}

fn cyl_bessel_je_c64(v: f64, z: C64) -> C64 {
    cyl_bessel_j_kode(v, z, 2)
}

fn cyl_bessel_je_f64(v: f64, x: f64) -> f64 {
    if v != v.floor() && x < 0.0 {
        return f64::NAN;
    }

    // Unlike in xsf, an expansion of J_v(x) for large |v| is used where one is accurate (for
    // real x, e^-|Im x| = 1).
    if let Some(j) = bessel_j_large_order(v, x) {
        return j;
    }

    let je = cyl_bessel_je_c64(v, c64(x, 0.0)).re;
    // Unlike in xsf (scipy/xsf#402), use the cephes fallback of J_v(x) if AMOS underflowed.
    if amos_underflow(je, v, x) {
        return cyl_bessel_j_f64(v, x);
    }
    cyl_bessel_overflow(cyl_bessel_je_c64, -1.0, v, x, je)
}

fn cyl_bessel_ye_c64(v: f64, z: C64) -> C64 {
    scaled_bessel(
        z,
        c64(-z.im.abs(), 0.0),
        |kode| cyl_bessel_y_kode(v, z, kode),
        || scaled_bessel_y_on_imaginary_axis(v, z).or_else(|| scaled_bessel_from_k(v, z, false)),
    )
}

// Unlike in xsf, recover integer-order Y on the imaginary axis when its
// scaled I component is below half the least subnormal (DLMF 10.27.8 and 10.27.11).
#[cold]
fn scaled_bessel_y_on_imaginary_axis(v: f64, z: C64) -> Option<C64> {
    let k = scaled_axis_k_when_i_rounds_to_zero(v, z)?;
    let quadrant = v.rem_euclid(4.0) as u32;
    let (re, im) = Quadrant::from_index(u64::from(quadrant)).rotate(0.0, neg_if_odd(quadrant, -k));
    let value = c64(re, im);
    Some(if z.im < 0.0 { value.conj() } else { value })
}

#[cold]
fn scaled_axis_k_when_i_rounds_to_zero(v: f64, z: C64) -> Option<f64> {
    let (n, y) = (v.abs(), z.im.abs());
    // Normal decay implies y < 710. The integral over [1, 2] alone gives
    // K_n(y) > exp(n - 4y)/2 (DLMF 10.32.9), so n >= 4096 overflows.
    if z.re != 0.0 || v.fract() != 0.0 || !(1.0..4096.0).contains(&n) {
        return None;
    }
    let decay = (-y).exp();
    if !decay.is_normal() {
        return None;
    }
    let (raw_k, k) = normal_scaled_axis_k(n, y, decay)?;
    // The Wronskian and K recurrence give 2n I_n(y) K_n(y) < 1 (DLMF 10.28.2, 10.29.1).
    // With n >= 1, a zero quotient leaves at least a factor of two for rounding in K.
    (decay / raw_k == 0.0).then_some(k)
}

fn normal_scaled_axis_k(n: f64, y: f64, decay: f64) -> Option<(f64, f64)> {
    let raw = cyl_bessel_k_f64(n, y);
    let scaled = (FRAC_2_PI * raw) * decay;
    scaled.is_normal().then_some((raw, scaled))
}

fn cyl_bessel_ye_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    // Unlike in xsf, an expansion of Y_v(x) for large |v| is used where one is accurate (for
    // real x, e^-|Im x| = 1).
    if let Some((_, y)) = bessel_jy_large_order(v, x) {
        return y;
    }

    let ye = cyl_bessel_ye_c64(v, c64(x, 0.0)).re;
    cyl_bessel_overflow(cyl_bessel_ye_c64, -1.0, v, x, ye)
}

fn scaled_bessel(
    z: C64,
    exponent: C64,
    evaluate: impl Fn(i32) -> C64,
    recover: impl FnOnce() -> Option<C64>,
) -> C64 {
    let value = evaluate(2);
    // Unlike in xsf, recover a failed scaled result from the unscaled kernel.
    if value.re.is_nan() && value.im.is_nan() && z.is_finite() && !z.is_real() {
        rescale_bessel(evaluate(1), exponent)
            .or_else(recover)
            .unwrap_or(value)
    } else {
        value
    }
}

#[cold]
fn rescale_bessel(raw: C64, exponent: C64) -> Option<C64> {
    // Halve before rotating so the sum of two components cannot overflow.
    if !raw.is_finite() {
        return None;
    }
    let raw = 0.5 * raw;
    let rotated = raw * C64::cis(exponent.im);
    // The rotation amplifies errors by at most sqrt(2) in max norm. Keep
    // components large enough to limit the subsequent relative error.
    let min = (raw.linf_norm() / 32.0).max(2.0 * f64::MIN_POSITIVE);
    if rotated.re.abs().min(rotated.im.abs()) < min {
        return None;
    }
    // Apply the exponential in two steps, sharing the factor between components.
    let factor = (0.5 * exponent.re).exp();
    let value = rotated * factor * factor * 2.0;
    // Leave headroom from the overflow and subnormal thresholds.
    let range = 2.0 * f64::MIN_POSITIVE..=0.5 * f64::MAX;
    (range.contains(&value.re.abs()) && range.contains(&value.im.abs())).then_some(value)
}

fn cyl_hankel_scaled_c64(v: f64, z: C64, m: i32) -> C64 {
    let sign = if m == 1 { 1.0 } else { -1.0 };
    if z.re == 0.0 && sign * z.im < 0.0 && v.fract() == 0.0 {
        // Unlike in xsf, evaluate both normal imaginary-axis projections separately.
        if let Some(value) = damped_hankel_axis_normal(v, z.im) {
            return value;
        }
    }
    scaled_bessel(
        z,
        c64(sign * z.im, -sign * z.re),
        |kode| cyl_hankel_c64(v, z, kode, m),
        || {
            if sign * z.im < 0.0 {
                // Unlike in xsf, recover a failed damped Hankel value on the imaginary axis.
                damped_hankel_on_imaginary_axis(v, z).or_else(|| scaled_bessel_from_k(v, z, true))
            } else {
                scaled_hankel_overflow(v, z)
            }
        },
    )
}

// The damped Hankel function on the imaginary axis has orthogonal I and K
// projections for integer order (DLMF 10.4.4 and 10.27.8).
fn damped_hankel_axis_normal(v: f64, y: f64) -> Option<C64> {
    let t = y.abs();
    let decay = (-t).exp();
    if !decay.is_normal() {
        return None;
    }
    let n = v.abs();
    let i = cyl_bessel_i_amos(n, c64(t, 0.0), false).re * decay;
    if !i.is_normal() {
        return None;
    }
    let (_, k) = normal_scaled_axis_k(n, t, decay)?;
    Some(damped_hankel_axis_projection(v, y, 2.0 * i, k))
}

// With 2 exp(-|Im z|) I_n rounded to zero (DLMF 10.4.4 and 10.27.8),
// the remaining projection is a quarter turn of the scaled Y value.
#[cold]
fn damped_hankel_on_imaginary_axis(v: f64, z: C64) -> Option<C64> {
    // 2 I_n K_n < 1/n. An order of at least two leaves the same
    // factor-of-two margin as the single I_n in the scaled Y check.
    if v.abs() < 2.0 {
        return None;
    }
    let k = scaled_axis_k_when_i_rounds_to_zero(v, z)?;
    Some(damped_hankel_axis_projection(v, z.im, 0.0, k))
}

fn damped_hankel_axis_projection(v: f64, y: f64, i: f64, k: f64) -> C64 {
    let quadrant = v.rem_euclid(4.0) as u32;
    let (re, im) = Quadrant::from_index(u64::from(quadrant)).rotate(neg_if_odd(quadrant, k), i);
    let value = c64(re, im);
    if y < 0.0 { value.conj() } else { value }
}

// Limit the relative error from cancellation in either projection.
const MAX_PROJECTION_RATIO: f64 = 8.0;

// Unlike in xsf, recover the damped Hankel and Y when their continued K term
// dominates I and the ordinary intermediate K values overflow. This
// uses the K continuations (DLMF 10.27.8 and 10.27.11) and
// H1 = 2J - H2 (DLMF 10.4.4).
#[cold]
fn scaled_bessel_from_k(v: f64, z: C64, hankel: bool) -> Option<C64> {
    let w = c64(z.im.abs(), z.re);
    let order = v.abs();
    let k = rotate(amos::k_with_decay(order, w)?, 0.5 * order) * FRAC_2_PI;
    // Unlike in xsf, apply the order reflection (DLMF 10.4.5 and 10.4.6)
    // to the recovered projection.
    let value = if hankel {
        let value = c64(0.0, -1.0) * k * C64::cis(-z.re);
        if v < 0.0 { rotate(value, order) } else { value }
    } else if v < 0.0 {
        let cos = cephes_cospi(order);
        // When the Y term vanishes, the existing J evaluation supplies the result.
        if cos == 0.0 {
            return None;
        }
        -cos * k
    } else {
        -k
    };
    let (re, im) = (value.re.abs(), value.im.abs());
    let value = if re.min(im) >= re.max(im) / MAX_PROJECTION_RATIO {
        value
    } else if hankel && v < 0.0 && order.fract() == 0.5 {
        half_integer_hankel_from_k(order, w)?
    } else {
        return None;
    };
    Some(if z.im < 0.0 { value } else { value.conj() })
}

// Unlike in xsf, use the exact half-integer K seed (DLMF 10.39.2) and
// recurrence (DLMF 10.29.1) when the continued Debye result loses a projection.
fn half_integer_hankel_from_k(order: f64, w: C64) -> Option<C64> {
    const STEP: i32 = 400;
    // Start at orders -1/2 and 1/2 with the 2/pi prefactor from the Hankel
    // connection formula already applied. Leave exp(-w) until the end.
    let seed = c64(FRAC_2_PI.sqrt(), 0.0) / w.sqrt();
    let (mut previous, mut current) = (seed, seed);
    let (mut previous_bound, mut bound) = (seed.norm(), seed.norm());
    let inverse = 1.0 / w;
    let inverse_norm = inverse.norm();
    let mut shift = 0.0;
    for n in 0..order as i32 {
        let coefficient = f64::from(2 * n + 1);
        (previous, current) = (current, (coefficient * inverse) * current + previous);
        (previous_bound, bound) =
            (bound, (coefficient * inverse_norm).mul_add(bound, previous_bound));
        if bound > const { pow2(STEP) } {
            let scale = const { pow2(-STEP) };
            previous *= scale;
            current *= scale;
            previous_bound *= scale;
            bound *= scale;
            shift += f64::from(STEP);
        }
    }

    // The two exp(-w) factors come from the K seed and the Hankel scaling.
    // The remaining phase comes from the connection and order reflection.
    let projected = rotate(current * C64::cis(-2.0 * w.im), 1.5 * order - 0.5);
    // The norm majorant tracks perturbations from all recurrence steps, including
    // two successive complex products in the projection and continuation.
    let error = 32.0 * f64::EPSILON * order * bound;
    if projected.re.abs().min(projected.im.abs()) < 1e9 * error {
        return None;
    }

    let exponent = LOG2 * shift - 2.0 * w.re;
    let value = projected.mul_exp(f64::from(exponent));
    let range = 2.0 * f64::MIN_POSITIVE..=0.5 * f64::MAX;
    (range.contains(&value.re.abs()) && range.contains(&value.im.abs())).then_some(value)
}

// Unlike in xsf, classify both components before evaluating the overflowing exponential.
#[cold]
fn scaled_hankel_overflow(v: f64, z: C64) -> Option<C64> {
    let w = c64(z.im.abs(), z.re);
    let (factor, exponent) = amos::k_debye_parts(v.abs(), w)?;
    let scaled_exponent = exponent + w;
    // Keep a factor-of-two margin beyond overflow for both projections.
    let limit = MAXLOG + LN_2;
    if exponent.re <= 0.0 || scaled_exponent.re <= limit {
        return None;
    }

    // H_2e in the lower half-plane has a factor 2i/pi (DLMF 10.27.8).
    // Conjugate for H_1e in the upper.
    let phase = factor * C64::cis(scaled_exponent.im);
    let value = rotate(C64::I * phase, 0.5 * v);
    let value = if z.im < 0.0 { value } else { value.conj() };
    let (re, im) = (value.re.abs(), value.im.abs());
    let minor = re.min(im);
    (value.is_finite()
        && minor >= re.max(im) / MAX_PROJECTION_RATIO
        && scaled_exponent.re + (FRAC_2_PI * minor).ln() > limit)
        .then(|| value * f64::INFINITY)
}

/// Unlike in xsf, the four Hankel functions share this implementation, with `kode = 2` for the
/// exponentially scaled ones and `m` the kind.
fn cyl_hankel_c64(mut v: f64, z: C64, kode: i32, m: i32) -> C64 {
    let n = 1;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy = C64::NAN;

    if v.is_nan() || z.is_nan() {
        return cy;
    }
    if kode == 1 && m == 2 && v == 0.0 && z == C64::ZERO {
        cy.re = f64::NAN;
        cy.im = f64::INFINITY;
        return cy;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let x = z.re.abs();
    let cut = z.re < 0.0;
    // Unlike in xsf, take integer orders on the real axis up to the turning point from J and Y,
    // as AMOS loses the small real part there. Both imaginary zeros take the upper lip, where
    // the real coefficients of H1 and H2 are -1 and 3 times the integer-order parity, so that
    // H = a J + i b Y.
    let axis = (kode == 1 && z.is_real() && x <= v && v == v.trunc()).then(|| {
        let parity = if (sign < 0) == cut {
            1.0
        } else {
            cephes_cospi(v)
        };
        let (a, b) = match (cut, m) {
            (false, 1) => (1.0, 1.0),
            (false, _) => (1.0, -1.0),
            (true, 1) => (-1.0, 1.0),
            (true, _) => (3.0, -1.0),
        };
        (parity * a, parity * b)
    });
    let jy = axis.and_then(|_| bessel_jy_large_order(v, x));
    if let (Some((a, b)), Some((j, y))) = (axis, jy) {
        // 3 J would round twice if J is subnormal
        if y.is_finite() && j.is_finite() && (j.is_normal() || a.abs() == 1.0) {
            return c64(a * j, b * y);
        }
    }
    let nz = amos::besh(z, v, kode, m, n, core::slice::from_mut(&mut cy), &mut ierr);
    // Unlike in xsf, recover overflowing unscaled Hankel components on the real axis.
    if ierr == 2 && kode == 1 && z.is_real() && x >= f64::MIN_POSITIVE && v >= 1.0 {
        let log_leading = bessel_y_axis_log_leading(v, x);
        // Strong saturation makes J round to zero and each nonzero Y component overflow.
        if log_leading > 2.0 * MAXLOG && x * x < 4.0 * (v + 1.0) {
            return hankel_axis_saturated(v, sign < 0, cut, m);
        }
        if x <= 1.0 && (sign > 0 && !cut || v == v.trunc()) {
            let order = if sign < 0 { -v } else { v };
            let y = cyl_bessel_y_f64(order, x);
            if (y.is_finite() && log_leading + 0.5 * x < MAXLOG - BESSEL_Y_OVERFLOW_LOG_MARGIN)
                || (y.is_infinite() && log_leading > MAXLOG + BESSEL_Y_OVERFLOW_LOG_MARGIN)
            {
                let j = cyl_bessel_j_f64(order, x);
                let re = if cut && m == 2 {
                    hankel_cut_triple_j(v, x, log_leading, j)
                } else {
                    (j.is_finite() && bessel_j_axis_rounding_is_clear(v, x, log_leading, j))
                        .then_some(j)
                };
                if let Some(j) = re {
                    let parity = if cut { cephes_cospi(v) } else { 1.0 };
                    let re = if cut && m == 1 { -j } else { j };
                    let im = if m == 1 { y } else { -y };
                    return c64(parity * re, parity * im);
                }
            }
        }
    }
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if sign == -1 {
        cy = rotate(cy, if m == 1 { v } else { -v });
    }
    let Some((a, _)) = axis else {
        return cy;
    };

    // Unlike in xsf, recover a small real component from J. When |Re| >= |Im|/2, the norm is at
    // most sqrt(5) times |Re|.
    if cy.is_finite() && cy.re.abs() < 0.5 * cy.im.abs() {
        let j = if v == 1.0 {
            bessel_j1(x)
        } else if let Some((j, _)) = jy.filter(|(j, _)| j.is_finite()) {
            j
        } else if let Some(j) = (v < LARGE_ORDER_MIN_V)
            .then(|| bessel_j_axis_wronskian(v, x, cy.im.abs()))
            // a subnormal J means a tiny x, where the J of AMOS is accurate
            .filter(|j| j.is_normal())
        {
            j
        } else {
            cyl_bessel_j_amos(v, x)
        };
        if j.is_finite() {
            cy.re = a * j;
        }
    }
    cy
}

/// `J_n(x)` for an integer `n >= 1` and `0 < x <= n`, from `|Y_n(x)|`
///
/// By the Wronskian `J_(n+1) Y_n - J_n Y_(n+1) = 2 / (pi x)` (DLMF 10.5.2), with `Y_(n+1) / Y_n`
/// recurred up from `Y_1 / Y_0` (DLMF 10.6.1) and `J_(n+1) / J_n` from its continued fraction
/// (DLMF 10.10.1). Up to the turning point, `Y_n < 0` and `J_(n+1) / J_n < 1 < Y_(n+1) / Y_n`, so
/// their difference does not cancel much.
fn bessel_j_axis_wronskian(n: f64, x: f64, y_abs: f64) -> f64 {
    let t = 2.0 / x;
    let mut s = bessel_y1(x) / bessel_y0(x);
    let mut k = 1.0;
    while k <= n {
        s = k * t - 1.0 / s;
        k += 1.0;
    }
    // J_n / J_(n+1) = b_1 - 1 / (b_2 - ...) with b_k = (n + k) t, by the modified Lentz method
    let mut f = (n + 1.0) * t;
    let (mut c, mut d) = (f, 0.0);
    let mut k = n + 2.0;
    loop {
        let b = k * t;
        d = 1.0 / (b - d);
        c = b - 1.0 / c;
        f *= c * d;
        // also on NaN, where `t` overflows
        if !((c * d - 1.0).abs() > f64::EPSILON) {
            break;
        }
        k += 1.0;
    }
    FRAC_2_PI / (x * (s - 1.0 / f)) / y_abs
}

// If J rounds to zero and every nonzero multiple of Y overflows, only
// the coefficients from DLMF 10.4.6 and 10.11.5 are needed.
fn hankel_axis_saturated(v: f64, neg_order: bool, cut: bool, m: i32) -> C64 {
    let component = |j, y: f64| {
        if y == 0.0 {
            0.0_f64.copysign(j)
        } else {
            f64::INFINITY.copysign(-y)
        }
    };
    let h = if m == 1 { 1.0 } else { -1.0 };
    match (cut, neg_order, m) {
        (false, false, _) => c64(0.0, f64::INFINITY.copysign(-h)),
        (true, true, 1) => c64(-0.0, f64::NEG_INFINITY),
        (true, true, _) => {
            let sin = cephes_sinpi(2.0 * v);
            let cos = cephes_cospi(2.0 * v);
            // The coefficient 2 + cos(2 pi v) is positive.
            c64(component(1.0, -sin), component(-sin, -cos))
        },
        (cut, _, _) => {
            let sin = cephes_sinpi(v);
            let cos = cephes_cospi(v);
            if cut {
                c64(component(-h * cos, h * sin), component(sin, h * cos))
            } else {
                c64(component(cos, -sin), component(h * sin, h * cos))
            }
        },
    }
}

// On the upper cut, the real part of H2 at an integer order is ±3J. Round
// the product directly if the scalar J is subnormal.
fn hankel_cut_triple_j(v: f64, x: f64, log_y_leading: f64, j: f64) -> Option<f64> {
    use crate::utils::{next_down_nonneg as down, next_up_nonneg as up};

    const INV_MIN_NORMAL: f64 = 1.0 / f64::MIN_POSITIVE;
    const ROUNDING_MARGIN: f64 = f64::EPSILON * f64::EPSILON * 16_777_216.0;
    let result = 3.0 * j;
    if j.is_normal()
        || (result == 0.0
            && -log_y_leading - LN_PI - v.ln() + 3.0_f64.ln()
                < -1075.0 * LN_2 - BESSEL_Y_OVERFLOW_LOG_MARGIN)
    {
        return Some(result);
    }
    if !result.is_finite() || v >= 200.0 || x < 2.0 * f64::MIN_POSITIVE {
        return None;
    }

    let xh = 0.5 * x;
    let mut leading = Dd::new(3.0 * INV_MIN_NORMAL);
    for k in 1..=v as u32 {
        leading = leading * xh / f64::from(k);
        if leading.hi < 0.25 * f64::EPSILON {
            return None;
        }
    }
    // In units of MIN_POSITIVE, bound the first correction by the alternating
    // series and allow 2^-80 relative error for normal double-double operations.
    let error = leading * ROUNDING_MARGIN;
    let correction = leading * up(up(xh * xh) / (v + 1.0));
    let lo = leading - error - correction;
    let hi = leading + error;
    let center = (leading - correction * 0.5).hi * f64::MIN_POSITIVE;
    let midpoint =
        |a: f64, b: f64| (Dd::new(a * INV_MIN_NORMAL) + b * INV_MIN_NORMAL).mul_pwr2(0.5);
    [center, down(center), up(center)]
        .into_iter()
        .find(|&y| lo > midpoint(down(y), y) && hi < midpoint(y, up(y)))
        .map(|y| y.copysign(j))
}

pub(crate) fn cyl_bessel_j_c64(v: f64, z: C64) -> C64 {
    cyl_bessel_j_kode(v, z, 1)
}

/// Unlike in xsf, `J_v` and `e^-|Im z| J_v` share this implementation (xsf duplicates it), with
/// `kode = 2` for exponential scaling.
fn cyl_bessel_j_kode(mut v: f64, z: C64, kode: i32) -> C64 {
    let mut sign = 1;

    if v.is_nan() || z.is_nan() {
        return C64::NAN;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let (mut cy_j, ierr) = amos_j(z, v, kode);
    if kode == 1 && ierr == 2 {
        /* overflow */
        cy_j = cyl_bessel_j_kode(v, z, 2);
        cy_j.re *= f64::INFINITY;
        cy_j.im *= f64::INFINITY;
    }

    if sign == -1 && !reflect_jy(&mut cy_j, v) {
        let (mut cy_y, ierr) = amos_y(z, v, kode);
        if overflowed_real(ierr, z) {
            cy_y = C64::NEG_INFINITY;
        }
        cy_j = rotate_jy(cy_j, cy_y, v);
    }
    cy_j
}

/// The `|v|` from which the expansions for large orders are used for `J_v(x)` and `Y_v(x)`; below
/// it AMOS is accurate to a few eps of the amplitude, and faster
const LARGE_ORDER_MIN_V: f64 = 10.0;

/// Whether Debye's expansion in `bessel_jy_debye` is accurate to a few eps
///
/// Measured: with `w^2 = x^2 - v^2`, its 11 terms are for `1/w + 3 v^2/w^3 < 1/70`.
fn debye_is_accurate(v: f64, x: f64) -> bool {
    let cos2 = (v / x).powi(2);
    let sin2 = 1.0 - cos2;
    v.abs() >= LARGE_ORDER_MIN_V && x * sin2 * sin2.sqrt() > 70.0 * (1.0 + 2.0 * cos2)
}

/// `(J_v(x), Y_v(x))` by Debye's expansions (DLMF 10.19.6), for `debye_is_accurate(v, x)`
///
/// With `x = |v| sec(beta)`, `w = |v| tan(beta)` and `S = sum_k U_k(i cot(beta)) / |v|^k`, this is
/// `J_v - i Y_v = sqrt(2 / (pi w)) S e^(-i xi)` for `xi = w - |v| beta - pi/4`. For `v < 0`,
/// `beta - pi` replaces `beta` in `xi` (DLMF 10.4.6). The phase is in double-double, as
/// its rounding error is the error of `J` and `Y` relative to their amplitude.
fn bessel_jy_debye(v: f64, x: f64) -> (f64, f64) {
    let (order, v) = (v, v.abs());
    // t = tan(beta / 2), u = tan(pi / 4 - beta / 2)
    let sum_xv = Dd::new(x) + v;
    let t = ((Dd::new(x) - v) / sum_xv).sqrt();
    let u = (1.0 - t) / (t + 1.0);
    let w = t.hi * sum_xv.hi;
    let mut beta = FRAC_PI_2_DD - u.atan().mul_pwr2(2.0);
    if order < 0.0 {
        beta = beta - FRAC_PI_2_DD.mul_pwr2(2.0);
    }
    let sum = ikv_asymptotic_uniform_sums(v, c64(0.0, v / w)).0;
    // xi - arg(S) = x - shift, as w = x - |v| u
    let shift = (u + beta) * v + FRAC_PI_2_DD.mul_pwr2(0.5) + sum.im.atan2(sum.re);
    let (hi, lo) = two_sum(x, -shift.hi);
    let (sin, cos) = Dd { hi, lo }.sin_cos();
    let (sin_lo, cos_lo) = shift.lo.sin_cos();
    let amp = (FRAC_2_PI / w).sqrt() * sum.norm();
    let (j, y) = (cos * cos_lo + sin * sin_lo, sin * cos_lo - cos * sin_lo);
    (amp * j, amp * y)
}

/// The `w` up to which `amos::unhj` sums power series (those are for `w^2 <= 1/4`)
const UNHJ_SERIES_MAX_W: f64 = 0.49;

/// The `eta` from which the sums of `amos::unhj` are accurate to a few eps beyond its power
/// series (measured), for `x <= v`; below it their terms cancel
const UNHJ_MIN_ETA: f64 = 4.0;

/// `UNHJ_MIN_ETA` for `x > v`, where fewer terms cancel
const UNHJ_MIN_ETA_OSC: f64 = 3.0;

/// The `w` up to which `eta` is summed by its series in `w`; the closed forms cancel to `~w^3`,
/// and beyond it `v` is small for any `eta < 2 MAXLOG`
const ETA_SERIES_MAX_W: f64 = 1e-3;

/// `(J_v(x), Y_v(x))` by the uniform expansions in Airy functions (DLMF 10.20.4, 10.20.5) with
/// `phi` and the sums of AMOS, or `None` where those are inaccurate
///
/// With `w = sqrt|1 - (x / v)^2|`, the Airy functions are at `+-t` for `eta = 2/3 t^(3/2)`:
/// `v (atanh(w) - w)` for `x <= v`, else `v (w - atan(w))` (DLMF 10.20.2, 10.20.3). `eta` is in
/// double-double, as it is the exponent or the phase of `J` and `Y`.
fn bessel_jy_airy(v: f64, x: f64) -> Option<(f64, f64)> {
    if !(v >= LARGE_ORDER_MIN_V && x > 0.0) {
        return None;
    }
    let u = (Dd::new(v) - x) / v * ((Dd::new(v) + x) / v);
    let w = if x == v { u } else { u.abs().sqrt() };
    let eta = if w.hi <= ETA_SERIES_MAX_W {
        // w^3 sum_k u^k / (2 k + 3), with u = 1 - (x / v)^2 = +-w^2
        let tail = u.hi * (0.2 + u.hi * (1.0 / 7.0 + u.hi / 9.0));
        u.abs() * w * (Dd::ONE / 3.0 + tail)
    } else if x <= v {
        ((w + 1.0) * v / x).log() - w
    } else {
        w - w.atan()
    } * v;
    let min_eta = if x <= v {
        UNHJ_MIN_ETA
    } else {
        UNHJ_MIN_ETA_OSC
    };
    // from eta = 2 MAXLOG on, J = 0 and Y = -inf
    if !(eta.hi < 2.0 * MAXLOG && (w.hi <= UNHJ_SERIES_MAX_W || eta.hi >= min_eta)) {
        return None;
    }
    let [mut phi, mut arg, mut zeta1, mut zeta2, mut asum, mut bsum] = [C64::ZERO; 6];
    let (z, (tol, ..)) = (c64(x, 0.0), amos::MACH);
    #[rustfmt::skip]
    amos::unhj(z, v, 0, tol, &mut phi, &mut arg, &mut zeta1, &mut zeta2, &mut asum, &mut bsum);
    let (phi, asum, bsum) = (phi.re, asum.re, bsum.re);
    // J (Y) from Ai (-Bi) and its derivative
    let sum = |f: f64, fp: f64| phi * (f * asum + fp * bsum);
    let t = (1.5 * eta.hi).cbrt().powi(2);
    if x <= v {
        // Ai(t) e^eta and Bi(t) e^-eta
        let [ai, aip, bi, bip] = airy_scaled(t);
        let j = sum(ai, aip) * (1.0 - eta.lo);
        let y = -sum(bi, bip) * (1.0 + eta.lo);
        return Some((j * (-eta.hi).exp(), mul_exp(y, eta.hi)));
    }
    // The Airy functions at dt - t, to first order in the rounding error dt of -t, as
    // f''(-t) = -t f(-t) for f = Ai, Bi.
    let s = Dd::new(t).sqrt();
    let dt = (s * t - eta * 1.5).hi / (1.5 * s.hi);
    let sum_shifted = |f: f64, fp: f64| sum(f + fp * dt, fp - t * f * dt);
    let [ai, aip, bi, bip] = airy(-t);
    Some((sum_shifted(ai, aip), -sum_shifted(bi, bip)))
}

/// `bessel_jy_airy`, recurred up in the order (DLMF 10.6.1) where its sums cancel
///
/// For `x > v`, `eta < UNHJ_MIN_ETA_OSC` there, and `eta` is convex in `v` with the slope
/// `-atan(w)`, so it is large enough from `m` orders below `v` on.
///
/// For `x < v`, `eta < UNHJ_MIN_ETA` there, and `w <= UNHJ_SERIES_MAX_W` from `m` orders below
/// `v` on. Only `Y` is stable upward, so `J` is `J_(v-m)` times the ratios `J_n / J_(n-1)` up to
/// `v`, by their continued fraction (DLMF 10.10.1).
fn bessel_jy_airy_recur(v: f64, x: f64) -> Option<(f64, f64)> {
    let jy = bessel_jy_airy(v, x);
    if jy.is_some() || !(v >= LARGE_ORDER_MIN_V && x > 0.0) {
        return jy;
    }
    // eta = v |w - a| (DLMF 10.20.2, 10.20.3)
    let w = (1.0 - (x / v).powi(2)).abs().sqrt();
    let a = if x > v { w.atan() } else { w.atanh() };
    let eta = v * (w - a).abs();
    let m = if x > v {
        ((UNHJ_MIN_ETA_OSC - eta) / a).ceil()
    } else if eta < MAXLOG {
        (v - x / (1.0 - UNHJ_SERIES_MAX_W.powi(2)).sqrt()).ceil()
    } else {
        // eta >= MAXLOG, where Y ~ e^eta is near overflow: left to AMOS
        return None;
    };
    if m < 1.0 {
        return None;
    }
    // one order lower where `m` rounds to the edge of the range of `bessel_jy_airy`
    let seeds = |m: f64| Some((m, bessel_jy_airy(v - m - 1.0, x)?, bessel_jy_airy(v - m, x)?));
    let (m, mut jy0, mut jy1) = seeds(m).or_else(|| seeds(m + 1.0))?;
    // exact, as are the orders up to v
    let mut n = v - m;
    let mut j = jy1.0;
    while n < v {
        let c = 2.0 * n / x;
        (jy0, jy1) = (jy1, (c * jy1.0 - jy0.0, c * jy1.1 - jy0.1));
        n += 1.0;
    }
    if x > v {
        return Some(jy1);
    }
    // r = J_n / J_(n-1) at n = v - m + k; the tail of its continued fraction after k terms
    // decays as ((1 - w) / (1 + w))^k
    let mut k = m + (-f64::EPSILON.ln() / (2.0 * a)).ceil();
    let mut r = 0.0;
    while k > 0.0 {
        r = 1.0 / (2.0 * (v - m + k) / x - r);
        if k <= m {
            j *= r;
        }
        k -= 1.0;
    }
    Some((j, jy1.1))
}

/// `(J_v(x), Y_v(x))` for `x > 0` by an expansion for large `|v|`, where one is accurate
fn bessel_jy_large_order(v: f64, x: f64) -> Option<(f64, f64)> {
    if debye_is_accurate(v, x) {
        return Some(bessel_jy_debye(v, x));
    }
    let (j, y) = bessel_jy_airy_recur(v.abs(), x)?;
    if v > 0.0 {
        return Some((j, y));
    }
    // DLMF 10.4.6; AMOS handles an overflow of Y
    let (s, c) = (cephes_sinpi(-v), cephes_cospi(-v));
    y.is_finite().then_some((c * j - s * y, s * j + c * y))
}

/// `J_v(x)` by `bessel_jy_large_order`, with `J_v(-x) = (-1)^v J_v(x)` for an integer `v`
/// (DLMF 10.11.1)
fn bessel_j_large_order(v: f64, x: f64) -> Option<f64> {
    let (j, _) = bessel_jy_large_order(v, x.abs())?;
    Some(if x < 0.0 && v % 2.0 != 0.0 { -j } else { j })
}

pub(crate) fn cyl_bessel_j_f64(v: f64, x: f64) -> f64 {
    if v != v.floor() && x < 0.0 {
        return f64::NAN;
    }
    // Unlike in xsf, an expansion for large |v| is used where one is accurate, as the phase and
    // the exponent of AMOS lose ~eps x there.
    if let Some(j) = bessel_j_large_order(v, x) {
        return j;
    }
    cyl_bessel_j_amos(v, x)
}

fn cyl_bessel_j_amos(v: f64, x: f64) -> f64 {
    let res = cyl_bessel_j_c64(v, c64(x, 0.0));
    // Unlike in xsf (scipy/xsf#402), also if AMOS underflowed.
    if !res.re.is_finite() || amos_underflow(res.re, v, x) {
        /* AMOS returned NaN or inf, possibly due to overflow */
        let r = crate::xsf::cephes::jv(v, x);
        if !r.is_nan() {
            return r;
        }
    }

    res.re
}

pub(crate) fn cyl_bessel_y_c64(v: f64, z: C64) -> C64 {
    let mut y = cyl_bessel_y_kode(v, z, 1);
    // Unlike in xsf, use the real overflow recurrence on the positive real axis.
    let order = v.abs();
    if y.re.is_infinite()
        && z.is_real()
        && z.re > 0.0
        && z.re <= 1.0
        && order >= 1.0
        && (v > 0.0 || order == order.trunc())
        && bessel_y_axis_log_leading(order, z.re) + 0.5 * z.re
            < MAXLOG - BESSEL_Y_OVERFLOW_LOG_MARGIN
    {
        y.re = cyl_bessel_y_overflow(v, z.re, y.re);
    }
    // Unlike in xsf, recover finite integer-order values on the other axes by continuation.
    if y.re.is_nan() && y.im.is_nan() && order >= 1.0 && order == order.trunc() {
        return bessel_y_axis_continuation(v, z).unwrap_or(y);
    }
    y
}

// Both functions have a leading multiple of z^(-v) with opposite signs.
#[cold]
fn bessel_small_argument_overflow(v: f64, z: C64, weight: f64) -> Option<C64> {
    let m = z.linf_norm();
    let fraction = v.fract();
    let supported =
        matches!(fraction, 0.0 | 0.5) || v >= 1.0 && (0.125..=0.875).contains(&fraction);
    if v <= 0.0 || !supported || !z.is_finite() || z.re == 0.0 || z.im == 0.0 || m > 0.125 {
        return None;
    }
    let (im, re) = (-v * z.im.atan2(z.re)).sin_cos();
    // For half-integers, the relative tails are at most exp(2|z|)-1 (K)
    // and exp(|z|²/2)-1 (Y). For integers, the Laurent and logarithmic
    // tails (DLMF 10.8.1 and 10.31.1) are also less than 8m.
    // At fractional parts in [1/8, 7/8], every |k - v| is at least 1/8.
    // This and |sin(pi*v)| >= sin(pi/8) bound the series and connection terms
    // (DLMF 10.2.2–3, 10.25.2 and 10.27.4) by the same relative 8m
    // when v >= 1 and |z|/2 <= 1/8.
    let clearance = re.abs().min(im.abs()) - 8.0 * m - 16.0 * f64::EPSILON * v;
    if clearance <= 0.0 {
        return None;
    }
    // Gamma(v) m^(-v) / pi is below either leading magnitude since r <= 2m.
    let mut log_component = crate::gammaln(v) - v * m.ln() - LN_PI + clearance.ln();
    if weight != 1.0 {
        log_component += weight.abs().ln();
    }
    (log_component > MAXLOG + 1.0).then(|| {
        let value = c64(f64::INFINITY.copysign(re), f64::INFINITY.copysign(im));
        if weight < 0.0 { -value } else { value }
    })
}

/// `J_v(z) = (z/2)^v / Gamma(v + 1) 0F1(; v + 1; -z^2/4)` (DLMF 10.16.9) where the AMOS `seri`
/// sums it (`|z| <= 2` or `|z|^2/4 <= v + 1`), with `J_-n = (-1)^n J_n` (DLMF 10.4.1)
///
/// A subnormal result is rounded once: for an integer order up to the overflow of its
/// `Gamma(v + 1)`, the power is a product of `z 2^(k-1)` with `|z 2^(k-1)|^v <= Gamma(v + 1)`,
/// rescaled by `2^(-k v)`, so the parts that vanish on the axes and diagonals stay exactly 0.
/// Otherwise the power is scaled by `1 / eps` as in `seri`.
#[cold]
fn bessel_j_series_underflow(v: f64, z: C64) -> Option<C64> {
    const SCALE: f64 = 1.0 / f64::EPSILON;
    let order = v.abs();
    let az = z.norm();
    if !(az <= 2.0 || 0.25 * az * az <= order + 1.0) {
        return None;
    }
    let half_z = 0.5 * z;
    let sum: C64 = crate::hyp0f1(order + 1.0, Complex::from(-(half_z * half_z))).into();
    let gamma = crate::gamma(order + 1.0);
    let mut j = if order >= 1.0 && order == order.trunc() && gamma.is_finite() {
        let n = order as i32;
        let k = ((gamma.ln() / order - az.ln()) / LN_2 + 1.0).floor() as i32;
        (z.ldexp(k - 1).powi(n) * sum / gamma).ldexp(-k * n)
    } else {
        // as in `besj`, the upper lip of the cut for `Im z = -0`; the phase in units of pi keeps
        // the parts that vanish on the axes exactly 0
        let w = c64(z.re, z.im + 0.0);
        let ln = (w.ln().re - LN_2) * order - crate::gammaln(order + 1.0);
        rotate(sum, order * w.im.atan2(w.re) / PI) * (ln + SCALE.ln()).exp() / SCALE
    };
    // a non-integer `v < 0` also needs `Y_|v|` (DLMF 10.2.3)
    if v < 0.0 && !reflect_jy(&mut j, order) {
        return None;
    }
    (j.is_finite() && j != C64::ZERO).then_some(j)
}

fn bessel_y_axis_continuation(v: f64, z: C64) -> Option<C64> {
    let order = v.abs();
    let x = if z.is_real() && z.re < 0.0 {
        -z.re
    } else if z.re == 0.0 {
        z.im.abs()
    } else {
        return None;
    };
    if !(x > 0.0
        && x <= 1.0
        && bessel_y_axis_log_leading(order, x) + 0.5 * x < MAXLOG - BESSEL_Y_OVERFLOW_LOG_MARGIN)
    {
        return None;
    }
    let odd = order % 2.0 == 1.0;
    let value = if z.is_real() {
        // AMOS takes the upper lip for either imaginary zero (DLMF 10.11.6).
        let y = cyl_bessel_y_f64(v, x);
        let j = cyl_bessel_j_f64(v, x);
        let upper = c64(y, 2.0 * j);
        if odd { -upper } else { upper }
    } else {
        // DLMF 10.27.8 with H1 = J + iY. Scale the K seeds before raising the order.
        let seed = |n, z| c64(cyl_bessel_k_c64(n, z).re * FRAC_2_PI, 0.0);
        let k = cyl_bessel_overflow(seed, 1.0, order, x, seed(order, c64(x, 0.0)).re);
        let i = crate::xsf::cephes::iv(order, x);
        let upper = rotate(c64(-k, if odd { -i } else { i }), -0.5 * order);
        let value = if z.im < 0.0 { upper.conj() } else { upper };
        if v < 0.0 && odd { -value } else { value }
    };
    // The finite Y bound keeps both true components above half the smallest subnormal.
    (value.is_finite() && value.re != 0.0 && value.im != 0.0).then_some(value)
}

// This margin covers arithmetic rounding and the additive terms in the bounds.
const BESSEL_Y_OVERFLOW_LOG_MARGIN: f64 = 1e-8;

/// The log of the leading factor Gamma(v) (2/x)^v / pi, for v >= 1 and 0 < x <= 1.
/// DLMF 10.9.7 with u = exp(t) bounds |Y_v(x)| by exp(x/2) times the factor plus 1 + 1/(pi v).
fn bessel_y_axis_log_leading(v: f64, x: f64) -> f64 {
    crate::gammaln(v) + v * (LN_2 - x.ln()) - LN_PI
}

fn bessel_j_axis_rounding_is_clear(v: f64, x: f64, log_y_leading: f64, j: f64) -> bool {
    use crate::utils::{next_down_nonneg as down, next_up_nonneg as up};

    const LOG_HALF_MIN_SUBNORMAL: f64 = -1075.0 * LN_2;
    if j.is_normal() {
        return true;
    }
    // The leading terms of J and |Y| multiply to 1/(pi v).
    if j == 0.0 {
        return -log_y_leading - LN_PI - v.ln()
            < LOG_HALF_MIN_SUBNORMAL - BESSEL_Y_OVERFLOW_LOG_MARGIN;
    }
    if v != v.trunc() || v >= 200.0 {
        // Only use the product at integer orders. For v >= 200 and x <= 1
        // the leading term is already below half the smallest subnormal.
        return false;
    }
    // Bound the first term in units of MIN_POSITIVE, then include the first
    // correction as a lower bound for the alternating series (DLMF 10.2.2).
    let (x_lo, x_hi) = (down(0.5 * x), up(0.5 * x));
    let (mut lo, mut hi) = (1.0 / f64::MIN_POSITIVE, 1.0 / f64::MIN_POSITIVE);
    for k in 1..=v as u32 {
        let k = f64::from(k);
        lo = down(down(lo * x_lo) / k);
        hi = up(up(hi * x_hi) / k);
    }
    let correction = up(up(x * x) / (4.0 * (v + 1.0)));
    lo = down(lo * down(1.0 - correction));
    let magnitude = j.abs() / f64::MIN_POSITIVE;
    lo > magnitude - 0.5 * f64::EPSILON && hi < magnitude + 0.5 * f64::EPSILON
}

fn cyl_bessel_y_overflow(v: f64, x: f64, y: f64) -> f64 {
    cyl_bessel_overflow(|v, z| cyl_bessel_y_kode(v, z, 1), -1.0, v, x, y)
}

/// Unlike in xsf, `Y_v` and `e^-|Im z| Y_v` share this implementation (xsf duplicates it), with
/// `kode = 2` for exponential scaling.
fn cyl_bessel_y_kode(mut v: f64, z: C64, kode: i32) -> C64 {
    let mut sign = 1;
    let mut cy_y = C64::NAN;

    if v.is_nan() || z.is_nan() {
        return cy_y;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }

    if kode == 1 && z == C64::ZERO {
        /* overflow */
        cy_y.re = f64::NEG_INFINITY;
        cy_y.im = 0.0;
    } else {
        let ierr;
        (cy_y, ierr) = amos_y(z, v, kode);
        if overflowed_real(ierr, z) {
            /* overflow */
            cy_y.re = f64::NEG_INFINITY;
            cy_y.im = 0.0;
        }
    }

    if sign == -1 && !reflect_jy(&mut cy_y, v) {
        let (cy_j, _) = amos_j(z, v, kode);
        cy_y = rotate_jy(cy_y, cy_j, -v);
    }
    cy_y
}

pub(crate) fn cyl_bessel_y_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    // Unlike in xsf, an expansion for large |v| is used where one is accurate, as for J.
    if let Some((_, y)) = bessel_jy_large_order(v, x) {
        return y;
    }

    // Unlike in xsf, apply the shared overflow recovery once, to the raw result.
    let res = cyl_bessel_y_kode(v, c64(x, 0.0), 1);
    if res.re.is_nan() {
        return crate::xsf::cephes::yv(v, x);
    }

    cyl_bessel_y_overflow(v, x, res.re)
}

// Bessel J

/// Bessel function of the first kind, $J_v(z)$
///
/// Corresponds to [`scipy.special.jv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html
///
/// # See also
/// - [`bessel_j0`]: faster version of $J_0$ for real argument
/// - [`bessel_j1`]: faster version of $J_1$ for real argument
/// - [`bessel_je`]: $J_v$ without the leading exponential factor
/// - [`bessel_j_prime`](crate::bessel_j_prime): derivative $J_v\'(z)$
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "jv", alias = "cyl_bessel_j")]
#[must_use]
#[inline]
pub fn bessel_j<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_j(v)
}

/// Bessel function of the first kind of order 0, $J_0(x)$
///
/// Corresponds to [`scipy.special.j0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j0.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$ of real order and complex argument
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "j0", alias = "cyl_bessel_j0")]
#[must_use]
#[inline]
pub fn bessel_j0(x: f64) -> f64 {
    crate::xsf::cephes::j0(x)
}

/// Bessel function of the first kind of order 1, $J_1(x)$
///
/// Corresponds to [`scipy.special.j1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j1.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$ of real order and complex argument
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "j1", alias = "cyl_bessel_j1")]
#[must_use]
#[inline]
pub fn bessel_j1(x: f64) -> f64 {
    crate::xsf::cephes::j1(x)
}

/// Exponentially scaled Bessel function of the first kind, $e^{-\abs{\mathrm{Im}(z)}} J_v(z)$
///
/// Corresponds to [`scipy.special.jve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jve.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$
#[doc(alias = "jve", alias = "cyl_bessel_je")]
#[must_use]
#[inline]
pub fn bessel_je<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_je(v)
}

// Bessel Y (Neumann)

/// Bessel function of the second kind of order 0, $Y_0(x)$
///
/// Corresponds to [`scipy.special.y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y0.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$ of real order and complex argument
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "y0", alias = "cyl_bessel_y0", alias = "cyl_neumann_0")]
#[must_use]
#[inline]
pub fn bessel_y0(x: f64) -> f64 {
    crate::xsf::cephes::y0(x)
}

/// Bessel function of the second kind of order 1, $Y_1(x)$
///
/// Corresponds to [`scipy.special.y1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$ of real order and complex argument
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "y1", alias = "cyl_bessel_y1", alias = "cyl_neumann_1")]
#[must_use]
#[inline]
pub fn bessel_y1(x: f64) -> f64 {
    crate::xsf::cephes::y1(x)
}

/// Bessel function of the second kind, $Y_v(z)$
///
/// Corresponds to [`scipy.special.yv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yv.html
///
/// # See also
/// - [`bessel_y0`]: faster version of $Y_0$ for real argument
/// - [`bessel_y1`]: faster version of $Y_1$ for real argument
/// - [`bessel_ye`]: $Y_v$ without the leading exponential factor
/// - [`bessel_y_prime`](crate::bessel_y_prime): derivative $Y_v\'(z)$
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "yv", alias = "cyl_neumann", alias = "cyl_bessel_y")]
#[must_use]
#[inline]
pub fn bessel_y<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_y(v)
}

/// Exponentially scaled Bessel function of the second kind, $e^{-\abs{\mathrm{Im}(z)}} Y_v(z)$
///
/// Corresponds to [`scipy.special.yve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yve.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$
#[doc(alias = "yve", alias = "cyl_neumann_e", alias = "cyl_bessel_ye")]
#[must_use]
#[inline]
pub fn bessel_ye<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ye(v)
}

// Bessel I

/// Modified Bessel function of the first kind of order 0, $I_0(x)$
///
/// Corresponds to [`scipy.special.i0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$ of real order and complex argument
/// - [`bessel_i0e`]: exponentially scaled $I_0$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "i0", alias = "cyl_bessel_i0")]
#[must_use]
#[inline]
pub fn bessel_i0(x: f64) -> f64 {
    crate::xsf::cephes::i0(x)
}

/// Exponentially scaled modified Bessel function of the first kind of order 0, $e^{-|x|} I_0(x)$
///
/// Corresponds to [`scipy.special.i0e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0e.html
///
/// # See also
/// - [`bessel_i0`]: modified Bessel function $I_0(x)$
/// - [`bessel_ie`]: exponentially scaled $I_v(z)$
#[doc(alias = "i0e", alias = "cyl_bessel_i0e")]
#[must_use]
#[inline]
pub fn bessel_i0e(x: f64) -> f64 {
    crate::xsf::cephes::i0e(x)
}

/// Modified Bessel function of the first kind of order 1, $I_1(x)$
///
/// Corresponds to [`scipy.special.i1`][i1].
///
/// [i1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$ of real order and complex argument
/// - [`bessel_i1e`]: exponentially scaled $I_1$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "i1", alias = "cyl_bessel_i1")]
#[must_use]
#[inline]
pub fn bessel_i1(x: f64) -> f64 {
    crate::xsf::cephes::i1(x)
}

/// Exponentially scaled modified Bessel function of the first kind of order 1, $e^{-|x|} I_1(x)$
///
/// Corresponds to [`scipy.special.i1e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1e.html
///
/// # See also
/// - [`bessel_i1`]: modified Bessel function $I_1(x)$
/// - [`bessel_ie`]: exponentially scaled $I_v(z)$
#[doc(alias = "i1e", alias = "cyl_bessel_i1e")]
#[must_use]
#[inline]
pub fn bessel_i1e(x: f64) -> f64 {
    crate::xsf::cephes::i1e(x)
}

/// Modified Bessel function of the first kind, $I_v(z)$
///
/// Corresponds to [`scipy.special.iv`][scipy].
///
/// Near overflow, also well below `f64::MAX`, a part of a complex result can be infinite, NaN or 0
/// although it is finite and nonzero. [`scipy.special.iv`][scipy] has the same cases, and more.
/// The scaled [`bessel_ie`] avoids most of them.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iv.html
///
/// # See also
/// - [`bessel_i0`]: faster version of $I_0$ for real argument
/// - [`bessel_i1`]: faster version of $I_1$ for real argument
/// - [`bessel_ie`]: $I_v$ without the leading exponential factor
/// - [`bessel_i_prime`](crate::bessel_i_prime): derivative $I_v\'(z)$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "iv", alias = "cyl_bessel_i")]
#[must_use]
#[inline]
pub fn bessel_i<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_i(v)
}

/// Exponentially scaled modified Bessel function of the first kind,
/// $e^{-\abs{\mathrm{Re}(z)}} I_v(z)$
///
/// Corresponds to [`scipy.special.ive`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ive.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$
#[doc(alias = "ive", alias = "cyl_bessel_ie")]
#[must_use]
#[inline]
pub fn bessel_ie<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ie(v)
}

// Bessel K

/// Modified Bessel function of the second kind of order 0, $K_0(x)$
///
/// Corresponds to [`scipy.special.k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$ of real order and complex argument
/// - [`bessel_k0e`]: exponentially scaled $K_0$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "k0", alias = "cyl_bessel_k0")]
#[must_use]
#[inline]
pub fn bessel_k0(x: f64) -> f64 {
    crate::xsf::cephes::k0(x)
}

/// Exponentially scaled modified Bessel function of the second kind of order 0, $e^x K_0(x)$
///
/// Corresponds to [`scipy.special.k0e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0e.html
///
/// # See also
/// - [`bessel_k0`]: modified Bessel function $K_0(x)$
/// - [`bessel_ke`]: exponentially scaled $K_v(z)$
#[doc(alias = "k0e", alias = "cyl_bessel_k0e")]
#[must_use]
#[inline]
pub fn bessel_k0e(x: f64) -> f64 {
    crate::xsf::cephes::k0e(x)
}

/// Modified Bessel function of the second kind of order 1, $K_1(x)$
///
/// Corresponds to [`scipy.special.k1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$ of real order and complex argument
/// - [`bessel_k1e`]: exponentially scaled $K_1$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "k1", alias = "cyl_bessel_k1")]
#[must_use]
#[inline]
pub fn bessel_k1(x: f64) -> f64 {
    crate::xsf::cephes::k1(x)
}

/// Exponentially scaled modified Bessel function of the second kind of order 1, $e^x K_1(x)$
///
/// Corresponds to [`scipy.special.k1e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1e.html
///
/// # See also
/// - [`bessel_k1`]: modified Bessel function $K_1(x)$
/// - [`bessel_ke`]: exponentially scaled $K_v(z)$
#[doc(alias = "k1e", alias = "cyl_bessel_k1e")]
#[must_use]
#[inline]
pub fn bessel_k1e(x: f64) -> f64 {
    crate::xsf::cephes::k1e(x)
}

/// Modified Bessel function of the second kind, $K_v(z)$
///
/// Corresponds to [`scipy.special.kv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kv.html
///
/// # See also
/// - [`bessel_k0`]: faster version of $K_0$ for real argument
/// - [`bessel_k1`]: faster version of $K_1$ for real argument
/// - [`bessel_ke`]: $K_v$ without the leading exponential factor
/// - [`bessel_k_prime`](crate::bessel_k_prime): derivative $K_v\'(z)$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "kv", alias = "cyl_bessel_k")]
#[must_use]
#[inline]
pub fn bessel_k<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_k(v)
}

/// Exponentially scaled modified Bessel function of the second kind, $e^z K_v(z)$
///
/// Corresponds to [`scipy.special.kve`][kve].
///
/// [kve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kve.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$
#[doc(alias = "kve", alias = "cyl_bessel_ke")]
#[must_use]
#[inline]
pub fn bessel_ke<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ke(v)
}

// Zeros

/// Compute $N$ zeros of Bessel functions $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$
///
/// Corresponds to [`scipy.special.jnyn_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnyn_zeros.html
///
/// # Panics
/// - Panics if `N` is greater than 1200
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$
/// - [`bessel_y`]: Bessel function $Y_v(z)$
#[doc(alias = "jnyn_zeros")]
#[must_use]
#[inline]
pub fn bessel_zeros<const N: usize>(v: u32) -> [[f64; N]; 4] {
    assert!(N <= 1200, "N must be at most 1200");

    let (mut j, mut jp) = ([f64::NAN; N], [f64::NAN; N]);
    let (mut y, mut yp) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 || v > 100_000_000 {
        // Like in SciPy (where `nt` must be positive), `jyzo` isn't called for `N == 0`, since it
        // writes at least 1 zero.
        //
        // Unlike in SciPy, the same is done for huge `v`, to avoid the `int` overflows in `msta1`
        // and `msta2`, where the zeros are NaN anyway (as for `v >~ 5000`).
        return [j, jp, y, yp];
    }
    specfun::jyzo(v.try_into().unwrap(), N.try_into().unwrap(), &mut j, &mut jp, &mut y, &mut yp);
    [j, jp, y, yp]
}

/// The first `N` real or complex zeros of `Y0`, `Y1` (`kf = 1`) or `Y1'` (`kf = 2`), and the
/// values of `Y1`, `Y0` or `Y1` there
fn cyzo<const N: usize>(kf: i32, complex: bool) -> [[Complex<f64>; N]; 2] {
    let (mut zo, mut zv) = ([C64::ZERO; N], [C64::ZERO; N]);
    specfun::cyzo(kf, (!complex).into(), &mut zo, &mut zv);
    [zo.map(Into::into), zv.map(Into::into)]
}

/// Compute $N$ zeros $z_k$ of the Bessel function $Y_0(z)$, and the values $Y_1(z_k)$
///
/// Corresponds to [`scipy.special.y0_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y0_zeros.html
///
/// # Returns
///
/// - The real zeros $z_k > 0$ if `complex` is `false`, else the complex zeros with
///   $\mathrm{Re}\\, z_k < 0$ and $\mathrm{Im}\\, z_k > 0$ (their conjugates are zeros as well)
/// - The values $Y_1(z_k) = -Y_0\'(z_k)$ (the SciPy docs call them $Y_0\'(z_k)$)
///
/// # See also
/// - [`bessel_y1_zeros`]: zeros of $Y_1(z)$
/// - [`bessel_y1_prime_zeros`]: zeros of $Y_1\'(z)$
/// - [`bessel_zeros`]: real zeros of $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$
/// - [`bessel_y0`]: Bessel function $Y_0(x)$
#[doc(alias = "y0_zeros", alias = "cyzo")]
#[must_use]
#[inline]
pub fn bessel_y0_zeros<const N: usize>(complex: bool) -> [[Complex<f64>; N]; 2] {
    cyzo(0, complex)
}

/// Compute $N$ zeros $z_k$ of the Bessel function $Y_1(z)$, and the values $Y_1\'(z_k)$
///
/// Corresponds to [`scipy.special.y1_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1_zeros.html
///
/// # Returns
///
/// - The real zeros $z_k > 0$ if `complex` is `false`, else the complex zeros with
///   $\mathrm{Re}\\, z_k < 0$ and $\mathrm{Im}\\, z_k > 0$ (their conjugates are zeros as well)
/// - The values $Y_1\'(z_k) = Y_0(z_k)$
///
/// # See also
/// - [`bessel_y0_zeros`]: zeros of $Y_0(z)$
/// - [`bessel_y1_prime_zeros`]: zeros of $Y_1\'(z)$
/// - [`bessel_zeros`]: real zeros of $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$
/// - [`bessel_y1`]: Bessel function $Y_1(x)$
#[doc(alias = "y1_zeros", alias = "cyzo")]
#[must_use]
#[inline]
pub fn bessel_y1_zeros<const N: usize>(complex: bool) -> [[Complex<f64>; N]; 2] {
    cyzo(1, complex)
}

/// Compute $N$ zeros $z_k$ of the derivative $Y_1\'(z)$, and the values $Y_1(z_k)$
///
/// Corresponds to [`scipy.special.y1p_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1p_zeros.html
///
/// # Returns
///
/// - The real zeros $z_k > 0$ if `complex` is `false`, else the complex zeros with
///   $\mathrm{Im}\\, z_k > 0$, all but the first with $\mathrm{Re}\\, z_k < 0$ (their conjugates
///   are zeros as well)
/// - The values $Y_1(z_k)$
///
/// # See also
/// - [`bessel_y0_zeros`]: zeros of $Y_0(z)$
/// - [`bessel_y1_zeros`]: zeros of $Y_1(z)$
/// - [`bessel_zeros`]: real zeros of $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$
/// - [`bessel_y_prime`](crate::bessel_y_prime): derivatives of $Y_v(z)$
#[doc(alias = "y1p_zeros", alias = "cyzo")]
#[must_use]
#[inline]
pub fn bessel_y1_prime_zeros<const N: usize>(complex: bool) -> [[Complex<f64>; N]; 2] {
    cyzo(2, complex)
}

/// Compute the first $N$ zeros of all $J_n(x)$ and $J_n\'(x)$ with integer $n \ge 0$, in ascending
/// order
///
/// Corresponds to [`scipy.special.jnjnp_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnjnp_zeros.html
///
/// # Returns
///
/// - The zeros $x_l$, starting with $x_1 = 0$ (of $J_0\'(x)$)
/// - The orders $n_l$
/// - The serial numbers $m_l$ of the zeros $x_l$ of $J_{n_l}(x)$ or $J_{n_l}\'(x)$ (from 0 for
///   $J_0\'(x)$, so that $x_1 = 0$ has $m_1 = 0$, and from 1 otherwise)
/// - Whether $x_l$ is a zero of $J_{n_l}\'(x)$ (`true`) or of $J_{n_l}(x)$ (`false`)
///
/// # Panics
/// - if `N` is greater than 1200
///
/// # See also
/// - [`bessel_zeros`]: zeros of $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$ for a given $v$
/// - [`bessel_j`]: Bessel function $J_v(z)$
#[doc(alias = "jnjnp_zeros", alias = "jdzo")]
#[must_use]
pub fn bessel_j_jprime_zeros<const N: usize>() -> ([f64; N], [u32; N], [u32; N], [bool; N]) {
    assert!(N <= 1200, "N must be at most 1200");

    // `jdzo` computes at most `N + 24` zeros (for `N <= 1200`), into the 1-based `zo`
    let len = N + 25;
    let (mut zo, mut n, mut m, mut p) =
        (vec![0.0; len], vec![0; len], vec![0; len], vec![false; len]);
    specfun::jdzo(N as i32, &mut zo, &mut n, &mut m, &mut p);
    (
        core::array::from_fn(|l| zo[l + 1]),
        core::array::from_fn(|l| n[l]),
        core::array::from_fn(|l| m[l]),
        core::array::from_fn(|l| p[l]),
    )
}

// Hankel 1

/// Hankel function of the first kind, $H_v^{(1)}(z)$
///
/// Corresponds to [`scipy.special.hankel1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1.html
///
/// # See also
/// - [`hankel_1e`]: $H_v^{(1)}$ without the leading exponential factor
/// - [`hankel_1_prime`](crate::hankel_1_prime): $n$th derivative, ${d^n\over dz^n}H_v^{(1)}(z)$
/// - [`hankel_2`]: Hankel function of the second kind $H_v^{(2)}(z)$
#[doc(alias = "h1v", alias = "hankel1", alias = "cyl_hankel_1")]
#[must_use]
#[inline]
pub fn hankel_1<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    let z: Complex<f64> = z.into();
    z.via(|z| cyl_hankel_c64(v, z, 1, 1))
}

/// Exponentially scaled Hankel function of the first kind, $e^{-iz} H_v^{(1)}(z)$
///
/// Corresponds to [`scipy.special.hankel1e`][h1ve].
///
/// [h1ve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1e.html
///
/// # See also
/// - [`hankel_1`]: Hankel function $H_v^{(1)}(z)$
#[doc(alias = "h1ve", alias = "hankel1e", alias = "cyl_hankel_1e")]
#[must_use]
#[inline]
pub fn hankel_1e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    let z: Complex<f64> = z.into();
    z.via(|z| cyl_hankel_scaled_c64(v, z, 1))
}

// Hankel 2

/// Hankel function of the second kind, $H_v^{(2)}(z)$
///
/// Corresponds to [`scipy.special.hankel2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2.html
///
/// # See also
/// - [`hankel_2e`]: $H_v^{(2)}$ without the leading exponential factor
/// - [`hankel_2_prime`](crate::hankel_2_prime): $n$th derivative, ${d^n\over dz^n}H_v^{(2)}(z)$
/// - [`hankel_1`]: Hankel function of the first kind $H_v^{(1)}(z)$
#[doc(alias = "h2v", alias = "hankel2", alias = "cyl_hankel_2")]
#[must_use]
#[inline]
pub fn hankel_2<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    let z: Complex<f64> = z.into();
    z.via(|z| cyl_hankel_c64(v, z, 1, 2))
}

/// Exponentially scaled Hankel function of the second kind, $e^{iz} H_v^{(2)}(z)$
///
/// Corresponds to [`scipy.special.hankel2e`][h2ve].
///
/// [h2ve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2e.html
///
/// # See also
/// - [`hankel_2`]: Hankel function $H_v^{(2)}(z)$
#[doc(alias = "h2ve", alias = "hankel2e", alias = "cyl_hankel_2e")]
#[must_use]
#[inline]
pub fn hankel_2e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    let z: Complex<f64> = z.into();
    z.via(|z| cyl_hankel_scaled_c64(v, z, 2))
}

// Integrals

/// `int_0^1 t^lambda J_nu(2a t) dt` as the Neumann series DLMF 10.22.10 with `x = 2a`
fn besselpoly_neumann(a: f64, lambda: f64, nu: f64) -> f64 {
    // Unlike in xsf, use AMOS where cephes `jv` recurs down from the power series at order
    // ~ v + x, which cancels for large x and can underflow: for a non-integer order 2x < v < 500
    // with x >= 3.6 sqrt(v).
    let jv = |v: f64, x: f64| {
        if v.fract() != 0.0 && v > 2.0 * x && v < 500.0 && x >= 3.6 * v.sqrt() {
            cyl_bessel_j_f64(v, x)
        } else {
            crate::xsf::cephes::jv(v, x)
        }
    };

    let b = 2.0 * a;
    let c0 = (lambda + nu + 1.0) / 2.0;
    let c1 = c0 - lambda;
    let (mut f, mut g, mut sum) = (1.0, 0.0, 0.0); // J_mu(b), J_(mu+1)(b) up to a common factor
    let k_max = (f64::from(specfun::msta2(b, (nu + 1.0).ceil() as i32, 15)) - nu) / 2.0;
    for k in (0..=k_max as i32).rev() {
        let k = f64::from(k);
        let mu = nu + 1.0 + 2.0 * k;
        for j in [2.0, 1.0] {
            // J_(mu+j-1) = 2(mu+j)/b J_(mu+j) - J_(mu+j+1)
            (f, g) = (2.0 * (mu + j) / b * f - g, f);
        }
        sum = mu * f + (c1 + k) / (c0 + k + 1.0) * sum;
    }
    sum * (jv(nu + 1.0, b) * f + jv(nu + 2.0, b) * g) / (f * f + g * g) / b / c0
}

/// Weighted integral of the Bessel function of the first kind,
/// $\int_0^1 x^\lambda \mathop{J}_v(2ax) \dd x$
///
/// Corresponds to [`scipy.special.besselpoly`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.besselpoly.html
#[must_use]
#[inline]
pub fn besselpoly(a: f64, lambda: f64, v: f64) -> f64 {
    // Unlike in xsf, J_v(-y) = (-1)^v J_v(y) for a < 0 and integer v, instead of NaN
    if a < 0.0 {
        let s = besselpoly(-a, lambda, v);
        return if v % 2.0 == 0.0 {
            s
        } else if v.fract() == 0.0 {
            -s
        } else {
            f64::NAN
        };
    }

    let n = if v < 0.0 && v.floor() == v { -v } else { v }; // J_(-n) = (-1)^n J_n
    let negate = n != v && n % 2.0 == 1.0;

    // Unlike in xsf, +-inf at the poles of the series (lambda + n + 1 = 0, -2, ...), with the sign
    // of the leading term, instead of a partial sum (only where `pole` is exact, so above -2^53)
    let (s, e0) = two_sum(lambda, n);
    let (pole, e1) = two_sum(s, 1.0);
    if a > 0.0 && e0 == 0.0 && e1 == 0.0 && pole <= 0.0 && pole % 2.0 == 0.0 {
        let s = f64::INFINITY * crate::xsf::cephes::gammasgn(n + 1.0);
        return if negate { -s } else { s };
    }

    // Use the Neumann series where the power series cancels badly, and where the Neumann series
    // is accurate:
    // - a >= 1.5 and n < a^2 - 1: otherwise the power series terms grow to at most ~e^3, or
    //   decrease from the start.
    // - n > -1: all orders n + 2k + 1 are positive, where Miller's backward recurrence gives J.
    // - lambda >= -1: then c0 > -1/2 and |c1 + k| <= c0 + k + 1, so the coefficients don't grow.
    // - a < 500: the recurrence isn't rescaled, and f * f overflows from about a = 800.
    // For lambda + n + 1 <= 0 the integral diverges, but the Neumann series (like the power
    // series) still converges, since J_(n+2k+1)(2a) decays rapidly in k, to the same analytic
    // continuation in lambda.
    if (1.5..500.0).contains(&a)
        && n > -1.0
        && n < a * a - 1.0
        && lambda >= -1.0
        && lambda.is_finite()
    {
        let s = besselpoly_neumann(a, lambda, n);
        return if negate { -s } else { s };
    }
    crate::xsf::cephes::besselpoly(a, lambda, v)
}

/* Integrals of bessel functions */

fn itika(x: f64) -> (f64, f64) {
    // =======================================================
    // Purpose: Integrate modified Bessel functions I0(t) and
    //          K0(t) with respect to t from 0 to x
    // Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
    // Output:  TI --- Integration of I0(t) from 0 to x
    //          TK --- Integration of K0(t) from 0 to x
    // =======================================================

    const A: [f64; 28] = [
        0.625,
        1.0078125,
        2.5927734375,
        9.1868591308594,
        4.1567974090576e+1,
        2.2919635891914e+2,
        1.491504060477e+3,
        1.1192354495579e+4,
        9.515939374212e+4,
        9.0412425769041e+5,
        // Unlike in xsf (scipy/xsf#379), also A(11:28), for the asymptotic expansion of `ti`:
        // A(k) = a(k) + (k - 1/2) A(k-1), with a(k) from DLMF 10.40.1
        9.493856041645449e6,
        1.0918238256943358e8,
        1.3647980398733943e9,
        1.842489237671708e10,
        2.6716177232170163e11,
        4.1410137239378687e12,
        6.832677651456437e13,
        1.195719014944093e15,
        2.21208056127209e16,
        4.313557459164584e17,
        8.842793156188478e18,
        1.9012005669158696e20,
        4.277701317750422e21,
        1.0052598145253639e23,
        2.462886551414386e24,
        6.280360713393542e25,
        1.6642955899969174e27,
        4.576812873771245e28,
    ];
    let el = 0.5772156649015329;

    let mut ti;
    let mut tk = 0.0;
    let mut r;
    let x2;

    if x == 0.0 {
        return (0.0, 0.0);
    } else if x < 40.0 {
        // Unlike in xsf (scipy/xsf#379), the series of `ti` is used for x < 40 (instead of 20)
        x2 = x * x;
        ti = 1.0;
        r = 1.0;
        for k in 1..=50 {
            r = 0.25 * r * (f64::from(2 * k) - 1.0) / (f64::from(2 * k) + 1.0) / f64::from(k * k)
                * x2;
            ti += r;
            // Unlike in xsf (scipy/xsf#379), stop at machine epsilon instead of 1e-12, which loses
            // up to 2e-13
            if (r / ti).abs() < f64::EPSILON {
                break;
            }
        }
        ti *= x;
    } else {
        x2 = 0.0;
        ti = 1.0;
        r = 1.0;
        for k in 1..=28 {
            r /= x;
            ti += A[k - 1] * r;
        }
        let rc1 = 1.0 / (TAU * x).sqrt();
        ti *= rc1 * x.exp();
    }

    if x < 0.5 {
        let e0 = el + (x / 2.0).ln();
        let mut b1 = 1.0 - e0;
        let mut b2 = 0.0;
        let mut rs = 0.0;
        r = 1.0;
        let mut tw = 0.0;
        for k in 1..=50 {
            r = 0.25 * r * (f64::from(2 * k) - 1.0) / (f64::from(2 * k) + 1.0) / f64::from(k * k)
                * x2;
            b1 += r * (1.0 / f64::from(2 * k + 1) - e0);
            rs += 1.0 / f64::from(k);
            b2 += r * rs;
            tk = b1 + b2;
            if ((tk - tw) / tk).abs() < 1.0e-12 {
                break;
            }
            tw = tk;
        }
        tk *= x;
    } else if x < 40.0 {
        // Unlike in xsf (scipy/xsf#379), pi/2 - int_x^inf K0 for 0.5 <= x < 40 (instead of the
        // series for x < 12), with int_x^inf K0 = int_0^inf exp(-x cosh(u)) / cosh(u) du
        // (DLMF 10.32.9)
        tk = FRAC_PI_2 - trapezoid_cosh(x, |w| (-w).exp() / (x + w)) * x * (-x).exp();
    } else {
        tk = 1.0;
        r = 1.0;
        for k in 1..=10 {
            r = -r / x;
            tk += A[k - 1] * r;
        }
        let rc2 = (PI / (2.0 * x)).sqrt();
        tk = FRAC_PI_2 - rc2 * tk * (-x).exp();
    }
    (ti, tk)
}

/// `h (f(0) / 2 + f(w_1) + f(w_2) + ...)` with `w_j = x (cosh(j h) - 1)`, i.e. the trapezoidal
/// rule for `int_0^inf f(x (cosh(u) - 1)) du`, for `f(w)` that decays like `exp(-w)` and has its
/// nearest singularity at `u = i pi/2` (`w = -x`).
///
/// Toward it, `|f|` grows like `exp(x)`, so the step `h = pi^2 / (x + 45)` gives a relative
/// discretization error of order `exp(x - pi^2 / h) = exp(-45)`.
fn trapezoid_cosh(x: f64, f: impl Fn(f64) -> f64) -> f64 {
    let h = PI * PI / (x + 45.0);
    let mut s = 0.5 * f(0.0);
    for j in 1..=40 {
        let v = (0.5 * h * f64::from(j)).sinh();
        let t = f(2.0 * x * v * v);
        s += t;
        if t < 1e-17 * s {
            break;
        }
    }
    h * s
}

/// Calls `f(n, J_n(x))` for `n = N, ..., 1` with the unnormalized `J_n` of Miller's backward
/// recurrence from `J_{N+1} = 0`, `J_N = 1`, and returns `J_0 + 2 (J_2 + J_4 + ...)`, i.e. 1 for
/// the exact `J_n` (DLMF 10.12.4). For `0.5 <= x < 40`, the odd `N = 2 floor(x/2) + 51` has
/// `|J_{N+1}(x) / Y_{N+1}(x)| < 1e-43`, so the error from `J_{N+1} = 0` is far below eps.
pub(crate) fn bessel_j_miller(x: f64, mut f: impl FnMut(u32, f64)) -> f64 {
    let (mut jp, mut j, mut je) = (0.0, 1.0, 0.0);
    for k in (1..=2 * (0.5 * x) as u32 + 51).rev() {
        f(k, j);
        if k % 2 == 0 {
            je += j;
        }
        (jp, j) = (j, 2.0 * f64::from(k) / x * j - jp);
    }
    j + 2.0 * je
}

/// `[int_0^x J0, int_0^x (1 - J0(t))/t dt, int_0^x Y0, int_x^inf Y0(t)/t dt]` for
/// `0.5 <= x < 40`, by integrating DLMF 10.23.16, `Y0 = 2/pi (L J0 - 2 sum_k (-1)^k J_2k / k)`
/// with `L = ln(x/2) + euler`, term by term:
///
/// - `int_0^x J_n = 2 (J_{n+1} + J_{n+3} + ...)` (DLMF 10.22.8)
/// - `int_0^x f_n(t)/t dt = (2 (f_n + f_{n+2} + ...) - f_n) / n`, for `f_n = J_n` (`w_n`), and
///   also for `f_n = w_n`
/// - `int_0^x (1 - J0(t))/t dt = 2 sum_k w_2k`, since `1 - J0 = 2 (J_2 + J_4 + ...)`
/// - `L J0` by parts, and for `int_x^inf` as `L/t - L (1 - J0(t))/t`, where `L/t` gives
///   `-L^2/pi`, plus the integration constant `pi/6`
fn itj0y0_neumann(x: f64) -> [f64; 4] {
    // the sums of J_n and w_n over the odd n, of w_n and of int_0^x w_n(t)/t dt over the even
    // n, and two alternating sums (`a = t_k - a` for k downwards gives `a = -sum_k (-1)^k t_k`)
    let (mut jo, mut wo, mut je, mut we, mut ww) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut ajo, mut awe) = (0.0, 0.0);
    let norm = bessel_j_miller(x, |k, j| {
        let n = f64::from(k);
        if k % 2 == 1 {
            jo += j;
            wo += (2.0 * jo - j) / n;
        } else {
            je += j;
            let w = (2.0 * je - j) / n;
            we += w;
            ww += (2.0 * we - w) / n;
            ajo = jo / n - ajo;
            awe = w / n - awe;
        }
    });
    let l = (0.5 * x).ln() + EGAMMA;
    let (tj, ttj) = (2.0 * jo / norm, 2.0 * we / norm);
    let ty = FRAC_2_PI * (l * tj - 2.0 * (wo - 4.0 * ajo) / norm);
    let tty = PI / 6.0 + FRAC_2_PI * (l * (ttj - 0.5 * l) - 2.0 * (ww + 2.0 * awe) / norm);
    [tj, ttj, ty, tty]
}

/// The coefficients `a_k` of the asymptotic expansions of the integrals of `J0`, `Y0`, `H0` and
/// `L0`, which `ITJYA`, `ITSH0` and `ITSL0` of specfun each compute with this recurrence
pub(crate) const IT0_ASYMP: [f64; 25] = {
    let mut a = [5.0 / 8.0; 25];
    let (mut a0, mut k) = (1.0, 1);
    while k < 25 {
        let (kf, a1) = (k as f64, a[k - 1]);
        a[k] = (1.5 * (kf + 0.5) * (kf + 5.0 / 6.0) * a1
            - 0.5 * (kf + 0.5) * (kf + 0.5) * (kf - 0.5) * a0)
            / (kf + 1.0);
        a0 = a1;
        k += 1;
    }
    a
};

/// `(int_0^x J0(t) dt, int_0^x Y0(t) dt)` by `n + 1` terms of each of the two sums of their
/// asymptotic expansion, as in `ITJYA` (also used by `ITSH0` of specfun)
pub(crate) fn itj0y0_asymp(x: f64, n: usize) -> (f64, f64) {
    let mut bf = 1.0;
    let mut r = 1.0;
    for k in 1..=n {
        r = -r / (x * x);
        bf += IT0_ASYMP[2 * k - 1] * r;
    }
    let mut bg = IT0_ASYMP[0] / x;
    r = 1.0 / x;
    for k in 1..=n {
        r = -r / (x * x);
        bg += IT0_ASYMP[2 * k] * r;
    }
    let xp = x + FRAC_PI_4;
    let rc = (2.0 / (PI * x)).sqrt();
    let tj = 1.0 - rc * (bf * xp.cos() + bg * xp.sin());
    (tj, rc * (bg * xp.cos() - bf * xp.sin()))
}

fn itjya(x: f64) -> (f64, f64) {
    // Unlike in xsf (scipy/xsf#446), use the Neumann series for 0.5 <= x < 40. That leaves the
    // power series below to x < 0.5, and the asymptotic expansion to x >= 40.
    if (0.5..40.0).contains(&x) {
        let [tj, _, ty, _] = itj0y0_neumann(x);
        return (tj, ty);
    }

    let el = 0.5772156649015329;
    let eps = 1.0e-12;

    let tj;
    let ty;
    if x == 0.0 {
        tj = 0.0;
        ty = 0.0;
    } else if x <= 20.0 {
        let x2 = x * x;
        let mut tj_ = x;
        let mut r = x;
        for ki in 1..61 {
            let k = f64::from(ki);
            r = -0.25 * r * (2.0 * k - 1.0) / (2.0 * k + 1.0) / f64::from(ki * ki) * x2;
            tj_ += r;
            if r.abs() < tj_.abs() * eps {
                break;
            }
        }
        tj = tj_;
        let ty1 = (el + (x / 2.0).ln()) * tj;
        let mut rs = 0.0;
        let mut ty2 = 1.0;
        r = 1.0;
        for ki in 1..61 {
            let k = f64::from(ki);
            r = -0.25 * r * (2.0 * k - 1.0) / (2.0 * k + 1.0) / f64::from(ki * ki) * x2;
            rs += 1.0 / k;
            let r2 = r * (rs + 1.0 / (2.0 * k + 1.0));
            ty2 += r2;
            if r2.abs() < ty2.abs() * eps {
                break;
            }
        }
        ty = (ty1 - x * ty2) * 2.0 / PI;
    } else {
        (tj, ty) = itj0y0_asymp(x, 8);
    }
    (tj, ty)
}

fn ittika(x: f64) -> (f64, f64) {
    // =========================================================
    // Purpose: Integrate [I0(t)-1]/t with respect to t from 0
    //          to x, and K0(t)/t with respect to t from x to ∞
    // Input :  x   --- Variable in the limits  ( x ≥ 0 )
    // Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
    //          TTK --- Integration of K0(t)/t from x to ∞
    // =========================================================

    const C: [f64; 28] = [
        1.625,
        4.1328125,
        1.45380859375e1,
        6.553353881835e1,
        3.6066157150269e2,
        2.3448727161884e3,
        1.7588273098916e4,
        1.4950639538279e5,
        // Unlike in xsf (scipy/xsf#379), also C(9:28), for both asymptotic expansions:
        // C(k) = a(k) + (k + 1/2) C(k-1), with a(k) from DLMF 10.40.1
        1.4203351366661638e6,
        1.491362895213499e7,
        1.715072842854485e8,
        2.143844091658617e9,
        2.8941913495146805e10,
        4.1965786451805493e11,
        6.504697732889155e12,
        1.0732751884562256e14,
        1.8782316298679842e15,
        3.4747285578496924e16,
        6.775720726169453e17,
        1.3890227525115777e19,
        2.986398921548903e20,
        6.719397577318566e21,
        1.5790584310917602e23,
        3.868693156660214e24,
        9.865167550066271e25,
        2.6142694008404304e27,
        7.189240852405946e28,
        2.048933642948492e30,
    ];

    let el = 0.5772156649015329;
    let mut tti;
    let mut r;

    if x == 0.0 {
        return (0.0, SPECFUN_INF);
    }
    if x < 40.0 {
        tti = 1.0;
        r = 1.0;
        for k in 2..=50 {
            r = 0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            tti += r;
            // Unlike in xsf (scipy/xsf#379), stop at machine epsilon instead of 1e-12, which loses
            // up to 2e-13
            if (r / tti).abs() < f64::EPSILON {
                break;
            }
        }
        tti *= 0.125 * x * x;
    } else {
        tti = 1.0;
        r = 1.0;
        for k in 1..=28 {
            r /= x;
            tti += C[k - 1] * r;
        }
        let rc = x * (TAU * x).sqrt();
        tti = tti * x.exp() / rc;
    }
    let ttk = if x <= 1.0 {
        let ln_half_x = (x / 2.0).ln();
        let e0 = (0.5 * ln_half_x + el) * ln_half_x + PI * PI / 24.0 + 0.5 * el * el;
        let mut b1 = 1.5 - (el + ln_half_x);
        let mut rs = 1.0;
        r = 1.0;
        for k in 2..=50 {
            r = 0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            rs += 1.0 / f64::from(k);
            let r2 = r * (rs + 1.0 / (2.0 * f64::from(k)) - (el + ln_half_x));
            b1 += r2;
            if (r2 / b1).abs() < 1.0e-12 {
                break;
            }
        }
        e0 - 0.125 * x * x * b1
    } else if x < 40.0 {
        // Unlike in xsf (scipy/xsf#379), the trapezoidal rule on int_0^inf E1(x cosh(u)) du
        // (DLMF 10.32.9)
        trapezoid_cosh(x, |w| (-w).exp() * scaled_exp1(x + w) / (x + w)) * (-x).exp()
    } else {
        let mut ttk_ = 1.0;
        r = 1.0;
        for k in 1..=28 {
            r = -r / x;
            ttk_ += C[k - 1] * r;
        }
        let rc = x * (FRAC_2_PI * x).sqrt();
        ttk_ * (-x).exp() / rc
    };
    (tti, ttk)
}

fn ittjya(x: f64) -> (f64, f64) {
    // =========================================================
    // Purpose: Integrate [1-J0(t)]/t with respect to t from 0
    //          to x, and Y0(t)/t with respect to t from x to ∞
    // Input :  x   --- Variable in the limits  ( x ≥ 0 )
    // Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
    //          TTY --- Integration of Y0(t)/t from x to ∞
    // =========================================================

    // Unlike in xsf (scipy/xsf#446), use the Neumann series for 0.5 <= x < 34 (from about 34, the
    // asymptotic expansion below is more accurate). That leaves the power series below to x < 0.5.
    if (0.5..34.0).contains(&x) {
        let [_, ttj, _, tty] = itj0y0_neumann(x);
        return (ttj, tty);
    }

    let el = 0.5772156649015329;
    let mut ttj;
    let tty;
    let mut r;

    if x == 0.0 {
        ttj = 0.0;
        tty = -SPECFUN_INF;
    } else if x <= 20.0 {
        ttj = 1.0;
        r = 1.0;
        for k in 2..=100 {
            r = -0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            ttj += r;
            if r.abs() < ttj.abs() * 1.0e-12 {
                break;
            }
        }
        ttj *= 0.125 * x * x;
        let ln_half_x = (x / 2.0).ln();
        let e0 = 0.5 * (FRAC_PI_SQ_6 - el * el) - (0.5 * ln_half_x + el) * ln_half_x;
        let mut b1 = el + ln_half_x - 1.5;
        let mut rs = 1.0;
        r = -1.0;
        for k in 2..=100 {
            r = -0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            rs += 1.0 / f64::from(k);
            let r2 = r * (rs + 1.0 / (2.0 * f64::from(k)) - (el + ln_half_x));
            b1 += r2;
            if r2.abs() < b1.abs() * 1.0e-12 {
                break;
            }
        }
        tty = FRAC_2_PI * (e0 + 0.125 * x * x * b1);
    } else {
        let a0 = (2.0 / (PI * x)).sqrt();
        let mut bj0 = 0.0;
        let mut by0 = 0.0;
        let mut bj1 = 0.0;
        let mut by1 = 0.0;
        for l in 0..=1 {
            let l = f64::from(l);
            let vt = 4.0 * l * l;
            let mut px = 1.0;
            r = 1.0;
            // `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM
            // only does when optimizing
            for k in 1..=14 {
                let k = f64::from(k);
                let (p0, p1) = (4.0 * k - 3.0, 4.0 * k - 1.0);
                r = -0.0078125 * r * (vt - p0 * p0) / (x * k) * (vt - p1 * p1)
                    / ((2.0 * k - 1.0) * x);
                px += r;
                if r.abs() < px.abs() * 1.0e-12 {
                    break;
                }
            }
            let mut qx = 1.0;
            r = 1.0;
            for k in 1..=14 {
                let k = f64::from(k);
                let (q0, q1) = (4.0 * k - 1.0, 4.0 * k + 1.0);
                r = -0.0078125 * r * (vt - q0 * q0) / (x * k) * (vt - q1 * q1)
                    / ((2.0 * k + 1.0) * x);
                qx += r;
                if r.abs() < qx.abs() * 1.0e-12 {
                    break;
                }
            }
            qx *= 0.125 * (vt - 1.0) / x;
            let xk = x - (0.25 + 0.5 * l) * PI;
            bj1 = a0 * (px * xk.cos() - qx * xk.sin());
            by1 = a0 * (px * xk.sin() + qx * xk.cos());
            if l == 0.0 {
                bj0 = bj1;
                by0 = by1;
            }
        }
        let t = 2.0 / x;
        let mut g0 = 1.0;
        let mut r0 = 1.0;
        for k in 1..=10 {
            let k = f64::from(k);
            r0 = -r0 * k * k * t * t;
            g0 += r0;
        }
        let mut g1 = 1.0;
        let mut r1 = 1.0;
        for k in 1..=10 {
            let k = f64::from(k);
            r1 = -r1 * k * (k + 1.0) * t * t;
            g1 += r1;
        }
        ttj = 2.0 * g1 * bj0 / (x * x) - g0 * bj1 / x + el + (x / 2.0).ln();
        tty = 2.0 * g1 * by0 / (x * x) - g0 * by1 / x;
    }
    (ttj, tty)
}

/// Integrals of Bessel functions of order 0
///
/// Corresponds to [`scipy.special.itj0y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itj0y0.html
///
/// # Returns
/// - $\int_0^x J_0(t) \dd t$
/// - $\int_0^x Y_0(t) \dd t$
///
/// # See also
/// - [`it2j0y0`]
/// - [`bessel_j0`]
/// - [`bessel_y0`]
#[doc(alias = "itj0y0")]
#[must_use]
#[inline]
pub fn it1j0y0(mut x: f64) -> (f64, f64) {
    /* int(j0(t),t=0..x) */
    /* int(y0(t),t=0..x) */
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (mut j0int, mut y0int) = itjya(x);
    if flag {
        j0int = -j0int;
        y0int = f64::NAN; /* domain error */
    }
    (j0int, y0int)
}

/// Integrals related to Bessel functions of order 0
///
/// Corresponds to [`scipy.special.it2j0y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2j0y0.html
///
/// # Returns
/// - $\int_0^x {1 - J_0(t) \over t} \dd t$
/// - $\int_x^\infty {Y_0(t) \over t} \dd t$
///
/// # See also
/// - [`it1j0y0`]
/// - [`bessel_j0`]
/// - [`bessel_y0`]
#[must_use]
#[inline]
pub fn it2j0y0(mut x: f64) -> (f64, f64) {
    /* int((1-j0(t))/t,t=0..x) */
    /* int(y0(t)/t,t=x..inf) */
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (j0int, mut y0int) = ittjya(x);
    if flag {
        y0int = f64::NAN; /* domain error */
    }
    (j0int, y0int)
}

/// Integrals of modified Bessel functions of order 0
///
/// Corresponds to [`scipy.special.iti0k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iti0k0.html
///
/// # Returns
/// - $\int_0^x I_0(t) \dd t$
/// - $\int_0^x K_0(t) \dd t$
///
/// # See also
/// - [`it2i0k0`]
/// - [`bessel_i0`]
/// - [`bessel_k0`]
#[doc(alias = "iti0k0")]
#[must_use]
#[inline]
pub fn it1i0k0(mut x: f64) -> (f64, f64) {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (mut i0int, mut k0int) = itika(x);
    if flag {
        i0int = -i0int;
        k0int = f64::NAN; /* domain error */
    }
    (i0int, k0int)
}

/// Integrals related to modified Bessel functions of order 0.
///
/// Corresponds to [`scipy.special.it2i0k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2i0k0.html
///
/// # Returns
/// - $\int_0^x {I_0(t) - 1 \over t} \dd t$
/// - $\int_x^\infty {K_0(t) \over t} \dd t$
///
/// # See also
/// - [`it1i0k0`]
/// - [`bessel_i0`]
/// - [`bessel_k0`]
#[must_use]
#[inline]
pub fn it2i0k0(mut x: f64) -> (f64, f64) {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (i0int, mut k0int) = ittika(x);
    if flag {
        k0int = f64::NAN; /* domain error */
    }
    (i0int, k0int)
}

// Riccati-Bessel

fn rctj(x: f64, nm: &mut i32, rj: &mut [f64], dj: &mut [f64]) {
    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the first
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of jn(x)  ( n = 0,1,2,... )
    // Output:  RJ(n) --- x·jn(x)
    //          DJ(n) --- [x·jn(x)]'
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    // ========================================================

    // NOTE: rj[0] and rj[1] are always written (for |x| >= 1e-100), so buffers of length 1
    // (n == 0) are out of bounds; riccati_j uses length 2 for N == 1.
    //
    // NOTE: for x NaN, +-inf, or |x| >= 2^31 / 1.1, msta2 hits its `(int)` UB (see msta1), and its
    // `+ 5` / `+ 10` can overflow (panic in debug builds) for x +-inf and |x| >= ~2^31.
    let n = rj.len() as i32 - 1;

    *nm = n;
    if x.abs() < 1.0e-100 {
        // Unlike in xsf, these are the leading terms `x^(k+1) / (2k+1)!!` (DLMF 10.53.1) instead
        // of 0 (except `dj[0] = 1`)
        let mut u = 1.0; // x^k / (2k+1)!!
        for k in 0..=n {
            rj[k as usize] = u * x;
            dj[k as usize] = f64::from(k + 1) * u;
            u *= x / f64::from(2 * k + 3);
        }
        return;
    }
    rj[0] = x.sin();
    rj[1] = rj[0] / x - x.cos();
    let rj0 = rj[0];
    let rj1 = rj[1];
    let mut f = 0.0;

    // Unlike in xsf, this also recurs for `n = 1` if `|x| < 1`, where `rj[1]` cancels. It always
    // starts at `msta2`, instead of leaving the orders above `msta1(x, 200)` unset. The
    // derivatives are of the unnormalized `f`, since `rj[k]` can underflow where `dj[k]` doesn't.
    if n >= 2 || x.abs() < 1.0 {
        let m = specfun::msta2(x, n, 15);

        let mut f0 = 0.0;
        let mut f1 = 1.0e-100;

        for k in (0..=m).rev() {
            f = (2.0 * f64::from(k) + 3.0) * f1 / x - f0;
            if k <= *nm {
                rj[k as usize] = f;
            }
            if k < *nm {
                dj[(k + 1) as usize] = f64::from(-(k + 1)) * f1 / x + f;
            }
            if f.abs() > SPECFUN_BIG {
                (f, f1) = (f / SPECFUN_BIG, f1 / SPECFUN_BIG);
                for (r, d) in rj.iter_mut().zip(&mut *dj).skip(k as usize) {
                    (*r, *d) = (*r / SPECFUN_BIG, *d / SPECFUN_BIG);
                }
            }
            f0 = f1;
            f1 = f;
        }
        let cs = if rj0.abs() > rj1.abs() {
            rj0 / f
        } else {
            rj1 / f0
        };
        for k in 0..=*nm {
            rj[k as usize] *= cs;
            dj[k as usize] *= cs;
        }
    } else {
        dj[1] = -rj1 / x + rj0;
    }
    dj[0] = x.cos();
}

// Unlike in xsf, reorder and halve overflowing products in a * b / x + c instead of using FMA.
// Return the product without c if it still overflows.
#[inline]
fn mul_div_add(a: f64, b: f64, x: f64, c: f64) -> f64 {
    #[cold]
    fn reordered(a: f64, b: f64, x: f64, c: f64) -> f64 {
        // halved, since `a * b / x` can be up to twice `f64::MAX` if `c` cancels it
        let (a, c) = (0.5 * a, 0.5 * c);
        let t = a * (b / x);
        let t = if t.is_finite() { t } else { a / x * b };
        if t.is_infinite() { t } else { 2.0 * (t + c) }
    }
    let t = a * b / x;
    if t.is_finite() {
        t + c
    } else {
        reordered(a, b, x, c)
    }
}

fn rcty(x: f64, nm: &mut i32, ry: &mut [f64], dy: &mut [f64]) {
    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the second
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of yn(x)
    // Output:  RY(n) --- x·yn(x)
    //          DY(n) --- [x·yn(x)]'
    //          NM --- Highest order computed
    // ========================================================

    let n = ry.len() as i32 - 1;

    *nm = n;
    ry[0] = -x.cos();
    if n <= 0 {
        dy[0] = x.sin();
        return;
    }
    ry[1] = ry[0] / x - x.sin();
    let mut rf0 = ry[0];
    let mut rf1 = ry[1];

    // Unlike in xsf, continue the same recurrence after overflow.
    for k in 2..=n {
        let rf2 = mul_div_add(2.0 * f64::from(k) - 1.0, rf1, x, -rf0);
        ry[k as usize] = rf2;
        rf0 = rf1;
        rf1 = rf2;
    }

    dy[0] = x.sin();
    for k in 1..=n {
        dy[k as usize] = mul_div_add(f64::from(-k), ry[k as usize], x, ry[(k - 1) as usize]);
    }
}

/// Compute Riccati-Bessel function of the first kind and derivatives for the first $N$ orders
///
/// The Riccati-Bessel function of the first kind is defined as $x \mathop{j}_n(x)$, where $j_n$ is
/// the spherical Bessel function of the first kind of order $n$.
///
/// This function computes the value and first derivative of the function for the orders `0..N`.
///
/// Corresponds to [`scipy.special.riccati_jn`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_jn.html
///
/// # Returns
/// - $\[x j_0(x), \ldots, x j_{N-1}(x)\]$
/// - Their derivatives with respect to $x$
///
/// # See also
/// - [`riccati_y`]
/// - [`bessel_j`]
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
#[doc(alias = "rctj", alias = "riccati_jn")]
#[must_use]
#[inline]
pub fn riccati_j<const N: usize>(x: f64) -> ([f64; N], [f64; N]) {
    let (mut rj, mut dj) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 {
        return (rj, dj);
    }
    if N == 1 {
        // like in scipy, since `rctj` writes at least 2 terms
        let (rj2, dj2) = riccati_j::<2>(x);
        rj.copy_from_slice(&rj2[..N]);
        dj.copy_from_slice(&dj2[..N]);
        return (rj, dj);
    }

    if N > 2 && (x.is_nan() || x.abs() > 1e9) {
        // Unlike in SciPy, this avoids the `int` casts of `|x|` in `msta2`, which are undefined for
        // NaN and `|x| >= 2**31 / 1.1`.
        return (rj, dj);
    }

    let mut nm = 0;
    rctj(x, &mut nm, &mut rj, &mut dj);
    (rj, dj)
}

/// Compute Riccati-Bessel function of the second kind and derivatives for the first $N$ orders
///
/// The Riccati-Bessel function of the second kind is defined here as $x \mathop{y}_n(x)$, where
/// $y_n$ is the spherical Bessel function of the second kind of order $n$.
/// *Note that this is in contrast to a common convention that includes a minus sign in the
/// definition.*
///
/// This function computes the value and first derivative of the function for the orders `0..N`.
///
/// Corresponds to [`scipy.special.riccati_yn`][scipy].
///
/// # Returns
/// - $\[x y_0(x), \ldots, x y_{N-1}(x)\]$
/// - Their derivatives with respect to $x$
///
/// # See also
/// - [`riccati_j`]
/// - [`bessel_y`]
/// - [`sph_bessel_y`](crate::sph_bessel_y)
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_yn.html
#[doc(alias = "rcty", alias = "riccati_yn")]
#[must_use]
#[inline]
pub fn riccati_y<const N: usize>(x: f64) -> ([f64; N], [f64; N]) {
    let (mut ry, mut dy) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 {
        return (ry, dy);
    }
    let mut nm = 0;
    rcty(x, &mut nm, &mut ry, &mut dy);
    (ry, dy)
}

/// Evaluate the Jahnke-Emden Lambda function $\Lambda_v(x)$ and its derivatives
///
/// The Jahnke-Emden Lambda function is defined as
///
/// $$\Lambda_v(x) = \Gamma(v+1) \left({2 \over x}\right)^v J_v(x)$$
///
/// where $\Gamma$ is the [Gamma function](crate::gamma) and $J_v$ is the
/// [Bessel function of the first kind](crate::bessel_j).
///
/// Corresponds to [`scipy.special.lmbda`][scipy] in SciPy.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lmbda.html
///
/// # Arguments
/// - `v`: Order of the Lambda function as a non-negative real number
/// - `x`: Value at which to evaluate the function and derivatives
///
/// # Returns
/// - `vl`: Values of $\Lambda_{v_i}(x)$ for
///   $v_i = v-\lfloor v\rfloor, v-\lfloor v\rfloor+1, \ldots, v$
/// - `dl`: Derivatives $\Lambda_{v_i}\'(x)$ for
///   $v_i = v-\lfloor v\rfloor, v-\lfloor v\rfloor+1, \ldots, v$
///
/// Both vectors have length $\lfloor v\rfloor+1$, i.e. `(v as usize) + 1`.
/// They are all NaN if `x` is NaN or if $|x| > 10^9$.
///
/// # See also
/// - [`bessel_j`]
/// - [`gamma`](crate::gamma)
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
/// # Panics
/// - Panics if `v` is negative
/// - Panics if `v` is too large to convert to `c_int`
#[doc(alias = "lmbda")]
#[must_use]
pub fn jahnke_emden_lambda<V: Into<f64>>(v: V, x: f64) -> (Vec<f64>, Vec<f64>) {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L2008-L2055
    let v: f64 = v.into();
    assert!(v >= 0.0, "v must be non-negative");

    let n = v.to_usize().unwrap();
    // checked before allocating, since `lamv` casts `v` to `int` as well
    assert!(core::ffi::c_int::try_from(n).is_ok(), "v is too large");
    if x.is_nan() || x.abs() > 1e9 {
        // Unlike in SciPy, this avoids out-of-bounds indexing in `lamn` and `lamv`, where `int`
        // casts of `|x|` (e.g. in `msta1` and `msta2`) are undefined for NaN and
        // `|x| >= 2**31 / 1.1`.
        return (vec![f64::NAN; n + 1], vec![f64::NAN; n + 1]);
    }
    let v0 = v - n.to_f64().unwrap();
    let n1 = if n < 1 { 1 } else { n };
    let v1 = v0 + n1.to_f64().unwrap();

    // one more order, for the derivative of the highest
    let size = v1.to_usize().unwrap() + 2;
    let (mut vl, mut dl) = (vec![f64::NAN; size], vec![f64::NAN; size]);
    if v == n.to_f64().unwrap() {
        let mut nm = 0;
        let n = v1.to_usize().unwrap().try_into().unwrap();
        specfun::lamn(n, x, &mut nm, &mut vl, &mut dl);
    } else {
        let mut vm = 0.0;
        specfun::lamv(v1, x, &mut vm, &mut vl, &mut dl);
    }

    vl.resize(n + 1, f64::NAN);
    dl.resize(n + 1, f64::NAN);
    (vl, dl)
}
