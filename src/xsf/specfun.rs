//! The specfun ports (`xsf/specfun/specfun.h`, in the submodules) and their wrappers
//! (`xsf/specfun.h`).
//!
//! Unlike in xsf, the specfun ports compute some `int` products that can overflow (e.g. `n * n`
//! for `n > 46340`) in `f64`, as `f64::from(n) * f64::from(n)`, which is exact for products
//! below 2^53.

mod bessel;
mod erf;
mod gamma;
mod hyp1f1_recurrence;
mod hyperg;
mod legendre;
mod mathieu;
mod numbers;
mod par_cyl;
mod sphd;

pub(crate) use bessel::{cyzo, jdzo, jynb, jyzo, lamn, lamv, msta2, sphj, sphy};
pub(crate) use erf::cerzo;
pub(crate) use gamma::{cgama, gamma2, psi_spec};
pub(crate) use hyperg::{
    cchg, chgm, chgm_neg_x_asymptotic, chgu, hyp1f1_neg_b_fallback, hyp1f1_small_b,
};
pub(crate) use legendre::{lpmns, lpmv, lqmns};
pub(crate) use mathieu::{characteristic_value, coefficient_cv, fcoef, mtu0, mtu12};
pub(crate) use numbers::{bernob, eulerb};
pub(crate) use par_cyl::cpbdn;
pub(crate) use sphd::{aswfa, rswfo, rswfp, segv};

// Newton step for g = f / prod(z - zi): g / g' = f / (df - s * f), where
// s = sum(1 / (z - zi)) is the logarithmic derivative of the product (no overflow).
pub(crate) fn deflated_newton_step(z: C64, f: C64, df: C64, zeros: &[C64]) -> C64 {
    let s = zeros.iter().fold(C64::ZERO, |s, &zi| s + 1.0 / (z - zi));
    f / (df - s * f)
}

/// The `1e300` that specfun returns in place of infinity, e.g. `gamma2` at the poles
pub(crate) const SPECFUN_INF: f64 = 1e300;

/// The magnitude above which some backward recurrences rescale (by this power of 2), since they
/// would overflow otherwise
pub(crate) const SPECFUN_BIG: f64 = pow2(SPECFUN_BIG_LOG2);
pub(crate) const SPECFUN_BIG_LOG2: i32 = 512;

/// `SPECFUN_CONVINF`: specfun's `+-1e300` sentinels as signed infinities
pub(crate) fn specfun_convinf(x: f64) -> f64 {
    if x.abs() == SPECFUN_INF {
        f64::INFINITY.copysign(x)
    } else {
        x
    }
}

/// `SPECFUN_ZCONVINF` (only the real part is checked)
pub(crate) fn specfun_zconvinf(z: &mut Complex<f64>) {
    z.re = specfun_convinf(z.re);
}

/// `specfun::Status`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    Ok,
    NoMemory,
    Other,
}

/// `new (std::nothrow) T[n]()`, i.e. `None` if the allocation fails
pub(crate) fn try_zeroed(n: usize) -> Option<Vec<f64>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).ok()?;
    v.resize(n, 0.0);
    Some(v)
}

use num_complex::Complex;

use crate::utils::{is_nonpos_int, pow2};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble, frexp_f64, ldexp_f64, two_sum};

pub trait Hyp1F1Arg: crate::sealed::Sealed {
    fn hyp1f1(self, a: f64, b: f64) -> Self;
}

impl Hyp1F1Arg for f64 {
    #[inline]
    fn hyp1f1(self, a: f64, b: f64) -> Self {
        // `xsf::hyp1f1`
        let hg = hyperg::chgm_long_poly(a, b, self).unwrap_or_else(|| chgm(self, a, b));
        let hg = hyp1f1_neg_b_fallback(a, b, self, hg);
        hyp1f1_small_b(a, b, c64(self, 0.0), c64(hg, 0.0)).re
    }
}

impl Hyp1F1Arg for Complex<f64> {
    #[inline]
    fn hyp1f1(self, a: f64, b: f64) -> Self {
        // `xsf::hyp1f1`
        self.via(|z| {
            if a.is_nan() || b.is_nan() || z.is_nan() {
                return C64::NAN;
            }
            if let Some(out) = hyperg::hyp1f1_pole_poly(a, b, z) {
                return out;
            }
            let mut outz = hyp1f1_small_b(a, b, z, cchg(a, b, z));
            if outz.re == SPECFUN_INF {
                outz.re = f64::INFINITY;
            }
            outz
        })
    }
}

/// Bernoulli numbers $B_0, \ldots, B_{N-1}$
///
/// Corresponds to [`scipy.special.bernoulli`][scipy] in SciPy, and uses a translation of
/// `xsf::specfun::bernob`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bernoulli.html
///
/// # Examples
/// ```
/// use xsf::bernoulli;
/// assert_eq!(bernoulli::<0>(), []);
/// assert_eq!(bernoulli::<1>(), [1.0]);
/// assert_eq!(bernoulli::<2>(), [1.0, -0.5]);
/// assert_eq!(bernoulli::<4>(), [1.0, -0.5, 1.0 / 6.0, 0.0]);
/// assert_eq!(bernoulli::<1000>()[999], 0.0);
/// ```
///
/// # Panics
/// - If `N` is so large that `N - 1` overflows `i32`.
///
/// # See also
/// - [`euler`]: Euler numbers $E_0,\dotsc,E_{N-1}$
#[must_use]
#[inline]
pub fn bernoulli<const N: usize>() -> [f64; N] {
    let mut bn = [0.0; N];
    if N < 3 {
        // `bernob` requires `N >= 3` (i.e. `n >= 2`), so handle small `N` manually
        if N >= 1 {
            bn[0] = 1.0;
            if N >= 2 {
                bn[1] = -0.5;
            }
        }
    } else {
        bernob((N - 1).try_into().unwrap(), &mut bn);
    }
    bn
}

/// Euler numbers $E_0, \ldots, E_{N-1}$
///
/// Corresponds to [`scipy.special.euler`][scipy] in SciPy, and uses a translation of
/// `xsf::specfun::eulerb`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.euler.html
///
/// # Examples
/// ```
/// use xsf::euler;
/// assert_eq!(euler::<0>(), []);
/// assert_eq!(euler::<1>(), [1.0]);
/// assert_eq!(euler::<2>(), [1.0, 0.0]);
/// assert_eq!(euler::<4>(), [1.0, 0.0, -1.0, 0.0]);
/// assert_eq!(euler::<1000>()[999], 0.0);
/// ```
///
/// # Panics
/// - If `N` is so large that `N - 1` overflows `i32`.
///
/// # See also
/// - [`bernoulli`]: Bernoulli numbers $B_0,\dotsc,B_{N-1}$
#[must_use]
#[inline]
pub fn euler<const N: usize>() -> [f64; N] {
    let mut en = [0.0; N];
    if N < 3 {
        // `eulerb` requires `N >= 3` (i.e. `n >= 2`), so handle small `N` manually
        if N >= 1 {
            en[0] = 1.0;
        }
    } else {
        eulerb((N - 1).try_into().unwrap(), &mut en);
    }
    en
}

/// Kummer's confluent hypergeometric function $_1F_1$
///
/// Corresponds to [`scipy.special.hyp1f1`][scipy] in SciPy, and accepts both `f64` and
/// `num_complex::Complex<f64>` inputs for `z`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp1f1.html
///
/// # Definition
///
/// $$
/// M(a,b,z) = \hyp 1 1 a b z = \sum\_{n=0}^\infty  {\rpow a n \over \rpow b n} {z^n \over n!}
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1].
///
/// For complex $z$ and integers $b \le a \le 0$, the polynomial extension has degree $-a$.
///
/// Results at $b = 0$ can be unreliable. The function may return NaN for ill-conditioned real
/// inputs or when evaluating a terminating series at complex $z$.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.2#E2>
///
/// # See also
/// - [`hypu`]: Tricomi's confluent hypergeometric function $U(a,b,x)$
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function,
///   $_0F_1\left[b\middle\| z\right]$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
#[inline]
pub fn hyp1f1<T: Hyp1F1Arg>(a: f64, b: f64, z: T) -> T {
    z.hyp1f1(a, b)
}

/// `xsf::hypu`
fn xsf_hypu(a: f64, b: f64, x: f64) -> f64 {
    xsf_hypu_digits(a, b, x).0
}

/// `xsf_hypu` with the numbers of digits that `chgu` estimates, as its methods compare them and
/// for the result (more if by the recurrence in a)
pub(crate) fn xsf_hypu_digits(a: f64, b: f64, x: f64) -> (f64, i32, i32) {
    let mut out;
    let mut md = 0; /* method code --- not returned */
    let mut isfer = 0;
    let (mut digits, mut digits_rec) = (0, i32::MIN);

    out = chgu(x, a, b, 0.0, &mut md, &mut isfer, &mut digits, &mut digits_rec);
    // Unlike in xsf, if a - b + 1 rounds to an integer n <= 0 that it isn't (only for b < 1, see
    // `aa_rounding`), U is also computed with the rounding error d. That result replaces the
    // original evaluation if it has enough digits to be returned, and the two disagree by the test
    // of step 3 of `chgu`, but allowing ten times as much disagreement, as the estimate isn't a
    // bound (see `digits_lacking`). This test is heuristic.
    let (rounded, d) = hyperg::aa_rounding(a, b);
    if rounded {
        let (mut md_d, mut isfer_d, mut digits_d, mut digits_rec_d) = (0, 0, 0, i32::MIN);
        let u = chgu(x, a, b, d, &mut md_d, &mut isfer_d, &mut digits_d, &mut digits_rec_d);
        if isfer_d == 0 && (u - out).abs() > 10f64.powi(2 - digits_d) * u.abs() {
            (out, isfer, digits, digits_rec) = (u, isfer_d, digits_d, digits_rec_d);
        }
    }
    if out == SPECFUN_INF {
        out = f64::INFINITY;
    }
    if isfer != 0 {
        out = f64::NAN;
    }
    (out, digits, digits.max(digits_rec))
}

#[cold]
fn hypu_scaled_rising(b: f64, m: u32) -> f64 {
    const SHIFT: i32 = 512;
    let mut r = DoubleDouble::new(pow2(-SHIFT));
    let b = DoubleDouble::new(b);
    for k in 0..m {
        r = r * (b + f64::from(k));
        if !r.is_finite() || r.hi >= pow2(SHIFT) {
            return f64::INFINITY;
        }
    }
    ldexp_f64(r.hi, SHIFT)
}

// Translated from `xsf/cpu/hyperu.h` (scipy/xsf@470575f).
/// Tricomi's confluent hypergeometric function $U(a,b,x)$
///
/// Corresponds to [`scipy.special.hyperu`][scipy] in SciPy.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyperu.html
///
/// # Notes
///
/// It is defined as the solution to the equation
///
/// $$
/// x \frac{\dd^2 w}{\dd x^2} + (b - x) \frac{\dd w}{\dd x} - aw = 0
/// $$
///
/// which satisfies the property $U(a,b,x) \sim x^{-a}$ as $x \to \infty$.
/// See [^DLMF] for more details.
///
/// # See also
/// - [`hyp1f1`]: Kummer's confluent hypergeometric function $M(a,b,z)$
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.2#E6>
#[doc(alias = "hyperu")]
#[must_use]
pub fn hypu(a: f64, b: f64, x: f64) -> f64 {
    if a.is_nan() || b.is_nan() || x.is_nan() {
        f64::NAN
    } else if x < 0.0 {
        // Domain error
        f64::NAN
    } else if x == 0.0 {
        if b > 1.0 {
            // Odd floor(a) gives (-1)^m for a = -m and the sign of 1 / Gamma(a) otherwise.
            let fa = a.floor();
            let finite = a.is_finite();
            let odd = a < 0.0 && finite && fa % 2.0 != 0.0;
            let mut r = f64::INFINITY;
            if a <= 0.0 && a == fa && finite {
                // DLMF 13.2.7: U(-m, b, 0) = (-1)^m (b)_m.
                // (b)_m > m! overflows for m >= 171.
                r = 1.0;
                let m = (-a).min(171.0) as u32;
                for k in 0..m {
                    r *= b + f64::from(k);
                    if !r.is_finite() {
                        break;
                    }
                }
                // Unlike in xsf, use a scaled double-double product near overflow. The guard
                // allows for up to m rounded factors and m rounded multiplications.
                let margin = 2.0 * f64::from(m) * f64::EPSILON;
                if m < 171 && b.is_finite() && r >= f64::MAX * (1.0 - margin) {
                    r = hypu_scaled_rising(b, m);
                }
            }
            return if odd { -r } else { r };
        }

        // Unlike in xsf, form (b)_m directly for b <= 1 to preserve small factors.
        if is_nonpos_int(a) {
            // (-1)^m (b)_m for a = -m (DLMF 13.2.7), i.e. (-1)^m b^m for infinite b
            let mut p = if b.is_infinite() {
                b.powf(-a)
            } else {
                crate::pow_rising(b, -a)
            };
            if (p == 0.0 || !p.is_finite()) && b.is_finite() && !is_nonpos_int(b) {
                // `pow_rising` multiplies from b + m - 1 down, so unlike the product from b up,
                // it can overflow before the factors near 0, e.g. for U(-172, 1/64, 0) = 2.1e307,
                // and for a subnormal b, it can give 0 or NaN
                p = 1.0;
                let mut j = 0.0;
                while j < -a && p.is_finite() && p != 0.0 {
                    p *= b + j;
                    j += 1.0;
                }
                // an overflow has the sign of the ceil(-b) negative factors, not of those so far
                if p.is_infinite() {
                    p = f64::INFINITY * crate::cospi((-b).ceil().clamp(0.0, -a));
                }
            }
            crate::cospi(a) * p
        } else if b >= 1.0 {
            // Unlike in xsf, use the sign of 1 / Gamma(a) at b = 1.
            f64::INFINITY.copysign(crate::gammasgn(a))
        } else {
            // Unlike in xsf, use the falling factorial to avoid rounding a - b + 1 onto a pole.
            crate::pow_falling(-b, -a)
        }
    } else {
        // Unlike in xsf, not by DLMF 13.3.7 for b = 1, small a, and x < 1, which cancels near the
        // zeros of U (`chgu` is accurate there)
        xsf_hypu(a, b, x)
    }
}

/// `U(a, b, x)` for `a < 1` and `b > a` by the backward recurrence in `a` (DLMF 13.3.7), from
/// `a + m` in `[1, 2)` (or `a + j` nearest 0, see below), with its estimated (fractional) number of
/// digits as `chgu` compares it and for the result; NaN if a (heuristic) check finds that the
/// recurrence is ill-conditioned; for small `a` by `hypu_sum_b`, unless that cancels
///
/// `chgu` calls this only for `a < 1`, and the starting values have `a >= 1`, so their `xsf_hypu`
/// (i.e. `chgu`) calls don't come back here (or only once, if `a + m` rounds to just below 1). A
/// start `a + j` in `(-0.5, 0.5]` comes back here, but doesn't start near 0 again.
///
/// Unlike in xsf, `a + a_lo` (for `b < 1`, `a - b + 1 + a_lo`) is the value before rounding.
fn hypu_rec(a: f64, a_lo: f64, b: f64, x: f64) -> (f64, f64, f64) {
    if b < 1.0 {
        // DLMF 13.2.40; the recurrence is better behaved for b >= 1
        let (u, id, id_result) = hypu_rec(a - b + 1.0, a_lo, 2.0 - b, x);
        return (x.powf(1.0 - b) * u, id, id_result);
    }
    // the digits for the result, at most `ID_MIN`, as these estimates are not bounds
    let capped = |id: f64| id.min(f64::from(hyperg::ID_MIN));
    // Unlike in xsf, DLMF 13.3.9 for small a, where the recurrence in a can cancel without the
    // check below noticing, e.g. for tiny a near an integer b
    if a.abs() < 0.5 && (1.5..1000.0).contains(&b) {
        let (u0, s, id_terms) = hypu_sum_b(a + a_lo, b, x);
        let u = u0 + s;
        if u.is_infinite() || u0.abs() + s.abs() < REC_MAX_CANCEL * u.abs() {
            // the terms are taken to have 15 digits, but for the result have `id_terms`
            let id = digits_left(15.0, u0, s);
            return (u, id, digits_lacking(id, id_terms));
        }
    }
    // U(c - 1, b, x) = s U(c, b, x) - t U(c + 1, b, x), with (s, t) = `coef(c)`
    let coef = |c: f64| (2.0 * c - b + x, c * (c - b + 1.0));
    // The recurrence from U(a + m, b, x) and U(a + m + 1, b, x) with `id` digits, or None if they
    // cancel. They are recurred separately, to detect that, and so is a rounding error of each step
    // (as in `hypu_fwd_b`), and for the digits of the result also with an estimate of the rounding
    // of `coef(c)`.
    let rec = |m: f64, u0: f64, u1: f64, id: f64| {
        let (mut p0, mut p1) = (u0, 0.0);
        let (mut q0, mut q1) = (0.0, u1);
        let (mut e0, mut e1) = (0.0, 0.0);
        let (mut f0, mut f1) = (0.0, 0.0);
        // with the sign of the error so far, which it could cancel otherwise (e.g. after a step
        // that cancels)
        let add =
            |e: f64, rounding: f64, u: f64| e + rounding.copysign(if e == 0.0 { u } else { e });
        for k in (1..=m as i32).rev() {
            let (c, c_err) = two_sum(a, f64::from(k));
            let (c, c_err) = two_sum(c, c_err + a_lo);
            let (s, t) = coef(c);
            let (v0, v1) = (p0 + q0, p1 + q1);
            let rounding = f64::EPSILON * (s * v0).abs().max((t * v1).abs());
            // (the rounding errors of s and t, with those from the rounding `c_err` of c)
            let ds = f64::EPSILON * ((2.0 * c - b).abs() + x.abs()) + 2.0 * c_err.abs();
            let dt = 2.0 * f64::EPSILON * (c * ((c - b).abs() + 1.0)).abs()
                + (2.0 * c - b + 1.0).abs() * c_err.abs();
            (p0, p1) = (s * p0 - t * p1, p0);
            (q0, q1) = (s * q0 - t * q1, q0);
            (e0, e1) = (add(s * e0 - t * e1, rounding, p0 + q0), e0);
            let coefs = ds * v0.abs() + dt * v1.abs();
            (f0, f1) = (add(s * f0 - t * f1, rounding + coefs, p0 + q0), f0);
        }
        let u = p0 + q0;
        let id = digits_left(id, p0, q0).min(error_digits(e0, u));
        (p0.abs() + q0.abs() < REC_MAX_CANCEL * u.abs())
            .then(|| (u, id, id.min(error_digits(f0, u))))
    };
    let m = (1.0 - a).ceil();
    let ((u0, id0, _), (u1, id1, _)) =
        (xsf_hypu_digits(a + m + a_lo, b, x), xsf_hypu_digits(a + m + 1.0 + a_lo, b, x));
    // the starting values are taken to have 15 digits, but for the result have `id0` and `id1`
    let mut best = rec(m, u0, u1, 15.0);
    let mut from_start = best.is_some();
    // Unlike in xsf (scipy/xsf#380), the recurrence in b is used if the one in a cancels (for
    // large b and a << 0), e.g. for U(-60.2, 30.6, 10)
    if best.is_none() && (3.0..1000.0).contains(&b) {
        let (u, id) = hypu_rec_b(a, b, x);
        best = Some((u, id, id));
    }
    if a < -0.5 && (1.5..1000.0).contains(&b) {
        // The step to a + j nearest 0 (the first or second) cancels for a near a negative integer
        // and small x, e.g. by 1e7 for U(-2 - 1e-7, 5 - 1e-10, 0.005). Then the recurrence also
        // starts from U(a + j, b, x) (which may use `hypu_sum_b`), and the result with more digits
        // is kept.
        let j = (-a).round();
        let (mut v0, mut v1, mut cancel) = (u0, u1, 1.0);
        for k in (j as i32 + 1..=m as i32).rev() {
            let (s, t) = coef(a + f64::from(k));
            let (f, g) = (s * v0, t * v1);
            (v0, v1) = (f - g, v0);
            cancel = f.abs().max(g.abs()) / v0.abs();
        }
        if cancel > REC_MAX_CANCEL {
            let (v, id, _) = xsf_hypu_digits(a + j + a_lo, b, x);
            let near_0 = rec(j, v, v1, f64::from(id.min(15)));
            if near_0.is_some_and(|(_, id, _)| best.is_none_or(|(_, id_best, _)| id > id_best)) {
                (best, from_start) = (near_0, false);
            }
        }
    }
    let (u, id, id_coefs) = best.unwrap_or((f64::NAN, 0.0, 0.0));
    if from_start {
        return (u, id, digits_lacking(id_coefs, f64::from(id0.min(id1))));
    }
    (u, id, capped(id))
}

/// The largest cancellation `(|p| + |q|) / |p + q|` of a sum that `hypu_rec` accepts (empirical)
const REC_MAX_CANCEL: f64 = 1e4;

/// The power of 2 by which the recurrences in `b` scale down (exactly) values that can overflow
const LOG2_SCALE: i32 = 512;
const SCALE: f64 = pow2(LOG2_SCALE);

/// The two terms of DLMF 13.3.9 for `|a| < 0.5` and `1.5 <= b < 1000`, i.e. of
/// `U(a, b, x) = U(a, b - n, x) + a sum_j U(a + 1, b - j, x)` for `j < n` with `b - n` nearest 1
/// (as `a + 1 > 0.5`, the `U(a + 1, c, x)` don't come back here), and the estimated digits of
/// the `U` that it starts from and of the sum
fn hypu_sum_b(a: f64, b: f64, x: f64) -> (f64, f64, f64) {
    let n = (b - 1.0).round();
    let (u0, id0, _) = xsf_hypu_digits(a, b - n, x);
    // U(a + 1, c, x) for c = b - n + 1, ..., b, by the forward recurrence in c (DLMF 13.3.8) from
    // c = b - n and b - n + 1 (near 1 and 2, where `chguit` integrates over ln(1 + t))
    let mut c = b - n + 1.0;
    let ((v0, id1, _), (v1, id2, _)) =
        (xsf_hypu_digits(a + 1.0, c - 1.0, x), xsf_hypu_digits(a + 1.0, c, x));
    let mut v = [v0, v1];
    // the sum is sum 2^e, as U(a + 1, c, x) can overflow where a U(a + 1, c, x) doesn't, e.g. for
    // U(0.05, 15.5, 2.5e-21) = 6.4e307
    let (mut sum, mut e) = (v[1], 0);
    // and the rounding errors of the steps, as in `hypu_fwd_b`
    let (mut r, mut r_sum) = ([0.0; 2], 0.0);
    while c < b {
        // U(a', c + 1, x) = ((c + x - 1) U(a', c, x) - (c - a' - 1) U(a', c - 1, x)) / x
        let step = |v: [f64; 2]| ((c + x - 1.0) * v[1] - (c - a - 2.0) * v[0]) / x;
        // scaled down (exactly) above SCALE, which keeps the sum finite, and while the step
        // overflows, e.g. for U(0.001, 3, 1e-155) = 1.0e307 (but not below 1, to not underflow)
        while v[1].is_finite() && (v[1].abs() > SCALE || v[1].abs() > 1.0 && step(v).is_infinite())
        {
            (v, sum, e) = (v.map(|t| t / SCALE), sum / SCALE, e + LOG2_SCALE);
            (r, r_sum) = (r.map(|t| t / SCALE), r_sum / SCALE);
        }
        let (f, g) = ((c + x - 1.0) * v[1], (c - a - 2.0) * v[0]);
        let rounding = f64::EPSILON * f.abs().max(g.abs()) / x;
        (v, r) = ([v[1], (f - g) / x], [r[1], step(r) + rounding]);
        sum += v[1];
        r_sum += r[1];
        c += 1.0;
    }
    let id = f64::from(id0.min(id1).min(id2)).min(error_digits(r_sum, sum));
    (u0, ldexp_f64(a * sum, e), id)
}

/// The digits of a result with `id` estimated digits from values taken to have 15 digits: one
/// less, as the estimate isn't a bound (with `id`, a few results with 13 are 1.02e-13 off), and
/// less those that the values lack by their estimated `id_values` digits, but for those not below
/// `ID_MIN` (unless `id - 1` is), as the estimates of `chgu` are often too low
fn digits_lacking(id: f64, id_values: f64) -> f64 {
    let (id, lack) = (id - 1.0, 15.0 - id_values.min(15.0));
    id.min((id - lack).max(f64::from(hyperg::ID_MIN)))
}

/// `U(a, b, x)` for `a < 1` and `3 <= b < 1000` by the forward recurrence in `b` (DLMF 13.3.8),
/// from `hypu_rec` at `b - n` and `b - n + 1` in `[1, 3)`, with its estimated number of digits
fn hypu_rec_b(a: f64, b: f64, x: f64) -> (f64, f64) {
    hypu_fwd_b(a, b, x, (b - 1.0).floor(), |c| {
        let (u, id, _) = hypu_rec(a, 0.0, c, x);
        (u, id)
    })
}

/// `U(a, b, x)` by the forward recurrence in `b` (DLMF 13.3.8) from `U(a, c, x)` and its
/// estimated number of digits by `start(c)` at `c = b - n` and `b - n + 1`
pub(super) fn hypu_fwd_b(
    a: f64,
    b: f64,
    x: f64,
    n: f64,
    start: impl Fn(f64) -> (f64, f64),
) -> (f64, f64) {
    let ((mut u0, id0), (mut u1, id1)) = (start(b - n), start(b - n + 1.0));
    // (fewer for a subnormal start, as in `chgu`)
    let id_start = subnormal_digits(u0, id0).min(subnormal_digits(u1, id1));
    // To estimate how much the starting values cancel, the recurrence also runs from (1, 0), and
    // its result is scaled by the first starting value at the end, as the part of U from that can
    // overflow where U doesn't, e.g. for U(-164, 9 - 0.016, 0.014) = 3.9e306, where they cancel by
    // 1e5. As the recurrence from (1, 0) can overflow too, e.g. for U(150, 71 + 1e-14, 0.0005) =
    // 5.3e68, both of its values are scaled down by `SCALE` after each step that takes it above.
    let (u_start, mut p0, mut p1, mut p_exp) = (u0, 1.0, 0.0, 0);
    // A rounding error of each step is recurred as well, as it can grow even if the starting values
    // don't cancel, e.g. past b = x for a near a negative integer, where U is nearly the minimal
    // solution. This is an estimate, not a bound, as the propagated errors can cancel; a bound by
    // their absolute values is far too pessimistic where U is dominant.
    let (mut e0, mut e1) = (0.0, 0.0);
    let mut c = b - n + 1.0;
    while c < b {
        // U(a, c + 1, x) = ((c + x - 1) U(a, c, x) - (c - a - 1) U(a, c - 1, x)) / x
        let next = |u0: f64, u1: f64| ((c + x - 1.0) * u1 - (c - a - 1.0) * u0) / x;
        let rounding =
            f64::EPSILON * ((c + x - 1.0) * u1).abs().max(((c - a - 1.0) * u0).abs()) / x;
        (u0, u1) = (u1, next(u0, u1));
        (p0, p1) = (p1, next(p0, p1));
        if p1.abs() > SCALE {
            (p0, p1, p_exp) = (p0 / SCALE, p1 / SCALE, p_exp + LOG2_SCALE);
        }
        (e0, e1) = (e1, next(e0, e1) + rounding);
        c += 1.0;
    }
    // (the fraction of U from the first starting value, p1 2^p_exp u_start / u1, where u_start and
    // u1 can be subnormal or large; no digits if it isn't finite)
    let ((m_start, e_start), (m_u, e_u)) = (frexp_f64(u_start), frexp_f64(u1));
    let w = ldexp_f64(p1 * m_start / m_u, p_exp + e_start - e_u);
    let id = if w.is_finite() {
        digits_left(id_start, w, 1.0 - w)
    } else {
        0.0
    };
    let id = id.min(error_digits(e1, u1));
    (if u1.is_finite() { u1 } else { f64::NAN }, id)
}

/// `id`, but at most the digits of a subnormal `u`, whose rounding error can be a few spacings
/// (up to about 7 in the sum of a subnormal integrand in `chguit`, empirical)
fn subnormal_digits(u: f64, id: f64) -> f64 {
    const SPACINGS: f64 = 8.0;
    if u == 0.0 || u.abs() >= f64::MIN_POSITIVE {
        return id;
    }
    id.min(error_digits(SPACINGS * f64::EPSILON * f64::MIN_POSITIVE, u))
}

/// The number of digits of `u` with an error `e`, or 0 if that's NaN (e.g. for `e = u = 0`)
fn error_digits(e: f64, u: f64) -> f64 {
    let id = -(e / u).abs().log10();
    if id.is_nan() { 0.0 } else { id }
}

/// The estimated number of digits of `p + q`, if `p` and `q` have `id` digits
fn digits_left(id: f64, p: f64, q: f64) -> f64 {
    // the ratio is at least 1, but 0 or NaN (inf / inf) if the sum overflows
    id - ((p.abs() + q.abs()) / (p + q).abs()).log10().max(0.0)
}

/// Associated Legendre function for $|x| \leq 1$
///
/// Corresponds to [`scipy.special.lpmv`][scipy] in SciPy.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lpmv.html
///
/// # See also
/// - [`assoc_legendre_p_all`](crate::assoc_legendre_p_all)
#[must_use]
#[inline]
pub fn pmv(m: i32, v: f64, x: f64) -> f64 {
    // Unlike in xsf, the scaled integer-degree path can return a finite 1e300.
    lpmv(x, m, v)
}
