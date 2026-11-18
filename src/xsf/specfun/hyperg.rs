//! Confluent hypergeometric functions, translated into pure Rust from `xsf/cpu/specfun/specfun.h`
//! (scipy/xsf@6406578).

use core::cell::LazyCell;
use core::f64::consts::{LN_2, PI};
use core::ops::{Div, Mul};

use super::hyp1f1_recurrence::{
    MAX_SERIES_ITERATIONS, hyp1f1_b_recurrence, hyp1f1_neg_a, hyp1f1_neg_b,
};
use super::{SPECFUN_INF, cgama, gamma2, psi_spec};
use crate::ext::{Scaled, add_exp, gamma_quot, ln_dd, ln_gamma, mul_gamma_ratio};
use crate::numpy::sinc;
use crate::utils::{EGAMMA, EPSILON_SQ, SQRT_EPSILON, is_nonpos_int, mul_exp, neg_if_odd, pow2};
use crate::xsf::bessel::cyl_bessel_j_f64;
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::consts::{MACHEP, MAXGAM};
use crate::xsf::cephes::dd_real::{DoubleDouble, frexp_f64, ldexp_f64, two_sum};
use crate::xsf::cephes::rgamma;
use crate::xsf::cexpm1;
use crate::xsf::tools::{int_x86, kahan_step};
use crate::xsf::trig::cephes_sinpi;
use crate::{digamma, exprel, gamma, gammaln, hyp0f1};

/// The cost bound of `hyp1f1_sum`, of the sum in `chgm` for non-integer `b < 0`, and of an
/// overflowing sum in `chgul`
const MAX_TERMS: i32 = 100_000;

/// Above this `|a|`, `hyp1f1_neg_b` is only tried if `chgm` gives NaN, as its initial values can
/// be slow there (e.g. a polynomial of degree `-a`)
const NEG_B_MAX_A: f64 = 1e5;

/// The largest relative error bound that `chgm_neg_x_asymptotic` accepts (its bounds are at least
/// 3x the actual errors, and up to about 100 eps for `|a|, |b| <= 100`)
const NEG_X_ASYMPTOTIC_TOL: f64 = 128.0 * f64::EPSILON;
/// The most terms of DLMF 13.7.2 that `chgm_neg_x_asymptotic` sums
const NEG_X_ASYMPTOTIC_TERMS: u32 = 64;
/// Covers the roundings in the error bound of `chgm_neg_x_asymptotic`. Its parts are sums and
/// products of positive terms, chains of at most 32 operations, where the only differences
/// `1 - 2 sigma^2 >= 1/2`, `1 - nu sigma >= 1 - 1/sqrt(2)` and `1 - sigma >= 1/2` have condition
/// numbers of at most 1, 2.5 and 1, so each part is within `128 u` (`u = eps / 2`). The arguments
/// of `exp` are at most `MAX_BOUND_EXP`, so it adds at most `2^13 u`, and `powf(nu, n)`
/// (`n <= 64`) at most `2^13 u` as well; all well below `2^-30`.
const BOUND_SLACK: f64 = 1.0 + 1.0 / pow2(30);
/// The largest argument of `exp` in the error bound of `chgm_neg_x_asymptotic`
const MAX_BOUND_EXP: f64 = 64.0;

/// The (exact) factor that keeps the terms of `chgm_poly`, `rescale_term` and `chgul` in range
const LOG2_SCALE: i32 = 512;
const SCALE: f64 = pow2(LOG2_SCALE);

/// From this `x0` down, `chgm` applies `e^x0` gradually, as it nears the subnormals (`e^-708.4`)
const GRADUAL_EXP_X0: f64 = -700.0;

/// The term `t_j` of the series DLMF 13.2.2 without its factor `z`, from `t_(j-1)`.
///
/// With `b + (j - 1)`, since `b + j - 1` loses the digits of `b` near 0, -1, -2, ...
fn kummer_term<T: Mul<f64, Output = T> + Div<f64, Output = T>>(t: T, a: f64, b: f64, j: i32) -> T {
    let n = f64::from(j - 1);
    t * (a + n) / (f64::from(j) * (b + n))
}

/// The ratio `q_n = (a + n) |z| / ((n + 1) (b + n))`, with `|q_n| = |t_n+1 / t_n|` for the terms
/// of DLMF 13.2.2, or `None` if the terms can grow again after `t_n`: if `|q_n| >= 1` or
/// `b + n < |z|`, e.g. for tiny `a`. Otherwise they can't, because `(n + 1) (b + n) +- (a + n) |z|`
/// increase for `n >= N`.
fn term_ratio(a: f64, b: f64, abs_z: f64, n: f64) -> Option<f64> {
    let q = (a + n) / (n + 1.0) * (abs_z / (b + n));
    (q.abs() < 1.0 && b + n >= abs_z).then_some(q)
}

/// The series term `t / SCALE^k` as `(t, term, |term|)`, where `t` with `norm(t) = |t|` is first
/// scaled by `SCALE` (raising `k`) if it is below `TINY`, and back once it is `TINY SCALE` again.
///
/// Unlike in xsf (scipy/xsf#544), for `b < 0` the terms of DLMF 13.2.2 can underflow before they
/// grow again near `n = -b`, and that part of the sum would be lost, or rebuilt from a few bits.
fn rescale_term<T: Copy + Default + Mul<f64, Output = T>>(
    mut t: T,
    norm: impl Fn(T) -> f64,
    k: &mut u32,
) -> (T, T, f64) {
    // far enough from the subnormals for the double-double products to be exact
    const TINY: f64 = pow2(-900);
    let mut abs_t = norm(t);
    if (f64::MIN_POSITIVE..TINY).contains(&abs_t) {
        (t, abs_t, *k) = (t * SCALE, abs_t * SCALE, *k + 1);
    } else if *k > 0 && abs_t >= TINY * SCALE {
        (t, abs_t, *k) = (t * (1.0 / SCALE), abs_t * (1.0 / SCALE), *k - 1);
    }
    match *k {
        0 => (t, t, abs_t),
        // below 2^-1075 the term rounds to 0, which is slow to get by an underflow
        1 if abs_t >= pow2(LOG2_SCALE - 1075) => {
            let term = t * (1.0 / SCALE);
            (t, term, norm(term))
        },
        _ => (t, T::default(), 0.0),
    }
}

/// Whether the terms `t_m` (`m > n`) of DLMF 13.2.2 are 0, as for `a = 0, -1, ...`, where it is a
/// polynomial of degree `-a`
fn series_ended(a: f64, n: f64) -> bool {
    is_nonpos_int(a) && n >= -a
}

/// Whether the terms `t_m` (`m > n`) of DLMF 13.2.2 after a negligible `t_n` are negligible: if
/// they are 0 (see `series_ended`) or can't grow again (see `term_ratio`), or, for non-integer
/// `b < 0` where they can near `m = -b`, if they sum to at most `2 |t_n|`.
///
/// Since `|a + m| <= |b + m| + d` with `d = |a - b|`, `|t_m+1 / t_m| <= s = 2 |z| / (n + 1)` if
/// `|b + m| >= d`. So if `s <= 1/2`, they only grow for the (at most `2d + 2`) `m` with
/// `|b + m| < d`, by at most `16^d (1 + d / e)^2` in all, with `e` the distance of `b` to the
/// integers. It thus suffices that the `k` terms before shrink by `(2d + 4) 16^d (1 + d / e)^2`.
///
/// Likewise, `|t_m+1 / t_m| <= h / |b + m|` with `h = c |z|`, if `c >= max(|a + j| / (j + 1), 1)`
/// for all `j >= m`. That is at most 1/2 if `|b + m| >= 2h`, and the (at most `4h + 1`) other
/// terms grow by at most `e^(2h) (1 + 2h / e)^2` in all. So it also suffices that the `k` terms
/// before those shrink by `(4h + 4) e^(2h) (1 + 2h / e)^2`, e.g. for `|z|` far below `-b`.
fn rest_is_negligible(a: f64, b: f64, abs_z: f64, n: f64) -> bool {
    let d = (a - b).abs();
    let s = 2.0 * abs_z / (n + 1.0);
    let k = (-b - d).ceil() - n;
    let e = (b - b.round()).abs();
    let growth = (2.0 * d + 4.0).ln() + d * 16_f64.ln() + 2.0 * (d / e).ln_1p();
    // c(j) bounds |a + m| / (m + 1) for m >= j; the h from c(n) gives a later j, for a smaller h
    let c = |j: f64| ((a + j).abs() / (j + 1.0)).max(1.0);
    let h = c(n.max(-b - 2.0 * c(n) * abs_z)) * abs_z;
    let k_h = (-b - 2.0 * h).ceil() - n;
    let growth_h = (4.0 * h + 4.0).ln() + 2.0 * h + 2.0 * (2.0 * h / e).ln_1p();
    series_ended(a, n)
        || term_ratio(a, b, abs_z, n).is_some()
        || (s <= 0.5 && k > 0.0 && -k * s.ln() >= growth)
        || (k_h > 0.0 && k_h * LN_2 >= growth_h)
}

pub(crate) fn cchg(mut a: f64, b: f64, mut z: C64) -> C64 {
    // ===================================================
    // Purpose: Compute confluent hypergeometric function
    //          M(a,b,z) with real parameters a, b and a
    //          complex argument z
    // Input :  a --- Parameter
    //          b --- Parameter
    //          z --- Complex argument
    // Output:  CHG --- M(a,b,z)
    // Routine called: CGAMA for computing complex ln[Г(x)]
    // ===================================================

    let ci = C64::I;
    let mut a0 = a;
    // NOTE: `a1 = a; z0 = z;` are omitted: they are only used for the dead stores
    // `a = a1; z = z0;` at the end.
    let mut cy0 = C64::ZERO;
    let mut cy1 = C64::ZERO;
    // NOTE: uninitialized in C++, but always written before it is read (the `n` loop below
    // runs at least once).
    let mut chg = C64::ZERO;
    // NOTE: `(int)fabs(b)` is UB for |b| >= 2^31 (see `int_x86`); INT_MIN negates to INT_MIN, so
    // b = -2^31 returns 1e300.
    if (b == 0.0) || (b == f64::from(int_x86(b.abs()).wrapping_neg())) {
        return c64(SPECFUN_INF, 0.0);
    }
    if (a == 0.0) || (z == C64::ZERO) {
        return C64::ONE;
    }
    if a == -1.0 {
        return 1.0 - z / b;
    }
    if a == b {
        return z.exp();
    }
    if a - b == 1.0 {
        return (1.0 + z / b) * z.exp();
    }
    if (a == 1.0) && (b == 2.0) {
        // Unlike in xsf, `expm1` near 0, where `exp(z) - 1` cancels (not beyond, where it can
        // overflow, or cancel in its `cos(y) - 1`)
        let em1 = if z.norm() < 1.0 {
            cexpm1(z)
        } else {
            z.exp() - 1.0
        };
        return em1 / z;
    }
    if (a == f64::from(a as i32)) && (a < 0.0) {
        // NOTE: for a = -2^31, `(int)(-a)` is UB (see `int_x86`), so the loop doesn't run and the
        // result is 1. The C++ `k < (m + 1)` is written as `1..=m`:
        // for m = INT_MAX (a = -(2^31 - 1)) `m + 1` overflows (UB), and g++ -O2 then also runs
        // ~2^31 iterations (i.e. treats it as `k <= m`), giving NaN.
        let m = int_x86(-a);
        let mut cr = C64::ONE;
        chg = C64::ONE;
        for k in 1..=m {
            let n = f64::from(k - 1);
            cr = cr * (a + n) / f64::from(k) / (b + n) * z;
            chg += cr;
        }
    } else {
        let x0 = z.re;
        if x0 < 0.0 {
            a = b - a;
            a0 = a;
            z = -z;
        }
        let mut nl = 0;
        let mut la = 0;
        if a >= 2.0 {
            nl = 1;
            // NOTE: `(int)a` is UB for a >= 2^31 (see `int_x86`). `la + 1` overflows (and wraps, as
            // with g++) for a in [2^31 - 1, 2^31).
            la = int_x86(a);
            a -= f64::from(la.wrapping_add(1));
        }
        // ns = 0;  (dead store)
        for n in 0..=nl {
            if a0 >= 2.0 {
                a += 1.0;
            }
            if (z.norm() < 20.0 + b.abs()) || (a < 0.0) {
                chg = C64::ONE;
                let mut chw = C64::ZERO;
                let mut crg = C64::ONE;
                // Unlike in xsf (scipy/xsf#454), a negligible term only ends the sum if the terms
                // can't grow again (see `term_ratio`). Otherwise the sum up to the first one is used
                // if the sum doesn't end, or if the later terms are below its rounding error.
                let mut early = None;
                let mut ended = false;
                // the sum of |crg|
                let mut rs = 1.0;
                for j in 1..501 {
                    crg = kummer_term(crg, a, b, j) * z;
                    chg += crg;
                    rs += crg.l1_norm();
                    if ((chg - chw) / chg).norm() < 1e-15 {
                        if term_ratio(a, b, z.norm(), f64::from(j)).is_some() {
                            ended = true;
                            break;
                        }
                        early.get_or_insert(chg);
                    }
                    chw = chg;
                }
                if let Some(s) = early {
                    if !ended || (chg - s).l1_norm() <= f64::EPSILON * rs {
                        chg = s;
                    }
                }
            } else {
                // y = 0.0;  (dead store)
                let cg1 = cgama(c64(a, 0.0), 0);
                let cg2 = cgama(c64(b, 0.0), 0);
                let cg3 = cgama(c64(b - a, 0.0), 0);
                let mut cs1 = C64::ONE;
                let mut cs2 = C64::ONE;
                let mut cr1 = C64::ONE;
                let mut cr2 = C64::ONE;
                for i in 1..=8 {
                    cr1 = -cr1 * (a + f64::from(i) - 1.0) * (a - b + f64::from(i))
                        / (z * f64::from(i));
                    cr2 = cr2 * (b - a + f64::from(i) - 1.0) * (f64::from(i) - a)
                        / (z * f64::from(i));
                    cs1 += cr1;
                    cs2 += cr2;
                }
                let y = z.im;
                // DLMF 13.7.2, with the sign of ph z (Re z >= 0 here)
                let ns = if y > 0.0 { 1 } else { -1 };
                let mut cfac = (f64::from(ns) * ci * PI * a).exp();
                if y == 0.0 {
                    cfac = c64((PI * a).cos(), 0.0);
                }
                let chg1 = (cg2 - cg3).exp() * z.powf(-a) * cfac * cs1;
                let chg2 = (cg2 - cg1 + z).exp() * z.powf(a - b) * cs2;
                chg = chg1 + chg2;
            }
            if n == 0 {
                cy0 = chg;
            }
            if n == 1 {
                cy1 = chg;
            }
        }
        if a0 >= 2.0 {
            for _i in 1..la {
                chg = ((2.0 * a - b + z) * cy1 + (b - a) * cy0) / a;
                cy0 = cy1;
                cy1 = chg;
                a += 1.0;
            }
            // Unlike in xsf (scipy/scipy#5349), the series (DLMF 13.2.2) replaces the forward
            // recurrence in a above if they differ by more than twice the series' error bound.
            if b > 0.0 {
                if let Some((s, _)) = hyp1f1_series(a0, b, z, C64::ONE, 0.0).filter(|&(s, err)| {
                    err.is_finite() && (!chg.is_finite() || (chg - s).norm() > 2.0 * err)
                }) {
                    chg = s;
                }
            }
        }
        if x0 < 0.0 {
            chg *= (-z).exp();
        }
    }
    // a = a1;
    // z = z0;
    chg
}

/// The terms of `hyp1f1_sum`: complex, or double-double for about twice the precision.
trait Term: Copy + Default + Mul<f64, Output = Self> {
    /// the relative precision `2 u` of the error bound
    const EPS: f64;
    fn norm(self) -> f64;
    /// `sum += self`, with the Kahan compensation `comp` for `C64`
    fn add_to(self, sum: &mut Self, comp: &mut Self);
}

impl Term for C64 {
    const EPS: f64 = f64::EPSILON;
    fn norm(self) -> f64 {
        C64::norm(self)
    }
    fn add_to(self, sum: &mut Self, comp: &mut Self) {
        kahan_step(sum, comp, self);
    }
}

impl Term for DoubleDouble {
    const EPS: f64 = EPSILON_SQ;
    fn norm(self) -> f64 {
        // `inf` on overflow, where the error terms of the products are NaN
        if self.is_finite() {
            self.hi.abs()
        } else {
            f64::INFINITY
        }
    }
    fn add_to(self, sum: &mut Self, _: &mut Self) {
        *sum = *sum + self;
    }
}

/// `hyp1f1_sum` in complex arithmetic.
fn hyp1f1_series(a: f64, b: f64, z: C64, scale: C64, log_scale: f64) -> Option<(C64, f64)> {
    let next_term = |t, n| t * ((a + n) / (n + 1.0) * (z / (b + n)));
    hyp1f1_sum(a, b, z.norm(), next_term, scale, log_scale)
}

/// `scale * M(a, b, z)` from the series DLMF 13.2.2, with the terms `t_n+1 = next_term(t_n, n)` and
/// `abs_z = |z|`, for `a > -1` or `b < 0`, with an error bound for `a >= 1` and `b > 0` (an
/// estimate for `b < 0`), or `None` if it doesn't converge within 100000 terms (a cost bound).
///
/// The error bound (to first order in `u = T::EPS / 2`) is the sum of:
/// - the tail after the last term `t_N`: at most `|t_N| q_N / (1 - q_N)`, since the term ratios
///   `q_n = (a + n) |z| / ((n + 1) (b + n))` are non-increasing;
/// - the rounding error: at most `8 u sum (n + 1) |t_n|` (`8 n u` per term, `2 u |t_n|` per Kahan
///   step of `C64`).
///
/// A negligible `t_N` is skipped if the terms can grow again (see `term_ratio`), unless the sum
/// overflowed. Until then `rescale_term` keeps them from underflowing.
///
/// For `a >= 0`, `b > 0` and `z >= 0` there is no cancellation, and the partial sums only
/// overflow if the result does. For `-1 < a < 0` only the leading 1 cancels.
///
/// The extra factor `e^log_scale` (for `log_scale < 0`) is applied in steps once the terms are
/// large, so that the result is finite even if `e^log_scale` underflows and `M` overflows:
/// - while `err > 1e250`, which leaves room for the next terms below `MAX ~ 1.8e308`;
/// - in steps of at least `e^-600`, which is normal, with an integer exponent so that `log_scale`
///   stays exact (if `|log_scale| < 2^53`);
/// - the rest in 4 steps, so that (since `|sum| <= MAX`) a step only underflows if the result
///   does, and so that the exponent of a step is exact.
fn hyp1f1_sum<T: Term>(
    a: f64,
    b: f64,
    abs_z: f64,
    next_term: impl Fn(T, f64) -> T,
    scale: T,
    mut log_scale: f64,
) -> Option<(T, f64)> {
    // the term `tn` as scaled `k` times by `rescale_term`
    let (mut t, mut k) = (scale, 0);
    let mut sum = scale;
    let mut comp = T::default();
    let mut err = t.norm();
    for n in 0..MAX_TERMS {
        let n = f64::from(n);
        let (mut tn, mut norm);
        (t, tn, norm) = rescale_term(next_term(t, n), T::norm, &mut k);
        tn.add_to(&mut sum, &mut comp);
        err += (n + 2.0) * norm;
        if log_scale < 0.0 && err > 1e250 {
            let s = log_scale.max(-600.0);
            log_scale -= s;
            let s = s.exp();
            (t, tn, sum, comp, err) = (t * s, tn * s, sum * s, comp * s, err * s);
            norm = tn.norm();
        }
        if norm <= T::EPS * sum.norm() {
            let q = term_ratio(a, b, abs_z, n + 1.0);
            if q.is_none() && sum.norm().is_finite() {
                continue;
            }
            let tail = q.map_or(f64::INFINITY, |q| norm * q / (1.0 - q));
            let bound = 4.0 * T::EPS * err + tail;
            if log_scale == 0.0 {
                return Some((sum, bound));
            }
            let s = (log_scale / 4.0).exp();
            return sum
                .norm()
                .is_finite()
                .then(|| (sum * s * s * s * s, bound * s * s * s * s));
        }
    }
    None
}

/// Unlike in xsf: `hg` if finite, else `M(a, b, z)` for `0 < |b| < 1` from `hyp1f1_series` with
/// the first term `|b|`, as the next ones are `~ 1/b` and can overflow first. For `Re z < 0` via
/// DLMF 13.2.39, unless rounding `b - a` moves it relatively far (over `sqrt(eps)`) from an
/// integer `<= 0`. Also `hg` if the sum loses over half of its digits or underflows.
pub(crate) fn hyp1f1_small_b(a: f64, b: f64, z: C64, hg: C64) -> C64 {
    if hg.is_finite() || !(a.is_finite() && z.is_finite() && b != 0.0 && b.abs() < 1.0) {
        return hg;
    }
    let (a, z, ln_scale, phase) = if z.re < 0.0 {
        let ak = DoubleDouble::new(b) - a;
        if rounding_ratio(ak) > SQRT_EPSILON {
            return hg;
        }
        (ak.hi, -z, z.re, z.im)
    } else {
        (a, z, 0.0, 0.0)
    };
    let scale = b.abs() * C64::cis(phase);
    hyp1f1_series(a, b, z, scale, ln_scale)
        .filter(|&(m, err)| m.norm().is_normal() && err <= SQRT_EPSILON * m.norm())
        .map_or(hg, |(m, _)| m / b.abs())
}

/// The polynomial `M(-n, b, x)` for `n >= 0` and `b > 0`, using the recurrence of
/// `xsf::detail::eval_genlaguerre_l`, since `M(-n, b, x) = L_n^(b-1)(x) / binom(n+b-1, n)` (DLMF
/// 13.6.19).
///
/// The recurrence is linear in `(d, p)`, so they are scaled up by `2^512` (exactly) whenever both
/// are below `2^-512`, far from the subnormals, where the arithmetic is inexact and slow, and down
/// whenever one is above `2^512` (e.g. `~ 1/b` for tiny `b`). Returns `scale M`, as `M` can
/// overflow when `scale M` doesn't.
fn chgm_poly(n: i32, b: f64, x: f64, scale: f64) -> f64 {
    let (mut d, mut p) = (0.0, 1.0);
    // `M = p 2^-e`
    let mut e = 0_i32;
    for k in 0..n {
        let k = f64::from(k);
        d = (d * k - p * x) / (k + b);
        p += d;
        let max = d.abs().max(p.abs());
        if max * SCALE < 1.0 {
            (d, p, e) = (d * SCALE, p * SCALE, e.saturating_add(LOG2_SCALE));
        } else if max > SCALE {
            (d, p, e) = (d / SCALE, p / SCALE, e.saturating_sub(LOG2_SCALE));
        }
    }
    ldexp_f64(p * scale, -e)
}

/// `hyp1f1_series` for `M(a, b, x) e^min(x0, 0)` (cf. DLMF 13.2.39), where `e^x0` is applied
/// gradually for `x0 <= -700`, to avoid underflow
fn chgm_series(a: f64, b: f64, x: f64, x0: f64) -> Option<(f64, f64)> {
    let s = if x0 > GRADUAL_EXP_X0 {
        x0.min(0.0)
    } else {
        0.0
    };
    hyp1f1_series(a, b, c64(x, 0.0), c64(s.exp(), 0.0), x0.min(0.0) - s)
        .map(|(hg, err)| (hg.re, err))
}

/// `chgm_series` in double-double, with `log_scale = min(x0, 0)`.
fn chgm_series_dd(a: DoubleDouble, b: DoubleDouble, x: f64, x0: f64) -> Option<(f64, f64)> {
    let one = DoubleDouble::new(1.0);
    let next_term = |t, n| t * (a + n) * x / ((b + n) * (n + 1.0));
    hyp1f1_sum(a.hi, b.hi, x, next_term, one, x0.min(0.0)).map(|(hg, err)| (hg.hi, err))
}

/// Whether the terms `t_m` (`m >= n`) of DLMF 13.2.2 with `z > 0` have one sign: `t_(k+1) / t_k`
/// is negative only for the integers `k` between `-a` and `-b` (and `t_(k+1) = 0` for `k >= -a`
/// if `a = 0, -1, ...`).
fn terms_keep_sign(a: f64, b: f64, n: i32) -> bool {
    let poly = is_nonpos_int(a);
    let end = if poly { -a } else { (-a).max(-b) };
    f64::from(n).max((-a).min(-b).floor() + 1.0) >= end
}

/// Whether the value of `(value, error bound)` is finite and keeps at least half of its digits
fn half_accurate(&(hg, err): &(f64, f64)) -> bool {
    hg.is_finite() && err <= SQRT_EPSILON * hg.abs()
}

/// The largest relative error of `a.hi + k` for `a + k` over the integers `0 <= k < MAX_TERMS`, as
/// in the factors `(a.hi)_n` of the series terms that `hyp1f1_sum` can reach; infinite if one is 0
/// but its `a + k` isn't
fn rounding_ratio(a: DoubleDouble) -> f64 {
    if a.lo == 0.0 {
        return 0.0;
    }
    let last = f64::from(MAX_TERMS - 1);
    let dist = if a.hi > 0.0 {
        a.hi
    } else if a.hi < -last {
        -a.hi - last
    } else {
        (a.hi - a.hi.round()).abs()
    };
    a.lo.abs() / dist
}

/// `(a, x)`, or `(b - a, -x)` of DLMF 13.2.39 for `x < 0`, with `b - a` also without rounding, e.g.
/// to the integer -100 for `hyp1f1(99.7, -0.3, -800)`
fn kummer_args(x: f64, a: f64, b: f64) -> (DoubleDouble, f64) {
    if x < 0.0 {
        (DoubleDouble::new(b) - a, -x)
    } else {
        (DoubleDouble::new(a), x)
    }
}

/// `e^min(x0, 0) M(a, b, x)` for `a, x > 0` and non-integer `b < 0` as `K + I` (DLMF 13.2.42), if
/// that loses at most half of its digits, with
/// - `K = Gamma(a - b + 1) / Gamma(1 - b) U(a, b, x)`,
/// - `I = -Gamma(a - b + 1) Gamma(b - 1) / (Gamma(1 - b) Gamma(a)) x^(1 - b) M2`,
/// - `M2 = M(a - b + 1, 2 - b, x)` and `M3 = M(a - b + 2, 3 - b, x)`, whose terms are positive.
///
/// `K` is the minimal solution of DLMF 13.3.1 for increasing `a` (DLMF 13.8.11, 13.8.12), so
/// `r = K(a + 1) / K(a)` is a continued fraction. The Wronskian of `K` and `x^(1 - b) M2` is that
/// of `M`, i.e. `(1 - b) x^-b e^x` (DLMF 13.2.33), so by DLMF 13.3.15 and 13.3.24,
/// `K = (1 - b) e^x / ((1 - b + a (1 - r)) M2 + x (a - b + 1) / (2 - b) M3)`, where `0 < r < 1`.
fn chgm_minimal(a: DoubleDouble, b: f64, x: f64, x0: f64) -> Option<f64> {
    // a (1 - r) = a / g, with g = 1 + (a_1 - b) / (x + a_1 / (1 + (a_2 - b) / (x + ...))) and
    // a_n = a + n, by the modified Lentz method. Its terms are positive, so nothing cancels, and
    // successive convergents enclose g.
    let (mut g, mut c, mut d) = (1.0, 1.0, 0.0);
    let n = (1..=MAX_TERMS).find(|&n| {
        let an = a.hi + f64::from(n);
        for (p, q) in [(an - b, x), (an, 1.0)] {
            (c, d) = (q + p / c, 1.0 / (q + p * d));
            g *= c * d;
        }
        (c * d - 1.0).abs() <= f64::EPSILON
    })?;
    // e^-x M2 and e^-x M3, with their absolute errors
    let (a2, b1) = (a - b + 1.0, DoubleDouble::new(1.0) - b);
    let (m2, err2) = chgm_series_dd(a2, b1 + 1.0, x, -x).filter(half_accurate)?;
    let (m3, err3) = chgm_series_dd(a2 + 1.0, b1 + 2.0, x, -x).filter(half_accurate)?;
    let ar = a.hi / g;
    let den = (b1.hi + ar) * m2 + x * a2.hi / (b1.hi + 1.0) * m3;
    // The error of g, which only counts for the share of its term in `den`: 1 eps for the
    // truncation, and at most 16 roundings per step, as the products of the c and of the d are
    // ratios of sums of products of the rounded terms, off by up to 5 and 6 roundings per step,
    // and multiplying them into g adds 4.
    let err_g = (1.0 + 8.0 * f64::from(n)) * f64::EPSILON * (ar * m2 / den);
    // (and about 10 roundings for the rest)
    let err_k = err2 / m2 + err3 / m3 + 5.0 * f64::EPSILON + err_g;
    let k = Scaled::new(b1.hi / den, 0, err_k)?;
    let i = Scaled::new(-m2, 0, err2 / m2)?.with_l(ln_dd(DoubleDouble::new(x)) * b1 + x0.max(0.0));
    let i = mul_gamma_ratio([a2, -b1], [b1, a], || Some(i))?;
    let (hg, err) = add_exp([k.with_l(DoubleDouble::new(x0.min(0.0))), i])?.to_f64();
    (err <= SQRT_EPSILON).then_some(hg)
}

/// Unlike in xsf: `hg` unless it is NaN, or else `M(a, b, x)` for non-integer `b < 0`
/// - from `chgm_minimal`, with `a` and `x` as in `kummer_args` if both are positive, where the
///   series of `M` cancels and DLMF 13.3.2 is unstable in both directions, e.g. for
///   `hyp1f1(316.3, -1000.3, 100)`;
/// - or else for finite `|a| > NEG_B_MAX_A` and `x` from `hyp1f1_neg_b`.
pub(crate) fn hyp1f1_neg_b_fallback(a: f64, b: f64, x: f64, hg: f64) -> f64 {
    if !(hg.is_nan() && b < 0.0 && b != b.trunc()) {
        return hg;
    }
    let (ak, xk) = kummer_args(x, a, b);
    // the sums of M2 and M3 need over x terms
    let minimal = ak > 0.0 && 0.0 < xk && xk < f64::from(MAX_TERMS);
    let huge_a = a.abs() > NEG_B_MAX_A && a.is_finite() && x.is_finite();
    let neg_b = || hyp1f1_neg_b(a, b, x, |b| chgm(x, a, b));
    (minimal.then(|| chgm_minimal(ak, b, xk, x)).flatten())
        .or_else(|| huge_a.then(neg_b).flatten())
        .unwrap_or(hg)
}

/// `M(a, b, x)` for `x < 0` from `chgm_neg_x_series`, or else (for `x <= -30 - |b|`) from
/// `chgm_neg_x_asymptotic`
fn chgm_neg_x(x: f64, a: f64, b: f64) -> Option<f64> {
    chgm_neg_x_series(x, a, b).or_else(|| chgm_neg_x_asymptotic(x, a, b))
}

/// `M(a, b, x)` for `x <= -30 - |b|`, or for `x < 0` and `b > 0` if rounding `b - a` errs by more
/// than a rounding (see `rounding_ratio`), from the series of DLMF 13.2.39 in double-double with
/// the exact `b - a`, if that is a non-integer below 0 and the error bound is below `eps`.
///
/// Near a pole of `Gamma(b - a)`, whose reciprocal scales the `e^-x` part, rounding `b - a` and the
/// recurrence in `a` lose all digits, and even the sign. The series needs over `-x` terms.
fn chgm_neg_x_series(x: f64, a: f64, b: f64) -> Option<f64> {
    let ak = DoubleDouble::new(b) - a;
    // (for `b < 0` see `chgm_neg_b`, where the series can need over `MAX_TERMS` terms)
    let near = 0.0 < -x && b > 0.0 && rounding_ratio(ak) > f64::EPSILON;
    if !(near || -x >= 30.0 + b.abs()) || -x >= f64::from(MAX_TERMS) {
        return None;
    }
    if ak >= 0.0 || (ak.hi.fract() == 0.0 && ak.lo.fract() == 0.0) {
        return None;
    }
    let (hg, err) = chgm_series_dd(ak, DoubleDouble::new(b), -x, x)?;
    (hg.is_finite() && err <= f64::EPSILON * hg.abs()).then_some(hg)
}

/// `M(a, b, x)` for `x <= -30 - |b|` and `-x > 2 |b - 2a|`, by the asymptotic expansion of
/// DLMF 13.7.2, if its error bound is at most `NEG_X_ASYMPTOTIC_TOL`. DLMF 13.2.41, averaged over
/// both signs, gives, with `z = -x` and `q = b - a`,
///
/// ```text
/// M(a, b, -z) = Gamma(b) / Gamma(q) Re[e^(-i pi a) U(a, b, z e^(-i pi))]
///     + Gamma(b) / Gamma(a) cos(pi q) e^-z U(q, b, z)
/// ```
///
/// The first term is `P S_n` plus the remainder of DLMF 13.7.4, with `P = Gamma(b) / Gamma(q) z^-a`
/// and `S_n = sum_(k<n) t_k`, `t_k = (a)_k (1 - q)_k / (k! z^k)`, and the second term is only
/// bounded. Both bounds are those of DLMF 13.7.5, with `sigma = |b - 2a| / z = |b - 2q| / z`: in
/// region R3 (at `ph z = -pi`) for the remainder, and in R1 with `n = 1` for the second term.
fn chgm_neg_x_asymptotic(x: f64, a: f64, b: f64) -> Option<f64> {
    const U: f64 = MACHEP;
    let z = -x;
    // (rounded up, as all bounds increase with it)
    let r = (b - 2.0 * a).abs() * BOUND_SLACK;
    if !(z >= 30.0 + b.abs() && z > 2.0 * r) {
        return None;
    }
    let q = DoubleDouble::new(b) - a;
    let (ln_gb, sign_b, err_gb) = ln_gamma(DoubleDouble::new(b))?;
    let (ln_gq, sign_q, err_gq) = ln_gamma(q)?;
    let (ln_ga, _, err_ga) = ln_gamma(DoubleDouble::new(a))?;
    let ln_z = ln_dd(DoubleDouble::new(z));

    // DLMF 13.7.6, 13.7.8 and 13.7.9 (in R3 with `nu sigma` for `sigma`), bounded from above
    // without cancellation (see `BOUND_SLACK`): `nu <= (1 - 2 sigma^2)^(-1/2)` as `sqrt(y) >= y`
    // for `0 <= y <= 1`, `|2a^2 - 2ab + b| = |2q^2 - 2qb + b| = |b - 2aq| <= |b| + 2 |a q|`, and
    // Gautschi's inequality `chi(n) < sqrt(pi (n + 2) / 2)` (DLMF 5.6.4).
    let sigma = r / z;
    let nu = 1.0 / (1.0 - 2.0 * sigma * sigma).sqrt();
    let rho_0 = (b.abs() + 2.0 * (a * q.hi).abs()) / 2.0;
    let c_n = |n: f64| ((PI * (n + 2.0) / 2.0).sqrt() + sigma * nu * nu * n) * nu.powf(n);
    let remainder = remainder_factor(nu * sigma, rho_0, c_n(1.0), z)?;

    // The terms, each with at most 7 roundings per step (`fc` is within `u (1 + 3u)`, as the
    // double-double plus binary64 sum errs by at most `2u^2` relative to its result), so that
    // `|t_k| <= |t^_k| / (1 - gamma_7k)`, and the error bound `err_s` of their sum `s`.
    let gamma_7 = |k: u32| 7.0 * f64::from(k) * U / (1.0 - 7.0 * f64::from(k) * U);
    let (mut t, mut tail, mut err_s) = (1.0_f64, 0.0_f64, 0.0);
    let (s, err_r) = (1..=NEG_X_ASYMPTOTIC_TERMS).find_map(|n| {
        let k = f64::from(n - 1);
        let (fa, fc) = (a + k, (-q + (k + 1.0)).hi);
        let t_n = t * (fa * fc / ((k + 1.0) * z));
        // the series ends at an exact zero, but an underflow loses the relative error bound
        let underflow = t_n.abs() < f64::MIN_POSITIVE && (t_n != 0.0 || fa != 0.0 && fc != 0.0);
        if underflow || t_n.abs() > t.abs() {
            return Some(None);
        }
        let g = gamma_7(n);
        let s = 1.0 + tail;
        let err_r = remainder * c_n(f64::from(n)) * t_n.abs() / (1.0 - g);
        // (leaving most of the tolerance to the roundings)
        if err_r <= NEG_X_ASYMPTOTIC_TOL * s.abs() / 8.0 {
            return Some(Some((s, err_r)));
        }
        tail += t_n;
        err_s += g / (1.0 - g) * t_n.abs() + U * tail.abs();
        t = t_n;
        None
    })??;
    let err_s = err_s + U * s.abs();

    // `P S^`, with the errors of `l` (and `fold_exp`'s), and the bound `e_p` of `exp(err) - 1`,
    // where `exp(e) - 1 <= e (1 + e)` holds for small `e`
    let l = ln_gb - ln_gq - ln_z * a;
    let err_l =
        err_gb + err_gq + 4.0 * U * U * (ln_gb.hi.abs() + ln_gq.hi.abs() + l.hi.abs() + 1.0);
    let (hg, err_y) = Scaled::new(sign_b * sign_q * s, 0, err_l)?
        .with_l(l)
        .to_f64();
    if !(hg.is_finite() && hg.abs() >= f64::MIN_POSITIVE && err_y <= 1e-10) {
        return None;
    }
    let e_p = err_y * (1.0 + err_y);

    // the bound of the second term relative to `hg`, with `|q (q - b + 1)| = |q (1 - a)|`
    let q1 = q.hi;
    let ln_u1 = (remainder_factor(sigma, rho_0, 1.0, z)? * (q1 * (1.0 - a)).abs() / z).ln_1p();
    let ln_e = [ln_gb.hi, -ln_ga.hi, -z, -q1 * ln_z.hi, ln_u1, -hg.abs().ln()];
    let err_ln_e = err_gb
        + err_ga
        + q.lo.abs() * (ln_z.hi.abs() + 1.0)
        + 8.0 * U * (ln_e.iter().map(|v| v.abs()).sum::<f64>() + 1.0);
    let e_rel = (ln_e.iter().sum::<f64>() + err_ln_e).exp();

    let rel = (e_p + (1.0 + e_p) * (err_s + err_r) / s.abs() + e_rel) * BOUND_SLACK;
    (rel < 1.0 && rel / (1.0 - rel) <= NEG_X_ASYMPTOTIC_TOL).then_some(hg)
}

/// The factor `2 alpha exp(2 alpha rho C_1 / z)` of the remainder bound of DLMF 13.7.5, with
/// `alpha` and `rho` of DLMF 13.7.9 for `s = sigma` (or `nu sigma` in R3) and
/// `rho_0 = |2a^2 - 2ab + b| / 2`, or `None` if the argument of `exp` exceeds `MAX_BOUND_EXP`
fn remainder_factor(s: f64, rho_0: f64, c_1: f64, z: f64) -> Option<f64> {
    let alpha = 1.0 / (1.0 - s);
    let rho = rho_0 + s * (1.0 + s / 4.0) * alpha * alpha;
    let w = 2.0 * alpha * rho * c_1 / z;
    (w <= MAX_BOUND_EXP).then(|| 2.0 * alpha * w.exp())
}

/// `M(a, b, x)` for non-integer `b < 0` from the series, or else from `hyp1f1_neg_b`, or else from
/// the series in double-double, if any loses at most half of its digits. If rounding `b - a` (for
/// `x < 0`) errs by more than a rounding (see `rounding_ratio`), the series in double-double comes
/// before `hyp1f1_neg_b`, which rounds it too, if its error bound is below `eps`. `hyp1f1_neg_b`
/// is skipped for `|a| > NEG_B_MAX_A`, where `chgm` is often accurate (see
/// `hyp1f1_neg_b_fallback`).
fn chgm_neg_b(x: f64, a: f64, b: f64) -> Option<f64> {
    /// The series if it ends, with an error bound below `tol` relative to its value
    fn checked(
        x: f64,
        ends: bool,
        tol: f64,
        ln_max: impl Fn() -> f64,
        series: impl Fn(f64) -> Option<(f64, f64)>,
    ) -> Option<f64> {
        let mut hg = if ends { series(x) } else { None };
        // if the partial sums overflow, `e^-s M` with `s = ceil(ln_max)` is summed instead (an
        // integer, as `exp` amplifies the rounding of a large argument), with the rounding error
        // `d` of `x - s` applied after, as it can be over `eps` relative to `x`
        let (mut ln_scale, mut d) = (0.0, 0.0);
        if hg.is_some_and(|(hg, _)| !hg.is_finite()) {
            ln_scale = ln_max().ceil();
            let x0;
            (x0, d) = two_sum(x.min(0.0), -ln_scale);
            hg = ln_scale.is_finite().then(|| series(x0)).flatten();
        }
        // (with the error of an underflow in the last 4 steps of `hyp1f1_sum`, each by up to half
        // the smallest subnormal `MIN_POSITIVE eps`, as `err` underflows too)
        let underflow = 2.0 * f64::MIN_POSITIVE * f64::EPSILON;
        hg.filter(|&(hg, err)| hg.is_finite() && err + underflow <= tol * hg.abs())
            .map(|(hg, _)| mul_exp(hg, ln_scale) * d.exp())
    }
    let (ak, xk) = kummer_args(x, a, b);
    // whether the series can end with a finite sum, i.e. within `MAX_TERMS` (see `term_ratio`)
    let ends = b + f64::from(MAX_TERMS) >= xk;
    let series_dd = |x0| chgm_series_dd(ak, DoubleDouble::new(b), xk, x0);
    // with an estimate of the error of rounding `b - a`, which changes the term `t_n` relatively by
    // up to about `n r`, while `err` is at least `4 eps sum (n + 1) |t_n|` (see `hyp1f1_sum`)
    let r = rounding_ratio(ak);
    let series = |x0| {
        chgm_series(ak.hi, b, xk, x0).map(|(hg, err)| (hg, err + r / (4.0 * f64::EPSILON) * err))
    };
    let small_rounding = r <= f64::EPSILON;
    // an estimate of `ln max_n |t_n|` for the terms `t_n` of `e^min(x, 0) M(ak, b, xk)`, summing
    // the logarithms of their factors until they shrink (as in `term_ratio`)
    let ln_max = LazyCell::new(|| {
        let (mut ln_t, mut ln_max) = (x.min(0.0), x.min(0.0));
        for n in (0..MAX_TERMS).map(f64::from) {
            let ln_q = (ak + n).hi.abs().ln() - (n + 1.0).ln() + xk.ln() - (b + n).abs().ln();
            if ln_q < 0.0 && b + n >= xk {
                break;
            }
            ln_t += ln_q;
            ln_max = ln_max.max(ln_t);
        }
        ln_max
    });
    let neg_b = || hyp1f1_neg_b(a, b, x, |b| chgm(x, a, b));
    checked(x, ends, SQRT_EPSILON, || *ln_max, series)
        .or_else(|| checked(x, ends && !small_rounding, f64::EPSILON, || *ln_max, series_dd))
        .or_else(|| (a.abs() <= NEG_B_MAX_A).then(neg_b).flatten())
        .or_else(|| checked(x, ends, SQRT_EPSILON, || *ln_max, series_dd))
}

/// `e^(x0/2 - x/2) M(a, b, x)` for `b, x > 0` and `c = b - 2a > 0` from Tricomi's expansion
/// `Gamma(b) e^(x0/2) (t/2)^(1-b) sum_n B_n J_(b-1+n)(t)` (DLMF 13.11.3, 13.11.4), with
/// `t = sqrt(2 c x)`, `B_n = A_n w^n` and `w = sqrt(x / (2c))`. If `b > t` and
/// `Gamma(b) (t/2)^(1-b) > 1`, its terms are `B_n p_n 0F1(; b + n; -t^2/4)` with
/// `p_n = (t/2)^n / (b)_n` (DLMF 10.16.9) instead, which underflow later.
///
/// If `q = ((b + 3) w^2 + c w^3) / 3 <= 1/2`, then for `n >= 2`,
/// `|B_(n+1)| <= q max(|B_n|, |B_(n-1)|, |B_(n-2)|)`, so as `|J_v|, |0F1| <= 1` (DLMF 10.14.1,
/// 10.14.4), the rest after `B_n` is at most `3 p_n` times that maximum, and `sum |B_n| < 5`.
/// `None` if `q > 1/2`, after 500 terms, if underflowing terms aren't negligible, or if the result
/// isn't finite.
fn chgm_tricomi(a: f64, b: f64, x: f64, x0: f64) -> Option<f64> {
    let c = b - 2.0 * a;
    let t = (2.0 * c * x).sqrt();
    let w2 = x / (2.0 * c);
    let cw3 = c * w2 * w2.sqrt();
    let q = ((b + 3.0) * w2 + cw3) / 3.0;
    if !(q <= 0.5) {
        return None;
    }
    let ln_scale = gammaln(b) + (1.0 - b) * (t / 2.0).ln();
    let scaled = b > t && ln_scale > 0.0;
    let j = |v| cyl_bessel_j_f64(v, t);
    let f = |n| {
        if scaled {
            hyp0f1(b + n, -t * t / 4.0)
        } else if n == 0.0 {
            // J_(b-1) = (2b/t) J_b - J_(b+1) (DLMF 10.6.1), since `b - 1` can round off a tiny `b`
            2.0 * b / t * j(b) - j(b + 1.0)
        } else {
            j(b - 1.0 + n)
        }
    };
    let mut sum = f(0.0);
    let mut abs_sum = sum.abs();
    let mut p = 1.0;
    // B_(n-2), B_(n-1) and B_n, with B_0 = 1 and B_1 = 0
    let [mut bn2, mut bn1, mut bn] = [0.0, 1.0, 0.0_f64];
    for n in (1..500).map(f64::from) {
        if scaled {
            p *= t / 2.0 / (b + n - 1.0);
        }
        let term = bn * p * f(n);
        sum += term;
        abs_sum += term.abs();
        if 3.0 * bn.abs().max(bn1.abs()).max(bn2.abs()) * p <= f64::EPSILON * abs_sum {
            // the terms that underflow add less than `5 MIN_POSITIVE`, i.e. `5 EPSILON abs_sum`
            if abs_sum < f64::MIN_POSITIVE / f64::EPSILON {
                return None;
            }
            let e = (x0 / 2.0).exp();
            if scaled {
                return Some(e * sum).filter(|hg| hg.is_finite());
            }
            // the product if in range, as `exp` amplifies the rounding of a large logarithm; if
            // only the prefactor overflows, `sum` first, or else in the exponent
            let [g, pw] = [gamma(b), (t / 2.0).powf(1.0 - b)];
            let ln_scale = ln_scale + x0 / 2.0;
            let scale = if [g, pw, e, g * pw, g * pw * e].iter().all(|v| v.is_normal()) {
                g * pw * e
            } else {
                ln_scale.exp()
            };
            let product = g * sum * pw * e;
            let hg = if scale.is_finite() {
                scale * sum
            } else if product.is_finite() {
                product
            } else {
                (ln_scale + sum.abs().ln()).exp().copysign(sum)
            };
            return Some(hg).filter(|hg| hg.is_finite());
        }
        // (n + 1) B_(n+1) = (n + b - 1) w^2 B_(n-1) - c w^3 B_(n-2)
        (bn2, bn1, bn) = (bn1, bn, ((n + b - 1.0) * w2 * bn1 - cw3 * bn2) / (n + 1.0));
    }
    None
}

/// The series DLMF 13.2.2 as summed by `chgm_sum`, with what `chgm` needs to judge it
struct KummerSum {
    /// the sum, without the factor of DLMF 13.2.39
    hg: f64,
    /// its last term
    rg: f64,
    /// the sum of |rg|, to measure the cancellation
    rs: f64,
    /// the number of terms
    terms: i32,
    /// the number of terms of the last finite partial sum
    last_finite: i32,
    /// whether a negligible term ended the sum
    converged: bool,
}

impl KummerSum {
    /// Whether the term `rg` is negligible in the sum `hg`
    fn negligible(rg: f64, hg: f64) -> bool {
        rg == 0.0 || (rg / hg).abs() < 1e-15
    }

    /// Whether the sum converged and loses at most half of its digits
    fn accurate(&self) -> bool {
        self.converged && self.rs * SQRT_EPSILON <= self.hg.abs()
    }

    /// By how much the sum cancels, or `inf` if it didn't converge
    fn cancellation(&self) -> f64 {
        if self.converged {
            self.rs / self.hg.abs()
        } else {
            f64::INFINITY
        }
    }

    /// Whether the sum of a series with a < -1 and b > 0 doesn't converge or is estimated to lose
    /// more than half of its digits (as in Boost's checked series)
    fn lossy(&self, a: f64, b: f64, x: f64) -> bool {
        a < -1.0 && b > 0.0 && !self.accurate() && x > 0.0
    }

    /// A finite sum is wrong if its last term (or after 500 terms, the rest) isn't negligible, or
    /// for b < 0 if it loses over half of its digits (as in `chgm_neg_b`). For b < 0 an overflow is
    /// wrong if the terms can change sign after it, since then so can the sum.
    fn wrong(&self, a: f64, b: f64, x: f64, neg_b: bool) -> bool {
        let (hg, terms) = (self.hg, self.terms);
        if hg.is_finite() {
            let rest_negligible = || rest_is_negligible(a, b, x, f64::from(terms));
            let truncated = !Self::negligible(self.rg, hg) || terms > 500 && !rest_negligible();
            truncated || b < 0.0 && self.rs * SQRT_EPSILON > hg.abs()
        } else {
            neg_b && !(hg.is_infinite() && terms_keep_sign(a, b, self.last_finite + 1))
        }
    }
}

/// The series of `chgm`: up to `max_terms` terms of DLMF 13.2.2
fn chgm_sum(a: f64, b: f64, x: f64, max_terms: i32) -> KummerSum {
    let mut hg = 1.0;
    let mut rg = 1.0;
    // the term `rg` as scaled `k` times by `rescale_term`
    let (mut t, mut k) = (1.0, 0);
    let mut rs = 1.0;
    let mut converged = false;
    let mut terms = 0;
    let mut last_finite = 0;
    for j in 1..=max_terms {
        terms = j;
        (t, rg, _) = rescale_term(kummer_term(t, a, b, j) * x, f64::abs, &mut k);
        hg += rg;
        if hg.is_finite() {
            last_finite = j;
        }
        rs += rg.abs();
        let rest_negligible = || rest_is_negligible(a, b, x, f64::from(j));
        let ended = || series_ended(a, f64::from(j));
        // Unlike in xsf (scipy/xsf#454), not if the terms can grow again (see `term_ratio`),
        // unless they are 0: all that follow (see `series_ended`), or within 500 terms exactly
        // (`t == 0`) or else with a negligible rest
        if (hg != 0.0)
            && ((rg / hg).abs() < 1e-15)
            && ((rg == 0.0 && (ended() || j <= 500 && (t == 0.0 || rest_negligible())))
                || term_ratio(a, b, x, f64::from(j)).is_some())
        {
            converged = true;
            break;
        }
        // Unlike in xsf (scipy/xsf#365), for non-integer b < 0 (`max_terms > 500`) the sum goes on
        // after 500 terms until the rest is negligible
        if j >= 500 && (!hg.is_finite() || KummerSum::negligible(rg, hg) && rest_negligible()) {
            break;
        }
    }
    KummerSum {
        hg,
        rg,
        rs,
        terms,
        last_finite,
        converged,
    }
}

/// What `chgm` tries in turn in place of its series `sum`, which gave `hg = e^min(x0, 0) sum.hg`
fn chgm_fallback(
    a: f64,
    a0: f64,
    b: f64,
    x: f64,
    x0: f64,
    hg: f64,
    sum: &KummerSum,
) -> Option<f64> {
    // for -1 < a < 2 and b > 0, the series in full if it doesn't converge within 500 terms
    if !sum.converged && hg.is_finite() && a > -1.0 && a0 < 2.0 && b > 0.0 {
        if let Some((hg, _)) = chgm_series(a, b, x, x0) {
            return Some(hg);
        }
    }
    // for a < -1 and b > 0, first Tricomi's expansion if the series cancels by over 1e3 or
    // doesn't converge; below that the series is usually more accurate
    let neg_a = a < -1.0 && b > 0.0;
    if neg_a && sum.cancellation() > 1e3 {
        if let Some(hg) = chgm_tricomi(a, b, x, x0) {
            return Some(hg);
        }
    }
    // then, if the series is `lossy`, the backward recurrences from Boost.Math for non-integer a
    if sum.lossy(a, b, x) && a != a.trunc() {
        // NOTE: here `a` is in (-1, 2); only the (unused) error bound needs `a >= 1`.
        let m = |a, b| chgm_series(a, b, x, 0.0).map(|(s, _)| s);
        if let Some(hg) = hyp1f1_neg_a(a, b, x, x0.min(0.0), m) {
            return Some(hg);
        }
    }
    // next, if the series cancels by over 1e5, DLMF 13.3.2 backward from b + n >= -a x / 2
    // (n <= 5000) if they differ by over twice its error
    let n = (-a * x / 2.0 - b).ceil();
    let cancels = sum.cancellation() > 1e5 && (1.0..=5e3).contains(&n);
    if neg_a && sum.accurate() && cancels && a != a.trunc() && x > 0.0 {
        let m = |b| chgm_series(a, b, x, x0).map_or(f64::NAN, |(hg, _)| hg);
        let rec = hyp1f1_b_recurrence(a, b, x, n as i32, m);
        if let Some((r, _)) = rec.filter(|&(r, err)| (r - hg).abs() > 2.0 * err * r.abs()) {
            return Some(r);
        }
    }
    None
}

/// Kummer's function `M(a, b, x)`.
///
/// Unlike in xsf, in code order (`x0` is the argument; from 4 on, `a` and `x` are `b - a` and
/// `-x0` if `x0 < 0`, DLMF 13.2.39):
///
/// 1. NaNs, and `a = 0` or `a = b` for finite non-pole `b`, return at once (scipy/xsf#365).
/// 2. Non-integer `b - a < 0`, `30 + |b| < -x0 < 100000`, or `b > 0` and `0 < -x0` if rounding
///    `b - a` errs by more than a rounding: the series of DLMF 13.2.39 in double-double with the
///    exact `b - a`, if its error bound is below `eps` (no issue).
/// 3. Non-integer `b < 0`: `chgm_neg_b`, if that is accurate (scipy/xsf#365).
/// 4. `a = -n >= -10^6`, `b > 0`, `n x > b`, `x0 > -700`: `chgm_poly`, if finite (scipy/xsf#365).
/// 5. `a >= 2`, `b > 0`, or `a > -1`, `x > 30 + |b|`: the series in full (scipy/xsf#365).
/// 6. `a >= 2`, `b = 0`: NaN at once (no issue).
/// 7. No recurrence in `a` over 10^6 steps, nor for the series with non-integer `b < 0`, but that
///    series at `a` itself, or the one at `x0 < 0`, if accurate, or else an overflow (see below), 0
///    or NaN (scipy/xsf#365).
/// 8. The terms don't underflow (see `rescale_term`, scipy/xsf#544), and a negligible nonzero term
///    only marks convergence if they can't grow again (scipy/xsf#454).
/// 9. Non-integer `b < 0`: up to 100000 terms, until the rest is negligible (scipy/xsf#365).
/// 10. An unconverged sum is scaled by `e^x0` (DLMF 13.2.39) too, and past 500 terms with
///     `mul_exp` for `x0 <= -700` (scipy/xsf#365).
/// 11. A sum that cancels or doesn't converge: see `chgm_fallback` (scipy/xsf#365).
/// 12. NaN if a sum is wrong or lossy (see `KummerSum`), for finite `a`, `x0` and non-pole `b`
///     (scipy/xsf#365).
#[allow(clippy::needless_late_init)]
pub(crate) fn chgm(mut x: f64, mut a: f64, b: f64) -> f64 {
    // ===================================================
    // Purpose: Compute confluent hypergeometric function
    //          M(a,b,x)
    // Input  : a  --- Parameter
    //          b  --- Parameter ( b <> 0,-1,-2,... )
    //          x  --- Argument
    // Output:  HG --- M(a,b,x)
    // Routine called: CGAMA for computing complex ln[Г(x)]
    // ===================================================

    let mut a0 = a;
    let a1 = a;
    let x0 = x;
    let mut hg = 0.0;

    // Unlike in xsf (scipy/xsf#365), NaN at once, and M(0, b, x) = 1 and M(b, b, x) = e^x
    // (DLMF 13.6.3, 13.6.1), without the O(a) recurrence below
    if a.is_nan() || b.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    let neg_b = b < 0.0 && b != b.trunc();
    // b is not a pole 0, -1, ...
    let regular = b > 0.0 || neg_b;
    if b.is_finite() && regular && (a == 0.0 || a == b) {
        return if a == 0.0 { 1.0 } else { x.exp() };
    }

    // Unlike in xsf, negative x has its own methods (see `chgm_neg_x`)
    if regular && a.is_finite() {
        if let Some(hg) = chgm_neg_x(x, a, b) {
            return hg;
        }
    }
    // Unlike in xsf (scipy/xsf#365), non-integer b < 0 has its own methods (see `chgm_neg_b`)
    if neg_b && a.is_finite() && x.is_finite() {
        if let Some(hg) = chgm_neg_b(x, a, b) {
            return hg;
        }
    }

    //  DLMF 13.2.39
    if x < 0.0 {
        a = b - a;
        a0 = a;
        x = x.abs();
    }
    // Unlike in xsf (scipy/xsf#365), the polynomial M(-n, b, x) from `chgm_poly` if finite, up to
    // n = 10^6 (the cost bound of `hyp1f1_neg_a`)
    let poly = a < 0.0 && a == a.trunc() && a >= -MAX_SERIES_ITERATIONS && b > 0.0 && -a * x > b;
    if poly && x0 > GRADUAL_EXP_X0 {
        let p = chgm_poly(-a as i32, b, x, x0.min(0.0).exp());
        if p.is_finite() {
            return p;
        }
    }
    // Unlike in xsf (scipy/xsf#365), all of the series DLMF 13.2.2, also for x <= -700; for b < 0
    // only if finite, since the alternating terms (n < -b) can overflow first
    let series = (a >= 2.0 && b > 0.0) || (a > -1.0 && x > 30.0 + b.abs());
    if series && a.is_finite() && b.is_finite() && x.is_finite() {
        if let Some((hg, _)) = chgm_series(a, b, x, x0).filter(|&(hg, _)| b > 0.0 || hg.is_finite())
        {
            return hg;
        }
    }
    // Unlike in xsf, NaN at the pole without the O(a) recurrence below
    if b == 0.0 && a >= 2.0 {
        return f64::NAN;
    }
    // Unlike in xsf (scipy/xsf#365), NaN if a sum below is wrong, except at the poles of b and for
    // infinite a or x
    let check = a1.is_finite() && x0.is_finite() && regular;
    let mut nl = 0;
    let mut la = 0;
    let mut direct = false;
    if a >= 2.0 {
        // preparing terms for DLMF 13.3.1
        nl = 1;
        // NOTE: `(int)a` is UB for a >= 2^31 (see `int_x86`).
        la = int_x86(a);
        // Unlike in xsf (scipy/xsf#365), not if that takes over 10^6 steps, nor for the series
        // below with b < 0, where it is unstable (see `direct` below)
        if f64::from(la) > MAX_SERIES_ITERATIONS || check && neg_b && x <= 30.0 + b.abs() {
            (direct, nl, la) = (true, 0, 0);
        }
        a -= f64::from(la + 1);
    }
    let mut y0 = 0.0;
    let mut y1 = 0.0;
    let mut accurate = false;
    let mut overflow = false;
    let mut wrong = false;
    for n in 0..=nl {
        if a0 >= 2.0 {
            a += 1.0;
        }
        if (x <= 30.0 + b.abs()) || (a < 0.0) {
            // Unlike in xsf (scipy/xsf#365), with more than 500 terms for non-integer b < 0, and
            // NaN in the end if the sum is wrong or lossy
            let sum = chgm_sum(a, b, x, if check && neg_b { MAX_TERMS } else { 500 });
            accurate = sum.accurate();
            wrong |= check && (sum.wrong(a, b, x, neg_b) || sum.lossy(a, b, x));
            // |M| > MAX if the sum overflows where its terms keep their sign (not `wrong`)
            overflow = check && sum.hg.is_infinite() && !wrong;
            hg = sum.hg;
            // DLMF 13.2.39 (cf. above); unlike in xsf (scipy/xsf#365), also if the series doesn't
            // converge within 500 terms, and after 500 terms (finite if not `wrong`) with `mul_exp`
            // for x <= -700, where e^x can be subnormal
            if x0 <= GRADUAL_EXP_X0 && sum.terms > 500 {
                hg = mul_exp(hg, x0);
            } else if x0 < 0.0 {
                hg *= x0.exp();
            }
            // Unlike in xsf (scipy/xsf#365), a sum that cancels or doesn't converge is replaced
            // if another method applies
            if let Some(hg) = chgm_fallback(a, a0, b, x, x0, hg, &sum) {
                return hg;
            }
        } else {
            // DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum
            let cta = cgama(c64(a, 0.0), 0);
            let ctb = cgama(c64(b, 0.0), 0);
            let xg = b - a;
            let ctba = cgama(c64(xg, 0.0), 0);
            let mut sum1 = 1.0;
            let mut sum2 = 1.0;
            let mut r1 = 1.0;
            let mut r2 = 1.0;
            for i in 1..9 {
                r1 = -r1 * (a + f64::from(i) - 1.0) * (a - b + f64::from(i)) / (x * f64::from(i));
                r2 = -r2 * (b - a + f64::from(i) - 1.0) * (a - f64::from(i)) / (x * f64::from(i));
                sum1 += r1;
                sum2 += r2;
            }
            let hg1;
            let hg2;
            if x0 >= 0.0 {
                hg1 = (ctb - ctba).exp().re * x.powf(-a) * (PI * a).cos() * sum1;
                hg2 = (ctb - cta + x).exp().re * x.powf(a - b) * sum2;
            } else {
                // DLMF 13.2.39 (cf. above)
                hg1 = (ctb - ctba + x0).exp().re * x.powf(-a) * (PI * a).cos() * sum1;
                hg2 = (ctb - cta).exp().re * x.powf(a - b) * sum2;
            }
            hg = hg1 + hg2;
        }
        /* 25 */
        if n == 0 {
            y0 = hg;
        }
        if n == 1 {
            y1 = hg;
        }
    }
    // Unlike in xsf (scipy/xsf#365), if `direct`, the series above (at a itself) if accurate, else
    // for x < 0 the one without DLMF 13.2.39, else an overflow
    if direct && !accurate {
        let hg = (x0 < 0.0).then(|| chgm_series(a1, b, x0, 0.0)).flatten();
        if let Some((hg, _)) = hg.filter(half_accurate) {
            return hg;
        }
        // and |e^x0 M| > MAX too if the series with that factor also overflows
        let hg = overflow.then(|| chgm_series(a, b, x, x0)).flatten();
        if let Some((hg, _)) = hg.filter(|(hg, _)| hg.is_infinite()) {
            return hg;
        }
        // DLMF 13.4.1: 0 < M(a, b, x) <= I_s(a, b - a) + e^(s x) for x < 0 < a < b, s in [0, 1]
        let s = a1 / (2.0 * b);
        let zero =
            x0 < 0.0 && 0.0 < a1 && a1 < b && crate::betainc(a1, b - a1, s) + (s * x0).exp() == 0.0;
        return if zero { 0.0 } else { f64::NAN };
    }
    if wrong {
        return f64::NAN;
    }
    if a0 >= 2.0 {
        // DLMF 13.3.1
        for _i in 1..la {
            hg = ((2.0 * a - b + x) * y1 + (b - a) * y0) / a;
            y0 = y1;
            y1 = hg;
            a += 1.0;
        }
    }
    // a = a1;
    // x = x0;
    hg
}

/// Tricomi's function `U(a, b, x)`.
///
/// Unlike in xsf (scipy/xsf#380), which returns the first of `chgus` and `chgul` with 9 digits,
/// and else its step 2:
///
/// 1. `chgus` (non-integer `b`, or 1), then `chgul` (any of `il1` to `il4`, where `il3` and `il4`
///    exclude `b > 1` with `x <= |b - 2a|`). The first with `ID_MIN` digits is returned (`chgus` not
///    if `chgus_deferred`); otherwise the better of the two is kept as `(hu0, id0)`.
/// 2. One method by region:
///    - `a >= 1`: `chgubi` (nonzero integer `b < 1001`, any of `bl1` to `bl3`); `chguit`
///      otherwise, or if `chgubi` has fewer than `ID_MIN` digits and both are finite
///      - an overflowing `chguit` is redone without `x^(1-b)` in the integral, unless step 1 has 9
///        digits
///      - then, if neither this nor step 1 has `ID_ACCEPT` digits, the recurrence in `b` (for
///        `LOG_MAP_B <= b < 1000` and `x < 1`)
///    - `a < 1`, `b <= a`: `chguit` for `U(a-b+1, 2-b, x)`, times `x^(1-b)`; for `-1e5 <= b <= 0`
///      with `x^(1-b)` in the integral if that product isn't normal or step 1 overflowed
///    - `a < 1`, `b > a`: `chgubi` (nonzero integer `b`); the recurrence in `a > -1000` otherwise,
///      or if neither `chgubi` (finite) nor step 1 has `ID_MIN` digits, unless a
///      finite `chgubi` (or `chgus_near_1` in step 1) has at least as many estimated digits as
///      the recurrence (by the cancellation of its two starting values)
/// 3. If `hu0` has 9 to 11 digits, or more than step 2 (and isn't NaN), the result of step 2 only
///    replaces it if the two agree within 10 times its error estimate (step 2 can be far off, e.g.
///    for small `x`). A NaN of step 2 doesn't replace a finite `hu0`. If neither has `ID_ACCEPT`
///    digits, `chgul` is tried if step 1 skipped it for `x <= |b - 2a|` and step 2 had no recurrence
///    in `b`, and kept if it has more.
///
/// `ID_MIN` and the factor 10 are empirical.
pub(crate) fn chgu(
    x: f64,
    mut a: f64,
    mut b: f64,
    md: &mut i32,
    isfer: &mut i32,
    digits: &mut i32,
) -> f64 {
    // =======================================================
    // Purpose: Compute the confluent hypergeometric function
    //          U(a,b,x)
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument  ( x > 0 )
    // Output:  HU --- U(a,b,x)
    //          MD --- Method code
    //          ISFER --- Error flag
    //          DIGITS --- Estimated number of digits (unlike in xsf)
    // Routines called:
    //      (1) CHGUS for small x ( MD=1 )
    //      (2) CHGUL for large x ( MD=2 )
    //      (3) CHGUBI for integer b ( MD=3 )
    //      (4) CHGUIT for numerical integration ( MD=4 )
    // =======================================================

    // Unlike in xsf (scipy/xsf#380), chgus and chgul need 12 digits (instead of 9) to be returned
    // at once (step 1); chgubi needs as many to be kept as is (step 2)
    const ID_MIN: i32 = 12;
    // results with fewer digits are flagged (and NaN in `hypu`)
    const ID_ACCEPT: i32 = 6;

    let mut id1: i32 = 0;
    let mut hu = 0.0;

    let aa = a - b + 1.0;
    *isfer = 0;
    // NOTE: `(int)a` (and `(int)aa` below) is UB in C++ if out of `int` range (incl.
    // inf/NaN; GCC/x86 gives INT_MIN); here it saturates. Either way an out-of-range value never
    // compares equal, so the results of these comparisons are the same.
    let il1 = (a == f64::from(a as i32)) && (a <= 0.0);
    let il2 = (aa == f64::from(aa as i32)) && (aa <= 0.0);
    let bl1 = (x <= 5.0) || (x <= 10.0 && a <= 2.0);
    let bl2 = (x > 5.0) && (x <= 12.5) && ((a >= 1.0) && (b >= a + 4.0));
    let bl3 = (x > 12.5) && (a >= 5.0) && (b >= a + 5.0);
    let bn = b.fract() == 0.0 && b != 0.0;
    // Unlike in xsf (scipy/xsf#380), conservatively not for b > 1 and x <= |b - 2a|: the error
    // bound of the series in chgul (DLMF 13.7.5) is only usable for x > |b - 2a|, and for tiny a
    // the series overestimates its digits (e.g. it is 1 for U(1e-15, 35 + 1e-14, 4) = 1830)
    let unbounded = b > 1.0 && x <= (b - 2.0 * a).abs();
    let il3 = (a * (a - b + 1.0)).abs() / x <= 2.0 && !unbounded;

    let mut id = -100;
    let mut hu1 = 0.0;
    // Unlike in xsf (scipy/xsf#380), also not for the integers |b| >= 2^31, but for b = 1 (by
    // `chgus_near_1`, which also holds there)
    if b != b.trunc() || b == 1.0 {
        hu = chgus(x, a, b, &mut id1);
        *md = 1;
        if id1 >= ID_MIN && !chgus_deferred(a, b, x, id1) {
            *digits = super::subnormal_digits(hu, id1.into()) as i32;
            return hu;
        }
        hu1 = hu;
    }
    // Unlike in xsf (scipy/xsf#380), chgul is also tried for a < 1 and b > a, with b zero or
    // non-integer, before the recurrence in a below (not for |a| < 1.1e-16, where `a + 1 - 1`
    // rounds to 0 in chgul)
    let il4 = a < 1.0 && b > a && !bn && a + 1.0 != 1.0 && !unbounded;
    if il1 || il2 || il3 || il4 {
        hu = chgul(x, a, b, &mut id);
        *md = 2;
        if id >= ID_MIN {
            *digits = super::subnormal_digits(hu, id.into()) as i32;
            return hu;
        }
        if id1 > id {
            *md = 1;
            id = id1;
            hu = hu1;
        }
    }
    // step 1's best
    let (hu0, id0) = (hu, id.max(id1));
    if a >= 1.0 {
        // `chguit` at `(a, b)` and its digits
        let chguit_id = |a: f64, b: f64, kummer: bool| {
            let mut id = 0;
            (chguit(x, a, b, kummer, &mut id), id)
        };
        // Unlike in xsf (scipy/xsf#380), chguit is used for b - 1 >= 1000, where chgubi gives up
        let ubi = bn && (bl1 || bl2 || bl3) && b - 1.0 < 1000.0;
        if ubi {
            hu = chgubi(x, a, b, &mut id);
            *md = 3;
        }
        // Unlike in xsf (scipy/xsf#380), chguit is also used if chgubi has fewer than `ID_MIN`
        // digits, unless chgubi overflows (e.g. for b < -170, where chguit is often wrong too) or
        // chguit does (e.g. U(150, 9, 0.1) = 2.23e-250)
        if !ubi || id < ID_MIN && hu.is_finite() {
            let (mut u, mut id_it) = chguit_id(a, b, false);
            // Unlike in xsf, an overflowing integral is redone without x^(1-b) (DLMF 13.2.40), as
            // Gamma(a) U, or t^(a-1) in it, can overflow if U doesn't, e.g. for
            // U(20, 7 - 2e-11, 6e-18) = 2.1e88. Then x^(1-b) is applied in halves, as it can
            // overflow too. Only a finite result replaces the infinite one (e.g. not inf * 0 if
            // U overflows and x^(1-b) underflows).
            //
            // This is skipped if step 1 has 9 digits, since step 3 then keeps it anyway. After
            // chgubi, it is also done for inf / inf, where Gamma(a) overflows as well, and then
            // only used if it has more digits than chgubi.
            if (u.is_infinite() || ubi && u.is_nan()) && id0 < 9 {
                let h = x.powf(0.5 * (1.0 - b));
                let (v, id_k) = chguit_id(a, b, true);
                let v = v * h * h;
                if v.is_finite() && (u.is_infinite() || id_k > id) {
                    (u, id_it) = (v, id_k);
                }
            }
            if !ubi || u.is_finite() {
                (hu, id) = (u, id_it);
                *md = 4;
            }
        }
        // Unlike in xsf, the forward recurrence in b if this and step 1 fail (e.g. near an integer
        // b), from b - n near 1 and 2, where chguit integrates over ln(1 + t)
        if (hu.is_nan() || id < ID_ACCEPT)
            && id0 < ID_ACCEPT
            && x < 1.0
            && (LOG_MAP_B..1000.0).contains(&b)
        {
            let (u, id_b) = super::hypu_fwd_b(a, b, x, (b - 1.0).round(), |c| {
                let (u, id) = chguit_id(a, c, false);
                (u, f64::from(id))
            });
            let id_b = id_b as i32;
            if !u.is_nan() && (hu.is_nan() || id_b > id) {
                (hu, id) = (u, id_b);
            }
        }
    } else if b <= a {
        let b00 = b;
        a -= b - 1.0;
        b = 2.0 - b;
        *md = 4;
        hu = chguit(x, a, b, false, &mut id);
        hu *= x.powf(1.0 - b00);
        // Unlike in xsf (scipy/xsf#380), if this or the sum in chgul over- or underflows, x^(1-b)
        // is moved into the integral, for -1e5 <= b <= 0 (where its integrand is bounded at 0; in
        // log space it loses about |b| ln|b| eps)
        if (!hu.is_normal() || hu0.is_infinite()) && (-1e5..=0.0).contains(&b00) {
            hu = chguit(x, a, b, true, &mut id);
            if hu0.is_infinite() {
                *digits = super::subnormal_digits(hu, id.into()) as i32;
                return hu;
            }
        }
    } else {
        // NOTE: xsf checks `bn && ~il1`, where `~il1` (-1 or -2) is always true
        if bn {
            hu = chgubi(x, a, b, &mut id);
            *md = 3;
        }
        // Unlike in xsf (scipy/xsf#380), the recurrence in a is used if no other method applies,
        // or if chgubi and chgul have fewer than `ID_MIN` digits (or chgubi overflows); a finite
        // chgubi is only replaced if the recurrence has more estimated digits, and so is
        // `chgus_near_1`, as its estimate is conservative
        let failed = id < ID_MIN || !hu.is_finite();
        let near_1 = (b - 1.0).abs() < NEAR_1_DD && id0 == id1 && hu0.is_finite();
        if (!bn || failed && id0 < ID_MIN) && a > -1000.0 {
            let (u, id_rec) = super::hypu_rec(a, b, x);
            let id_rec = id_rec as i32;
            if !u.is_nan() && (!bn || !hu.is_finite() || id_rec > id) && (!near_1 || id_rec > id0) {
                // (at most `ID_MIN`, as the estimate of the recurrence in a is not a bound)
                (hu, id) = (u, id_rec.min(ID_MIN));
            }
        }
    }
    // NOTE: if none of the branches is taken (e.g. b = 0 and a < 0), `*md` is only set if
    // chgus/chgul was called (`hypu` passes an uninitialized `md`, but never reads it).

    // Unlike in xsf, a NaN here doesn't replace a finite result of step 1, e.g. of chguit for
    // subnormal x, where 1/x overflows
    if hu.is_nan() && hu0.is_finite() {
        (hu, id) = (hu0, id0);
    }
    // Unlike in xsf, chgul is the last resort where it was skipped as `unbounded`, e.g. for
    // U(172, 173 + 3e-14, 50), where chgus and chguit overflow, and U(171, 172 + 3e-14, 50), where
    // chgus has too few digits. Not where step 2 has the recurrence in b (a < 1 and 3 <= b < 1000),
    // as where it fails (near a negative integer a) the smallest term underestimates the error,
    // e.g. it gives -9.8e7 for U(-3 - 5.6e-12, 904 + 1.2e-9, 442.7) = 1.8e69
    let rec_b = a < 1.0 && (3.0..1000.0).contains(&b);
    if unbounded && !rec_b && (hu.is_nan() || id < ID_ACCEPT) {
        let mut id_l = 0;
        let u = chgul(x, a, b, &mut id_l);
        if hu.is_nan() || id_l > id {
            (hu, id) = (u, id_l);
            *md = 2;
        }
    }
    // Step 3. Unlike in xsf, also if step 1 has more digits than step 2, e.g. if chguit doesn't
    // converge
    if (id0 >= 9 || id0 > id) && !hu0.is_nan() {
        let tol = 10f64.powi(1 - id0) * hu0.abs();
        let agrees = hu == hu0 || hu0.is_finite() && (hu - hu0).abs() <= tol;
        if !agrees {
            (hu, id) = (hu0, id0);
        }
        id = id.max(id0);
    }
    if id < ID_ACCEPT {
        *isfer = 6;
    }
    // Unlike in xsf, a subnormal U has fewer digits (also at the returns above), which doesn't
    // flag it (as no method does better), and isn't compared above (where it would favor other
    // methods that cancel to a normal U, e.g. `chgus` near U(167.7, -4, 1.22) = 4.9e-318)
    *digits = super::subnormal_digits(hu, id.into()) as i32;
    hu
}

fn chgubi(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
    // Unlike in xsf (scipy/xsf#380), for a < 1 with `digamma` and `gamma` instead of `psi_spec` and
    // `gamma2`, which lose digits near the poles (e.g. both are 8% off at -2 + 1e-15), unless Gamma
    // isn't normal or the result isn't finite, e.g. U(-169 - 3e-14, 14, 4) = inf (not inf - inf)
    if a < 1.0 {
        let hu = chgubi_with(x, a, b, id, false);
        if hu.is_finite() {
            return hu;
        }
    }
    chgubi_with(x, a, b, id, true)
}

/// `chgubi` with `psi_spec` and `gamma2` if `xsf`, else with `digamma` and `gamma`
fn chgubi_with(x: f64, a: f64, b: f64, id: &mut i32, xsf: bool) -> f64 {
    // ======================================================
    // Purpose: Compute confluent hypergeometric function
    //          U(a,b,x) with integer b ( b = ±1,±2,... )
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument
    // Output:  HU --- U(a,b,x)
    //          ID --- Estimated number of significant digits
    // Routines called:
    //      (1) GAMMA2 for computing gamma function Г(x)
    //      (2) PSI_SPEC for computing psi function
    // ======================================================

    const EL: f64 = EGAMMA;

    /// `y f(1) f(2) ... f(m)`, with `y` and the result as `(p, e)` for `p 2^e`
    fn prod(y: (f64, i32), m: i32, f: impl Fn(f64) -> f64) -> (f64, i32) {
        (1..=m).fold(y, |(p, e), j| {
            let (p, k) = frexp_f64(p * f(f64::from(j)));
            (p, e + k)
        })
    }

    // Unlike in xsf (scipy/xsf#380), NaN for |b - 1| >= 1000, without ~150 n steps
    if (b - 1.0).abs() >= 1000.0 {
        return f64::NAN;
    }
    *id = -100;
    let n = (b - 1.0).abs() as i32;
    let ps = if xsf { psi_spec(a) } else { digamma(a) };
    // Unlike in xsf, the prefactors `ua` and `ub` are products of ratios with their exponents
    // apart (`p 2^e`), as n!, x^n and Gamma(a) over- or underflow where U(a, b, x) doesn't (e.g.
    // for U(0.936, 172, 197)), and `a -+ n` rounds.
    //
    // 1 / Gamma(a) = rg, by Gamma(a) = Gamma(a0) (a - 1) ... (a - m) for a > 171, where it
    // overflows, and Gamma(a0) = Gamma(a) (a0 - 1) ... (a0 + m) for a < -170, where it underflows.
    // Beyond |m| = 1e5, rg is 0 or infinite, like U(a, b, x).
    let m = (a - 171.0).clamp(0.0, 1e5).ceil() + (a + 170.0).clamp(-1e5, 0.0).floor();
    let a0 = a - m;
    let ga = if xsf { gamma2(a0) } else { gamma(a0) };
    // `chgubi` then retries with `psi_spec` and `gamma2`
    if !xsf && !ga.is_normal() {
        return f64::NAN;
    }
    let (pg, eg) = frexp_f64(ga);
    let shift = |j| if m > 0.0 { 1.0 / (a - j) } else { a0 - j };
    let rg = prod((1.0 / pg, -eg), m.abs() as i32, shift);
    // x^n = pn 2^(n ex)
    let (px, ex) = frexp_f64(x);
    let pn = px.powf(f64::from(n));
    let (n0, n2, ua, mut ub);
    if b > 0.0 {
        (n0, n2) = (0, -n);
        // (-1)^(n+1) / (n! Gamma(a - n)), by `ln_1p` for |a| < 1, as j - a rounds off a tiny a
        ua = if a.abs() < 1.0 {
            let l: f64 = (1..=n).map(|j| (-a / f64::from(j)).ln_1p()).sum();
            (-rg.0 * l.exp(), rg.1)
        } else {
            prod((-rg.0, rg.1), n, |j| (j - a) / j)
        };
        // (n - 1)! / (x^n Gamma(a))
        ub = prod((rg.0 / pn, rg.1 - n * ex), n - 1, |j| j);
    } else {
        (n0, n2) = (n, 0);
        // (-1)^(n-1) x^n / (n! Gamma(a)) and (n - 1)! / Gamma(a + n)
        ua = prod((-rg.0 * pn, rg.1 + n * ex), n, |j| -1.0 / j);
        ub = prod((rg.0 / a, rg.1), n - 1, |j| j / (a + j));
    }
    // Unlike in xsf, the factors of the series are a + (n0 + k - 1) (xsf has a0 + k - 1 with
    // a0 = a + n0, and likewise for a2 = a + n2), which round once, as a +- n loses the low digits
    // of a for large n, e.g. for U(0.936, 150, 170)
    let mut hm1 = 1.0;
    let mut r = 1.0;
    let mut hmax: f64 = 0.0;
    let mut h0 = 0.0;
    for k in 1..=150 {
        r = r * (a + f64::from(n0 + k - 1)) * x / f64::from((n + k) * k);
        hm1 += r;
        let hu1 = hm1.abs();

        if hu1 > hmax {
            hmax = hu1;
        }
        if (hm1 - h0).abs() < hm1.abs() * 1.0e-15 {
            break;
        }
        h0 = hm1;
    }

    let (hmax1, t1) = (hmax, r);
    hm1 *= x.ln();
    let mut s0 = 0.0;

    for m in 1..=n {
        if b >= 0.0 {
            s0 -= 1.0 / f64::from(m);
        }
        if b < 0.0 {
            s0 += (1.0 - a) / (f64::from(m) * (a + f64::from(m) - 1.0));
        }
    }

    let mut hm2 = ps + 2.0 * EL + s0;
    let mut t2 = 0.0;
    r = 1.0;
    hmax = 0.0;

    for k in 1..=150_i32 {
        let mut s1 = 0.0;
        let mut s2 = 0.0;

        if b > 0.0 {
            for m in 1..=k {
                s1 -= (f64::from(m) + 2.0 * a - 2.0) / (f64::from(m) * (f64::from(m) + a - 1.0));
            }
            for m in 1..=n {
                s2 += 1.0 / f64::from(k + m);
            }
        } else {
            for m in 1..=k + n {
                s1 += (1.0 - a) / (f64::from(m) * (f64::from(m) + a - 1.0));
            }
            for m in 1..=k {
                s2 += 1.0 / f64::from(m);
            }
        }

        let hw = 2.0 * EL + ps + s1 - s2;
        r = r * (a + f64::from(n0 + k - 1)) * x / f64::from((n + k) * k);
        t2 = r * hw;
        hm2 += t2;
        let hu2 = hm2.abs();

        if hu2 > hmax {
            hmax = hu2;
        }

        // NOTE: `h0` is not reset: in the first iteration it still holds the last `hm1` of the
        // first series (before `hm1 *= log(x)`).
        if ((hm2 - h0) / hm2).abs() < 1.0e-15 {
            break;
        }

        h0 = hm2;
    }

    let mut hm3 = 1.0;
    if n == 0 {
        hm3 = 0.0;
    }
    r = 1.0;
    for k in 1..n {
        r = r * (a + f64::from(n2 + k - 1)) / f64::from((k - n) * k) * x;
        hm3 += r;
        // Unlike in xsf, scaled down (exactly) above SCALE, as hm3 overflows where ub underflows
        if hm3.abs() > SCALE {
            (r, hm3, ub.1) = (r / SCALE, hm3 / SCALE, ub.1 + LOG2_SCALE);
        }
    }
    // Unlike in xsf, the terms are added with the larger exponent of the nonzero ones split off
    let (sa, sb) = (ua.0 * (hm1 + hm2), ub.0 * hm3);
    let by_a = sb == 0.0 || sa != 0.0 && ua.1 > ub.1;
    let e = if by_a { ua.1 } else { ub.1 };
    let [ua, sa, sb] = [ua, (sa, ua.1), (sb, ub.1)].map(|(p, k)| ldexp_f64(p, k - e));
    let hu = ldexp_f64(sa + sb, e);
    // (as in xsf, inf - inf is NaN if both overflow)
    let t = ldexp_f64(sa, e) + ldexp_f64(sb, e);
    let hu = if hu.is_infinite() { t } else { hu };
    // Unlike in xsf (scipy/xsf#380), the digits are estimated from the rounding errors (of the
    // largest partial sums) and the truncation errors (the last terms) of the series. If `hu`
    // overflows (inf / inf) or underflows to 0 (0 / 0), `min` turns the NaN into 15 digits.
    //
    // Unlike in xsf, if n! overflows or Gamma(a) is shifted, an infinite `hu` has the digits of the
    // sum without its exponent, as its sign can be a rounding error (e.g. for U(-160.5, 400, 300)).
    let unscaled = hu.is_infinite() && (m != 0.0 || f64::from(n) > MAXGAM - 1.0);
    let k = if unscaled { 0 } else { e };
    let lnx = x.ln().abs();
    let round = 1e-15 * (ua.abs() * (hmax1 * lnx + hmax) + sb.abs());
    let trunc = ua.abs() * (t1.abs() * lnx + t2.abs());
    *id = (-(ldexp_f64(round + trunc, k) / ldexp_f64(sa + sb, k).abs()).log10()).min(15.0) as i32;
    hu
}

/// `chguit` integrates over `ln(1 + t)` for `b` below this (and `x < 1`), where the integrand at
/// `t` of order 1 isn't negligible (empirical; above it, it can overflow for tiny `x`)
const LOG_MAP_B: f64 = 2.5;

/// Whether `chgu` goes on to its other methods at `(a, b, x)` even if `chgus` claims `ID_MIN`
/// digits (`id`), for `b` near 1
fn chgus_deferred(a: f64, b: f64, x: f64, id: i32) -> bool {
    (a >= 1.0 || x >= 1.0) && (b - 1.0).abs() < NEAR_1_DD && id < NEAR_1_ID_DEFER
}

#[allow(clippy::let_and_return)]
fn chguit(x: f64, a: f64, b: f64, kummer: bool, id: &mut i32) -> f64 {
    // ======================================================
    // Purpose: Compute hypergeometric function U(a,b,x) by
    //          using Gaussian-Legendre integration (n=60)
    // Input  : a  --- Parameter ( a > 0 )
    //          b  --- Parameter
    //          x  --- Argument ( x > 0 )
    // Output:  HU --- U(a,b,z)
    //          ID --- Estimated number of significant digits
    // Routine called: GAMMA2 for computing Г(x)
    // ======================================================

    const T: [f64; 30] = [
        0.259597723012478e-01,
        0.778093339495366e-01,
        0.129449135396945e+00,
        0.180739964873425e+00,
        0.231543551376029e+00,
        0.281722937423262e+00,
        0.331142848268448e+00,
        0.379670056576798e+00,
        0.427173741583078e+00,
        0.473525841761707e+00,
        0.518601400058570e+00,
        0.562278900753945e+00,
        0.604440597048510e+00,
        0.644972828489477e+00,
        0.683766327381356e+00,
        0.720716513355730e+00,
        0.755723775306586e+00,
        0.788693739932264e+00,
        0.819537526162146e+00,
        0.848171984785930e+00,
        0.874519922646898e+00,
        0.898510310810046e+00,
        0.920078476177628e+00,
        0.939166276116423e+00,
        0.955722255839996e+00,
        0.969701788765053e+00,
        0.981067201752598e+00,
        0.989787895222222e+00,
        0.995840525118838e+00,
        0.999210123227436e+00,
    ];
    const W: [f64; 30] = [
        0.519078776312206e-01,
        0.517679431749102e-01,
        0.514884515009810e-01,
        0.510701560698557e-01,
        0.505141845325094e-01,
        0.498220356905502e-01,
        0.489955754557568e-01,
        0.480370318199712e-01,
        0.469489888489122e-01,
        0.457343797161145e-01,
        0.443964787957872e-01,
        0.429388928359356e-01,
        0.413655512355848e-01,
        0.396806954523808e-01,
        0.378888675692434e-01,
        0.359948980510845e-01,
        0.340038927249464e-01,
        0.319212190192963e-01,
        0.297524915007890e-01,
        0.275035567499248e-01,
        0.251804776215213e-01,
        0.227895169439978e-01,
        0.203371207294572e-01,
        0.178299010142074e-01,
        0.152746185967848e-01,
        0.126781664768159e-01,
        0.100475571822880e-01,
        0.738993116334531e-02,
        0.471272992695363e-02,
        0.202681196887362e-02,
    ];
    /// The composite 60-point Gauss-Legendre rule for `f` over `m` panels of width `2 g` from 0
    fn panels(m: i32, g: f64, f: impl Fn(f64) -> f64) -> f64 {
        let (mut sum, mut d) = (0.0, g);
        for _j in 1..=m {
            let mut s = 0.0;
            for k in 1..=30 {
                let t1 = d + g * T[k - 1];
                let t2 = d - g * T[k - 1];
                s += W[k - 1] * (f(t1) + f(t2));
            }
            sum += s * g;
            d += 2.0 * g;
        }
        sum
    }
    // DLMF 13.4.4, integration up to C (12/X, or as below)
    let a1 = a - 1.0;
    let b1 = b - a - 1.0;
    // Unlike in xsf (scipy/xsf#380), with `kummer` this is x^(b-1) U(a, b, x) = U(a-b+1, 2-b, x)
    // (DLMF 13.2.40), integrated over s = x t in log space, which can be finite where U or
    // x^(b-1) over- or underflows. The integrand peaks at s < b, so up to C = 12 + b.
    //
    // Unlike in xsf, for b < LOG_MAP_B (e.g. near 0, 1 or 2, where chgus cancels) and x < 1, up
    // to C = 1/x and over ln(1 + t): the integrand is then ~ t^(b-2) from t ~ 1 up to t ~ 1/x,
    // which uniform panels in t can't resolve for small x
    let log_map = !kummer && b < LOG_MAP_B && x < 1.0;
    let c = if kummer {
        12.0 + b
    } else if log_map {
        1.0 / x
    } else {
        12.0 / x
    };
    let lg = if kummer { gammaln(a) } else { 0.0 };
    let f = |t: f64| {
        if kummer {
            // s^(a-1) (s + x)^(b-a-1) e^(-s) / Gamma(a)
            (a1 * t.ln() + b1 * (t + x).ln() - t - lg).exp()
        } else if log_map {
            // without overflow of t^(a-1)
            (-x * t).exp() * (t / (1.0 + t)).powf(a1) * (1.0 + t).powf(b - 2.0)
        } else {
            (-x * t).exp() * t.powf(a1) * (1.0 + t).powf(b1)
        }
    };
    // Unlike in xsf (scipy/xsf#380), t = c u^3 for non-integer a < 4, which smooths the endpoint
    // singularity of t^(a-1), where the quadrature converges slowly (or not at all for a near 1)
    let int_a = a == a.trunc();
    let cubic = a < 4.0 && !int_a;
    let lc = c.ln_1p();
    // the integrand, times dt/du = 3 c u^2 if `cubic`; if `log_map`, t = (1 + c)^v - 1 with
    // v = u^3 for any non-integer a (as for `cubic`), else v = u
    let f_head = |u: f64| {
        if log_map {
            let (v, dv) = if int_a {
                (u, 1.0)
            } else {
                (u * u * u, 3.0 * u * u)
            };
            let t = (lc * v).exp_m1();
            dv * lc * (1.0 + t) * f(t)
        } else if cubic {
            let v = c * u * u;
            3.0 * v * f(v * u)
        } else {
            f(u)
        }
    };
    let mut hu0 = 0.0;
    let mut hu1 = 0.0;
    // Unlike in xsf, a refinement m (by dm) of an integral is only final if its last two relative
    // changes (if there are two) are small, as one can be small by chance where the error
    // oscillates. Their max times m / dm, about the error if it is ~ m^-p with p >= 1, gives the
    // id (in xsf it is 9).
    let (mut d_head, mut d_tail, mut d0) = (0.0, 0.0, 0.0);
    for m in (10..=100).step_by(5) {
        let g = if cubic || log_map { 0.5 } else { 0.5 * c } / f64::from(m);
        hu1 = panels(m, g, f_head);
        let d = (1.0 - hu0 / hu1).abs();
        d_head = d.max(d0) * f64::from(m / 5);
        if d < 1.0e-9 && d0 < 1.0e-9 {
            break;
        }
        // (the first "change" is from 0)
        d0 = if m == 10 { 0.0 } else { d };
        hu0 = hu1;
    }
    let ga = if kummer { 1.0 } else { gamma2(a) };
    hu1 /= ga;
    // DLMF 13.4.4 with substitution t=C/(1-u)
    // integration u from 0 to 1, i.e. t from C to infinity
    let mut hu2 = 0.0;
    // Unlike in xsf (scipy/xsf#380), a non-finite integrand (0 * inf) is recomputed in log space,
    // unless b is infinite or Gamma(a) overflows (a > 171)
    let log_space = b1.is_finite() && a < 171.0;
    // the integrand times dt/du = t^2 / c
    let f_tail = |t: f64| {
        // t^2 overflows for x < ~1e-150
        if log_map {
            return t / c * (t * f(t));
        }
        let dt = t * t / c;
        if kummer {
            return dt * f(t);
        }
        let ft = dt * (-x * t).exp() * t.powf(a1) * (1.0 + t).powf(b1);
        if log_space && !ft.is_finite() {
            dt * (a1 * t.ln() + b1 * t.ln_1p() - x * t).exp()
        } else {
            ft
        }
    };
    d0 = 0.0;
    for m in (2..=10).step_by(2) {
        hu2 = panels(m, 0.5 / f64::from(m), |u| f_tail(c / (1.0 - u)));
        // NOTE: in the first iteration `hu0` still holds a `hu1` of the first integral (the
        // second-to-last one if that loop converged, else the last one), not divided by `ga`.
        let d = (1.0 - hu0 / hu2).abs();
        d_tail = d.max(d0) * f64::from(m / 2);
        if d < 1.0e-9 && d0 < 1.0e-9 {
            break;
        }
        d0 = if m == 2 { 0.0 } else { d };
        hu0 = hu2;
    }
    hu2 /= ga;
    let hu = hu1 + hu2;
    // An integral that is 0 (e.g. if its integrand underflows) adds no error. If U is 0 or
    // infinite (e.g. if it under- or overflows) it has 9 digits as in xsf, else none if NaN.
    // (relative to U, as d |h| can underflow)
    let err = |d: f64, h: f64| if h == 0.0 { 0.0 } else { d * (h / hu).abs() };
    let rel = err(d_head, hu1) + err(d_tail, hu2);
    // the rounding error, which grows with the exponents in the integrand and ln Gamma(a)
    // (empirical)
    let rounding = 2.0 * f64::EPSILON * (1.0 + a.abs() + b.abs() + lg.abs());
    *id = if hu == 0.0 || hu.is_infinite() {
        9
    } else if rel < 1.0 {
        (-rel.max(rounding).log10()) as i32
    } else {
        0
    };
    hu
}

fn chgul(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
    // =======================================================
    // Purpose: Compute the confluent hypergeometric function
    //          U(a,b,x) for large argument x
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument
    // Output:  HU --- U(a,b,x)
    //          ID --- Estimated number of significant digits
    // =======================================================

    let mut hu;
    let mut r;
    let mut r0 = 0.0;
    let mut ra: f64 = 0.0;

    *id = -100;
    let aa = a - b + 1.0;
    // Unlike in xsf (scipy/xsf#380), a (or aa) = -2^31 is excluded
    let is_nonpos_i32 = |z: f64| z == z.trunc() && f64::from(i32::MIN) < z && z <= 0.0;
    let il1 = is_nonpos_i32(a);
    let il2 = is_nonpos_i32(aa);
    let mut nm = 0;
    if il1 {
        nm = -a as i32;
    }
    if il2 {
        nm = -aa as i32;
    }
    // IL1: DLMF 13.2.7 with k=-s-a
    // IL2: DLMF 13.2.8
    if il1 || il2 {
        hu = 1.0;
        r = 1.0;
        let mut rmax: f64 = 1.0;
        // Unlike in xsf (scipy/xsf#380), the sum is hu 2^e_sum, where hu is scaled down (exactly)
        // before it overflows, so that an infinite U has its sign
        let mut e_sum = 0_i32;
        // (a + k - 1) (a - b + k) < 0 for k between 1 - a and b - a, so the terms keep their sign
        // there, and a sum that is +-inf at k >= k_inf stays +-inf
        let (lo, hi) = ((1.0 - a).min(b - a), (1.0 - a).max(b - a));
        let k_inf = if f64::from(nm) < hi {
            lo.floor()
        } else {
            f64::INFINITY
        };
        // whether the sum stopped early, as a lower bound of |sum|
        let mut partial = false;
        for k in 1..=nm {
            let k = f64::from(k);
            let prev = r;
            r = -r * (a + (k - 1.0)) * (a - b + k) / (k * x);
            hu += r;
            rmax = rmax.max(r.abs());
            if rmax > SCALE {
                (r, hu, rmax) = (r / SCALE, hu / SCALE, rmax / SCALE);
                e_sum = e_sum.saturating_add(LOG2_SCALE);
            }
            // Unlike in xsf (scipy/xsf#380), the loop stops once
            // - the sum is NaN, or +-inf for good (k >= k_inf)
            // - the terms decrease below `MIN_POSITIVE`
            // - the sum overflows (unscaled) and has more than `MAX_TERMS` terms: it is then
            //   `partial` if it has the sign of the terms that are left (k >= k_inf), and else NaN
            if e_sum > LOG2_SCALE && k > f64::from(MAX_TERMS) && k < f64::from(nm) {
                partial = k >= k_inf && hu * r > 0.0;
                hu = if partial { hu } else { f64::NAN };
            }
            if hu.is_nan()
                || partial
                || hu.is_infinite() && k >= k_inf
                || r.abs() < f64::MIN_POSITIVE && r.abs() <= prev.abs()
            {
                break;
            }
        }
        *id = 10;
        // Unlike in xsf (scipy/xsf#380), if the sum cancels or overflows, U(-n, c, x) with
        // c + x >= 0 is computed by the recurrence in a (DLMF 13.3.7, found to be stable then), and
        // otherwise the id reflects the cancellation, so that other methods are tried
        let cancels = rmax > 1e6 * hu.abs();
        // (unscaled, i.e. if its largest term does)
        let overflows = e_sum > LOG2_SCALE || !hu.is_finite();
        // U(a, b, x) = x^(1-b) U(a-b+1, 2-b, x) (DLMF 13.2.40)
        let kummer = !il1 || (il2 && b < 1.0);
        let (n, c) = if kummer { (aa, 2.0 - b) } else { (a, b) };
        // U(n, c, x) = u 2^e, where u is scaled down (exactly) before it overflows, so that an
        // infinite U has its sign
        let (mut u, mut e) = (f64::NAN, 0_i32);
        if (cancels || overflows) && c + x >= 0.0 {
            // from U(0, c, x) = 1 (U(1, c, x) has coefficient 0)
            let mut u1 = 0.0;
            u = 1.0;
            for k in 0..-n as i32 {
                // U(k - 1, c, x) = (2k - c + x) U(k, c, x) - k (k - c + 1) U(k + 1, c, x)
                let k = -f64::from(k);
                (u, u1) = ((2.0 * k - c + x) * u - k * (k - c + 1.0) * u1, u);
                if u.abs() > SCALE {
                    (u, u1, e) = (u / SCALE, u1 / SCALE, e.saturating_add(LOG2_SCALE));
                }
            }
        }
        // m 2^e x^p, for finite m and x in log space if x^p, or its product with m, over- or
        // underflows
        let pow_scaled = |m: f64, e: i32, p: f64| {
            let xp = x.powf(p);
            if xp.is_normal() && (m * xp).is_normal() || !(m.is_finite() && x.is_finite()) {
                ldexp_f64(m * xp, e)
            } else {
                let l = ln_dd(DoubleDouble::new(x)) * p;
                Scaled::new(m, e, 0.0).map_or(0.0, |y| y.with_l(l).to_f64().0)
            }
        };
        if u.is_finite() {
            hu = pow_scaled(u, e, if kummer { 1.0 - b } else { 0.0 });
            // the recurrence is stable here, so this has (nearly) all digits
            *id = 15;
        } else {
            let sum = hu;
            hu = pow_scaled(sum, e_sum, -a);
            // a partial sum only shows that U overflows
            if partial && hu.is_finite() {
                hu = f64::NAN;
            } else if sum.is_finite() {
                *id = (15.0 - (rmax / sum.abs()).log10()) as i32;
            }
        }
    } else {
        // DLMF 13.7.3
        hu = 1.0;
        r = 1.0;
        let mut rmax: f64 = 1.0;
        for k in 1..=25 {
            r = -r * (a + f64::from(k) - 1.0) * (a - b + f64::from(k)) / (f64::from(k) * x);
            ra = r.abs();
            if ((k > 5) && (ra >= r0)) || (ra < 1e-15) {
                break;
            }
            r0 = ra;
            hu += r;
            rmax = rmax.max(ra);
        }
        // Unlike in xsf (scipy/xsf#380), the error estimate is relative to min(|sum|, 1), and
        // includes the rounding of the largest term, so a last term of 0 needs no special case
        // (unlike the one of scipy/xsf#430)
        //
        // NOTE: the `int` conversion is UB in C++ if the estimate is not finite; here it
        // saturates (0 for NaN).
        *id = (-(ra.max(rmax * 1e-15) / hu.abs().min(1.0)).log10()) as i32;
        hu *= x.powf(-a);
    }
    hu
}

fn chgus(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
    // ======================================================
    // Purpose: Compute confluent hypergeometric function
    //          U(a,b,x) for small argument x
    // Input  : a  --- Parameter
    //          b  --- Parameter ( b <> 0,-1,-2,...)
    //          x  --- Argument
    // Output:  HU --- U(a,b,x)
    //          ID --- Estimated number of significant digits
    // Routine called: GAMMA2 for computing gamma function
    // ======================================================

    // DLMF 13.2.42 with prefactors rewritten according to
    // DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2

    if (b - 1.0).abs() < NEAR_1_DD {
        // Unlike in xsf, in divided-difference form near b = 1, where the two series cancel
        if let Some(hu) = chgus_near_1(x, a, b - 1.0, id) {
            return hu;
        } else if b == 1.0 {
            // DLMF 13.2.42 has no limit form here
            *id = 0;
            return f64::NAN;
        }
    }
    *id = 100;
    let ga = gamma2(a);
    let gab = gamma2(1.0 + a - b);
    // Unlike in xsf, with sinpi, as sin(pi b) loses about eps / |b - n| near an integer n, and
    // reflected (DLMF 5.5.3) where gamma2(b) (b < -1) or gamma2(2 - b) (b > 3) loses as much
    let hu0 = PI / cephes_sinpi(b);
    let mut r1 = if b < -1.0 {
        gamma2(1.0 - b) / gab
    } else {
        hu0 / (gab * gamma2(b))
    };
    // r2 without x^(1-b), which it is multiplied by last (as Gamma(b - 1) x^(1-b) can overflow)
    let c = if b > 3.0 {
        -gamma2(b - 1.0) / ga
    } else {
        hu0 / (ga * gamma2(2.0 - b))
    };
    let xp = x.powf(1.0 - b);
    let mut r2 = c * xp;
    // Unlike in xsf (scipy/xsf#380), 1/Gamma is 0 at the poles (where gamma2 returns 1e300)
    if ga == SPECFUN_INF {
        r2 = 0.0;
    }
    if gab == SPECFUN_INF {
        r1 = 0.0;
    }
    // Unlike in xsf, no digits if r1 or r2 is 0 but not at a pole (or as x^(1-b) underflows, see
    // below), e.g. for U(115.7, -62.2, 4.86) = 5.3e-243, where r1 = Gamma(63.2) / inf and
    // hu = 1.5e-226
    let lost = r1 == 0.0 && gab != SPECFUN_INF || c == 0.0 && ga != SPECFUN_INF && xp != 0.0;
    // Unlike in xsf (scipy/xsf#380), if r2 overflows, or x^(1-b) or r2 underflows for b > 1, the
    // terms are scaled by 1/|r2| = exp(-lr2) (and hu back), unless r1 then overflows. Its terms
    // grow by up to about e^x, so a tiny r2 isn't negligible, e.g. for U(0.94, 165.43, 139.13) =
    // 1.83, nor is the rounding error of a subnormal x^(1-b), e.g. 1e-9 for U(100.5, 103.5, 1200).
    // For b < 1, x^(1-b) only underflows for tiny x, where they don't grow. lr2 is a
    // double-double, as an ulp of |lr2| > 700 is 1e-13 of exp(lr2).
    let mut lr2 = None;
    // m exp(l), and 0 for m = 0 (e.g. r1 at a pole, where r1 exp(-lr2) can be 0 inf)
    let mul_exp_dd = |m: f64, l| Scaled::new(m, 0, 0.0).map_or(0.0, |y| y.with_l(l).to_f64().0);
    let tiny = b > 1.0 && !(xp.is_normal() && r2.is_normal()) && c != 0.0 && ga != SPECFUN_INF;
    if (r2.is_infinite() || tiny) && c.is_finite() {
        let l = ln_dd(DoubleDouble::new(c.abs())) + ln_dd(DoubleDouble::new(x)) * (1.0 - b);
        let r1s = mul_exp_dd(r1, -l);
        if r1s.is_finite() {
            (r1, r2, lr2) = (r1s, c.signum(), Some(l));
        }
    }
    let (mut s1, mut s2) = (r1, r2);
    let mut hu = r1 - r2;
    let mut hmax: f64 = 0.0;
    let mut hmin: f64 = SPECFUN_INF;
    // Unlike in xsf, the id also has the rounding errors of the terms, about j eps for term j
    let mut terms = r1.abs() + r2.abs();
    // Unlike in xsf (scipy/xsf#380), the series only converge once the terms of both can't grow
    // again (see `term_ratio`), e.g. the r2 terms of U(123.77, 128.58, 42.47) after j = b - 1, and
    // the id of unconverged series is capped by the remaining terms. These are kept from
    // underflowing before that (see `rescale_term`), e.g. the r2 terms of U(164.9, 168.90001, 55)
    // from j = 64.
    let mut converged = false;
    // the terms r1 and r2 as scaled k1 and k2 times
    let (mut t1, mut t2, mut k1, mut k2) = (r1, r2, 0, 0);
    let r2_next = |r2: f64, j: f64| r2 * (a - b + j) / (j * (1.0 - b + j)) * x;
    for j in 1..151 {
        let n = f64::from(j);
        (t1, r1, _) = rescale_term(kummer_term(t1, a, b, j) * x, f64::abs, &mut k1);
        (t2, r2, _) = rescale_term(r2_next(t2, n), f64::abs, &mut k2);
        hu += r1 - r2;
        s1 += r1;
        s2 += r2;
        terms += n * (r1.abs() + r2.abs());
        let hua = hu.abs();
        if hua > hmax {
            hmax = hua;
        }
        if hua < hmin {
            hmin = hua;
        }
        // Unlike in xsf, the terms must be small, not their difference (which is also small if they
        // cancel, e.g. at j = 2 for U(10, 1e-17, 0.5))
        if r1.abs() + r2.abs() < hu.abs() * 1e-15
            && term_ratio(a, b, x, n).is_some()
            && term_ratio(a - b + 1.0, 2.0 - b, x, n).is_some()
        {
            converged = true;
            break;
        }
    }
    let d1 = hmax.log10();
    let mut d2 = 0.0;
    if hmin != 0.0 {
        d2 = hmin.log10();
    }
    // NOTE: `15 - (...)` is evaluated in double and then truncated to int. It is -inf (UB in C++,
    // GCC/x86: INT_MIN, same as the saturated i32::MIN here) if `hmax` is 0 or inf.
    *id = (15.0 - (if d1 - d2 < 0.0 { d2 - d1 } else { d1 - d2 })) as i32;
    if !converged {
        // the last term, and the sum of the next 150 (which can be about sqrt(x) times their
        // largest, e.g. for r2 terms that peak at j ~ x)
        let mut rest = (r1 - r2).abs();
        for j in 151..=300 {
            (t1, r1, _) = rescale_term(kummer_term(t1, a, b, j) * x, f64::abs, &mut k1);
            (t2, r2, _) = rescale_term(r2_next(t2, f64::from(j)), f64::abs, &mut k2);
            rest += (r1 - r2).abs();
        }
        *id = (*id).min((-(rest / hu.abs()).log10()) as i32);
    }
    // Unlike in xsf (scipy/xsf#380), the id includes the errors of the sums s1 and s2, incl. those
    // of gamma2 and hu0 near their poles (z < 0), amplified if s1 and s2 cancel
    let k = |z: f64| {
        if z < 0.0 && z.fract() != 0.0 {
            f64::EPSILON * (z / cephes_sinpi(z)).abs()
        } else {
            0.0
        }
    };
    // Unlike in xsf, e1 and e2 also have the rounding errors d of the arguments z of gamma2, times
    // |digamma(z)|, which is about 1 / |z - n| near a pole n (e.g. for 1 + a - b) and ln(z) for
    // large z (e.g. for 1 - b = 128.1, with d = 1.4e-14), and that of 1 - b times ln(x) for x^(1-b)
    let psi_err = |z: f64, d: f64| {
        if d == 0.0 {
            0.0
        } else {
            (digamma(z) * d).abs()
        }
    };
    let (s, d_s) = two_sum(1.0, a);
    let (aab, d_aab) = two_sum(s, -b);
    let (bm, d_bm) = two_sum(1.0, -b);
    let e_r1 = if b < -1.0 { psi_err(bm, d_bm) } else { 0.0 };
    let e_c = if b > 3.0 {
        psi_err(-bm, d_bm)
    } else {
        let (b2, d_b2) = two_sum(2.0, -b);
        psi_err(b2, d_b2)
    };
    let e1 = 1e-15 + k(1.0 + a - b) + psi_err(aab, d_s + d_aab) + e_r1;
    let e2 = 1e-15 + k(a) + k(b) + k(2.0 - b) + e_c + (d_bm * x.ln()).abs();
    let err = e1 * s1.abs() + e2 * s2.abs() + f64::EPSILON * terms;
    *id = if lost {
        0
    } else {
        (*id).min((-(err / hu.abs()).log10()) as i32)
    };
    lr2.map_or(hu, |l| mul_exp_dd(hu, l))
}

/// `chgus` uses `chgus_near_1` for `|b - 1|` below this, the largest its `gamma_quot(.., 2e)` allows
const NEAR_1_DD: f64 = 2e-2;

/// `chgu` returns a `chgus_near_1` with fewer digits than this only after its other methods, for
/// `a >= 1` or `x >= 1` (empirical; e.g. it has 12 at `U(16.5, 0.994, 0.1)`, and an error of
/// 5.6e-14, vs 9.2e-16 in `chguit`)
const NEAR_1_ID_DEFER: i32 = 14;

/// The relative error that `chgus_near_1` assumes for each of its terms
const NEAR_1_ROUND: f64 = 8.0 * f64::EPSILON;

/// `U(a, 1 + e, x)` by DLMF 13.2.42 (DLMF 13.2.9 for `e = 0`), as in `chgus`, but with the `k`th
/// terms of its two series, both of order `1 / e`, combined into `e pi / sin(pi e) d_k` without
/// cancellation; `None` if `1 / Gamma(a)` under- or overflows (unless at a pole), at a pole of
/// `Gamma(a - e)`, or if the sum isn't finite
fn chgus_near_1(x: f64, a: f64, e: f64, id: &mut i32) -> Option<f64> {
    if a.abs() > MAXGAM {
        return None;
    }
    let (ra, rp, pole) = (rgamma(a), rgamma(1.0 + e), is_nonpos_int(a));
    if !ra.is_normal() && !pole {
        return None;
    }
    // Gamma(1 + e) / Gamma(1 - e) = 1 + 2e q_1, and s = (1 / Gamma(a) - 1 / Gamma(a - e)) / e
    // from Gamma(a) / Gamma(a - e) = 1 + e q_2, or at a pole a = -n by DLMF 5.5.3 as
    // (-1)^n Gamma(1 - a + e) sin(pi e) / (pi e)
    let (dd_1, dd_e) = (DoubleDouble::new(1.0), DoubleDouble::new(e));
    let q_1 = gamma_quot(dd_1 - e, dd_e * 2.0)?.0;
    let rm = rp * (1.0 + 2.0 * e * q_1);
    let s = if pole {
        let q = gamma_quot(dd_1 - a, dd_e)?.0;
        neg_if_odd(a, gamma(1.0 - a) * (1.0 + e * q) * sinc(e))
    } else {
        -ra * gamma_quot(DoubleDouble::new(a) - e, dd_e)?.0
    };
    let lx = x.ln();
    // d_k = (p_k - q_k) / e, with p_k = x^(k-e) (a - e)_k / (k! Gamma(a) Gamma(1 - e + k)) and
    // q_k = x^k (a)_k / (k! Gamma(a - e) Gamma(1 + e + k)), so d_0 is the sum of
    // (x^-e - 1) / (e Gamma(a) Gamma(1 - e)), (1 / Gamma(1 - e) - 1 / Gamma(1 + e)) / (e Gamma(a))
    // and s / Gamma(1 + e)
    let parts = [-lx * exprel(-e * lx) * ra * rm, 2.0 * q_1 * rp * ra, s * rp];
    let mut d = parts[0] + parts[1] + parts[2];
    let mut p = (-e * lx).exp() * ra * rm;
    // `m` bounds the error of d_k, with each rounding error (`NEAR_1_ROUND`) taken as positive
    let mut m = NEAR_1_ROUND * parts.iter().map(|t| t.abs()).sum::<f64>();
    let (mut sum, mut err) = (d, m);
    *id = 0;
    // as many terms as `chgus` sums at most
    for k in 0..300 {
        let k = f64::from(k);
        // p_(k+1) = p_k u x / (k + 1) and q_(k+1) = q_k v x / (k + 1), with
        // u = (a - e + k) / (k + 1 - e), v = (a + k) / (k + 1 + e) and
        // u - v = e (2a + k - 1 - e) / ((k + 1 - e) (k + 1 + e))
        let v = (a + k) / (k + 1.0 + e);
        let w = (2.0 * a + k - 1.0 - e) / ((k + 1.0 - e) * (k + 1.0 + e));
        let t = x / (k + 1.0);
        d = (d * v + p * w) * t;
        m = (m * v.abs() + NEAR_1_ROUND * (p * w).abs()) * t;
        p *= (a - e + k) / (k + 1.0 - e) * t;
        sum += d;
        err += m;
        if k > x && d.abs() + m / NEAR_1_ROUND <= 0.25 * f64::EPSILON * sum.abs() {
            *id = (-(err / sum.abs()).log10()) as i32;
            break;
        }
    }
    sum.is_finite().then_some(sum / sinc(e))
}

#[cfg(test)]
mod tests {
    use super::{chgm_neg_x_asymptotic, chgu, chguit, chgus};

    /// Whether the relative error is within the estimated digits
    fn digits_ok(u: f64, expected: f64, id: i32) -> bool {
        (u / expected - 1.0).abs() < 10f64.powi(-id)
    }

    #[test]
    fn test_chgm_neg_x_asymptotic_poles() {
        // `b - a` is -2 (a pole of `Gamma`), and -1 - 8.3e-17 (next to one)
        assert_eq!(chgm_neg_x_asymptotic(-1e5, 2.5, 0.5), None);
        assert!(chgm_neg_x_asymptotic(-1e5, 1.1, 0.1).is_some());
        assert!(chgm_neg_x_asymptotic(-1e5, 1.5, 1.0).is_some());
    }

    #[test]
    fn test_chguit_digits() {
        let mut id = 0;
        // U(2.25, 0, 0.5) (Arb), with 15 digits (9 in xsf)
        let u = chguit(0.5, 2.25, 0.0, false, &mut id);
        assert!(id >= 13 && digits_ok(u, 0.123_369_748_576_127_63, id));
        // x^(b-1) U(a, b, x) (Arb), where the error oscillates over the refinements, so that one
        // change is small by chance (9.548_667_857_677_493e-20 with 9 digits if that is final)
        let (x, a, b) = (6.294_217_108_153_153e-4, 21.564_919_933_139_418, 1.579_409_465_206_605_7);
        let u = chguit(x, a, b, true, &mut id);
        assert!(digits_ok(u, 9.548_667_779_662_119_8e-20, id));
    }

    #[test]
    fn test_chgus_digits() {
        // U(a, b, x) = 6.67e-286 (Arb), where the r2 terms underflow from j = 64 before they grow
        // again after j = b - 1, so that the series can't converge in time
        let mut id = 0;
        chgus(55.0, 164.9, 168.900_01, &mut id);
        assert!(id <= 0);

        // U(a, b, x) (Arb), where the rounding error 1.4e-14 of 1 - b = 128.1 makes Gamma(1 - b)
        // 6.9e-14 off (13.2 digits)
        let (x, a, b) = (1.922_426_784_912_072_6, 16.014_115_853_033_63, -127.117_796_766_930_57);
        let mut id = 0;
        let u = chgus(x, a, b, &mut id);
        assert!(digits_ok(u, 5.652_129_176_700_298e-35, id));

        // U(a, b, x) (Arb), near b = 1, where 1 / Gamma(a) is near the underflow threshold
        let u = chgus(1e-3, 171.3, 1.001, &mut id);
        assert!(id >= 13 && digits_ok(u, 3.223_069_881_767_204_4e-308, id));
    }

    #[test]
    fn test_chgu_subnormal_digits() {
        // U(170, 1 + 1e-11, 1) (Arb), by `chguit`, is subnormal, with 7.7 digits
        let (mut md, mut isfer, mut digits) = (0, 0, 0);
        let u = chgu(1.0, 170.0, 1.000_000_000_01, &mut md, &mut isfer, &mut digits);
        assert_eq!((md, isfer), (4, 0));
        assert!(digits_ok(u, 9.220_520_059_940_690_4e-317, digits));
        // the same for an early return, here of `chgul` (Arb)
        let u = chgu(1500.0, 100.0, 101.0, &mut md, &mut isfer, &mut digits);
        assert_eq!(md, 2);
        assert!(digits_ok(u, 2.459_654_426_579_829_3e-318, digits));
    }
}
