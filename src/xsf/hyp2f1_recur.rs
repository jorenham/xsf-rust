//! The Gauss hypergeometric function `F(a, b; c; x)` for real `x < 1` and large real parameters,
//! by backward recurrence in `c`.
//!
//! This is not a port; the real `hyp2f1` (in `cephes/hyp2f1.rs`) tries it first for large
//! parameters, since its series and transformations can lose all precision there.
//!
//! # Algorithm
//!
//! 1. Write `F(a, b; c; x)` as `(1 - x)^s (1 - z)^p F(p, q; c; z)` with `z <= 0` (DLMF 15.8.1),
//!    where
//!
//!    - for `x < 0`, `z = x`, `(p, q)` is `(a, b)` or `(c - a, c - b)`, and `s` is `-a` or `-b`
//!    - for `0 < x < 1`, `z = x / (x - 1)`, `(p, q)` is `(a, c - b)` or `(c - a, b)`, and `s` is
//!      `0` or `c - a - b`
//!
//!    and likewise with `p` and `q` swapped, i.e. 4 combinations.
//!
//! 2. With `w = z / (z - 1)` (so `0 <= w < 1`, and `w = x` if `0 < x < 1`), and any `C`,
//!
//!    `(1 - z)^p F(p, q; C; z) = F(p, C - q; C; w)`        (DLMF 15.8.1)
//!
//!    For `k > -p`, the terms of this series have the same sign as long as `C + k` and `C - q + k`
//!    have the same sign, e.g. for all `k` if `C > max(0, q)`. So there is little cancellation
//!    if `p >= 0`, or if `p < 0` is small. We sum it at `C = c`.
//!
//! 3. If `c <= max(0, q)` and that loses too much precision, we also sum it at `C = c + n` and
//!    `C = c + n + 1`, with the smallest `n` for which `c + n > max(0, q)`, and recur down to
//!    `C = c` with the contiguous relation in `c` (DLMF 15.5.18). For `Re z < 1/2`,
//!    `F(p, q; c + n; z)` is the minimal solution of this recurrence as `n -> inf` (Gil, Segura and
//!    Temme 2007, sec. 3.2), so recurring towards smaller `c` is the stable direction. But this is
//!    only asymptotic, for `c + n` much larger than `|p|` and `|q|`, and the recurrence can be
//!    unstable in between (Deano and Segura 2007). So instead of predicting this, we bound the
//!    error of the result to first order, see `recur_down`.
//!
//! 4. For `x < -1`, `w` is close to 1 if `|x|` is large, and the series take about `|p| |x|`
//!    terms. So we first try the 1/x transformation (DLMF 15.8.2), whose two hypergeometric
//!    functions of `1/x` are computed as above, with `w = 1 / (1 - x)`, and only continue with the
//!    4 combinations if its error bound (incl. the cancellation of its two terms) exceeds
//!    `HYP2F1_ETHRESH`.
//!
//! 5. If `x >= -1` and all error bounds so far exceed `HYP2F1_RECUR_TOL`, also try the same
//!    transformation of `F(a, b; c; x)` if `x < 0` (i.e. DLMF 15.8.3), or of the Pfaff form
//!    `(1 - x)^-b F(c - a, b; c; x / (x - 1))` if `0 < x < 1` (i.e. the 1 - x transformation,
//!    DLMF 15.8.4). Its series are in `w = 1 / (1 - x)` or `w = 1 - x`. The 4 combinations can
//!    all cancel if `a` and `b` are large with opposite signs and `0 < x < 1` (or with the same
//!    sign and `x < 0`), e.g. `F(-975.2, 241.4; 0.00135; 0.39)`, while these two terms rarely do.
//!
//! 6. Return the result with the smallest error bound, or `None` if that exceeds
//!    `HYP2F1_RECUR_TOL`, for the caller to fall back to another method.
//!
//! Some details:
//!
//! - The parameters are double-double, so that e.g. `C - q + k` is accurate if `q = c - b` and
//!   `c` is tiny, and `(1 - x)^s` is computed in double-double, since `|s ln(1 - x)|` can be large.
//!   For the same reason, the `ln Gamma` of steps 4 and 5 are computed in double-double for large
//!   arguments.
//! - The series and the recurrence keep their values as `m * 2^e`, to avoid overflow and
//!   underflow.
//! - The cost is bounded by `HYP2F1_RECUR_MAXITER` series terms and recurrence steps per
//!   combination. It's about `(|p| w + 40) / (1 - w)` terms, plus `max(0, q) - c` steps.
//!
//! # References
//!
//! - A. Gil, J. Segura and N.M. Temme, "Numerically satisfactory solutions of hypergeometric
//!   recursions", Math. Comp. 76 (2007), 1449-1468. <https://doi.org/10.1090/S0025-5718-07-01918-7>
//! - A. Deano and J. Segura, "Transitory minimal solutions of hypergeometric recursions and
//!   pseudoconvergence of associated continued fractions", Math. Comp. 76 (2007), 879-901.
//!   <https://doi.org/10.1090/S0025-5718-07-01934-5>
//! - NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/5.5.E3>,
//!   <https://dlmf.nist.gov/5.11.E1>, <https://dlmf.nist.gov/15.5>, <https://dlmf.nist.gov/15.8>

use core::f64::consts::{LN_2, PI};

use crate::utils::MACHEP;
use crate::xsf::cephes::dd_real::{
    DoubleDouble as Dd, LN_PI_DD, LN_SQRT_2PI_DD, LOG2, frexp_f64, ldexp_f64,
};
use crate::xsf::cephes::{HYP2F1_ETHRESH, lgam_sgn, lgam_stirling_corr, psi};
use crate::xsf::trig::{cospi, sinpi};

/// The largest accepted error bound. These bounds are pessimistic (the actual errors are usually
/// 10 to 300 times smaller), and the methods of the caller are often far worse if this fails.
const HYP2F1_RECUR_TOL: f64 = 1e-10;

/// The maximum number of series terms plus recurrence steps, per combination
const HYP2F1_RECUR_MAXITER: f64 = 1e6;

/// `m * 2^e`, with a bound of its relative error
#[derive(Clone, Copy, Debug)]
struct Scaled {
    m: f64,
    e: i32,
    err: f64,
}

/// The exponent `k` that scales `|x|` to about 1 by `2^-k`, if `|x|` is far from 0 and 1
#[allow(clippy::cast_possible_truncation)]
fn rescale_exp(x: f64) -> Option<i32> {
    // (|x| is finite and nonzero, so log2(|x|) is in [-1075, 1024])
    let big = 2.0_f64.powi(500);
    let x = x.abs();
    (x != 0.0 && x.is_finite() && !(1.0 / big..=big).contains(&x)).then(|| x.log2().round() as i32)
}

/// `v + k`, rounded to `f64` (if `v` is a double, `v.hi + k` is the same, and faster)
fn add_hi(v: Dd, k: f64) -> f64 {
    if v.lo == 0.0 { v.hi + k } else { (v + k).hi }
}

/// The one of `u` and `v` with the smaller error bound `err`, or `u` if they're equal
fn min_err<T>(u: Option<T>, v: Option<T>, err: impl Fn(&T) -> f64) -> Option<T> {
    match (u, v) {
        (Some(u), Some(v)) => Some(if err(&u) <= err(&v) { u } else { v }),
        (u, v) => u.or(v),
    }
}

/// `F(p, r; c; w)` for `0 <= w < 1`, by summing its series, with an estimate of its rounding
/// errors; `budget` is the remaining number of terms.
///
/// - The arguments are double-double, so that e.g. `r + k` is accurate if `r = c - b` and `c` is
///   tiny.
/// - The terms are evaluated at `w.hi`, and the sum is corrected to first order in `w.lo`, since it
///   can be much more sensitive to `w` than `F(a, b; c; x)` itself.
/// - Each term `t_k` has a relative rounding error of at most `8 eps` more than `t_(k-1)` (the 3
///   roundings of `p + k`, `r + k` and `c + k`, and 5 operations), so at most `8 k eps`, and the
///   sum of `k` terms adds at most `k eps` of the sum of `|t_j|`.
#[allow(clippy::cast_precision_loss)]
fn series_w(p: Dd, r: Dd, c: Dd, w: Dd, budget: &mut f64) -> Option<Scaled> {
    let (w, w_lo) = (w.hi, w.lo);
    // from here on, the factors of the term ratios have constant signs
    let kmin = (-p.hi).max(-r.hi).max(-c.hi).max(0.0).ceil();
    // The terms are `t 2^te`, and the sums are scaled by `2^-e`, since the terms can underflow
    // before `k = -c`, and be huge after it.
    let (mut sum, mut sum_abs, mut k_sum, mut k_sum_abs) = (1.0, 1.0_f64, 0.0, 0.0_f64);
    let (mut t, mut te, mut e, mut k) = (1.0_f64, 0, 0, 0.0);
    // the two factors of the next term ratio, `(p + k) / (k + 1)` and `(r + k) / (c + k)`, and
    // their numerators
    let (mut pk, mut rk) = (p.hi, r.hi);
    let (mut f0, mut f1) = (pk, rk / c.hi);
    loop {
        // Multiply by the two factors separately, with `|t|` in `[2^-500, 2^500]`, or in `[1, 2)`
        // if a factor is outside that range, so that the product can't underflow (e.g. for tiny p
        // and r), and fail if a factor underflows, or a factor or the product overflows.
        for (f, num) in [(f0, pk), (f1 * w, rk)] {
            if !f.is_finite() || f.abs() < f64::MIN_POSITIVE && num != 0.0 {
                return None;
            }
            if t != 0.0 && rescale_exp(f).is_some() {
                let (m, s) = frexp_f64(t);
                (t, te) = (2.0 * m, te + s - 1);
            }
            t *= f;
            if !t.is_finite() {
                return None;
            }
            if let Some(s) = rescale_exp(t) {
                (t, te) = (ldexp_f64(t, -s), te + s);
            }
        }
        let term = if te == 0 { t } else { ldexp_f64(t, te) };
        k += 1.0;
        sum += term;
        sum_abs += term.abs();
        k_sum += k * term;
        k_sum_abs += k * term.abs();
        if !k_sum_abs.is_finite() {
            return None;
        }
        (pk, rk) = (add_hi(p, k), add_hi(r, k));
        (f0, f1) = (pk / (k + 1.0), rk / add_hi(c, k));
        // For j >= k > kmin, the next term ratios are w (p + j) / (j + 1) * (r + j) / (c + j),
        // where both factors are positive and tend monotonically to 1, so they are at most `rho`,
        // which bounds the rest of the series geometrically.
        if k > kmin {
            let rho = w * f0.max(1.0) * f1.max(1.0);
            if rho < 1.0 && term.abs() * rho / (1.0 - rho) <= sum.abs() * MACHEP / 2.0 {
                break;
            }
        }
        *budget -= 1.0;
        if *budget < 0.0 {
            return None;
        }
        if let Some(s) = rescale_exp(sum_abs) {
            (sum, sum_abs) = (ldexp_f64(sum, -s), ldexp_f64(sum_abs, -s));
            (k_sum, k_sum_abs) = (ldexp_f64(k_sum, -s), ldexp_f64(k_sum_abs, -s));
            (e, te) = (e + s, te - s);
        }
    }
    // dF/dw = (sum of k t_k) / w
    if w != 0.0 {
        sum += k_sum * (w_lo / w);
    }
    let err = MACHEP * (8.0 * k_sum_abs + (k + 1.0) * sum_abs) / sum.abs();
    (sum != 0.0 && err.is_finite()).then_some(Scaled { m: sum, e, err })
}

/// `F(p, q; c; z)` by recurring down from `y_n = F(p, q; c + n; z)` and `y_(n+1)`, using
/// DLMF 15.5.18 (with `C = c + k`):
///
/// ```text
/// A_k y_(k-1) + B_k y_k + C_k y_(k+1) = 0
/// A_k = C (C - 1) (z - 1),  B_k = C (C - 1 - (2 C - p - q - 1) z),  C_k = (C - p) (C - q) z
/// ```
///
/// The relative error estimate is a first-order bound:
///
/// - Let `g` be another solution, and `d_k = g_(k+1) / y_(k+1) - g_k / y_k`. A relative error
///   `t_m` of `y_m` (with `y_(m+1)` exact) changes `y_0` by the relative error
///   `t_m (d_0 + ... + d_m) / d_m`. This is small if `y` is the dominant solution in the
///   direction of the recurrence, and large if it is the minimal one.
/// - The Casoratian `W_k = y_k g_(k+1) - y_(k+1) g_k` satisfies `C_k W_k = A_k W_(k-1)`, so
///   `d_(k-1) / d_k = C_k / A_k * y_(k+1) / y_(k-1)`.
/// - So `u_m = (d_0 + ... + d_m) / d_m` satisfies `u_0 = 1` and `u_m = 1 + u_(m-1) d_(m-1) / d_m`,
///   which we evaluate on the way back up. The bound is the sum of `|t_m u_m|`, plus the errors
///   of `y_n` and `y_(n+1)` times `|u_n|` and `|1 - u_n|`.
/// - Where `y` oscillates in `k`, the `d_k` change sign, and `|u_m|` can stay small, while
///   `(|d_0| + ... + |d_m|) / |d_m|` (which could be accumulated on the way down) grows like `m`.
#[allow(clippy::many_single_char_names)]
fn recur_down(p: Dd, q: Dd, c: Dd, z: f64, n: u32, y_n: Scaled, y_n1: Scaled) -> Option<Scaled> {
    // y_(k+1) and y_k, scaled by 2^-e
    let (mut y1, mut y0, mut e) = (ldexp_f64(y_n1.m, y_n1.e - y_n.e), y_n.m, y_n.e);
    // `t_(k-1)` and `d_(k-1) / d_k`, for the error bound
    let mut steps = Vec::with_capacity(n as usize);
    for k in (1..=n).rev() {
        // (C - 1, C - p and C - q are computed from c, since they can be tiny)
        let (cc, cc1) = (add_hi(c, f64::from(k)), add_hi(c, f64::from(k - 1)));
        let cc_p = (-p + c.hi + f64::from(k) + c.lo).hi;
        let cc_q = (-q + c.hi + f64::from(k) + c.lo).hi;
        let coef_a = cc * cc1 * (z - 1.0);
        let coef_b = cc * (cc1 - (cc_p + cc1 - q.hi) * z);
        let coef_c = cc_p * cc_q * z;
        let (t0, t1) = (coef_b * y0, coef_c * y1);
        let ym = -(t0 + t1) / coef_a;
        if ym == 0.0 || !(ym.is_finite() && t0.is_finite() && t1.is_finite()) || coef_a == 0.0 {
            return None;
        }
        // the rounding errors of the step, and of `coef_b` (incl. its cancellation)
        let err_b = cc * (cc1.abs() + (cc_p.abs() + cc1.abs() + q.hi.abs()) * z.abs());
        // (the roundings of t0 and t1 and their sum (6), of the parts of `coef_b` (5), and of
        // `coef_a` (C, C - 1, z - 1, 2 products) and the division (6))
        let t_m = MACHEP * (6.0 * (t0.abs() + t1.abs()) + 5.0 * (err_b * y0).abs())
            / (coef_a * ym).abs()
            + 6.0 * MACHEP;
        // d_(k-1) / d_k, which is 0 if C_k = 0, where the Casoratian vanishes
        let rho = coef_c / coef_a * (y1 / ym);
        if rho == 0.0 || !rho.is_finite() {
            return None;
        }
        steps.push((t_m, rho));
        (y1, y0) = (y0, ym);
        if let Some(s) = rescale_exp(y0) {
            (y0, y1) = (ldexp_f64(y0, -s), ldexp_f64(y1, -s));
            e += s;
        }
    }
    let (mut u, mut err) = (1.0_f64, 0.0);
    for (t_m, rho) in steps.into_iter().rev() {
        err += t_m * u.abs();
        u = 1.0 + u * rho;
    }
    err += y_n.err * u.abs() + y_n1.err * (1.0 - u).abs();
    err.is_finite().then_some(Scaled { m: y0, e, err })
}

/// `(1 - z)^p F(p, q; c; z)` for `z <= 0` (steps 2 and 3), with `w = z / (z - 1)`
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn hyp2f1_recur_form(p: Dd, q: Dd, c: Dd, z: f64, w: Dd) -> Option<Scaled> {
    // The series takes about (|p| w + 40) / (1 - w) terms, so give up early if that's too many:
    // the terms peak near k = |p| w / (1 - w), and then decay like w^k, which takes about
    // ln(2^53) / (1 - w) ~ 37 / (1 - w) more terms.
    if (p.hi.abs() * w.hi + 40.0) / (1.0 - w.hi) > HYP2F1_RECUR_MAXITER {
        return None;
    }
    // the number of series terms plus recurrence steps left
    let mut budget = HYP2F1_RECUR_MAXITER;
    let r = -q + c.hi + c.lo;
    let direct = series_w(p, r, c, w, &mut budget);
    if c.hi > q.hi.max(0.0) || direct.is_some_and(|y| y.err <= HYP2F1_ETHRESH) {
        return direct;
    }
    // the smallest n with c + n > max(0, q)
    let n = (q.hi.max(0.0) - c.hi).floor() + 1.0;
    budget -= n;
    if budget < 0.0 {
        return direct;
    }
    let (c_n, r_n) = (c + n, r + n);
    let recur = series_w(p, r_n, c_n, w, &mut budget).and_then(|y_n| {
        let y_n1 = series_w(p, r_n + 1.0, c_n + 1.0, w, &mut budget)?;
        // (the recurrence is symmetric in p and q)
        recur_down(p, q, c, z, n as u32, y_n, y_n1)
    });
    min_err(direct, recur, |y| y.err)
}

/// The 4 combinations `(p, q, s)` of step 1 for `x = z < 0`
fn neg_z_combinations(a: Dd, b: Dd, c: Dd) -> [(Dd, Dd, Dd); 4] {
    let (ca, cb) = (c - a, c - b);
    [(a, b, -a), (b, a, -b), (ca, cb, -b), (cb, ca, -a)]
}

/// The one of `combinations` `(p, q, s)` (step 1) with the smallest error bound (steps 2, 3 and
/// 6), as `(y, s ln(1 - x))`; `ln_1mx` is `ln(1 - x)`, and `z`, `w` are as in the module docs
fn best_combination(
    mut combinations: [(Dd, Dd, Dd); 4],
    c: Dd,
    z: f64,
    w: Dd,
    ln_1mx: Dd,
) -> Option<(Scaled, Dd)> {
    // First try those with p >= 0 and small |p|, whose series cancel less and have fewer terms, and
    // stop at the first with an error bound below `HYP2F1_ETHRESH`.
    combinations.sort_by(|u, v| {
        (u.0.hi < 0.0, u.0.hi.abs())
            .partial_cmp(&(v.0.hi < 0.0, v.0.hi.abs()))
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let mut best: Option<(Scaled, Dd)> = None;
    for (p, q, s) in combinations {
        let Some(y) = hyp2f1_recur_form(p, q, c, z, w) else {
            continue;
        };
        best = min_err(best, Some((y, s * ln_1mx)), |y| y.0.err);
        if y.err <= HYP2F1_ETHRESH {
            break;
        }
    }
    best
}

/// `ln(v)` for `v > 0`, after scaling `v` to `[1, 2)` by a power of 2, since `DoubleDouble::log`
/// fails for `v` beyond about `e^709` (as it computes `exp(-ln(v))`)
fn ln_dd(v: Dd) -> Dd {
    let (_, k) = frexp_f64(v.hi);
    v.ldexp(1 - k).log() + LOG2 * f64::from(k - 1)
}

/// `ln|Gamma(g)|`, the sign of `Gamma(g)`, and a bound of the absolute error of the former, for
/// a double-double `g`, or `None` at (or too close to) a pole:
///
/// - for `g >= 13`, by Stirling's formula in double-double, with the correction of cephes' `lgam`
///   (<https://dlmf.nist.gov/5.11.E1>)
/// - for `g <= -13`, by the reflection formula `Gamma(g) Gamma(1 - g) = pi / sin(pi g)`
///   (<https://dlmf.nist.gov/5.5.E3>), to first order in `g.lo`
/// - else by `lgam` at `g.hi`, to first order in `g.lo`
///
/// The second-order terms `f'' g.lo^2 / 2` are bounded on the interval between `g.hi` and
/// `g.hi + g.lo`: with `d` the smallest distance to a pole there, `f'' < 1/d^2 + pi^2` for
/// `f = ln Gamma` and `g < 1/2`, `f'' < 1/g + 1/g^2` for `g > 0`, and `|f''| < (pi / 2d)^2` for
/// `f = ln|sin(pi g)|` (since `|sin(pi g)| >= 2d`).
fn ln_gamma(g: Dd) -> Option<(Dd, f64, f64)> {
    if g.hi >= 13.0 {
        let v = (g - 0.5) * ln_dd(g) - g + LN_SQRT_2PI_DD + lgam_stirling_corr(g.hi);
        // (the correction is below 0.01, and its error below eps)
        return Some((v, 1.0, 2.0 * MACHEP + 16.0 * MACHEP * MACHEP * v.hi.abs()));
    }
    let d = (g.hi - g.hi.round()).abs() - g.lo.abs();
    if g.hi < 0.5 && d <= 0.0 {
        return None;
    }
    if g.hi <= -13.0 {
        let s = sinpi(g.hi);
        let (v, _, err) = ln_gamma(-g + 1.0)?;
        // d/dg ln|sin(pi g)| = pi cot(pi g)
        let dsin = PI * cospi(g.hi) / s * g.lo;
        let ln_sin = s.abs().ln();
        let v = LN_PI_DD - v - (ln_sin + dsin);
        let err = err
            + 4.0 * MACHEP * (1.0 + ln_sin.abs() + dsin.abs())
            + (PI / (2.0 * d) * g.lo).powi(2) / 2.0;
        return Some((v, s.signum(), err));
    }
    let (v, sign) = lgam_sgn(g.hi);
    if !v.is_finite() {
        return None;
    }
    let dpsi = psi(g.hi) * g.lo;
    let f2 = if g.hi < 0.5 {
        d.powi(-2) + PI * PI
    } else {
        let g = g.hi - g.lo.abs();
        (1.0 + 1.0 / g) / g
    };
    let err = 4.0 * MACHEP * (v.abs() + dpsi.abs() + 1.0) + f2 * g.lo * g.lo / 2.0;
    Some((Dd::new(v) + dpsi, f64::from(sign), err))
}

/// `F(a, b; c; t)` for `t < 0` by the 1/t transformation (steps 4 and 5), as `(y, l)` with the
/// value `y exp(l)`, where `t = x` if `x < 0`, and `t = x / (x - 1)` if `0 < x < 1`; see
/// <https://dlmf.nist.gov/15.8.E2>:
///
/// ```text
/// F(a, b; c; t) = G(a, b) (-t)^-a F(a, a - c + 1; a - b + 1; 1/t) + G(b, a) (-t)^-b F(b, ...)
/// G(a, b) = Gamma(c) Gamma(b - a) / (Gamma(b) Gamma(c - a))
/// ```
///
/// The error bound includes those of the two `F` and of the `ln Gamma` (of their double-double
/// arguments), and the cancellation of the two terms. It fails if `b - a` is an integer (a pole
/// of `G`), and its bound is large if `b - a` is close to one.
fn hyp2f1_recip(a: Dd, b: Dd, c: Dd, x: f64) -> Option<(Scaled, Dd)> {
    let one = Dd::new(1.0);
    let (ln_x, ln_1mx) = (ln_dd(Dd::new(x.abs())), ln_dd(one - x));
    // 1/t, w = (1/t) / (1/t - 1) = 1 / (1 - t), ln(-t), and ln(1 - 1/t), computed from x, since
    // e.g. w = 1 - x if 0 < x < 1
    let (z, w, ln_mt, ln_1mz) = if x < 0.0 {
        (1.0 / x, one / (one - x), ln_x, ln_1mx - ln_x)
    } else {
        ((x - 1.0) / x, one - x, ln_x - ln_1mx, -ln_x)
    };
    // (sign, y, l) of a term `sign y exp(l)`, and the error bound of its `ln Gamma`
    let term = |a: Dd, b: Dd| {
        // F(p, q; cc; 1/t)
        let (p, q, cc) = (a, a - c + 1.0, a - b + 1.0);
        let (y, l) = best_combination(neg_z_combinations(p, q, cc), cc, z, w, ln_1mz)?;
        let (mut sign, mut ln_g, mut ln_g_err) = (1.0, Dd::new(0.0), 0.0);
        for (g, power) in [(c, 1.0), (b - a, 1.0), (b, -1.0), (c - a, -1.0)] {
            let (v, sg, err) = ln_gamma(g)?;
            (ln_g, sign, ln_g_err) = (ln_g + v * power, sign * sg, ln_g_err + err);
        }
        Some((sign, y, l - a * ln_mt + ln_g, ln_g_err))
    };
    let ((s1, y1, l1, g1), (s2, y2, l2, g2)) = (term(a, b)?, term(b, a)?);
    // the two terms, scaled by `exp(-lmax)`, with `lmax` the larger exponent
    let (l1, l2) = (l1 + LOG2 * f64::from(y1.e), l2 + LOG2 * f64::from(y2.e));
    let lmax = if l1.hi >= l2.hi { l1 } else { l2 };
    let exp_rel = |l: Dd| {
        let r = l - lmax;
        r.hi.exp() * (1.0 + r.lo)
    };
    let (t1, t2) = (s1 * y1.m * exp_rel(l1), s2 * y2.m * exp_rel(l2));
    let m = t1 + t2;
    // (3 eps for `exp_rel` and the products, and 1 eps for the sum)
    let err = (t1.abs() * (y1.err + g1 + 3.0 * MACHEP) + t2.abs() * (y2.err + g2 + 3.0 * MACHEP))
        / m.abs()
        + MACHEP;
    (m != 0.0 && err.is_finite()).then_some((Scaled { m, e: 0, err }, lmax))
}

/// `F(a, b; c; x)` for real `x < 1` with `x != 0` (see the module docs), or `None` if the
/// estimated relative error is too large.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub(crate) fn hyp2f1_recur(a: f64, b: f64, c: f64, x: f64) -> Option<f64> {
    // (beyond this, the coefficients of the recurrence can overflow)
    if !(a.abs().max(b.abs()).max(c.abs()) <= 1e15 && x < 1.0 && x != 0.0) {
        return None;
    }
    let (a, b, c) = (Dd::new(a), Dd::new(b), Dd::new(c));
    // the combinations `(p, q, s)` of step 1, with `F(a, b; c; x)` equal to
    // `(1 - x)^s (1 - z)^p F(p, q; c; z)`
    let (z, w, combinations) = if x < 0.0 {
        // (if x is huge, w rounds to 1, and the series fails)
        let w = Dd::new(x) / (Dd::new(x) - 1.0);
        (x, w, neg_z_combinations(a, b, c))
    } else {
        let (ca, cb, zero) = (c - a, c - b, Dd::new(0.0));
        let cab = ca - b;
        let combinations = [(a, cb, zero), (cb, a, cab), (ca, b, cab), (b, ca, zero)];
        (x / (x - 1.0), Dd::new(x), combinations)
    };
    // For x < -1, first try the 1/x transformation, whose series are much shorter if w is near 1.
    let recip = if x < -1.0 {
        hyp2f1_recip(a, b, c, x)
    } else {
        None
    };
    // (`log1p` of `DoubleDouble` is only as accurate as that of `f64`, but `log` is accurate)
    let ln_1mx = (Dd::new(1.0) - x).log();
    let best = if recip.is_some_and(|(y, _)| y.err <= HYP2F1_ETHRESH) {
        recip
    } else {
        min_err(
            best_combination(combinations, c, z, w, ln_1mx),
            recip,
            |y| y.0.err,
        )
    };
    // If that fails for x >= -1, also try the 1/t transformation (step 5).
    let best = if x >= -1.0 && !best.is_some_and(|(y, _)| y.err <= HYP2F1_RECUR_TOL) {
        let recip = if x < 0.0 {
            hyp2f1_recip(a, b, c, x)
        } else {
            // F(a, b; c; x) = (1 - x)^-b F(c - a, b; c; x / (x - 1))
            hyp2f1_recip(c - a, b, c, x).map(|(y, l)| (y, l - b * ln_1mx))
        };
        min_err(best, recip, |y| y.0.err)
    } else {
        best
    };
    let (y, l) = best?;
    if !(y.err <= HYP2F1_RECUR_TOL) {
        return None;
    }
    // `y 2^e exp(l)`, with `exp(l) = 2^k exp(r)` and `|r| <= ln(2) / 2`. If `|k|` is so large that
    // it's clamped, then `e + k` is clamped as well, since `|e|` is at most about 1e9 (2^1000 per
    // term), and the result overflows or underflows.
    let k = (l.hi / LN_2).round().clamp(-1e12, 1e12);
    let r = l - LOG2 * k;
    #[allow(clippy::cast_possible_truncation)]
    let e = (f64::from(y.e) + k).clamp(-1e4, 1e4) as i32;
    Some(ldexp_f64(y.m * r.hi.exp() * (1.0 + r.lo), e))
}

#[cfg(test)]
mod tests {
    use super::{Dd, Scaled, hyp2f1_recur, recur_down};

    #[test]
    fn test_recur_down_start_errors() {
        // y_0 = y_2 / 8 (as B_1 = 0), so the error of y_2 carries over, although u_1 = 0
        let y = |m, err| Scaled { m, e: 0, err };
        let (p, q, c) = (Dd::new(3.5), Dd::new(1.5), Dd::new(1.0));
        let y_0 = recur_down(p, q, c, -0.5, 1, y(0.5, 0.0), y(0.8, 1e-6)).unwrap();
        assert!((y_0.m - 0.1).abs() < 1e-15 && y_0.err >= 1e-6, "{y_0:?}");
    }

    #[test]
    fn test_hyp2f1_recur_overflow() {
        // the terms or their ratios over- or underflow: fail instead of returning e.g. NaN or 1.0
        // (direct sums of the series in high precision)
        for (a, b, c, x, desired) in [
            (10.1, 0.2, 1.4e-309, -0.18, -9.670_216_976_660_568e307),
            (5e-324, 2000.3, 10.7, 0.5, 6.216_729_240_191_533e249),
            (10.0, 9e306, 1e308, 0.5, 1.584_770_156_297_082),
            (0.125, 1.135, -20.2, -1e308, 7.103_393_504_675_01e-39),
            (1e-150, 1e-200, -1000.5, 0.8, -5.767_096_442_064_845e249),
            (1e-140, 1e-180, -1000.5, 0.8, -5.767_096_442_064_845e279),
        ] {
            let actual = hyp2f1_recur(a, b, c, x);
            assert!(
                actual.is_none_or(|y| (y / desired - 1.0).abs() < 1e-10),
                "{actual:?} != {desired}"
            );
        }
    }
}
