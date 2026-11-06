//! The Gauss hypergeometric function `F(a, b; c; x)` for real `x < 1` and large real parameters,
//! by backward recurrence in `c`. Not a port: `hyp2f1` (in `cephes/hyp2f1.rs`) tries it first for
//! large parameters, where its series and transformations can lose all precision.
//!
//! # Algorithm
//!
//! 1. Write `F(a, b; c; x) = (1 - x)^s (1 - z)^p F(p, q; c; z)` with `z <= 0` (DLMF 15.8.1), in 4
//!    combinations (these and the same with `a` and `b` swapped):
//!    - for `x < 0`: `z = x`, and `(p, q, s)` is `(a, b, -a)` or `(c - a, c - b, -b)`
//!    - for `0 < x < 1`: `z = x / (x - 1)`, and `(p, q, s)` is `(a, c - b, 0)` or
//!      `(c - a, b, c - a - b)`
//!
//! 2. With `w = z / (z - 1)` (so `0 <= w < 1`, and `w = x` if `x > 0`), for any `C`:
//!
//!    `(1 - z)^p F(p, q; C; z) = F(p, C - q; C; w)` (DLMF 15.8.1)
//!
//!    For `k > -p`, its terms have one sign while `C + k` and `C - q + k` do, e.g. for all `k` if
//!    `C > max(0, q)`. So it cancels little if `p >= 0` or `|p|` is small. We sum it at `C = c`.
//!
//! 3. If `c <= max(0, q)` and that loses too much precision, we also sum it at `C = c + n` and
//!    `C = c + n + 1`, with the smallest `n` with `c + n > max(0, q)`, and recur down to `C = c`
//!    (DLMF 15.5.18). For `Re z < 1/2`, `F(p, q; c + n; z)` is minimal as `n -> inf` (Gil, Segura
//!    and Temme 2007, sec. 3.2), so this is stable for `c + n` much larger than `|p|` and `|q|`,
//!    but can be unstable in between (Deano and Segura 2007), so `recur_down` bounds its error.
//!
//! 4. For `x < -1`, `w` is near 1 if `|x|` is large, and the series take about `|p| |x|` terms.
//!    So we first try the 1/x transformation (DLMF 15.8.2), with its two `F(..; 1/x)` computed as
//!    above (with `w = 1 / (1 - x)`), and only continue with the 4 combinations if its error bound
//!    (incl. the cancellation of its two terms) exceeds `HYP2F1_ETHRESH`.
//!
//! 5. If `x >= -1` and all error bounds so far exceed `HYP2F1_RECUR_TOL`, also try the same
//!    transformation of `F(a, b; c; x)` if `x < 0` (i.e. DLMF 15.8.3), or of the Pfaff form
//!    `(1 - x)^-b F(c - a, b; c; x / (x - 1))` if `0 < x < 1` (i.e. the 1 - x transformation,
//!    DLMF 15.8.4). Its series are in `w = 1 / (1 - x)` or `w = 1 - x`. The 4 combinations can
//!    all cancel if `a` and `b` are large with opposite signs and `0 < x < 1` (or with the same
//!    sign and `x < 0`), e.g. `F(-975.2, 241.4; 0.00135; 0.39)`, while these two terms rarely do.
//!
//! Finally, return the result with the smallest error bound, if that is at most `HYP2F1_RECUR_TOL`.
//!
//! The parameters are double-double, so that e.g. `C - q + k` is accurate if `q = c - b` and `c`
//! is tiny, and so are `s ln(1 - x)` and the `ln Gamma` of the transformations. The values are
//! kept as `m 2^e`, to avoid over- and underflow. The cost is at most `HYP2F1_RECUR_MAXITER`
//! series terms and recurrence steps per combination: about `(|p| w + 40) / (1 - w)` terms, plus
//! `max(0, q) - c` steps.
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

use core::cmp::Ordering;
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

/// The combinations `(p, q, s)` of step 1
type Combinations = [(Dd, Dd, Dd); 4];

/// `m 2^e exp(l)`, with a bound of its relative error
#[derive(Clone, Copy, Debug)]
struct Scaled {
    m: f64,
    e: i32,
    l: Dd,
    err: f64,
}

impl Scaled {
    /// `Some(m 2^e)` if `m` is nonzero and `err` is finite
    fn new(m: f64, e: i32, err: f64) -> Option<Self> {
        let l = Dd::new(0.0);
        (m != 0.0 && err.is_finite()).then_some(Self { m, e, l, err })
    }

    fn with_l(self, l: Dd) -> Self {
        Self { l, ..self }
    }
}

/// The error bound of `y`, or infinity for `None`
fn err_of(y: Option<Scaled>) -> f64 {
    y.map_or(f64::INFINITY, |y| y.err)
}

/// The one of `u` and `v` with the smaller error bound, or `u` if they're equal
fn min_err(u: Option<Scaled>, v: Option<Scaled>) -> Option<Scaled> {
    if err_of(v) < err_of(u) { v } else { u }
}

/// The first of `ys` with an error bound of at most `HYP2F1_ETHRESH`, or else the one with the
/// smallest error bound
fn best_of(ys: impl IntoIterator<Item = Option<Scaled>>) -> Option<Scaled> {
    let mut best = None;
    for y in ys {
        best = min_err(best, y);
        if err_of(best) <= HYP2F1_ETHRESH {
            break;
        }
    }
    best
}

/// The exponent `k` that scales `|x|` to about 1 by `2^-k`, if `|x|` is outside `[2^-500, 2^500]`
fn rescale_exp(x: f64) -> Option<i32> {
    // (|x| is finite and nonzero, so log2(|x|) is in [-1075, 1024])
    let big = 2.0_f64.powi(500);
    let x = x.abs();
    (x != 0.0 && x.is_finite() && !(1.0 / big..=big).contains(&x)).then(|| x.log2().round() as i32)
}

/// `t 2^e`, with `t` scaled to about 1 if `|t|` is outside `[2^-500, 2^500]`
fn rescaled(t: f64, e: i32) -> (f64, i32) {
    rescale_exp(t).map_or((t, e), |s| (ldexp_f64(t, -s), e + s))
}

/// `t 2^e` times `f`, or `None` if `f` underflows (unless its numerator `num` is 0), or `f` or the
/// product overflows. If `f` is far from 1, `|t|` is first scaled to `[1, 2)`, so that `t f` can't
/// underflow (e.g. for tiny `p` and `r` in `series_w`).
fn mul_scaled((t, e): (f64, i32), f: f64, num: f64) -> Option<(f64, i32)> {
    if !f.is_finite() || num != 0.0 && !f.is_normal() {
        return None;
    }
    let (t, e) = if t != 0.0 && rescale_exp(f).is_some() {
        let (m, s) = frexp_f64(t);
        (2.0 * m, e + s - 1)
    } else {
        (t, e)
    };
    let t = t * f;
    t.is_finite().then(|| rescaled(t, e))
}

/// `v + k`, rounded to `f64` (if `v` is a double, `v.hi + k` is the same, and faster)
fn add_hi(v: Dd, k: f64) -> f64 {
    if v.lo == 0.0 { v.hi + k } else { (v + k).hi }
}

/// `F(p, r; c; w)` for `0 <= w < 1` by its series, with a bound of its rounding errors; `budget`
/// is the number of terms left. The terms are evaluated at `w.hi`, and the sum is corrected to
/// first order in `w.lo`, since it can be much more sensitive to `w` than `F(a, b; c; x)` itself.
/// Each term `t_k` has at most `8 eps` more relative error than `t_(k-1)` (from `p + k`, `r + k`,
/// `c + k` and 5 operations), and summing `k` terms adds at most `k eps` of `sum_j |t_j|`.
fn series_w(p: Dd, r: Dd, c: Dd, w: Dd, budget: &mut f64) -> Option<Scaled> {
    let (w, w_lo) = (w.hi, w.lo);
    // from here on, the factors of the term ratios have constant signs
    let kmin = [-p.hi, -r.hi, -c.hi].into_iter().fold(0.0, f64::max).ceil();
    // The terms are `t 2^te`, and the sums are scaled by `2^-e`, since the terms can underflow
    // before `k = -c`, and be huge after it.
    let (mut sum, mut sum_abs, mut k_sum, mut k_sum_abs) = (1.0, 1.0_f64, 0.0, 0.0_f64);
    let (mut t, mut te, mut e, mut k) = (1.0, 0, 0, 0.0);
    // the factors `(p + k) / (k + 1)` and `(r + k) / (c + k)` of the next term ratio, and their
    // numerators
    let (mut pk, mut rk) = (p.hi, r.hi);
    let (mut f0, mut f1) = (pk, rk / c.hi);
    loop {
        (t, te) = mul_scaled((t, te), f0, pk)?;
        (t, te) = mul_scaled((t, te), f1 * w, rk)?;
        let term = if te == 0 { t } else { ldexp_f64(t, te) };
        let term_abs = term.abs();
        k += 1.0;
        sum += term;
        sum_abs += term_abs;
        k_sum += k * term;
        k_sum_abs += k * term_abs;
        if !k_sum_abs.is_finite() {
            return None;
        }
        (pk, rk) = (add_hi(p, k), add_hi(r, k));
        (f0, f1) = (pk / (k + 1.0), rk / add_hi(c, k));
        // For j >= k > kmin, both factors of the term ratio are positive and tend monotonically
        // to 1, so the ratio is at most `rho`, which bounds the rest of the series geometrically.
        if k > kmin {
            let rho = w * f0.max(1.0) * f1.max(1.0);
            if rho < 1.0 && term_abs * rho / (1.0 - rho) <= sum.abs() * MACHEP / 2.0 {
                break;
            }
        }
        *budget -= 1.0;
        if *budget < 0.0 {
            return None;
        }
        if let Some(s) = rescale_exp(sum_abs) {
            [sum, sum_abs, k_sum, k_sum_abs] =
                [sum, sum_abs, k_sum, k_sum_abs].map(|v| ldexp_f64(v, -s));
            (e, te) = (e + s, te - s);
        }
    }
    // dF/dw = (sum of k t_k) / w
    if w != 0.0 {
        sum += k_sum * (w_lo / w);
    }
    let err = MACHEP * (8.0 * k_sum_abs + (k + 1.0) * sum_abs) / sum.abs();
    Scaled::new(sum, e, err)
}

/// `F(p, q; c; z)` by recurring down from `y_n = F(p, q; c + n; z)` and `y_(n+1)` with
/// DLMF 15.5.18, where `C = c + k`:
///
/// ```text
/// A_k y_(k-1) + B_k y_k + C_k y_(k+1) = 0
/// A_k = C (C - 1) (z - 1)
/// B_k = C (C - 1 - (2 C - p - q - 1) z)
/// C_k = (C - p) (C - q) z
/// ```
///
/// The relative error bound is to first order, with `g` another solution and
/// `d_k = g_(k+1) / y_(k+1) - g_k / y_k`:
///
/// - A relative error `t_m` of `y_m` (with `y_(m+1)` exact) gives `y_0` the relative error
///   `t_m u_m`, with `u_m = (d_0 + ... + d_m) / d_m`, which is small if `y` is dominant in the
///   direction of the recurrence, and large if it is minimal.
/// - The Casoratian `W_k = y_k g_(k+1) - y_(k+1) g_k` satisfies `C_k W_k = A_k W_(k-1)`, so
///   `d_(k-1) / d_k = C_k / A_k * y_(k+1) / y_(k-1)`, and on the way back up we evaluate
///   `u_m = 1 + u_(m-1) d_(m-1) / d_m`.
/// - The bound is `sum_m |t_m u_m|`, plus the errors of `y_n` and `y_(n+1)` times `|u_n|` and
///   `|1 - u_n|`. Where `y` oscillates in `k`, the `d_k` change sign, and `|u_m|` can stay small,
///   while `(|d_0| + ... + |d_m|) / |d_m|` grows like `m`.
fn recur_down(p: Dd, q: Dd, c: Dd, z: f64, n: f64, y_n: Scaled, y_n1: Scaled) -> Option<Scaled> {
    let n = n as u32;
    let finite_nonzero = |v: f64| (v != 0.0 && v.is_finite()).then_some(v);
    // y_(k+1) and y_k, scaled by 2^-e
    let (mut y1, mut y0, mut e) = (ldexp_f64(y_n1.m, y_n1.e - y_n.e), y_n.m, y_n.e);
    let (q_abs, z_abs) = (q.hi.abs(), z.abs());
    // `t_(k-1)` and `d_(k-1) / d_k`, for the error bound
    let mut steps = Vec::with_capacity(n as usize);
    for k in (1..=n).rev().map(f64::from) {
        // (C - 1, C - p and C - q are computed from c, since they can be tiny)
        let (cc, cc1) = (add_hi(c, k), add_hi(c, k - 1.0));
        let cc_minus = |v: Dd| (-v + c.hi + k + c.lo).hi;
        let (cc_p, cc_q) = (cc_minus(p), cc_minus(q));
        let coef_a = cc * cc1 * (z - 1.0);
        let coef_b = cc * (cc1 - (cc_p + cc1 - q.hi) * z);
        let coef_c = cc_p * cc_q * z;
        let (t0, t1) = (coef_b * y0, coef_c * y1);
        let ym = finite_nonzero(-(t0 + t1) / coef_a)?;
        // The rounding errors of the step: of t0, t1 and their sum (6), of the parts of `coef_b`
        // (5), whose absolute values sum to `|b_abs|`, and of `coef_a` and the division (6).
        let cc1_abs = cc1.abs();
        let b_abs = cc * (cc1_abs + (cc_p.abs() + cc1_abs + q_abs) * z_abs);
        let t_m = MACHEP * (6.0 * (t0.abs() + t1.abs()) + 5.0 * (b_abs * y0).abs())
            / (coef_a * ym).abs()
            + 6.0 * MACHEP;
        // d_(k-1) / d_k, failing if it's 0 (if C_k = 0, where the Casoratian vanishes)
        let rho = finite_nonzero(coef_c / coef_a * (y1 / ym))?;
        steps.push((t_m, rho));
        (y1, y0) = (y0, ym);
        if let Some(s) = rescale_exp(y0) {
            (y0, y1) = (ldexp_f64(y0, -s), ldexp_f64(y1, -s));
            e += s;
        }
    }
    let (u, err) = steps
        .into_iter()
        .rfold((1.0_f64, 0.0), |(u, err), (t_m, rho)| {
            (1.0 + u * rho, err + t_m * u.abs())
        });
    let err = err + (y_n.err * u.abs() + y_n1.err * (1.0 - u).abs());
    Scaled::new(y0, e, err)
}

/// `(1 - z)^p F(p, q; c; z)` for `z <= 0` (steps 2 and 3), with `w = z / (z - 1)`
fn hyp2f1_recur_form(p: Dd, q: Dd, c: Dd, z: f64, w: Dd) -> Option<Scaled> {
    // Give up early if the series takes too many terms: they peak near k = |p| w / (1 - w), then
    // decay like w^k, for about ln(2^53) / (1 - w) ~ 40 / (1 - w) terms.
    if (p.hi.abs() * w.hi + 40.0) / (1.0 - w.hi) > HYP2F1_RECUR_MAXITER {
        return None;
    }
    // the number of series terms plus recurrence steps left
    let mut budget = HYP2F1_RECUR_MAXITER;
    let r = -q + c.hi + c.lo;
    let direct = series_w(p, r, c, w, &mut budget);
    // the smallest n with c + n > max(0, q)
    let n = (q.hi.max(0.0) - c.hi).floor() + 1.0;
    if c.hi > q.hi.max(0.0) || err_of(direct) <= HYP2F1_ETHRESH || budget < n {
        return direct;
    }
    budget -= n;
    let (c_n, r_n) = (c + n, r + n);
    let recur = series_w(p, r_n, c_n, w, &mut budget).and_then(|y_n| {
        let y_n1 = series_w(p, r_n + 1.0, c_n + 1.0, w, &mut budget)?;
        // (the recurrence is symmetric in p and q)
        recur_down(p, q, c, z, n, y_n, y_n1)
    });
    min_err(direct, recur)
}

/// The 4 combinations `(p, q, s)` of step 1 for `x = z < 0`
fn neg_z_combinations(a: Dd, b: Dd, c: Dd) -> Combinations {
    let (ca, cb) = (c - a, c - b);
    [(a, b, -a), (b, a, -b), (ca, cb, -b), (cb, ca, -a)]
}

/// `best_of` the forms of `combs` (steps 1 to 3); `ln_1mx` is `ln(1 - x)`
fn best_combination(mut combs: Combinations, c: Dd, z: f64, w: Dd, ln_1mx: Dd) -> Option<Scaled> {
    // First try those with p >= 0 and small |p|, whose series cancel less and are shorter.
    let key = |(p, ..): &(Dd, Dd, Dd)| (p.hi < 0.0, p.hi.abs());
    combs.sort_by(|u, v| key(u).partial_cmp(&key(v)).unwrap_or(Ordering::Equal));
    best_of(combs.into_iter().map(|(p, q, s)| {
        let y = hyp2f1_recur_form(p, q, c, z, w)?;
        Some(y.with_l(s * ln_1mx))
    }))
}

/// `ln(v)` for `v > 0`, after scaling `v` to `[1, 2)` by a power of 2, since `DoubleDouble::log`
/// fails for `v` beyond about `e^709` (as it computes `exp(-ln(v))`)
fn ln_dd(v: Dd) -> Dd {
    let (_, k) = frexp_f64(v.hi);
    v.ldexp(1 - k).log() + LOG2 * f64::from(k - 1)
}

/// `ln|Gamma(g)|`, the sign of `Gamma(g)`, and a bound of the absolute error of the former, or
/// `None` at (or too close to) a pole:
///
/// - for `g >= 13`, by Stirling's formula in double-double, with the correction of cephes' `lgam`
///   (DLMF 5.11.1)
/// - for `g <= -13`, by the reflection formula `Gamma(g) Gamma(1 - g) = pi / sin(pi g)`
///   (DLMF 5.5.3), to first order in `g.lo`
/// - else by `lgam` at `g.hi`, to first order in `g.lo`
///
/// The second-order terms `f'' g.lo^2 / 2` are bounded between `g.hi` and `g.hi + g.lo`, with `d`
/// the smallest distance to a pole there: for `f = ln Gamma`, `f'' < 1/d^2 + pi^2` if `g < 1/2`,
/// and `f'' < 1/g + 1/g^2` if `g > 0`, and for `f = ln|sin(pi g)|`, `|f''| < (pi / 2d)^2` (since
/// `|sin(pi g)| >= 2d`).
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
    let err = 4.0 * MACHEP * (v.abs() + dpsi.abs() + 1.0) + trigamma_lo(g, d) * g.lo.abs() / 2.0;
    Some((Dd::new(v) + dpsi, f64::from(sign), err))
}

/// `|g.lo|` times the bound of `(ln Gamma)''` in `ln_gamma` (0 if `g.lo = 0`, even if it overflows)
fn trigamma_lo(g: Dd, d: f64) -> f64 {
    let lo = g.lo.abs();
    if lo == 0.0 {
        0.0
    } else if g.hi < 0.5 {
        (d.powi(-2) + PI * PI) * lo
    } else {
        let g = g.hi - lo;
        (1.0 + 1.0 / g) / g * lo
    }
}

/// `y()` times `Gamma(num_0) Gamma(num_1) ... / (Gamma(den_0) Gamma(den_1) ...)`, with the error
/// bounds of their `ln_gamma` added; `y` is only evaluated if no argument is at a pole
fn mul_gamma_ratio(num: &[Dd], den: &[Dd], y: impl FnOnce() -> Option<Scaled>) -> Option<Scaled> {
    let num = num.iter().map(|&g| (g, 1.0));
    let mut gammas = num.chain(den.iter().map(|&g| (g, -1.0)));
    let (ln_g, sign, err) =
        gammas.try_fold((Dd::new(0.0), 1.0, 0.0), |(ln_g, sign, err), (g, power)| {
            let (v, v_sign, v_err) = ln_gamma(g)?;
            Some((ln_g + v * power, sign * v_sign, err + v_err))
        })?;
    let mut y = y()?;
    (y.m, y.l, y.err) = (sign * y.m, y.l + ln_g, y.err + err);
    Some(y)
}

/// The sum of the 1 or 2 `terms`, with an error bound that includes their cancellation
fn add_exp<const N: usize>(terms: [Scaled; N]) -> Option<Scaled> {
    // the terms with `l` incl. `2^e`, scaled by `exp(-lmax)`, with `lmax` the largest `l`
    let terms = terms.map(|y| (y, y.l + LOG2 * f64::from(y.e)));
    let (_, lmax) = terms
        .into_iter()
        .reduce(|u, v| if u.1.hi >= v.1.hi { u } else { v })?;
    let (m, err_abs) = terms.into_iter().fold((0.0, 0.0), |(m, err_abs), (y, l)| {
        let r = l - lmax;
        let t = y.m * (r.hi.exp() * (1.0 + r.lo));
        // (3 eps for `exp(r)` and the products)
        (m + t, err_abs + t.abs() * (y.err + 3.0 * MACHEP))
    });
    // (and 1 eps for the sum)
    Some(Scaled::new(m, 0, err_abs / m.abs() + MACHEP)?.with_l(lmax))
}

/// `F(a, b; c; t)` for `t < 0` by the 1/t transformation (steps 4 and 5, DLMF 15.8.2), where
/// `t = x` if `x < 0`, and `t = x / (x - 1)` if `0 < x < 1`:
///
/// ```text
/// F(a, b; c; t) = G(a, b) (-t)^-a F(a, a - c + 1; a - b + 1; 1/t) + G(b, a) (-t)^-b F(b, ...)
/// G(a, b) = Gamma(c) Gamma(b - a) / (Gamma(b) Gamma(c - a))
/// ```
///
/// It fails if `b - a` is an integer (a pole of `G`), and its error bound is large if `b - a` is
/// close to one.
fn hyp2f1_recip(a: Dd, b: Dd, c: Dd, x: f64) -> Option<Scaled> {
    let one = Dd::new(1.0);
    let (ln_x, ln_1mx) = (ln_dd(Dd::new(x.abs())), ln_dd(one - x));
    // z = 1/t, w = 1 / (1 - t), ln(-t) and ln(1 - z), from x (e.g. w = 1 - x if 0 < x < 1)
    let (z, w, ln_mt, ln_1mz) = if x < 0.0 {
        (1.0 / x, one / (one - x), ln_x, ln_1mx - ln_x)
    } else {
        ((x - 1.0) / x, one - x, ln_x - ln_1mx, -ln_x)
    };
    // the term `G(a, b) (-t)^-a F(a, a - c + 1; a - b + 1; 1/t)`
    let term = |a: Dd, b: Dd| {
        let (p, q, cc) = (a, a - c + 1.0, a - b + 1.0);
        let y = best_combination(neg_z_combinations(p, q, cc), cc, z, w, ln_1mz)?;
        mul_gamma_ratio(&[c, b - a], &[b, c - a], || Some(y.with_l(y.l - a * ln_mt)))
    };
    add_exp([term(a, b)?, term(b, a)?])
}

/// `F(a, b; c; x)` by steps 1 to 5, with `ln_1mx = ln(1 - x)`
fn hyp2f1_recur_dd(a: Dd, b: Dd, c: Dd, x: f64, ln_1mx: Dd) -> Option<Scaled> {
    // `z`, `w`, and the combinations `(p, q, s)` of step 1
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
    let recip = (x < -1.0).then(|| hyp2f1_recip(a, b, c, x)).flatten();
    if err_of(recip) <= HYP2F1_ETHRESH {
        return recip;
    }
    let best = min_err(best_combination(combinations, c, z, w, ln_1mx), recip);
    // If that fails for x >= -1, also try the 1/t transformation (step 5).
    if x >= -1.0 && !(err_of(best) <= HYP2F1_RECUR_TOL) {
        let recip = if x < 0.0 {
            hyp2f1_recip(a, b, c, x)
        } else {
            // F(a, b; c; x) = (1 - x)^-b F(c - a, b; c; x / (x - 1))
            hyp2f1_recip(c - a, b, c, x).map(|y| y.with_l(y.l - b * ln_1mx))
        };
        min_err(best, recip)
    } else {
        best
    }
}

/// `F(a, b; c; x)` for real `x < 1` with `x != 0` (see the module docs), or `None` if the
/// estimated relative error is too large.
pub(crate) fn hyp2f1_recur(a: f64, b: f64, c: f64, x: f64) -> Option<f64> {
    // (beyond 1e15, the coefficients of the recurrence can overflow)
    if !(x < 1.0 && x != 0.0 && [a, b, c].iter().all(|v| (-1e15..=1e15).contains(v))) {
        return None;
    }
    let (a, b, c) = (Dd::new(a), Dd::new(b), Dd::new(c));
    // (`log1p` of `DoubleDouble` is only as accurate as that of `f64`, but `log` is accurate)
    let ln_1mx = (Dd::new(1.0) - x).log();
    let best = hyp2f1_recur_dd(a, b, c, x, ln_1mx);
    let y = best.filter(|y| y.err <= HYP2F1_RECUR_TOL)?;
    // `y = m 2^e exp(l)`, with `exp(l) = 2^k exp(r)` and `|r| <= ln(2) / 2`. If `|k|` is so large
    // that it's clamped, then so is `e + k`, since `|e|` is at most about 1e9 (2^1000 per term),
    // and the result overflows or underflows.
    let k = (y.l.hi / LN_2).round().clamp(-1e12, 1e12);
    let r = y.l - LOG2 * k;
    let e = (f64::from(y.e) + k).clamp(-1e4, 1e4) as i32;
    Some(ldexp_f64(y.m * r.hi.exp() * (1.0 + r.lo), e))
}

#[cfg(test)]
mod tests {
    use super::{Dd, Scaled, hyp2f1_recur, recur_down};

    #[test]
    fn test_recur_down_start_errors() {
        // y_0 = y_2 / 8 (as B_1 = 0), so the error of y_2 carries over, although u_1 = 0
        let y = |m, err| Scaled::new(m, 0, err).unwrap();
        let (p, q, c) = (Dd::new(3.5), Dd::new(1.5), Dd::new(1.0));
        let y_0 = recur_down(p, q, c, -0.5, 1.0, y(0.5, 0.0), y(0.8, 1e-6)).unwrap();
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
