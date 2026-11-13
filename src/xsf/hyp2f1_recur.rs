//! The Gauss hypergeometric function `F(a, b; c; x)` for real `x < 1` and large real parameters,
//! by backward recurrence in `c`. Not a port: `hyp2f1` (in `cephes/hyp2f1.rs`) tries it first for
//! large parameters, where its series and transformations can lose all precision, and for a tiny
//! negative `a` or `b`.
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
//!    The twin combinations `(p, q)` and `(c - q, c - p)` share this series, as `F` is symmetric
//!    in its first two parameters, so the twin reuses its sum (if it succeeded).
//!
//!    But if `p = -n` is a negative integer (with `n <= HYP2F1_RECUR_MAXITER`), its terms can
//!    alternate in sign and cancel, e.g. for `F(-522, -123.5; 533.5; -0.9)`. So unless its sum or
//!    that of `F(-n, q; c; z)` is within `HYP2F1_ETHRESH`, we also recur up in the degree, from
//!    `F(0) = 1` and `F(-1)` to `F(-n)` (DLMF 15.5.11), with the error bound of step 3, for
//!    `F(-n, c - q; c; w)` and the form (DLMF 15.8.7), often more stable if `c < 0`,
//!
//!    `F(-n, q; c; z) = (c - q)_n / (c)_n F(-n, q; q - c - n + 1; 1 - z)`
//!
//!    and return the best of these if its error bound is at most `HYP2F1_RECUR_TOL`.
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
//!    (incl. the cancellation of its two terms) exceeds `HYP2F1_ETHRESH`. If `b - a` is an integer
//!    (a pole of its Gamma factors), we use its limit instead (DLMF 15.8.9), whose series in `w`
//!    has `psi` weights, so it can't be transformed or recurred like the others, and it can cancel
//!    whatever the signs of the parameters (see `recip_log`), e.g. `F(-900.5, -800.5; 50; -40)`.
//!
//! 5. If `-1 <= x < 0.9` and all error bounds so far exceed `HYP2F1_RECUR_TOL`, also try the same
//!    transformation of `F(a, b; c; x)` if `x < 0` (i.e. DLMF 15.8.3), or of the Pfaff form
//!    `(1 - x)^-b F(c - a, b; c; x / (x - 1))` if `0 < x < 1` (i.e. the 1 - x transformation,
//!    DLMF 15.8.4). Its series are in `w = 1 / (1 - x)` or `w = 1 - x`. The 4 combinations can
//!    all cancel if `a` and `b` are large with opposite signs and `0 < x < 1` (or with the same
//!    sign and `x < 0`), e.g. `F(-975.2, 241.4; 0.00135; 0.39)`, while these two terms rarely do.
//!    For `x >= 0.9`, where `w = x` is close to 1, we try it first instead, as in step 4. Its log
//!    case (step 4) applies here too: if `b - a` is an integer for `x < 0`, or `c - a - b` for
//!    `0 < x < 1` (DLMF 15.8.10).
//!
//! 6. If all error bounds still exceed `HYP2F1_RECUR_TOL`, recur down in `c` as in step 3, but
//!    with other start values. Those of step 3 are sums of the same series, and can lose precision
//!    for every `n`, e.g. if `p` is large and negative, while other forms at `c + n` don't. Here
//!    `F(p, q; c + k; z)` is a multiple of (DLMF 15.8.1)
//!    - `F(a, b; c + k; x)` for `(p, q, z) = (a, b, x)`
//!    - `(1 - x)^k F(a + k, b + k; c + k; x)` for `(c - a, c - b, x)`
//!    - `F(a, b + k; c + k; x)` for `(a, c - b, x / (x - 1))`
//!    - `F(a + k, b; c + k; x)` for `(c - a, b, x / (x - 1))`
//!
//!    So for each of these, we compute `y_n` and `y_(n+1)` (with `n` as in step 3, but at least 1)
//!    by steps 1 to 5, and recur down to `k = 0`, e.g. `F(758.1, 7.7; -151.2; -0.3)` from
//!    `F(758.1, 159.7; 0.8; -0.3)` and `F(758.1, 160.7; 1.8; -0.3)`. If `z > 0`, the recurrence
//!    need not be stable, but its error bound still holds. We skip those with `c > max(0, q)` (so
//!    `n = 1`) unless steps 1 to 5 got within `HYP2F1_RECUR_SHIFT_TOL`: their start values at
//!    `c + 1` and `c + 2` hardly differ from `F(a, b; c; x)`, and in our tests improved the bound
//!    at most 15-fold.
//!
//! Finally, return the result with the smallest error bound, which is too large if it exceeds
//! `HYP2F1_RECUR_TOL`.
//!
//! The parameters are double-double, so that e.g. `C - q + k` is accurate if `q = c - b` and `c`
//! is tiny, and so are `s ln(1 - x)` and the `ln Gamma` of the transformations. The values are
//! kept as `m 2^e`, to avoid over- and underflow. The cost is at most `HYP2F1_RECUR_MAXITER`
//! series terms and recurrence steps per combination, which is skipped (or only its step 3) if it
//! would exceed that: about `(|p| w + 40) / (1 - w)` terms (but at least
//! `(w (c - q) - c) / (1 - w)`), plus `max(0, q) - c` steps. Step 6 costs up to 8 times as much as
//! steps 1 to 5, but only if those fail.
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
use crate::xsf::cephes::{HYP2F1_ETHRESH, HYP2F1_RECUR_NEAR1, lgam_sgn, lgam_stirling_corr, psi};
use crate::xsf::trig::{cospi, sinpi};

/// The largest accepted error bound. These bounds are pessimistic (the actual errors are usually
/// 10 to 300 times smaller), and the methods of the caller are often far worse if this fails.
pub(crate) const HYP2F1_RECUR_TOL: f64 = 1e-10;

/// The largest error bound of steps 1 to 5 for which step 6 also tries the forms with
/// `c > max(0, q)` (empirical)
const HYP2F1_RECUR_SHIFT_TOL: f64 = 100.0 * HYP2F1_RECUR_TOL;

/// The maximum number of series terms plus recurrence steps, per combination
const HYP2F1_RECUR_MAXITER: f64 = 1e6;

/// The combinations `(p, q, s)` of step 1
type Combinations = [(Dd, Dd, Dd); 4];

/// `(A_k, B_k, C_k, b_abs)` of `recur_down`
type Coefs = (f64, f64, f64, f64);

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
///
/// With `weights = Some((beta_0, err))` (step 4), it sums `t_k beta_k` instead, with `err` bounding
/// the error of `beta_0`, and
///
/// `beta_(k+1) = beta_k + 1/(k + 1) + 1/(c + k) - 1/(p + k) - 1/(r + k)`
///
/// It sums through at most `t_kmax`. It's inlined, since the weights make the unweighted series
/// about 5% slower otherwise.
#[allow(clippy::inline_always)]
#[inline(always)]
fn series_weighted(
    p: Dd,
    r: Dd,
    c: Dd,
    w: Dd,
    budget: &mut f64,
    weights: Option<(f64, f64)>,
    kmax: f64,
) -> Option<Scaled> {
    let (w, w_lo) = (w.hi, w.lo);
    // from here on, the factors of the term ratios have constant signs
    let kmin = [-p.hi, -r.hi, -c.hi].into_iter().fold(0.0, f64::max).ceil();
    // `beta_k` and a bound of its absolute error, the sum of `|t_k|` times that bound, and the
    // sum of `|parts|` of `beta_k - beta_(k-1)`
    let (mut beta, mut beta_err) = weights.unwrap_or((1.0, 0.0));
    let (mut sum_err, mut d_abs) = (beta_err, 0.0);
    // The terms are `t 2^te`, and the sums are scaled by `2^-e`, since the terms can underflow
    // before `k = -c`, and be huge after it.
    let (mut sum, mut sum_abs, mut k_sum, mut k_sum_abs) = (beta, beta.abs(), 0.0, 0.0_f64);
    let (mut t, mut te, mut e, mut k) = (1.0, 0, 0, 0.0);
    // the factors `(p + k) / (k + 1)` and `(r + k) / (c + k)` of the next term ratio, and `p + k`,
    // `r + k` and `c + k`
    let (mut pk, mut rk, mut ck) = (p.hi, r.hi, c.hi);
    let (mut f0, mut f1) = (pk, rk / ck);
    let tame = 2.0_f64.powi(-250)..=2.0_f64.powi(250);
    while k < kmax {
        let fw = f1 * w;
        if [f0, fw].into_iter().all(|f| tame.contains(&f.abs())) {
            // (usually: then `|t f0 fw|` is 0 or in `[2^-1000, 2^1000]`, and one rescale equals two
            // `mul_scaled`, as scaling by powers of 2 is exact)
            (t, te) = rescaled(t * f0 * fw, te);
        } else {
            (t, te) = mul_scaled((t, te), f0, pk)?;
            (t, te) = mul_scaled((t, te), fw, rk)?;
        }
        let term = if te == 0 { t } else { ldexp_f64(t, te) };
        let term_abs = term.abs();
        k += 1.0;
        if weights.is_some() {
            let parts = [1.0 / k, 1.0 / ck, -1.0 / pk, -1.0 / rk];
            d_abs = parts.into_iter().map(f64::abs).sum();
            beta += parts.into_iter().sum::<f64>();
            // (2 roundings per part, and 4 additions)
            beta_err += MACHEP * (5.0 * d_abs + beta.abs());
            sum_err += term_abs * beta_err;
        }
        // (the rounding of `term * beta` is within the `(k + 1) eps` of `sum_abs` below)
        let term_beta = term * beta;
        let term_beta_abs = term_beta.abs();
        sum += term_beta;
        sum_abs += term_beta_abs;
        k_sum += k * term_beta;
        k_sum_abs += k * term_beta_abs;
        if !k_sum_abs.is_finite() {
            return None;
        }
        (pk, rk, ck) = (add_hi(p, k), add_hi(r, k), add_hi(c, k));
        (f0, f1) = (pk / (k + 1.0), rk / ck);
        // For j >= k > kmin, both factors of the term ratio are positive and tend monotonically
        // to 1, so the ratio is at most `rho`, which bounds the rest of the series geometrically.
        // With weights, the parts of `beta_j - beta_(j-1)` decrease as well, so
        // `|beta_j| <= |beta_k| + (j - k) d_abs`.
        if k > kmin {
            let rho = w * f0.max(1.0) * f1.max(1.0);
            let beta_max = beta.abs() + d_abs / (1.0 - rho);
            if rho < 1.0 && term_abs * beta_max * rho / (1.0 - rho) <= sum.abs() * MACHEP / 2.0 {
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
            // (separately, since a 5-element array makes the series about 4% slower)
            sum_err = ldexp_f64(sum_err, -s);
            (e, te) = (e + s, te - s);
        }
    }
    // dF/dw = (sum of k t_k) / w
    if w != 0.0 {
        sum += k_sum * (w_lo / w);
    }
    let err = (MACHEP * (8.0 * k_sum_abs + (k + 1.0) * sum_abs) + sum_err) / sum.abs();
    Scaled::new(sum, e, err)
}

/// `series_weighted` without weights
fn series_w(p: Dd, r: Dd, c: Dd, w: Dd, budget: &mut f64) -> Option<Scaled> {
    series_weighted(p, r, c, w, budget, None, f64::INFINITY)
}

/// `y_0` from `y_n` and `y_(n+1)` by the recurrence `A_k y_(k-1) + B_k y_k + C_k y_(k+1) = 0`,
/// where `coefs(k)` is `(A_k, B_k, C_k, b_abs)`, with `|b_abs|` the sum of `|parts|` of `B_k` (at
/// most 5 roundings in `A_k` and `B_k`, and 4 in `C_k`).
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
fn recur_down(n: f64, y_n: Scaled, y_n1: Scaled, coefs: impl Fn(f64) -> Coefs) -> Option<Scaled> {
    let n = n as u32;
    let finite_nonzero = |v: f64| (v != 0.0 && v.is_finite()).then_some(v);
    // y_(k+1) and y_k, scaled by 2^-e
    let (mut y1, mut y0, mut e) = (ldexp_f64(y_n1.m, y_n1.e - y_n.e), y_n.m, y_n.e);
    // `t_(k-1)` and `d_(k-1) / d_k`, for the error bound
    let mut steps = Vec::with_capacity(n as usize);
    for k in (1..=n).rev().map(f64::from) {
        let (coef_a, coef_b, coef_c, b_abs) = coefs(k);
        let (t0, t1) = (coef_b * y0, coef_c * y1);
        let ym = finite_nonzero(-(t0 + t1) / coef_a)?;
        // The rounding errors of the step: of t0, t1 (incl. `coef_c`) and their sum (6), of the
        // parts of `coef_b` (5), and of `coef_a` and the division (6).
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

/// The coefficients of `recur_down` for `y_k = F(p, q; c + k; z)` (DLMF 15.5.18, `C = c + k`):
///
/// ```text
/// A_k = C (C - 1) (z - 1)
/// B_k = C (C - 1 - (2 C - p - q - 1) z)
/// C_k = (C - p) (C - q) z
/// ```
fn coefs_c(p: Dd, q: Dd, c: Dd, z: Dd) -> impl Fn(f64) -> Coefs {
    // (z - 1 is computed from z, since it can be tiny)
    let (z, z1) = (z.hi, add_hi(z, -1.0));
    let (q_abs, z_abs) = (q.hi.abs(), z.abs());
    move |k| {
        // (C - 1, C - p and C - q are computed from c, since they can be tiny)
        let (cc, cc1) = (add_hi(c, k), add_hi(c, k - 1.0));
        let cc_minus = |v: Dd| (-v + c.hi + k + c.lo).hi;
        let (cc_p, cc_q) = (cc_minus(p), cc_minus(q));
        let coef_b = cc * (cc1 - (cc_p + cc1 - q.hi) * z);
        let cc1_abs = cc1.abs();
        let b_abs = cc * (cc1_abs + (cc_p.abs() + cc1_abs + q_abs) * z_abs);
        (cc * cc1 * z1, coef_b, cc_p * cc_q * z, b_abs)
    }
}

/// The coefficients of `recur_down` for `y_k = F(k - n, r; c; t)` (DLMF 15.5.11, `a = k - n`):
///
/// ```text
/// A_k = c - a
/// B_k = 2 a - c + (r - a) t
/// C_k = a (t - 1)
/// ```
fn coefs_a(n: f64, r: Dd, c: Dd, t: Dd) -> impl Fn(f64) -> Coefs {
    // (t - 1 is computed from t, since it can be tiny)
    let t1 = add_hi(t, -1.0);
    move |k| {
        let a = k - n;
        let (c_2a, r_a_t) = (add_hi(c, -2.0 * a), add_hi(r, -a) * t.hi);
        let b_abs = c_2a.abs() + r_a_t.abs();
        (add_hi(c, -a), r_a_t - c_2a, a * t1, b_abs)
    }
}

/// `F(-n, r; c; t)` for `n >= 1` (step 2), by recurring up from `F(0) = 1` and
/// `F(-1) = 1 - r t / c`, i.e. down in `k` for `y_k = F(k - n, r; c; t)`
fn poly_up(n: f64, r: Dd, c: Dd, t: Dd) -> Option<Scaled> {
    // (r / c first, since r t can underflow if c is tiny, and fail if r / c is subnormal)
    let r_c = Some(r.hi / c.hi).filter(|&r_c| r_c == 0.0 || r_c.is_normal())?;
    let rtc = r_c * t.hi;
    // (the roundings of r, t, c, the quotient and the product (5), and of the difference)
    let err = MACHEP * (5.0 * rtc.abs() + 1.0) / (1.0 - rtc).abs();
    let (y_1, y_0) = (Scaled::new(1.0 - rtc, 0, err)?, Scaled::new(1.0, 0, 0.0)?);
    recur_down(n - 1.0, y_1, y_0, coefs_a(n, r, c, t))
}

/// The polynomial forms of step 2 of `combs` with `p = -n`; `ln_1mx` is `ln(1 - x)`
fn hyp2f1_poly(combs: Combinations, c: Dd, x: f64, w: Dd, ln_1mx: Dd) -> Option<Scaled> {
    // 1 - z, which is 1 - x if x < 0, and 1 / (1 - x) else, and its log
    let one = Dd::new(1.0);
    let (one_mz, ln_1mz) = if x < 0.0 {
        (one - x, ln_1mx)
    } else {
        (1.0 / (one - x), -ln_1mx)
    };
    let is_poly = |p: Dd| p.is_integer() && (1.0..=HYP2F1_RECUR_MAXITER).contains(&-p.hi);
    combs
        .into_iter()
        .filter(|&(p, ..)| is_poly(p))
        .map(|(p, q, s)| {
            let (n, r, l) = (-p.hi, c - q, s * ln_1mx);
            // (with the factor `(1 - z)^p` of the forms in `z` and `1 - z`)
            let l_z = l + p * ln_1mz;
            // The sums of `F(p, q; c; z)` and `F(p, c - q; c; w)` come first: if the terms of one
            // have one sign, it is more accurate than the recurrences.
            let sum = |q: Dd, t: Dd, l: Dd| {
                let mut budget = HYP2F1_RECUR_MAXITER;
                Some(series_weighted(p, q, c, t, &mut budget, None, n)?.with_l(l))
            };
            let direct = min_err(sum(q, one - one_mz, l_z), sum(r, w, l));
            if err_of(direct) <= HYP2F1_ETHRESH {
                return direct;
            }
            // `(1 - z)^p F(p, q; c; z) = F(p, c - q; c; w)`, as in step 2
            let y_w = poly_up(n, r, c, w).map(|y| y.with_l(l));
            // `(c - q)_n / (c)_n = Gamma(c - q + n) Gamma(c) / (Gamma(c - q) Gamma(c + n))`
            let y_1mz = mul_gamma_ratio(&[r + n, c], &[r, c + n], || {
                let y = poly_up(n, q, q - c + (1.0 - n), one_mz)?;
                Some(y.with_l(l_z))
            });
            min_err(direct, min_err(y_w, y_1mz))
        })
        .fold(None, min_err)
}

/// The smallest `n >= 1` with `c + n > max(0, q)`
fn recur_start(q: Dd, c: Dd) -> f64 {
    (q.hi.max(0.0) - c.hi).max(0.0).floor() + 1.0
}

/// Whether the series `F(p, r; c; w)` likely takes more than `HYP2F1_RECUR_MAXITER` terms
fn series_too_long(p: Dd, r: Dd, c: Dd, w: Dd) -> bool {
    // The terms peak near k = |p| w / (1 - w), then decay like w^k, for about
    // ln(2^53) / (1 - w) ~ 40 / (1 - w) terms. It only stops once w (r + k) / (c + k) < 1, i.e.
    // after k = (w r - c) / (1 - w).
    (p.hi.abs() * w.hi + 40.0).max(w.hi * r.hi - c.hi) / (1.0 - w.hi) > HYP2F1_RECUR_MAXITER
}

/// `(1 - z)^p F(p, q; c; z)` for `z <= 0` (steps 2 and 3), with `w = z / (z - 1)` and `r = c - q`.
/// `shared` holds the step-2 sum and the budget left after it, for the twin `(c - q, c - p)`.
fn hyp2f1_recur_form(
    p: Dd,
    q: Dd,
    r: Dd,
    c: Dd,
    z: f64,
    w: Dd,
    shared: &mut Option<(Scaled, f64)>,
) -> Option<Scaled> {
    if series_too_long(p, r, c, w) {
        return None;
    }
    // the step-2 sum, and the number of series terms plus recurrence steps left
    let (direct, mut budget) = if let Some((y, budget)) = *shared {
        (Some(y), budget)
    } else {
        let mut budget = HYP2F1_RECUR_MAXITER;
        (series_w(p, r, c, w, &mut budget), budget)
    };
    // (if it failed, the twin sums it again: e.g. for subnormal `p`, `p / (k + 1)` underflows,
    // but `p / (c + k)` may not)
    *shared = direct.map(|y| (y, budget));
    let n = recur_start(q, c);
    // Give up if too little is left for the n steps and the series at c + n and c + n + 1, each
    // estimated at n terms fewer than this one (without its first -c terms).
    let terms = (HYP2F1_RECUR_MAXITER - budget - n).max(0.0);
    if c.hi > q.hi.max(0.0) || err_of(direct) <= HYP2F1_ETHRESH || budget - n < 2.0 * terms {
        return direct;
    }
    budget -= n;
    let (c_n, r_n) = (c + n, r + n);
    let recur = series_w(p, r_n, c_n, w, &mut budget).and_then(|y_n| {
        let y_n1 = series_w(p, r_n + 1.0, c_n + 1.0, w, &mut budget)?;
        // (the recurrence is symmetric in p and q)
        recur_down(n, y_n, y_n1, coefs_c(p, q, c, Dd::new(z)))
    });
    min_err(direct, recur)
}

/// The 4 combinations `(p, q, s)` of step 1 for `x = z < 0`
fn neg_z_combinations(a: Dd, b: Dd, c: Dd) -> Combinations {
    let (ca, cb) = (c - a, c - b);
    [(a, b, -a), (b, a, -b), (ca, cb, -b), (cb, ca, -a)]
}

/// `best_of` the forms of `combs` (steps 1 to 3); `ln_1mx` is `ln(1 - x)`. `combs[i]` and
/// `combs[3 - i]` must be twins, `(p, q)` and `(c - q, c - p)`.
fn best_combination(combs: Combinations, c: Dd, z: f64, w: Dd, ln_1mx: Dd) -> Option<Scaled> {
    // First try those with p >= 0 and small |p|, whose series cancel less and are shorter.
    let key = |i: usize| (combs[i].0.hi < 0.0, combs[i].0.hi.abs());
    let mut order = [0, 1, 2, 3];
    order.sort_by(|&i, &j| key(i).partial_cmp(&key(j)).unwrap_or(Ordering::Equal));
    // the step-2 sums of the twins
    let mut shared = [None; 2];
    best_of(order.into_iter().map(|i| {
        let (p, q, s) = combs[i];
        // `c - q` is the `p` of the twin, while `c - (c - p)` would round a tiny `p`
        let r = combs[3 - i].0;
        let y = hyp2f1_recur_form(p, q, r, c, z, w, &mut shared[i.min(3 - i)])?;
        Some(y.with_l(s * ln_1mx))
    }))
}

/// `ln(v)` for `v > 0`, after scaling `v` to `[1, 2)` by a power of 2, since `DoubleDouble::log`
/// fails for `v` beyond about `e^709` (as it computes `exp(-ln(v))`)
fn ln_dd(v: Dd) -> Dd {
    let (_, k) = frexp_f64(v.hi);
    v.ldexp(1 - k).log() + LOG2 * f64::from(k - 1)
}

/// A lower bound `d` of the distance of `g.hi + [-|g.lo|, |g.lo|]` to the integers (negative if it
/// contains one), or `None` if `g < 1/2` and `d <= 0`, i.e. at (or too close to) a pole of `Gamma`
fn pole_distance(g: Dd) -> Option<f64> {
    let d = (g.hi - g.hi.round()).abs() - g.lo.abs();
    (!(g.hi < 0.5 && d <= 0.0)).then_some(d)
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
    let d = pole_distance(g)?;
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

/// `psi(v)` at `v.hi`, and a bound of its absolute error, or `None` at (or too close to) a pole.
/// For `v < 0`, `psi` uses `psi(v) = psi(1 - v) - pi cot(pi v)`, with `|pi cot(pi v)| <= 1/d`
/// (`d` as in `pole_distance`), and its error is relative to both terms.
fn psi_dd(v: Dd) -> Option<(f64, f64)> {
    let d = pole_distance(v)?;
    let y = psi(v.hi);
    let cot = if v.hi < 0.0 { 1.0 / d } else { 0.0 };
    let err = 8.0 * MACHEP * (y.abs() + 2.0 * cot + 1.0) + trigamma_lo(v, d);
    Some((y, err))
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

/// `hyp2f1_recip` if `m = b - a` is a non-negative integer, by the limit `b - a -> m` of its
/// transformation (DLMF 15.8.9), with `w = 1 / (1 - t)`:
///
/// ```text
/// F(a, b; c; t) / Gamma(c) = (1 - t)^-a Gamma(m) / (Gamma(b) Gamma(c - a)) S_1
///     + (-1)^m (1 - t)^-b / (Gamma(a) Gamma(c - b) m!) S_2
/// S_1 = sum_(k < m) (a)_k (c - b)_k / ((1 - m)_k k!) w^k
/// S_2 = sum_k (b)_k (c - a)_k / ((m + 1)_k k!) w^k beta_k
/// beta_k = ln(1 - t) + psi(k + 1) + psi(m + k + 1) - psi(b + k) - psi(c - a + k)
/// ```
///
/// `S_2` is the series of step 2 of the second term of `hyp2f1_recip`, but with the weights
/// `beta_k` of `series_weighted`. As `ln(1 - t) = -ln(w)`, its terms are
///
/// ```text
/// t_k beta_k = -g'(k) / g(0)
/// g(s) = Gamma(b + s) Gamma(c - a + s) / (Gamma(1 + s) Gamma(m + 1 + s)) w^s
/// ```
///
/// so they change sign where `|g|` peaks, and `S_2` cancels if that peak is far above `|g(0)|`
/// (e.g. if `b (c - a) w / (m + 1)` is large), whatever the signs of `b` and `c - a`.
fn recip_log(a: Dd, b: Dd, c: Dd, w: Dd, ln_1mt: Dd) -> Option<Scaled> {
    let (m, ca, cb) = ((b - a).hi, c - a, c - b);
    let m1 = Dd::new(m + 1.0);
    if series_too_long(b, ca, m1, w) || m > HYP2F1_RECUR_MAXITER {
        return None;
    }
    // beta_0 and a bound of its absolute error
    let (beta, beta_err) = [(Dd::new(1.0), 1.0), (m1, 1.0), (b, -1.0), (ca, -1.0)]
        .into_iter()
        .try_fold((ln_1mt, 0.0), |(beta, beta_err), (v, sign)| {
            let (y, err) = psi_dd(v)?;
            Some((beta + sign * y, beta_err + err))
        })?;
    let weights = Some((beta.hi, beta_err + MACHEP * beta.hi.abs()));
    let mut budget = HYP2F1_RECUR_MAXITER;
    let term2 = mul_gamma_ratio(&[c], &[a, cb, m1], || {
        let mut s2 = series_weighted(b, ca, m1, w, &mut budget, weights, f64::INFINITY)?;
        // (-1)^m
        s2.m *= if m % 2.0 == 0.0 { 1.0 } else { -1.0 };
        Some(s2.with_l(-b * ln_1mt))
    })?;
    if m == 0.0 {
        // (S_1 is an empty sum)
        return add_exp([term2]);
    }
    let s1 = series_weighted(a, cb, Dd::new(1.0 - m), w, &mut budget, None, m - 1.0)?;
    let term1 = mul_gamma_ratio(&[c, Dd::new(m)], &[b, ca], || Some(s1.with_l(-a * ln_1mt)))?;
    add_exp([term1, term2])
}

/// `F(a, b; c; t)` for `t < 0` by the 1/t transformation (steps 4 and 5, DLMF 15.8.2), where
/// `t = x` if `x < 0`, and `t = x / (x - 1)` if `0 < x < 1`:
///
/// ```text
/// F(a, b; c; t) = G(a, b) (-t)^-a F(a, a - c + 1; a - b + 1; 1/t) + G(b, a) (-t)^-b F(b, ...)
/// G(a, b) = Gamma(c) Gamma(b - a) / (Gamma(b) Gamma(c - a))
/// ```
///
/// If `b - a` is an integer (a pole of `G`; in step 5 with `0 < x < 1`, that's an integer
/// `c - a - b`), it uses `recip_log` instead. Its error bound is large if `b - a` is close to an
/// integer. `ln_1mx` is `ln(1 - x)`.
fn hyp2f1_recip(a: Dd, b: Dd, c: Dd, x: f64, ln_1mx: Dd) -> Option<Scaled> {
    let one = Dd::new(1.0);
    let ln_x = ln_dd(Dd::new(x.abs()));
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
    let m = b - a;
    if m.is_integer() {
        let (a, b) = if m.hi < 0.0 { (b, a) } else { (a, b) };
        // ln(-t) + ln(1 - 1/t) = ln(1 - t)
        recip_log(a, b, c, w, ln_mt + ln_1mz)
    } else {
        add_exp([term(a, b)?, term(b, a)?])
    }
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
    // the polynomial cases of step 2
    let poly = hyp2f1_poly(combinations, c, x, w, ln_1mx);
    if err_of(poly) <= HYP2F1_RECUR_TOL {
        return poly;
    }
    // the 1/x transformation (step 4), or the 1/t transformation (step 5)
    let try_recip = || {
        if x < 0.0 {
            hyp2f1_recip(a, b, c, x, ln_1mx)
        } else {
            // F(a, b; c; x) = (1 - x)^-b F(c - a, b; c; x / (x - 1))
            hyp2f1_recip(c - a, b, c, x, ln_1mx).map(|y| y.with_l(y.l - b * ln_1mx))
        }
    };
    // For x < -1 or x >= 0.9, first try that, as its series are much shorter if w is near 1.
    let recip_first = !(-1.0..HYP2F1_RECUR_NEAR1).contains(&x);
    let recip = recip_first.then(try_recip).flatten();
    if err_of(recip) <= HYP2F1_ETHRESH {
        return recip;
    }
    let best = min_err(best_combination(combinations, c, z, w, ln_1mx), recip);
    // Else, if that fails, try it last (step 5).
    if !recip_first && !(err_of(best) <= HYP2F1_RECUR_TOL) {
        min_err(best, try_recip())
    } else {
        best
    }
}

/// `y` with `exp(l)` folded into `m 2^e` (so `l = 0`), and its roundings into the error bound
fn fold_exp(y: Scaled) -> Scaled {
    // `exp(l) = 2^k exp(r)` with `|r| <= ln(2) / 2`. If `|k|` is so large that it's clamped, then
    // so is `e + k`, since `|e|` is at most about 1e9 (2^1000 per term), and the value overflows
    // or underflows.
    let k = (y.l.hi / LN_2).round().clamp(-1e12, 1e12);
    let r = y.l - LOG2 * k;
    let e = (f64::from(y.e) + k).clamp(-1e9, 1e9) as i32;
    let (m, l) = (y.m * r.hi.exp() * (1.0 + r.lo), Dd::new(0.0));
    // (the roundings of `exp`, `1 + r.lo`, and the 2 products)
    let err = y.err + 4.0 * MACHEP;
    Scaled { m, e, l, err }
}

/// Step 6: `F(a, b; c; x)` by recurring down in `c` from start values by `hyp2f1_recur_dd`, with
/// `ln_1mx = ln(1 - x)`, where `err` is the error bound of steps 1 to 5
fn hyp2f1_recur_shifted(a: Dd, b: Dd, c: Dd, x: f64, ln_1mx: Dd, err: f64) -> Option<Scaled> {
    let x_dd = Dd::new(x);
    let z_pfaff = x_dd / (x_dd - 1.0);
    // `(p, q, z)` of `F(p, q; c + k; z)`, and the direction `(i, j)` of `(a, b)` along which it's a
    // multiple of `F(a + i k, b + j k; c + k; x) (1 - x)^(i j k)` (DLMF 15.8.1)
    let forms = [
        (a, b, x_dd, 0.0, 0.0),
        (c - a, c - b, x_dd, 1.0, 1.0),
        (a, c - b, z_pfaff, 0.0, 1.0),
        (c - a, b, z_pfaff, 1.0, 0.0),
    ];
    best_of(forms.into_iter().map(|(p, q, z, i, j)| {
        // (if c > max(0, q), only if steps 1 to 5 nearly succeeded, see the module docs)
        let n = recur_start(q, c);
        if c.hi > q.hi.max(0.0) && err > HYP2F1_RECUR_SHIFT_TOL || n > HYP2F1_RECUR_MAXITER {
            return None;
        }
        // (only the start values that we would return)
        let start = |k: f64| {
            let y = hyp2f1_recur_dd(a + i * k, b + j * k, c + k, x, ln_1mx);
            let y = y.filter(|y| y.err <= HYP2F1_RECUR_TOL)?;
            Some(y.with_l(y.l + ln_1mx * (i * j * k)))
        };
        let (y_n, y_n1) = (start(n)?, start(n + 1.0)?);
        // (with `y_(n+1)` relative to `exp(y_n.l)`)
        let y_n1 = fold_exp(y_n1.with_l(y_n1.l - y_n.l));
        Some(recur_down(n, y_n, y_n1, coefs_c(p, q, c, z))?.with_l(y_n.l))
    }))
}

/// `F(a, b; c; x)` for real `x < 1` with `x != 0` (see the module docs) and its error bound, which
/// is too large if it exceeds `HYP2F1_RECUR_TOL`.
pub(crate) fn hyp2f1_recur(a: f64, b: f64, c: f64, x: f64) -> Option<(f64, f64)> {
    // (beyond 1e15, the coefficients of the recurrence can overflow)
    if !(x < 1.0 && x != 0.0 && [a, b, c].iter().all(|v| (-1e15..=1e15).contains(v))) {
        return None;
    }
    let (a, b, c) = (Dd::new(a), Dd::new(b), Dd::new(c));
    // (`log1p` of `DoubleDouble` is only as accurate as that of `f64`, but `ln_dd` is accurate)
    let ln_1mx = ln_dd(Dd::new(1.0) - x);
    let mut best = hyp2f1_recur_dd(a, b, c, x, ln_1mx);
    let err = err_of(best);
    if !(err <= HYP2F1_RECUR_TOL) {
        best = min_err(best, hyp2f1_recur_shifted(a, b, c, x, ln_1mx, err));
    }
    let y = fold_exp(best?);
    Some((ldexp_f64(y.m, y.e), y.err))
}

#[cfg(test)]
mod tests {
    use super::{Dd, HYP2F1_RECUR_TOL, Scaled, coefs_c, hyp2f1_recur, recur_down};

    #[test]
    fn test_recur_down_start_errors() {
        // y_0 = y_2 / 8 (as B_1 = 0), so the error of y_2 carries over, although u_1 = 0
        let y = |m, err| Scaled::new(m, 0, err).unwrap();
        let (p, q, c, z) = (Dd::new(3.5), Dd::new(1.5), Dd::new(1.0), Dd::new(-0.5));
        let y_0 = recur_down(1.0, y(0.5, 0.0), y(0.8, 1e-6), coefs_c(p, q, c, z)).unwrap();
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
            let actual = hyp2f1_recur(a, b, c, x).filter(|&(_, err)| err <= HYP2F1_RECUR_TOL);
            assert!(
                actual.is_none_or(|(y, _)| (y / desired - 1.0).abs() < 1e-10),
                "{actual:?} != {desired}"
            );
        }
    }
}
