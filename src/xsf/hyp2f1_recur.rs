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
//!    if `p >= 0`, or if `p < 0` is small.
//!
//! 3. So for each combination, we sum this series at `C = c`. If `c <= max(0, q)`, and that loses
//!    too much precision, we also sum it at `C = c + n` and `C = c + n + 1`, with the smallest `n`
//!    for which `c + n > max(0, q)`, and recur down to `C = c` with the contiguous relation in `c`
//!    (DLMF 15.5.18). For `Re z < 1/2`, `F(p, q; c + n; z)` is the minimal solution of this
//!    recurrence as `n -> inf` (Gil, Segura and Temme 2007, sec. 3.2), so recurring towards
//!    smaller `c` is the stable direction. But this is only asymptotic, for `c + n` much larger
//!    than `|p|` and `|q|`, and the recurrence can be unstable in between (Deano and Segura 2007).
//!    So instead of predicting this, we bound the error of the result to first order, see
//!    `recur_down`.
//!
//! 4. Return the result with the smallest error bound, or `None` if that exceeds
//!    `HYP2F1_RECUR_TOL`, for the caller to fall back to another method.
//!
//! Some details:
//!
//! - The parameters are double-double, so that e.g. `C - q + k` is accurate if `q = c - b` and
//!   `c` is tiny, and `(1 - x)^s` is computed in double-double, since `|s ln(1 - x)|` can be large.
//! - The series and the recurrence keep their values as `m * 2^e`, to avoid overflow and
//!   underflow.
//! - The cost is bounded by `HYP2F1_RECUR_MAXITER` series terms and recurrence steps per
//!   combination. It's about `(|p| w + 40) / (1 - w)` terms, plus `max(0, q) - c` steps.
//!
//! # References
//!
//! - A. Gil, J. Segura and N.M. Temme, "Numerically satisfactory solutions of hypergeometric
//!   recursions", Math. Comp. 76 (2007), 1449-1468. <https://doi.org/10.1090/S0025-5718-07-01918-7>
//! - A. Gil, J. Segura and N.M. Temme, "The ABC of hyper recursions", J. Comput. Appl. Math. 190
//!   (2006), 270-286. <https://doi.org/10.1016/j.cam.2005.01.041>
//! - A. Deano and J. Segura, "Transitory minimal solutions of hypergeometric recursions and
//!   pseudoconvergence of associated continued fractions", Math. Comp. 76 (2007), 879-901.
//!   <https://doi.org/10.1090/S0025-5718-07-01934-5>
//! - NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/15.5>,
//!   <https://dlmf.nist.gov/15.8>

use core::f64::consts::LN_2;

use crate::utils::MACHEP;
use crate::xsf::cephes::HYP2F1_ETHRESH;
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, LOG2, frexp_f64, ldexp_f64};

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
    loop {
        let (pk, rk, ck) = ((p + k).hi, (r + k).hi, (c + k).hi);
        // Multiply by the two factors separately, after scaling `|t|` to `[1, 2)`, so that e.g.
        // tiny p and r don't underflow, and fail if a factor or product over- or underflows.
        for (f, num) in [(pk / (k + 1.0), pk), (rk / ck * w, rk)] {
            if !f.is_finite() || f.abs() < f64::MIN_POSITIVE && num != 0.0 {
                return None;
            }
            if t != 0.0 {
                let (m, s) = frexp_f64(t);
                (t, te) = (2.0 * m, te + s - 1);
            }
            t *= f;
            if !t.is_finite() {
                return None;
            }
        }
        let term = ldexp_f64(t, te);
        k += 1.0;
        sum += term;
        sum_abs += term.abs();
        k_sum += k * term;
        k_sum_abs += k * term.abs();
        if !k_sum_abs.is_finite() {
            return None;
        }
        // For j >= k > kmin, the next term ratios are w (p + j) / (j + 1) * (r + j) / (c + j),
        // where both factors are positive and tend monotonically to 1, so they are at most `rho`,
        // which bounds the rest of the series geometrically.
        if k > kmin {
            let rho = w * ((p + k).hi / (k + 1.0)).max(1.0) * ((r + k).hi / (c + k).hi).max(1.0);
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
/// - The sum over `m` of `|t_m| |d_0 + ... + d_m| / |d_m|` is at most the sum over `j` of
///   `|d_j| T_j`, with `T_j` the sum over `m >= j` of `|t_m| / |d_m|`, which we accumulate on the
///   way down. The errors of `y_n` and `y_(n+1)` act like `t_n`.
#[allow(clippy::many_single_char_names)]
fn recur_down(p: Dd, q: Dd, c: f64, z: f64, n: u32, y_n: Scaled, y_n1: Scaled) -> Option<Scaled> {
    // y_(k+1) and y_k, scaled by 2^-e
    let (mut y1, mut y0, mut e) = (ldexp_f64(y_n1.m, y_n1.e - y_n.e), y_n.m, y_n.e);
    // d_k (scaled), T_k (scaled inversely), and the error estimate
    let (mut d, mut t_sum) = (1.0_f64, y_n.err + y_n1.err);
    let mut err = t_sum;
    for k in (1..=n).rev() {
        // (C - 1, C - p and C - q are computed from c, since they can be tiny)
        let (cc, cc1) = (c + f64::from(k), c + f64::from(k - 1));
        let cc_p = (-p + c + f64::from(k)).hi;
        let cc_q = (-q + c + f64::from(k)).hi;
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
        d *= coef_c / coef_a * (y1 / ym);
        t_sum += t_m / d.abs();
        err += d.abs() * t_sum;
        (y1, y0) = (y0, ym);
        if let Some(s) = rescale_exp(y0) {
            (y0, y1) = (ldexp_f64(y0, -s), ldexp_f64(y1, -s));
            e += s;
        }
        if let Some(s) = rescale_exp(d) {
            // (this leaves |d| T, and so `err`, unchanged)
            (d, t_sum) = (ldexp_f64(d, -s), ldexp_f64(t_sum, s));
        }
        if !err.is_finite() {
            return None;
        }
    }
    Some(Scaled { m: y0, e, err })
}

/// `(1 - z)^p F(p, q; c; z)` for `z <= 0` (steps 2 and 3), with `w = z / (z - 1)`
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn hyp2f1_recur_form(p: Dd, q: Dd, c: f64, z: f64, w: Dd) -> Option<Scaled> {
    // The series takes about (|p| w + 40) / (1 - w) terms, so give up early if that's too many:
    // the terms peak near k = |p| w / (1 - w), and then decay like w^k, which takes about
    // ln(2^53) / (1 - w) ~ 37 / (1 - w) more terms.
    if (p.hi.abs() * w.hi + 40.0) / (1.0 - w.hi) > HYP2F1_RECUR_MAXITER {
        return None;
    }
    // the number of series terms plus recurrence steps left
    let mut budget = HYP2F1_RECUR_MAXITER;
    let r = -q + c;
    let direct = series_w(p, r, Dd::new(c), w, &mut budget);
    if c > q.hi.max(0.0) || direct.is_some_and(|y| y.err <= HYP2F1_ETHRESH) {
        return direct;
    }
    // the smallest n with c + n > max(0, q)
    let n = (q.hi.max(0.0) - c).floor() + 1.0;
    budget -= n;
    if budget < 0.0 {
        return direct;
    }
    let (c_n, r_n) = (Dd::new(c) + n, r + n);
    let recur = series_w(p, r_n, c_n, w, &mut budget).and_then(|y_n| {
        let y_n1 = series_w(p, r_n + 1.0, c_n + 1.0, w, &mut budget)?;
        // (the recurrence is symmetric in p and q)
        recur_down(p, q, c, z, n as u32, y_n, y_n1)
    });
    match (direct, recur) {
        (Some(d), Some(r)) => Some(if d.err <= r.err { d } else { r }),
        (d, r) => d.or(r),
    }
}

/// `F(a, b; c; x)` for real `x < 1` with `x != 0` (see the module docs), or `None` if the
/// estimated relative error is too large.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub(crate) fn hyp2f1_recur(a: f64, b: f64, c: f64, x: f64) -> Option<f64> {
    // (beyond this, the coefficients of the recurrence can overflow)
    if !(a.abs().max(b.abs()).max(c.abs()) <= 1e15 && x < 1.0 && x != 0.0) {
        return None;
    }
    let (a, b, c1) = (Dd::new(a), Dd::new(b), Dd::new(c));
    let (ca, cb) = (c1 - a, c1 - b);
    // the combinations `(p, q, s)` of steps 1 and 2, with `F(a, b; c; x)` equal to
    // `(1 - x)^s (1 - z)^p F(p, q; c; z)`
    let (z, w, mut combinations) = if x < 0.0 {
        let combinations = [(a, b, -a), (b, a, -b), (ca, cb, -b), (cb, ca, -a)];
        // (if x is huge, w rounds to 1, and the series fails)
        (x, Dd::new(x) / (Dd::new(x) - 1.0), combinations)
    } else {
        let (zero, cab) = (Dd::new(0.0), ca - b);
        let combinations = [(a, cb, zero), (cb, a, cab), (ca, b, cab), (b, ca, zero)];
        (x / (x - 1.0), Dd::new(x), combinations)
    };
    // (`log1p` of `DoubleDouble` is only as accurate as that of `f64`, but `log` is accurate)
    let ln_1mx = (Dd::new(1.0) - x).log();
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
        if best.is_none_or(|(best, _)| y.err < best.err) {
            best = Some((y, s * ln_1mx));
        }
        if y.err <= HYP2F1_ETHRESH {
            break;
        }
    }
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
    use super::hyp2f1_recur;

    #[test]
    fn test_hyp2f1_recur_overflow() {
        // the terms or their ratios over- or underflow: fail instead of returning e.g. NaN or 1.0
        // (direct sums of the series in high precision)
        for (a, b, c, x, desired) in [
            (10.1, 0.2, 1.4e-309, -0.18, -9.670_216_976_660_568e307),
            (5e-324, 2000.3, 10.7, 0.5, 6.216_729_240_191_533e249),
            (10.0, 9e306, 1e308, 0.5, 1.584_770_156_297_082),
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
