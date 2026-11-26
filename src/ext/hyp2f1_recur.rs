//! The Gauss hypergeometric function `F(a, b; c; x)` for real `x < 1` and large real parameters,
//! by backward recurrence in `c`. Not a port: `hyp2f1` (in `cephes/hyp2f1.rs`) tries it first for
//! large parameters, where its series and transformations can lose all precision, for a tiny
//! negative `a` or `b`, near a negative integer `a`, `b`, `c - a` or `c - b`, and near an integer
//! `b - a` or `c - a - b` (the log case of steps 4 and 5).
//!
//! # Algorithm
//!
//! 1. Write `F(a, b; c; x) = (1 - x)^s (1 - z)^p F(p, q; c; z)` with `z <= 0` (DLMF 15.8.1), in 4
//!    combinations (these and the same with `a` and `b` swapped):
//!    - for `x < 0`: `z = x`, and `(p, q, s)` is `(a, b, -a)` or `(c - a, c - b, -b)`
//!    - for `0 < x < 1`: `z = x / (x - 1)`, and `(p, q, s)` is `(a, c - b, 0)` or
//!      `(c - a, b, c - a - b)`
//!
//! 2. Sum `(1 - z)^p F(p, q; C; z) = F(p, C - q; C; w)` (DLMF 15.8.1) at `C = c`, with
//!    `w = z / (z - 1)` (so `0 <= w < 1`, and `w = x` if `x > 0`).
//!    - For `k > -p`, its terms have one sign while `C + k` and `C - q + k` do, e.g. for all `k` if
//!      `C > max(0, q)`. So it cancels little if `p >= 0` or `|p|` is small.
//!    - The twins `(p, q)` and `(c - q, c - p)` share this series, as `F` is symmetric in its
//!      first two parameters, so the twin reuses its sum (if it succeeded).
//!    - If `p = -n` is a negative integer (with `n <= HYP2F1_RECUR_MAXITER`), its terms can
//!      alternate in sign and cancel. So unless its sum or that of `F(-n, q; c; z)` is within
//!      `HYP2F1_ETHRESH`, we also recur up in the degree, from `F(0) = 1` and `F(-1)` to `F(-n)`
//!      (DLMF 15.5.11), with the error bound of step 3, for `F(-n, c - q; c; w)` and the form
//!      (DLMF 15.8.7) `F(-n, q; c; z) = (c - q)_n / (c)_n F(-n, q; q - c - n + 1; 1 - z)`, often
//!      more stable if `c < 0`. If neither is within `HYP2F1_ETHRESH` either, we also sum
//!      `F(-n, c - q; c; w)` in double-double, as near a root both recurrences can lose about as
//!      many digits as the sums. The best of these is returned if its error bound is at most
//!      `HYP2F1_RECUR_TOL`.
//!
//! 3. If `c <= max(0, q)` and that loses too much precision, also sum it at `C = c + n` and
//!    `C = c + n + 1`, with the smallest `n` with `c + n > max(0, q)`, and recur down to `C = c`
//!    (DLMF 15.5.18).
//!    - For `Re z < 1/2`, `F(p, q; c + n; z)` is minimal as `n -> inf` (Gil, Segura and Temme
//!      2007, sec. 3.2), so this is stable for `c + n` much larger than `|p|` and `|q|`.
//!    - In between, it can be unstable (Deano and Segura 2007), so `recur_down` bounds its error.
//!
//! 4. For `x < -1`, first try the 1/x transformation (DLMF 15.8.2), as `w` is near 1 for large
//!    `|x|`, where the series take about `|p| |x|` terms.
//!    - Its two `F(..; 1/x)` are computed as above, with `w = 1 / (1 - x)`.
//!    - If `b - a` is not an integer and `a` or `c - b` is a non-positive integer, the term of `b`
//!      is 0 (at a pole of its Gamma factors), which leaves that of `a`: for `a = -n`, the
//!      polynomial in `1/x`.
//!    - The 4 combinations follow only if its error bound (incl. the cancellation of its two
//!      terms) exceeds `HYP2F1_ETHRESH`.
//!    - If `b - a` is near an integer (a pole of its Gamma factors), its two terms combined come
//!      first (DLMF 15.8.9 at an integer). Their series in `w` has weights, so it can't be
//!      transformed or recurred like the others, and it can cancel whatever the signs of the
//!      parameters (see `recip_log`).
//!
//! 5. If `-1 <= x < 0.9` and all error bounds so far exceed `HYP2F1_RECUR_TOL`, also try the same
//!    transformation of `F(a, b; c; x)` if `x < 0` (i.e. DLMF 15.8.3), or of the Pfaff form
//!    `(1 - x)^-b F(c - a, b; c; x / (x - 1))` if `0 < x < 1` (i.e. the 1 - x transformation,
//!    DLMF 15.8.4), with series in `w = 1 / (1 - x)` or `w = 1 - x`.
//!    - The 4 combinations can all cancel if `a` and `b` are large with opposite signs and
//!      `0 < x < 1` (or with the same sign and `x < 0`), while these two terms rarely do.
//!    - For `x >= 0.9`, where `w = x` is close to 1, it is tried first instead, as in step 4.
//!    - The log case of step 4 applies too: if `b - a` is near an integer for `x < 0`, or
//!      `c - a - b` for `0 < x < 1` (DLMF 15.8.10).
//!
//! 6. If `c < 0` and all error bounds so far exceed `HYP2F1_RECUR_TOL`, also try `F` from two
//!    other solutions of its differential equation (DLMF 15.10.19): the second solution at 0,
//!    `|x|^(1-c) F(a - c + 1, b - c + 1; 2 - c; x)`, and one term of the transformation of step 4
//!    (or of step 5 if `0 < x < 1`).
//!    - That replaces its other term, the one with a third parameter of at most 1
//!      (`1 - |b - a|` if `x < 0`), by a solution with the third parameter `2 - c > 2`.
//!    - For `c > 0`, it rarely had the smallest error bound in our tests.
//!
//! 7. If all error bounds still exceed `HYP2F1_RECUR_TOL`, recur down in `c` as in step 3, but
//!    from other start values: those of step 3 are sums of the same series, which can lose
//!    precision for every `n` (e.g. if `p` is large and negative), while other forms at `c + n`
//!    don't.
//!    - `F(p, q; c + k; z)` is a multiple of (DLMF 15.8.1)
//!      - `F(a, b; c + k; x)` for `(p, q, z) = (a, b, x)`
//!      - `(1 - x)^k F(a + k, b + k; c + k; x)` for `(c - a, c - b, x)`
//!      - `F(a, b + k; c + k; x)` for `(a, c - b, x / (x - 1))`
//!      - `F(a + k, b; c + k; x)` for `(c - a, b, x / (x - 1))`
//!    - For each, `y_n` and `y_(n+1)` (with `n` as in step 3, but at least 1) come from steps 1
//!      to 5, then recur down to `k = 0`. If `z > 0`, the recurrence need not be stable, but its
//!      error bound still holds.
//!    - Those with `c > max(0, q)` (so `n = 1`) are skipped unless steps 1 to 6 got within
//!      `HYP2F1_RECUR_SHIFT_TOL`: their start values at `c + 1` and `c + 2` hardly differ from
//!      `F(a, b; c; x)`, and in our tests improved the bound at most 15-fold.
//!
//! 8. If they still do in the log case of steps 4 and 5, interpolate that transformation in `b`,
//!    between Chebyshev nodes around the integer `b - a` (or `c - a - b`).
//!    - `F` is entire in `b`, and off the integer the two terms lose only about the inverse of
//!      the distance, while the series of `recip_log` can cancel completely.
//!    - Its error bound is an estimate (see `recip_interp`), unlike those of the other steps.
//!
//! Finally, return the result with the smallest error bound, which is too large if it exceeds
//! `HYP2F1_RECUR_TOL`.
//!
//! For an integer `c <= 0` at or below an integer `a` or `b < 0`, where the steps above pass
//! through the poles of `1 / (c)_k`, only the polynomial is computed (DLMF 15.2.5), for any finite
//! `x`.
//!
//! For example, the sum of step 2 cancels for `F(-522, -123.5; 533.5; -0.9)`, the log case of
//! step 4 for `F(-900.5, -800.5; 50; -40)`, and all 4 combinations for
//! `F(-975.2, 241.4; 0.00135; 0.39)`. Step 7 gets `F(758.1, 7.7; -151.2; -0.3)` from
//! `F(758.1, 159.7; 0.8; -0.3)` and `F(758.1, 160.7; 1.8; -0.3)`.
//!
//! The parameters are double-double, so that e.g. `C - q + k` is accurate if `q = c - b` and `c`
//! is tiny, and so are `s ln(1 - x)` and the `ln Gamma` of the transformations. The values are
//! kept as `m 2^e`, to avoid over- and underflow. The cost is at most `HYP2F1_RECUR_MAXITER`
//! series terms and recurrence steps per combination, which is skipped (or only its step 3) if it
//! would exceed that: about `(|p| w + 40) / (1 - w)` terms (but at least
//! `(w (c - q) - c) / (1 - w)`), plus `max(0, q) - c` steps. Step 6 costs at most as much as steps
//! 1 to 4, and step 7 up to 8 times as much as steps 1 to 5, but only if those fail, and step 8
//! `HYP2F1_RECUR_NODES` times the transformation of step 4.
//!
//! # References
//!
//! - A. Gil, J. Segura and N.M. Temme, "Numerically satisfactory solutions of hypergeometric
//!   recursions", Math. Comp. 76 (2007), 1449-1468. <https://doi.org/10.1090/S0025-5718-07-01918-7>
//! - A. Deano and J. Segura, "Transitory minimal solutions of hypergeometric recursions and
//!   pseudoconvergence of associated continued fractions", Math. Comp. 76 (2007), 879-901.
//!   <https://doi.org/10.1090/S0025-5718-07-01934-5>
//! - R.C. Forrey, "Computing the hypergeometric function", J. Comput. Phys. 137 (1997), 79-100.
//!   <https://doi.org/10.1006/jcph.1997.5794>
//! - N. Michel and M.V. Stoitsov, "Fast computation of the Gauss hypergeometric function with all
//!   its parameters complex with application to the Poschl-Teller-Ginocchio potential wave
//!   functions", Comput. Phys. Commun. 178 (2008), 535-551.
//!   <https://doi.org/10.1016/j.cpc.2007.11.007>
//! - NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/5.5.E3>,
//!   <https://dlmf.nist.gov/5.11.E1>, <https://dlmf.nist.gov/5.15.E9>,
//!   <https://dlmf.nist.gov/15.5>, <https://dlmf.nist.gov/15.8>, <https://dlmf.nist.gov/15.10.E19>

use core::cmp::Ordering;
use core::f64::consts::{LN_2, LOG2_E, PI};

use crate::numpy::sinc;
use crate::utils::{DENORM_MIN, MIN_SUBNORMAL_EXP, neg_one_pow, pow2};
use crate::xsf::cephes::consts::MACHEP;
use crate::xsf::cephes::dd_real::{
    DoubleDouble as Dd, LN_PI_DD, LN_SQRT_2PI_DD, LOG2, frexp_f64, ldexp_f64,
};
use crate::xsf::cephes::hyp2f1::{HYP2F1_ETHRESH, HYP2F1_RECUR_NEAR1};
use crate::xsf::cephes::{lgam_sgn, lgam_stirling_corr, psi};
use crate::xsf::{Real, cospi, sinpi};

/// The largest accepted error bound. These bounds are pessimistic (the actual errors are usually
/// 10 to 300 times smaller), and the methods of the caller are often far worse if this fails.
pub(crate) const HYP2F1_RECUR_TOL: f64 = 1e-10;

/// The largest error bound of steps 1 to 6 for which step 7 also tries the forms with
/// `c > max(0, q)` (empirical)
const HYP2F1_RECUR_SHIFT_TOL: f64 = 100.0 * HYP2F1_RECUR_TOL;

/// The largest distance of `b - a` to an integer for which `hyp2f1_recip` first tries `recip_log`,
/// as its two terms lose about the inverse of that distance (and for which `hyp2f1` first tries
/// `hyp2f1_recur`)
pub(crate) const HYP2F1_RECUR_NEAR_INT: f64 = 1e-3;

/// The number of nodes of `recip_interp` (even, so that none is at the integer), and half the width
/// of their interval in `b` (empirical: closer nodes cancel more, wider ones interpolate worse)
const HYP2F1_RECUR_NODES: usize = 16;
const HYP2F1_RECUR_NODE_WIDTH: f64 = 0.05;

/// The maximum number of series terms plus recurrence steps, per combination
pub(crate) const HYP2F1_RECUR_MAXITER: f64 = 1e6;

/// The combinations `(p, q, s)` of step 1
type Combinations = [(Dd, Dd, Dd); 4];

/// `(A_k, B_k, C_k, b_abs)` of `recur_down`
type Coefs = (f64, f64, f64, f64);

/// `m 2^e exp(l)`, with a bound of its relative error
#[derive(Clone, Copy, Debug)]
pub(crate) struct Scaled {
    pub(crate) m: f64,
    pub(crate) e: i32,
    pub(crate) l: Dd,
    pub(crate) err: f64,
}

impl Scaled {
    /// `Some(m 2^e)` if `m` is nonzero and `err` is finite
    pub(crate) fn new(m: f64, e: i32, err: f64) -> Option<Self> {
        let l = Dd::ZERO;
        (m != 0.0 && err.is_finite()).then_some(Self { m, e, l, err })
    }

    pub(crate) fn with_l(self, l: Dd) -> Self {
        Self { l, ..self }
    }

    /// This times `exp(l)`
    pub(crate) fn mul_exp(self, l: Dd) -> Self {
        self.with_l(self.l + l)
    }

    /// A zero value whose scaled absolute error `|self| (1 + self.err)` is stored in `err`
    ///
    /// Only valid as an input to `add_exp_in` / `add_exp`; `fold_exp` and `to_f64` interpret
    /// `err` as a relative error and must not be called on the returned value.
    pub(crate) fn into_zero_err(self) -> Self {
        let (m, k) = frexp_f64(self.m.abs());
        Self {
            m: 0.0,
            e: self.e + k,
            l: self.l,
            err: m * (1.0 + self.err),
        }
    }

    /// The value as an `f64` (which can overflow or underflow), and its error bound (not counting
    /// the rounding of an underflow)
    pub(crate) fn to_f64(self) -> (f64, f64) {
        let y = fold_exp(self);
        (ldexp_f64(y.m, y.e), y.err)
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

const RESCALE_LIMIT_EXP: i32 = 500;

/// The exponent `k` that scales `|x|` to about 1 by `2^-k`, if `|x|` is outside `[2^-500, 2^500]`
fn rescale_exp(x: f64) -> Option<i32> {
    const MIN: u64 = pow2(-RESCALE_LIMIT_EXP).to_bits();
    const MAX: u64 = pow2(RESCALE_LIMIT_EXP).to_bits();
    // (|x| is finite and nonzero, so log2(|x|) is in [-1075, 1024])
    let x = x.abs();
    // For positive finite values the bit order matches the numerical order.
    let bits = x.to_bits();
    (!(MIN..=MAX).contains(&bits) && x != 0.0 && x.is_finite()).then(|| x.log2().round() as i32)
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
/// With `weights = Some((beta_0, err, eps))` (step 4), it sums `t_k beta_k` instead, with `err`
/// bounding the error of `beta_0`, and
///
/// ```text
/// beta_(k+1) = beta_k + (1 + eps beta_k) d_k
/// 1 + eps d_k = (1 - eps/(p + k)) (1 - eps/(r + k)) / ((1 - eps/(k + 1)) (1 - eps/(c + k)))
/// ```
///
/// i.e. `d_k = 1/(k + 1) + 1/(c + k) - 1/(p + k) - 1/(r + k)` if `eps = 0`.
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
    weights: Option<(f64, f64, f64)>,
    kmax: f64,
) -> Option<Scaled> {
    let (w, w_lo) = (w.hi, w.lo);
    // from here on, the factors of the term ratios have constant signs
    let kmin = [-p.hi, -r.hi, -c.hi].into_iter().fold(0.0, f64::max).ceil();
    // `beta_k` and a bound of its absolute error, the sum of `|t_k|` times that bound, and the
    // sum of `|parts|` of `d_(k-1)`
    let (mut beta, mut beta_err, eps) = weights.unwrap_or((1.0, 0.0, 0.0));
    let (mut sum_err, mut d_abs) = (beta_err, 0.0);
    // The terms are `t 2^te`, and the sums are scaled by `2^-e`, since the terms can underflow
    // before `k = -c`, and be huge after it.
    let (mut sum, mut sum_abs, mut k_sum, mut k_sum_abs) = (beta, beta.abs(), 0.0, 0.0_f64);
    let (mut t, mut te, mut e, mut k) = (1.0, 0, 0, 0.0);
    // the factors `(p + k) / (k + 1)` and `(r + k) / (c + k)` of the next term ratio, and `p + k`,
    // `r + k` and `c + k`
    let (mut pk, mut rk, mut ck) = (p.hi, r.hi, c.hi);
    let (mut f0, mut f1) = (pk, rk / ck);
    // As in rescale_exp, absolute binary64 encodings preserve the in-range test.
    let tame = pow2(-250).to_bits()..=pow2(250).to_bits();
    while k < kmax {
        let fw = f1 * w;
        if [f0, fw]
            .into_iter()
            .all(|f| tame.contains(&f.abs().to_bits()))
        {
            // (usually: then `|t f0 fw|` is 0 or in `[2^-1000, 2^1000]`, and one rescale equals two
            // `mul_scaled`, as scaling by powers of 2 is exact)
            (t, te) = rescaled(t * f0 * fw, te);
        } else {
            (t, te) = mul_scaled((t, te), f0, pk)?;
            (t, te) = mul_scaled((t, te), fw, rk)?;
        }
        let term = if te == 0 {
            t
        } else if te < -1075 - RESCALE_LIMIT_EXP {
            // |t| <= 2^RESCALE_LIMIT_EXP, so this term rounds to signed zero.
            0.0_f64.copysign(t)
        } else {
            ldexp_f64(t, te)
        };
        let term_abs = term.abs();
        k += 1.0;
        if weights.is_some() {
            let parts = [1.0 / k, 1.0 / ck, -1.0 / pk, -1.0 / rk];
            // the second-order terms of the denominator and the numerator of `1 + eps d`, and
            // the inverse of that denominator
            let (den, num) = (eps * parts[0] * parts[1], eps * parts[2] * parts[3]);
            let scale = 1.0 / (1.0 - eps * (parts[0] + parts[1] - den));
            let d = (parts.into_iter().sum::<f64>() + (num - den)) * scale;
            let parts_abs = parts.into_iter().map(f64::abs).sum::<f64>();
            d_abs = (parts_abs + num.abs() + den.abs()) * scale.abs();
            let ratio = 1.0 + eps * beta;
            beta += ratio * d;
            // (the error of `beta_(k-1)`, and 2 roundings per part, and 4 additions)
            beta_err = beta_err * (1.0 + (eps * d).abs())
                + MACHEP * (5.0 * d_abs * ratio.abs() + beta.abs());
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
        // With weights, the parts of `d_j` decrease as well, and `|1 + eps beta_j|` grows by at
        // most `u = 1 + |eps| d_abs` per term, which is in `rho`, so
        // `|beta_j| <= |beta_k| + (j - k) u^(j-k) |1 + eps beta_k| d_abs`.
        if k > kmin {
            let rho = w * f0.max(1.0) * f1.max(1.0) * (1.0 + eps.abs() * d_abs);
            let beta_max = beta.abs() + (1.0 + eps * beta).abs() * d_abs / (1.0 - rho);
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
        .rfold((1.0_f64, 0.0), |(u, err), (t_m, rho)| (1.0 + u * rho, err + t_m * u.abs()));
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

/// The polynomial forms of step 2 of `combs` with `p = -n`; `ln_1mx` is `ln(1 - x)`, and
/// `(one_mz, ln_1mz)` is `1 - z` and its log
fn hyp2f1_poly(
    combs: Combinations,
    c: Dd,
    w: Dd,
    ln_1mx: Dd,
    (one_mz, ln_1mz): (Dd, Dd),
) -> Option<Scaled> {
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
            let direct = min_err(sum(q, Dd::ONE - one_mz, l_z), sum(r, w, l));
            if err_of(direct) <= HYP2F1_ETHRESH {
                return direct;
            }
            // `(1 - z)^p F(p, q; c; z) = F(p, c - q; c; w)`, as in step 2
            let y_w = poly_up(n, r, c, w).map(|y| y.with_l(l));
            // `(c - q)_n / (c)_n = Gamma(c - q + n) Gamma(c) / (Gamma(c - q) Gamma(c + n))`
            let y_1mz = mul_gamma_ratio([r + n, c], [r, c + n], || {
                let y = poly_up(n, q, q - c + (1.0 - n), one_mz)?;
                Some(y.with_l(l_z))
            });
            let best = min_err(direct, min_err(y_w, y_1mz));
            if err_of(best) <= HYP2F1_ETHRESH {
                return best;
            }
            // (e.g. near a root, where the sum of `F(-8, 18.02; 10.81; 0.909)` cancels by 1e12 and
            // the better recurrence has an error bound of 2e-11), only if accurate, as the
            // rounding of `w` can make it less so than the paths of the caller
            let y_dd = poly_dd(n, r, c, w, false).filter(|y| y.err <= HYP2F1_ETHRESH);
            min_err(best, y_dd.map(|y| y.with_l(l)))
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
    // At exact Euler balance, both twin pairs can have the same parameters.
    // Keep one pair so each still uses the same first series and budget.
    let bits = |(p, q, s): (Dd, Dd, Dd)| [p, q, s].map(|d| (d.hi.to_bits(), d.lo.to_bits()));
    let duplicate = bits(combs[0]) == bits(combs[1]) && bits(combs[2]) == bits(combs[3]);
    // the step-2 sums of the twins
    let mut shared = [None; 2];
    let order = order
        .into_iter()
        .filter(|&i| !duplicate || i == 0 || i == 3);
    best_of(order.map(|i| {
        let (p, q, s) = combs[i];
        // `c - q` is the `p` of the twin, while `c - (c - p)` would round a tiny `p`
        let r = combs[3 - i].0;
        let y = hyp2f1_recur_form(p, q, r, c, z, w, &mut shared[i.min(3 - i)])?;
        Some(y.with_l(s * ln_1mx))
    }))
}

/// `ln(v)` for `v > 0`, after scaling `v` to `[1, 2)` by a power of 2, since `DoubleDouble::log`
/// fails for `v` beyond about `e^709` (as it computes `exp(-ln(v))`)
pub(crate) fn ln_dd(v: impl Into<Dd>) -> Dd {
    let v = v.into();
    let (_, k) = frexp_f64(v.hi);
    v.ldexp(1 - k).log() + LOG2 * f64::from(k - 1)
}

/// A lower bound `d` of the distance of `g.hi + [-|g.lo|, |g.lo|]` to the integers (negative if it
/// contains one), or `None` if `g < 1/2` and `d <= 0`, i.e. at (or too close to) a pole of `Gamma`
pub(crate) fn pole_distance(g: Dd) -> Option<f64> {
    let d = (g.hi - g.hi.round()).abs() - g.lo.abs();
    (!(g.hi < 0.5 && d <= 0.0)).then_some(d)
}

/// The relative rounding error of a few `f64` operations, in the error bounds of `ln_gamma`
const F64_ROUNDING: f64 = 4.0 * MACHEP;
/// The relative rounding error of a few double-double operations, in the same
const DD_ROUNDING: f64 = 4.0 * MACHEP * MACHEP;

/// `ln|Gamma(g)|`, the sign of `Gamma(g)`, and a bound of the absolute error of the former, or
/// `None` at a pole or for a non-finite `g`:
///
/// - for `g >= 13`, by Stirling's formula in double-double, with the correction of cephes' `lgam`
///   (DLMF 5.11.1)
/// - for `g <= -13` or within `|g.lo|` of a pole, by the reflection formula
///   `Gamma(g) Gamma(1 - g) = pi / sin(pi g)` (DLMF 5.5.3), with `ln_sinpi`
/// - else by `lgam` at `g.hi`, to first order in `g.lo`
///
/// The second-order term `f'' g.lo^2 / 2` of `f = ln Gamma` is bounded between `g.hi` and
/// `g.hi + g.lo`, with `d` the smallest distance to a pole there: `f'' < 1/d^2 + pi^2` if
/// `g < 1/2`, and `f'' < 1/g + 1/g^2` if `g > 0`.
pub(crate) fn ln_gamma(g: impl Into<Dd>) -> Option<(Dd, f64, f64)> {
    let g = g.into();
    // (a non-finite `g` would recurse forever in the reflection)
    if !g.hi.is_finite() || !g.lo.is_finite() {
        return None;
    }
    if g.hi >= 13.0 {
        let v = (g - 0.5) * ln_dd(g) - g + LN_SQRT_2PI_DD + lgam_stirling_corr(g.hi);
        // (the correction is below 0.01, and its error below eps)
        return Some((v, 1.0, 2.0 * MACHEP + 16.0 * MACHEP * MACHEP * v.hi.abs()));
    }
    match pole_distance(g) {
        Some(d) if g.hi > -13.0 => {
            let (v, sign) = lgam_sgn(g.hi);
            if !v.is_finite() {
                return None;
            }
            let dpsi = psi(g.hi) * g.lo;
            let err =
                F64_ROUNDING * (v.abs() + dpsi.abs() + 1.0) + trigamma_lo(g, d) * g.lo.abs() / 2.0;
            Some((Dd::new(v) + dpsi, f64::from(sign), err))
        },
        d => {
            let (ln_sin, sign, sin_err) = ln_sinpi(g, d)?;
            let (v, _, err) = ln_gamma(-g + 1.0)?;
            let err = err + sin_err + DD_ROUNDING * (v.hi.abs() + ln_sin.hi.abs() + 1.0);
            Some((LN_PI_DD - v - ln_sin, sign, err))
        },
    }
}

/// `ln|sin(pi g)|`, its sign, and a bound of its absolute error, for `g` at a distance `d` from
/// the integers (`None` if within `|g.lo|` of one):
///
/// - to first order in `g.lo`, with the second-order term bounded by `|f''| < (pi / 2d)^2` for
///   `f = ln|sin(pi g)|` (since `|sin(pi g)| >= 2d`)
/// - else for `g = k + e`, by `ln|sin(pi e)| = ln(pi |e|) - x^2/6 - x^4/180 - ...` with `x = pi e`,
///   whose tail is below `x^4/100` for `|x| <= pi/1000`; `None` if `e = 0`, or if `|e| > 1/1000`
///   (only for `|g| > 2^42`)
fn ln_sinpi(g: Dd, d: Option<f64>) -> Option<(Dd, f64, f64)> {
    if let Some(d) = d {
        let s = sinpi(g.hi);
        // d/dg ln|sin(pi g)| = pi cot(pi g)
        let dsin = PI * cospi(g.hi) / s * g.lo;
        let ln_sin = s.abs().ln();
        let err = F64_ROUNDING * (1.0 + ln_sin.abs() + dsin.abs())
            + (PI / (2.0 * d) * g.lo).powi(2) / 2.0;
        return Some((Dd::new(ln_sin + dsin), s.signum(), err));
    }
    let k = g.hi.round();
    let e = g - k;
    if e.hi == 0.0 || e.hi.abs() > 1e-3 {
        return None;
    }
    let ln_e = ln_dd(if e.hi < 0.0 { -e } else { e });
    let x2 = (PI * e.hi).powi(2);
    let err = DD_ROUNDING * (ln_e.hi.abs() + 2.0) + F64_ROUNDING * x2 + x2 * x2 / 100.0;
    Some((LN_PI_DD + ln_e - x2 / 6.0, cospi(k) * e.hi.signum(), err))
}

/// `|g.lo|` times the bound of `(ln Gamma)''` in `ln_gamma` (0 if `g.lo = 0`, even if it overflows)
pub(crate) fn trigamma_lo(g: Dd, d: f64) -> f64 {
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

/// `(Gamma(v + e) / Gamma(v) - 1) / e` for `|e| <= 0.04`, which is `psi(v)` at
/// `e = 0`, and a bound of its absolute error, or `None` at a pole of either:
///
/// - for `v < 1/4`, by the reflection formula (DLMF 5.5.3), with `u = 1 - v - e` and
///   `s = sin(pi v) / sin(pi (v + e))`:
///
///   ```text
///   Gamma(v + e) / Gamma(v) = s Gamma(u + e) / Gamma(u)
///   (s - 1) / e = -pi cos(pi (v + e/2)) sinc(e/2) / sin(pi (v + e))
///   ```
///
/// - else as `expm1(e l) / e`, with `e l = ln Gamma(v + e) - ln Gamma(v)` from `v + n >= 13`,
///   where it is `e psi(x) + e^3 psi''(x) / 24 + e^5 psi''''(x) / 1920` at the midpoint `x`
///   (DLMF 5.15.9 for the derivatives, up to about `e^2 / (80 x^10)` in `l`), minus the
///   `ln(1 + e / (v + i))`, `i < n`.
pub(crate) fn gamma_quot(v: Dd, e: Dd) -> Option<(f64, f64)> {
    if v.hi < 0.25 && e.hi != 0.0 {
        let (ve, mid) = (v + e, v + e * 0.5);
        if v.is_integer() || ve.is_integer() {
            return None;
        }
        // sin(pi g) or cos(pi g), from the distance of g to the nearest integer k
        let reduced = |g: Dd, f: fn(f64) -> f64| {
            let k = g.hi.round();
            neg_one_pow(k) * f((g - k).hi)
        };
        let sin_ve = reduced(ve, sinpi);
        let q_sin = -PI * reduced(mid, cospi) * sinc(0.5 * e.hi) / sin_ve;
        // (a few roundings each of the sine and the cosine, which can be near 0)
        let sin_err = 8.0 * MACHEP * (q_sin.abs() + PI / sin_ve.abs());
        let (q, q_err) = gamma_quot(-ve + 1.0, e)?;
        let y = q_sin + q * (1.0 + e.hi * q_sin);
        let err = sin_err * (1.0 + (e.hi * q).abs()) + q_err * (1.0 + e.hi * q_sin).abs();
        return Some((y, err + 2.0 * MACHEP * y.abs()));
    }
    let e = e.hi;
    // `psi(v)` at e = 0, and for a subnormal `e`, where `e / (v + i)` and `e l` underflow, but
    // this is within about `e (psi'(v) + psi(v)^2) / 2` of it
    if !e.is_normal() {
        return psi_dd(v);
    }
    let n = (13.0 - v.hi).ceil().max(0.0);
    let (mut l, l_abs) = (0..n as u32)
        .map(|i| (e / (v.hi + f64::from(i))).ln_1p() / e)
        .fold((0.0, 0.0), |(l, l_abs), t| (l - t, l_abs + t.abs()));
    let x = v.hi + n + 0.5 * e;
    let x2 = 1.0 / (x * x);
    let psi2 = -x2 * (1.0 + 1.0 / x + x2 * (0.5 - x2 / 6.0 * (1.0 - x2)));
    let psi4 = -6.0 * x2 * x2 * (1.0 + 2.0 / x + 5.0 / 3.0 * x2);
    let psi_x = psi(x);
    l += psi_x + e * e / 24.0 * (psi2 + e * e / 80.0 * psi4);
    let y = (e * l).exp_m1() / e;
    // (the roundings of the terms, and of `psi` as in `psi_dd`)
    let l_err = MACHEP * (4.0 * l_abs + 8.0 * (psi_x.abs() + 1.0));
    Some((y, l_err * (1.0 + (e * y).abs()) + 2.0 * MACHEP * y.abs()))
}

/// `y()` times `Gamma(num_0) Gamma(num_1) ... / (Gamma(den_0) Gamma(den_1) ...)`, with the error
/// bounds of their `ln_gamma` added; `y` is only evaluated if no argument is at a pole
pub(crate) fn mul_gamma_ratio(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    y: impl FnOnce() -> Option<Scaled>,
) -> Option<Scaled> {
    let num = num.into_iter().map(|g| (g, 1.0));
    let mut gammas = num.chain(den.into_iter().map(|g| (g, -1.0)));
    let (ln_g, sign, err) =
        gammas.try_fold((Dd::ZERO, 1.0, 0.0), |(ln_g, sign, err), (g, power)| {
            let (v, v_sign, v_err) = ln_gamma(g)?;
            Some((ln_g + v * power, sign * v_sign, err + v_err))
        })?;
    let mut y = y()?;
    (y.m, y.l, y.err) = (sign * y.m, y.l + ln_g, y.err + err);
    Some(y)
}

/// The sum of the `terms`, with an error bound that includes their cancellation
pub(crate) fn add_exp(terms: impl IntoIterator<Item = Scaled, IntoIter: Clone>) -> Option<Scaled> {
    add_exp_in(terms, |r| r.hi.exp() * (1.0 + r.lo), MACHEP)
}

/// `add_exp` in the arithmetic of `T`, whose operations and `exp` have the relative error `eps`
pub(crate) fn add_exp_in<T: Real>(
    terms: impl IntoIterator<Item = Scaled, IntoIter: Clone>,
    exp: impl Fn(Dd) -> T,
    eps: f64,
) -> Option<Scaled> {
    // normalize each value and error independently into `l` so neither a subnormal value nor a
    // tiny value with a large relative error can underflow and discard its uncertainty
    let terms = terms.into_iter().map(|y| {
        let (m, km) = frexp_f64(y.m);
        let l_val = y.l + LOG2 * (i64::from(y.e) + i64::from(km)) as f64;
        // (3 eps for `exp(r)` and the products; a zero term stores its scaled absolute error
        // in `y.err`)
        let (err_m, ke) = if y.m == 0.0 {
            let (em, ke) = frexp_f64(y.err * (1.0 + 3.0 * eps));
            (em, i64::from(ke))
        } else {
            let (em, ke) = frexp_f64(y.err + 3.0 * eps);
            (m.abs() * em, i64::from(km) + i64::from(ke))
        };
        let l_err = y.l + LOG2 * (i64::from(y.e) + ke) as f64;
        let l_ref = if y.m == 0.0 && y.err == 0.0 {
            Dd::new(f64::NEG_INFINITY)
        } else if y.m == 0.0 {
            l_err
        } else {
            l_val
        };
        (m, l_val, err_m, l_err, l_ref)
    });
    let lmax = terms
        .clone()
        .map(|(_, _, _, _, l_ref)| l_ref)
        .reduce(|u, v| if u.hi >= v.hi { u } else { v })?;
    let (m, err_abs) =
        terms.fold((T::from_f64(0.0), 0.0), |(m, err_abs), (ym, l_val, err_m, l_err, _)| {
            let t = if ym == 0.0 {
                T::from_f64(0.0)
            } else {
                T::from_f64(ym) * exp(l_val - lmax)
            };
            let err = if err_m == 0.0 {
                0.0
            } else if !err_m.is_finite() {
                err_m
            } else {
                let d = l_err - lmax;
                // range-reduce the full double-double exponent before `exp` so its low word is
                // preserved (e.g. `LOG2 * 1024` for `y.err ~ 1e308`)
                if d.hi > 750.0 {
                    f64::INFINITY
                } else if d.hi < -750.0 {
                    0.0
                } else {
                    let q = (d.hi * LOG2_E).round() as i32;
                    let r = d - LOG2 * f64::from(q);
                    let scale = r.hi.exp() * (1.0 + r.lo);
                    ldexp_f64(err_m * (1.0 + 8.0 * MACHEP) * scale, q)
                }
            };
            (m + t, err_abs + err)
        });
    // (and 1 eps for the sum, as an `f64`)
    let m = m.to_f64();
    Some(Scaled::new(m, 0, err_abs / m.abs() + MACHEP)?.with_l(lmax))
}

/// `hyp2f1_recip` if `m = b - a = n + eps >= 0` is within `HYP2F1_RECUR_NEAR_INT` of an integer
/// `n`, where its two terms cancel, with the terms `k >= n` of the series of the first added to
/// the second analytically (as Forrey 1997 and Michel and Stoitsov 2008 do). At `eps = 0`, that is
/// the limit of the transformation (DLMF 15.8.9). With `w = 1 / (1 - t)`:
///
/// ```text
/// F(a, b; c; t) / Gamma(c) = (1 - t)^-a Gamma(m) / (Gamma(b) Gamma(c - a)) S_1
///     + (-1)^n (1 - t)^-b / (Gamma(a) Gamma(c - b) Gamma(m + 1)) S_2 / sinc(eps)
/// S_1 = sum_(k < n) (a)_k (c - b)_k / ((1 - m)_k k!) w^k
/// S_2 = sum_k (b)_k (c - a)_k / ((m + 1)_k k!) w^k beta_k
/// beta_k = (g(k - eps) / g(k) - 1) / eps
///        = ln(1 - t) + psi(k + 1) + psi(m + k + 1) - psi(b + k) - psi(c - a + k) if eps = 0
/// ```
///
/// `S_2` is the series of step 2 of the second term of `hyp2f1_recip`, but with the weights
/// `beta_k` of `series_weighted` (and `beta_0` by `gamma_quot`). Its terms are
///
/// ```text
/// t_k beta_k = (g(k - eps) - g(k)) / (eps g(0))
/// g(s) = Gamma(b + s) Gamma(c - a + s) / (Gamma(1 + s) Gamma(m + 1 + s)) w^s
/// ```
///
/// so they change sign where `|g|` peaks, and `S_2` cancels if that peak is far above `|g(0)|`
/// (e.g. if `b (c - a) w / (m + 1)` is large), whatever the signs of `b` and `c - a`.
fn recip_log(a: Dd, b: Dd, c: Dd, w: Dd, ln_1mt: Dd) -> Option<Scaled> {
    let (one, m, ca, cb) = (Dd::ONE, b - a, c - a, c - b);
    let n = m.hi.round();
    let eps = m - n;
    let (e, m1) = (eps.hi, m + 1.0);
    let near_int = e.abs() <= HYP2F1_RECUR_NEAR_INT;
    if !near_int || series_too_long(b, ca, m1, w) || n > HYP2F1_RECUR_MAXITER {
        return None;
    }
    // `(w^-eps - 1) / eps`, and a bound of its absolute error (`ln(1 - t)` if `eps` is 0 or
    // subnormal, where `eps ln(1 - t)` underflows)
    let ln_quot = if e.is_normal() {
        let q = (e * ln_1mt.hi).exp_m1() / e;
        (Dd::new(q), F64_ROUNDING * q.abs())
    } else {
        (ln_1mt, 0.0)
    };
    // `beta_0` and a bound of its absolute error, from the factors `1 + eps q` of
    // `g(-eps) / g(0)`, with `q = sign * gamma_quot(v, sign * eps)`
    let (beta, beta_err) = [(one - eps, 1.0), (one + n, 1.0), (b, -1.0), (ca, -1.0)]
        .into_iter()
        .try_fold(ln_quot, |(beta, beta_err), (v, sign)| {
            let (y, err) = gamma_quot(v, eps * sign)?;
            let ratio = 1.0 + e * beta.hi;
            let beta_err = beta_err * (1.0 + (e * y).abs()) + err * ratio.abs();
            Some((beta + sign * y * ratio, beta_err))
        })?;
    let weights = Some((beta.hi, beta_err + MACHEP * beta.hi.abs(), e));
    let mut budget = HYP2F1_RECUR_MAXITER;
    let term2 = mul_gamma_ratio([c], [a, cb, m1], || {
        let mut s2 = series_weighted(b, ca, m1, w, &mut budget, weights, f64::INFINITY)?;
        // (-1)^n / sinc(eps)
        s2.m *= neg_one_pow(n) / sinc(e);
        Some(s2.with_l(-b * ln_1mt))
    })?;
    if n == 0.0 {
        // (S_1 is an empty sum)
        return add_exp([term2]);
    }
    let s1 = series_weighted(a, cb, one - m, w, &mut budget, None, n - 1.0)?;
    let term1 = mul_gamma_ratio([c, m], [b, ca], || Some(s1.with_l(-a * ln_1mt)))?;
    add_exp([term1, term2])
}

/// Step 6: `F(p, q; c; t)` for `p <= q` from `w5 = (-t)^-q F(q, q - c + 1; q - p + 1; 1/t)`, a term
/// of `hyp2f1_recip` (with its `t < 0`), and the second solution at 0, by steps 1 to 3 (DLMF
/// 15.10.19, solved for `F`):
///
/// ```text
/// F(p, q; c; t) = Gamma(q - c + 1) Gamma(1 - p) / (Gamma(1 - c) Gamma(q - p + 1)) w5
///     - Gamma(c - 1) Gamma(q - c + 1) Gamma(1 - p) / (Gamma(1 - c) Gamma(q) Gamma(c - p)) w2
/// w2 = (-t)^(1-c) F(q - c + 1, p - c + 1; 2 - c; t)
/// ```
///
/// `ln_mt` and `ln_1mt` are `ln(-t)` and `ln(1 - t)`.
fn recip_connect(
    (p, q, c): (Dd, Dd, Dd),
    x: f64,
    (ln_mt, ln_1mt): (Dd, Dd),
    w5: impl FnOnce() -> Option<Scaled>,
) -> Option<Scaled> {
    // t and t / (t - 1), from x
    let (t, w) = if x < 0.0 {
        (x, Dd::new(x) / (Dd::new(x) - 1.0))
    } else {
        (x / (x - 1.0), Dd::new(x))
    };
    let (q2, p1, c1) = (q - c + 1.0, -p + 1.0, -c + 1.0);
    let term5 = mul_gamma_ratio([q2, p1], [c1, q - p + 1.0], w5)?;
    let term2 = mul_gamma_ratio([c - 1.0, q2, p1], [c1, q, c - p], || {
        let combs = neg_z_combinations(q2, p - c + 1.0, c1 + 1.0);
        let y = best_combination(combs, c1 + 1.0, t, w, ln_1mt)?;
        Some(Scaled { m: -y.m, ..y }.mul_exp(c1 * ln_mt))
    })?;
    add_exp([term5, term2])
}

/// `F(a, b; c; t)` for `t < 0` by the 1/t transformation (steps 4 and 5, DLMF 15.8.2), where
/// `t = x` if `x < 0`, and `t = x / (x - 1)` if `0 < x < 1`:
///
/// ```text
/// F(a, b; c; t) = G(a, b) (-t)^-a F(a, a - c + 1; a - b + 1; 1/t) + G(b, a) (-t)^-b F(b, ...)
/// G(a, b) = Gamma(c) Gamma(b - a) / (Gamma(b) Gamma(c - a))
/// ```
///
/// If `b - a` is near an integer (a pole of `G`; in step 5 with `0 < x < 1`, that's `c - a - b`),
/// it first tries `recip_log`. `ln_1mx` is `ln(1 - x)`. With `connect`, it is by `recip_connect`
/// instead.
fn hyp2f1_recip(a: Dd, b: Dd, c: Dd, x: f64, ln_1mx: Dd, connect: bool) -> Option<Scaled> {
    let one = Dd::ONE;
    let ln_x = ln_dd(x.abs());
    // z = 1/t, w = 1 / (1 - t), 1 - z, ln(-t) and ln(1 - z), from x (e.g. w = 1 - x if
    // 0 < x < 1)
    let (z, w, one_mz, ln_mt, ln_1mz) = if x < 0.0 {
        let one_mz = (Dd::new(x) - 1.0) / x;
        (1.0 / x, one / (one - x), one_mz, ln_x, ln_1mx - ln_x)
    } else {
        ((x - 1.0) / x, one - x, one / x, ln_x - ln_1mx, -ln_x)
    };
    // `(-t)^-a F(a, a - c + 1; a - b + 1; 1/t)`, also by its polynomial forms, e.g. for
    // `F(c - a, 1 - a; ..)` in step 5 with an integer `a > 0`. Unlike in `hyp2f1_recur_dd`, this
    // keeps the smaller bound even if both are within tolerance, since the two terms can cancel.
    let series = |a: Dd, b: Dd| {
        let (p, q, cc) = (a, a - c + 1.0, a - b + 1.0);
        let combs = neg_z_combinations(p, q, cc);
        let poly = hyp2f1_poly(combs, cc, w, ln_1mz, (one_mz, ln_1mz));
        let y = min_err(poly, best_combination(combs, cc, z, w, ln_1mz))?;
        Some(y.with_l(y.l - a * ln_mt))
    };
    // the term `G(a, b) (-t)^-a F(a, a - c + 1; a - b + 1; 1/t)`
    let term = |a: Dd, b: Dd| mul_gamma_ratio([c, b - a], [b, c - a], || series(a, b));
    let m = b - a;
    let (p, q) = if m.hi < 0.0 { (b, a) } else { (a, b) };
    // ln(-t) + ln(1 - 1/t) = ln(1 - t)
    let ln_1mt = ln_mt + ln_1mz;
    if connect {
        return recip_connect((p, q, c), x, (ln_mt, ln_1mt), || series(q, p));
    }
    let log = recip_log(p, q, c, w, ln_1mt);
    // (at an integer, the two terms have poles)
    if m.is_integer() || err_of(log) <= HYP2F1_ETHRESH {
        return log;
    }
    // (a term is 0 at a pole of its `Gamma(b) Gamma(c - a)`)
    let is_zero = |&(a, b): &(Dd, Dd)| [b, c - a].iter().any(|g| g.hi <= 0.0 && g.is_integer());
    let terms = [(a, b), (b, a)].into_iter().filter(|t| !is_zero(t));
    let sum = terms.map(|(a, b)| term(a, b)).collect::<Option<Vec<_>>>();
    min_err(log, sum.and_then(add_exp))
}

/// `F(a, b; c; x)` by `hyp2f1_recip`: the 1/x transformation (step 4) if `x < 0`, or else the 1/t
/// transformation of its Pfaff form (step 5), with `ln_1mx = ln(1 - x)`; with `connect`, step 6
fn recip_of_x(a: Dd, b: Dd, c: Dd, x: f64, ln_1mx: Dd, connect: bool) -> Option<Scaled> {
    if x < 0.0 {
        hyp2f1_recip(a, b, c, x, ln_1mx, connect)
    } else {
        // F(a, b; c; x) = (1 - x)^-b F(c - a, b; c; x / (x - 1))
        hyp2f1_recip(c - a, b, c, x, ln_1mx, connect).map(|y| y.with_l(y.l - b * ln_1mx))
    }
}

/// Step 8: `recip_of_x` by polynomial interpolation in `b`, between Chebyshev nodes around
/// `b - shift`, where `b - a` (or `c - a - b` if `x > 0`) is an integer.
///
/// Its error bound is heuristic: the bounds of the nodes, plus the distance to the interpolants
/// through every other node, which estimates the error of the interpolation but doesn't bound it.
fn recip_interp(a: Dd, b: Dd, c: Dd, x: f64, ln_1mx: Dd, shift: f64) -> Option<Scaled> {
    const N: usize = HYP2F1_RECUR_NODES;
    // the nodes relative to `b` (Chebyshev points of the first kind)
    let offsets: [f64; N] = core::array::from_fn(|j| {
        HYP2F1_RECUR_NODE_WIDTH * cospi((2 * j + 1) as f64 / (2 * N) as f64) - shift
    });
    let mut ys = [recip_of_x(a, b + offsets[0], c, x, ln_1mx, false)?; N];
    for (y, o) in ys.iter_mut().zip(offsets).skip(1) {
        *y = recip_of_x(a, b + o, c, x, ln_1mx, false)?;
    }
    // the interpolant at `b` through every `step`-th node from `first`
    let interp = |first: usize, step: usize| {
        let nodes = || ys.iter().zip(offsets).skip(first).step_by(step);
        add_exp(nodes().map(|(&y, o)| {
            // the Lagrange basis polynomial of this node at 0, with 3 roundings per factor
            let basis = nodes().filter(|&(_, k)| k != o).map(|(_, k)| k / (k - o));
            let m = y.m * basis.product::<f64>();
            let err = y.err + (3 * N) as f64 * MACHEP;
            Scaled { m, err, ..y }
        }))
    };
    let full = interp(0, 1)?;
    // the distance of `y` to `full`, relative to that (`add_exp` returns them as `m exp(l)`)
    let dist = |y: Scaled| (y.m / full.m * (y.l - full.l).hi.exp() - 1.0).abs();
    let err = full.err + dist(interp(0, 2)?).max(dist(interp(1, 2)?));
    Some(Scaled { err, ..full })
}

/// `F(a, b; c; x)` for finite `x != 0` and an integer `a` or `b = -n < 0`, with `c` not an integer
/// in `(-n, 0]`: the polynomial (DLMF 15.2.4, or 15.2.5 for an integer `c <= -n`) in
/// double-double, as its terms can cancel by `1 / eps` or more, e.g. by 7e15 for
/// `F(-186, -26.75; -333; -9)`; if that is inaccurate, its Pfaff transformation (DLMF 15.8.1)
/// for `x < 1`, and its reversal about 1 (DLMF 15.8.7); at `x = 1`, Chu-Vandermonde
pub(crate) fn hyp2f1_poly_dd(
    a: impl Into<Dd>,
    b: impl Into<Dd>,
    c: impl Into<Dd>,
    x: f64,
) -> Option<Scaled> {
    if x == 0.0 || !x.is_finite() {
        return None;
    }
    let (a, b, c) = (a.into(), b.into(), c.into());
    // (`a >= b`, as in `hyp2f1_recur`, so the lower degree first)
    let (a, b) = if b.hi > a.hi { (b, a) } else { (a, b) };
    // (no factor of `(c)_n` is 0)
    let no_pole = |p: Dd, c: Dd| p.hi >= c.hi || !(c.hi <= 0.0 && c.is_integer());
    let is_poly =
        |p: Dd| p.is_integer() && no_pole(p, c) && (1.0..=HYP2F1_RECUR_MAXITER).contains(&-p.hi);
    let (p, q) = if is_poly(a) {
        (a, b)
    } else if is_poly(b) {
        (b, a)
    } else {
        return None;
    };
    let n = -p.hi;
    let y = if x == 1.0 {
        chu_vandermonde(n, q, c, true)
    } else {
        let direct = || poly_dd(n, q, c, Dd::new(x), true);
        let pfaff = || {
            if !(x < 1.0) {
                return None;
            }
            let w = Dd::new(x) / (Dd::new(x) - Dd::ONE);
            let y = poly_dd(n, c - q, c, w, true)?;
            Some(y.with_l(ln_dd(Dd::ONE - x) * n))
        };
        // `F(-n, q; d; 1 - x)` times `F(-n, q; c; 1)`, with `1 - x` exact in double-double
        let reversal = || {
            let d = q - c + (1.0 - n);
            let y = no_pole(p, d).then(|| poly_dd(n, q, d, Dd::ONE - x, true))??;
            let at_one = chu_vandermonde(n, q, c, false)?;
            Scaled::new(at_one.m * y.m, at_one.e + y.e, at_one.err + y.err + MACHEP)
        };
        let methods: [&dyn Fn() -> Option<Scaled>; 3] = [&direct, &pfaff, &reversal];
        best_of(methods.into_iter().map(|method| method()))
    };
    // (not if its mantissa is NaN or infinite, which `Scaled::new` accepts)
    y.filter(|y| y.m.is_finite())
}

/// `F(-n, b; c; x)` for `n >= 1` by its terms in double-double, with an error bound as in
/// `series_weighted`, plus the terms it skips if `c <= -n`. If `scale_c`, also with the mantissa
/// of `c + k`, so that the term ratio can't overflow for a subnormal or huge `c`.
fn poly_dd(n: f64, b: Dd, c: Dd, x: Dd, scale_c: bool) -> Option<Scaled> {
    // the term `t 2^te`, and the sums scaled by `2^-e`
    let (mut t, mut te, mut e) = (Dd::ONE, 0, 0);
    let (mut sum, mut sum_abs, mut k_sum_abs, mut k, mut tail) = (Dd::ONE, 1.0, 0.0, 0.0, 0.0);
    // (the mantissas of `b + k` and `x`, which avoids subnormal rounding in the update)
    let (xm, xe) = x.frexp();
    // The relative error of `x` if it was computed, as its low word can be subnormal, e.g.
    // for `1 / (1 - x)` with `x = -1.3e308`. A term `t_k` has `k` factors `x`.
    let x_err = (DENORM_MIN / x.hi.abs()).min(1.0);
    let rounding = |k: f64, sum_abs: f64, k_sum_abs: f64| {
        DD_ROUNDING * (8.0 * k_sum_abs + (k + 1.0) * sum_abs) + x_err * k_sum_abs
    };
    while k < n {
        let (bm, be) = (b + k).frexp();
        let (cm, ce) = if scale_c { (c + k).frexp() } else { (c + k, 0) };
        t = t * (bm * (k - n)) * xm / (cm * (k + 1.0));
        te += be + xe - ce;
        k += 1.0;
        if let Some(s) = rescale_exp(t.hi) {
            (t, te) = (t.ldexp(-s), te + s);
        }
        // If k times the term overflows at the scale of the sums (e.g. the first term for
        // x = 1e308), restart them at its own: `sum_abs` is at most 2^500 and `k_sum_abs` at most
        // k 2^500 (see below), so they are below 2^-503 of `|term|` and `k |term|` (as
        // k <= `HYP2F1_RECUR_MAXITER`), negligible for the error bound and the stop rule.
        if !(k * t.ldexp(te - e).hi).is_finite() {
            (sum, sum_abs, k_sum_abs, e) = (Dd::ZERO, 0.0, 0.0, te);
        }
        let term = t.ldexp(te - e);
        let term_abs = term.hi.abs();
        sum = sum + term;
        sum_abs += term_abs;
        k_sum_abs += k * term_abs;
        if !k_sum_abs.is_finite() {
            return None;
        }
        // Stop once the rest is below `thr`: for c <= -n and b + k > 0, the ratios of the next
        // terms are at most `r`, as |c + k| >= n - k and (b + k) / (k + 1) tends monotonically to
        // 1, so the rest is at most `|t_k| r / (1 - r)` (with slack for rounding). Only if the
        // rounding of `sum` is at most `|sum| / 4`, so that the full sum is at least `|sum| / 2`
        // (hence `2 tail` below).
        let bk = (b + k).hi;
        let r = x.hi.abs() * (bk / (k + 1.0)).max(1.0) * (1.0 + 4.0 * MACHEP);
        let abs = sum.hi.abs();
        let thr = DD_ROUNDING * abs;
        let negligible = term_abs.max(f64::MIN_POSITIVE)
            <= thr * (1.0 - r) / r * (1.0 - 4.0 * MACHEP)
            && rounding(k, sum_abs, k_sum_abs) <= abs / 4.0;
        if c.hi <= -n && bk > 0.0 && r < 1.0 && negligible {
            tail = thr;
            break;
        }
        if let Some(s) = rescale_exp(sum_abs) {
            sum = sum.ldexp(-s);
            [sum_abs, k_sum_abs] = [sum_abs, k_sum_abs].map(|v| ldexp_f64(v, -s));
            e += s;
        }
    }
    let err = MACHEP + (rounding(k, sum_abs, k_sum_abs) + 2.0 * tail) / sum.hi.abs();
    Scaled::new(sum.hi, e, err)
}

/// `F(-n, b; c; 1) = (c - b)_n / (c)_n` (DLMF 15.4.24) for `n >= 1` and `c` not an integer in
/// `(-n, 0]`, in double-double. If `round_subnormal`, a subnormal value is rounded with its low
/// word, so that `Scaled::to_f64` returns it without rounding the high word again.
fn chu_vandermonde(n: f64, b: Dd, c: Dd, round_subnormal: bool) -> Option<Scaled> {
    let cb = c - b;
    // the product `t 2^e`, with the mantissas of the factors as in `poly_dd`, or `None` if a
    // numerator is 0
    let product = (0..n as u32)
        .map(f64::from)
        .try_fold((Dd::ONE, 0), |(t, e), k| {
            let ((num, num_e), (den, den_e)) = ((cb + k).frexp(), (c + k).frexp());
            let t = (num.hi != 0.0).then(|| t * num / den)?;
            let e = e + num_e - den_e;
            Some(rescale_exp(t.hi).map_or((t, e), |s| (t.ldexp(-s), e + s)))
        });
    match product {
        Some((t, e)) => {
            let units =
                round_subnormal.then(|| t.ldexp(e - MIN_SUBNORMAL_EXP).round_subnormal_units());
            let m = units
                .filter(|v| v.is_subnormal())
                .map_or(t.hi, |v| ldexp_f64(v, -e));
            // (an addition, a product, and a division per factor)
            Scaled::new(m, e, MACHEP + 3.0 * n * DD_ROUNDING)
        },
        None => Some(Scaled {
            m: 0.0,
            e: 0,
            l: Dd::ZERO,
            err: 0.0,
        }),
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
        let (ca, cb, zero) = (c - a, c - b, Dd::ZERO);
        let cab = ca - b;
        let combinations = [(a, cb, zero), (cb, a, cab), (ca, b, cab), (b, ca, zero)];
        (x / (x - 1.0), Dd::new(x), combinations)
    };
    // 1 - z, which is 1 - x if x < 0, and 1 / (1 - x) else, and its log
    let (one_mz, ln_1mz) = if x < 0.0 {
        (Dd::ONE - x, ln_1mx)
    } else {
        (1.0 / (Dd::ONE - x), -ln_1mx)
    };
    // the polynomial cases of step 2
    let poly = hyp2f1_poly(combinations, c, w, ln_1mx, (one_mz, ln_1mz));
    if err_of(poly) <= HYP2F1_RECUR_TOL {
        return poly;
    }
    // the 1/x transformation (step 4), or the 1/t transformation (step 5)
    let try_recip = || recip_of_x(a, b, c, x, ln_1mx, false);
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

/// `(r, e + k)` with `exp(l) = 2^k exp(r)` and `|r| <= ln(2) / 2`
fn split_exp(l: Dd, e: i32) -> (Dd, i32) {
    // If `|k|` is so large that it's clamped, then so is `e + k`, since `|e|` is at most about 1e9
    // (2^1000 per term), and the value overflows or underflows.
    let k = (l.hi / LN_2).round().clamp(-1e12, 1e12);
    let e = (f64::from(e) + k).clamp(-1e9, 1e9) as i32;
    (l - LOG2 * k, e)
}

/// `exp(l)` rounded once from double-double, also if subnormal
pub(crate) fn exp_rounded(l: Dd) -> f64 {
    let (r, e) = split_exp(l, 0);
    let y = r.exp();
    let v = ldexp_f64(y.hi, e);
    if v.abs() <= f64::MIN_POSITIVE {
        // (in units of the smallest subnormal, as `ldexp` would round `y.hi` a second time)
        y.ldexp(e - MIN_SUBNORMAL_EXP).round_subnormal_units()
    } else {
        v
    }
}

/// `y` with `exp(l)` folded into `m 2^e` (so `l = 0`), and its roundings into the error bound
pub(crate) fn fold_exp(y: Scaled) -> Scaled {
    let (r, e) = split_exp(y.l, y.e);
    // For an unclamped `k`, `exp(r.lo)` rounds to 1. It is left out, as `r` is huge for a clamped
    // `k`, where `m` must keep its sign.
    let (m, l) = (y.m * r.hi.exp(), Dd::ZERO);
    // (the roundings of `exp` and the product, and `exp(r.lo)`)
    let err = y.err + 4.0 * MACHEP;
    Scaled { m, e, l, err }
}

/// Step 7: `F(a, b; c; x)` by recurring down in `c` from start values by `hyp2f1_recur_dd`, with
/// `ln_1mx = ln(1 - x)`, where `err` is the error bound of steps 1 to 6
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
        // (if c > max(0, q), only if steps 1 to 6 nearly succeeded, see the module docs)
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

/// `F(a, b; c; x)` for real `x < 1` with `x != 0`, or any finite `x` for a non-positive integer
/// `c` (see the module docs). Its error bound (that of `Scaled::to_f64`) is too large if it
/// exceeds `HYP2F1_RECUR_TOL`.
pub(crate) fn hyp2f1_recur(a: f64, b: f64, c: f64, x: f64) -> Option<Scaled> {
    // (beyond 1e15, the coefficients of the recurrence can overflow)
    let in_range =
        x != 0.0 && x.is_finite() && [a, b, c].iter().all(|v| (-1e15..=1e15).contains(v));
    if !in_range {
        return None;
    }
    // In a fixed order, as step 7 is not symmetric in `a` and `b` (its `n` is from `q`): the larger
    // first, which was more often accurate in our tests, e.g. for
    // `F(-104 - 9e-7, 36; 0.0074; 0.66)`, whose form `(c - a, c - b)` then has `n = 105`, not 1.
    let (a, b) = if b > a { (b, a) } else { (a, b) };
    let (a, b, c) = (Dd::new(a), Dd::new(b), Dd::new(c));
    if c.hi <= 0.0 && c.is_integer() {
        return hyp2f1_poly_dd(a, b, c, x);
    }
    if x >= 1.0 {
        return None;
    }
    // (`log1p` of `DoubleDouble` is only as accurate as that of `f64`, but `ln_dd` is accurate)
    let ln_1mx = ln_dd(Dd::ONE - x);
    let mut best = hyp2f1_recur_dd(a, b, c, x, ln_1mx);
    // step 6
    if !(err_of(best) <= HYP2F1_RECUR_TOL) && c.hi < 0.0 {
        best = min_err(best, recip_of_x(a, b, c, x, ln_1mx, true));
    }
    let err = err_of(best);
    if err <= HYP2F1_RECUR_TOL {
        return best;
    }
    best = min_err(best, hyp2f1_recur_shifted(a, b, c, x, ln_1mx, err));
    if err_of(best) <= HYP2F1_RECUR_TOL {
        return best;
    }
    // Step 8, with `eps` the distance of `b - a` (or `c - a - b` if `x > 0`) to an integer. Last,
    // as its error bound is heuristic.
    let m = if x < 0.0 { b - a } else { c - a - b };
    let eps = (m - m.hi.round()).hi;
    if eps.abs() <= HYP2F1_RECUR_NEAR_INT {
        let shift = if x < 0.0 { eps } else { -eps };
        best = min_err(best, recip_interp(a, b, c, x, ln_1mx, shift));
    }
    best
}

#[cfg(test)]
mod tests {
    use super::{
        Dd, HYP2F1_RECUR_TOL, Scaled, coefs_c, exp_rounded, gamma_quot, hyp2f1_recur, ln_gamma,
        psi_dd, recur_down,
    };
    use crate::utils::DENORM_MIN;

    #[test]
    fn test_ln_gamma_near_pole() {
        // Arb; `-n + e` with `e = +-1e-17` only in the low word
        let ln_gamma_at =
            |n: f64, e| ln_gamma(Dd { hi: -n, lo: e }).map(|(v, sign, _)| (v.hi, sign));
        assert_eq!(ln_gamma_at(1.0, 1e-17), Some((39.143_946_580_898_78, -1.0)));
        assert_eq!(ln_gamma_at(1.0, -1e-17), Some((39.143_946_580_898_78, 1.0)));
        assert_eq!(ln_gamma_at(2.0, 1e-17), Some((38.450_799_400_338_83, 1.0)));
        assert_eq!(ln_gamma_at(1.0, 0.0), None);
    }

    #[test]
    fn test_ln_gamma_non_finite() {
        for (hi, lo) in [(1.0, f64::NAN), (1.0, f64::INFINITY), (f64::NAN, 0.0)] {
            assert!(ln_gamma(Dd { hi, lo }).is_none(), "{hi:e} + {lo:e}");
        }
    }

    #[test]
    fn test_gamma_quot() {
        // Arb; the last two are next to a pole
        for (v, e, desired) in [
            (2.5, 1e-3, 0.703_649_225_211_809_1_f64),
            (1.0, -2.5e-4, -0.577_462_985_633_150_5),
            (1000.0, 2.5e-4, 6.913_222_526_976_623),
            (0.24, 1e-9, -4.406_160_326_799_55),
            (-7.3, -3.3e-5, 4.336_750_268_419_202),
            (0.98, 0.04, -0.570_759_743_970_275_3),
            (13.5, -0.04, 2.436_589_945_228_180_4),
            (-7.3, 0.02, 4.699_610_668_760_759),
            (-135.999_997_964, -2.0344e-6, -624_989_001.899_548_8),
            (-399.999_999_998_5, 1e-12, -666_229_553.282_955_2),
        ] {
            let (y, err) = gamma_quot(Dd::new(v), Dd::new(e)).unwrap();
            let tol = 1e-14 * desired.abs().max(1.0);
            assert!((y - desired).abs() <= err && err <= tol, "{y} {err}");
        }
        // psi(v) at a subnormal e, whatever the branch, and to within the bound at the smallest
        // normal ones, where e / (v + i) is subnormal
        let tiny = DENORM_MIN;
        for (v, e, desired) in [
            (1.0, tiny, -0.577_215_664_901_532_9),
            (-7.3, -tiny, 4.337_307_305_510_05),
            (0.3, 4.0 * f64::MIN_POSITIVE, -3.502_524_222_200_133),
            (12.5, -f64::MIN_POSITIVE, 2.485_195_651_274_912_3),
        ] {
            let (y, err) = gamma_quot(Dd::new(v), Dd::new(e)).unwrap();
            assert!((y - desired).abs() <= err && err <= 1e-14, "{y} {err}");
        }
        // psi at e = 0, and None at a pole of Gamma(v) or Gamma(v + e)
        let (one, pole, e) = (Dd::ONE, Dd::new(-3.0), Dd::new(1.0 / 4096.0));
        assert_eq!(gamma_quot(one, Dd::ZERO), psi_dd(one));
        assert_eq!(gamma_quot(pole, e), None);
        assert_eq!(gamma_quot(pole + e, -e), None);
    }

    #[test]
    fn test_recur_down_start_errors() {
        // y_0 = y_2 / 8 (as B_1 = 0), so the error of y_2 carries over, although u_1 = 0
        let y = |m, err| Scaled::new(m, 0, err).unwrap();
        let (p, q, c, z) = (Dd::new(3.5), Dd::new(1.5), Dd::ONE, Dd::new(-0.5));
        let y_0 = recur_down(1.0, y(0.5, 0.0), y(0.8, 1e-6), coefs_c(p, q, c, z)).unwrap();
        assert!((y_0.m - 0.1).abs() < 1e-15 && y_0.err >= 1e-6, "{y_0:?}");
    }

    #[test]
    fn test_hyp2f1_recur_interp() {
        // Arb; step 8 at an integer `b - a` and `c - a - b`, and at 5e-4 from one: within its error
        // bound
        for (a, b, c, x, desired) in [
            (-130.5, -106.5, 88.75, -56.25, -7.673_134_033_185_444e147),
            (-209.25, 294.875, 145.625, 0.955, -1.377_551_721_883_881_6e-135),
            (-129.875, -111.875_5, 121.5, -71.875, -2.123_783_401_591_176_3e153),
            (-294.75, 181.0, -146.750_5, 0.984, -9.089_494_319_571_015e122),
        ] {
            let (y, err) = hyp2f1_recur(a, b, c, x).map(Scaled::to_f64).unwrap();
            assert!((y / desired - 1.0).abs() <= err && err <= HYP2F1_RECUR_TOL, "{y} {err}");
        }
    }

    #[test]
    fn test_exp_rounded() {
        // (2^52 - 1/2 - 2^-8) times the smallest subnormal (Arb), which `ldexp` alone rounds up to
        // the smallest normal
        let l = Dd {
            hi: -708.396_418_532_264_1,
            lo: -2.758_730_638_543_521_8e-14,
        };
        assert_eq!(exp_rounded(l), f64::MIN_POSITIVE - DENORM_MIN);
        assert_eq!(exp_rounded(Dd::new(-1e300)), 0.0);
        assert_eq!(exp_rounded(Dd::new(1e300)), f64::INFINITY);
    }

    #[test]
    fn test_hyp2f1_recur_connect() {
        // Arb; step 6 for x < 0 and 0 < x < 1: within its error bound
        for (a, b, c, x, desired) in [
            (-430.375, -193.75, -651.5, -18.5, 9.891_276_940_572_099e225),
            (786.875, 195.75, -474.75, 0.05, 2.450_524_293_296_754_5e-7),
        ] {
            let (y, err) = hyp2f1_recur(a, b, c, x).map(Scaled::to_f64).unwrap();
            assert!((y / desired - 1.0).abs() <= err && err <= HYP2F1_RECUR_TOL, "{y} {err}");
        }
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
            let actual = hyp2f1_recur(a, b, c, x).map(Scaled::to_f64);
            let actual = actual.filter(|&(_, err)| err <= HYP2F1_RECUR_TOL);
            assert!(
                actual.is_none_or(|(y, _)| (y / desired - 1.0).abs() < 1e-10),
                "{actual:?} != {desired}"
            );
        }
    }
}
