//! Translated into pure Rust from `xsf/cephes/hyp2f1.h` (scipy/xsf@9a1dfd1).
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
//!
//! Unlike in xsf:
//!
//! - The `loss` estimates are returned with the results, instead of written to output references.
//! - `set_error` is a no-op in our build of xsf, so its calls are omitted, `goto hypdon` and
//!   `goto hypdiv` become plain returns, and the never-read assignments to `err` and `q` in
//!   `hyp2f1` are dropped.
//!
//! Our `psi` is more accurate than xsf's just above negative integers, so the psi function
//! expansion (for `x > 0.9` and integer `c - a - b`) can differ from xsf when `a`, `b`, `c - a` or
//! `c - b` is close to a negative integer.

use crate::ext::{
    HYP2F1_RECUR_MAXITER, HYP2F1_RECUR_NEAR_INT, HYP2F1_RECUR_TOL, Scaled, hyp2f1_recur, ln_dd,
};
use crate::utils::{is_nonpos_int, pow2};
use crate::xsf::cephes::consts::MACHEP;
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, frexp_f64, ldexp_f64, two_prod, two_sum};
use crate::xsf::cephes::{gamma, gammasgn, lgam_ratio, lgam_sgn, psi, rgamma};

pub(crate) const HYP2F1_EPS: f64 = 1.0e-13;

pub(crate) const HYP2F1_ETHRESH: f64 = 1.0e-12;
/// The largest heuristic error estimate with which the 1/x transformation (DLMF 15.8.2) replaces
/// an accurate series of the Pfaff transformations, as that is often more accurate than its own
/// estimate (empirical)
const HYP2F1_RECIP_MAX_ERR: f64 = 1.0e-13;
/// The estimated relative error from which `hyp2f1` returns NaN: if a sum is only rounding error,
/// then its estimate is about 1/2 to 1.
const HYP2F1_LOST: f64 = 0.5;
const HYP2F1_MAXITER: u64 = 10000;
/// About `sqrt(f64::MAX)`, to rescale tiny terms exactly, far from overflow and underflow
const HYP2F1_SCALE: f64 = pow2(512);

/// Whether `t` is a non-positive integer, up to `HYP2F1_EPS` (as in cephes)
fn near_neg_int(t: f64) -> bool {
    t <= 0.0 && (t - t.round()).abs() < HYP2F1_EPS
}

/// Unlike in xsf (scipy/xsf#371), `hys2f1`, or the power series and `true` if the recurrence in
/// `a` has a small `loss` but a large error estimate, and the series is accurate and confirms the
/// discrepancy (see `less_accurate`). Not if `a` or `b` is near a non-positive integer, where the
/// series can stop too early.
fn hys2f1_checked(a: f64, b: f64, c: f64, x: f64) -> ((f64, f64, f64), bool) {
    let (y, loss, err) = hys2f1(a, b, c, x);
    if loss < HYP2F1_ETHRESH && err >= HYP2F1_ETHRESH && !(near_neg_int(a) || near_neg_int(b)) {
        let series = hys2f1_series(a, b, c, x);
        if series.1 < HYP2F1_ETHRESH && less_accurate(y, series.0, series.1) {
            return (series, true);
        }
    }
    ((y, loss, err), false)
}

/// Whether `y` is less accurate than `alt` with the estimated relative error `err`, by at least
/// `HYP2F1_ETHRESH`: if they differ by more than `2 err`, then the error of `y` exceeds `err`.
fn less_accurate(y: f64, alt: f64, err: f64) -> bool {
    let diff = (y - alt).abs();
    diff > 2.0 * err * alt.abs() && diff > HYP2F1_ETHRESH * alt.abs()
}

/// Whether `t` is `n` up to `HYP2F1_EPS`. Unlike in xsf, relative to `n` if `|n| < 1`, so that a
/// tiny `c`, `c - a` or `c - b` isn't 0.
fn near(t: f64, n: f64) -> bool {
    (t - n).abs() <= HYP2F1_EPS * n.abs().min(1.0)
}

/// `hys2f1` and its error estimate, or if that is large, `hys2f1_series` and its `loss` if
/// smaller (not if `a` or `b` is near a non-positive integer, as in `hys2f1_checked`)
fn hys2f1_min(a: f64, b: f64, c: f64, x: f64) -> (f64, f64) {
    let (y, _, err) = hys2f1(a, b, c, x);
    if err < HYP2F1_ETHRESH || near_neg_int(a) || near_neg_int(b) {
        return (y, err);
    }
    let (series, loss, _) = hys2f1_series(a, b, c, x);
    if loss < err { (series, loss) } else { (y, err) }
}

/// Whether `c` is near, but not at, a non-positive integer: unlike in xsf, 2F1 has no pole there
fn near_pole(c: f64) -> bool {
    c < 0.0 && c != c.round() && near(c, c.round())
}

/// The series of the Pfaff transformations (DLMF 15.8.1) in `x / (x - 1)`, by `hys2f1_min`, with
/// the smaller error estimate. Unlike in xsf, for c `near_pole`, that includes the rounding of
/// `c - b` (or `c - a`) relative to its distance to an integer, which the pole of `1 / (c)_k`
/// near it amplifies (e.g. for a tiny `b`).
fn hys2f1_pfaff(a: f64, b: f64, c: f64, x: f64) -> (f64, f64) {
    let rounding = |t: f64| {
        let (d, e) = two_sum(c, -t);
        if e == 0.0 {
            0.0
        } else {
            (e / (d - d.round())).abs()
        }
    };
    let w = x / (x - 1.0);
    let (ya, mut ea) = hys2f1_min(a, c - b, c, w);
    let (yb, mut eb) = hys2f1_min(b, c - a, c, w);
    if near_pole(c) {
        (ea, eb) = (ea + rounding(b), eb + rounding(a));
    }
    // if one of these overflows, its error estimate is NaN, so use the other
    let (y, err, p) = if ea <= eb || eb.is_nan() {
        (ya, ea, a)
    } else {
        (yb, eb, b)
    };
    // Unlike in xsf, for c `near_pole`, split the power factor if it isn't normal, as the product
    // can still be (as for the Euler transformation in `hyp2f1_impl`)
    let s = 1.0 - x;
    let f = s.powf(-p);
    if f.is_normal() || !near_pole(c) {
        (f * y, err)
    } else {
        let h = s.powf(-p / 2.0);
        (h * y * h, err)
    }
}

/// Defining power series expansion of Gauss hypergeometric function.
/// The `loss` (the second return value) estimates loss of significance. Unlike in xsf, the third
/// is an estimate of the relative error that also includes the error of the recurrence in `a`.
fn hys2f1(mut a: f64, mut b: f64, c: f64, x: f64) -> (f64, f64, f64) {
    let mut f;
    let mut intflag = false;

    if b.abs() > a.abs() {
        /* Ensure that |a| > |b| ... */
        f = b;
        b = a;
        a = f;
    }

    let ib = b.round() as i32;

    if (b - f64::from(ib)).abs() < HYP2F1_EPS
        && ib <= 0
        && (b.abs() < a.abs() || (b == -a && a <= HYP2F1_MAXITER as f64))
    {
        /* .. except when `b` is a negative integer with |b| <= |a| */
        f = b;
        b = a;
        a = f;
        intflag = a.abs() != b.abs(); /* not for the tie `b == -a` */
    }

    if (a.abs() > c.abs() + 1.0 || intflag) && (c - a).abs() > 2.0 && a.abs() > 2.0 {
        /* |a| >> |c| implies that large cancellation error is to be expected.
         *
         * We try to reduce it with the recurrence relations
         */
        return hyp2f1ra(a, b, c, x);
    }
    hys2f1_series(a, b, c, x)
}

/// The power series of `hys2f1`, without the recurrence in `a`
fn hys2f1_series(a: f64, b: f64, c: f64, x: f64) -> (f64, f64, f64) {
    let mut i: u64 = 0;
    let mut umax = 0.0;
    let f = a;
    let g = b;
    let h = c;
    let mut s = 1.0;
    let mut u = 1.0;
    let mut k = 0.0;
    // Unlike in xsf, if x > 0 and c < 0 is not an integer (even if near one), don't stop before
    // the terms have passed k = -c and stopped growing, since they can be tiny before and huge
    // after it.
    let pole_ahead = x > 0.0 && h < 0.0 && h.fract() != 0.0;
    // up to `HYP2F1_MAXITER` terms before and after k = -c
    let maxiter = HYP2F1_MAXITER * if pole_ahead { 2 } else { 1 };
    // the terms are scaled by `HYP2F1_SCALE^n_scaled` so that they don't underflow, and are only
    // summed (and counted in `loss`) if `n_scaled == 0`, as they are negligible otherwise
    let (mut n_scaled, mut i_scaled) = (0, 0);
    // Unlike in xsf, for c `near_pole`, `loss` weighs each term by its index instead of counting
    // the largest one once, as its rounding accumulates over the term ratios, and a term after
    // the pole can be the largest (times MACHEP before summing, so that it doesn't overflow)
    let (weigh, mut weighted) = (near_pole(h), 0.0);
    loop {
        // Unlike in xsf, only if c = 0, since the series is finite for a tiny c != 0
        if h == 0.0 {
            return (f64::INFINITY, 1.0, 1.0);
        }
        let m = k + 1.0;
        let r = (f + k) * (g + k) * x / ((h + k) * m);
        u *= r;
        // fail if the term ratio underflows, e.g. for tiny a and b
        if pole_ahead && r.abs() < f64::MIN_POSITIVE && f + k != 0.0 && g + k != 0.0 {
            return (s, 1.0, 1.0);
        }
        if pole_ahead && u != 0.0 && u.abs() < 1.0 / HYP2F1_SCALE {
            (u, n_scaled) = (u * HYP2F1_SCALE, n_scaled + 1);
        } else if n_scaled > 0 && u.abs() > 1.0 {
            (u, n_scaled) = (u / HYP2F1_SCALE, n_scaled - 1);
        }
        let t = if n_scaled == 0 { u } else { 0.0 };
        s += t;
        if weigh {
            weighted += MACHEP * m * t.abs();
        }
        k = t.abs(); /* remember largest term summed */
        if k > umax {
            umax = k;
        }
        k = m;
        i += 1;
        i_scaled += u64::from(n_scaled > 0);
        if i > maxiter {
            /* should never happen */
            return (s, 1.0, 1.0);
        }
        let growing = pole_ahead && u != 0.0 && (k < 1.0 - h || r.abs() > 1.0);
        // Unlike in xsf (scipy/xsf#470), a negligible term doesn't end the sum if a later one
        // isn't, e.g. a tiny first term from a tiny `a`, after which the terms can grow again.
        let grows_back = || (t / s).abs() * ln_max_growth(f, g, h, x, k).exp() > MACHEP;
        if !(s == 0.0 || (t / s).abs() > MACHEP || growing || grows_back()) {
            break;
        }
    }

    /* return estimated relative error */
    let lost = if weigh { weighted } else { MACHEP * umax };
    let loss = lost / s.abs() + (MACHEP * (i - i_scaled) as f64);

    (s, loss, loss)
}

/// The largest `ln(t_m / t_n)` (`m >= n`) of the terms `t_m` of the power series, if their ratios
/// `x (a + m) (b + m) / ((c + m) (m + 1))` are positive, from `0 < x < 1` and `a + n`, `b + n`,
/// `c + n > 0`, else 0. The ratios are at least 1 between the roots of a quadratic in `m`, so the
/// largest term is at its larger root (rounded up). Not for alternating ratios, as summing
/// growing terms then cancels.
fn ln_max_growth(a: f64, b: f64, c: f64, x: f64, n: f64) -> f64 {
    if !(x > 0.0 && x < 1.0 && a + n > 0.0 && b + n > 0.0 && c + n > 0.0) {
        return 0.0;
    }
    // x (a + m) (b + m) - (c + m) (m + 1) = p m^2 + q m + r, with p < 0
    let (p, q, r) = (x - 1.0, x * (a + b) - c - 1.0, x * a * b - c);
    let m = ((-q - (q * q - 4.0 * p * r).sqrt()) / (2.0 * p)).ceil();
    if !(m > n) {
        return 0.0;
    }
    // the sum of `ln((u + j) / (v + j))` for `n <= j < m`, with the smaller of `m - n` and
    // `|u - v|` as the (non-negative) second argument of `lgam_ratio`, to bound the cancellation
    let ln_ratio = |u: f64, v: f64| {
        let (u, v, sign) = if u >= v { (u, v, 1.0) } else { (v, u, -1.0) };
        sign * if m - n < u - v {
            lgam_ratio(u + n, m - n) - lgam_ratio(v + n, m - n)
        } else {
            lgam_ratio(v + m, u - v) - lgam_ratio(v + n, u - v)
        }
    };
    (m - n) * x.ln() + ln_ratio(a, 1.0) + ln_ratio(b, c)
}

/// Apply transformations for |x| near 1 then call the power series. Unlike in xsf, the third
/// return value is an error estimate as in `hys2f1`.
fn hyt2f1(a: f64, b: f64, c: f64, x: f64) -> (f64, f64, f64) {
    let (mut p, mut q, mut r, mut t, mut w, d);
    let (d1, d2, e, mut y1);
    let mut sign;

    let mut neg_int_a = false;
    let mut neg_int_b = false;

    let ia = a.round() as i32;
    let ib = b.round() as i32;

    if a <= 0.0 && (a - f64::from(ia)).abs() < HYP2F1_EPS {
        /* a is a negative integer */
        neg_int_a = true;
    }

    if b <= 0.0 && (b - f64::from(ib)).abs() < HYP2F1_EPS {
        /* b is a negative integer */
        neg_int_b = true;
    }

    let (mut err, err_full);
    let s = 1.0 - x;
    let y = 'done: {
        if x < -0.5 && !(neg_int_a || neg_int_b) {
            let y;
            if b > a {
                (y, err, err_full) = hys2f1(a, c - b, c, -x / s);
                break 'done s.powf(-a) * y;
            }
            (y, err, err_full) = hys2f1(c - a, b, c, -x / s);
            break 'done s.powf(-b) * y;
        }

        d = c - a - b;
        let id = d.round(); /* nearest integer to d */

        if x > 0.9 && !(neg_int_a || neg_int_b) {
            /* Try the power series first */
            // Unlike in xsf, also if c-a-b is an integer, as the psi function expansion below
            // can be inexact or fail.
            //
            // Unlike in xsf (scipy/xsf#371), if its `loss` is small but its error estimate incl.
            // the recurrence in `a` isn't, it is only replaced by AMS55 #15.3.6 if that shows it
            // to be less accurate.
            let ((y0, loss0, full0), _) = hys2f1_checked(a, b, c, x);
            if loss0 < HYP2F1_ETHRESH && full0 < HYP2F1_ETHRESH {
                (err, err_full) = (loss0, full0);
                break 'done y0;
            }

            // Unlike in xsf (scipy/xsf#461), use the psi function expansion if c-a-b is an
            // integer up to `HYP2F1_EPS`, not `MACHEP`, as AMS55 #15.3.6 cancels near one, and
            // only up to `HYP2F1_MAXITER`, beyond which its finite sum can hang
            if (d - id).abs() > HYP2F1_EPS || id.abs() > HYP2F1_MAXITER as f64 {
                let (mut lg, mut sgngam);
                let (err1, q_full, r_full);

                /* test for integer c-a-b */
                /* If power series fails, then apply AMS55 #15.3.6 */
                ((q, err, q_full), _) = hys2f1_checked(a, b, 1.0 - d, s);
                sign = 1;
                (w, sgngam) = lgam_sgn(d);
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(c - a);
                w -= lg;
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(c - b);
                w -= lg;
                sign *= sgngam;
                let (w_q, sign_q) = (w, f64::from(sign));
                ((r, err1, r_full), _) = hys2f1_checked(c - a, c - b, d + 1.0, s);
                sign = 1;
                (w, sgngam) = lgam_sgn(-d);
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(a);
                w -= lg;
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(b);
                w -= lg;
                sign *= sgngam;
                // Unlike in xsf (scipy/xsf#473), if a factor leaves the normal range, combine
                // `Gamma(c)` and `s^d` with the logarithms of the others, as the terms can still be
                // normal (e.g. for a tiny `b`). Not at a pole of `Gamma(+-d)`, nor if an inner
                // series is inaccurate.
                let (mut exp_q, mut s_d, mut exp_r, mut gamma_c) =
                    (w_q.exp(), s.powf(d), w.exp(), gamma(c));
                let in_range = [exp_q, s_d, exp_r, gamma_c].into_iter().all(f64::is_normal);
                let accurate = q_full < HYP2F1_ETHRESH && r_full < HYP2F1_ETHRESH;
                if !in_range && w_q.is_finite() && w.is_finite() && accurate {
                    // `s^d` and `|Gamma(c)|` are moved into the exponents
                    let (lg_c, sign_c) = lgam_sgn(c);
                    exp_q = (w_q + lg_c).exp();
                    (s_d, exp_r) = (1.0, (w + d * s.ln() + lg_c).exp());
                    gamma_c = f64::from(sign_c);
                }
                q *= sign_q * exp_q;
                r = r * s_d * (f64::from(sign) * exp_r);
                let mut y = q + r;

                // the error of `y` from `q` and `r` (incl. their recurrences in `a`), their sum,
                // and `Gamma(+-d)`: its poles amplify the rounding of `d` by `|a|+|b|+|c|` to
                // `|d - id|`
                let near_pole = 1.0 + (a.abs() + b.abs() + c.abs()) / (d - id).abs();
                let full = (q.abs() * q_full
                    + r.abs() * r_full
                    + MACHEP * near_pole * q.abs().max(r.abs()))
                    / y.abs();

                q = q.abs(); /* estimate cancellation error */
                r = r.abs();
                if q > r {
                    r = q;
                }
                err += err1 + (MACHEP * r) / y;

                y *= gamma_c;
                if loss0 < HYP2F1_ETHRESH && !less_accurate(y0, y, full) {
                    (err, err_full) = (loss0, full0);
                    break 'done y0;
                }
                err_full = full;
                break 'done y;
            }
            /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12
             *
             * Although AMS55 does not explicitly state it, this expansion fails
             * for negative integer a or b, since the psi and Gamma functions
             * involved have poles.
             */

            let aid;
            if id >= 0.0 {
                e = d;
                d1 = d;
                d2 = 0.0;
                aid = id as i32;
            } else {
                e = -d;
                d1 = 0.0;
                d2 = d;
                aid = -id as i32;
            }

            let ax = s.ln();

            /* sum for t = 0 */
            let mut y = psi(1.0) + psi(1.0 + e) - psi(a + d1) - psi(b + d1) - ax;
            y *= rgamma(e + 1.0);
            // Unlike in xsf, estimate the errors of both sums and their cancellation (scaled
            // as the sums), as the `t`-th term is a product of about `t` roundings
            let (mut y_err, mut y1_err) = (y.abs(), 1.0_f64);
            // the relative error of `y` from the absolute error `e`, but no digit
            // (`HYP2F1_LOST`) if `y` is 0 or either isn't finite, e.g. as the factors below
            // over- or underflow
            let rel_err = |e: f64, y: f64| {
                let err = MACHEP * e / y.abs();
                if err.is_finite() && y.is_finite() && y != 0.0 {
                    err
                } else {
                    HYP2F1_LOST
                }
            };

            p = (a + d1) * (b + d1) * s * rgamma(e + 2.0); /* Poch for t=1 */
            t = 1.0;
            loop {
                r = psi(1.0 + t) + psi(1.0 + t + e) - psi(a + t + d1) - psi(b + t + d1) - ax;
                q = p * r;
                y += q;
                y_err += t * q.abs();
                p *= s * (a + t + d1) / (t + 1.0);
                p *= (b + t + d1) / (t + 1.0 + e);
                t += 1.0;
                if t > HYP2F1_MAXITER as f64 {
                    /* should never happen */
                    (err, err_full) = (1.0, 1.0);
                    break 'done f64::NAN;
                }
                if !(y == 0.0 || (q / y).abs() > HYP2F1_EPS) {
                    break;
                }
            }

            if id == 0.0 {
                let f = gamma(c) / (gamma(a) * gamma(b));
                y *= f;
                err = rel_err(y_err * f.abs(), y);
                err_full = err;
                break 'done y; // goto psidon
            }

            y1 = 1.0;

            'nosum: {
                if aid == 1 {
                    break 'nosum;
                }

                t = 0.0;
                p = 1.0;
                for _ in 1..aid {
                    r = 1.0 - e + t;
                    p *= s * (a + t + d2) * (b + t + d2) / r;
                    t += 1.0;
                    p /= t;
                    y1 += p;
                    y1_err += t * p.abs();
                }
            }
            // nosum:
            p = gamma(c);
            let f1 = gamma(e) * p * (rgamma(a + d1) * rgamma(b + d1));
            let f = p * (rgamma(a + d2) * rgamma(b + d2));
            // `y` or `y1` times `s^id`
            let q = (s.powf(id), id * ax); /* s to the id power */
            let one = (1.0, 0.0);
            let (q, q1) = if id > 0.0 { (q, one) } else { (one, q) };
            (y1, y1_err) = times_gamma_ratio((y1, y1_err), f1, q1, [e, c], [a + d1, b + d1]);
            (y, y_err) = times_gamma_ratio((y, y_err), f, q, [1.0, c], [a + d2, b + d2]);
            if (aid & 1) != 0 {
                y = -y;
            }

            y += y1;
            err = rel_err(y_err + y1_err, y);
            err_full = err;
            // psidon:
            break 'done y;
        }

        /* Use defining power series if no special cases */
        // Unlike in xsf (scipy/xsf#371), check the recurrence in `a` here too
        let y;
        ((y, err, err_full), _) = hys2f1_checked(a, b, c, x);
        y
    };
    // done:
    (y, err, err_full)
}

/// `y f q` and its error in `MACHEP`, from `y` and its error, the ratio
/// `f = Gamma(num_0) Gamma(num_1) / (Gamma(den_0) Gamma(den_1))`, and `(q, ln(q))`
fn times_gamma_ratio(
    (y, y_err): (f64, f64),
    f: f64,
    (q, ln_q): (f64, f64),
    num: [f64; 2],
    den: [f64; 2],
) -> (f64, f64) {
    let yf = y * f;
    let fy = yf * q;
    if f.is_normal() && q.is_normal() && (y == 0.0 || (yf.is_normal() && fy.is_normal())) {
        return (fy, y_err * f.abs() * q);
    }
    // Unlike in xsf, by logarithms if a factor or a product leaves the normal range: e.g. `f` can
    // underflow to 0 while `y f` is the larger part of the result.
    let (mut ln_f, mut ln_f_err, mut sign) = (ln_q, ln_q.abs(), y.signum());
    for (x, power) in [(num[0], 1.0), (num[1], 1.0), (den[0], -1.0), (den[1], -1.0)] {
        let (lg, lg_sign) = lgam_sgn(x);
        (ln_f, ln_f_err, sign) =
            (ln_f + power * lg, ln_f_err + lg.abs(), sign * f64::from(lg_sign));
    }
    if ln_f == f64::NEG_INFINITY {
        // a pole of the denominator: exactly 0
        return (0.0, 0.0);
    }
    let ln_y = y.abs().ln();
    let fy = sign * (ln_f + ln_y).exp();
    // the rounding of the exponent, which also covers that of `fy` if subnormal (as it is > 700)
    let rounding = (ln_f_err + ln_y.abs()) * fy.abs().max(f64::MIN_POSITIVE);
    (fy, (ln_f + y_err.ln()).exp() + rounding)
}

/*
  15.4.2 Abramowitz & Stegun.
*/
fn hyp2f1_neg_c_equal_bc(a: f64, b: f64, x: f64) -> f64 {
    let mut collector = 1.0;
    let mut sum = 1.0;
    let mut collector_max = 1.0_f64;

    if !(b.abs() < 1e5) {
        return f64::NAN;
    }

    let mut k = 1.0;
    while k <= -b {
        collector *= (a + k - 1.0) * x / k;
        collector_max = collector.abs().max(collector_max);
        sum += collector;
        k += 1.0;
    }

    if 1e-16 * (1.0 + collector_max / sum.abs()) > 1e-7 {
        return f64::NAN;
    }

    sum
}

/// `(1 - x)^m`, by `ln(1 - x)` in double-double for x < 1, where it's finite, as `1 - x` can round
fn pow_1mx(x: f64, m: f64) -> f64 {
    // (for a tiny x by its series, as `ln_dd` is accurate only to about 2^-106 absolute)
    let ln_1mx = if x.abs() < pow2(-27) {
        let x = Dd::new(x);
        -(x + x * x * (x / 3.0 + 0.5))
    } else {
        ln_dd(Dd::new(1.0) - x)
    };
    let l = ln_1mx * m;
    let y = Scaled::new(1.0, 0, 0.0).filter(|_| l.is_finite());
    y.map_or_else(|| (1.0 - x).powf(m), |y| y.with_l(l).to_f64().0)
}

/// Whether `sum_{k <= m} (b)_k x^k / k!` (DLMF 15.2.5) is exactly 0, as are all its terms and
/// partial sums in `f64` (for `m <= HYP2F1_RECUR_MAXITER`)
fn sums_to_zero(b: f64, m: f64, x: f64) -> bool {
    // `p q` if exact: `two_prod` gives its rounding error exactly if `|p q| >= 2^-969`
    // (2^53 times the smallest normal), else only if `p` or `q` is 0
    let exact_prod = |p: f64, q: f64| {
        let (r, err) = two_prod(p, q);
        (err == 0.0 && (r.abs() >= const { pow2(53 - 1022) } || p == 0.0 || q == 0.0)).then_some(r)
    };
    let (mut t, mut sum, mut k) = (1.0, 1.0, 0.0);
    while k < m && m <= HYP2F1_RECUR_MAXITER {
        let (u, u_err) = two_sum(b, k);
        let Some(q) = exact_prod(t, u).and_then(|p| exact_prod(p, x)) else {
            return false;
        };
        let next = q / (k + 1.0);
        let (s, s_err) = two_sum(sum, next);
        // (`next` is exact iff `next (k + 1) = q` exactly)
        if u_err != 0.0 || s_err != 0.0 || exact_prod(next, k + 1.0) != Some(q) {
            return false;
        }
        (t, sum, k) = (next, s, k + 1.0);
    }
    sum == 0.0
}

/// Evaluate hypergeometric function by two-term recurrence in `a`.
///
/// This avoids some of the loss of precision in the strongly alternating
/// hypergeometric series, and can be used to reduce the `a` and `b` parameters
/// to smaller values.
///
/// AMS55 #15.2.10
fn hyp2f1ra(a: f64, b: f64, c: f64, x: f64) -> (f64, f64, f64) {
    let mut f2;
    // the coefficients `(A, B)` of a step `f0 = A f1 + B f2`, and its rounding error in MACHEP
    let mut steps = Vec::new();
    let step = |t: f64, q: f64, den: f64, f1: f64, f2: f64| {
        let p = 2.0 * t - c - t * x + b * x;
        let p_abs = (2.0 * t).abs() + c.abs() + (t * x).abs() + (b * x).abs();
        let e = (p_abs * f1.abs() + (q * f2).abs()) / den.abs();
        (-p / den, -q / den, e)
    };

    /* Don't cross c or zero */
    let da = if (c < 0.0 && a <= c) || (c >= 0.0 && a >= c) {
        (a - c).round()
    } else {
        a.round()
    };
    let mut t = a - da;

    debug_assert!(da != 0.0);

    if da.abs() > HYP2F1_MAXITER as f64 {
        /* Too expensive to compute this value, so give up */
        return (f64::NAN, 1.0, 1.0);
    }

    // Unlike in xsf (scipy/xsf#371), check the starting values, but only replace both, as their
    // errors correlate
    let t0 = if da < 0.0 { t - 1.0 } else { t + 1.0 };
    let ((mut start_f1, checked1), (mut start_f0, checked0)) =
        (hys2f1_checked(t, b, c, x), hys2f1_checked(t0, b, c, x));
    if checked1 != checked0 {
        (start_f1, start_f0) = (hys2f1(t, b, c, x), hys2f1(t0, b, c, x));
    }
    let ((mut f1, loss1, _), (mut f0, loss0, _)) = (start_f1, start_f0);
    let loss = loss1 + loss0;

    if da < 0.0 {
        /* Recurse down */
        t -= 1.0;
        let mut n = 1.0;
        while n < -da {
            f2 = f1;
            f1 = f0;
            f0 = -(2.0 * t - c - t * x + b * x) / (c - t) * f1 - t * (x - 1.0) / (c - t) * f2;
            steps.push(step(t, t * (x - 1.0), c - t, f1, f2));
            t -= 1.0;
            n += 1.0;
        }
    } else {
        /* Recurse up */
        t += 1.0;
        let mut n = 1.0;
        while n < da {
            f2 = f1;
            f1 = f0;
            f0 = -((2.0 * t - c - t * x + b * x) * f1 + (c - t) * f2) / (t * (x - 1.0));
            steps.push(step(t, c - t, t * (x - 1.0), f1, f2));
            t += 1.0;
            n += 1.0;
        }
    }

    // Unlike in xsf, also estimate the error of the result, which is large if the recurrence is
    // unstable:
    //
    // - the error of each starting value (incl. that of its own recurrence) and of each step is
    //   weighted by its effect on `f0`, separately, since `f0` can have a large relative but a
    //   tiny absolute error
    // - that is `(g, h) = (d f0 / d f1, d f0 / d f2)`, from the last step backwards
    let (mut g, mut h, mut rounding) = (1.0_f64, 0.0, 0.0);
    for &(a_k, b_k, e) in steps.iter().rev() {
        rounding += e * g.abs();
        (g, h) = (a_k * g + h, b_k * g);
    }
    // `e |w s|` for a weight `w`, but inf (not NaN) if `e` is inf and `s` is 0
    let weigh = |w: f64, (s, _, e): (f64, f64, f64)| {
        if e.is_infinite() {
            e
        } else {
            e * (w * s).abs()
        }
    };
    let start = weigh(g, start_f0) + weigh(h, start_f1);
    (f0, loss, (start + MACHEP * rounding) / f0.abs())
}

/// ```text
///      Gauss hypergeometric function   F
///                                     2 1
/// ```
///
/// ```text
///  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
///                           2 1
///
///           inf.
///            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
///   =  1 +   >   -----------------------------  x   .
///            -         c(c+1)...(c+k) (k+1)!
///          k = 0
/// ```
///
/// ```text
///  Cases addressed are
///      Tests and escapes for negative integer a, b, or c
///      Linear transformation if c - a or c - b negative integer
///      Special case c = a or c = b
///      Linear transformation for  x near +1
///      Transformation for x < -0.5
///      Psi function expansion if x > 0.9 and c - a - b integer
///      Conditionally, a recurrence on c to make c-a-b > 0
///
///      x < -1  AMS 15.3.7 transformation applied (Travis Oliphant)
///         valid for b,a,c,(b-a) != integer and (c-a),(c-b) != negative integer
/// ```
///
/// x >= 1 is rejected (unless special cases are present)
///
/// The parameters a, b, c are considered to be integer
/// valued if they are within 1.0e-14 of the nearest integer
/// (1.0e-13 for IEEE arithmetic).
///
/// ACCURACY:
///
/// ```text
///               Relative error (-1 < x < 1):
/// arithmetic   domain     # trials      peak         rms
///    IEEE      -1,7         20000      1.6e-12     1.6e-14
/// ```
///
/// Several special cases also tested with a, b, c in
/// the range -7 to 7.
///
/// ERROR MESSAGES:
///
/// A "partial loss of precision" message is printed if
/// the internally estimated relative error exceeds 1^-12.
/// A "singularity" message is printed on overflow or
/// in cases not addressed (such as x < -1).
pub(crate) fn hyp2f1(a: f64, b: f64, c: f64, x: f64) -> f64 {
    hyp2f1_impl(a, b, c, x, in_recur_regime(a, b, c, x), true)
}

/// Unlike in xsf, `hyp2f1_recur` is tried first if `max(|a|, |b|, |c|)` is at least this
const HYP2F1_RECUR_MIN: f64 = 10.0;

/// The `x` from which `hyp2f1_recur` tries the 1 - x transformation first, as `hyt2f1` does
pub(crate) const HYP2F1_RECUR_NEAR1: f64 = 0.9;

fn is_large(a: f64, b: f64, c: f64) -> bool {
    [a, b, c].into_iter().any(|v| v.abs() >= HYP2F1_RECUR_MIN)
}

/// Whether `t` is near an integer: for `b - a` or `c - a - b`, the two terms of the 1/x or the
/// 1 - x transformation (DLMF 15.8.2, 15.8.4) then cancel, and `hyp2f1_recur` combines them
fn near_int(t: f64) -> bool {
    (t - t.round()).abs() <= HYP2F1_RECUR_NEAR_INT
}

/// Whether `b - a` (for x <= -9) or `c - a - b` (for x >= 0.9) is near an integer, and `a` and `b`
/// aren't non-positive integers (the polynomials)
fn near_log(a: f64, b: f64, c: f64, x: f64) -> bool {
    let (t, w) = if x < 0.0 {
        (b - a, x / (x - 1.0))
    } else {
        (c - a - b, x)
    };
    w >= HYP2F1_RECUR_NEAR1 && near_int(t) && !(near_neg_int(a) || near_neg_int(b))
}

/// `a`, `b`, `c - a` and `c - b`: if one is a non-positive integer, 2F1 is a polynomial, or one
/// times a power of `1 - x` (DLMF 15.8.1)
fn poly_params(a: f64, b: f64, c: f64) -> [Dd; 4] {
    [Dd::new(a), Dd::new(b), Dd::new(c) - a, Dd::new(c) - b]
}

/// Whether one of `poly_params` is a negative integer up to `HYP2F1_EPS`, but not exactly
fn near_poly(a: f64, b: f64, c: f64) -> bool {
    let near = |t: Dd| t.hi < -0.5 && near_neg_int(t.hi) && !t.is_integer();
    poly_params(a, b, c).into_iter().any(near)
}

/// Whether `hyp2f1_impl` first tries `hyp2f1_recur`, if this also holds for the outermost call
fn in_recur_regime(a: f64, b: f64, c: f64, x: f64) -> bool {
    is_large(a, b, c) || near_log(a, b, c, x) || near_poly(a, b, c)
}

/// `hyp2f1`, where `recur` is whether to try `hyp2f1_recur` in its regime: if the outermost call
/// is `in_recur_regime` (and not in the Pfaff transformation below for x <= -9), and `outer`
/// whether this is the outermost call
fn hyp2f1_impl(a: f64, b: f64, c: f64, x: f64, recur: bool, outer: bool) -> f64 {
    let (mut p, mut q, mut r, mut y);
    let mut neg_int_ca_or_cb = false;

    let ax = x.abs();
    let mut s = 1.0 - x;
    let mut ia = a.round(); /* nearest integer to a */
    let mut ib = b.round();

    if x == 0.0 {
        return 1.0;
    }

    let d = c - a - b;
    let mut id = d.round();

    // Unlike in xsf, also for c = 0 (DLMF 15.2.5 for m = 0, see below)
    if a == 0.0 || b == 0.0 {
        return 1.0;
    }

    // Unlike in xsf, 2F1 -> 1 as c -> inf (DLMF 15.12.2), also at x = 1 (DLMF 15.4.20)
    if c == f64::INFINITY && a.is_finite() && b.is_finite() && x.is_finite() && x <= 1.0 {
        return 1.0;
    }

    // Unlike in xsf, also for b = -1, where the paths below can give NaN. And as `(c - b x) / c`
    // in double-double, instead of rounding results below 1e-15 to 0.
    if (a == -1.0 || b == -1.0) && c != 0.0 {
        let b = if a == -1.0 { b } else { a };
        /* hyp2f1(-1, b; c; x) = 1 - (b/c)*x, c != 0 */
        let rounded = 1.0 - (b / c) * x;
        // As `(c' - b' x') / c'` for the significands `b'` and `c'` of `b` and `c`, and
        // `x' = x 2^(e_b - e_c)`, so that the product (and its rounding error) is subnormal only
        // if it is negligible next to `c'`, and NaN if it overflows
        let ((b, e_b), (c, e_c)) = (frexp_f64(b), frexp_f64(c));
        let x = ldexp_f64(x, e_b - e_c);
        let exact = f64::from((Dd::new(c) - Dd::new(b) * x) / c);
        return if exact.is_nan() { rounded } else { exact };
    }

    let accurate = |&(_, err): &(f64, f64)| err <= HYP2F1_RECUR_TOL;
    // Unlike in xsf, the polynomial of degree m of DLMF 15.2.5 for c = a = -m or c = b = -m, not
    // (1 - x)^(-b) of DLMF 15.2.6 or inf.
    //
    // Unlike in xsf, also at x = 1 for a or b a negative integer above c, the Chu-Vandermonde
    // product of `hyp2f1_recur` (DLMF 15.4.24), as the paths below can give NaN there.
    //
    // (only in the outermost call, as in the transformations c - a can round to c)
    let (c_eq, ends) = (a == c || b == c, |t: f64| is_nonpos_int(t) && t > c);
    if outer && c.is_finite() && is_nonpos_int(c) && (c_eq || x == 1.0 && (ends(a) || ends(b))) {
        // By `hyp2f1_recur` in double-double, else for a = b = c by (1 - x)^m (as both are that
        // then), or 0 if exactly that, else by the paths below
        let y = hyp2f1_recur(a, b, c, x).map(Scaled::to_f64);
        if let Some((y, _)) = y.filter(accurate) {
            return y;
        }
        if a == b && b == c {
            return pow_1mx(x, -c);
        }
        if c_eq && sums_to_zero(if a == c { b } else { a }, -c, x) {
            return 0.0;
        }
    }

    let neg_int_a = near_neg_int(a);
    let neg_int_b = near_neg_int(b);

    // Unlike in xsf, for large parameters first try the backward recurrence in c of
    // `hyp2f1_recur`, since the series and transformations below can lose all precision then. This
    // excludes c close to (but not at) a non-positive integer (unless tiny), and a zero or
    // subnormal negative a, b, c - a or c - b. The transformations below also use it, but only if
    // the outermost call is `in_recur_regime`, so that the other results are unchanged.
    //
    // Unlike in xsf, a tiny negative `c` is not the integer 0 (see `near`), so it isn't excluded:
    // the paths below lose as much precision for it as for any other `c`.
    //
    // Unlike in xsf (scipy/xsf#470), also try it for a tiny negative `a` or `b`. That counts as the
    // integer 0 below, but the paths for polynomials lose it: their power series stops after its
    // tiny first term, or diverges for x < -1. Not for a subnormal one, where the series of
    // `hyp2f1_recur` underflow.
    //
    // Unlike in xsf, also try it in the outermost call for an integer c <= 0 below an integer
    // a or b < 0 (DLMF 15.2.5): it sums that polynomial in double-double, as the paths below lose
    // all digits when its terms cancel. Not in the transformations, where c - a can round to one.
    //
    // Unlike in xsf, also try it for `near_log`, where it combines the two terms of the 1/x or the
    // 1 - x transformation, up to their limit at an integer (DLMF 15.8.9, 15.8.10), while the paths
    // below:
    //
    // - add the two terms (for x < -2, and AMS55 #15.3.6 in `hyt2f1`), which loses about the
    //   inverse of the distance to the integer
    // - at an integer `b - a` (up to `HYP2F1_EPS`), use the Pfaff transformation (DLMF 15.8.1),
    //   where `hyt2f1` needs `1 - w` for its `w = x / (x - 1)` near 1, whose rounding loses about
    //   `log10(-x)` digits of it
    //
    // Unlike in xsf, also try it for `near_poly`, where the paths below take 2F1 for the
    // polynomial at that integer, whose series doesn't terminate: it diverges for x < -1, and
    // can lose all digits elsewhere (it cancels by 1e43 for `F(65.9, -57.6; 7.9; 0.861)`, where
    // c - a is -58 - 5e-15).
    let tiny = |t: f64| t < 0.0 && -t < HYP2F1_EPS && t.is_normal();
    // (0 or a negative subnormal)
    let degenerate = |t: Dd| t.hi <= 0.0 && !t.hi.is_normal();
    let recur_first = |a: f64, b: f64, c: f64| {
        let first = (recur && in_recur_regime(a, b, c, x) || tiny(a) || tiny(b))
            && (tiny(c) || !near_neg_int(c) || outer && c == c.round())
            && !poly_params(a, b, c).into_iter().any(degenerate);
        first.then(|| hyp2f1_recur(a, b, c, x)).flatten()
    };
    // the result of `hyp2f1_recur` and its error bound, if that exceeds its tolerance
    let inexact = recur_first(a, b, c).map(Scaled::to_f64);
    if let Some((y, _)) = inexact.filter(accurate) {
        return y;
    }
    // Unlike in xsf, the result `y` of a transformation below is replaced by that of
    // `hyp2f1_recur` if this shows `y` to be less accurate, or `y` is not finite
    let or_recur = |y: f64| {
        let worse = |&(alt, bound): &(f64, f64)| {
            bound < HYP2F1_LOST && (!y.is_finite() || less_accurate(y, alt, bound))
        };
        inexact.filter(worse).map_or(y, |(alt, _)| alt)
    };

    // Unlike in xsf, the divergent 2F1 at x = 1 has the sign of its limit from below, that of
    // `Gamma(c) / (Gamma(a) Gamma(b))` (DLMF 15.4.21, 15.4.23). Also for c - a - b <= -1, where
    // the Euler transformation below is `inf * 0` if its gamma functions overflow.
    if d <= 0.0 && x == 1.0 && !(neg_int_a || neg_int_b) {
        // (+inf at a pole of `Gamma(c)`, as for x < 1)
        let pole = c <= 0.0 && near(c, c.round());
        let sign = gammasgn(a) * gammasgn(b) * gammasgn(c);
        return f64::INFINITY * if pole { 1.0 } else { sign }; // goto hypdiv
    }

    // Unlike in xsf, for c `near_pole`, sum the power series (or for x < 0 that of the Pfaff
    // transformations, as it can stop early), or NaN if inaccurate: at c = -n + e its terms can
    // jump by about `1 / e` at the term n + 1, which the paths below miss. Not if a or b is an
    // integer that terminates it before that term, nor if x / (x - 1) rounds to 1 or exceeds it
    // (x >= 1).
    let ic = c.round(); /* nearest integer to c */
    // (whether a or b is an integer `is_int` that terminates the series before the term 1 - c)
    let terminates = |is_int: fn(f64) -> bool| is_int(a) && ia > ic || is_int(b) && ib > ic;
    // Unlike in xsf, +inf at a pole of c (as below) for any finite x, where xsf can return -inf,
    // NaN or a finite value. Not for c = a or c = b (see above).
    let pole = c.is_finite() && is_nonpos_int(c) && !terminates(is_nonpos_int);
    if pole && a != c && b != c && x.is_finite() {
        return f64::INFINITY; // goto hypdiv
    }
    if near_pole(c) && !terminates(is_nonpos_int) && x / (x - 1.0) < 1.0 {
        let (y, err) = if x > 0.0 {
            hys2f1_min(a, b, c, x)
        } else {
            hys2f1_pfaff(a, b, c, x)
        };
        return if err < HYP2F1_ETHRESH { y } else { f64::NAN };
    }

    // Unlike in xsf, for x > 1 and c = a or c = b, the power `(1 - x)^-b` or `(1 - x)^-a` of
    // DLMF 15.4.6 if its exponent is an integer, else NaN (see below), as the Euler transformation
    // rounds `c - a - b`, and skips it if `c - a` or `c - b` rounds to `c`. Not for an integer
    // c <= 0, where c = a = -m is the polynomial of DLMF 15.2.5 (see above), nor for an exponent
    // `near` an integer, which the paths below take as one, nor in the transformations, whose
    // `c - a` or `c - b` can round to `c`.
    if outer && x > 1.0 && (a == c || b == c) && !is_nonpos_int(c) {
        let e = if a == c { b } else { a };
        if e == e.round() {
            return s.powf(-e);
        } else if !near(e, e.round()) {
            return f64::NAN;
        }
    }

    // Unlike in xsf (scipy/xsf#470), not if `c - a` or `c - b` rounds to `c`, where the Euler
    // transformation (DLMF 15.8.1) loses `a` or `b` and returns about 1, as its power factor
    // cancels `F(c, c - b; c; x) = (1 - x)^(b - c)` (DLMF 15.4.6). Only if c > 0, as for c < 0
    // the other paths can be less accurate (in the recurrence in `a`).
    let lost = c > 0.0 && (c - a == c || c - b == c);
    if d <= -1.0 && !((d - id).abs() > HYP2F1_EPS && s < 0.0) && !(neg_int_a || neg_int_b || lost) {
        // Unlike in xsf (scipy/xsf#464), if c < 0 is more than `HYP2F1_EPS` from an integer,
        // where this transformation tends to overflow or cancel, and -1 < x < 0, first try the
        // series of the Pfaff transformations (as for x < -2 below)
        let c_pole = c < 0.0 && (c - c.round()).abs() > HYP2F1_EPS;
        if c_pole && x < 0.0 && x > -1.0 {
            let (y, err) = hys2f1_pfaff(a, b, c, x);
            if err < HYP2F1_ETHRESH {
                return y;
            }
        }
        // Unlike in xsf (scipy/xsf#464), if 0 < x < 1, first try the power series (see
        // `hys2f1_min`) for such a c < 0, or for c > 0 if the series of this transformation is
        // estimated to lose precision
        let lossy = || c > 0.0 && hys2f1(c - a, c - b, c, x).1 >= HYP2F1_ETHRESH;
        if x > 0.0 && x < 1.0 && (c_pole || lossy()) {
            let (y, err) = hys2f1_min(a, b, c, x);
            // for c > 0, also keep it if it didn't fail (estimate 1) and its estimate is smaller
            // than that of `hyt2f1`, which the call below uses unless c - (c - a) or c - (c - b)
            // rounds to a non-positive integer
            let special = [c - (c - a), c - (c - b)].into_iter().any(near_neg_int);
            let euler_worse =
                || c > 0.0 && err < 1.0 && !special && err < hyt2f1(c - a, c - b, c, x).1;
            if err < HYP2F1_ETHRESH || euler_worse() {
                return y;
            }
        }
        let p = s.powf(d);
        let y = hyp2f1_impl(c - a, c - b, c, x, recur, false);
        // Unlike in xsf, if `y` overflows for a tiny `c` and is by `hyp2f1_recur`, multiply its
        // scaled value by `s^d`, as their product can still be finite.
        //
        // Only if that is at least `1 / HYP2F1_ETHRESH`: the rounding of `c - a` and `c - b` (e.g.
        // to a negative integer, where `y` is a polynomial) loses a part of the result that
        // isn't proportional to `1 / c`.
        if y.is_infinite() && c.abs() < HYP2F1_EPS {
            // (in double-double, as a large `d` amplifies the rounding of `d` and `s`)
            let ln_p = (Dd::new(c) - a - b) * ln_dd(Dd::new(1.0) - x);
            let scaled = recur_first(c - a, c - b, c).map(|y| y.mul_exp(ln_p).to_f64());
            let large = |&(y, _): &(f64, f64)| y.abs() * HYP2F1_ETHRESH >= 1.0;
            if let Some((y, _)) = scaled.filter(accurate).filter(large) {
                return y;
            }
        }
        let mut euler = p * y;
        if !p.is_normal() {
            // Unlike in xsf (scipy/xsf#471), split `s^d` in two if s > 0 and `y` is normal, as
            // their product can still be normal (e.g. for x << -1, where `y` is huge)
            if s > 0.0 && y.is_normal() {
                let h = s.powf(d / 2.0);
                euler = h * y * h;
            } else if x < 0.0 && x / (x - 1.0) < 1.0 {
                // Unlike in xsf (scipy/xsf#471), else for x < 0 try the Pfaff series (as for
                // x < -2 below) and keep a normal result (e.g. about 1 for a tiny `b` and c < 0);
                // not if their argument x / (x - 1) rounds to 1, which their error estimate misses
                let (y1, e1) = hys2f1_pfaff(a, b, c, x);
                if e1 < HYP2F1_ETHRESH && y1.is_normal() {
                    return y1;
                }
            }
        }
        return or_recur(euler);
    }
    if ax < 1.0 || x == -1.0 {
        /* 2F1(a,b;b;x) = (1-x)**(-a) */
        if near(b, c) {
            /* b = c */
            if neg_int_b {
                y = hyp2f1_neg_c_equal_bc(a, b, x);
            } else {
                y = s.powf(-a); /* s to the -a power */
            }
            return y; // goto hypdon
        }
        if near(a, c) {
            /* a = c */
            y = s.powf(-b); /* s to the -b power */
            return y; // goto hypdon
        }
    }

    'hypf: {
        'hypok: {
            if c <= 0.0 && near(c, ic) {
                /* c is a negative integer */
                /* check if termination before explosion */
                if terminates(near_neg_int) {
                    break 'hypok;
                }
                return f64::INFINITY; // goto hypdiv
            }

            if neg_int_a || neg_int_b {
                /* function is a polynomial */
                break 'hypok;
            }

            let mut t1 = (b - a).abs();
            if x < -2.0 && (t1 - t1.round()).abs() > HYP2F1_EPS {
                // Unlike in xsf (scipy/xsf#465), also sum the series of the Pfaff transformations
                // as for x < -1, since the inner functions of the 1/x transformation can be
                // inaccurate
                let (y1, e1) = hys2f1_pfaff(a, b, c, x);
                let y1_ok = e1 < HYP2F1_ETHRESH;
                // The rounding of `t` relative to its distance to the nearest pole of `Gamma(t)`.
                // Unlike in xsf, at a pole of `Gamma(c - a)` or `Gamma(c - b)`, keep the series if
                // it is accurate.
                let pole = |t: f64| MACHEP * (t / (t - t.round().min(0.0))).abs();
                let (pa, pb) = (pole(c - a), pole(c - b));
                if y1_ok && !(pa + pb).is_finite() {
                    return y1;
                }
                /* This transform has a pole for b-a integer, and
                 * may produce large cancellation errors for |1/x| close 1
                 */
                // Unlike in xsf, round `b - a` only once, since the poles of `Gamma(b - a)` and
                // `Gamma(a - b)` only cancel those of `p` and `q` at the same `b - a`
                let e = b - a;
                // Unlike in xsf, near an integer `b - a`, where `s p` and `y q` below cancel, first
                // try `hyp2f1_recur` (as for `near_log`), unless it was tried above or the series
                // is accurate
                if inexact.is_none() && near_int(e) && !y1_ok {
                    let log = hyp2f1_recur(a, b, c, x).map(Scaled::to_f64);
                    if let Some((y, _)) = log.filter(accurate) {
                        return y;
                    }
                }
                // Unlike in xsf, if the series is accurate, the inner functions are by `hyt2f1`,
                // for its error estimate
                let f = |a, b, c| {
                    if y1_ok {
                        let (y, _, err) = hyt2f1(a, b, c, 1.0 / x);
                        (y, err)
                    } else {
                        (hyp2f1_impl(a, b, c, 1.0 / x, recur, false), f64::NAN)
                    }
                };
                let (ep, eq);
                (p, ep) = f(a, 1.0 - c + a, 1.0 - e);
                (q, eq) = f(b, 1.0 - c + b, 1.0 + e);
                p *= (-x).powf(-a);
                q *= (-x).powf(-b);
                t1 = gamma(c);
                s = t1 * gamma(e) * (rgamma(b) * rgamma(c - a));
                y = t1 * gamma(-e) * (rgamma(a) * rgamma(c - b));
                // Unlike in xsf, keep the series if its error estimate is smaller than the
                // cancellation of `s p + y q`, or if it is accurate, unless a heuristic error
                // estimate of `s p + y q` is smaller and at most `HYP2F1_RECIP_MAX_ERR`, and `p`,
                // `q` and the sum are normal (the estimate misses the rounding of a subnormal).
                //
                // That adds, for `p` and `q` by their share of the sum, their error estimates and
                // the rounding of `b - a` and of `c - a` or `c - b`, each as in `pole`.
                let (u, v) = (s * p, y * q);
                let canc = MACHEP * u.abs().max(v.abs()) / (u + v).abs();
                let weigh = |t: f64, err: f64| (t / (u + v)).abs() * (err + pole(-e.abs()));
                let full = canc + weigh(u, ep + pa) + weigh(v, eq + pb);
                let normal = [p, q, u + v].into_iter().all(f64::is_normal);
                if e1 < canc || y1_ok && !(normal && full < e1.min(HYP2F1_RECIP_MAX_ERR)) {
                    return y1;
                }
                return or_recur(u + v);
            } else if x < -1.0 || lost && near(x, -1.0) {
                // Unlike in xsf (scipy/xsf#470), also near x = -1 if `a` or `b` is lost (see
                // above), as the power series diverges there if c - a - b <= -1.
                //
                // (not `hyp2f1_recur` at x / (x - 1) >= 0.9, i.e. x <= -9, since its 1 - x
                // transformation and combinations there are those that it tried for x)
                let w = x / (x - 1.0);
                let recur = recur && w < HYP2F1_RECUR_NEAR1;
                if a.abs() < b.abs() {
                    return or_recur(s.powf(-a) * hyp2f1_impl(a, c - b, c, w, recur, false));
                }
                return or_recur(s.powf(-b) * hyp2f1_impl(b, c - a, c, w, recur, false));
            }

            // Unlike in xsf, NaN instead of inf for x > 1 (x < -1 returned above), where 2F1 is
            // generally complex, as it jumps across its branch cut (DLMF 15.2.3). (Most of its
            // real cases, the polynomials and their Euler transformations, are handled above.)
            if ax > 1.0 {
                /* series diverges  */
                return f64::NAN;
            }

            p = c - a;
            ia = p.round(); /* nearest integer to c-a */
            if (ia <= 0.0) && near(p, ia) {
                /* negative int c - a */
                neg_int_ca_or_cb = true;
            }

            r = c - b;
            ib = r.round(); /* nearest integer to c-b */
            if (ib <= 0.0) && near(r, ib) {
                /* negative int c - b */
                neg_int_ca_or_cb = true;
            }

            id = d.round(); /* nearest integer to d */

            /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
             * for reporting a bug here.  */
            if (ax - 1.0).abs() < HYP2F1_EPS {
                /* |x| == 1.0   */
                // Unlike in xsf, only at x = 1 itself: up to `HYP2F1_EPS` below it, 2F1 is finite
                // for d <= 0, and differs from its value at 1 by about `(1 - x)^min(d, 1)` (times
                // `ln(1 - x)` if d = 1) for d > 0, while the transformations of `hyt2f1` converge.
                if x == 1.0 {
                    if neg_int_ca_or_cb {
                        if d >= 0.0 {
                            break 'hypf;
                        }
                        return f64::INFINITY; // goto hypdiv
                    }
                    if d <= 0.0 {
                        return f64::INFINITY; // goto hypdiv
                    }
                    y = gamma(c) * gamma(d) * (rgamma(p) * rgamma(r));
                    return y; // goto hypdon
                }
                if x < 0.0 && d <= -1.0 {
                    return f64::INFINITY; // goto hypdiv
                }
            }

            /* Conditionally make d > 0 by recurrence on c
             * AMS55 #15.2.27
             */
            if d < 0.0 {
                /* Try the power series first */
                let (err, err_full);
                (y, err, err_full) = hyt2f1(a, b, c, x);
                // Unlike in xsf, a negative integer c - a or c - b uses the polynomial
                // transformation (as for d >= 0) if its error estimate is smaller than that of `y`
                // including the recurrence in `a`, and it shows `y` to be less accurate: `err`
                // misses that error, and the recurrence in `c` below can lose all digits (e.g.
                // across c = 0).
                if neg_int_ca_or_cb && !(err_full < HYP2F1_ETHRESH) {
                    let (f, _, f_full) = hys2f1(c - a, c - b, c, x);
                    let z = s.powf(d) * f;
                    if !(err_full <= f_full) && less_accurate(y, z, f_full) {
                        return z;
                    }
                }
                if err < HYP2F1_ETHRESH {
                    return y; // goto hypdon
                }
                /* Apply the recurrence if power series fails */
                let aid = (2.0 - id) as i32;
                let mut e = c + f64::from(aid);
                // Unlike in xsf (scipy/xsf#466, scipy/xsf#470), use Euler (DLMF 15.8.1), or Pfaff
                // for x < 0.5 if its error estimate is smaller: for |e| >= 2^53, where the unit
                // shifts can round and recurse forever, or if `a` or `b` is lost in `c - a` or
                // `c - b` (see above), as the recurrence can then amplify the rounding of its
                // starting values
                if e.abs() * MACHEP >= 1.0 || lost {
                    let (z, err, _) = hyt2f1(c - a, c - b, c, x);
                    if x < 0.5 {
                        let (y, pfaff_err) = hys2f1_pfaff(a, b, c, x);
                        if pfaff_err < err {
                            return y;
                        }
                    }
                    return s.powf(d) * z;
                }
                let mut d2 = hyp2f1_impl(a, b, e, x, recur, false);
                let mut d1 = hyp2f1_impl(a, b, e + 1.0, x, recur, false);
                q = a + b + 1.0;
                for _ in 0..aid {
                    r = e - 1.0;
                    y = (e * (r - (2.0 * e - q) * x) * d2 + (e - a) * (e - b) * x * d1)
                        / (e * r * s);
                    e = r;
                    d1 = d2;
                    d2 = y;
                }
                return y; // goto hypdon
            }

            if neg_int_ca_or_cb {
                break 'hypf; /* negative integer c-a or c-b */
            }
        }
        // hypok:
        let (y, err, err_full) = hyt2f1(a, b, c, x);

        // Unlike in xsf, the result of `hyp2f1_recur` instead if its error bound is below the
        // estimate of `y`, or shows `y` to be less accurate.
        if let Some((alt, bound)) = inexact.filter(|&(_, e)| e < HYP2F1_LOST) {
            if !(err_full <= bound) || less_accurate(y, alt, bound) {
                return alt;
            }
        }

        // hypdon:
        // Unlike in xsf (scipy/xsf#371), NaN instead of a finite `y` if its error estimates
        // (without and with the recurrences in `a`) say that no digit is correct.
        //
        // Not 1 itself: a series that fails to converge is estimated as 1, while its value can
        // be accurate.
        let no_digit = |e| (HYP2F1_LOST..1.0).contains(&e);
        // The estimate is inf if it overflows or `y` is 0. Not for a polynomial: near a root its
        // estimate also exceeds 1, while `y` is accurate in absolute terms.
        let wrong = err_full > 1.0 && err_full.is_finite() && !(neg_int_a || neg_int_b);
        if y.is_finite() && (no_digit(err) && no_digit(err_full) || wrong) {
            return f64::NAN;
        }
        return y;
    }

    /* The transformation for c-a or c-b negative integer
     * AMS55 #15.3.3
     */
    // hypf:
    (y, _, _) = hys2f1(c - a, c - b, c, x);
    s.powf(d) * y // goto hypdon
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `times_gamma_ratio` with `q = s^k`
    fn ratio(y: (f64, f64), s: f64, k: f64, num: [f64; 2], den: [f64; 2]) -> (f64, f64) {
        let [n0, n1] = num.map(gamma);
        let f = n0 * n1 * (rgamma(den[0]) * rgamma(den[1]));
        times_gamma_ratio(y, f, (s.powf(k), k * s.ln()), num, den)
    }

    #[test]
    fn test_times_gamma_ratio_out_of_range() {
        // Arb; `s^k` or the gamma factors underflow, but their product with `y` doesn't
        for (s, k, den, desired) in [
            (0.5, 2000.0, [1.0, 1.0], 8.709_809_816_217_217e-303),
            (1.0, 0.0, [150.0, 150.0], 6.892_800_301_195_469e-222),
        ] {
            let (y, err) = ratio((1e300, 0.0), s, k, [1.0, 1.0], den);
            let rel_err = err * MACHEP / desired;
            assert!((y / desired - 1.0).abs() <= rel_err && rel_err < 1e-12, "{y} {err}");
        }
        // normal factors, but `y f` underflows: 2^-1074 / Gamma(3) 2^1000 = 2^-75
        let tiny = f64::from_bits(1);
        let (y, err) = ratio((tiny, tiny), 0.5, -1000.0, [1.0, 1.0], [3.0, 1.0]);
        let rel_err = err * MACHEP / 2f64.powi(-75);
        assert!((y / 2f64.powi(-75) - 1.0).abs() <= rel_err && rel_err < 1e-12, "{y} {err}");
        // subnormal products (0.1 2^-1022 and 7.9e-324), and a pole of the denominator
        let (y, err) = ratio((0.1, 0.1), 0.5, 1022.0, [1.0, 1.0], [1.0, 1.0]);
        assert!((y - 2.225_073_858_507_201_4e-309).abs() <= err * MACHEP, "{y} {err}");
        let (y, err) = ratio((0.1, 0.0), 0.5, 1070.0, [1.0, 1.0], [1.0, 1.0]);
        assert!((y - 7.905e-324).abs() <= err * MACHEP, "{y} {err}");
        assert_eq!(ratio((1e300, 1.0), 0.5, 2000.0, [1.0, 1.0], [-3.0, 1.0]), (0.0, 0.0));
    }
}
