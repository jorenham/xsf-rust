//! Translated into pure Rust from `xsf/cephes/hyp2f1.h` (scipy/xsf@f2a68ec).
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
//! - The threshold of the psi function expansion is corrected in the comments.
//!
//! Our `psi` is more accurate than xsf's just above negative integers, so the psi function
//! expansion (for `x > 0.9` and integer `c - a - b`) can differ from xsf when `a`, `b`, `c - a` or
//! `c - b` is close to a negative integer.

use crate::utils::MACHEP;
use crate::xsf::cephes::{gamma, lgam_ratio, lgam_sgn, psi, rgamma};

pub(crate) const HYP2F1_EPS: f64 = 1.0e-13;

pub(crate) const HYP2F1_ETHRESH: f64 = 1.0e-12;
/// The estimated relative error from which `hyp2f1` returns NaN: if a sum is only rounding error,
/// then its estimate is about 1/2 to 1.
const HYP2F1_LOST: f64 = 0.5;
const HYP2F1_MAXITER: u64 = 10000;
/// 2^512 (about `sqrt(f64::MAX)`), to rescale tiny terms exactly, far from overflow and underflow
const HYP2F1_SCALE: f64 = 1.340_780_792_994_259_7e154;

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

/// The series of the Pfaff transformations (DLMF 15.8.1) in `x / (x - 1)`, by `hys2f1_min`, with
/// the smaller error estimate
fn hys2f1_pfaff(a: f64, b: f64, c: f64, x: f64) -> (f64, f64) {
    let w = x / (x - 1.0);
    let (ya, ea) = hys2f1_min(a, c - b, c, w);
    let (yb, eb) = hys2f1_min(b, c - a, c, w);
    // if one of these overflows, its error estimate is NaN, so use the other
    let (y, err, p) = if ea <= eb || eb.is_nan() {
        (ya, ea, a)
    } else {
        (yb, eb, b)
    };
    ((1.0 - x).powf(-p) * y, err)
}

/// Defining power series expansion of Gauss hypergeometric function.
/// The `loss` (the second return value) estimates loss of significance. Unlike in xsf, the third
/// is an estimate of the relative error that also includes the error of the recurrence in `a`.
fn hys2f1(mut a: f64, mut b: f64, c: f64, x: f64) -> (f64, f64, f64) {
    let mut f;
    let mut intflag = false;

    // Unlike in xsf (scipy/xsf#460), also for `b = -a` a negative integer, so that the recurrence
    // is in `b`
    let neg_int_tie = b == -a && near_neg_int(b) && a <= HYP2F1_MAXITER as f64;
    if b.abs() > a.abs() || neg_int_tie {
        /* Ensure that |a| > |b| ... */
        f = b;
        b = a;
        a = f;
    }

    let ib = b.round() as i32;

    if (b - f64::from(ib)).abs() < HYP2F1_EPS && ib <= 0 && b.abs() < a.abs() {
        /* .. except when `b` is a smaller negative integer */
        f = b;
        b = a;
        a = f;
        intflag = true;
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
    // after it:
    //
    // - so this may take up to twice `HYP2F1_MAXITER` terms (before and after k = -c)
    // - terms are scaled by `HYP2F1_SCALE^n_scaled` so that they don't underflow, and only
    //   summed (and counted in `loss`) if `n_scaled == 0` (else they are negligible)
    // - it fails if their ratio underflows, e.g. for tiny a and b
    let pole_ahead = x > 0.0 && h < 0.0 && h.fract() != 0.0;
    let maxiter = HYP2F1_MAXITER * if pole_ahead { 2 } else { 1 };
    let (mut n_scaled, mut i_scaled) = (0, 0);
    loop {
        // Unlike in xsf, only if c = 0, since the series is finite for a tiny c != 0
        if h == 0.0 {
            return (f64::INFINITY, 1.0, 1.0);
        }
        let m = k + 1.0;
        let r = (f + k) * (g + k) * x / ((h + k) * m);
        u *= r;
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
    let loss = (MACHEP * umax) / s.abs() + (MACHEP * (i - i_scaled) as f64);

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

/// Apply transformations for |x| near 1 then call the power series. Unlike in xsf, also returns
/// an error estimate that includes the recurrences in `a`.
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

    let (mut err, mut err_full) = (0.0, 0.0);
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
            // can be inexact or fail. Unlike in xsf (scipy/xsf#371), if its `loss` is small but
            // its error estimate incl. the recurrence in `a` isn't, only replace it by AMS55
            // #15.3.6 if that shows it to be less accurate.
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

            p = (a + d1) * (b + d1) * s * rgamma(e + 2.0); /* Poch for t=1 */
            t = 1.0;
            loop {
                r = psi(1.0 + t) + psi(1.0 + t + e) - psi(a + t + d1) - psi(b + t + d1) - ax;
                q = p * r;
                y += q;
                p *= s * (a + t + d1) / (t + 1.0);
                p *= (b + t + d1) / (t + 1.0 + e);
                t += 1.0;
                if t > HYP2F1_MAXITER as f64 {
                    /* should never happen */
                    err = 1.0;
                    break 'done f64::NAN;
                }
                if !(y == 0.0 || (q / y).abs() > HYP2F1_EPS) {
                    break;
                }
            }

            if id == 0.0 {
                y *= gamma(c) / (gamma(a) * gamma(b));
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
                }
            }
            // nosum:
            p = gamma(c);
            y1 *= gamma(e) * p * (rgamma(a + d1) * rgamma(b + d1));

            y *= p * (rgamma(a + d2) * rgamma(b + d2));
            if (aid & 1) != 0 {
                y = -y;
            }

            q = s.powf(id); /* s to the id power */
            if id > 0.0 {
                y *= q;
            } else {
                y1 *= q;
            }

            y += y1;
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
    let (mut p, mut q, mut r, mut y);
    let mut neg_int_a = false;
    let mut neg_int_b = false;
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

    if (a == 0.0 || b == 0.0) && c != 0.0 {
        return 1.0;
    }

    // Unlike in xsf, 2F1 -> 1 as c -> inf (DLMF 15.12.2), also at x = 1 (DLMF 15.4.20)
    if c == f64::INFINITY && a.is_finite() && b.is_finite() && x.is_finite() && x <= 1.0 {
        return 1.0;
    }

    if a == -1.0 && c != 0.0 {
        let result = 1.0 - (b / c) * x; /* hyp2f1(-1, b; c; x) = 1 - (b/c)*x, c != 0 */
        if result.abs() < 1e-15 {
            return 0.0;
        }
        return result;
    }

    if a <= 0.0 && (a - ia).abs() < HYP2F1_EPS {
        /* a is a negative integer */
        neg_int_a = true;
    }

    if b <= 0.0 && (b - ib).abs() < HYP2F1_EPS {
        /* b is a negative integer */
        neg_int_b = true;
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
        let y = hyp2f1(c - a, c - b, c, x);
        if !p.is_normal() {
            // Unlike in xsf (scipy/xsf#471), split `s^d` in two if s > 0 and `y` is normal, as
            // their product can still be normal (e.g. for x << -1, where `y` is huge)
            if s > 0.0 && y.is_normal() {
                let h = s.powf(d / 2.0);
                return h * y * h;
            }
            // Unlike in xsf (scipy/xsf#471), else for x < 0 try the Pfaff series (as for x < -2
            // below) and keep a normal result (e.g. about 1 for a tiny `b` and c < 0); not if
            // their argument x / (x - 1) rounds to 1, as their error estimate misses that
            if x < 0.0 && x / (x - 1.0) < 1.0 {
                let (y1, e1) = hys2f1_pfaff(a, b, c, x);
                if e1 < HYP2F1_ETHRESH && y1.is_normal() {
                    return y1;
                }
            }
        }
        return p * y;
    }
    if d <= 0.0 && x == 1.0 && !(neg_int_a || neg_int_b) {
        return f64::INFINITY; // goto hypdiv
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
            if c <= 0.0 {
                let ic = c.round(); /* nearest integer to c */
                if near(c, ic) {
                    /* c is a negative integer */
                    /* check if termination before explosion */
                    if neg_int_a && (ia > ic) {
                        break 'hypok;
                    }
                    if neg_int_b && (ib > ic) {
                        break 'hypok;
                    }
                    return f64::INFINITY; // goto hypdiv
                }
            }

            if neg_int_a || neg_int_b {
                /* function is a polynomial */
                break 'hypok;
            }

            let mut t1 = (b - a).abs();
            if x < -2.0 && (t1 - t1.round()).abs() > HYP2F1_EPS {
                // Unlike in xsf (scipy/xsf#465), first try the series of the Pfaff transformations
                // as for x < -1, since the inner functions of the 1/x transformation can be
                // inaccurate
                let (y1, e1) = hys2f1_pfaff(a, b, c, x);
                if e1 < HYP2F1_ETHRESH {
                    return y1;
                }
                /* This transform has a pole for b-a integer, and
                 * may produce large cancellation errors for |1/x| close 1
                 */
                // Unlike in xsf, round `b - a` only once, since the poles of `Gamma(b - a)` and
                // `Gamma(a - b)` only cancel those of `p` and `q` at the same `b - a`
                let e = b - a;
                p = hyp2f1(a, 1.0 - c + a, 1.0 - e, 1.0 / x);
                q = hyp2f1(b, 1.0 - c + b, 1.0 + e, 1.0 / x);
                p *= (-x).powf(-a);
                q *= (-x).powf(-b);
                t1 = gamma(c);
                s = t1 * gamma(e) * (rgamma(b) * rgamma(c - a));
                y = t1 * gamma(-e) * (rgamma(a) * rgamma(c - b));
                // Unlike in xsf, keep the series if its error estimate is smaller than the
                // cancellation of `s p + y q`
                let (u, v) = (s * p, y * q);
                if e1 < MACHEP * u.abs().max(v.abs()) / (u + v).abs() {
                    return y1;
                }
                return u + v;
            } else if x < -1.0 || lost && near(x, -1.0) {
                // Unlike in xsf (scipy/xsf#470), also near x = -1 if `a` or `b` is lost (see
                // above), as the power series diverges there if c - a - b <= -1.
                if a.abs() < b.abs() {
                    return s.powf(-a) * hyp2f1(a, c - b, c, x / (x - 1.0));
                }
                return s.powf(-b) * hyp2f1(b, c - a, c, x / (x - 1.0));
            }

            if ax > 1.0 {
                /* series diverges  */
                return f64::INFINITY; // goto hypdiv
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
                if x > 0.0 {
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
                if d <= -1.0 {
                    return f64::INFINITY; // goto hypdiv
                }
            }

            /* Conditionally make d > 0 by recurrence on c
             * AMS55 #15.2.27
             */
            if d < 0.0 {
                /* Try the power series first */
                let err;
                (y, err, _) = hyt2f1(a, b, c, x);
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
                let mut d2 = hyp2f1(a, b, e, x);
                let mut d1 = hyp2f1(a, b, e + 1.0, x);
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

        // hypdon:
        // Unlike in xsf (scipy/xsf#371), NaN instead of a finite `y` if its error estimates say
        // that no digit is correct:
        //
        // - if both, without and with the recurrences in `a`, are in [`HYP2F1_LOST`, 1); not 1
        //   itself, as a series that fails (doesn't converge) is estimated as 1, while its value
        //   can be accurate
        // - or if the latter is above 1 and finite (it is inf if it overflows or `y` is 0); not for
        //   a polynomial, as near a root its estimate is also above 1, while `y` is accurate in
        //   absolute terms
        let no_digit = |e| (HYP2F1_LOST..1.0).contains(&e);
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
