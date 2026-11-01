//! Translated into pure Rust from `xsf/cephes/hyp2f1.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
//!
//! Unlike in xsf, the `loss` estimates are returned together with the results, instead of written
//! to output references. And because `set_error` is a no-op in our build of xsf, its calls are
//! omitted, so that `goto hypdon` and `goto hypdiv` become plain returns, and the assignments to
//! `err` and `q` in `hyp2f1` that are never read are dropped. The threshold of the psi function
//! expansion is corrected in the comments.
//!
//! Note that our `psi` is more accurate than the one in xsf just above negative integers. So the
//! results of the psi function expansion (for `x > 0.9` and integer `c - a - b`) can differ from
//! those in xsf when `a`, `b`, `c - a` or `c - b` is close to a negative integer.

use crate::utils::MACHEP;
use crate::xsf::cephes::{gamma, lgam_sgn, psi, rgamma};

const HYP2F1_EPS: f64 = 1.0e-13;

const HYP2F1_ETHRESH: f64 = 1.0e-12;
const HYP2F1_MAXITER: u64 = 10000;

/// Defining power series expansion of Gauss hypergeometric function.
/// The `loss` (the second return value) estimates loss of significance.
#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]
fn hys2f1(mut a: f64, mut b: f64, c: f64, x: f64) -> (f64, f64) {
    let mut f;
    let mut intflag = false;

    if b.abs() > a.abs() {
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

    let mut i: u64 = 0;
    let mut umax = 0.0;
    f = a;
    let g = b;
    let h = c;
    let mut s = 1.0;
    let mut u = 1.0;
    let mut k = 0.0;
    loop {
        if h.abs() < HYP2F1_EPS {
            return (f64::INFINITY, 1.0);
        }
        let m = k + 1.0;
        u *= (f + k) * (g + k) * x / ((h + k) * m);
        s += u;
        k = u.abs(); /* remember largest term summed */
        if k > umax {
            umax = k;
        }
        k = m;
        i += 1;
        if i > HYP2F1_MAXITER {
            /* should never happen */
            return (s, 1.0);
        }
        if !(s == 0.0 || (u / s).abs() > MACHEP) {
            break;
        }
    }

    /* return estimated relative error */
    let loss = (MACHEP * umax) / s.abs() + (MACHEP * i as f64);

    (s, loss)
}

/// Apply transformations for |x| near 1 then call the power series
#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]
fn hyt2f1(a: f64, b: f64, c: f64, x: f64) -> (f64, f64) {
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

    let mut err = 0.0;
    let s = 1.0 - x;
    let y = 'done: {
        if x < -0.5 && !(neg_int_a || neg_int_b) {
            let y;
            if b > a {
                (y, err) = hys2f1(a, c - b, c, -x / s);
                break 'done s.powf(-a) * y;
            }
            (y, err) = hys2f1(c - a, b, c, -x / s);
            break 'done s.powf(-b) * y;
        }

        d = c - a - b;
        let id = d.round(); /* nearest integer to d */

        if x > 0.9 && !(neg_int_a || neg_int_b) {
            if (d - id).abs() > MACHEP {
                let (mut lg, mut sgngam);
                let err1;
                let mut y;

                /* test for integer c-a-b */
                /* Try the power series first */
                (y, err) = hys2f1(a, b, c, x);
                if err < HYP2F1_ETHRESH {
                    break 'done y;
                }
                /* If power series fails, then apply AMS55 #15.3.6 */
                (q, err) = hys2f1(a, b, 1.0 - d, s);
                sign = 1;
                (w, sgngam) = lgam_sgn(d);
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(c - a);
                w -= lg;
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(c - b);
                w -= lg;
                sign *= sgngam;
                q *= f64::from(sign) * w.exp();
                (r, err1) = hys2f1(c - a, c - b, d + 1.0, s);
                r *= s.powf(d);
                sign = 1;
                (w, sgngam) = lgam_sgn(-d);
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(a);
                w -= lg;
                sign *= sgngam;
                (lg, sgngam) = lgam_sgn(b);
                w -= lg;
                sign *= sgngam;
                r *= f64::from(sign) * w.exp();
                y = q + r;

                q = q.abs(); /* estimate cancellation error */
                r = r.abs();
                if q > r {
                    r = q;
                }
                err += err1 + (MACHEP * r) / y;

                y *= gamma(c);
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
                    return (f64::NAN, 1.0);
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
        let y;
        (y, err) = hys2f1(a, b, c, x);
        y
    };
    // done:
    (y, err)
}

/*
  15.4.2 Abramowitz & Stegun.
*/
#[allow(clippy::neg_cmp_op_on_partial_ord)]
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
#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
fn hyp2f1ra(a: f64, b: f64, c: f64, x: f64) -> (f64, f64) {
    let (mut f2, mut f1, mut f0);
    let mut err;

    /* Don't cross c or zero */
    let da = if (c < 0.0 && a <= c) || (c >= 0.0 && a >= c) {
        (a - c).round()
    } else {
        a.round()
    };
    let mut t = a - da;

    let mut loss = 0.0;

    debug_assert!(da != 0.0);

    if da.abs() > HYP2F1_MAXITER as f64 {
        /* Too expensive to compute this value, so give up */
        return (f64::NAN, 1.0);
    }

    if da < 0.0 {
        /* Recurse down */
        (f1, err) = hys2f1(t, b, c, x);
        loss += err;
        (f0, err) = hys2f1(t - 1.0, b, c, x);
        loss += err;
        t -= 1.0;
        let mut n = 1.0;
        while n < -da {
            f2 = f1;
            f1 = f0;
            f0 = -(2.0 * t - c - t * x + b * x) / (c - t) * f1 - t * (x - 1.0) / (c - t) * f2;
            t -= 1.0;
            n += 1.0;
        }
    } else {
        /* Recurse up */
        (f1, err) = hys2f1(t, b, c, x);
        loss += err;
        (f0, err) = hys2f1(t + 1.0, b, c, x);
        loss += err;
        t += 1.0;
        let mut n = 1.0;
        while n < da {
            f2 = f1;
            f1 = f0;
            f0 = -((2.0 * t - c - t * x + b * x) * f1 + (c - t) * f2) / (t * (x - 1.0));
            t += 1.0;
            n += 1.0;
        }
    }

    (f0, loss)
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
#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::cast_possible_truncation,
    clippy::too_many_lines
)]
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

    if d <= -1.0 && !((d - id).abs() > HYP2F1_EPS && s < 0.0) && !(neg_int_a || neg_int_b) {
        return s.powf(d) * hyp2f1(c - a, c - b, c, x);
    }
    if d <= 0.0 && x == 1.0 && !(neg_int_a || neg_int_b) {
        return f64::INFINITY; // goto hypdiv
    }

    if ax < 1.0 || x == -1.0 {
        /* 2F1(a,b;b;x) = (1-x)**(-a) */
        if (b - c).abs() < HYP2F1_EPS {
            /* b = c */
            if neg_int_b {
                y = hyp2f1_neg_c_equal_bc(a, b, x);
            } else {
                y = s.powf(-a); /* s to the -a power */
            }
            return y; // goto hypdon
        }
        if (a - c).abs() < HYP2F1_EPS {
            /* a = c */
            y = s.powf(-b); /* s to the -b power */
            return y; // goto hypdon
        }
    }

    'hypf: {
        'hypok: {
            if c <= 0.0 {
                let ic = c.round(); /* nearest integer to c */
                if (c - ic).abs() < HYP2F1_EPS {
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
                /* This transform has a pole for b-a integer, and
                 * may produce large cancellation errors for |1/x| close 1
                 */
                p = hyp2f1(a, 1.0 - c + a, 1.0 - b + a, 1.0 / x);
                q = hyp2f1(b, 1.0 - c + b, 1.0 - a + b, 1.0 / x);
                p *= (-x).powf(-a);
                q *= (-x).powf(-b);
                t1 = gamma(c);
                s = t1 * gamma(b - a) * (rgamma(b) * rgamma(c - a));
                y = t1 * gamma(a - b) * (rgamma(a) * rgamma(c - b));
                return s * p + y * q;
            } else if x < -1.0 {
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
            if (ia <= 0.0) && ((p - ia).abs() < HYP2F1_EPS) {
                /* negative int c - a */
                neg_int_ca_or_cb = true;
            }

            r = c - b;
            ib = r.round(); /* nearest integer to c-b */
            if (ib <= 0.0) && ((r - ib).abs() < HYP2F1_EPS) {
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
                (y, err) = hyt2f1(a, b, c, x);
                if err < HYP2F1_ETHRESH {
                    return y; // goto hypdon
                }
                /* Apply the recurrence if power series fails */
                let aid = (2.0 - id) as i32;
                let mut e = c + f64::from(aid);
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
        (y, _) = hyt2f1(a, b, c, x);

        // hypdon:
        return y;
    }

    /* The transformation for c-a or c-b negative integer
     * AMS55 #15.3.3
     */
    // hypf:
    (y, _) = hys2f1(c - a, c - b, c, x);
    s.powf(d) * y // goto hypdon
}
