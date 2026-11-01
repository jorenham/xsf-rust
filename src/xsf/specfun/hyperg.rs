//! Confluent hypergeometric functions, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.
//! For example, `chgul` then returns `U(a, b, x) ~ 1` for tiny `a` (e.g. `1e-17`), instead of
//! NaN or a less accurate value (scipy/xsf#350).

use core::f64::consts::PI;

use num_complex::{Complex64, c64};

use super::hyp1f1_recurrence::{hyp1f1_b_recurrence, hyp1f1_neg_a, hyp1f1_neg_b};
use super::{cgama, gamma2, psi_spec};
use crate::utils::{SCIPY_EULER, SQRT_EPS};
use crate::xsf::complex::{cdiv, cmul, cpow, rsub};
use crate::xsf::exp::cexp;
use crate::xsf::tools::{int_x86, kahan_step};

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines
)]
pub(crate) fn cchg(mut a: f64, b: f64, mut z: Complex64) -> Complex64 {
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

    let ci = c64(0.0, 1.0);
    let mut a0 = a;
    // NOTE: `a1 = a; z0 = z;` are omitted: they are only used for the dead stores
    // `a = a1; z = z0;` at the end.
    let mut cy0 = c64(0.0, 0.0);
    let mut cy1 = c64(0.0, 0.0);
    // NOTE: uninitialized in C++, but always written before it is read (the `n` loop below
    // runs at least once).
    let mut chg = c64(0.0, 0.0);
    // NOTE: `(int)fabs(b)` is UB for |b| >= 2^31; like g++ on x86-64, `int_x86` gives INT_MIN,
    // which negates to INT_MIN (so b = -2^31 returns 1e300).
    if (b == 0.0) || (b == f64::from(int_x86(b.abs()).wrapping_neg())) {
        return c64(1e300, 0.0);
    }
    if (a == 0.0) || (z == c64(0.0, 0.0)) {
        return c64(1.0, 0.0);
    }
    if a == -1.0 {
        return rsub(1.0, z / b);
    }
    if a == b {
        return cexp(z);
    }
    if a - b == 1.0 {
        return cmul(1.0 + z / b, cexp(z));
    }
    if (a == 1.0) && (b == 2.0) {
        return cdiv(cexp(z) - 1.0, z);
    }
    if (a == f64::from(a as i32)) && (a < 0.0) {
        // NOTE: for a = -2^31, `(int)(-a)` is UB; like g++ on x86-64, `int_x86` gives INT_MIN, so
        // the loop doesn't run and the result is 1. The C++ `k < (m + 1)` is written as `1..=m`:
        // for m = INT_MAX (a = -(2^31 - 1)) `m + 1` overflows (UB), and g++ -O2 then also runs
        // ~2^31 iterations (i.e. treats it as `k <= m`), giving NaN.
        let m = int_x86(-a);
        let mut cr = c64(1.0, 0.0);
        chg = c64(1.0, 0.0);
        for k in 1..=m {
            cr = cmul(
                cr * (a + f64::from(k) - 1.0) / f64::from(k) / (b + f64::from(k) - 1.0),
                z,
            );
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
            // NOTE: `(int)a` is UB for a >= 2^31; like g++ on x86-64, `int_x86` gives INT_MIN.
            // `la + 1` overflows (and wraps, as with g++) for a in [2^31 - 1, 2^31).
            la = int_x86(a);
            a -= f64::from(la.wrapping_add(1));
        }
        // ns = 0;  (dead store)
        for n in 0..=nl {
            if a0 >= 2.0 {
                a += 1.0;
            }
            if (z.norm() < 20.0 + b.abs()) || (a < 0.0) {
                chg = c64(1.0, 0.0);
                let mut chw = c64(0.0, 0.0);
                let mut crg = c64(1.0, 0.0);
                for j in 1..501 {
                    crg = cmul(
                        crg * (a + f64::from(j) - 1.0) / (f64::from(j) * (b + f64::from(j) - 1.0)),
                        z,
                    );
                    chg += crg;
                    if cdiv(chg - chw, chg).norm() < 1e-15 {
                        break;
                    }
                    chw = chg;
                }
            } else {
                // y = 0.0;  (dead store)
                let cg1 = cgama(c64(a, 0.0), 0);
                let cg2 = cgama(c64(b, 0.0), 0);
                let cg3 = cgama(c64(b - a, 0.0), 0);
                let mut cs1 = c64(1.0, 0.0);
                let mut cs2 = c64(1.0, 0.0);
                let mut cr1 = c64(1.0, 0.0);
                let mut cr2 = c64(1.0, 0.0);
                for i in 1..=8 {
                    cr1 = cdiv(
                        -cr1 * (a + f64::from(i) - 1.0) * (a - b + f64::from(i)),
                        z * f64::from(i),
                    );
                    cr2 = cdiv(
                        cr2 * (b - a + f64::from(i) - 1.0) * (f64::from(i) - a),
                        z * f64::from(i),
                    );
                    cs1 += cr1;
                    cs2 += cr2;
                }
                let x = z.re;
                let y = z.im;
                let phi = y.atan2(x);
                let ns = if phi > -0.5 * PI { 1 } else { -1 };
                let mut cfac = cexp(f64::from(ns) * ci * PI * a);
                if y == 0.0 {
                    cfac = c64((PI * a).cos(), 0.0);
                }
                let chg1 = cmul(cmul(cmul(cexp(cg2 - cg3), cpow(z, -a)), cfac), cs1);
                let chg2 = cmul(cmul(cexp(cg2 - cg1 + z), cpow(z, a - b)), cs2);
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
                chg = (cmul(2.0 * a - b + z, cy1) + (b - a) * cy0) / a;
                cy0 = cy1;
                cy1 = chg;
                a += 1.0;
            }
            // Unlike in xsf (scipy/scipy#5349), the series (DLMF 13.2.2) replaces the forward
            // recurrence in a above if they differ by more than twice the series' error bound,
            // e.g. `hyp1f1(30, 70, 20j)` was -1748.2+51.4j instead of -0.32066+0.38180j.
            if b > 0.0 {
                if let Some((s, _)) =
                    hyp1f1_series(a0, b, z, Complex64::ONE, 0.0).filter(|&(s, err)| {
                        err.is_finite() && (!chg.is_finite() || (chg - s).norm() > 2.0 * err)
                    })
                {
                    chg = s;
                }
            }
        }
        if x0 < 0.0 {
            chg = cmul(chg, cexp(-z));
        }
    }
    // a = a1;
    // z = z0;
    chg
}

/// Unlike in xsf (scipy/xsf#365): `scale * M(a, b, z)` from the series DLMF 13.2.2 for `a > -1`
/// or `b < 0`, with an error bound for `a >= 1` and `b > 0` (an estimate for `b < 0`), or `None`
/// if it doesn't converge within 100000 terms.
///
/// The error bound (to first order in `u = EPSILON / 2`) is the sum of:
/// - the tail after the last term `t_N`: at most `|t_N| q_N / (1 - q_N)`, since the term ratios
///   `q_n = (a + n) |z| / ((n + 1) (b + n))` are non-increasing;
/// - the rounding error: at most `8 u sum (n + 1) |t_n|` (`8 n u` per term, `2 u |t_n|` per Kahan
///   step).
///
/// A negligible `t_N` is skipped if `q_N >= 1` or `b + N < |z|` (unless the sum overflowed), since
/// then the terms can grow again, e.g. for tiny `a`. Otherwise they can't, because
/// `(n + 1) (b + n) - (a + n) |z|` increases for `n >= N`.
///
/// For `a >= 0`, `b > 0` and `z >= 0` there is no cancellation, and the partial sums only
/// overflow if the result does. For `-1 < a < 0` only the leading 1 cancels.
///
/// The extra factor `e^log_scale` (for `log_scale < 0`) is applied in steps once the terms are
/// large, so that the result is finite even if `e^log_scale` underflows and `M` overflows.
fn hyp1f1_series(
    a: f64,
    b: f64,
    z: Complex64,
    scale: Complex64,
    mut log_scale: f64,
) -> Option<(Complex64, f64)> {
    let mut t = scale;
    let mut sum = scale;
    let mut comp = Complex64::ZERO;
    let mut err = t.norm();
    for n in 0..100_000 {
        let n = f64::from(n);
        t *= (a + n) / (n + 1.0) * (z / (b + n));
        kahan_step(&mut sum, &mut comp, t);
        err += (n + 2.0) * t.norm();
        if log_scale < 0.0 && err > 1e250 {
            let s = log_scale.max(-600.0);
            log_scale -= s;
            let s = s.exp();
            (t, sum, comp, err) = (t * s, sum * s, comp * s, err * s);
        }
        if t.norm() <= f64::EPSILON * sum.norm() {
            let q = (a + n + 1.0) / (n + 2.0) * (z.norm() / (b + n + 1.0));
            if (q >= 1.0 || b + n + 1.0 < z.norm()) && sum.norm().is_finite() {
                continue;
            }
            let tail = if q < 1.0 {
                t.norm() * q / (1.0 - q)
            } else {
                f64::INFINITY
            };
            let bound = 4.0 * f64::EPSILON * err + tail;
            if log_scale == 0.0 {
                return Some((sum, bound));
            }
            // in 3 steps, since `|sum| < e^710`
            let s = (log_scale.max(-1500.0) / 3.0).exp();
            return sum
                .norm()
                .is_finite()
                .then(|| (sum * s * s * s, bound * s * s * s));
        }
    }
    None
}

/// Unlike in xsf (scipy/xsf#365): the polynomial `M(-n, b, x)` for `n >= 0` and `b > 0`, using the
/// recurrence of `eval_genlaguerre_l` in `scipy/special/orthogonal_eval.pxd`, since
/// `M(-n, b, x) = L_n^(b-1)(x) / binom(n+b-1, n)` (DLMF 13.6.19).
fn chgm_poly(n: i32, b: f64, x: f64) -> f64 {
    let (mut d, mut p) = (0.0, 1.0);
    for k in 0..n {
        let k = f64::from(k);
        d = (d * k - p * x) / (k + b);
        p += d;
    }
    p
}

/// Unlike in xsf (scipy/xsf#365): `hyp1f1_series` for `M(a, b, x) e^min(x0, 0)` with `x >= 0`
/// (cf. DLMF 13.2.39), where `e^x0` is applied gradually for `x0 <= -700`, since it underflows
fn chgm_series(a: f64, b: f64, x: f64, x0: f64) -> Option<(f64, f64)> {
    let s = if x0 > -700.0 { x0.min(0.0) } else { 0.0 };
    hyp1f1_series(a, b, c64(x, 0.0), c64(s.exp(), 0.0), x0.min(0.0) - s)
        .map(|(hg, err)| (hg.re, err))
}

/// Unlike in xsf (scipy/xsf#365): `M(a, b, x)` for non-integer `b < 0` from the series, or else
/// from `hyp1f1_neg_b`, if either loses at most half of its digits.
fn chgm_neg_b(x: f64, a: f64, b: f64) -> Option<f64> {
    let (ak, xk) = if x < 0.0 { (b - a, -x) } else { (a, x) };
    let mut series = chgm_series(ak, b, xk, x);
    // if the partial sums overflow, `e^-1400 M` is summed instead
    let mut scale = 1.0;
    if series.is_some_and(|(hg, _)| !hg.is_finite()) {
        series = chgm_series(ak, b, xk, x.min(0.0) - 1400.0);
        scale = 700_f64.exp();
    }
    series
        .filter(|&(hg, err)| hg.is_finite() && err <= SQRT_EPS * hg.abs())
        .map(|(hg, _)| hg * scale * scale)
        // not for huge |a|, where the initial values are too slow (e.g. a polynomial of degree -a)
        .or_else(|| {
            (a.abs() <= 1e5)
                .then(|| hyp1f1_neg_b(a, b, x, SQRT_EPS, |b| chgm(x, a, b)))
                .flatten()
        })
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::needless_late_init,
    clippy::similar_names,
    clippy::too_many_lines
)]
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
    // NOTE: `a1 = a` is omitted: it is only used for the dead store `a = a1;` at the end.
    let x0 = x;
    let mut hg = 0.0;

    // Unlike in xsf (scipy/xsf#365), for non-integer b < 0 the series or a recurrence in b is used
    // if accurate (see `chgm_neg_b`).
    if b < 0.0 && b != b.trunc() && a.is_finite() && x.is_finite() {
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
    // Unlike in xsf (scipy/xsf#365), for a = -n (n a positive integer), b > 0 and n x > b, the
    // polynomial M(-n, b, x) is evaluated with `chgm_poly` (unless that overflows), e.g.
    // `hyp1f1(-1000, 0.5, 3)` was 5.4e29 instead of -4.1088.
    //
    // Unlike in xsf (scipy/xsf#365), this includes b - a = -n >= -1e5 for -700 < x < 0.
    let poly = x0 >= 0.0 || (x0 > -700.0 && a >= -1e5);
    if a < 0.0 && a == a.trunc() && a > f64::from(i32::MIN) && b > 0.0 && -a * x > b && poly {
        let p = chgm_poly(-a as i32, b, x) * x0.min(0.0).exp();
        if p.is_finite() {
            return p;
        }
    }
    // Unlike in xsf (scipy/xsf#365), the series DLMF 13.2.2 is summed directly for a >= 2 and
    // b > 0, and for a > -1 and x > 30 + |b|, e.g. `hyp1f1(300, 400, 1e-10)` was 1.2e87 instead
    // of 1.000000000075.
    //
    // Unlike in xsf (scipy/xsf#365), this includes x <= -700, where e^x underflows.
    //
    // For b < 0 an infinite sum isn't used, since the alternating terms (n < -b) can overflow
    // first.
    let series = (a >= 2.0 && b > 0.0) || (a > -1.0 && x > 30.0 + b.abs());
    if series && a.is_finite() && b.is_finite() && x.is_finite() {
        if let Some((hg, _)) = chgm_series(a, b, x, x0).filter(|&(hg, _)| b > 0.0 || hg.is_finite())
        {
            return hg;
        }
    }
    let mut nl = 0;
    let mut la = 0;
    if a >= 2.0 {
        // preparing terms for DLMF 13.3.1
        nl = 1;
        // NOTE: `(int)a` is UB for a >= 2^31; like g++ on x86-64, `int_x86` gives INT_MIN.
        // `la + 1` overflows (and wraps, as with g++) for a in [2^31 - 1, 2^31).
        la = int_x86(a);
        a -= f64::from(la.wrapping_add(1));
    }
    let mut y0 = 0.0;
    let mut y1 = 0.0;
    for n in 0..=nl {
        if a0 >= 2.0 {
            a += 1.0;
        }
        if (x <= 30.0 + b.abs()) || (a < 0.0) {
            hg = 1.0;
            let mut rg = 1.0;
            // Unlike in xsf (scipy/xsf#365): the sum of |rg|, to measure the cancellation
            let mut rs = 1.0;
            let mut cancellation = f64::INFINITY;
            let mut accurate = false;
            let mut converged = false;
            for j in 1..501 {
                rg = rg * (a + f64::from(j) - 1.0) / (f64::from(j) * (b + f64::from(j) - 1.0)) * x;
                hg += rg;
                rs += rg.abs();
                if (hg != 0.0) && ((rg / hg).abs() < 1e-15) {
                    accurate = rs * SQRT_EPS <= hg.abs();
                    cancellation = rs / hg.abs();
                    converged = true;
                    // DLMF 13.2.39 (cf. above)
                    if x0 < 0.0 {
                        hg *= x0.exp();
                    }
                    break;
                }
            }
            // Unlike in xsf (scipy/xsf#365), for -1 < a < 2 and b > 0 the series is summed in full
            // if it doesn't converge within 500 terms.
            if !converged && hg.is_finite() && a > -1.0 && a0 < 2.0 && b > 0.0 {
                if let Some((hg, _)) = chgm_series(a, b, x, x0) {
                    return hg;
                }
            }
            // Unlike in xsf (scipy/xsf#365), for non-integer a < -1 and b > 0, the backward
            // recurrences from Boost.Math (`hyp1f1_recurrence`) are used if the series doesn't
            // converge or loses more than half of its digits (as in Boost's checked series).
            //
            // E.g. `hyp1f1(-1000.5, 1, 1)` was 5.5e8 instead of 0.15565.
            if !accurate && a < -1.0 && a != a.trunc() && b > 0.0 && x > 0.0 {
                // NOTE: here `a` is in (-1, 2); only the (unused) error bound needs `a >= 1`.
                let m =
                    |a, b| hyp1f1_series(a, b, c64(x, 0.0), Complex64::ONE, 0.0).map(|(s, _)| s.re);
                if let Some(hg) = hyp1f1_neg_a(a, b, x, x0.min(0.0), m) {
                    return hg;
                }
            }
            // Unlike in xsf (scipy/xsf#365), if the series cancels by over 1e5, DLMF 13.3.2 is used
            // backward from b + n >= -a x / 2 (n <= 5000) if they differ by over twice its error.
            let n = (-a * x / 2.0 - b).ceil();
            let cancels = cancellation > 1e5 && (1.0..=5e3).contains(&n);
            if accurate && cancels && a < -1.0 && a != a.trunc() && b > 0.0 && x > 0.0 {
                let m = |b| chgm_series(a, b, x, x0).map_or(f64::NAN, |(hg, _)| hg);
                let rec = hyp1f1_b_recurrence(a, b, x, n as i32, true, m);
                if let Some((r, _)) = rec.filter(|&(r, err)| (r - hg).abs() > 2.0 * err * r.abs()) {
                    return r;
                }
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
                hg1 = cexp(ctb - ctba).re * x.powf(-a) * (PI * a).cos() * sum1;
                hg2 = cexp(ctb - cta + x).re * x.powf(a - b) * sum2;
            } else {
                // DLMF 13.2.39 (cf. above)
                hg1 = cexp(ctb - ctba + x0).re * x.powf(-a) * (PI * a).cos() * sum1;
                hg2 = cexp(ctb - cta).re * x.powf(a - b) * sum2;
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

#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::similar_names,
    clippy::collapsible_else_if
)]
pub(crate) fn chgu(x: f64, mut a: f64, mut b: f64, md: &mut i32, isfer: &mut i32) -> f64 {
    // =======================================================
    // Purpose: Compute the confluent hypergeometric function
    //          U(a,b,x)
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument  ( x > 0 )
    // Output:  HU --- U(a,b,x)
    //          MD --- Method code
    //          ISFER --- Error flag
    // Routines called:
    //      (1) CHGUS for small x ( MD=1 )
    //      (2) CHGUL for large x ( MD=2 )
    //      (3) CHGUBI for integer b ( MD=3 )
    //      (4) CHGUIT for numerical integration ( MD=4 )
    // =======================================================

    let mut id1: i32 = 0;
    let mut hu = 0.0;

    let aa = a - b + 1.0;
    *isfer = 0;
    // NOTE: `(int)a` (and `(int)aa`, `(int)b` below) is UB in C++ if out of `int` range (incl.
    // inf/NaN; GCC/x86 gives INT_MIN); here it saturates. Either way an out-of-range value never
    // compares equal, so the results of these comparisons are the same.
    let il1 = i32::from((a == f64::from(a as i32)) && (a <= 0.0));
    let il2 = i32::from((aa == f64::from(aa as i32)) && (aa <= 0.0));
    let il3 = i32::from((a * (a - b + 1.0)).abs() / x <= 2.0);
    let bl1 = i32::from((x <= 5.0) || (x <= 10.0 && a <= 2.0));
    let bl2 = i32::from((x > 5.0) && (x <= 12.5) && ((a >= 1.0) && (b >= a + 4.0)));
    let bl3 = i32::from((x > 12.5) && (a >= 5.0) && (b >= a + 5.0));
    let bn = i32::from((b == f64::from(b as i32)) && (b != 0.0));

    let mut id = -100;
    let mut hu1 = 0.0;
    if b != f64::from(b as i32) {
        hu = chgus(x, a, b, &mut id1);
        *md = 1;
        if id1 >= 9 {
            return hu;
        }
        hu1 = hu;
    }
    // Unlike in xsf (scipy/xsf#380), chgul is also tried if no method below applies (not for
    // |a| < 1.1e-16, where `a + 1 - 1` rounds to 0 in chgul)
    let il4 = a < 1.0 && b > a && bn == 0 && a + 1.0 != 1.0;
    if il1 != 0 || il2 != 0 || il3 != 0 || il4 {
        hu = chgul(x, a, b, &mut id);
        *md = 2;
        if id >= 9 {
            return hu;
        }
        if id1 > id {
            *md = 1;
            id = id1;
            hu = hu1;
        }
    }
    if a >= 1.0 {
        if bn != 0 && (bl1 != 0 || bl2 != 0 || bl3 != 0) {
            hu = chgubi(x, a, b, &mut id);
            *md = 3;
        } else {
            hu = chguit(x, a, b, &mut id);
            *md = 4;
        }
    } else {
        if b <= a {
            let b00 = b;
            a -= b - 1.0;
            b = 2.0 - b;
            hu = chguit(x, a, b, &mut id);
            hu *= x.powf(1.0 - b00);
            *md = 4;
        } else if bn != 0 && (!il1) != 0 {
            // NOTE: `~il1` is a bitwise NOT of 0 or 1, i.e. -1 or -2, which is always truthy, so this
            // condition is just `bn`. Reproduced as is.
            hu = chgubi(x, a, b, &mut id);
            *md = 3;
        }
        // NOTE: if none of the branches is taken (e.g. b = 0 and a < 0), `*md` is only set if
        // chgus/chgul was called (`hypu` passes an uninitialized `md`, but never reads it).
    }
    if id < 6 {
        *isfer = 6;
    }
    hu
}

#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::approx_constant
)]
pub(crate) fn chgubi(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
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

    const EL: f64 = SCIPY_EULER;

    let (a0, a1, a2, ga1, ua, ub);

    *id = -100;
    // NOTE: `(int)fabs(b - 1)` is UB in C++ if |b - 1| >= 2^31; like g++ on x86-64, `int_x86` gives
    // INT_MIN, so that the loops are skipped, instead of running ~150 * 2^31 iterations. The `int`
    // arithmetic with `n` below wraps (as with g++), since it overflows for |n| close to 2^31.
    let n = int_x86((b - 1.0).abs());
    let mut rn1 = 1.0;
    let mut rn = 1.0;
    for j in 1..=n {
        rn *= f64::from(j);
        if j == n - 1 {
            rn1 = rn;
        }
    }
    let ps = psi_spec(a);
    let ga = gamma2(a);
    if b > 0.0 {
        a0 = a;
        a1 = a - f64::from(n);
        a2 = a1;
        ga1 = gamma2(a1);
        ua = (-1.0_f64).powf(f64::from(n.wrapping_sub(1))) / (rn * ga1);
        ub = rn1 / ga * x.powf(f64::from(n.wrapping_neg()));
    } else {
        a0 = a + f64::from(n);
        a1 = a0;
        a2 = a;
        ga1 = gamma2(a1);
        ua = (-1.0_f64).powf(f64::from(n.wrapping_sub(1))) / (rn * ga) * x.powf(f64::from(n));
        ub = rn1 / ga1;
    }
    let mut hm1 = 1.0;
    let mut r = 1.0;
    let mut hmax: f64 = 0.0;
    let mut hmin: f64 = 1e300;
    let mut h0 = 0.0;
    for k in 1..=150 {
        // NOTE: `(n + k) * k` overflows `int` (C++ UB) for |b| >~ 1.43e7; `wrapping_mul` does what
        // GCC's code does (and avoids a panic in debug builds).
        r = r * (a0 + f64::from(k) - 1.0) * x / f64::from(n.wrapping_add(k).wrapping_mul(k));
        hm1 += r;
        let hu1 = hm1.abs();

        if hu1 > hmax {
            hmax = hu1;
        }
        if hu1 < hmin {
            hmin = hu1;
        }
        if (hm1 - h0).abs() < hm1.abs() * 1.0e-15 {
            break;
        }
        h0 = hm1;
    }

    let da1 = hmax.log10();
    let mut da2 = 0.0;

    if hmin != 0.0 {
        da2 = hmin.log10();
    }

    // NOTE: if `hmax` is inf (or `hm1` is NaN, so that `hmax == 0` and `da1 = -inf`), `(int)inf`
    // is UB in C++ (GCC/x86: INT_MIN, so `15 - INT_MIN` overflows to -2147483633); here it
    // saturates to i32::MAX, giving -2147483632. Same for the `(int)fabs(db1 - db2)` below.
    *id = 15 - (da1 - da2).abs() as i32;
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
    r = 1.0;
    hmax = 0.0;
    hmin = 1.0e+300;

    for k in 1..=150_i32 {
        let mut s1 = 0.0;
        let mut s2 = 0.0;

        if b > 0.0 {
            for m in 1..=k {
                s1 -= (f64::from(m) + 2.0 * a - 2.0) / (f64::from(m) * (f64::from(m) + a - 1.0));
            }
            for m in 1..=n {
                s2 += 1.0 / f64::from(k.wrapping_add(m));
            }
        } else {
            for m in 1..=k.wrapping_add(n) {
                s1 += (1.0 - a) / (f64::from(m) * (f64::from(m) + a - 1.0));
            }
            for m in 1..=k {
                s2 += 1.0 / f64::from(m);
            }
        }

        let hw = 2.0 * EL + ps + s1 - s2;
        r = r * (a0 + f64::from(k) - 1.0) * x / f64::from(n.wrapping_add(k).wrapping_mul(k));
        hm2 += r * hw;
        let hu2 = hm2.abs();

        if hu2 > hmax {
            hmax = hu2;
        }

        if hu2 < hmin {
            hmin = hu2;
        }

        // NOTE: `h0` is not reset: in the first iteration it still holds the last `hm1` of the
        // first series (before `hm1 *= log(x)`).
        if ((hm2 - h0) / hm2).abs() < 1.0e-15 {
            break;
        }

        h0 = hm2;
    }

    let db1 = hmax.log10();
    let mut db2 = 0.0;
    if hmin != 0.0 {
        db2 = hmin.log10();
    }
    let mut id1 = 15 - (db1 - db2).abs() as i32;
    if id1 < *id {
        *id = id1;
    }
    let mut hm3 = 1.0;
    if n == 0 {
        hm3 = 0.0;
    }
    r = 1.0;
    for k in 1..n {
        // NOTE: `(k - n) * k` overflows `int` (C++ UB) for |b| >~ 92682; `wrapping_mul` does what
        // GCC's code does (and avoids a panic in debug builds).
        r = r * (a2 + f64::from(k) - 1.0) / f64::from(k.wrapping_sub(n).wrapping_mul(k)) * x;
        hm3 += r;
    }
    let sa = ua * (hm1 + hm2);
    let sb = ub * hm3;
    let hu = sa + sb;
    let mut id2 = 0;
    // NOTE: `(int)log10(fabs(sa))` (and of `hu`) is UB in C++ if `sa` (`hu`) is inf or NaN
    // (GCC/x86: INT_MIN; here: i32::MAX for inf, 0 for NaN). The wrapping arithmetic below avoids
    // panics in debug builds for those saturated values.
    if sa != 0.0 {
        id1 = sa.abs().log10() as i32;
    }
    if hu != 0.0 {
        id2 = hu.abs().log10() as i32;
    }
    if sa * sb < 0.0 {
        *id = (*id).wrapping_sub(id1.wrapping_sub(id2).wrapping_abs());
    }
    hu
}

#[allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::let_and_return
)]
pub(crate) fn chguit(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
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
    *id = 9;
    // DLMF 13.4.4, integration up to C=12/X
    let a1 = a - 1.0;
    let b1 = b - a - 1.0;
    let c = 12.0 / x;
    let mut hu0 = 0.0;
    let mut hu1 = 0.0;
    for m in (10..=100).step_by(5) {
        hu1 = 0.0;
        let g = 0.5 * c / f64::from(m);
        let mut d = g;
        for _j in 1..=m {
            let mut s = 0.0;
            for k in 1..=30 {
                let t1 = d + g * T[k - 1];
                let t2 = d - g * T[k - 1];
                let f1 = (-x * t1).exp() * t1.powf(a1) * (1.0 + t1).powf(b1);
                let f2 = (-x * t2).exp() * t2.powf(a1) * (1.0 + t2).powf(b1);
                s += W[k - 1] * (f1 + f2);
            }
            hu1 += s * g;
            d += 2.0 * g;
        }
        if (1.0 - hu0 / hu1).abs() < 1.0e-9 {
            break;
        }
        hu0 = hu1;
    }
    let mut ga = gamma2(a);
    hu1 /= ga;
    // DLMF 13.4.4 with substitution t=C/(1-u)
    // integration u from 0 to 1, i.e. t from C=12/X to infinity
    let mut hu2 = 0.0;
    // Unlike in xsf (scipy/xsf#380), a non-finite integrand (0 * inf) is recomputed in log space,
    // unless b is infinite or Gamma(a) overflows (a > 171)
    let log_space = b1.is_finite() && a < 171.0;
    for m in (2..=10).step_by(2) {
        hu2 = 0.0;
        let g = 0.5 / f64::from(m);
        let mut d = g;
        for _j in 1..=m {
            let mut s = 0.0;
            for k in 1..=30 {
                let t1 = d + g * T[k - 1];
                let t2 = d - g * T[k - 1];
                let t3 = c / (1.0 - t1);
                let t4 = c / (1.0 - t2);
                let mut f1 = t3 * t3 / c * (-x * t3).exp() * t3.powf(a1) * (1.0 + t3).powf(b1);
                let mut f2 = t4 * t4 / c * (-x * t4).exp() * t4.powf(a1) * (1.0 + t4).powf(b1);
                if log_space && !f1.is_finite() {
                    f1 = t3 * t3 / c * (a1 * t3.ln() + b1 * t3.ln_1p() - x * t3).exp();
                }
                if log_space && !f2.is_finite() {
                    f2 = t4 * t4 / c * (a1 * t4.ln() + b1 * t4.ln_1p() - x * t4).exp();
                }
                s += W[k - 1] * (f1 + f2);
            }
            hu2 += s * g;
            d += 2.0 * g;
        }
        // NOTE: in the first iteration `hu0` still holds a `hu1` of the first integral (the
        // second-to-last one if that loop converged, else the last one), not divided by `ga`.
        if (1.0 - hu0 / hu2).abs() < 1.0e-9 {
            break;
        }
        hu0 = hu2;
    }
    ga = gamma2(a);
    hu2 /= ga;
    let hu = hu1 + hu2;
    hu
}

#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn chgul(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
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
    let il1 = i32::from((a == f64::from(a as i32)) && (a <= 0.0));
    let il2 = i32::from((aa == f64::from(aa as i32)) && (aa <= 0.0));
    let mut nm = 0;
    // NOTE: for a (or aa) = -2^31, `(int)fabs(a)` is UB in C++; like g++ on x86-64, `int_x86`
    // gives INT_MIN, so that the loop below doesn't run.
    if il1 != 0 {
        nm = int_x86(a.abs());
    }
    if il2 != 0 {
        nm = int_x86(aa.abs());
    }
    // IL1: DLMF 13.2.7 with k=-s-a
    // IL2: DLMF 13.2.8
    if il1 != 0 || il2 != 0 {
        hu = 1.0;
        r = 1.0;
        for k in 1..=nm {
            r = -r * (a + f64::from(k) - 1.0) * (a - b + f64::from(k)) / (f64::from(k) * x);
            hu += r;
        }
        hu *= x.powf(-a);
        *id = 10;
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
        // includes the rounding of the largest term
        //
        // NOTE: the `int` conversion is UB in C++ if the estimate is not finite; here it
        // saturates (0 for NaN).
        *id = (-(ra.max(rmax * 1e-15) / hu.abs().min(1.0)).log10()) as i32;
        hu *= x.powf(-a);
    }
    hu
}

#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::similar_names
)]
pub(crate) fn chgus(x: f64, a: f64, b: f64, id: &mut i32) -> f64 {
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

    *id = 100;
    let ga = gamma2(a);
    let gb = gamma2(b);
    let gab = gamma2(1.0 + a - b);
    let gb2 = gamma2(2.0 - b);
    let hu0 = PI / (PI * b).sin();
    let mut r1 = hu0 / (gab * gb);
    let mut r2 = hu0 * x.powf(1.0 - b) / (ga * gb2);
    let mut hu = r1 - r2;
    let mut hmax: f64 = 0.0;
    let mut hmin: f64 = 1e300;
    let mut h0 = 0.0;
    // Unlike in xsf (scipy/xsf#380), the series can't converge before j > |b - 1| (if < 150),
    // or j > 1 - b + x if b < 1 (the r1 terms grow until about then), and the id of unconverged
    // series is capped by the remaining terms
    let mut converged = false;
    let jmin = (b - 1.0).abs();
    let jconv = if b < 1.0 { jmin + x } else { jmin };
    for j in 1..151 {
        r1 = r1 * (a + f64::from(j) - 1.0) / (f64::from(j) * (b + f64::from(j) - 1.0)) * x;
        r2 = r2 * (a - b + f64::from(j)) / (f64::from(j) * (1.0 - b + f64::from(j))) * x;
        hu += r1 - r2;
        let hua = hu.abs();
        if hua > hmax {
            hmax = hua;
        }
        if hua < hmin {
            hmin = hua;
        }
        if (f64::from(j) > jconv || jmin >= 150.0) && (hu - h0).abs() < hu.abs() * 1e-15 {
            converged = true;
            break;
        }
        h0 = hu;
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
        // the last terms, and the r1 terms that still grow (up to j = 300)
        let mut rmax = (r1 - r2).abs();
        let mut j = 151.0;
        while j <= jconv && j <= 300.0 {
            r1 = r1 * (a + j - 1.0) / (j * (b + j - 1.0)) * x;
            rmax = rmax.max(r1.abs());
            j += 1.0;
        }
        *id = (*id).min((-(rmax / hu.abs()).log10()) as i32);
    }
    hu
}
