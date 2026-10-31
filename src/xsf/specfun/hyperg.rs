//! Confluent hypergeometric functions, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.
//! For example, `chgul` then returns `U(a, b, x) ~ 1` for tiny `a` (e.g. `1e-17`), instead of
//! NaN or a less accurate value.

use core::f64::consts::PI;

use num_complex::{Complex64, c64};

use super::{cgama, gamma2, psi_spec};
use crate::xsf::complex::{cdiv, cmul, cpow, rsub};
use crate::xsf::exp::cexp;
use crate::xsf::tools::{int_x86, series_eval_kahan};

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
        }
        if x0 < 0.0 {
            chg = cmul(chg, cexp(-z));
        }
    }
    // a = a1;
    // z = z0;
    chg
}

/// Unlike in xsf: `scale * M(a, b, x)` from the series DLMF 13.2.2, for `a >= 0`, `b > 0` and
/// `x >= 0`, or `None` if it doesn't converge within 100000 terms. All terms are positive, so
/// there is no cancellation, and the partial sums don't overflow unless the result does.
fn chgm_series_pos(a: f64, b: f64, x: f64, scale: f64) -> Option<f64> {
    let mut n = 0.0;
    let mut t = scale;
    let terms = || {
        t *= (a + n) / (n + 1.0) * (x / (b + n));
        n += 1.0;
        t
    };
    let (hg, k) = series_eval_kahan(terms, f64::EPSILON, 100_000, scale);
    (k > 0).then_some(hg)
}

/// Unlike in xsf: `M(a, b, x)` for `a <= -1`, `b > 0` and `x > 0`, using the recurrence
/// DLMF 13.3.1 downward in `a`, starting from `M(a0, b, x)` and `M(a0 + 1, b, x)` with
/// `a0 = a - floor(a)` (exact, because `|a| >= 1`), together with an estimate of its relative
/// error in units of eps: `steps^2 max(1, |M(a + 1, b, x) / M(a, b, x)|)` (as observed; the
/// second factor accounts for a nearby zero of `M` in `a`). The steps from `a0` and
/// `a0 - 1` divide by `b - a0` and `b - a0 + 1`, which can be tiny (the recurrence is singular at
/// `a = b`), so `None` is returned if their local rounding error can exceed 1e4 eps times the
/// previous value. It is also `None` if the result isn't finite, or (partly) underflows.
fn chgm_recur_neg_a(a: f64, b: f64, x: f64) -> Option<(f64, f64)> {
    let mut ak = a - a.floor();
    let steps = ak - a;
    let mut y1 = chgm_series_pos(ak, b, x, 1.0)?;
    let mut y2 = chgm_series_pos(ak + 1.0, b, x, 1.0)?;
    while ak > a {
        // (b - a) M(a - 1) = a M(a + 1) - (2a - b + x) M(a)
        let p = ak * y2;
        let q = (2.0 * ak - b + x) * y1;
        if ak > -1.0 && p.abs() + q.abs() > 1e4 * ((b - ak) * y1).abs() {
            return None;
        }
        let y0 = (p - q) / (b - ak);
        y2 = y1;
        y1 = y0;
        ak -= 1.0;
    }
    let err = steps * steps * (y2 / y1).abs().max(1.0);
    y1.is_normal().then_some((y1, err))
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::needless_late_init,
    clippy::similar_names
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

    //  DLMF 13.2.39
    if x < 0.0 {
        a = b - a;
        a0 = a;
        x = x.abs();
    }
    // Unlike in xsf, for a >= 2 and b > 0 the series DLMF 13.2.2 is summed directly, instead of
    // using the forward recurrence in a (DLMF 13.3.1) below, which is unstable for b > 2a + x
    // (roughly). E.g. `hyp1f1(300, 400, 1e-10)` was 1.2e87 instead of 1.000000000075.
    if a >= 2.0 && b > 0.0 && a.is_finite() && b.is_finite() && x.is_finite() && x0 > -700.0 {
        if let Some(hg) = chgm_series_pos(a, b, x, x0.min(0.0).exp()) {
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
            // Unlike in xsf, the cancellation in the series is measured (see below)
            let mut ha = 1.0;
            let mut hg_raw = hg;
            let mut hg_tol = f64::INFINITY;
            for j in 1..501 {
                rg = rg * (a + f64::from(j) - 1.0) / (f64::from(j) * (b + f64::from(j) - 1.0)) * x;
                hg += rg;
                ha += rg.abs();
                if (hg != 0.0) && ((rg / hg).abs() < 1e-15) {
                    hg_raw = hg;
                    // bound on the rounding error: 8 roundings per term ratio, 1 per addition
                    hg_tol = 8.0 * f64::from(j) * f64::EPSILON * ha;
                    // DLMF 13.2.39 (cf. above)
                    if x0 < 0.0 {
                        hg *= x0.exp();
                    }
                    break;
                }
            }
            // Unlike in xsf, for a <= -1 and b > 0, if the series did not converge, or lost more
            // than 4 digits, then the downward recurrence in a (DLMF 13.3.1) is used instead if
            // its estimated error is smaller, and it agrees with the series to within the bound
            // on its rounding error. The recurrence is unstable if x is much larger than b.
            // E.g. `hyp1f1(-1000.5, 1, 1)` was 5.5e8 instead of 0.15565.
            let canc = if hg_tol.is_finite() {
                ha / hg_raw.abs()
            } else {
                f64::INFINITY
            };
            if a <= -1.0 && b > 0.0 && x > 0.0 && x <= b + 5.0 && -a < 1e5 && canc > 1e4 {
                if let Some((h, _)) = chgm_recur_neg_a(a, b, x)
                    .filter(|&(h, err)| err < canc && (h - hg_raw).abs() <= hg_tol)
                {
                    hg = if x0 < 0.0 { h * x0.exp() } else { h };
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
    // NOTE: `(int)a` (and `(int)aa`, `(int)b` below) is UB in C++ if out of `int` range (incl. inf/NaN;
    // GCC/x86 gives INT_MIN); here it saturates. Either way an out-of-range value never compares equal,
    // so the results of these comparisons are the same.
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
    if il1 != 0 || il2 != 0 || il3 != 0 {
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
        // NOTE: if none of the branches is taken (a < 1, b > a, b not a nonzero integer, e.g. b = 0 and
        // a < 0), `*md` is only set if chgus/chgul was called, and otherwise left untouched (`hypu` passes
        // an uninitialized `md`, but never reads it).
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

    const EL: f64 = 0.5772156649015329;

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
        // NOTE: `(n + k) * k` overflows `int` (C++ UB) if n > about 1.43e7, i.e. |b| > about 1.43e7;
        // `wrapping_mul` reproduces what GCC's code does (and avoids a panic in debug builds).
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

    // NOTE: if `hmax` is inf (or `hm1` is NaN, so that `hmax == 0` and `da1 = -inf`), `(int)inf` is UB
    // in C++ (GCC/x86: INT_MIN, so `15 - INT_MIN` overflows to -2147483633); here it saturates to
    // i32::MAX, giving -2147483632. Same for the `(int)fabs(db1 - db2)` below.
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

        // NOTE: `h0` is not reset: in the first iteration it still holds the last `hm1` of the first
        // series (before `hm1 *= log(x)`).
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
        // NOTE: `(k - n) * k` overflows `int` (C++ UB) if n > about 92681, i.e. |b| > about 92682;
        // `wrapping_mul` reproduces what GCC's code does (and avoids a panic in debug builds).
        r = r * (a2 + f64::from(k) - 1.0) / f64::from(k.wrapping_sub(n).wrapping_mul(k)) * x;
        hm3 += r;
    }
    let sa = ua * (hm1 + hm2);
    let sb = ub * hm3;
    let hu = sa + sb;
    let mut id2 = 0;
    // NOTE: `(int)log10(fabs(sa))` (and of `hu`) is UB in C++ if `sa` (`hu`) is inf or NaN (GCC/x86:
    // INT_MIN; here: i32::MAX for inf, 0 for NaN). The wrapping arithmetic below avoids panics in debug
    // builds for those saturated values.
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
                let f1 = t3 * t3 / c * (-x * t3).exp() * t3.powf(a1) * (1.0 + t3).powf(b1);
                let f2 = t4 * t4 / c * (-x * t4).exp() * t4.powf(a1) * (1.0 + t4).powf(b1);
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
    // NOTE: for a (or aa) = -2147483648.0 exactly, `(int)fabs(a)` is UB in C++; like g++ on x86-64,
    // `int_x86` gives INT_MIN, so that the loop below doesn't run. For very negative integer a (or aa)
    // the loop runs |a| (|aa|) times.
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
        for k in 1..=25 {
            r = -r * (a + f64::from(k) - 1.0) * (a - b + f64::from(k)) / (f64::from(k) * x);
            ra = r.abs();
            if ((k > 5) && (ra >= r0)) || (ra < 1e-15) {
                break;
            }
            r0 = ra;
            hu += r;
        }
        // NOTE: `(int)fabs(log10(ra))` is UB in C++ if `ra` is 0 or inf (log10 -> -inf/inf), or NaN
        // (GCC/x86: INT_MIN); here it saturates to i32::MAX (0 for NaN), so `chgu` then returns this
        // result (`id >= 9`) where the C++ continues. `ra == 0` happens e.g. if `a + k - 1.0` rounds to
        // zero for a non-integer `a` close to `1 - k` (e.g. a = 1e-17).
        *id = ra.log10().abs() as i32;
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
        if (hu - h0).abs() < hu.abs() * 1e-15 {
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
    hu
}
