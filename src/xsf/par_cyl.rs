//! Translated into pure Rust from `xsf/par_cyl.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! `pow(2.0, y)` is written as `y.exp2()`, which is what LLVM optimizes it to anyway.
//!
//! Unlike in xsf, `pbdv` and `pbvv` return NaN for `|v| >= 2^31`, where `(int) v` is undefined
//! behavior in C++ (in practice the allocation fails there, which also returns NaN), and the
//! other `double` -> `int` conversions saturate. Like in xsf, they also return NaN if the
//! allocation of their work buffer fails.

#![allow(clippy::excessive_precision, clippy::unreadable_literal)]

use core::f64::consts::{LN_2, PI, SQRT_2};

use num_complex::c64;

use crate::utils::{SQRT_PI, SQRT2OPI};
use crate::xsf::cephes::polevl::polevl;
use crate::xsf::cephes::{lgam, lgam_stirling_corr};
use crate::xsf::specfun::{cgama, gamma2};
use crate::{airy, airy_scaled, cospi, log1pmx, sinpi};

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn dvsa(x: f64, va: f64) -> f64 {
    // ===================================================
    // Purpose: Compute parabolic cylinder function Dv(x)
    //          for small argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PD --- Dv(x)
    // Routine called: GAMMA2 for computing Г(x)
    // ===================================================

    let eps = 1.0e-15;
    let sq2 = SQRT_2;

    let mut pd;
    let ep = (-0.25 * x * x).exp();
    let va0 = 0.5 * (1.0 - va);
    if va == 0.0 {
        pd = ep;
    } else if x == 0.0 {
        if (va0 <= 0.0) && (va0 == f64::from(va0 as i32)) {
            pd = 0.0;
        } else {
            let ga0 = gamma2(va0);
            pd = SQRT_PI / ((-0.5 * va).exp2() * ga0);
        }
    } else {
        let g1 = gamma2(-va);
        let a0 = (-0.5 * va - 1.0).exp2() * ep / g1;
        let vt = -0.5 * va;
        let g0 = gamma2(vt);
        pd = g0;
        let mut r = 1.0;
        for m in 1..=250 {
            let vm = 0.5 * (f64::from(m) - va);
            let gm = gamma2(vm);
            r = -r * sq2 * x / f64::from(m);
            let r1 = gm * r;
            pd += r1;
            if r1.abs() < pd.abs() * eps {
                break;
            }
        }
        pd *= a0;
    }
    pd
}

fn dvla(x: f64, va: f64) -> f64 {
    // ====================================================
    // Purpose: Compute parabolic cylinder functions Dv(x)
    //          for large argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PD --- Dv(x)
    // Routines called:
    //       (1) VVLA for computing Vv(x) for large |x|
    //       (2) GAMMA2 for computing Г(x)
    // ====================================================

    let eps = 1.0e-12;
    let ep = (-0.25 * x * x).exp();
    let a0 = x.abs().powf(va) * ep;
    let mut r = 1.0;
    let mut pd = 1.0;
    for k in 1..=16 {
        let kf = f64::from(k);
        r = -0.5 * r * (2.0 * kf - va - 1.0) * (2.0 * kf - va - 2.0) / (kf * x * x);
        pd += r;
        if (r / pd).abs() < eps {
            break;
        }
    }
    pd *= a0;
    if x < 0.0 {
        let x1 = -x;
        let vl = vvla(x1, va);
        let gl = gamma2(-va);
        pd = PI * vl / gl + (PI * va).cos() * pd;
    }
    pd
}

fn vvla(x: f64, va: f64) -> f64 {
    // ===================================================
    // Purpose: Compute parabolic cylinder function Vv(x)
    //          for large argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PV --- Vv(x)
    // Routines called:
    //       (1) DVLA for computing Dv(x) for large |x|
    //       (2) GAMMA2 for computing Г(x)
    // ===================================================

    let eps = 1e-12;

    let qe = (0.25 * x * x).exp();
    let a0 = x.abs().powf(-va - 1.0) * SQRT2OPI * qe;
    let mut r = 1.0;
    let mut pv = 1.0;
    for k in 1..=18 {
        let kf = f64::from(k);
        r = 0.5 * r * (2.0 * kf + va - 1.0) * (2.0 * kf + va) / (kf * x * x);
        pv += r;
        if (r / pv).abs() < eps {
            break;
        }
    }
    pv *= a0;
    if x < 0.0 {
        let x1 = -x;
        let pdl = dvla(x1, va);
        let gl = gamma2(-va);
        let dsl = (PI * va).sin() * (PI * va).sin();
        pv = dsl * gl / PI * pdl - (PI * va).cos() * pv;
    }
    pv
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn vvsa(x: f64, va: f64) -> f64 {
    // ===================================================
    // Purpose: Compute parabolic cylinder function Vv(x)
    //          for small argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PV --- Vv(x)
    // Routine called : GAMMA2 for computing Г(x)
    // ===================================================

    let eps = 1.0e-15;
    let ep = (-0.25 * x * x).exp();
    let va0 = 1.0 + 0.5 * va;

    let mut pv;
    if x == 0.0 {
        if ((va0 <= 0.0) && (va0 == f64::from(va0 as i32))) || va == 0.0 {
            pv = 0.0;
        } else {
            let vb0 = -0.5 * va;
            let sv0 = (va0 * PI).sin();
            let ga0 = gamma2(va0);
            pv = vb0.exp2() * sv0 / ga0;
        }
    } else {
        let sq2 = SQRT_2;
        let a0 = (-0.5 * va).exp2() * ep / (2.0 * PI);
        let sv = (-(va + 0.5) * PI).sin();
        let v1 = -0.5 * va;
        let g1 = gamma2(v1);
        pv = (sv + 1.0) * g1;
        let mut r = 1.0;
        let mut fac = 1.0;

        for m in 1..=250 {
            let vm = 0.5 * (f64::from(m) - va);
            let gm = gamma2(vm);
            r = r * sq2 * x / f64::from(m);
            fac = -fac;
            let gw = fac * sv + 1.0;
            let r1 = gw * r * gm;
            pv += r1;
            // Unlike in xsf (scipy/xsf#374), the convergence test ignores the weight `gw`
            // (instead of only skipping `gw == 0`), e.g. `pbvv(2.00000001, 7.0251687981058035).0`
            // was 1.8e-4 instead of 601.28.
            if (r * gm / pv).abs() < eps {
                break;
            }
        }
        pv *= a0;
    }
    pv
}

/// `pbdv` with `dv` and `dp` as `buf[..num]` and `buf[num..]`, like the C++ `malloc` block
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]
fn pbdv_impl(x: f64, mut v: f64, buf: &mut [f64], num: usize) -> (f64, f64) {
    // ====================================================
    // Purpose: Compute parabolic cylinder functions Dv(x)
    //          and their derivatives
    // Input:   x --- Argument of Dv(x)
    //          v --- Order of Dv(x)
    // Output:  DV(na) --- Dn+v0(x)
    //          DP(na) --- Dn+v0'(x)
    //          ( na = |n|, v0 = v-n, |v0| < 1,
    //            n = 0,±1,±2,… )
    //          PDF --- Dv(x)
    //          PDD --- Dv'(x)
    // Routines called:
    //       (1) DVSA for computing Dv(x) for small |x|
    //       (2) DVLA for computing Dv(x) for large |x|
    // ====================================================

    let dv = |k: i32| k as usize;
    let dp = |k: i32| num + k as usize;

    let mut pdf;
    let mut pd0 = 0.0;
    let mut pd1 = 0.0;

    let xa = x.abs();
    v += 1.0_f64.copysign(v);
    let nv = v as i32;
    let v0 = v - f64::from(nv);
    let na = nv.abs();
    let ep = (-0.25 * x * x).exp();
    let mut ja = 0;
    if na >= 1 {
        ja = 1;
    }
    if v >= 0.0 {
        if v0 == 0.0 {
            pd0 = ep;
            pd1 = x * ep;
        } else {
            for l in 0..=ja {
                let v1 = v0 + f64::from(l);
                if xa <= 5.8 {
                    pd1 = dvsa(x, v1);
                } else {
                    pd1 = dvla(x, v1);
                }
                if l == 0 {
                    pd0 = pd1;
                }
            }
        }
        buf[dv(0)] = pd0;
        buf[dv(1)] = pd1;
        for k in 2..=na {
            pdf = x * pd1 - (f64::from(k) + v0 - 1.0) * pd0;
            buf[dv(k)] = pdf;
            pd0 = pd1;
            pd1 = pdf;
        }
    } else if x <= 0.0 {
        if xa <= 5.8 {
            pd0 = dvsa(x, v0);
            let v1 = v0 - 1.0;
            pd1 = dvsa(x, v1);
        } else {
            pd0 = dvla(x, v0);
            let v1 = v0 - 1.0;
            pd1 = dvla(x, v1);
        }
        buf[dv(0)] = pd0;
        buf[dv(1)] = pd1;
        for k in 2..=na {
            let pd = (-x * pd1 + pd0) / (f64::from(k) - 1.0 - v0);
            buf[dv(k)] = pd;
            pd0 = pd1;
            pd1 = pd;
        }
    } else if x <= 2.0 {
        let mut v2 = f64::from(nv) + v0;
        if nv == 0 {
            v2 -= 1.0;
        }
        let nk = (-v2) as i32;
        let mut f1 = dvsa(x, v2);
        let v1 = v2 + 1.0;
        let mut f0 = dvsa(x, v1);
        buf[dv(nk)] = f1;
        buf[dv(nk - 1)] = f0;
        for k in (0..=nk - 2).rev() {
            let f = x * f0 + (f64::from(k) - v0 + 1.0) * f1;
            buf[dv(k)] = f;
            f1 = f0;
            f0 = f;
        }
    } else {
        if xa <= 5.8 {
            pd0 = dvsa(x, v0);
        } else {
            pd0 = dvla(x, v0);
        }
        buf[dv(0)] = pd0;
        let m = 100 + na;
        let mut f1 = 0.0;
        let mut f0 = 1e-30;
        let mut f = 0.0;
        for k in (0..=m).rev() {
            f = x * f0 + (f64::from(k) - v0 + 1.0) * f1;
            if k <= na {
                buf[dv(k)] = f;
            }
            f1 = f0;
            f0 = f;
        }
        let s0 = pd0 / f;
        for k in 0..=na {
            buf[dv(k)] *= s0;
        }
    }
    for k in 0..na {
        let v1 = v0.abs() + f64::from(k);
        if v >= 0.0 {
            buf[dp(k)] = 0.5 * x * buf[dv(k)] - buf[dv(k + 1)];
        } else {
            buf[dp(k)] = -0.5 * x * buf[dv(k)] - v1 * buf[dv(k + 1)];
        }
    }
    pdf = buf[dv(na - 1)];
    let pdd = buf[dp(na - 1)];
    (pdf, pdd)
}

/// `pbvv` with `vv` and `vp` as `buf[..num]` and `buf[num..]`, like the C++ `malloc` block
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::manual_range_contains,
    clippy::too_many_lines
)]
fn pbvv_impl(x: f64, mut v: f64, buf: &mut [f64], num: usize) -> (f64, f64) {
    // ===================================================
    // Purpose: Compute parabolic cylinder functions Vv(x)
    //          and their derivatives
    // Input:   x --- Argument of Vv(x)
    //          v --- Order of Vv(x)
    // Output:  VV(na) --- Vv(x)
    //          VP(na) --- Vv'(x)
    //          ( na = |n|, v = n+v0, |v0| < 1
    //            n = 0,±1,±2,… )
    //          PVF --- Vv(x)
    //          PVD --- Vv'(x)
    // Routines called:
    //       (1) VVSA for computing Vv(x) for small |x|
    //       (2) VVLA for computing Vv(x) for large |x|
    // ===================================================

    let vv = |k: i32| k as usize;
    let vp = |k: i32| num + k as usize;

    let mut f0 = 0.0;
    let mut f1 = 0.0;

    let xa = x.abs();
    v += 1.0_f64.copysign(v);
    let nv = v as i32;
    let v0 = v - f64::from(nv);
    let na = nv.abs();
    let qe = (0.25 * x * x).exp();
    let q2p = SQRT2OPI;
    let mut ja = 0;
    if na >= 1 {
        ja = 1;
    }
    let mut f = 0.0;
    if v <= 0.0 {
        if v0 == 0.0 {
            let pv0 = if xa <= 7.5 { vvsa(x, v0) } else { vvla(x, v0) };
            f0 = q2p * qe;
            f1 = x * f0;
            buf[vv(0)] = pv0;
            buf[vv(1)] = f0;
            buf[vv(2)] = f1;
        } else {
            for l in 0..=ja {
                let v1 = v0 - f64::from(l);
                if xa <= 7.5 {
                    f1 = vvsa(x, v1);
                } else {
                    f1 = vvla(x, v1);
                }
                if l == 0 {
                    f0 = f1;
                }
            }
            buf[vv(0)] = f0;
            buf[vv(1)] = f1;
        }
        let mut kv = 2;
        if v0 == 0.0 {
            kv = 3;
        }
        for k in kv..=na {
            f = x * f1 + (f64::from(k) - v0 - 2.0) * f0;
            buf[vv(k)] = f;
            f0 = f1;
            f1 = f;
        }
    } else if ((x >= 0.0) && (x <= 7.5)) || (xa > 7.5 && xa <= 10.0 && v <= 0.5 * x * x) {
        // Unlike in xsf (scipy/xsf#374), this branch is also used for 7.5 < |x| <= 10 and
        // v <= x^2 / 2 (here v is the order + 1), e.g. `pbvv(6.0, -8.0).0` was -1.5e263 instead
        // of -5.6633.
        let mut v2 = v;
        if v2 < 1.0 {
            v2 += 1.0;
        }
        f1 = vvsa(x, v2);
        let v1 = v2 - 1.0;
        let kv = v2 as i32;
        f0 = vvsa(x, v1);
        buf[vv(kv)] = f1;
        buf[vv(kv - 1)] = f0;
        for k in (0..=kv - 2).rev() {
            f = x * f0 - (f64::from(k) + v0 + 2.0) * f1;
            if k <= na {
                buf[vv(k)] = f;
            }
            f1 = f0;
            f0 = f;
        }
    } else if x > 7.5 {
        let pv0 = vvla(x, v0);
        let m = 100 + na.abs();
        buf[vv(1)] = pv0;
        f1 = 0.0;
        f0 = 1.0e-40;
        for k in (0..=m).rev() {
            f = x * f0 - (f64::from(k) + v0 + 2.0) * f1;
            if k <= na {
                buf[vv(k)] = f;
            }
            f1 = f0;
            f0 = f;
        }
        let s0 = pv0 / f;
        for k in 0..=na {
            buf[vv(k)] *= s0;
        }
    } else {
        if xa <= 7.5 {
            f0 = vvsa(x, v0);
            let v1 = v0 + 1.0;
            f1 = vvsa(x, v1);
        } else {
            f0 = vvla(x, v0);
            let v1 = v0 + 1.0;
            f1 = vvla(x, v1);
        }
        buf[vv(0)] = f0;
        buf[vv(1)] = f1;
        for k in 2..=na {
            f = (x * f1 - f0) / (f64::from(k) + v0);
            buf[vv(k)] = f;
            f0 = f1;
            f1 = f;
        }
    }
    for k in 0..na {
        let v1 = v0 + f64::from(k);
        if v >= 0.0 {
            buf[vp(k)] = 0.5 * x * buf[vv(k)] - (v1 + 1.0) * buf[vv(k + 1)];
        } else {
            buf[vp(k)] = -0.5 * x * buf[vv(k)] + buf[vv(k + 1)];
        }
    }
    (buf[vv(na - 1)], buf[vp(na - 1)])
}

#[allow(clippy::cast_sign_loss, clippy::float_cmp, clippy::similar_names)]
fn pbwa_impl(a: f64, x: f64) -> [f64; 4] {
    // ======================================================
    // Purpose: Compute parabolic cylinder functions W(a,±x)
    //          and their derivatives
    // Input  : a --- Parameter  ( 0 ≤ |a| ≤ 5 )
    //          x --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 )
    // Output : W1F --- W(a,x)
    //          W1D --- W'(a,x)
    //          W2F --- W(a,-x)
    //          W2D --- W'(a,-x)
    // Routine called:
    //         CGAMA for computing complex gamma function
    // ======================================================

    let mut d = [0.0; 80];
    let mut h = [0.0; 100];
    let eps = 1e-15;
    let p0 = 0.59460355750136;

    let g1;
    let g2;
    if a == 0.0 {
        g1 = 3.625609908222;
        g2 = 1.225416702465;
    } else {
        let ug = cgama(c64(0.25, 0.5 * a), 1);
        g1 = ug.norm();
        let vg = cgama(c64(0.75, 0.5 * a), 1);
        g2 = vg.norm();
    }
    let f1 = (g1 / g2).sqrt();
    let f2 = (2.0 * g2 / g1).sqrt();
    let mut h0 = 1.0;
    let mut h1 = a;
    h[0] = a;
    for l1 in 2..=100 {
        let hl = a * h1 - 0.25 * (f64::from(2 * l1) - 2.0) * (f64::from(2 * l1) - 3.0) * h0;
        h[l1 as usize - 1] = hl;
        h0 = h1;
        h1 = hl;
    }
    let mut y1f = 1.0;
    let mut r = 1.0;
    for k in 1..=100 {
        let kf = f64::from(k);
        r = 0.5 * r * x * x / (kf * (2.0 * kf - 1.0));
        let r1 = h[k as usize - 1] * r;
        y1f += r1;
        if (r1.abs() <= eps * y1f.abs()) && (k > 30) {
            break;
        }
    }
    let mut y1d = a;
    r = 1.0;
    for k in 1..100 {
        let kf = f64::from(k);
        r = 0.5 * r * x * x / (kf * (2.0 * kf + 1.0));
        let r1 = h[k as usize] * r;
        y1d += r1;
        if (r1.abs() <= eps * y1d.abs()) && (k > 30) {
            break;
        }
    }
    y1d *= x;
    let mut d1 = 1.0;
    let mut d2 = a;
    d[0] = 1.0;
    d[1] = a;
    for l2 in 3..=80 {
        let dl = a * d2 - 0.25 * (f64::from(2 * l2 - 1) - 2.0) * (f64::from(2 * l2 - 1) - 3.0) * d1;
        d[l2 as usize - 1] = dl;
        d1 = d2;
        d2 = dl;
    }
    let mut y2f = 1.0;
    r = 1.0;
    for k in 1..80 {
        let kf = f64::from(k);
        r = 0.5 * r * x * x / (kf * (2.0 * kf + 1.0));
        let r1 = d[k as usize] * r;
        y2f += r1;
        if (r1.abs() <= eps * y2f.abs()) && (k > 30) {
            break;
        }
    }
    y2f *= x;
    let mut y2d = 1.0;
    r = 1.0;
    for k in 1..80 {
        let kf = f64::from(k);
        r = 0.5 * r * x * x / (kf * (2.0 * kf - 1.0));
        let r1 = d[k as usize] * r;
        y2d += r1;
        if (r1.abs() <= eps * y2d.abs()) && (k > 30) {
            break;
        }
    }
    [
        p0 * (f1 * y1f - f2 * y2f),
        p0 * (f1 * y1d - f2 * y2d),
        p0 * (f1 * y1f + f2 * y2f),
        p0 * (f1 * y1d + f2 * y2d),
    ]
}

/// The `dv`/`dp` (or `vv`/`vp`) buffer of length `2 * num` with `num = std::abs((int) v) + 2`,
/// or `None` where `(int) v` is undefined behavior in C++, or if the allocation fails (like the
/// C++ `malloc`)
#[allow(clippy::cast_possible_truncation)]
fn alloc_buf(v: f64) -> Option<(Vec<f64>, usize)> {
    if v.abs() >= 2_147_483_648.0 {
        return None;
    }
    let num = (v as i32).unsigned_abs() as usize + 2;
    let len = num.checked_mul(2)?;
    let mut buf = Vec::new();
    buf.try_reserve_exact(len).ok()?;
    buf.resize(len, 0.0);
    Some((buf, num))
}

// Uniform asymptotic expansions of D_v(x) = U(-v - 1/2, x) for large |v| (not in xsf), following
// Gil, Segura & Temme, "Computing the real parabolic cylinder functions U(a, x), V(a, x)", ACM TOMS
// 32(1), 2006, and DLMF 12.10.

/// `u_s(t)`, `v_s(t)` (DLMF 12.10.9-13) for `s <= 10` as polynomials in `t^2` (highest degree
/// first), times `t` for odd `s`
#[rustfmt::skip]
const U_POLY: [&[f64]; 11] = [
    &[1.0],
    &[0.041666666666666664, -0.25],
    &[-0.0078125, 0.21614583333333334, 0.12586805555555555],
    &[-0.00974633487654321, 0.043858506944444445, -0.06820746527777778, -0.36650028935185186, -0.6252170138888888],
    &[0.0018274377893518518, -0.008071183569637345, -0.0038927368666409464, 1.2794328854407793, 3.079461293161651, 0.3208267274707433],
    &[0.012318482904510828, -0.0923886217838312, 0.29826426071534784, -0.542734630000888, 0.662854925623226, -5.587159132466885, -17.878599503337927, -5.084628339853797],
    &[-0.00230971554459578, 0.01713039028908537, -0.05076538089737496, 0.08601473529456106, -0.31866684067606715, 25.84962798338368, 121.71794608208658, 62.6892089202155, 2.5819438650657536],
    &[-0.03782293391770554, 0.3971408061359082, -1.8838925777686621, 5.32748236358745, -9.984394948940045, 12.980388793646037, -10.763440603293342, -146.84264579558297, -935.8851202729279, -746.193167477272, -86.57218967207372],
    &[0.007091800109569789, -0.07387291780801863, 0.34238630912780404, -0.9419515992002291, 1.7788360915042798, -2.2979420815767244, -4.788210716342768, 1078.6792540360673, 8021.524976664292, 9074.99256761248, 1991.4138217756224, 43.654373991872866],
    &[0.21514326767209896, -2.904434113573336, 18.14495655346351, -69.50879442192759, 182.34147186902612, -346.28371993293325, 490.1463712256148, -525.2468475857339, 475.3016575782134, -8889.45239528473, -76544.54715265075, -115048.85360754158, -40132.360646338035, -2526.179903176065],
    &[-0.04033936268851856, 0.541219782737624, -3.3428141089826, 12.62556001235141, -32.70922102548275, 61.506660792247686, -84.5530900042066, 92.36856692836709, -462.26787698400136, 77596.02387797361, 806760.5279375684, 1535600.7752102665, 770397.6550705286, 93261.43808375314, 1269.3042259173505],
];

#[rustfmt::skip]
const V_POLY: [&[f64]; 11] = [
    &[1.0],
    &[0.041666666666666664, 0.25],
    &[0.013020833333333334, -0.2838541666666667, -0.12413194444444445],
    &[-0.00974633487654321, 0.043858506944444445, -0.08773871527777778, 0.5749059606481481, 0.6252170138888888],
    &[-0.003045729648919753, 0.013858069902584877, -0.0038167820055298356, -1.4438563217351466, -3.33438419294946, -0.3043902864181456],
    &[0.012318482904510828, -0.0923886217838312, 0.30283285518872743, -0.5607805781707375, 0.5729826674329611, 5.264663972579894, 19.573475616622527, 5.084628339853797],
    &[0.0038495259076596336, -0.029063920602830234, 0.0903717091318846, -0.16085349184238465, 0.4309883571133622, -28.175558428156034, -129.70034337197194, -64.67370051767834, -2.5026844747880435],
    &[-0.03782293391770554, 0.3971408061359082, -1.8896668666301517, 5.3676136711748015, -10.074544460725786, 13.102329472881307, -13.249295077446435, 173.7283692713045, 984.1228658176296, 775.05056403463, 86.57218967207372],
    &[-0.011819666849282982, 0.12469748525993544, -0.5894549742489205, 1.660509942050939, -3.176285520674347, 4.199005915486087, 3.4530584725968723, -1147.4100754295257, -8448.478028408099, -9409.006632081693, -2021.895568933705, -42.91781568020086],
    &[0.21514326767209896, -2.904434113573336, 18.16268605373743, -69.68520294965315, 183.08952452219646, -348.1273549618347, 493.7529245992061, -529.0674183243831, 378.2972972123226, 8860.87318914695, 80483.0989479272, 119117.64388850002, 40974.42061406339, 2526.179903176065],
    &[0.06723227114753093, -0.9109972740490441, 5.698637542269786, -21.84760300051726, 57.50860700075974, -109.9844329920553, 155.45331145611354, -171.51515468055277, 552.0004544419786, -81630.52720143358, -842407.7602493344, -1588810.2610279743, -788089.4488024611, -94079.4112894267, -1256.8756772587144],
];

/// Taylor series in `t - 1` (highest degree first) of `zeta / (2^(1/3) (t - 1))`, `phi(zeta)`,
/// `A_1..3`, `B_0..3`, `C_0..3`, and `D_1..3` (DLMF 12.10.39-46), near the turning point `t = 1`
#[rustfmt::skip]
const TURNING_POINT: [[f64; 12]; 16] = [
    [4.122029906918741e-7, -1.0446999487501535e-6, 2.7122497343725224e-6, -7.251362109543162e-6, 2.0109234639206627e-5, -5.844030094370231e-5, 0.00018073862645291217, -0.000609853638425067, 0.002349206349206349, -0.011428571428571429, 0.1, 1.0],
    [-1.7603058148888408e-5, 3.768165559857385e-5, -8.120928973083745e-5, 0.00017646983367092497, -0.000387465244610215, 0.0008622023436729332, -0.001953542999683359, 0.004541969615373768, -0.01099482148181133, 0.02863603022593948, -0.08908987181403392, 0.8908987181403393],
    [0.00021621432944234035, -0.0003684168733650523, 0.0006190218606284738, -0.0010228417370837278, 0.0016561154879748186, -0.002614736703610845, 0.003997656268450358, -0.005857619657808733, 0.008090311998407236, -0.01022661088911089, 0.011118506493506493, -0.008645833333333333],
    [-0.0035864894779143227, 0.005070715974858761, -0.0069823771333452425, 0.00932252492324728, -0.011999636254880337, 0.014778094121550652, -0.017234566168666303, 0.01875554257408971, -0.018627521769155247, 0.016274206365427406, -0.01164729185566272, 0.005626622430831725],
    [0.06850951765253314, -0.08358431365651951, 0.09856853367739032, -0.11176315753861145, 0.12103188168616304, -0.12409849490076756, 0.11908018926746923, -0.10520923089101666, 0.08355275075441959, -0.05737133102096258, 0.0316618870574479, -0.011580906180050927],
    [3.765345480833658e-5, -7.220721948166304e-5, 0.00013796313404702427, -0.00026246489488324657, 0.0004967681809161184, -0.0009344827722691205, 0.0017449535036504898, -0.003229773415524117, 0.0059184303337401906, -0.010746859039899003, 0.019598771887253584, -0.040497462318049494],
    [-0.0007156487654069487, 0.001114491043223882, -0.0017028218574990867, 0.0025446254816655716, -0.0037045212360145886, 0.0052275940049471455, -0.0071035099177698445, 0.009214337195030392, -0.011279683854800649, 0.012851244720281004, -0.01352990226943992, 0.014635700675355903],
    [0.013287509432718636, -0.017631523436721152, 0.02274979801728682, -0.028424955166245974, 0.03421461865921327, -0.03941666392749874, 0.043099158278634614, -0.044243440711976195, 0.042040804437945556, -0.03634804543470931, 0.028298757152742497, -0.022199545918327485],
    [-0.2968177906846001, 0.3457754388127816, -0.38941339246402795, 0.4219882685412556, -0.43747386736451005, 0.4307527650287771, -0.3991549861960925, 0.34399761780737725, -0.27154671168661254, 0.19275485283032512, -0.12163465137120104, 0.07674149272374023],
    [8.743196765054662e-5, -0.00016651173062195388, 0.0003154438112200702, -0.0005936697065709325, 0.001107958505650478, -0.002044931594127975, 0.0037163840403811187, -0.006598836468245009, 0.01125611655031996, -0.017448813150602826, 0.007937005259840998, -0.07937005259840997],
    [-0.0011816571586011687, 0.0018035408082183709, -0.0026890853037427613, 0.0038987732568562293, -0.005462706553500596, 0.00733303456810808, -0.009311062666164658, 0.010955621609779732, -0.011498796954374599, 0.009751005926637022, -0.0009510984698477485, 0.009510984698477486],
    [0.018623631003904302, -0.023995351139765336, 0.029894554395694422, -0.03580012769366694, 0.04089012793287656, -0.044072926605810255, 0.04414661390145871, -0.040126484276724765, 0.031703783340347354, -0.019537081610309787, 0.0009691046990416329, -0.009691046990416329],
    [-0.36725985493893504, 0.4127352626244065, -0.4455853158679972, 0.45908216624632525, -0.44753912868448986, 0.40807289075239456, -0.34236410935742323, 0.25769967550324385, -0.16636802864431094, 0.08229673672622904, -0.002604923354379345, 0.026049233543793453],
    [-0.0002527431991470788, 0.00043053351551110183, -0.0007231793205165157, 0.001194604697365748, -0.0019337049379520828, 0.0030522984785932177, -0.004665858744642401, 0.006836303683533575, -0.0094431754665088, 0.011941777864992151, -0.01299549062049062, 0.010124007936507936],
    [0.0038828193069043717, -0.005490929386098703, 0.0075633500449262675, -0.010102459405287527, 0.013011034513137808, -0.016036486431530286, 0.018723268805834862, -0.020409037198499923, 0.02031967935018286, -0.017821446057856783, 0.01283732454368127, -0.006273704449353786],
    [-0.07230977746280635, 0.08824887926666419, -0.10411118040258356, 0.11810853717739263, -0.1279892604151968, 0.1313493268931241, -0.12619202490024678, 0.11168675305445423, -0.08892784364278655, 0.06131610707697965, -0.03408415881385055, 0.012641811383975566],
];

/// `u_s(t)` (or `v_s(t)`) from its `U_POLY` (or `V_POLY`) entry, or `i^s u_s(-i t)` if `bar`
/// (DLMF 12.10.27)
fn uv_poly(p: &[f64], s: usize, t: f64, bar: bool) -> f64 {
    let r = if bar {
        polevl(-t * t, p) * if s & 2 == 0 { 1.0 } else { -1.0 }
    } else {
        polevl(t * t, p)
    };
    if s % 2 == 1 { r * t } else { r }
}

/// `sinh(u) - u` (`hyperbolic`) or `u - sin(u)`, without cancellation for small `u >= 0`
fn sin_tail(u: f64, hyperbolic: bool) -> f64 {
    if u >= 2.0 {
        return if hyperbolic {
            u.sinh() - u
        } else {
            u - u.sin()
        };
    }
    let sign = if hyperbolic { 1.0 } else { -1.0 };
    let (mut term, mut sum) = (u, 0.0);
    for k in 1..20 {
        let k = f64::from(k);
        term *= sign * u * u / ((2.0 * k) * (2.0 * k + 1.0));
        sum += term;
    }
    sign * sum
}

/// `ln(Gamma(z + 1/2) / (sqrt(2 pi) z^z e^-z))`, for `z >= 20`; cf. DLMF 12.10.16
fn ln_gamma_ratio(z: f64) -> f64 {
    z * log1pmx(0.5 / z) + lgam_stirling_corr(z + 0.5)
}

/// `m * e^e`, without spurious over- or underflow of `e^e`
fn mul_exp(m: f64, e: f64) -> f64 {
    if e.abs() > 1400.0 {
        return m.signum() * (m.abs().ln() + e).exp();
    }
    let h = (0.5 * e).exp();
    m * h * h
}

/// `zeta`, `phi(zeta)`, and the coefficients `A_s`, `B_s`, `C_s`, `D_s` (DLMF 12.10.39-46) for
/// `s <= 4`
#[allow(clippy::needless_range_loop)]
fn airy_coefs(t: f64) -> (f64, f64, [[f64; 5]; 4]) {
    let eps = t - 1.0;
    let mut c = [[0.0; 5]; 4];
    if eps.abs() < 0.08 {
        c[0][0] = 1.0;
        c[3][0] = 1.0;
        let zeta = 2.0_f64.cbrt() * eps * polevl(eps, &TURNING_POINT[0]);
        let phi = polevl(eps, &TURNING_POINT[1]);
        let mut series = TURNING_POINT[2..].iter().map(|p| polevl(eps, p));
        for (i, row) in c.iter_mut().enumerate() {
            // A_1..3, B_0..3, C_0..3, D_1..3
            for cs in &mut row[usize::from(i % 3 == 0)..4] {
                *cs = series.next().unwrap();
            }
        }
        return (zeta, phi, c);
    }

    // 12.10.7, 12.10.23, 12.10.39, with t = cosh(theta) or cos(theta)
    let zeta = if eps > 0.0 {
        let theta = (eps + (eps * (t + 1.0)).sqrt()).ln_1p();
        (0.375 * sin_tail(2.0 * theta, true)).powf(2.0 / 3.0)
    } else {
        -(0.375 * sin_tail(2.0 * t.acos(), false)).powf(2.0 / 3.0)
    };
    let r = zeta / (eps * (t + 1.0));
    let p6 = r * r.sqrt(); // phi^6, 12.10.40

    // Fewer terms near the turning point, where 12.10.42 and 12.10.46 lose precision
    let n = if eps.abs() < 0.3 { 3 } else { 4 };
    let u: [f64; 10] = core::array::from_fn(|k| uv_poly(U_POLY[k], k, t, false));
    let v: [f64; 10] = core::array::from_fn(|k| uv_poly(V_POLY[k], k, t, false));
    // alpha_m and beta_m of 12.10.43, and phi^6m
    let (mut alpha, mut beta, mut p6k) = ([1.0; 10], [1.0; 10], [1.0; 10]);
    let mut mf = 0.0;
    for m in 1..10 {
        mf += 1.0;
        alpha[m] = alpha[m - 1] * (6.0 * mf - 5.0) * (6.0 * mf - 3.0) * (6.0 * mf - 1.0)
            / ((2.0 * mf - 1.0) * mf * 144.0);
        beta[m] = -(6.0 * mf + 1.0) / (6.0 * mf - 1.0) * alpha[m];
        p6k[m] = p6k[m - 1] * p6;
    }
    // 12.10.42, 12.10.46
    let z3 = zeta.powi(-3);
    let mut z3s = 1.0;
    for s in 0..=n {
        let mut sum = [0.0; 4];
        for m in 0..=2 * s + 1 {
            let k = 2 * s + 1 - m;
            sum[1] -= alpha[m] * p6k[k] * u[k];
            sum[2] -= beta[m] * p6k[k] * v[k];
            if m <= 2 * s {
                sum[0] += beta[m] * p6k[k - 1] * u[k - 1];
                sum[3] += alpha[m] * p6k[k - 1] * v[k - 1];
            }
        }
        c[0][s] = z3s * sum[0];
        c[1][s] = z3s * sum[1] / (zeta * zeta);
        c[2][s] = z3s * sum[2] / zeta;
        c[3][s] = z3s * sum[3];
        z3s *= z3;
    }
    (zeta, r.sqrt().sqrt(), c)
}

/// `D_v(y)`, `D_v'(y)`, `Gamma(v + 1) V(-v - 1/2, y)` and its derivative, for `v >= 20`, `y >= 0`,
/// as mantissas (with exponents `ea` for `D_v`, and `eb`), by 12.10.35-38
fn pbdv_airy(v: f64, y: f64) -> ([f64; 4], f64, f64) {
    let m = v + 0.5; // -a = mu^2 / 2
    let mu2 = 2.0 * m;
    let rm = m.sqrt();
    let y = y.min(1e100); // far beyond the turning point, where D_v(y) underflows anyway
    let t = 0.5 * y / rm;
    // ln h(mu) (12.10.15), and for t > 1 also ln h -+ mu^2 xi, with
    // mu^2 xi = (y / 2) s - m ln((y / 2 + s) / sqrt(m)) and s = sqrt(y^2 / 4 - m) (12.10.7)
    let l4 = 0.25 * (2.0 * mu2).ln();
    let lh = 0.5 * m * (m.ln() - 1.0) - l4;
    let (ea, eb) = if t > 1.0 {
        let s = ((0.5 * y - rm) * (0.5 * y + rm)).sqrt();
        let ea = -0.5 * y * s + m * (0.5 * y + s).ln() - 0.5 * m - l4;
        (ea, 2.0 * lh - ea)
    } else {
        (lh, lh)
    };
    if t > 1.0 && ea.abs() > 1000.0 {
        // D_v(y) and Ubar_v(y) under- or overflow, so only their signs matter
        return ([1.0, -1.0, 1.0, 1.0], ea, eb);
    }

    let (zeta, phi, c) = airy_coefs(t);
    let m4 = 1.0 / (mu2 * mu2);
    let [sa, sb, sc, sd] = c.map(|cs| cs.iter().rev().fold(0.0, |acc, &c| acc * m4 + c));
    let mu23 = mu2.cbrt(); // mu^(2/3)
    let mu43 = mu23 * mu23;
    let z = mu43 * zeta;
    let [ai, aip, bi, bip] = if t > 1.0 { airy_scaled(z) } else { airy(z) };
    // g(mu) / h(mu) = (1 + r) / 2 (12.10.14, 12.10.16)
    let g = 1.0 + 0.5 * ln_gamma_ratio(m).exp_m1();
    let pre = 2.0 * SQRT_PI * mu23.sqrt() * g * phi;
    let pre_d = SQRT_2 * SQRT_PI * mu23 * g / phi;
    let mu83 = mu43 * mu23 * mu23;
    let parts = [
        pre * (ai * sa + aip / mu83 * sb),
        pre_d * (ai / mu43 * sc + aip * sd),
        pre * (bi * sa + bip / mu83 * sb),
        pre_d * (bi / mu43 * sc + bip * sd),
    ];
    (parts, ea, eb)
}

/// `D_v(x)`, `D_v'(x)` for `v <= -20`, by 12.10.25 and 12.10.29, with `a = -v - 1/2`
fn pbdv_elementary(v: f64, x: f64) -> (f64, f64) {
    let a = -v - 0.5;
    let mu2 = 2.0 * a;
    let x = x.clamp(-1e100, 1e100); // the result under/overflows far beyond
    let s = (0.25 * x * x + a).sqrt();
    // gbar(mu) = (1 + 1/r) / (2 sqrt(2) mu h(mu)) (12.10.28, with r as in `ln_gamma_ratio`), and
    // e^(-mu^2 xibar) / (sqrt(2) mu h(mu) (t^2 + 1)^(1/4)) = e^l / sqrt(2 s) (12.10.15, 12.10.26)
    let ln_xs = if x >= 0.0 {
        (0.5 * x + s).ln()
    } else {
        a.ln() - (s - 0.5 * x).ln()
    };
    let l = -0.5 * x * s - a * ln_xs + 0.5 * a;
    if l.abs() > 1000.0 {
        return (mul_exp(1.0, l), mul_exp(-1.0, l));
    }

    let t = 0.5 * x / a.sqrt();
    let q = t.mul_add(t, 1.0);
    let w = 1.0 / (q * q.sqrt() * mu2);
    let (mut su, mut sv) = (0.0, 0.0);
    for k in (0..=10).rev() {
        su = su * w + uv_poly(U_POLY[k], k, t, true);
        sv = sv * w + uv_poly(V_POLY[k], k, t, true);
    }
    let g = 1.0 + 0.5 * (-ln_gamma_ratio(a)).exp_m1();
    (
        mul_exp(g * su / (2.0 * s).sqrt(), l),
        mul_exp(-g * sv * (0.5 * s).sqrt(), l),
    )
}

/// `D_v(x)`, `D_v'(x)` for `v >= 300` and `|x| sqrt(v) < 1`, by the Maclaurin series, with the
/// values at `x = 0` (12.2.6, 12.2.7, and 5.5.3), so that its zeros are exact, and
/// `w'' = (x^2 / 4 - v - 1/2) w` (12.2.2)
fn pbdv_maclaurin(v: f64, x: f64) -> (f64, f64) {
    let a = -v - 0.5;
    let lg = lgam(0.5 * v + 1.0);
    let e = 0.5 * v * LN_2 + lg; // with mantissas c0 and c1
    let c0 = cospi(0.5 * v) / SQRT_PI * (lgam(0.5 * v + 0.5) - lg).exp();
    let c1 = SQRT_2 * sinpi(0.5 * v) / SQRT_PI;
    if x.abs() < 1e-150 {
        // Only the linear terms matter; scale by 2^600 to keep x out of the subnormals
        let (k, ex) = (2.0_f64.powi(600), e - 600.0 * LN_2);
        return (
            mul_exp(c0 * k + c1 * (x * k), ex),
            mul_exp(c1 * k + a * c0 * (x * k), ex),
        );
    }
    // w = sum c_n x^n with c_{n+2} = (a c_n + c_{n-2} / 4) / ((n + 1) (n + 2))
    let mut c = [0.0, 0.0, c0, c1]; // c_{n-2}, c_{n-1}, c_n, c_{n+1}
    let (mut d, mut dp, mut xn) = (c0 + c1 * x, c1, x); // xn = x^(n+1)
    let mut prev = f64::INFINITY;
    for n in 0..60 {
        let nf = f64::from(n);
        let cn2 = (a * c[2] + 0.25 * c[0]) / ((nf + 1.0) * (nf + 2.0));
        c = [c[1], c[2], c[3], cn2];
        let term = (nf + 2.0) * cn2 * xn; // of dp
        dp += term;
        xn *= x;
        d += cn2 * xn;
        // the odd or even terms vanish for integer v
        let size = term.abs().max(prev);
        if size <= 1e-17 * dp.abs().min(d.abs() / x.abs()) {
            break;
        }
        prev = term.abs();
    }
    (mul_exp(d, e), mul_exp(dp, e))
}

/// `D_v(x)`, `D_v'(x)` by uniform asymptotic expansions, or `None` outside their regime
#[allow(clippy::cast_possible_truncation)]
fn pbdv_uniform(v: f64, x: f64) -> Option<(f64, f64)> {
    if !(v.abs() < 2_147_483_648.0 && x.is_finite()) {
        return None;
    }
    if v <= -20.0 {
        return (x >= 0.0 || x < -50.0).then(|| pbdv_elementary(v, x));
    }
    if v >= 300.0 && x.abs() * (v + 0.5).sqrt() < 1.0 {
        return Some(pbdv_maclaurin(v, x));
    }
    if x >= 0.0 {
        return (v >= 20.0 && (x >= 53.0 || v >= 300.0)).then(|| {
            let ([d, dp, ..], ea, _) = pbdv_airy(v, x);
            (mul_exp(d, ea), mul_exp(dp, ea))
        });
    }
    if v < 2.0 || (x >= -5.0 && v < 300.0) {
        return None;
    }

    // D_v(-y) = cos(pi v) D_v(y) - sin(pi v) Ubar_v(y), with Ubar_v(y) = Gamma(v + 1) V(-v - 1/2, y)
    // (12.2.15, 12.2.21)
    let y = -x;
    let (c, s) = (cospi(v), sinpi(v));
    if s == 0.0 {
        let (d, dp) = pbdv(v, y);
        return Some((c * d, -c * dp));
    }
    if v >= 20.0 {
        let ([d, dp, ub, ubp], ea, eb) = pbdv_airy(v, y);
        let r = (ea - eb).exp(); // <= 1
        return Some((
            mul_exp(c * d * r - s * ub, eb),
            mul_exp(s * ubp - c * dp * r, eb),
        ));
    }

    let (dy, dpy) = pbdv(v, y);
    let n = (20.0 - v).ceil();
    let ([.., mut u0, _], ea, eb) = pbdv_airy(v + n, y);
    if ea < -1000.0 {
        // Ubar_v(y) overflows
        let ub = mul_exp(s, eb);
        return Some((c * dy - ub, ub - c * dpy));
    }
    // Ubar_{w-1}(y) = (y Ubar_w(y) - Ubar_{w+1}(y)) / w (12.8.5), stable for decreasing w
    let ([.., u1, _], _, eb1) = pbdv_airy(v + n + 1.0, y);
    let mut u1 = mul_exp(u1, eb1 - eb);
    for k in (1..=n as i32).rev() {
        (u0, u1) = ((y * u0 - u1) / (v + f64::from(k)), u0);
    }
    // Ubar_v'(y) = y / 2 Ubar_v(y) - Ubar_{v+1}(y) (12.8.6)
    Some((
        c * dy - mul_exp(s * u0, eb),
        mul_exp(s * (0.5 * y * u0 - u1), eb) - c * dpy,
    ))
}

/// Parabolic cylinder function $W$
///
/// The function is a particular solution to the differential equation:
///
/// $$
/// y\'\' + \left( {1 \over 4} x^2 - a \right) y = 0
/// $$
///
/// Corresponds to [`scipy.special.pbwa`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbwa.html
///
/// # Arguments
/// - `a` - Real parameter
/// - `x` - Real argument
///
/// # Returns
/// - *w*: Value of the function
/// - *wp*: Value of the derivative in `x`
///
/// # See also
/// - [`pbdv`]
/// - [`pbvv`]
#[allow(clippy::manual_range_contains)]
#[must_use]
#[inline]
pub fn pbwa(a: f64, mut x: f64) -> (f64, f64) {
    /*
     * If x > 0 return w1f and w1d. Otherwise set x = abs(x) and return
     * w2f and -w2d.
     */
    if x < -5.0 || x > 5.0 || a < -5.0 || a > 5.0 {
        /*
         * The Zhang and Jin implementation only uses Taylor series;
         * return NaN outside of the range which they are accurate.
         */
        return (f64::NAN, f64::NAN);
    }
    let mut flag = false;
    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let [w1f, w1d, w2f, w2d] = pbwa_impl(a, x);
    if flag { (w2f, -w2d) } else { (w1f, w1d) }
}

/// Parabolic cylinder function $D$
///
/// Returns `(d, dp)` the parabolic cylinder function $Dv(x)$ in `d` and the derivative,
/// $Dv\'(x)$ in `dp`.
///
/// Corresponds to [`scipy.special.pbdv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbdv.html
///
/// # Arguments
/// - `v` - Real parameter
/// - `x` - Real argument
///
/// # Returns
/// - *d*: Value of the function
/// - *dp*: Value of the derivative in `x`
///
/// # See also
/// - [`pbwa`]
/// - [`pbvv`]
#[doc(alias = "pbdv_seq")]
#[must_use]
#[inline]
pub fn pbdv(v: f64, x: f64) -> (f64, f64) {
    if v.is_nan() || x.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    // Unlike in xsf (scipy/xsf#374), D_v(x) and D_{v+1}(x) are computed directly with the series
    // `dvsa` for -10 <= x < -5.8 and 0 < v < min(x^2 / 3, 20) (not near an integer v, where
    // `gamma2` is inaccurate), e.g. `pbdv(7.8219197547668315, -6.126185993244382).0` was 152.38
    // instead of 142.56.
    if v > 0.0
        && v < (x * x / 3.0).min(20.0)
        && (v - v.round()).abs() > 1e-5
        && (-10.0..-5.8).contains(&x)
    {
        let pdf = dvsa(x, v);
        return (pdf, 0.5 * x * pdf - dvsa(x, v + 1.0));
    }
    // Unlike in xsf, uniform asymptotic expansions are used for large |v| where the recurrences
    // are unstable or overflow, e.g. `pbdv(150.3, -30.0).0` was 3.1e293 instead of -2.45e144.
    if let Some(r) = pbdv_uniform(v, x) {
        return r;
    }
    /* NB. Indexing of DV/DP in specfun.f:PBDV starts from 0, hence +2 */
    let Some((mut buf, num)) = alloc_buf(v) else {
        return (f64::NAN, f64::NAN);
    };
    pbdv_impl(x, v, &mut buf, num)
}

/// Parabolic cylinder function $V$
///
/// Returns `(v, vp)` the parabolic cylinder function $Vv(x)$ in `v` and the derivative,
/// $Vv\'(x)$ in `vp`.
///
/// Corresponds to [`scipy.special.pbvv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbvv.html
///
/// # Arguments
/// - `v` - Real parameter
/// - `x` - Real argument
///
/// # Returns
/// - *v*: Value of the function
/// - *vp*: Value of the derivative in `x`
///
/// # See also
/// - [`pbwa`]
/// - [`pbdv`]
#[doc(alias = "pbvv_seq")]
#[must_use]
#[inline]
pub fn pbvv(v: f64, x: f64) -> (f64, f64) {
    if v.is_nan() || x.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    /* NB. Indexing of DV/DP in specfun.f:PBVV starts from 0, hence +2 */
    let Some((mut buf, num)) = alloc_buf(v) else {
        return (f64::NAN, f64::NAN);
    };
    pbvv_impl(x, v, &mut buf, num)
}
