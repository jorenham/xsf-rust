//! Translated into pure Rust from `xsf/par_cyl.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! `pow(2.0, y)` is written as `y.exp2()`, which is what LLVM optimizes it to anyway.
//!
//! Unlike in xsf, `pbdv` and `pbvv` return NaN for `|v| >= 2^31`, where `(int) v` is undefined
//! behavior in C++ (in practice the allocation fails there, which also returns NaN), and the
//! other `double` → `int` conversions saturate. Like in xsf, they also return NaN if the
//! allocation of their work buffer fails.

#![allow(clippy::excessive_precision, clippy::unreadable_literal)]

use core::f64::consts::{PI, SQRT_2};

use num_complex::c64;

use crate::utils::{SQRT_PI, SQRT2OPI};
use crate::xsf::specfun::{cgama, gamma2};

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
            if ((r1 / pv).abs() < eps) && (gw != 0.0) {
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
    } else if (x >= 0.0) && (x <= 7.5) {
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

#[cfg(test)]
mod tests {
    #[test]
    fn test_pbwa() {
        xsref::test("pbwa", "d_d-d_d", |x| crate::pbwa(x[0], x[1]));
    }

    #[test]
    fn test_pbdv() {
        xsref::test("pbdv", "d_d-d_d", |x| crate::pbdv(x[0], x[1]));
    }

    #[test]
    fn test_pbvv() {
        xsref::test("pbvv", "d_d-d_d", |x| crate::pbvv(x[0], x[1]));
    }
}
