//! `cpbdn`, `cpdla`, `cpdsa`, and `gaih`, translated into pure Rust from `xsf/specfun/specfun.h`
//! (scipy/xsf@7a006dc).
//!
//! `pow(2.0, y)` is written as `y.exp2()`, which is what LLVM optimizes it to anyway.

use core::f64::consts::{LN_2, SQRT_2};

use crate::utils::{SQRT_PI, SQRT_TAU_DOWN, pow2};
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::ldexp_f64;
use crate::xsf::cephes::lgam;

/// The minimum order `K` of the seeds `D_{-K}` and `D_{-K-1}` of the backward recurrence
const SEED_MIN_ORDER: usize = 40;
/// For `|Im z| > Re z` and `|z|^2 < SEED_NEAR * K`, `K` is raised to `SEED_FAR * |z|^2`, away from
/// the turning points `z^2 = -4K`, where the two saddle points of `seeds` merge
const SEED_NEAR: f64 = 10.0;
/// See [`SEED_NEAR`]
const SEED_FAR: f64 = 2.0;
/// The step of `seeds` in units of the width `1 / sqrt|K + t_s^2|` of the integrand
const SEED_STEP: f64 = 0.5;
/// The relative size of the last node of `seeds` in each direction
const SEED_TOL: f64 = 1e-18;
/// The maximum number of nodes of `seeds` in each direction
const SEED_MAX_NODES: u32 = 100;
/// The binary exponent of the factor by which the backward recurrence is rescaled
const RESCALE_EXP: i32 = 512;

pub(crate) fn cpbdn(n: i32, z: C64, cpb: &mut [C64], cpd: &mut [C64]) {
    // ==================================================
    // Purpose: Compute the parabolic cylinder functions
    //           Dn(z) and Dn'(z) for a complex argument
    // Input:   z --- Complex argument of Dn(z)
    //          n --- Order of Dn(z)  ( n=0,±1,±2,… )
    // Output:  CPB(|n|) --- Dn(z)
    //          CPD(|n|) --- Dn'(z)
    // Routines called:
    //      (1) CPDSA for computing Dn(z) for a small |z|
    //      (2) CPDLA for computing Dn(z) for a large |z|
    // ==================================================

    let x = z.re;
    let a0 = z.norm();
    let ca0 = (-0.25 * z * z).exp();
    let mut n0 = 0;

    if n >= 0 {
        let mut cf0 = ca0;
        let mut cf1 = z * ca0;

        cpb[0] = cf0;
        cpb[1] = cf1;

        for k in 2..=n {
            let cf = z * cf1 - (f64::from(k) - 1.0) * cf0;
            cpb[k as usize] = cf;
            cf0 = cf1;
            cf1 = cf;
        }
    } else {
        n0 = -n;

        if x < 0.0 {
            let mut cf0 = ca0;
            cpb[0] = cf0;

            let z1 = -z;

            if a0 <= 7.0 {
                cpb[1] = cpdsa(-1, z1);
            } else {
                cpb[1] = cpdla(-1, z1);
            }

            let mut cf1 = SQRT_TAU_DOWN / ca0 - cpb[1];
            cpb[1] = cf1;

            for k in 2..=n0 {
                let cf = (-z * cf1 + cf0) / (f64::from(k) - 1.0);
                cpb[k as usize] = cf;
                cf0 = cf1;
                cf1 = cf;
            }
        } else {
            // Unlike in xsf, a seeded backward recurrence, instead of `cpdsa` for |z| <= 3 and
            // else Miller's algorithm from order n0 + 100, which can be inaccurate or overflow
            pbdn_neg(z, ca0, &mut cpb[..=n0 as usize]);
        }
    }

    cpd[0] = -0.5 * z * cpb[0];

    if n >= 0 {
        for k in 1..=n as usize {
            cpd[k] = -0.5 * z * cpb[k] + k as f64 * cpb[k - 1];
        }
    } else {
        for k in 1..=n0 as usize {
            cpd[k] = 0.5 * z * cpb[k] - cpb[k - 1];
        }
    }
}

/// `D_{-k}(z)` for `k = 0, ..., cpb.len() - 1` and `Re z >= 0` from the backward recurrence
/// `D_{-k+1}(z) = z D_{-k}(z) + k D_{-k-1}(z)` (<https://dlmf.nist.gov/12.8.E1>), seeded by `seeds`
/// and normalized by `D_0(z) = ca0 = exp(-z^2/4)` to cancel the errors that they share
fn pbdn_neg(z: C64, ca0: C64, cpb: &mut [C64]) {
    let n0 = cpb.len() - 1;
    let r2 = z.norm_sqr();
    let mut kk = n0.max(SEED_MIN_ORDER);
    if z.im.abs() > z.re && r2 < SEED_NEAR * kk as f64 {
        kk = kk.max((SEED_FAR * r2).ceil() as usize);
    }
    let ldexp = |c: C64, e: f64| c64(ldexp_f64(c.re, e as i32), ldexp_f64(c.im, e as i32));
    let (big, small) = (pow2(RESCALE_EXP), pow2(-RESCALE_EXP));

    // D_{-k} and D_{-k-1} are `cfb 2^e` and `cfa 2^e`
    let (mut cfb, mut cfa, mut e) = seeds(kk as f64, z);
    if kk == n0 {
        cpb[n0] = ldexp(cfb, e);
    }
    for k in (0..kk).rev() {
        let mut cf = z * cfb + (k as f64 + 1.0) * cfa;
        if cf.l1_norm() > big {
            (cf, cfb) = (cf * small, cfb * small);
            e += f64::from(RESCALE_EXP);
        }
        if k <= n0 {
            cpb[k] = ldexp(cf, e);
        }
        (cfa, cfb) = (cfb, cf);
    }

    if cpb[0].norm().is_normal() && ca0.norm().is_normal() {
        let cs = ca0 / cpb[0];
        for c in cpb {
            *c = cs * *c;
        }
    }
}

/// `D_{-K}(z)` and `D_{-K-1}(z)` for `Re z >= 0`, as `(m0, m1, e)` for `m0 2^e` and `m1 2^e`
///
/// This is the trapezoidal rule for the integral <https://dlmf.nist.gov/12.5.E1> on the ray
/// `t = t_s e^v` through the saddle point `t_s^2 + z t_s = K` of `t^K exp(-t^2/2 - z t)`,
/// truncated where the integrand drops below `SEED_TOL` times the sum. Near the imaginary axis
/// the integrand grows again farther out on the ray (e.g. for K = 40 and z = 20i), so there the
/// truncated sum stands for a contour that leaves the ray for a sector where the integrand
/// decays. That error is not bounded here; on a test map of Re z >= 0, the normalized results
/// are within 5e-14 of mpmath.
fn seeds(kk: f64, z: C64) -> (C64, C64, f64) {
    // t_s = (sqrt(z^2 + 4K) - z) / 2 without cancellation, where Im sqrt has the sign of Im z,
    // also on the imaginary axis past the turning points (the limit from Re z > 0)
    let z2 = z * z;
    let s = (z2 + 4.0 * kk).sqrt();
    let ts = 2.0 * kk / (z + c64(s.re, s.im.copysign(z.im)));
    let a = ts * ts;
    let a_half = 0.5 * a;
    let h = SEED_STEP / (kk + a).norm().sqrt();

    // the integrands of D_{-K} and D_{-K-1} relative to their values at v = 0
    let (mut s0, mut s1) = (C64::ONE, C64::ONE);
    for dir in [1.0, -1.0] {
        for j in 1..=SEED_MAX_NODES {
            let v = dir * f64::from(j) * h;
            let em = v.exp_m1();
            let w = (-kk * (em - v) - a_half * em * em).exp();
            s0 += w;
            s1 += w * (em + 1.0);
            if !(w.norm() >= SEED_TOL * s0.norm()) {
                break;
            }
        }
    }
    // the integrand at v = 0, and 1 / Gamma(K), as 2^e times a number of order 1
    let lp = -0.25 * z2 + kk * ts.ln() - a_half - z * ts - lgam(kk);
    let e = (lp.re / LN_2).round();
    let cs = (lp - e * LN_2).exp() * h;
    (cs * s0, cs * s1 * ts / kk, e)
}

fn cpdla(n: i32, z: C64) -> C64 {
    // ===========================================================
    // Purpose: Compute complex parabolic cylinder function Dn(z)
    //          for large argument
    // Input:   z   --- Complex argument of Dn(z)
    //          n   --- Order of Dn(z) (n = 0,±1,±2,…)
    // Output:  CDN --- Dn(z)
    // ===========================================================

    let nf = f64::from(n);
    let cb0 = z.powi(n) * (-0.25 * z * z).exp();
    let mut cr = C64::ONE;
    let mut cdn = C64::ONE;
    for k in 1..=16 {
        let kf = f64::from(k);
        cr = -0.5 * cr * (2.0 * kf - nf - 1.0) * (2.0 * kf - nf - 2.0) / (kf * z * z);
        cdn += cr;
        if cr.norm() < cdn.norm() * 1e-12 {
            break;
        }
    }
    cdn * cb0
}

fn cpdsa(n: i32, z: C64) -> C64 {
    // ===========================================================
    // Purpose: Compute complex parabolic cylinder function Dn(z)
    //          for small argument
    // Input:   z   --- Complex argument of D(z)
    //          n   --- Order of D(z) (n = 0,-1,-2,...)
    // Output:  CDN --- Dn(z)
    // Routine called: GAIH for computing Г(x), x=n/2 (n=1,2,...)
    // ===========================================================

    let eps = 1.0e-15;
    let nf = f64::from(n);

    let mut cdn;
    let ca0 = (-0.25 * z * z).exp();
    let va0 = 0.5 * (1.0 - nf);
    if n == 0 {
        cdn = ca0;
    } else if z.norm() == 0.0 {
        if (va0 <= 0.0) && (va0 == f64::from(va0 as i32)) {
            cdn = C64::ZERO;
        } else {
            let ga0 = gaih(va0);
            let pd = SQRT_PI / ((-0.5 * nf).exp2() * ga0);
            cdn = C64::from(pd);
        }
    } else {
        let xn = -nf;
        let g1 = C64::from(gaih(xn));
        let cb0 = (-0.5 * nf - 1.0).exp2() * ca0 / g1;
        let vt = -0.5 * nf;
        let g0 = C64::from(gaih(vt));
        cdn = g0;
        let mut cr = C64::ONE;

        for m in 1..=250 {
            let vm = 0.5 * (f64::from(m) - nf);
            let gm = C64::from(gaih(vm));
            cr = -cr * SQRT_2 * z / f64::from(m);
            let cdw = gm * cr;
            cdn += cdw;
            if cdw.norm() < cdn.norm() * eps {
                break;
            }
        }
        cdn *= cb0;
    }
    cdn
}

fn gaih(x: f64) -> f64 {
    // =====================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x), x = n/2, n=1,2,…
    // Output:  GA --- Г(x)
    // =====================================================

    let mut ga;
    if (x == f64::from(x as i32)) && (x > 0.0) {
        ga = 1.0;
        let m = (x - 1.0) as i32;
        for k in 2..=m {
            ga *= f64::from(k);
        }
    } else if ((x + 0.5) == f64::from((x + 0.5) as i32)) && (x > 0.0) {
        let m = x as i32;
        ga = SQRT_PI;
        for k in 1..=m {
            ga *= 0.5 * (2.0 * f64::from(k) - 1.0);
        }
    } else {
        ga = f64::NAN;
    }
    ga
}
