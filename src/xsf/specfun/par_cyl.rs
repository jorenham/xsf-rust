//! `cpbdn`, `cpdla`, `cpdsa`, and `gaih`, translated into pure Rust from `xsf/specfun/specfun.h`
//! (scipy/xsf@370171f).
//!
//! `pow(2.0, y)` is written as `y.exp2()`, which is what LLVM optimizes it to anyway.

use core::f64::consts::SQRT_2;

use crate::utils::{SQRT_PI, SQRT_TAU};
use crate::xsf::c64::C64;

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
    // Unlike in xsf (scipy/xsf#537), e^(-z^2 / 4) instead of e^(-|z|^2 / 4), like in specfun.f
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

        if x <= 0.0 || a0 == 0.0 {
            let mut cf0 = ca0;
            cpb[0] = cf0;

            let z1 = -z;

            if a0 <= 7.0 {
                cpb[1] = cpdsa(-1, z1);
            } else {
                cpb[1] = cpdla(-1, z1);
            }

            let mut cf1 = SQRT_TAU / ca0 - cpb[1];
            cpb[1] = cf1;

            // Unlike in xsf (scipy/xsf#537), up to D_{-n0}(z), like in specfun.f
            for k in 2..=n0 {
                let cf = (-z * cf1 + cf0) / (f64::from(k) - 1.0);
                cpb[k as usize] = cf;
                cf0 = cf1;
                cf1 = cf;
            }
        } else if a0 <= 3.0 {
            // Unlike in xsf (scipy/xsf#537), the recurrence does not overwrite D_{-n0}(z) and
            // D_{-n0-1}(z) with D_0(z) and D_1(z), like in specfun.f
            let mut cfa = cpdsa(-n0, z);
            cpb[n0 as usize] = cfa;
            let n1 = n0 + 1;
            let mut cfb = cpdsa(-n1, z);
            cpb[n1 as usize] = cfb;

            let nm1 = n0 - 1;
            for k in (0..=nm1).rev() {
                let cf = z * cfa + (f64::from(k) + 1.0) * cfb;
                cpb[k as usize] = cf;
                cfb = cfa;
                cfa = cf;
            }
        } else {
            let m = 100 + n.abs();
            let mut cfa = C64::ZERO;
            let mut cfb = C64::from(1.0e-30);

            for k in (0..=m).rev() {
                let cf = z * cfb + (f64::from(k) + 1.0) * cfa;

                if k <= n0 {
                    cpb[k as usize] = cf;
                }

                cfa = cfb;
                cfb = cf;
            }

            let cs0 = ca0 / cfb;

            for c in &mut cpb[..=n0 as usize] {
                *c = cs0 * *c;
            }
        }
    }

    cpd[0] = -0.5 * z * cpb[0];

    if n >= 0 {
        for k in 1..=n as usize {
            cpd[k] = -0.5 * z * cpb[k] + k as f64 * cpb[k - 1];
        }
    } else {
        // Unlike in xsf (scipy/xsf#537), up to D_{-n0}'(z), like in specfun.f
        for k in 1..=n0 as usize {
            cpd[k] = 0.5 * z * cpb[k] - cpb[k - 1];
        }
    }
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
