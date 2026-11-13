//! `cerf` and `cerzo`, translated into pure Rust from `xsf/specfun/specfun.h` (scipy/xsf@65d313a).

use core::f64::consts::{FRAC_2_SQRT_PI, PI, TAU};

use crate::utils::SQRT_PI;
use crate::xsf::c64::{C64, c64};

fn cerf(z: C64, cer: &mut C64, cder: &mut C64) {
    // ==========================================================
    // Purpose: Compute complex Error function erf(z) & erf'(z)
    // Input:   z   --- Complex argument of erf(z)
    //          x   --- Real part of z
    //          y   --- Imaginary part of z
    // Output:  CER --- erf(z)
    //          CDER --- erf'(z)
    // ==========================================================

    const EPS: f64 = 1.0e-12;

    let mut er: f64;
    let mut r: f64;
    let mut w: f64;
    let mut c0: f64;
    let er0: f64;

    let x = z.re;
    let y = z.im;
    let x2 = x * x;

    if x <= 3.5 {
        er = 1.0;
        r = 1.0;
        w = 0.0;

        for k in 1..=100 {
            r = r * x2 / (f64::from(k) + 0.5);
            er += r;
            if (er - w).abs() <= EPS * er.abs() {
                break;
            }
            w = er;
        }

        c0 = FRAC_2_SQRT_PI * x * (-x2).exp();
        er0 = c0 * er;
        *cer = c64(er0, 0.0);
    } else {
        er = 1.0;
        r = 1.0;

        for k in 1..=12 {
            r = -r * (f64::from(k) - 0.5) / x2;
            er += r;
        }

        c0 = (-x2).exp() / (x * SQRT_PI);
        er0 = 1.0 - c0 * er;
        *cer = c64(er0, 0.0);
    }

    if y == 0.0 {
        let err = cer.re;
        let eri = 0.0;
        *cer = c64(err, eri);
    } else {
        // NOTE: for `x == 0` (and `y != 0`), `er1` and `ei1` are 0 / 0 = NaN, so the result is
        // NaN + NaN i (in C++ as well).
        let cs = (2.0 * x * y).cos();
        let ss = (2.0 * x * y).sin();
        let er1 = (-x2).exp() * (1.0 - cs) / (TAU * x);
        let ei1 = (-x2).exp() * ss / (TAU * x);
        let mut er2 = 0.0;
        let mut w1 = 0.0;

        for n in 1..=100 {
            er2 += (-0.25 * f64::from(n) * f64::from(n)).exp() / (f64::from(n * n) + 4.0 * x2)
                * (2.0 * x - 2.0 * x * (f64::from(n) * y).cosh() * cs
                    + f64::from(n) * (f64::from(n) * y).sinh() * ss);
            if ((er2 - w1) / er2).abs() < EPS {
                break;
            }
            w1 = er2;
        }

        c0 = 2.0 * (-x2).exp() / PI;
        let err = cer.re + er1 + c0 * er2;
        let mut ei2 = 0.0;
        let mut w2 = 0.0;

        for n in 1..=100 {
            ei2 += (-0.25 * f64::from(n) * f64::from(n)).exp() / (f64::from(n * n) + 4.0 * x2)
                * (2.0 * x * (f64::from(n) * y).cosh() * ss
                    + f64::from(n) * (f64::from(n) * y).sinh() * cs);
            if ((ei2 - w2) / ei2).abs() < EPS {
                break;
            }
            w2 = ei2;
        }
        *cer = c64(err, ei1 + c0 * ei2);
    }
    *cder = FRAC_2_SQRT_PI * (-z * z).exp();
}

pub(crate) fn cerzo(nt: i32, zo: &mut [C64]) {
    // ===============================================================
    // Purpose : Evaluate the complex zeros of error function erf(z)
    //           using the modified Newton's iteration method
    // Input :   NT --- Total number of zeros
    // Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
    // Routine called: CERF for computing erf(z) and erf'(z)
    // ===============================================================

    let mut it: i32;
    let mut z: C64;
    let mut zf = C64::ZERO;
    let mut zd = C64::ZERO;
    let mut zp: C64;
    let mut zw: C64;
    let mut zq: C64;
    let mut zfd: C64;
    let mut zgd: C64;
    let mut w0: f64;
    let mut w = 0.0;

    for nr in 1..=nt {
        let pu = (PI * (4.0 * f64::from(nr) - 0.5)).sqrt();
        let pv = PI * (2.0 * f64::from(nr) - 0.25).sqrt();
        let px = 0.5 * pu - 0.5 * pv.ln() / pu;
        let py = 0.5 * pu + 0.5 * pv.ln() / pu;
        z = c64(px, py);
        it = 0;

        loop {
            it += 1;
            cerf(z, &mut zf, &mut zd);
            zp = C64::ONE;

            for i in 1..nr {
                zp *= z - zo[(i - 1) as usize];
            }
            zfd = zf / zp;
            zq = C64::ZERO;
            for i in 1..nr {
                zw = C64::ONE;
                for j in 1..nr {
                    if j == i {
                        continue;
                    }
                    zw *= z - zo[(j - 1) as usize];
                }
                zq += zw;
            }
            zgd = (zd - zq * zfd) / zp;
            z -= zfd / zgd;
            w0 = w;
            w = z.norm();
            if !((it <= 50) && (((w - w0) / w).abs() > 1.0e-11)) {
                break;
            }
        }
        zo[(nr - 1) as usize] = z;
    }
}
