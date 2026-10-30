//! `cerf` and `cerzo`, translated into pure Rust from `xsf/specfun/specfun.h` (xsf v0.2.2).

use core::f64::consts::PI;

use num_complex::{Complex64, c64};

use crate::xsf::complex::{cdiv, cmul};

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines
)]
pub(crate) fn cerf(z: Complex64, cer: &mut Complex64, cder: &mut Complex64) {
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

        c0 = 2.0 / PI.sqrt() * x * (-x2).exp();
        er0 = c0 * er;
        *cer = c64(er0, 0.0);
    } else {
        er = 1.0;
        r = 1.0;

        for k in 1..=12 {
            r = -r * (f64::from(k) - 0.5) / x2;
            er += r;
        }

        c0 = (-x2).exp() / (x * PI.sqrt());
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
        let er1 = (-x2).exp() * (1.0 - cs) / (2.0 * PI * x);
        let ei1 = (-x2).exp() * ss / (2.0 * PI * x);
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
    *cder = 2.0 / PI.sqrt() * cmul(-z, z).exp();
}

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cast_sign_loss
)]
pub(crate) fn cerzo(nt: i32, zo: &mut [Complex64]) {
    // ===============================================================
    // Purpose : Evaluate the complex zeros of error function erf(z)
    //           using the modified Newton's iteration method
    // Input :   NT --- Total number of zeros
    // Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
    // Routine called: CERF for computing erf(z) and erf'(z)
    // ===============================================================

    let mut it: i32;
    let mut z: Complex64;
    let mut zf = c64(0.0, 0.0);
    let mut zd = c64(0.0, 0.0);
    let mut zp: Complex64;
    let mut zw: Complex64;
    let mut zq: Complex64;
    let mut zfd: Complex64;
    let mut zgd: Complex64;
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
            zp = c64(1.0, 0.0);

            for i in 1..nr {
                zp = cmul(zp, z - zo[(i - 1) as usize]);
            }
            zfd = cdiv(zf, zp);
            zq = c64(0.0, 0.0);
            for i in 1..nr {
                zw = c64(1.0, 0.0);
                for j in 1..nr {
                    if j == i {
                        continue;
                    }
                    zw = cmul(zw, z - zo[(j - 1) as usize]);
                }
                zq += zw;
            }
            zgd = cdiv(zd - cmul(zq, zfd), zp);
            z -= cdiv(zfd, zgd);
            w0 = w;
            w = z.norm();
            if !((it <= 50) && (((w - w0) / w).abs() > 1.0e-11)) {
                break;
            }
        }
        zo[(nr - 1) as usize] = z;
    }
}
