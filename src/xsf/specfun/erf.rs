//! `cerf` and `cerzo`, translated into pure Rust from `xsf/cpu/specfun/specfun.h` (scipy/xsf@13bb46d).

use core::f64::consts::{FRAC_2_PI, FRAC_2_SQRT_PI, PI, TAU};

use super::deflated_newton_step;
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
        let (mut er2, mut ei2, mut w1, mut w2) = (0.0, 0.0, 0.0, 0.0);

        // Fold exp(-x^2) into the terms because cosh(n y) overflows for large y.
        // The factor exp(-n^2 / 4 + n |y| - x^2) peaks at n = 2 |y|. Target 20 terms past
        // the peak, subject to the limits below.
        // With a = -n^2 / 4 - x^2, ea = exp(a + n |y|) and em = expm1(-2 n |y|):
        // exp(a) cosh(n y) = ea (1 + em / 2) and exp(a) sinh(n y) = -sgn(y) ea em / 2.
        let ay = y.abs();
        let nlimit = 1.0e6; // Bounds the loop for huge or infinite y.
        let nmax = (2.0 * ay + 20.0).max(100.0);
        let nmax = nmax.min(nlimit) as u32;
        for n in 1..=nmax {
            let n = f64::from(n);
            let a = -0.25 * n * n - x2;
            let d = n * n + 4.0 * x2;
            let ea = (a + n * ay).exp();
            let em = (-2.0 * n * ay).exp_m1();
            let ch = ea * (1.0 + 0.5 * em);
            let sh = (-0.5 * ea * em).copysign(y);
            er2 += (2.0 * x * (a.exp() - ch * cs) + n * sh * ss) / d;
            ei2 += (2.0 * x * ch * ss + n * sh * cs) / d;
            if ((er2 - w1) / er2).abs() < EPS && ((ei2 - w2) / ei2).abs() < EPS {
                break;
            }
            (w1, w2) = (er2, ei2);
        }
        c0 = FRAC_2_PI;
        *cer = c64(cer.re + er1 + c0 * er2, ei1 + c0 * ei2);
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
            z -= deflated_newton_step(z, zf, zd, &zo[..(nr - 1) as usize]);
            w0 = w;
            w = z.norm();
            if !((it <= 50) && (((w - w0) / w).abs() > 1.0e-11)) {
                break;
            }
        }
        zo[(nr - 1) as usize] = z;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cerf_large_imaginary_part() {
        for (z, expected) in [
            (c64(20.0, 20.2), c64(11.439_694_678_497_254, -60.688_458_125_933_05)),
            (c64(50.0, 50.0), c64(0.993_553_843_558_217_2, -0.004_702_025_221_298_015)),
            (
                c64(0.510_836_925_418_312_3, 2.0),
                c64(13.614_680_673_133_714, -1.508_413_136_495_032_6),
            ),
        ] {
            let (mut actual, mut derivative) = (C64::ZERO, C64::ZERO);
            cerf(z, &mut actual, &mut derivative);
            assert!((actual - expected).norm() <= 1e-12 * expected.norm(), "{z:?}: {actual:?}");
        }
    }
}
