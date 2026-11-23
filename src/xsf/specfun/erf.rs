//! `cerf` and `cerzo`, translated into pure Rust from `xsf/specfun/specfun.h` (scipy/xsf@65d313a).

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

    /// The sum over `n >= 1` of `term(n, ch, sh) / (n^2 + 4 x^2)`, with `ch` and `sh` the
    /// `e^(-x^2 - n^2 / 4) cosh(n y)` and `.. sinh(n y)`, until it converges
    ///
    /// Unlike in xsf, the factor `e^(-x^2)` is in the terms, and they are sums of exponentials, as
    /// `cosh(n y)` overflows (NaN) from `n y > 710`, e.g. at the 77th zero. There are more than
    /// 100 terms where needed, as they peak at `n = 2 |y|`.
    fn series(x: f64, y: f64, term: impl Fn(f64, f64, f64) -> f64) -> f64 {
        let x2 = x * x;
        let mut sum = 0.0;
        let mut w = 0.0;
        for n in 1..=(100.0 + 2.0 * y.abs()) as i32 {
            let n = f64::from(n);
            let e = |t: f64| (t - 0.25 * n * n - x2).exp();
            let (ep, em) = (e(n * y), e(-n * y));
            sum += term(n, 0.5 * (ep + em), 0.5 * (ep - em)) / (n * n + 4.0 * x2);
            if ((sum - w) / sum).abs() < EPS {
                break;
            }
            w = sum;
        }
        sum
    }

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
        let er2 = series(x, y, |n, ch, sh| {
            let e0 = (-0.25 * n * n - x2).exp();
            2.0 * x * e0 - 2.0 * x * ch * cs + n * sh * ss
        });
        c0 = FRAC_2_PI;
        let err = cer.re + er1 + c0 * er2;
        let ei2 = series(x, y, |n, ch, sh| 2.0 * x * ch * ss + n * sh * cs);
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
            // Unlike in xsf, cerzo uses the logarithmic derivative of the deflation product.
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
