//! `gamma2`, `cgama`, and `psi_spec`, translated into pure Rust from `xsf/specfun/specfun.h` (xsf v0.2.2).
//!
//! `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM only does
//! when optimizing.
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior,
//! except for the loop count in `gamma2`, which is INT_MIN like with g++ on x86-64.

#![allow(clippy::excessive_precision, clippy::unreadable_literal)]

use core::f64::consts::PI;

use num_complex::{Complex64, c64};

use crate::xsf::tools::int_x86;

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
pub(crate) fn gamma2(x: f64) -> f64 {
    // ==================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x)
    //                 ( x is not equal to 0,-1,-2,…)
    // Output:  GA --- Г(x)
    // ==================================================

    const G: [f64; 26] = [
        1.0000000000000000e+00,
        0.5772156649015329e+00,
        -0.6558780715202538e+00,
        -0.4200263503409520e-01,
        0.1665386113822915e+00,
        -0.4219773455554430e-01,
        -0.9621971527877000e-02,
        0.7218943246663000e-02,
        -0.1165167591859100e-02,
        -0.2152416741149000e-03,
        0.1280502823882000e-03,
        -0.2013485478070000e-04,
        -0.1250493482100000e-05,
        0.1133027232000000e-05,
        -0.2056338417000000e-06,
        0.6116095000000000e-08,
        0.5002007500000000e-08,
        -0.1181274600000000e-08,
        0.1043427000000000e-09,
        0.7782300000000000e-11,
        -0.3696800000000000e-11,
        0.5100000000000000e-12,
        -0.2060000000000000e-13,
        -0.5400000000000000e-14,
        0.1400000000000000e-14,
        0.1000000000000000e-15,
    ];

    let mut ga;
    if x == f64::from(x as i32) {
        if x > 0.0 {
            ga = 1.0;
            let m = (x - 1.0) as i32;
            for k in 2..=m {
                ga *= f64::from(k);
            }
        } else {
            ga = 1e300;
        }
    } else {
        let mut r = 1.0;
        let mut z;
        if x.abs() > 1.0 {
            z = x.abs();
            // like g++ on x86-64, INT_MIN for |x| >= 2^31 (UB in C++), instead of saturating and
            // running ~2^31 iterations
            let m = int_x86(z);
            for k in 1..=m {
                r *= z - f64::from(k);
            }
            z -= f64::from(m);
        } else {
            z = x;
        }
        let mut gr = G[25];
        for k in (1..=25).rev() {
            gr *= z;
            gr += G[k - 1];
        }
        ga = 1.0 / (gr * z);
        if x.abs() > 1.0 {
            ga *= r;
            if x < 0.0 {
                ga = -PI / (x * ga * (PI * x).sin());
            }
        }
    }
    ga
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::similar_names
)]
pub(crate) fn cgama(mut z: Complex64, kf: i32) -> Complex64 {
    // =========================================================
    // Purpose: Compute the gamma function Г(z) or ln[Г(z)]
    //          for a complex argument
    // Input :  z  --- Complex argument
    //          kf --- Function code
    //                 kf=0 for ln[Г(z)]
    //                 kf=1 for Г(z)
    // Output:  g  --- ln[Г(z)] or Г(z)
    // ========================================================

    const A: [f64; 10] = [
        8.333333333333333e-02,
        -2.777777777777778e-03,
        7.936507936507937e-04,
        -5.952380952380952e-04,
        8.417508417508418e-04,
        -1.917526917526918e-03,
        6.410256410256410e-03,
        -2.955065359477124e-02,
        1.796443723688307e-01,
        -1.392432216905900e+00,
    ];
    let mut xx = z.re;
    let mut yy = z.im;
    let z1;
    if (yy == 0.0) && (xx <= 0.0) && (xx == f64::from(xx as i32)) {
        return c64(1e300, 0.0);
    } else if xx < 0.0 {
        z1 = z;
        z = -z;
        xx = -xx;
        yy = -yy;
    } else {
        z1 = c64(xx, 0.0);
    }
    let mut x0 = xx;
    let mut na = 0;
    if xx <= 7.0 {
        na = (7.0 - xx) as i32;
        x0 = xx + f64::from(na);
    }
    let mut az0 = c64(x0, yy).norm();
    let th = yy.atan2(x0);
    let mut gr = (x0 - 0.5) * az0.ln() - th * yy - x0 + 0.5 * (2.0 * PI).ln();
    let mut gi = th * (x0 - 0.5) + yy * az0.ln() - yy;
    for k in 1..11 {
        let t = az0.powf(f64::from(1 - 2 * k));
        gr += A[k as usize - 1] * t * ((2.0 * f64::from(k) - 1.0) * th).cos();
        gi += -A[k as usize - 1] * t * ((2.0 * f64::from(k) - 1.0) * th).sin();
    }
    if xx <= 7.0 {
        let mut gr1 = 0.0;
        let mut gi1 = 0.0;
        for j in 0..na {
            let xj = xx + f64::from(j);
            gr1 += 0.5 * (xj * xj + yy * yy).ln();
            gi1 += yy.atan2(xx + f64::from(j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if z1.re < 0.0 {
        az0 = z.norm();
        let th1 = yy.atan2(xx);
        let sr = -(PI * xx).sin() * (PI * yy).cosh();
        let si = -(PI * xx).cos() * (PI * yy).sinh();
        let az1 = c64(sr, si).norm();
        let th2 = si.atan2(sr);
        gr = (PI / (az0 * az1)).ln() - gr;
        gi = -th1 - th2 - gi;
    }
    if kf == 1 {
        let e = gr.exp();
        c64(e * gi.cos(), e * gi.sin())
    } else {
        c64(gr, gi)
    }
}

#[allow(
    clippy::approx_constant,
    clippy::cast_possible_truncation,
    clippy::excessive_precision,
    clippy::float_cmp,
    clippy::unreadable_literal
)]
pub(crate) fn psi_spec(x: f64) -> f64 {
    // ======================================
    // Purpose: Compute Psi function
    // Input :  x  --- Argument of psi(x)
    // Output:  PS --- psi(x)
    // ======================================

    const A: [f64; 8] = [
        -0.8333333333333e-01,
        0.83333333333333333e-02,
        -0.39682539682539683e-02,
        0.41666666666666667e-02,
        -0.75757575757575758e-02,
        0.21092796092796093e-01,
        -0.83333333333333333e-01,
        0.4432598039215686,
    ];
    let el = 0.5772156649015329;
    let mut ps;
    let mut s = 0.0;
    let mut xa = x.abs();

    if (x == f64::from(x as i32)) && (x <= 0.0) {
        return 1e300;
    } else if xa == f64::from(xa as i32) {
        let n = xa as i32;
        for k in 1..n {
            s += 1.0 / f64::from(k);
        }
        ps = -el + s;
    } else if (xa + 0.5) == f64::from((xa + 0.5) as i32) {
        let n = (xa - 0.5) as i32;
        for k in 1..=n {
            s += 1.0 / (2.0 * f64::from(k) - 1.0);
        }
        ps = -el + 2.0 * s - 1.386294361119891; /* 2*log(2) */
    } else {
        if xa < 10.0 {
            let n = (10.0 - f64::from(xa as i32)) as i32;
            for k in 0..n {
                s += 1.0 / (xa + f64::from(k));
            }
            xa += f64::from(n);
        }
        let x2 = 1.0 / (xa * xa);
        ps = xa.ln() - 0.5 / xa;
        ps += x2
            * (((((((A[7] * x2 + A[6]) * x2 + A[5]) * x2 + A[4]) * x2 + A[3]) * x2 + A[2]) * x2
                + A[1])
                * x2
                + A[0]);
        ps -= s;
    }
    if x < 0.0 {
        ps -= PI * (PI * x).cos() / (PI * x).sin() + 1.0 / x;
    }
    ps
}
