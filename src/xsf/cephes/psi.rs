//! Translated into pure Rust from `xsf/cephes/psi.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
//!
//! Code for the rational approximation on [1, 2] is:
//!
//! (C) Copyright John Maddock 2006.
//! Use, modification and distribution are subject to the
//! Boost Software License, Version 1.0. (See accompanying file
//! LICENSE_1_0.txt or copy at <https://www.boost.org/LICENSE_1_0.txt>)
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf (scipy/xsf#311), the reflection for `x < 0` reduces `x` to the nearest integer
//! before evaluating `tan(pi x)`, as in the original Cephes. Otherwise, the relative error is up to
//! `~eps / (x - n)` for `x` just above a negative integer `n` (e.g. 5e-12 for `-30 < x < 0`,
//! instead of 1.3e-15).

use core::f64::consts::PI;

use crate::xsf::cephes::polevl::polevl;

// cephes/const.h
const SCIPY_EULER: f64 = 0.577_215_664_901_532_9;

const PSI_A: [f64; 7] = [
    8.333_333_333_333_333e-2,
    -2.109_279_609_279_609_4e-2,
    7.575_757_575_757_576e-3,
    -4.166_666_666_666_667e-3,
    3.968_253_968_253_968e-3,
    -8.333_333_333_333_333e-3,
    8.333_333_333_333_333e-2,
];

const PSI_Y: f32 = 0.995_581_6;

const PSI_ROOT1: f64 = 1_569_415_565.0 / 1_073_741_824.0;
const PSI_ROOT2: f64 = (381_566_830.0 / 1_073_741_824.0) / 1_073_741_824.0;
const PSI_ROOT3: f64 = 9.016_312_093_258_695e-20;

const PSI_P: [f64; 6] = [
    -2.071_332_116_774_595e-3,
    -4.525_132_144_873_906e-2,
    -2.891_912_644_477_478_4e-1,
    -6.503_185_377_089_651e-1,
    -3.255_503_118_680_449e-1,
    2.547_985_106_113_155_3e-1,
];
const PSI_Q: [f64; 7] = [
    -5.578_984_132_167_551e-7,
    2.128_498_701_782_114_6e-3,
    5.415_179_724_567_422_6e-2,
    4.359_352_969_266_596_7e-1,
    1.460_624_290_976_351_6,
    2.076_711_702_373_047,
    1.0,
];

fn digamma_imp_1_2(x: f64) -> f64 {
    /*
     * Rational approximation on [1, 2] taken from Boost.
     *
     * Now for the approximation, we use the form:
     *
     * digamma(x) = (x - root) * (Y + R(x-1))
     *
     * Where root is the location of the positive root of digamma,
     * Y is a constant, and R is optimised for low absolute error
     * compared to Y.
     *
     * Maximum Deviation Found:               1.466e-18
     * At double precision, max error found:  2.452e-17
     */
    let mut g = x - PSI_ROOT1;
    g -= PSI_ROOT2;
    g -= PSI_ROOT3;
    let r = polevl(x - 1.0, &PSI_P) / polevl(x - 1.0, &PSI_Q);

    g * f64::from(PSI_Y) + g * r
}

fn psi_asy(x: f64) -> f64 {
    let y = if x < 1.0e17 {
        let z = 1.0 / (x * x);
        z * polevl(z, &PSI_A)
    } else {
        0.0
    };

    x.ln() - (0.5 / x) - y
}

/// Psi (digamma) function
///
/// ```text
///              d      -
///   psi(x)  =  -- ln | (x)
///              dx
/// ```
///
/// is the logarithmic derivative of the gamma function.
/// For integer x,
///
/// ```text
///                   n-1
///                    -
/// psi(n) = -EUL  +   >  1/k.
///                    -
///                   k=1
/// ```
///
/// This formula is used for 0 < n <= 10.  If x is negative, it
/// is transformed to a positive argument by the reflection
/// formula  psi(1-x) = psi(x) + pi cot(pi x).
/// For general positive x, the argument is made greater than 10
/// using the recurrence  psi(x+1) = psi(x) + 1/x.
/// Then the following asymptotic expansion is applied:
///
/// ```text
///                           inf.   B
///                            -      2k
/// psi(x) = log(x) - 1/2x -   >   -------
///                            -        2k
///                           k=1   2k x
/// ```
///
/// where the B2k are Bernoulli numbers.
///
/// ACCURACY:
///
/// ```text
///    Relative error (except absolute when |psi| < 1):
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0,30        30000       3.7e-16     7.8e-17
///    IEEE      -30,0       40000       1.3e-15     1.9e-16
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///     message         condition      value returned
/// psi singularity    x integer <=0      INFINITY
/// ```
#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::if_same_then_else
)]
pub(crate) fn psi(mut x: f64) -> f64 {
    let mut y = 0.0;

    if x.is_nan() {
        return x;
    } else if x == f64::INFINITY {
        return x;
    } else if x == f64::NEG_INFINITY {
        return f64::NAN;
    } else if x == 0.0 {
        return f64::INFINITY.copysign(-x);
    } else if x < 0.0 {
        /* argument reduction before evaluating tan(pi * x) */
        let mut r = x.fract(); // std::modf(x, &q)
        if r == 0.0 {
            return f64::NAN;
        }
        // (unlike in xsf, scipy/xsf#311) reduce r to [-1/2, 1/2), where tan(pi r) doesn't lose
        // precision
        if r < -0.5 {
            r += 1.0;
        }
        y = -PI / (PI * r).tan();
        x = 1.0 - x;
    }

    /* check for positive integer up to 10 */
    if (x <= 10.0) && (x == x.floor()) {
        let n = x as i32;
        for i in 1..n {
            y += 1.0 / f64::from(i);
        }
        y -= SCIPY_EULER;
        return y;
    }

    /* use the recurrence relation to move x into [1, 2] */
    if x < 1.0 {
        y -= 1.0 / x;
        x += 1.0;
    } else if x < 10.0 {
        while x > 2.0 {
            x -= 1.0;
            y += 1.0 / x;
        }
    }
    if (1.0..=2.0).contains(&x) {
        y += digamma_imp_1_2(x);
        return y;
    }

    /* x is large, use the asymptotic series */
    y += psi_asy(x);
    y
}
