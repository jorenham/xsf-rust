//! `itstruve0`, `it2struve0`, `itmodstruve0`, `struve_h`, and `struve_l` are translated into pure
//! Rust from `xsf/cpu/struve.h` (scipy/xsf@e9eb557).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM only does
//! when optimizing.

#![allow(clippy::needless_late_init)]

use core::f64::consts::{FRAC_2_PI, FRAC_PI_2, PI, TAU};

use crate::xsf::bessel::{IT0_ASYMP, bessel_j_miller, itj0y0_asymp};
use crate::xsf::cephes;
use crate::xsf::specfun::specfun_convinf;

/// `(int_0^x H0(t) dt, int_0^x H0(t)/t dt)` for `0.5 <= x < 40`, by integrating DLMF 11.4.21,
/// `H0 = 4/pi sum_k J_{2k+1} / (2k+1)`, term by term, with Miller's backward recurrence for `J_n`.
fn itsh0_neumann(x: f64) -> (f64, f64) {
    // sum_{2j > n} J_{2j}, sum_{2j+1 >= n} J_{2j+1}, and the two series
    let (mut te, mut to, mut se, mut so) = (0.0, 0.0, 0.0, 0.0);
    let norm = bessel_j_miller(x, |k, j| {
        let n = f64::from(k);
        if k % 2 == 1 {
            to += j;
            so += (2.0 * to - j) / (n * n);
            se += te / n;
        } else {
            te += j;
        }
    });
    (4.0 * FRAC_2_PI * se / norm, 2.0 * FRAC_2_PI * so / norm)
}

fn itth0(x: f64) -> f64 {
    // ===========================================================
    // Purpose: Evaluate the integral H0(t)/t with respect to t
    //          from x to infinity
    // Input :  x   --- Lower limit  ( x ≥ 0 )
    // Output:  TTH --- Integration of H0(t)/t from x to infinity
    // ===========================================================

    // Unlike in xsf (scipy/xsf#399), use the Neumann series for 0.5 <= x < 40. That leaves the
    // power series below to x < 0.5, and the asymptotic expansion to x >= 40.
    if (0.5..40.0).contains(&x) {
        return FRAC_PI_2 - itsh0_neumann(x).1;
    }

    let tth;
    let mut s = 1.0;
    let mut r = 1.0;
    if x < 24.5 {
        for k in 1..61 {
            let k = f64::from(k);
            r = -r * x * x * (2.0 * k - 1.0) / (2.0 * k + 1.0).powf(3.0);
            s += r;
            if r.abs() < s.abs() * 1.0e-12 {
                break;
            }
        }
        tth = FRAC_PI_2 - FRAC_2_PI * x * s;
    } else {
        // Unlike in xsf, sum up to eps (instead of 1e-12), and integrate `Y0(t)/t` by parts with
        // `(t Y1)' = t Y0` and `Y0' = -Y1`, instead of using fits with a relative error of 1e-8.
        for k in 1..11 {
            let k = f64::from(k);
            r = -r * (2.0 * k - 1.0).powf(3.0) / ((2.0 * k + 1.0) * x * x);
            s += r;
            if r.abs() < s.abs() * f64::EPSILON {
                break;
            }
        }
        let tts = 2.0 / (PI * x) * s;
        // int_x^inf Y0(t)/t dt = p Y0(x) - q Y1(x), up to eps or the smallest term
        let (mut p, mut q) = (0.0, 0.0);
        let mut c = 1.0 / x;
        for k in 1..31 {
            let k2 = f64::from(2 * k);
            p += c * k2 / x;
            q += c;
            if c.abs() < tts * f64::EPSILON || k2 >= x {
                break;
            }
            c *= -(k2 * k2) / (x * x);
        }
        tth = tts + p * cephes::y0(x) - q * cephes::y1(x);
    }
    tth
}

fn itsh0(x: f64) -> f64 {
    // ===================================================
    // Purpose: Evaluate the integral of Struve function
    //          H0(t) with respect to t from 0 and x
    // Input :  x   --- Upper limit  ( x ≥ 0 )
    // Output:  TH0 --- Integration of H0(t) from 0 and x
    // ===================================================

    let mut s;
    let th0;
    let el = 0.57721566490153;

    // Unlike in xsf (scipy/xsf#399), use the Neumann series for 0.5 <= x < 40. That leaves the
    // power series below to x < 0.5, and the asymptotic expansion to x >= 40.
    if (0.5..40.0).contains(&x) {
        return itsh0_neumann(x).0;
    }

    let mut r = 1.0;
    if x <= 30.0 {
        s = 0.5;
        for k in 1..101 {
            let rd = if k == 1 { 0.5 } else { 1.0 };
            let k = f64::from(k);
            let t = x / (2.0 * k + 1.0);
            r = -r * rd * k / (k + 1.0) * (t * t);
            s += r;
            if r.abs() < s.abs() * 1.0e-12 {
                break;
            }
        }
        th0 = FRAC_2_PI * x * x * s;
    } else {
        s = 1.0;
        for k in 1..13 {
            let k = f64::from(k);
            let t = (2.0 * k + 1.0) / x;
            r = -r * k / (k + 1.0) * (t * t);
            s += r;
            if r.abs() < s.abs() * 1.0e-12 {
                break;
            }
        }
        let s0 = s / (PI * x * x) + FRAC_2_PI * ((2.0 * x).ln() + el);
        th0 = itj0y0_asymp(x, 10).1 + s0;
    }
    th0
}

fn itsl0(x: f64) -> f64 {
    // ===========================================================
    // Purpose: Evaluate the integral of modified Struve function
    //          L0(t) with respect to t from 0 to x
    // Input :  x   --- Upper limit  ( x ≥ 0 )
    // Output:  TL0 --- Integration of L0(t) from 0 to x
    // ===========================================================

    let mut s;
    let tl0;
    let el = 0.57721566490153;
    let mut r = 1.0;
    if x <= 40.0 {
        s = 0.5;
        for k in 1..101 {
            let rd = if k == 1 { 0.5 } else { 1.0 };
            let k = f64::from(k);
            let t = x / (2.0 * k + 1.0);
            r = r * rd * k / (k + 1.0) * (t * t);
            s += r;
            if (r / s).abs() < f64::EPSILON {
                break;
            }
        }
        tl0 = FRAC_2_PI * x * x * s;
    } else {
        s = 1.0;
        for k in 1..11 {
            let k = f64::from(k);
            let t = (2.0 * k + 1.0) / x;
            r = r * k / (k + 1.0) * (t * t);
            s += r;
            if (r / s).abs() < 1.0e-12 {
                break;
            }
        }
        let s0 = -s / (PI * x * x) + FRAC_2_PI * ((2.0 * x).ln() + el);
        let mut ti = 1.0;
        r = 1.0;
        // all 25 terms for x > 40, where the last one is at most 2.2e-16
        for ak in IT0_ASYMP {
            r /= x;
            ti += ak * r;
        }
        tl0 = ti / (TAU * x).sqrt() * x.exp() + s0;
    }
    tl0
}

/// Integral of the Struve $H_0$ function
///
/// $$ \int_0^x H_0(t) \dd t $$
///
/// Corresponds to [`scipy.special.itstruve0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itstruve0.html
///
/// # See also
/// - [`struve_h`]: Struve $H_v$ function
/// - [`it2struve0`]
#[doc(alias = "itstruve_h0", alias = "it1struve0")]
#[must_use]
#[inline]
pub fn itstruve0(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }

    specfun_convinf(itsh0(x))
}

/// Integral related to the Struve $H_0$ function
///
/// $$ \int_x^\infty {H_0(t) \over t} \dd t $$
///
/// Corresponds to [`scipy.special.it2struve0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2struve0.html
///
/// # See also
/// - [`struve_h`]: Struve $H_v$ function
/// - [`itstruve0`]
#[doc(alias = "it2struve_l0")]
#[must_use]
#[inline]
pub fn it2struve0(mut x: f64) -> f64 {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }

    let mut out = specfun_convinf(itth0(x));
    if flag {
        out = PI - out;
    }
    out
}

/// Integral of the modified Struve $L_0$ function
///
/// $$ \int_0^x L_0(t) \dd t $$
///
/// Corresponds to [`scipy.special.itmodstruve0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itmodstruve0.html
///
/// # See also
/// - [`struve_l`]: Struve $L_v$ function
#[doc(alias = "itstruve_l0")]
#[must_use]
#[inline]
pub fn itmodstruve0(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }

    specfun_convinf(itsl0(x))
}

/// Struve $H_v$ function
///
/// Corresponds to [`scipy.special.struve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.struve.html
///
/// # Definition
///
/// $$
/// H_v(x) =
/// (x/2)^{v+1}
/// \sum_{n=0}^\infty {(-1)^n (x/2)^{2n} \over \Gamma(n+3/2) \\ \Gamma(n+v+3/2)} \\ ,
/// $$
///
/// where $\Gamma$ is the Gamma function.
///
/// # See also
/// - [`struve_l`]: Struve $L_v$ function
/// - [`itstruve0`]: Integral of $H_0$
/// - [`it2struve0`]: Integral related to $H_0$
#[doc(alias = "struve")]
#[must_use]
#[inline]
pub fn struve_h(v: f64, x: f64) -> f64 {
    crate::xsf::cephes::struve_h(v, x)
}

/// Modified Struve $L_v$ function
///
/// Corresponds to [`scipy.special.modstruve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modstruve.html
///
/// # Definition
///
/// $$
/// \begin{align*}
/// L_v(x)
/// &= (x/2)^{v+1} \sum_{n=0}^\infty {(x/2)^{2n} \over \Gamma(n+3/2) \\ \Gamma(n+v+3/2)} \\\\
/// &= -i e^{-\pi iv/2} H_v(ix) \\ ,
/// \end{align*}
/// $$
///
/// where $\Gamma$ is the Gamma function, and $H_v$ the Struve function of order $v$ ([`struve_h`]).
///
/// # See also
/// - [`struve_h`]: Struve $H_v$ function
/// - [`itmodstruve0`]: Integral of $L_0$
#[doc(alias = "modstruve")]
#[must_use]
#[inline]
pub fn struve_l(v: f64, x: f64) -> f64 {
    crate::xsf::cephes::struve_l(v, x)
}
