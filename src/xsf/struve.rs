//! `itstruve0`, `it2struve0`, `itmodstruve0`, `struve_h`, and `struve_l` are translated into pure
//! Rust from `xsf/struve.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM only does
//! when optimizing.

#![allow(clippy::unreadable_literal)]

use core::f64::consts::{FRAC_2_PI, FRAC_PI_2, PI};

use crate::xsf::cephes;

#[allow(clippy::float_cmp)]
fn specfun_convinf(x: f64) -> f64 {
    if x == 1.0e300 {
        return f64::INFINITY;
    }
    if x == -1.0e300 {
        return f64::NEG_INFINITY;
    }
    x
}

/// `(int_0^x H0(t) dt, int_0^x H0(t)/t dt)` for `x >= 0.5`, by integrating
/// `H0 = 4/pi sum_k J_{2k+1} / (2k+1)` (DLMF 11.4.21) term by term:
///
/// - `int_0^x J_{2k+1}(t) dt = 2 sum_{j>k} J_{2j}(x)`
/// - `int_0^x J_{2k+1}(t)/t dt = (2 sum_{j>=k} J_{2j+1}(x) - J_{2k+1}(x)) / (2k+1)`
///
/// with the `J_n(x)` from Miller's backward recurrence, normalized by `J_0 + 2 sum_j J_{2j} = 1`.
fn itsh0_neumann(x: f64) -> (f64, f64) {
    let mut n = 2.0 * (0.5 * x).floor() + 51.0;
    // unnormalized J_{n+1}(x) and J_n(x)
    let (mut jp, mut j) = (0.0, 1.0);
    // sum_{2j > n} J_{2j}, sum_{2j+1 >= n} J_{2j+1}, and the two series
    let (mut te, mut to, mut se, mut so) = (0.0, 0.0, 0.0, 0.0);
    loop {
        to += j;
        so += (2.0 * to - j) / (n * n);
        se += te / n;
        (jp, j) = (j, 2.0 * n / x * j - jp);
        if n < 2.0 {
            break;
        }
        te += j;
        (jp, j) = (j, 2.0 * (n - 1.0) / x * j - jp);
        n -= 2.0;
    }
    let norm = j + 2.0 * te;
    (4.0 * FRAC_2_PI * se / norm, 2.0 * FRAC_2_PI * so / norm)
}

#[allow(clippy::needless_late_init)]
fn itth0(x: f64) -> f64 {
    // ===========================================================
    // Purpose: Evaluate the integral H0(t)/t with respect to t
    //          from x to infinity
    // Input :  x   --- Lower limit  ( x ≥ 0 )
    // Output:  TTH --- Integration of H0(t)/t from x to infinity
    // ===========================================================

    // Unlike in xsf, the Neumann series for 0.5 <= x < 40, where the series cancels and the
    // asymptotic expansion is inaccurate, e.g. `it2struve0(24.450718139337205)` was 0.0262689004
    // instead of 0.0262688934.
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
        tth = PI / 2.0 - FRAC_2_PI * x * s;
    } else {
        // Unlike in xsf, sum up to eps (instead of 1e-12), and integrate `Y0(t)/t` by parts with
        // `(t Y1)' = t Y0` and `Y0' = -Y1`, instead of using fits with a relative error of 1e-8.
        for k in 1..31 {
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

#[allow(clippy::manual_midpoint, clippy::needless_late_init)]
fn itsh0(x: f64) -> f64 {
    // ===================================================
    // Purpose: Evaluate the integral of Struve function
    //          H0(t) with respect to t from 0 and x
    // Input :  x   --- Upper limit  ( x ≥ 0 )
    // Output:  TH0 --- Integration of H0(t) from 0 and x
    // ===================================================

    let mut a = [0.0; 25];
    let mut s;
    let th0;
    let el = 0.57721566490153;

    // Unlike in xsf, the Neumann series for 0.5 <= x < 40, where the series cancels, e.g.
    // `itstruve0(29.746969855043805)` was 3.08382 instead of 3.08387.
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
        let mut a0 = 1.0;
        let mut a1 = 5.0 / 8.0;
        a[0] = a1;
        for k in 1_u32..21 {
            let kf = f64::from(k);
            let af = (1.5 * (kf + 0.5) * (kf + 5.0 / 6.0) * a1
                - 0.5 * (kf + 0.5) * (kf + 0.5) * (kf - 0.5) * a0)
                / (kf + 1.0);
            a[k as usize] = af;
            a0 = a1;
            a1 = af;
        }
        let mut bf = 1.0;
        r = 1.0;
        for k in 1..11 {
            r = -r / (x * x);
            bf += a[2 * k - 1] * r;
        }
        let mut bg = a[0] / x;
        r = 1.0 / x;
        for k in 1..=10 {
            r = -r / (x * x);
            bg += a[2 * k] * r;
        }
        let xp = x + 0.25 * PI;
        let ty = (2.0 / (PI * x)).sqrt() * (bg * xp.cos() - bf * xp.sin());
        th0 = ty + s0;
    }
    th0
}

#[allow(clippy::needless_late_init)]
fn itsl0(x: f64) -> f64 {
    // ===========================================================
    // Purpose: Evaluate the integral of modified Struve function
    //          L0(t) with respect to t from 0 to x
    // Input :  x   --- Upper limit  ( x ≥ 0 )
    // Output:  TL0 --- Integration of L0(t) from 0 to x
    // ===========================================================

    let mut a = [0.0; 25];
    let mut s;
    let tl0;
    let el = 0.57721566490153;
    let mut r = 1.0;
    // Unlike in xsf, sum the series up to eps and use it up to x = 40 (instead of 20), beyond
    // which 25 (instead of 10) asymptotic terms reach eps, e.g.
    // `itmodstruve0(21.281112146816426)` was 155893841.987 instead of 155893845.553.
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
        let mut a0 = 1.0;
        let mut a1 = 5.0 / 8.0;
        a[0] = a1;
        for k in 1_u32..25 {
            let kf = f64::from(k);
            let af = (1.5 * (kf + 0.50) * (kf + 5.0 / 6.0) * a1
                - 0.5 * ((kf + 0.5) * (kf + 0.5)) * (kf - 0.5) * a0)
                / (kf + 1.0);
            a[k as usize] = af;
            a0 = a1;
            a1 = af;
        }
        let mut ti = 1.0;
        r = 1.0;
        for k in 1..26 {
            r /= x;
            ti += a[k - 1] * r;
        }
        tl0 = ti / (2.0 * PI * x).sqrt() * x.exp() + s0;
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
/// &= (z/2)^{v+1} \sum_{n=0}^\infty {(-z/2)^{2n} \over \Gamma(n+3/2) \\ \Gamma(n+v+3/2)} \\\\
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
