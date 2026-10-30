//! `itstruve0`, `it2struve0`, `itmodstruve0`, `struve_h`, and `struve_l` are translated into pure
//! Rust from `xsf/struve.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM only does
//! when optimizing.

#![allow(clippy::unreadable_literal)]

use core::f64::consts::PI;

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

#[allow(clippy::needless_late_init)]
fn itth0(x: f64) -> f64 {
    // ===========================================================
    // Purpose: Evaluate the integral H0(t)/t with respect to t
    //          from x to infinity
    // Input :  x   --- Lower limit  ( x ≥ 0 )
    // Output:  TTH --- Integration of H0(t)/t from x to infinity
    // ===========================================================

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
        tth = PI / 2.0 - 2.0 / PI * x * s;
    } else {
        for k in 1..11 {
            let k = f64::from(k);
            r = -r * (2.0 * k - 1.0).powf(3.0) / ((2.0 * k + 1.0) * x * x);
            s += r;
            if r.abs() < s.abs() * 1.0e-12 {
                break;
            }
        }
        let t = 8.0 / x;
        let xt = x + 0.25 * PI;
        let f0 =
            (((((0.18118e-2 * t - 0.91909e-2) * t + 0.017033) * t - 0.9394e-3) * t - 0.051445) * t
                - 0.11e-5)
                * t
                + 0.7978846;
        let g0 = (((((-0.23731e-2 * t + 0.59842e-2) * t + 0.24437e-2) * t - 0.0233178) * t
            + 0.595e-4)
            * t
            + 0.1620695)
            * t;
        let tty = (f0 * xt.sin() - g0 * xt.cos()) / (x.sqrt() * x);
        tth = 2.0 / (PI * x) * s + tty;
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
        th0 = 2.0 / PI * x * x * s;
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
        let s0 = s / (PI * x * x) + 2.0 / PI * ((2.0 * x).ln() + el);
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
        let mut bg = a[0] * x;
        r = 1.0 / x;
        for k in 1..10 {
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

    let mut a = [0.0; 18];
    let mut s;
    let tl0;
    let el = 0.57721566490153;
    let mut r = 1.0;
    if x <= 20.0 {
        s = 0.5;
        for k in 1..101 {
            let rd = if k == 1 { 0.5 } else { 1.0 };
            let k = f64::from(k);
            let t = x / (2.0 * k + 1.0);
            r = r * rd * k / (k + 1.0) * (t * t);
            s += r;
            if (r / s).abs() < 1.0e-12 {
                break;
            }
        }
        tl0 = 2.0 / PI * x * x * s;
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
        let s0 = -s / (PI * x * x) + 2.0 / PI * ((2.0 * x).ln() + el);
        let mut a0 = 1.0;
        let mut a1 = 5.0 / 8.0;
        a[0] = a1;
        for k in 1_u32..11 {
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
        for k in 1..11 {
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

#[cfg(test)]
mod tests {
    #[test]
    fn test_itstruve0_f64() {
        xsref::test("itstruve0", "d-d", |x| crate::itstruve0(x[0]));
    }

    #[test]
    fn test_it2struve0_f64() {
        xsref::test("it2struve0", "d-d", |x| crate::it2struve0(x[0]));
    }

    #[test]
    fn test_itmodstruve0_f64() {
        xsref::test("itmodstruve0", "d-d", |x| crate::itmodstruve0(x[0]));
    }

    #[test]
    fn test_struve_h_f64() {
        xsref::test("struve_h", "d_d-d", |x| crate::struve_h(x[0], x[1]));
    }

    #[test]
    fn test_struve_l_f64() {
        xsref::test("struve_l", "d_d-d", |x| crate::struve_l(x[0], x[1]));
    }
}
