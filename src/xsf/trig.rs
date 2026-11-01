//! Translated into pure Rust from xsf v0.2.2:
//!
//! - `xsf/trig.h`: translated from Cython into C++ by SciPy developers in 2023, original author:
//!   Josh Wilson, 2016.
//! - `xsf/cephes/trig.h`: translated into C++ by SciPy developers in 2024, original author: Josh
//!   Wilson, 2020.
//! - `xsf/cephes/sindg.h`, `xsf/cephes/tandg.h`, and `xsf/cephes/unity.h`: translated into C++ by
//!   SciPy developers in 2024 from the Cephes Math Library Release 2.0 (April, 1987),
//!   Copyright 1984, 1987 by Stephen L. Moshier.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::cephes::polevl::polevl;

use core::f64::consts::{FRAC_PI_4, PI};
use num_complex::Complex;

// cephes/const.h

/// pi/180
const PI180: f64 = 1.745_329_251_994_329_5e-2;

// cephes/sindg.h

const SINCOF: [f64; 6] = [
    1.589_623_015_722_184_4e-10,
    -2.505_074_776_285_035_5e-8,
    2.755_731_362_138_567_6e-6,
    -1.984_126_982_958_954e-4,
    8.333_333_333_322_118e-3,
    -1.666_666_666_666_663e-1,
];

const COSCOF: [f64; 7] = [
    1.136_781_713_820_445_5e-11,
    -2.087_588_337_576_836_3e-9,
    2.755_731_554_298_166_3e-7,
    -2.480_158_729_361_863e-5,
    1.388_888_888_888_066_7e-3,
    -4.166_666_666_666_664e-2,
    5.0e-1,
];

const SINDG_LOSSTH: f64 = 1.0e14;

/* 1 arc second, in radians = 4.848136811095359935899141023579479759563533023727e-6 */
const SINDG_P64800: f64 = 4.848_136_811_095_36e-6;

// cephes/tandg.h

const TANDG_LOSSTH: f64 = 1.0e14;

// cephes/unity.h

const UNITY_COSCOF: [f64; 7] = [
    4.737_750_796_424_621e-14,
    -1.147_028_484_342_536e-11,
    2.087_675_428_708_152e-9,
    -2.755_731_921_499_979e-7,
    2.480_158_730_157_055e-5,
    -1.388_888_888_888_887_2e-3,
    4.166_666_666_666_666_4e-2,
];

/*
 * Implement sin(pi * x) and cos(pi * x) for real x. Since the periods
 * of these functions are integral (and thus representable in double
 * precision), it's possible to compute them with greater accuracy
 * than sin(x) and cos(x).
 */

/// `cephes::sinpi`: compute sin(pi * x).
pub(crate) fn cephes_sinpi(mut x: f64) -> f64 {
    let mut s = 1.0;

    if x.is_sign_negative() {
        x = -x;
        s = -1.0;
    }

    let r = x % 2.0; // std::fmod(x, 2.0)
    if r < 0.5 {
        s * (PI * r).sin()
    } else if r > 1.5 {
        s * (PI * (r - 2.0)).sin()
    } else {
        -s * (PI * (r - 1.0)).sin()
    }
}

/// `cephes::cospi`: compute cos(pi * x)
#[allow(clippy::float_cmp)]
pub(crate) fn cephes_cospi(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }

    let r = x % 2.0; // std::fmod(x, 2.0)
    if r == 0.5 {
        // We don't want to return -0.0
        return 0.0;
    }
    if r < 1.0 {
        -(PI * (r - 0.5)).sin()
    } else {
        (PI * (r - 1.5)).sin()
    }
}

/// Sine of complex `z`
///
/// This replaces `std::sin(std::complex<double>)` from C++. For finite `z`, it is computed like
/// glibc's `csin`, which, unlike `num_complex::Complex::sin`, doesn't overflow prematurely (to
/// infinity or NaN) when only `cosh(z.im)` or `sinh(z.im)` overflows.
pub(crate) fn csin(z: Complex<f64>) -> Complex<f64> {
    // (int) ((DBL_MAX_EXP - 1) * ln(2))
    const T: f64 = 709.0;

    if !z.re.is_finite() || !z.im.is_finite() {
        return z.sin();
    }

    let rx = z.re.abs();
    let (mut sinix, mut cosix) = if rx > f64::MIN_POSITIVE {
        rx.sin_cos()
    } else {
        (rx, 1.0)
    };
    if z.re.is_sign_negative() {
        sinix = -sinix;
    }
    if z.im.abs() > T {
        let exp_t = T.exp();
        let mut ix = z.im.abs();
        if z.im.is_sign_negative() {
            cosix = -cosix;
        }
        ix -= T;
        sinix *= exp_t / 2.0;
        cosix *= exp_t / 2.0;
        if ix > T {
            ix -= T;
            sinix *= exp_t;
            cosix *= exp_t;
        }
        if ix > T {
            // overflow (the original |z.im| > 3 T)
            Complex::new(f64::MAX * sinix, f64::MAX * cosix)
        } else {
            let ev = ix.exp();
            Complex::new(ev * sinix, ev * cosix)
        }
    } else {
        Complex::new(z.im.cosh() * sinix, z.im.sinh() * cosix)
    }
}

/// Cosine of complex `z`
///
/// This replaces `std::cos(std::complex<double>)` from C++. For finite `z`, it is computed like
/// glibc's `ccos`, i.e. as `ccosh(i z)`, which, unlike `num_complex::Complex::cos`, doesn't
/// overflow prematurely (to infinity or NaN) when only `cosh(z.im)` or `sinh(z.im)` overflows.
pub(crate) fn ccos(z: Complex<f64>) -> Complex<f64> {
    // (int) ((DBL_MAX_EXP - 1) * ln(2))
    const T: f64 = 709.0;

    if !z.re.is_finite() || !z.im.is_finite() {
        return z.cos();
    }

    // ccosh(x + i y) with x = -z.im and y = z.re
    let (x, y) = (-z.im, z.re);
    let (mut sinix, mut cosix) = if y.abs() > f64::MIN_POSITIVE {
        y.sin_cos()
    } else {
        (y, 1.0)
    };
    if x.abs() > T {
        let exp_t = T.exp();
        let mut rx = x.abs();
        if x.is_sign_negative() {
            sinix = -sinix;
        }
        rx -= T;
        sinix *= exp_t / 2.0;
        cosix *= exp_t / 2.0;
        if rx > T {
            rx -= T;
            sinix *= exp_t;
            cosix *= exp_t;
        }
        if rx > T {
            // overflow (the original |z.im| > 3 T)
            Complex::new(f64::MAX * cosix, f64::MAX * sinix)
        } else {
            let ev = rx.exp();
            Complex::new(ev * cosix, ev * sinix)
        }
    } else {
        Complex::new(x.cosh() * cosix, x.sinh() * sinix)
    }
}

/* Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
 * of these functions are integral (and thus better representable in
 * floating point), it's possible to compute them with greater accuracy
 * than sin(z), cos(z).
 */

/// Hyperbolic sine of complex `z`
///
/// This replaces `std::sinh(std::complex<double>)` from C++, as `sinh(z) = -i sin(i z)`.
pub(crate) fn csinh(z: Complex<f64>) -> Complex<f64> {
    let w = csin(Complex::new(-z.im, z.re));
    Complex::new(w.im, -w.re)
}

/// Hyperbolic cosine of complex `z`
///
/// This replaces `std::cosh(std::complex<double>)` from C++, as `cosh(z) = cos(-i z)`.
pub(crate) fn ccosh(z: Complex<f64>) -> Complex<f64> {
    ccos(Complex::new(z.im, -z.re))
}

/// `xsf::sinpi` for complex `z`
#[allow(clippy::float_cmp)]
fn csinpi(z: Complex<f64>) -> Complex<f64> {
    let x = z.re;
    let piy = PI * z.im;
    let abspiy = piy.abs();
    let sinpix = cephes_sinpi(x);
    let cospix = cephes_cospi(x);

    if abspiy < 700.0 {
        return Complex::new(sinpix * piy.cosh(), cospix * piy.sinh());
    }

    /* Have to be careful--sinh/cosh could overflow while cos/sin are small.
     * At this large of values
     *
     * cosh(y) ~ exp(y)/2
     * sinh(y) ~ sgn(y)*exp(y)/2
     *
     * so we can compute exp(y/2), scale by the right factor of sin/cos
     * and then multiply by exp(y/2) to avoid overflow. */
    // NOTE(xsf-rust): unlike xsf v0.2.2 (scipy/xsf#308), this accounts for sgn(y).
    //
    // Unlike xsf main, this also propagates NaN, and avoids premature rounding of tiny sin/cos
    // factors (scipy/xsf#308 (comment)).
    //
    // It still overflows prematurely if exp(|pi y|/2) overflows while the result wouldn't (only
    // for |x| < ~1e-308).
    let cospix = 1.0_f64.copysign(piy) * cospix;
    let exphpiy = (abspiy / 2.0).exp();
    let coshfac;
    let sinhfac;
    if exphpiy == f64::INFINITY {
        if sinpix == 0.0 {
            // Preserve the sign of zero.
            coshfac = 0.0_f64.copysign(sinpix);
        } else {
            coshfac = sinpix * f64::INFINITY;
        }
        if cospix == 0.0 {
            // Preserve the sign of zero.
            sinhfac = 0.0_f64.copysign(cospix);
        } else {
            sinhfac = cospix * f64::INFINITY;
        }
        return Complex::new(coshfac, sinhfac);
    }

    coshfac = sinpix * exphpiy * 0.5;
    sinhfac = cospix * exphpiy * 0.5;
    Complex::new(coshfac * exphpiy, sinhfac * exphpiy)
}

/// `xsf::cospi` for complex `z`
#[allow(clippy::float_cmp)]
fn ccospi(z: Complex<f64>) -> Complex<f64> {
    let x = z.re;
    let piy = PI * z.im;
    let abspiy = piy.abs();
    let sinpix = cephes_sinpi(x);
    let cospix = cephes_cospi(x);

    if abspiy < 700.0 {
        return Complex::new(cospix * piy.cosh(), -sinpix * piy.sinh());
    }

    // See csinpi(z) for an idea of what's going on here.
    // NOTE(xsf-rust): unlike xsf v0.2.2 (scipy/xsf#308), this uses the correct sign of the
    // imaginary part, accounts for sgn(y), and checks the right factor for zero.
    //
    // Unlike xsf main, this also propagates NaN, and avoids premature rounding of tiny sin/cos
    // factors (scipy/xsf#308 (comment), and see csinpi(z) for the remaining limitation).
    let sinpix = -1.0_f64.copysign(piy) * sinpix;
    let exphpiy = (abspiy / 2.0).exp();
    let coshfac;
    let sinhfac;
    if exphpiy == f64::INFINITY {
        if cospix == 0.0 {
            // Preserve the sign of zero.
            coshfac = 0.0_f64.copysign(cospix);
        } else {
            coshfac = cospix * f64::INFINITY;
        }
        if sinpix == 0.0 {
            // Preserve the sign of zero.
            sinhfac = 0.0_f64.copysign(sinpix);
        } else {
            sinhfac = sinpix * f64::INFINITY;
        }
        return Complex::new(coshfac, sinhfac);
    }

    coshfac = cospix * exphpiy * 0.5;
    sinhfac = sinpix * exphpiy * 0.5;
    Complex::new(coshfac * exphpiy, sinhfac * exphpiy)
}

/// `cephes::detail::tancot`
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
fn tancot(xx: f64, cotflg: bool) -> f64 {
    let mut x;
    let mut sign: i32;

    /* make argument positive but save the sign */
    if xx.is_sign_negative() {
        x = -xx;
        sign = -1;
    } else {
        x = xx;
        sign = 1;
    }

    if x > TANDG_LOSSTH {
        return 0.0;
    }

    /* modulo 180 */
    let k = (x / 180.0).floor();
    x -= 180.0 * k;
    if cotflg {
        if x <= 90.0 {
            // Unlike in xsf (scipy/xsf#363), cot(x) = 1 / tan(x) is used for 0 < x < 45, since
            // 90 - x rounds, which is amplified near the pole at x = 0.
            //
            // E.g. cotdg(1e-10) was 572956949790.9288 instead of 572957795130.8232, and
            // cotdg(1e-15) was inf.
            if x != 0.0 && x < 45.0 {
                return f64::from(sign) / (x * PI180).tan();
            }
            x = 90.0 - x;
        } else {
            x -= 90.0;
            sign *= -1;
        }
    } else if x > 90.0 {
        x = 180.0 - x;
        sign *= -1;
    }
    if x == 0.0 {
        if (k as i64) & 1 != 0 {
            sign *= -1;
        }
        return f64::from(sign) * 0.0;
    } else if x == 45.0 {
        return f64::from(sign) * 1.0;
    } else if x == 90.0 {
        if (k as i64) & 1 != 0 {
            sign *= -1;
        }
        return f64::from(sign) * f64::INFINITY;
    }
    /* x is now transformed into [0, 90) */
    // Unlike in xsf (scipy/xsf#363), tan(x) = 1 / tan(90 - x) is used for x > 45, since x * PI180
    // rounds, which is amplified near the pole at x = 90 (whereas 90 - x is exact).
    //
    // E.g. tandg(89.9999999999) was 572956949790.9288 instead of 572947570131.5804.
    if x > 45.0 {
        return f64::from(sign) / ((90.0 - x) * PI180).tan();
    }
    f64::from(sign) * (x * PI180).tan()
}

pub trait TrigArg: crate::sealed::Sealed {
    fn sinpi(self) -> Self;
    fn cospi(self) -> Self;
}

impl TrigArg for f64 {
    #[inline]
    fn sinpi(self) -> Self {
        cephes_sinpi(self)
    }

    #[inline]
    fn cospi(self) -> Self {
        cephes_cospi(self)
    }
}

impl TrigArg for Complex<f64> {
    #[inline]
    fn sinpi(self) -> Self {
        csinpi(self)
    }

    #[inline]
    fn cospi(self) -> Self {
        ccospi(self)
    }
}

/// $\sin(\pi z)$ for real or complex $z$
///
/// Has no corresponding function in `scipy.special`.
///
/// # See also
/// - [`cospi`]
#[must_use]
#[inline]
pub fn sinpi<T: TrigArg>(z: T) -> T {
    z.sinpi()
}

/// $\cos(\pi z)$ for real or complex $z$
///
/// Has no corresponding function in `scipy.special`.
///
/// # See also
/// - [`sinpi`]
#[must_use]
#[inline]
pub fn cospi<T: TrigArg>(z: T) -> T {
    z.cospi()
}

/// Sine of the angle x given in degrees.
///
/// Corresponds to [`scipy.special.sindg`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sindg.html
///
/// # See also
/// - [`cosdg`]
/// - [`tandg`]
/// - [`cotdg`]
#[must_use]
#[inline]
#[allow(clippy::cast_possible_truncation)]
pub fn sindg(mut x: f64) -> f64 {
    let mut y;
    let mut z;
    let mut j: i32;
    let mut sign: i32;

    /* make argument positive but save the sign */
    sign = 1;
    if x.is_sign_negative() {
        x = -x;
        sign = -1;
    }

    if x > SINDG_LOSSTH {
        return 0.0;
    }

    y = (x / 45.0).floor(); /* integer part of x/M_PI_4 */

    /* strip high bits of integer part to prevent integer overflow */
    z = y * 0.0625; // std::ldexp(y, -4), and `z * 16.0` below is std::ldexp(z, 4)
    z = z.floor(); /* integer part of y/8 */
    z = y - z * 16.0; /* y - 16 * (y/16) */

    j = z as i32; /* convert to integer for tests on the phase angle */
    /* map zeros to origin */
    if j & 1 != 0 {
        j += 1;
        y += 1.0;
    }
    j &= 0o7; /* octant modulo 360 degrees */
    /* reflect in x axis */
    if j > 3 {
        sign = -sign;
        j -= 4;
    }

    z = x - y * 45.0; /* x mod 45 degrees */
    z *= PI180; /* multiply by pi/180 to convert to radians */
    let zz = z * z;

    if j == 1 || j == 2 {
        y = 1.0 - zz * polevl(zz, &COSCOF);
    } else {
        y = z + z * (zz * polevl(zz, &SINCOF));
    }

    if sign < 0 {
        y = -y;
    }

    y
}

/// Cosine of angle in degrees
///
/// Corresponds to [`scipy.special.cosdg`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cosdg.html
///
/// # See also
/// - [`sindg`]
/// - [`tandg`]
/// - [`cotdg`]
#[must_use]
#[inline]
#[allow(clippy::cast_possible_truncation)]
pub fn cosdg(mut x: f64) -> f64 {
    let mut y;
    let mut z;
    let mut j: i32;
    let mut sign: i32;

    /* make argument positive */
    sign = 1;
    if x < 0.0 {
        x = -x;
    }

    if x > SINDG_LOSSTH {
        return 0.0;
    }

    y = (x / 45.0).floor();
    z = y * 0.0625; // std::ldexp(y, -4), and `z * 16.0` below is std::ldexp(z, 4)
    z = z.floor(); /* integer part of y/8 */
    z = y - z * 16.0; /* y - 16 * (y/16) */

    /* integer and fractional part modulo one octant */
    j = z as i32;
    if j & 1 != 0 {
        /* map zeros to origin */
        j += 1;
        y += 1.0;
    }
    j &= 0o7;
    if j > 3 {
        j -= 4;
        sign = -sign;
    }

    if j > 1 {
        sign = -sign;
    }

    z = x - y * 45.0; /* x mod 45 degrees */
    z *= PI180; /* multiply by pi/180 to convert to radians */

    let zz = z * z;

    if j == 1 || j == 2 {
        y = z + z * (zz * polevl(zz, &SINCOF));
    } else {
        y = 1.0 - zz * polevl(zz, &COSCOF);
    }

    if sign < 0 {
        y = -y;
    }

    y
}

/// Tangent of angle x given in degrees
///
/// Corresponds to [`scipy.special.tandg`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.tandg.html
///
/// # See also
/// - [`sindg`]
/// - [`cosdg`]
/// - [`cotdg`]
#[must_use]
#[inline]
pub fn tandg(x: f64) -> f64 {
    tancot(x, false)
}

/// Cotangent of the angle x given in degrees
///
/// Corresponds to [`scipy.special.cotdg`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cotdg.html
///
/// # See also
/// - [`sindg`]
/// - [`cosdg`]
/// - [`tandg`]
#[must_use]
#[inline]
pub fn cotdg(x: f64) -> f64 {
    tancot(x, true)
}

/// $\cos(x) - 1$ for use when $x$ is near zero
///
/// Corresponds to [`scipy.special.cosm1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cosm1.html
///
/// # See also
/// - [`expm1`](crate::expm1)
/// - [`log1p`](crate::log1p)
#[must_use]
#[inline]
pub fn cosm1(x: f64) -> f64 {
    #[allow(clippy::manual_range_contains)]
    if x < -FRAC_PI_4 || x > FRAC_PI_4 {
        return x.cos() - 1.0;
    }
    let mut xx = x * x;
    xx = -0.5 * xx + xx * xx * polevl(xx, &UNITY_COSCOF);
    xx
}

/// Convert from degrees to radians.
///
/// Returns the angle given in (d)egrees, (m)inutes, and (s)econds in radians.
///
/// Corresponds to [`scipy.special.radian`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.radian.html
#[must_use]
#[inline]
pub fn radian(d: f64, m: f64, s: f64) -> f64 {
    /* Degrees, minutes, seconds to radians: */
    ((d * 60.0 + m) * 60.0 + s) * SINDG_P64800
}
