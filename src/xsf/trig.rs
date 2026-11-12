//! Translated into pure Rust from xsf (scipy/xsf@d350b94):
//!
//! - `xsf/trig.h`: original author: Josh Wilson, 2016.
//! - `xsf/cephes/trig.h`: original author: Josh Wilson, 2020.
//! - `xsf/cephes/sindg.h`, `xsf/cephes/tandg.h`, and `xsf/cephes/unity.h`: from the Cephes Math
//!   Library Release 2.0 (April, 1987), Copyright 1984, 1987 by Stephen L. Moshier.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::utils::{mul_cosh_sinh, sin_cos_glibc};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::two_sum;
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

/// sin(pi (x - m)), without the rounding error of x - m near an integer
pub(crate) fn sinpi_sub(x: f64, m: f64) -> f64 {
    let (s, e) = two_sum(x, -m);
    let n = s.round();
    let t = cephes_sinpi((s - n) + e);
    if n % 2.0 == 0.0 { t } else { -t }
}

/// Sine of complex `z`
///
/// This replaces `std::sin(std::complex<double>)` from C++. For finite `z`, it is computed like
/// glibc's `csin`, which, unlike `num_complex::Complex::sin`, doesn't overflow prematurely (to
/// infinity or NaN) when only `cosh(z.im)` or `sinh(z.im)` overflows.
pub(crate) fn csin(z: C64) -> C64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        return Complex::from(z).sin().into();
    }

    let (mut sinix, cosix) = sin_cos_glibc(z.re.abs());
    if z.re.is_sign_negative() {
        sinix = -sinix;
    }
    mul_cosh_sinh(z.im, sinix, cosix)
}

/// Cosine of complex `z`
///
/// This replaces `std::cos(std::complex<double>)` from C++. For finite `z`, it is computed like
/// glibc's `ccos`, i.e. as `ccosh(i z)`, which, unlike `num_complex::Complex::cos`, doesn't
/// overflow prematurely (to infinity or NaN) when only `cosh(z.im)` or `sinh(z.im)` overflows.
pub(crate) fn ccos(z: C64) -> C64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        return Complex::from(z).cos().into();
    }

    // ccosh(x + i y) with x = -z.im and y = z.re
    let (sinix, cosix) = sin_cos_glibc(z.re);
    mul_cosh_sinh(-z.im, cosix, sinix)
}

/* Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
 * of these functions are integral (and thus better representable in
 * floating point), it's possible to compute them with greater accuracy
 * than sin(z), cos(z).
 */

/// Hyperbolic sine of complex `z`
///
/// This replaces `std::sinh(std::complex<double>)` from C++, as `sinh(z) = -i sin(i z)`.
pub(crate) fn csinh(z: C64) -> C64 {
    let w = csin(c64(-z.im, z.re));
    c64(w.im, -w.re)
}

/// Hyperbolic cosine of complex `z`
///
/// This replaces `std::cosh(std::complex<double>)` from C++, as `cosh(z) = cos(-i z)`.
pub(crate) fn ccosh(z: C64) -> C64 {
    ccos(c64(z.im, -z.re))
}

/// `xsf::sinpi` for complex `z`
pub(crate) fn csinpi(z: C64) -> C64 {
    let x = z.re;
    let piy = PI * z.im;
    let sinpix = cephes_sinpi(x);
    let cospix = cephes_cospi(x);

    // Unlike in xsf, `cosh` and `sinh` of large |pi y| are scaled as in glibc's `csin`, so that
    // this overflows only if the result does.
    mul_cosh_sinh(piy, sinpix, cospix)
}

/// `xsf::cospi` for complex `z`
pub(crate) fn ccospi(z: C64) -> C64 {
    let x = z.re;
    let piy = PI * z.im;
    let sinpix = cephes_sinpi(x);
    let cospix = cephes_cospi(x);

    // Unlike in xsf, this is scaled as in glibc's `ccos`; see csinpi(z).
    mul_cosh_sinh(-piy, cospix, sinpix)
}

/// `cephes::detail::tancot`
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
            // Unlike in xsf (scipy/xsf#363), cot(x) = 1 / tan(x) is used for 0 < x < 45
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
    // Unlike in xsf (scipy/xsf#363), tan(x) = 1 / tan(90 - x) is used for x > 45
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
        self.via(csinpi)
    }

    #[inline]
    fn cospi(self) -> Self {
        self.via(ccospi)
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
