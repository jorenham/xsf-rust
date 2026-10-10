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

    if x < 0.0 {
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
fn cephes_cospi(mut x: f64) -> f64 {
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

/* Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
 * of these functions are integral (and thus better representable in
 * floating point), it's possible to compute them with greater accuracy
 * than sin(z), cos(z).
 */

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
    // NOTE(xsf-rust): unlike xsf v0.2.2, this accounts for sgn(y), propagates NaN, and avoids
    // premature rounding of tiny sin/cos factors. It still overflows prematurely if
    // exp(|pi y|/2) overflows while the result wouldn't (only for |x| < ~1e-308).
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
    // NOTE(xsf-rust): unlike xsf v0.2.2, this uses the correct sign of the imaginary part,
    // accounts for sgn(y), checks the right factor for zero, propagates NaN, and avoids premature
    // rounding of tiny sin/cos factors (see csinpi(z) for the remaining limitation).
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
#[allow(clippy::float_cmp)]
fn tancot(xx: f64, cotflg: bool) -> f64 {
    let mut x;
    let mut sign: i32;

    /* make argument positive but save the sign */
    if xx < 0.0 {
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
    x -= 180.0 * (x / 180.0).floor();
    if cotflg {
        if x <= 90.0 {
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
        return 0.0;
    } else if x == 45.0 {
        return f64::from(sign) * 1.0;
    } else if x == 90.0 {
        return f64::INFINITY;
    }
    /* x is now transformed into [0, 90) */
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
    if x < 0.0 {
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

#[cfg(test)]
mod tests {
    use core::f64::consts::{FRAC_PI_4, PI};

    use num_complex::c64;

    #[test]
    fn test_sinpi_f64() {
        xsref::test("sinpi", "d-d", |x| crate::sinpi(x[0]));
    }

    #[test]
    fn test_sinpi_c64() {
        xsref::test("sinpi", "cd-cd", |x| crate::sinpi(c64(x[0], x[1])));
    }

    #[test]
    fn test_cospi_f64() {
        xsref::test("cospi", "d-d", |x| crate::cospi(x[0]));
    }

    #[test]
    fn test_cospi_c64() {
        xsref::test("cospi", "cd-cd", |x| crate::cospi(c64(x[0], x[1])));
    }

    #[test]
    fn test_sindg() {
        xsref::test("sindg", "d-d", |x| crate::sindg(x[0]));
    }

    #[test]
    fn test_cosdg() {
        xsref::test("cosdg", "d-d", |x| crate::cosdg(x[0]));
    }

    #[test]
    fn test_tandg() {
        xsref::test("tandg", "d-d", |x| crate::tandg(x[0]));
    }

    #[test]
    fn test_cotdg() {
        xsref::test("cotdg", "d-d", |x| crate::cotdg(x[0]));
    }

    #[test]
    fn test_cosm1() {
        xsref::test("cosm1", "d-d", |x| crate::cosm1(x[0]));
    }

    #[test]
    fn test_radian() {
        xsref::test("radian", "d_d_d-d", |x| crate::radian(x[0], x[1], x[2]));
    }

    // Edge cases that the xsref tables don't cover (they skip large inputs and outputs, and most
    // NaNs, and don't distinguish signed zeros).

    fn assert_close(actual: f64, expected: f64) {
        let ok = if expected.is_finite() {
            (actual - expected).abs() <= 1e-14 * expected.abs()
        } else {
            actual.to_bits() == expected.to_bits() || (actual.is_nan() && expected.is_nan())
        };
        assert!(ok, "actual: {actual:e}, expected: {expected:e}");
    }

    fn next_up(x: f64) -> f64 {
        f64::from_bits(x.to_bits() + 1)
    }

    #[test]
    fn test_nan() {
        for f in [
            crate::sinpi::<f64>,
            crate::cospi::<f64>,
            crate::sindg,
            crate::cosdg,
            crate::tandg,
            crate::cotdg,
            crate::cosm1,
        ] {
            assert!(f(f64::NAN).is_nan());
        }
        assert!(crate::radian(f64::NAN, 0.0, 0.0).is_nan());
        assert!(crate::sinpi(f64::INFINITY).is_nan());
        assert!(crate::cospi(f64::NEG_INFINITY).is_nan());
        assert!(crate::cosm1(f64::INFINITY).is_nan());
    }

    #[test]
    fn test_sinpi_cospi_exact() {
        assert_eq!(crate::sinpi(0.0_f64).to_bits(), 0.0_f64.to_bits());
        assert_eq!(crate::sinpi(-0.0_f64).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(crate::sinpi(0.5_f64), 1.0);
        assert_eq!(crate::sinpi(-0.5_f64), -1.0);
        assert_eq!(crate::sinpi(1.5_f64), -1.0);
        assert_eq!(crate::sinpi(1e300_f64), 0.0);
        assert_eq!(crate::cospi(0.0_f64), 1.0);
        assert_eq!(crate::cospi(-0.0_f64), 1.0);
        assert_eq!(crate::cospi(1.0_f64), -1.0);
        assert_eq!(crate::cospi(-3.0_f64), -1.0);
        // never -0.0
        for x in [0.5, -0.5, 2.5, -2.5, 1e15 + 0.5] {
            assert_eq!(crate::cospi(x).to_bits(), 0.0_f64.to_bits(), "cospi({x})");
        }
    }

    #[test]
    fn test_sindg_cosdg_exact() {
        for (x, s, c) in [
            (0.0, 0.0, 1.0),
            (90.0, 1.0, 0.0),
            (180.0, 0.0, -1.0),
            (270.0, -1.0, 0.0),
            (360.0, 0.0, 1.0),
            (-90.0, -1.0, 0.0),
            (-180.0, 0.0, -1.0),
            (450.0, 1.0, 0.0),
            (-720.0, 0.0, 1.0),
        ] {
            assert_eq!(crate::sindg(x), s, "sindg({x})");
            assert_eq!(crate::cosdg(x), c, "cosdg({x})");
        }
        assert_close(crate::sindg(30.0), 0.5);
        assert_close(crate::cosdg(60.0), 0.5);
        assert_close(crate::sindg(-45.0), -FRAC_PI_4.sin());
        assert_close(crate::cosdg(135.0), -FRAC_PI_4.cos());
    }

    #[test]
    fn test_tandg_cotdg_exact() {
        for (x, t) in [
            (0.0, 0.0),
            (45.0, 1.0),
            (135.0, -1.0),
            (-45.0, -1.0),
            (180.0, 0.0),
            (225.0, 1.0),
        ] {
            assert_eq!(crate::tandg(x), t, "tandg({x})");
        }
        for (x, c) in [
            (45.0, 1.0),
            (90.0, 0.0),
            (135.0, -1.0),
            (-45.0, -1.0),
            (270.0, 0.0),
        ] {
            assert_eq!(crate::cotdg(x), c, "cotdg({x})");
        }
        for x in [90.0, 270.0, -90.0] {
            assert!(crate::tandg(x).is_infinite(), "tandg({x})");
        }
        for x in [0.0, 180.0, -180.0] {
            assert!(crate::cotdg(x).is_infinite(), "cotdg({x})");
        }
        assert_close(crate::tandg(30.0), (PI / 6.0).tan());
        assert_close(crate::cotdg(60.0), 1.0 / (PI / 3.0).tan());
    }

    #[test]
    fn test_degrees_loss_threshold() {
        let x = 1e14;
        assert!(crate::sindg(x).abs() <= 1.0 && crate::sindg(x) != 0.0);
        assert!(crate::cosdg(x).abs() <= 1.0 && crate::cosdg(x) != 0.0);
        assert!(crate::tandg(x).is_finite() && crate::tandg(x) != 0.0);
        assert!(crate::cotdg(x).is_finite() && crate::cotdg(x) != 0.0);
        // no result: total loss of precision
        for x in [next_up(1e14), -next_up(1e14), 1e300] {
            assert_eq!(crate::sindg(x), 0.0);
            assert_eq!(crate::cosdg(x), 0.0);
            assert_eq!(crate::tandg(x), 0.0);
            assert_eq!(crate::cotdg(x), 0.0);
        }
    }

    #[test]
    fn test_cosm1_edges() {
        assert_eq!(crate::cosm1(0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(crate::cosm1(-0.0).to_bits(), 0.0_f64.to_bits());
        // the polynomial is used for |x| <= pi/4, and cos(x) - 1 beyond
        for x in [
            FRAC_PI_4,
            next_up(FRAC_PI_4),
            -FRAC_PI_4,
            -next_up(FRAC_PI_4),
        ] {
            assert_close(crate::cosm1(x), x.cos() - 1.0);
        }
        assert_close(crate::cosm1(1e-5), -5e-11 + 1e-20 / 24.0);
    }

    #[test]
    fn test_radian_edges() {
        assert_eq!(crate::radian(0.0, 0.0, 0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            crate::radian(-0.0, -0.0, -0.0).to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(crate::radian(1.0, -60.0, 0.0), 0.0);
        assert_eq!(crate::radian(-1.0, 59.0, 60.0), 0.0);
        assert_close(crate::radian(180.0, 0.0, 0.0), PI);
        assert_close(crate::radian(-90.0, 0.0, 0.0), -PI / 2.0);
        assert_close(crate::radian(0.0, 0.0, 3600.0), PI / 180.0);
        assert_close(crate::radian(0.0, -30.0, 0.0), -PI / 360.0);
    }

    #[test]
    fn test_sinpi_cospi_complex_small_imag() {
        // |pi y| < 700, both signs of x and y, just below the threshold
        for (x, y) in [
            (0.25, 222.8),
            (-0.25, -222.8),
            (0.75, -100.0),
            (-1.25, 3.0),
            (0.1, 0.0),
        ] {
            let (s, c) = (crate::sinpi(c64(x, y)), crate::cospi(c64(x, y)));
            let (sx, cx) = ((PI * x).sin(), (PI * x).cos());
            let (chy, shy) = ((PI * y).cosh(), (PI * y).sinh());
            assert_close(s.re, sx * chy);
            assert_close(s.im, cx * shy);
            assert_close(c.re, cx * chy);
            assert_close(c.im, -sx * shy);
        }
    }

    /// `a * cosh(t)` (`sgn = 1`) or `a * sinh(t)` (`sgn = sgn(t)`) for `|t| >= 700`, in the log
    /// domain, so that it's independent of the `exp(|t|/2)` scaling in `csinpi` and `ccospi`
    fn mul_cosh_sinh_large(a: f64, t: f64, sgn: f64) -> f64 {
        if a == 0.0 {
            0.0_f64.copysign(a * sgn)
        } else {
            a.signum() * sgn * (a.abs().ln() + t.abs() - core::f64::consts::LN_2).exp()
        }
    }

    fn assert_same_or_close(actual: f64, expected: f64, what: &str) {
        let ok = if expected.is_finite() && expected != 0.0 {
            (actual - expected).abs() <= 1e-12 * expected.abs()
        } else {
            actual.to_bits() == expected.to_bits() || (actual.is_nan() && expected.is_nan())
        };
        assert!(ok, "{what}: actual: {actual:e}, expected: {expected:e}");
    }

    #[test]
    fn test_sinpi_cospi_complex_large_imag() {
        // 700 <= |pi y|: below and beyond the overflow of exp(|pi y| / 2) at |y| ~ 451.9
        let xs = [
            0.0, -0.0, 0.25, -0.25, 0.5, -0.5, 0.75, 1.0, -1.0, 1.75, 2.5, 1e-10, -3e-200,
        ];
        let ys = [
            222.9,
            225.0,
            300.0,
            451.0,
            452.0,
            1000.0,
            1e300,
            f64::INFINITY,
        ];
        for x in xs {
            for y in ys.into_iter().flat_map(|y| [y, -y]) {
                let z = c64(x, y);
                let (sx, cx) = (crate::sinpi(x), crate::cospi(x));
                let (t, sgn) = (PI * y, y.signum());
                let (s, c) = (crate::sinpi(z), crate::cospi(z));
                assert_same_or_close(
                    s.re,
                    mul_cosh_sinh_large(sx, t, 1.0),
                    &format!("sinpi({z}).re"),
                );
                assert_same_or_close(
                    s.im,
                    mul_cosh_sinh_large(cx, t, sgn),
                    &format!("sinpi({z}).im"),
                );
                assert_same_or_close(
                    c.re,
                    mul_cosh_sinh_large(cx, t, 1.0),
                    &format!("cospi({z}).re"),
                );
                assert_same_or_close(
                    c.im,
                    mul_cosh_sinh_large(-sx, t, sgn),
                    &format!("cospi({z}).im"),
                );
            }
        }

        // the cases that xsf v0.2.2 got wrong
        let z = crate::cospi(c64(0.0, 1000.0));
        assert_eq!(
            (z.re, z.im.to_bits()),
            (f64::INFINITY, (-0.0_f64).to_bits())
        );
        let z = crate::sinpi(c64(0.0, -300.0));
        assert_eq!((z.re, z.im), (0.0, f64::NEG_INFINITY));
        let z = crate::cospi(c64(0.5, 1000.0));
        assert_eq!(
            (z.re.to_bits(), z.im),
            (0.0_f64.to_bits(), f64::NEG_INFINITY)
        );
        let z = crate::cospi(c64(0.25, 225.0));
        assert!(z.re > 0.0 && z.im < 0.0);

        // subnormal sin(pi x)
        let x = f64::from_bits(1);
        let (sx, z) = (crate::sinpi(x), c64(x, 222.9));
        let expected = mul_cosh_sinh_large(sx, PI * z.im, 1.0);
        assert_same_or_close(crate::sinpi(z).re, expected, "sinpi(5e-324+222.9i).re");
        assert_same_or_close(
            crate::cospi(z.conj()).im,
            expected,
            "cospi(5e-324-222.9i).im",
        );

        // NaN or infinite real part, or NaN imaginary part
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for y in [300.0, -1000.0, f64::INFINITY] {
                let (s, c) = (crate::sinpi(c64(x, y)), crate::cospi(c64(x, y)));
                assert!(s.re.is_nan() && s.im.is_nan(), "sinpi({x}+{y}i) = {s}");
                assert!(c.re.is_nan() && c.im.is_nan(), "cospi({x}+{y}i) = {c}");
            }
        }
        for x in [0.0, -0.5, 0.25, f64::NAN] {
            let (s, c) = (
                crate::sinpi(c64(x, f64::NAN)),
                crate::cospi(c64(x, f64::NAN)),
            );
            assert!(s.re.is_nan() && s.im.is_nan(), "sinpi({x}+NaNi) = {s}");
            assert!(c.re.is_nan() && c.im.is_nan(), "cospi({x}+NaNi) = {c}");
        }
    }
}
