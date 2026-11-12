use core::f64::consts::{FRAC_PI_2, LN_2};
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::utils::pow2;
use crate::xsf::c64::{C64, c64};
use crate::xsf::log::clog;

/// Square root of complex `z`
///
/// This replaces `std::sqrt(std::complex<double>)` in C++. Unlike `num_complex::Complex::sqrt`,
/// which uses the polar form, it avoids the cancellation in `|z| ± Re z`, and it scales `z` to
/// avoid premature overflow and underflow.
///
/// Translated from FreeBSD's `lib/msun/src/s_csqrt.c`:
///
/// > Copyright (c) 2007 David Schultz <das@FreeBSD.ORG>
/// > All rights reserved.
/// >
/// > Redistribution and use in source and binary forms, with or without
/// > modification, are permitted provided that the following conditions
/// > are met:
/// > 1. Redistributions of source code must retain the above copyright
/// >    notice, this list of conditions and the following disclaimer.
/// > 2. Redistributions in binary form must reproduce the above copyright
/// >    notice, this list of conditions and the following disclaimer in the
/// >    documentation and/or other materials provided with the distribution.
/// >
/// > THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
/// > ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
/// > IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
/// > ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
/// > FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
/// > DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
/// > OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
/// > HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
/// > LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
/// > OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
/// > SUCH DAMAGE.
#[allow(clippy::eq_op)] // `b - b` is NaN for infinite `b`
pub(crate) fn csqrt(z: C64) -> C64 {
    /* For avoiding overflow for components >= DBL_MAX / (1 + sqrt(2)). */
    const THRESH: f64 = 7.446_288_774_449_766e307; // 0x1.a827999fcef32p+1022
    const TWO_M1020: f64 = pow2(-1020);
    const TWO_54: f64 = pow2(54);
    const TWO_M27: f64 = pow2(-27);

    let mut a = z.re;
    let mut b = z.im;

    /* Handle special cases. */
    if a == 0.0 && b == 0.0 {
        return c64(0.0, b);
    }
    if b.is_infinite() {
        return c64(f64::INFINITY, b);
    }
    if a.is_nan() {
        let t = (b - b) / (b - b); /* raise invalid if b is not a NaN */
        return c64(a + 0.0 + t, a + 0.0 + t); /* NaN + NaN i */
    }
    if a.is_infinite() {
        /*
         * csqrt(inf + NaN i)  = inf +  NaN i
         * csqrt(inf + y i)    = inf +  0 i
         * csqrt(-inf + NaN i) = NaN +- inf i
         * csqrt(-inf + y i)   = 0   +  inf i
         */
        if a.is_sign_negative() {
            return c64((b - b).abs(), a.copysign(b));
        }
        return c64(a, (b - b).copysign(b));
    }
    if b.is_nan() {
        let t = (a - a) / (a - a); /* raise invalid */
        return c64(b + 0.0 + t, b + 0.0 + t); /* NaN + NaN i */
    }

    /* Scale to avoid overflow. */
    let mut scale = if a.abs() >= THRESH || b.abs() >= THRESH {
        /*
         * Don't scale a or b if this might give (spurious)
         * underflow.  Then the unscaled value is an equivalent
         * infinitesmal (or 0).
         */
        if a.abs() >= TWO_M1020 {
            a *= 0.25;
        }
        if b.abs() >= TWO_M1020 {
            b *= 0.25;
        }
        2.0
    } else {
        1.0
    };

    /* Scale to reduce inaccuracies when both components are denormal. */
    if a.abs() < f64::MIN_POSITIVE && b.abs() < f64::MIN_POSITIVE {
        a *= TWO_54;
        b *= TWO_54;
        scale = TWO_M27;
    }

    /* Algorithm 312, CACM vol 10, Oct 1967. */
    let (rx, ry);
    if a >= 0.0 {
        let t = ((a + a.hypot(b)) * 0.5).sqrt();
        rx = scale * t;
        // (unlike in FreeBSD) if `a == 0`, then `b / (2 t) = t sign(b)` in exact arithmetic, and
        // the latter is correctly rounded (unless `0.5 |b|` is subnormal and inexact), whereas the
        // former can be 1 ulp off
        ry = if a == 0.0 {
            (scale * t).copysign(b)
        } else {
            scale * b / (2.0 * t)
        };
    } else {
        let t = ((-a + a.hypot(b)) * 0.5).sqrt();
        rx = scale * b.abs() / (2.0 * t);
        ry = (scale * t).copysign(b);
    }

    c64(rx, ry)
}

/// `std::pow(z, y)` for complex `z` and real `y`, as in libstdc++
///
/// Unlike in xsf, `z^0` is 1 for all `z`, and the modulus is `|z|^y` instead of `exp(y log |z|)`,
/// which amplifies the relative rounding error of `log |z|` by `|y log |z||`. Only for a rounded
/// `|z|` in `[0.5, 2]`, where [`clog`] computes `log |z|` in double-double precision, or outside
/// the normal range, the modulus is as in xsf.
pub(crate) fn cpow(z: C64, y: f64) -> C64 {
    if y == 0.0 {
        C64::ONE
    } else if z.im == 0.0 && z.re > 0.0 {
        c64(z.re.powf(y), 0.0)
    } else {
        let t = clog(z);
        let r = z.norm();
        // `r` is exact on the axes
        let exact = z.re == 0.0 || z.im == 0.0;
        let rho = if exact || !(0.5..=2.0).contains(&r) && r.is_normal() {
            r.powf(y)
        } else {
            (y * t.re).exp()
        };
        // std::polar(rho, y * t.imag())
        let (sin, cos) = (y * t.im).sin_cos();
        c64(rho * cos, rho * sin)
    }
}

/// Inverse hyperbolic tangent of real `x`
///
/// This replaces `std::atanh(double)` in C++. Unlike `f64::atanh`, which computes
/// `0.5 * log1p(2x / (1 - x))` for all `x`, it is as accurate as glibc's `atanh`.
///
/// Translated from FreeBSD's `lib/msun/src/e_atanh.c`:
///
/// > Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
/// >
/// > Developed at SunSoft, a Sun Microsystems, Inc. business.
/// > Permission to use, copy, modify, and distribute this
/// > software is freely granted, provided that this notice
/// > is preserved.
#[allow(clippy::eq_op)] // `x - x` is NaN for infinite `x`
pub(crate) fn atanh(x: f64) -> f64 {
    /* atanh(x)
     * Method :
     *    1.Reduced x to positive by atanh(-x) = -atanh(x)
     *    2.For x>=0.5
     *                  1              2x                          x
     *      atanh(x) = --- * log(1 + -------) = 0.5 * log1p(2 * --------)
     *                  2             1 - x                      1 - x
     *
     *      For x<0.5
     *      atanh(x) = 0.5*log1p(2x+2x*x/(1-x))
     *
     * Special cases:
     *      atanh(x) is NaN if |x| > 1 with signal;
     *      atanh(NaN) is that NaN with no signal;
     *      atanh(+-1) is +-INF with signal.
     */
    const TWO_M28: f64 = pow2(-28);

    let xa = x.abs();
    if xa > 1.0 || xa.is_nan() {
        /* |x|>1 */
        return (x - x) / (x - x);
    }
    if xa == 1.0 {
        return x / 0.0;
    }
    if xa < TWO_M28 {
        return x; /* x<2**-28 */
    }
    let t = if xa < 0.5 {
        /* x < 0.5 */
        let t = xa + xa;
        0.5 * (t + t * xa / (1.0 - xa)).ln_1p()
    } else {
        0.5 * ((xa + xa) / (1.0 - xa)).ln_1p()
    };
    if x.is_sign_negative() { -t } else { t }
}

/// Arc tangent of complex `z`
///
/// This replaces `std::atan(std::complex<double>)` in C++. Unlike `num_complex::Complex::atan`,
/// which computes the difference of two logarithms, it doesn't lose precision near zero or near the
/// branch points `z = +-i`.
///
/// Translated from FreeBSD's `lib/msun/src/catrig.c`:
///
/// > Copyright (c) 2012 Stephen Montgomery-Smith <stephen@FreeBSD.ORG>
/// > All rights reserved.
/// >
/// > Redistribution and use in source and binary forms, with or without
/// > modification, are permitted provided that the following conditions
/// > are met:
/// > 1. Redistributions of source code must retain the above copyright
/// >    notice, this list of conditions and the following disclaimer.
/// > 2. Redistributions in binary form must reproduce the above copyright
/// >    notice, this list of conditions and the following disclaimer in the
/// >    documentation and/or other materials provided with the distribution.
/// >
/// > THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
/// > ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
/// > IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
/// > ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
/// > FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
/// > DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
/// > OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
/// > HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
/// > LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
/// > OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
/// > SUCH DAMAGE.
pub(crate) fn catan(z: C64) -> C64 {
    /*
     * catan(z) = reverse(catanh(reverse(z)))
     * where reverse(x + I*y) = y + I*x = I*conj(z).
     */
    let w = catanh(c64(z.im, z.re));
    c64(w.im, w.re)
}

/*
 * catanh(z) = log((1+z)/(1-z)) / 2
 *           = log1p(4*x / |z-1|^2) / 4
 *             + I * atan2(2*y, (1-x)*(1+x)-y*y) / 2
 *
 * catanh(z) = z + O(z^3)   as z -> 0
 *
 * catanh(z) = 1/z + sign(y)*I*PI/2 + O(1/z^3)   as z -> infinity
 * The above formula works for the real part as well, because
 * Re(catanh(z)) = x/|z|^2 + O(x/z^4)
 *    as z -> infinity, uniformly in x
 */
fn catanh(z: C64) -> C64 {
    const RECIP_EPSILON: f64 = 1.0 / f64::EPSILON;
    const SQRT_3_EPSILON: f64 = 2.580_956_827_951_785e-8; /*  0x1bb67ae8584caa.0p-78 */

    let x = z.re;
    let y = z.im;
    let ax = x.abs();
    let ay = y.abs();

    /* This helps handle many cases. */
    if y == 0.0 && ax <= 1.0 {
        return c64(atanh(x), y);
    }

    /* To ensure the same accuracy as atan(), and to filter out z = 0. */
    if x == 0.0 {
        return c64(x, y.atan());
    }

    if x.is_nan() || y.is_nan() {
        /* catanh(+-Inf + I*NaN) = +-0 + I*NaN */
        if x.is_infinite() {
            return c64(0.0_f64.copysign(x), y + y);
        }
        /* catanh(NaN + I*+-Inf) = sign(NaN)0 + I*+-PI/2 */
        if y.is_infinite() {
            return c64(0.0_f64.copysign(x), FRAC_PI_2.copysign(y));
        }
        /*
         * All other cases involving NaN return NaN + I*NaN.
         * C99 leaves it optional whether to raise invalid if one of
         * the arguments is not NaN, so we opt not to raise it.
         */
        return c64(x + y, x + y);
    }

    if ax > RECIP_EPSILON || ay > RECIP_EPSILON {
        return c64(real_part_reciprocal(x, y), FRAC_PI_2.copysign(y));
    }

    if ax < SQRT_3_EPSILON / 2.0 && ay < SQRT_3_EPSILON / 2.0 {
        /*
         * z = 0 was filtered out above.  All other cases must raise
         * inexact, but this is the only case that needs to do it
         * explicitly.
         */
        return z;
    }

    let rx = if ax == 1.0 && ay < f64::EPSILON {
        (LN_2 - ay.ln()) / 2.0
    } else {
        (4.0 * ax / sum_squares(ax - 1.0, ay)).ln_1p() / 4.0
    };

    let ry = if ax == 1.0 {
        2.0_f64.atan2(-ay) / 2.0
    } else if ay < f64::EPSILON {
        (2.0 * ay).atan2((1.0 - ax) * (1.0 + ax)) / 2.0
    } else {
        (2.0 * ay).atan2((1.0 - ax) * (1.0 + ax) - ay * ay) / 2.0
    };

    c64(rx.copysign(x), ry.copysign(y))
}

/*
 * sum_squares(x,y) = x*x + y*y (or just x*x if y*y would underflow).
 * Assumes x*x and y*y will not overflow.
 * Assumes x and y are finite.
 * Assumes y is non-negative.
 * Assumes fabs(x) >= DBL_EPSILON.
 */
#[inline]
fn sum_squares(x: f64, y: f64) -> f64 {
    const SQRT_MIN: f64 = pow2(-511); /* >= sqrt(DBL_MIN) */

    /* Avoid underflow when y is small. */
    if y < SQRT_MIN {
        return x * x;
    }

    x * x + y * y
}

/*
 * real_part_reciprocal(x, y) = Re(1/(x+I*y)) = x/(x*x + y*y).
 * Assumes x and y are not NaN, and one of x and y is larger than
 * RECIP_EPSILON.  We avoid unwarranted underflow.  It is important to not use
 * the code creal(1/z), because the imaginary part may produce an unwanted
 * underflow.
 */
#[allow(clippy::cast_possible_wrap)]
fn real_part_reciprocal(mut x: f64, mut y: f64) -> f64 {
    const BIAS: i32 = 1023; /* DBL_MAX_EXP - 1 */
    const CUTOFF: i32 = 53 / 2 + 1; /* DBL_MANT_DIG / 2 + 1, just half or 1 guard digit */

    /*
     * This code is inspired by the C99 document n1124.pdf, Section G.5.1,
     * example 2.
     */
    let ix = ((x.to_bits() >> 32) as u32 & 0x7ff0_0000) as i32;
    let iy = ((y.to_bits() >> 32) as u32 & 0x7ff0_0000) as i32;
    if ix - iy >= CUTOFF << 20 || x.is_infinite() {
        return 1.0 / x; /* +-Inf -> +-0 is special */
    }
    if iy - ix >= CUTOFF << 20 {
        return x / y / y; /* should avoid double div, but hard */
    }
    if ix <= (BIAS + 1024 / 2 - CUTOFF) << 20 {
        return x / (x * x + y * y);
    }
    let scale = f64::from_bits(u64::from((0x7ff0_0000 - ix) as u32) << 32); /* 2**(1-ilogb(x)) */
    x *= scale;
    y *= scale;
    x / (x * x + y * y) * scale
}

/// The arithmetic of `T` in C++, i.e. `double` or `std::complex<double>`
pub(crate) trait Num:
    Copy
    + From<f64>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    /// `std::sqrt(x)`
    fn sqrt(self) -> Self;
    /// `std::sin(x)`
    fn sin(self) -> Self;
    /// `std::cos(x)`
    fn cos(self) -> Self;
    /// `std::abs(x)`
    fn abs(self) -> f64;
    /// `std::real(x)`
    fn real(self) -> f64;
    /// `std::imag(x)`
    fn imag(self) -> f64;
    /// Whether `x` (each part, if complex) is finite
    fn is_finite(self) -> bool;
}

impl Num for f64 {
    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    #[inline]
    fn sin(self) -> Self {
        f64::sin(self)
    }
    #[inline]
    fn cos(self) -> Self {
        f64::cos(self)
    }
    #[inline]
    fn abs(self) -> f64 {
        f64::abs(self)
    }
    #[inline]
    fn real(self) -> f64 {
        self
    }
    #[inline]
    fn imag(self) -> f64 {
        0.0
    }
    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
}

impl Num for C64 {
    #[inline]
    fn sqrt(self) -> Self {
        C64::sqrt(self)
    }
    #[inline]
    fn sin(self) -> Self {
        C64::sin(self)
    }
    #[inline]
    fn cos(self) -> Self {
        C64::cos(self)
    }
    #[inline]
    fn abs(self) -> f64 {
        self.norm()
    }
    #[inline]
    fn real(self) -> f64 {
        self.re
    }
    #[inline]
    fn imag(self) -> f64 {
        self.im
    }
    #[inline]
    fn is_finite(self) -> bool {
        C64::is_finite(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::pow2;
    use crate::xsf::c64::{C64, c64};
    use core::f64::consts::FRAC_1_SQRT_2;
    use num_complex::Complex;

    #[test]
    fn test_csqrt_imaginary() {
        // correctly rounded, unlike in FreeBSD
        let r = super::csqrt(c64(0.0, -1.0));
        assert_eq!(
            [r.re.to_bits(), r.im.to_bits()],
            [FRAC_1_SQRT_2.to_bits(), (-FRAC_1_SQRT_2).to_bits()]
        );
    }

    #[test]
    fn test_cpow() {
        let rel = |w: C64, r: f64| (w.norm() / r - 1.0).abs();
        // `|z|^y` on an axis, where `exp(y log |z|)` loses 3 digits
        let w = super::cpow(c64(-3.0, 0.0), 600.0);
        assert!(rel(w, 3.0_f64.powf(600.0)) < 4.0 * f64::EPSILON);
        // `exp(y log |z|)` where `|z|` rounds to 1
        let w = super::cpow(c64(1.0, pow2(-27)), pow2(52));
        assert!(rel(w, 0.125_f64.exp()) < 4.0 * f64::EPSILON);
        // and where `|z|` overflows
        let w = super::cpow(c64(f64::MAX, f64::MAX), -0.5);
        assert!(rel(w, pow2(-512) / 2.0_f64.powf(0.25)) < 1e-12);

        // `z^0`, and both sides of the branch cut
        assert_eq!(super::cpow(c64(f64::NAN, 0.0), 0.0), c64(1.0, 0.0));
        assert_eq!(super::cpow(c64(-4.0, 0.0), 0.5).im, 2.0);
        assert_eq!(super::cpow(c64(-4.0, -0.0), 0.5).im, -2.0);
    }

    #[test]
    fn test_atanh() {
        // mpmath
        let actual = [0.5, 1e-10, -0.9].map(super::atanh);
        let expected = [0.549_306_144_334_054_9, 1e-10, -1.472_219_489_583_220_4];
        crate::np_assert_allclose!(&actual, &expected, rtol = 4e-16);

        assert_eq!(super::atanh(-0.0).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(super::atanh(-1.0), f64::NEG_INFINITY);
        assert!(super::atanh(1.5).is_nan());
    }

    #[test]
    fn test_catan() {
        // mpmath
        let actual = [
            c64(0.5, 0.5),
            c64(1e-5, 1.0),
            c64(2.0, -3.0),
            c64(1e300, 1.0),
        ]
        .map(|z| Complex::from(super::catan(z)));
        let expected = [
            c64(0.553_574_358_897_045_2, 0.402_359_478_108_525_07),
            c64(0.785_400_663_397_448_3, 6.103_036_322_771_337),
            c64(1.409_921_049_596_575_5, -0.229_072_682_968_538_78),
            c64(core::f64::consts::FRAC_PI_2, 0.0),
        ]
        .map(Complex::from);
        crate::np_assert_allclose!(&actual, &expected, rtol = 4e-16);

        // signed zeros and the branch points
        let r = super::catan(c64(-0.0, 0.0));
        assert_eq!([r.re.to_bits(), r.im.to_bits()], [(-0.0_f64).to_bits(), 0]);
        assert_eq!(super::catan(c64(0.0, 1.0)), c64(0.0, f64::INFINITY));
        assert_eq!(super::catan(c64(0.0, -1.0)), c64(0.0, f64::NEG_INFINITY));
    }
}
