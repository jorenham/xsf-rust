use core::f64::consts::{FRAC_PI_2, LN_2};

use num_complex::{Complex, c64};

use crate::xsf::log::clog;

/// Complex division `a / b` using Smith's algorithm, which avoids the overflow and underflow of
/// `|b|^2` in the naive formula. Tiny operands are first scaled by the same power of 2, so that
/// the intermediate results don't lose precision by becoming subnormal.
///
/// This replaces the complex division in C++, i.e. GCC's `__divdc3`. Like `__divdc3`, it recovers
/// the infinities and zeros that would otherwise be computed as NaN + NaN i (C99 Annex G), but
/// unlike `__divdc3`, it doesn't scale huge operands (e.g. `1 / (MAX + i MAX)` underflows to zero).
///
/// R. L. Smith, "Algorithm 116: Complex division", Communications of the ACM 5(8), 1962.
#[inline]
pub(crate) fn cdiv(mut a: Complex<f64>, mut b: Complex<f64>) -> Complex<f64> {
    const TWO_M400: f64 = f64::from_bits((1023 - 400) << 52); // 2^-400
    const TWO_600: f64 = f64::from_bits((1023 + 600) << 52); // 2^600

    if b.re.abs().max(b.im.abs()) < TWO_M400 {
        a *= TWO_600;
        b *= TWO_600;
    }
    let w = if b.re.abs() >= b.im.abs() {
        let r = b.im / b.re;
        let den = b.re + b.im * r;
        Complex::new((a.re + a.im * r) / den, (a.im - a.re * r) / den)
    } else {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        Complex::new((a.re * r + a.im) / den, (a.im * r - a.re) / den)
    };

    // Recover infinities and zeros that computed as NaN+iNaN; the only
    // cases are nonzero/zero, infinite/finite, and finite/infinite.
    if w.re.is_nan() && w.im.is_nan() {
        let (Complex { re: a, im: b }, Complex { re: c, im: d }) = (a, b);
        if c == 0.0 && d == 0.0 && (!a.is_nan() || !b.is_nan()) {
            let inf = f64::INFINITY.copysign(c);
            return Complex::new(inf * a, inf * b);
        } else if (a.is_infinite() || b.is_infinite()) && c.is_finite() && d.is_finite() {
            let a = box_inf(a);
            let b = box_inf(b);
            return Complex::new(
                f64::INFINITY * (a * c + b * d),
                f64::INFINITY * (b * c - a * d),
            );
        } else if (c.is_infinite() || d.is_infinite()) && a.is_finite() && b.is_finite() {
            let c = box_inf(c);
            let d = box_inf(d);
            return Complex::new(0.0 * (a * c + b * d), 0.0 * (b * c - a * d));
        }
    }
    w
}

/// Complex multiplication `a * b`, like GCC's `__muldc3` (C99 Annex G)
///
/// This replaces the complex multiplication in C++. For finite results it is the same as the naive
/// formula (i.e. `num_complex`'s `*`), but it recovers the infinities that would otherwise be
/// computed as NaN + NaN i, e.g. `(inf + inf i) * (1 + 0i)`.
#[inline]
pub(crate) fn cmul(a: Complex<f64>, b: Complex<f64>) -> Complex<f64> {
    let (
        Complex {
            re: mut a,
            im: mut b,
        },
        Complex {
            re: mut c,
            im: mut d,
        },
    ) = (a, b);
    let (ac, bd, ad, bc) = (a * c, b * d, a * d, b * c);
    let (mut x, mut y) = (ac - bd, ad + bc);

    if x.is_nan() && y.is_nan() {
        // Recover infinities that computed as NaN+iNaN.
        let mut recalc = false;
        if a.is_infinite() || b.is_infinite() {
            // z is infinite.  "Box" the infinity and change NaNs in
            // the other factor to 0.
            a = box_inf(a);
            b = box_inf(b);
            c = nan_to_zero(c);
            d = nan_to_zero(d);
            recalc = true;
        }
        if c.is_infinite() || d.is_infinite() {
            // w is infinite.  "Box" the infinity and change NaNs in
            // the other factor to 0.
            c = box_inf(c);
            d = box_inf(d);
            a = nan_to_zero(a);
            b = nan_to_zero(b);
            recalc = true;
        }
        if !recalc && (ac.is_infinite() || bd.is_infinite() || ad.is_infinite() || bc.is_infinite())
        {
            // Recover infinities from overflow by changing NaNs to 0.
            a = nan_to_zero(a);
            b = nan_to_zero(b);
            c = nan_to_zero(c);
            d = nan_to_zero(d);
            recalc = true;
        }
        if recalc {
            x = f64::INFINITY * (a * c - b * d);
            y = f64::INFINITY * (a * d + b * c);
        }
    }
    Complex::new(x, y)
}

/// `isnan(x) ? copysign(0, x) : x`, used by [`cmul`]
fn nan_to_zero(x: f64) -> f64 {
    if x.is_nan() { 0.0_f64.copysign(x) } else { x }
}

/// `copysign(isinf(x) ? 1 : 0, x)`, used to "box" infinities in C99 Annex G
fn box_inf(x: f64) -> f64 {
    (if x.is_infinite() { 1.0_f64 } else { 0.0 }).copysign(x)
}

/// `x - z` for real `x` as in C++, where the imaginary part is `-z.im` instead of `0.0 - z.im`
pub(crate) fn rsub(x: f64, z: Complex<f64>) -> Complex<f64> {
    Complex::new(x - z.re, -z.im)
}

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
#[allow(clippy::eq_op, clippy::manual_midpoint)] // `b - b` is NaN for infinite `b`
pub(crate) fn csqrt(z: Complex<f64>) -> Complex<f64> {
    /* For avoiding overflow for components >= DBL_MAX / (1 + sqrt(2)). */
    const THRESH: f64 = 7.446_288_774_449_766e307; // 0x1.a827999fcef32p+1022
    const TWO_M1020: f64 = f64::from_bits((1023 - 1020) << 52); // 0x1p-1020
    const TWO_54: f64 = f64::from_bits((1023 + 54) << 52); // 0x1p54
    const TWO_M27: f64 = f64::from_bits((1023 - 27) << 52); // 0x1p-27

    let mut a = z.re;
    let mut b = z.im;

    /* Handle special cases. */
    if a == 0.0 && b == 0.0 {
        return Complex::new(0.0, b);
    }
    if b.is_infinite() {
        return Complex::new(f64::INFINITY, b);
    }
    if a.is_nan() {
        let t = (b - b) / (b - b); /* raise invalid if b is not a NaN */
        return Complex::new(a + 0.0 + t, a + 0.0 + t); /* NaN + NaN i */
    }
    if a.is_infinite() {
        /*
         * csqrt(inf + NaN i)  = inf +  NaN i
         * csqrt(inf + y i)    = inf +  0 i
         * csqrt(-inf + NaN i) = NaN +- inf i
         * csqrt(-inf + y i)   = 0   +  inf i
         */
        if a.is_sign_negative() {
            return Complex::new((b - b).abs(), a.copysign(b));
        }
        return Complex::new(a, (b - b).copysign(b));
    }
    if b.is_nan() {
        let t = (a - a) / (a - a); /* raise invalid */
        return Complex::new(b + 0.0 + t, b + 0.0 + t); /* NaN + NaN i */
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

    Complex::new(rx, ry)
}

/// `std::pow(z, y)` for complex `z` and real `y`, as in libstdc++
pub(crate) fn cpow(z: Complex<f64>, y: f64) -> Complex<f64> {
    if z.im == 0.0 && z.re > 0.0 {
        return c64(z.re.powf(y), 0.0);
    }
    let t = clog(z);
    // std::polar(exp(y * t.real()), y * t.imag())
    let rho = (y * t.re).exp();
    let theta = y * t.im;
    c64(rho * theta.cos(), rho * theta.sin())
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
#[allow(clippy::eq_op, clippy::float_cmp)] // `x - x` is NaN for infinite `x`
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
    const TWO_M28: f64 = f64::from_bits((1023 - 28) << 52); // 2^-28

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
pub(crate) fn catan(z: Complex<f64>) -> Complex<f64> {
    /*
     * catan(z) = reverse(catanh(reverse(z)))
     * where reverse(x + I*y) = y + I*x = I*conj(z).
     */
    let w = catanh(Complex::new(z.im, z.re));
    Complex::new(w.im, w.re)
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
#[allow(clippy::float_cmp)]
fn catanh(z: Complex<f64>) -> Complex<f64> {
    const RECIP_EPSILON: f64 = 1.0 / f64::EPSILON;
    const SQRT_3_EPSILON: f64 = 2.580_956_827_951_785e-8; /*  0x1bb67ae8584caa.0p-78 */

    let x = z.re;
    let y = z.im;
    let ax = x.abs();
    let ay = y.abs();

    /* This helps handle many cases. */
    if y == 0.0 && ax <= 1.0 {
        return Complex::new(atanh(x), y);
    }

    /* To ensure the same accuracy as atan(), and to filter out z = 0. */
    if x == 0.0 {
        return Complex::new(x, y.atan());
    }

    if x.is_nan() || y.is_nan() {
        /* catanh(+-Inf + I*NaN) = +-0 + I*NaN */
        if x.is_infinite() {
            return Complex::new(0.0_f64.copysign(x), y + y);
        }
        /* catanh(NaN + I*+-Inf) = sign(NaN)0 + I*+-PI/2 */
        if y.is_infinite() {
            return Complex::new(0.0_f64.copysign(x), FRAC_PI_2.copysign(y));
        }
        /*
         * All other cases involving NaN return NaN + I*NaN.
         * C99 leaves it optional whether to raise invalid if one of
         * the arguments is not NaN, so we opt not to raise it.
         */
        return Complex::new(x + y, x + y);
    }

    if ax > RECIP_EPSILON || ay > RECIP_EPSILON {
        return Complex::new(real_part_reciprocal(x, y), FRAC_PI_2.copysign(y));
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

    Complex::new(rx.copysign(x), ry.copysign(y))
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
    const SQRT_MIN: f64 = f64::from_bits((1023 - 511) << 52); /* 0x1p-511 >= sqrt(DBL_MIN) */

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
#[allow(clippy::cast_possible_wrap, clippy::cast_sign_loss)]
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

#[cfg(test)]
mod tests {
    use core::f64::consts::FRAC_1_SQRT_2;
    use num_complex::c64;

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
        .map(super::catan);
        let expected = [
            c64(0.553_574_358_897_045_2, 0.402_359_478_108_525_07),
            c64(0.785_400_663_397_448_3, 6.103_036_322_771_337),
            c64(1.409_921_049_596_575_5, -0.229_072_682_968_538_78),
            c64(core::f64::consts::FRAC_PI_2, 0.0),
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 4e-16);

        // signed zeros and the branch points
        let r = super::catan(c64(-0.0, 0.0));
        assert_eq!([r.re.to_bits(), r.im.to_bits()], [(-0.0_f64).to_bits(), 0]);
        assert_eq!(super::catan(c64(0.0, 1.0)), c64(0.0, f64::INFINITY));
        assert_eq!(super::catan(c64(0.0, -1.0)), c64(0.0, f64::NEG_INFINITY));
    }
}
