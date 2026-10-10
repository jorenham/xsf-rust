//! Translated into pure Rust from `xsf/cephes/gamma.h` (xsf v0.2.2).
//!
//! Translated into C++ by SciPy developers in 2024 from the Cephes Math Library Release 2.2
//! (July, 1992), Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the parity of `floor(|x|)` is determined with an `i64` instead of an `int`, which
//! overflows (UB) for `|x| > 2^31`.

use crate::xsf::cephes::polevl::{p1evl, polevl};
use crate::xsf::trig::cephes_sinpi;

use core::f64::consts::PI;

// cephes/const.h

/// sqrt(2*pi)
const SQRT2PI: f64 = 2.506_628_274_631_000_7;
/// log(pi)
const LOGPI: f64 = 1.144_729_885_849_400_2;
/// Largest x such that Gamma(x) is finite
const MAXGAM: f64 = 171.624_376_956_302_7;

// cephes/gamma.h

const GAMMA_P: [f64; 7] = [
    1.601_195_224_767_518_5e-4,
    1.191_351_470_065_863_8e-3,
    1.042_137_975_617_615_8e-2,
    4.763_678_004_571_372e-2,
    0.207_448_227_648_435_98,
    0.494_214_826_801_497_1,
    1.0,
];

const GAMMA_Q: [f64; 8] = [
    -2.315_818_733_241_201_4e-5,
    5.396_055_804_933_034e-4,
    -4.456_419_138_517_973e-3,
    1.181_397_852_220_604_3e-2,
    3.582_363_986_054_986_5e-2,
    -0.234_591_795_718_243_35,
    7.143_049_170_302_73e-2,
    1.0,
];

/* Stirling's formula for the Gamma function */
const GAMMA_STIR: [f64; 5] = [
    7.873_113_957_930_937e-4,
    -2.295_499_616_133_781_3e-4,
    -2.681_326_178_057_812_4e-3,
    3.472_222_216_054_586_6e-3,
    8.333_333_333_334_822e-2,
];

const MAXSTIR: f64 = 143.016_08;

/* Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
fn stirf(x: f64) -> f64 {
    if x >= MAXGAM {
        return f64::INFINITY;
    }
    let mut w = 1.0 / x;
    w = 1.0 + w * polevl(w, &GAMMA_STIR);
    let mut y = x.exp();
    if x > MAXSTIR {
        /* Avoid overflow in pow() */
        let v = x.powf(0.5 * x - 0.25);
        y = v * (v / y);
    } else {
        y = x.powf(x - 0.5) / y;
    }
    SQRT2PI * y * w
}

/// Gamma function
///
/// Returns Gamma function of the argument.  The result is
/// correctly signed.
///
/// Arguments |x| <= 34 are reduced by recurrence and the function
/// approximated by a rational function of degree 6/7 in the
/// interval (2,3).  Large arguments are handled by Stirling's
/// formula. Large negative arguments are made positive using
/// a reflection formula.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE    -170, -33     20000       7.5e-16     1.9e-16
///    IEEE     -33,  33     20000       1.1e-15     2.1e-16
///    IEEE      33, 171.6   20000       5.7e-16     1.5e-16
/// ```
///
/// Error for arguments outside the test range will be larger
/// owing to error amplification by the exponential function.
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn gamma(mut x: f64) -> f64 {
    let mut sgngam = 1;

    if !x.is_finite() {
        if x > 0.0 {
            // gamma(+inf) = +inf
            return x;
        }
        // gamma(NaN) and gamma(-inf) both should equal NaN.
        return f64::NAN;
    }

    if x == 0.0 {
        /* For pole at zero, value depends on sign of zero.
         * +inf when approaching from right, -inf when approaching
         * from left. */
        return f64::INFINITY.copysign(x);
    }

    let q = x.abs();

    if q > 33.0 {
        let mut z;
        if x < 0.0 {
            let mut p = q.floor();
            if p == q {
                // x is a negative integer. This is a pole.
                return f64::NAN;
            }
            let i = p as i64;
            if (i & 1) == 0 {
                sgngam = -1;
            }
            z = q - p;
            if z > 0.5 {
                p += 1.0;
                z = q - p;
            }
            z = q * cephes_sinpi(z);
            if z == 0.0 {
                return f64::from(sgngam) * f64::INFINITY;
            }
            z = z.abs();
            z = PI / (z * stirf(q));
        } else {
            z = stirf(x);
        }
        return f64::from(sgngam) * z;
    }

    let mut z = 1.0;
    while x >= 3.0 {
        x -= 1.0;
        z *= x;
    }

    while x < 0.0 {
        if x > -1.0e-9 {
            return gamma_small(x, z);
        }
        z /= x;
        x += 1.0;
    }

    while x < 2.0 {
        if x < 1.0e-9 {
            return gamma_small(x, z);
        }
        z /= x;
        x += 1.0;
    }

    if x == 2.0 {
        return z;
    }

    x -= 2.0;
    let p = polevl(x, &GAMMA_P);
    let q = polevl(x, &GAMMA_Q);
    z * p / q
}

/// The `small:` label of `Gamma`
#[allow(clippy::float_cmp)]
fn gamma_small(x: f64, z: f64) -> f64 {
    if x == 0.0 {
        /* For this to have happened, x must have started as a negative integer. */
        f64::NAN
    } else {
        z / ((1.0 + 0.577_215_664_901_532_9 * x) * x)
    }
}

/* A[]: Stirling's formula expansion of log Gamma
 * B[], C[]: log Gamma function between 2 and 3
 */
const GAMMA_A: [f64; 5] = [
    8.116_141_674_705_085e-4,
    -5.950_619_042_843_014e-4,
    7.936_503_404_577_169e-4,
    -2.777_777_777_300_997e-3,
    8.333_333_333_333_319e-2,
];

const GAMMA_B: [f64; 6] = [
    -1_378.251_525_691_208_6,
    -38_801.631_513_463_784,
    -331_612.992_738_871_2,
    -1_162_370.974_927_623,
    -1_721_737.008_208_396_6,
    -853_555.664_245_765_4,
];

const GAMMA_C: [f64; 6] = [
    /* 1.00000000000000000000E0, */
    -351.815_701_436_523_45,
    -17_064.210_665_188_115,
    -220_528.590_553_854_45,
    -1_139_334.443_679_825_2,
    -2_532_523.071_775_829_4,
    -2_018_891.414_335_327_7,
];

/* log( sqrt( 2*pi ) ) */
const LS2PI: f64 = 0.918_938_533_204_672_8;

const MAXLGM: f64 = 2.556_348e305;

/* In xsf, optimizations for this function are disabled on 32 bit systems when compiling with GCC,
 * because they can result in degraded precision for this asymptotic approximation. Rust never
 * reassociates or contracts floating point operations, so that isn't needed here. */
fn lgam_large_x(x: f64) -> f64 {
    let q = (x - 0.5) * x.ln() - x + LS2PI;
    if x > 1.0e8 {
        return q;
    }
    let mut p = 1.0 / (x * x);
    p = ((7.936_507_936_507_937e-4 * p - 2.777_777_777_777_778e-3) * p + 8.333_333_333_333_333e-2)
        / x;
    q + p
}

/// [`lgam`] and the sign of the Gamma function, returned as `(lgam(x), sign)` instead of through
/// the `int *sign` out-parameter
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn lgam_sgn(mut x: f64) -> (f64, i32) {
    let mut sign = 1;

    if !x.is_finite() {
        return (x, sign);
    }

    if x < -34.0 {
        let q = -x;
        let w;
        (w, sign) = lgam_sgn(q);
        let mut p = q.floor();
        if p == q {
            // lgsing:
            return (f64::INFINITY, sign);
        }
        let i = p as i64;
        if (i & 1) == 0 {
            sign = -1;
        } else {
            sign = 1;
        }
        let mut z = q - p;
        if z > 0.5 {
            p += 1.0;
            z = p - q;
        }
        z = q * cephes_sinpi(z);
        if z == 0.0 {
            // goto lgsing;
            return (f64::INFINITY, sign);
        }
        /*     z = log(M_PI) - log( z ) - w; */
        z = LOGPI - z.ln() - w;
        return (z, sign);
    }

    if x < 13.0 {
        let mut z = 1.0;
        let mut p = 0.0;
        let mut u = x;
        while u >= 3.0 {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while u < 2.0 {
            if u == 0.0 {
                // goto lgsing;
                return (f64::INFINITY, sign);
            }
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if z < 0.0 {
            sign = -1;
            z = -z;
        } else {
            sign = 1;
        }
        if u == 2.0 {
            return (z.ln(), sign);
        }
        p -= 2.0;
        x += p;
        p = x * polevl(x, &GAMMA_B) / p1evl(x, &GAMMA_C);
        return (z.ln() + p, sign);
    }

    if x > MAXLGM {
        return (f64::from(sign) * f64::INFINITY, sign);
    }

    if x >= 1000.0 {
        return (lgam_large_x(x), sign);
    }

    let q = (x - 0.5) * x.ln() - x + LS2PI;
    let p = 1.0 / (x * x);
    (q + polevl(p, &GAMMA_A) / x, sign)
}

/// Natural logarithm of Gamma function
///
/// Returns the base e (2.718...) logarithm of the absolute
/// value of the Gamma function of the argument.
///
/// For arguments greater than 13, the logarithm of the Gamma
/// function is approximated by the logarithmic version of
/// Stirling's formula using a polynomial approximation of
/// degree 4. Arguments between -33 and +33 are reduced by
/// recurrence to the interval \[2,3\] of a rational approximation.
/// The cosecant reflection formula is employed for arguments
/// less than -33.
///
/// Arguments greater than MAXLGM return INFINITY and an error
/// message.  MAXLGM = 2.556348e305 for IEEE arithmetic.
///
/// ACCURACY:
///
/// ```text
/// arithmetic      domain        # trials     peak         rms
///    IEEE    0, 3                 28000     4.2e-16     7.7e-17
///    IEEE    2.718, 2.556e305     40000     3.4e-16     6.8e-17
/// ```
///
/// The error criterion was relative when the function magnitude
/// was greater than one but absolute when it was less than one.
///
/// The following test used the relative error criterion, though
/// at certain points the relative error could be much higher than
/// indicated.
///
/// ```text
///    IEEE    -200, -4             10000     7.7e-16     1.3e-16
/// ```
#[inline]
pub(crate) fn lgam(x: f64) -> f64 {
    lgam_sgn(x).0
}

/// Sign of the Gamma function
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn gammasgn(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x > 0.0 {
        return 1.0;
    }
    if x == 0.0 {
        return 1.0_f64.copysign(x);
    }
    if x.is_infinite() {
        // x > 0 case handled, so x must be negative infinity.
        return f64::NAN;
    }
    let fx = x.floor();
    if x - fx == 0.0 {
        return f64::NAN;
    }
    // sign of gamma for x in (-n, -n+1) for positive integer n is (-1)^n.
    if (fx as i64) % 2 != 0 {
        return -1.0;
    }
    1.0
}

#[cfg(test)]
mod tests {
    /// Unlike in xsf, the sign is also correct for `x < -2^31`
    #[test]
    #[allow(clippy::float_cmp)]
    fn test_sign_large_negative() {
        // floor(x) = -2^31 - 1 is odd, so Gamma(x) < 0
        assert_eq!(super::gammasgn(-2_147_483_648.5), -1.0);
        assert_eq!(super::lgam_sgn(-2_147_483_648.5).1, -1);
        assert!(super::gamma(-2_147_483_648.5).is_sign_negative());
        // floor(x) = -2^31 - 2 is even, so Gamma(x) > 0
        assert_eq!(super::gammasgn(-2_147_483_649.5), 1.0);
        assert_eq!(super::lgam_sgn(-2_147_483_649.5).1, 1);
        assert!(super::gamma(-2_147_483_649.5).is_sign_positive());
    }
}
