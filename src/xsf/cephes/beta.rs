//! Translated into pure Rust from `xsf/cephes/beta.h` (scipy/xsf@8e32110).
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `1 - a` in the negative integer helpers is computed in floating point instead of
//! in `int`, where it overflows (UB) for `a = -2^31` and `a = -2^31 + 1`.

use crate::utils::{LN_SQRT_2PI, MAXGAM, MAXLOG, is_nonpos_int};
use crate::xsf::cephes::incbet::incbet_corr;
use crate::xsf::cephes::{gamma, gammasgn, lgam, lgam_ratio, lgam_sgn, rgamma};
use crate::xsf::trig::{cephes_sinpi, sinpi_sub};

use core::f64::consts::PI;

const BETA_ASYMP_FACTOR: f64 = 1e6;

/// Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
///
/// Returns the sign as well, instead of through the `int *sgn` out-parameter.
fn lbeta_asymp(a: f64, b: f64) -> (f64, i32) {
    let (mut r, sgn) = lgam_sgn(b);
    r -= b * a.ln();

    /*
     * Evaluate the correction terms through the small ratios b/a and (1-b)/a
     * instead of the products b*(1-b), b*b*(1-b)*(1-b) and a*a*a, which overflow
     * for |b| > ~1e77 or a > ~5.6e102 and turned the result into NaN (inf/inf)
     * for e.g. lbeta(1e100, 1e107) although log|B| is an ordinary finite number.
     */
    let t = b / a;
    let u = (1.0 - b) / a;

    r += t * (1.0 - b) / 2.0;
    r += t * u * (1.0 - 2.0 * b) / 12.0;
    r += -t * t * u * (1.0 - b) / 12.0;

    (r, sgn)
}

/// `ln(|B(a, b)|)` and the sign of `B(a, b)` for `|b| < |a| / 8`, or `a > 0` and `b <= -a / 2`,
/// or `a, b > 0`, or `None` otherwise
///
/// Unlike in xsf (scipy/scipy#26423, scipy/xsf#357), `lgam_ratio` or Stirling's formula avoids
/// the cancellation in `lgam(a) + lgam(b) - lgam(a + b)`. The negative one of the two largest
/// arguments is reflected.
fn lbeta_ratio(a: f64, b: f64) -> Option<(f64, f64)> {
    if a > 0.0 && b.abs() < 0.125 * a {
        let (y, sign) = lgam_sgn(b);
        Some((y - lgam_ratio(a, b), f64::from(sign)))
    } else if lbeta_reflects_a(a, b) {
        // Gamma(a) / Gamma(a + b) = sin(pi (a + b)) / sin(pi a) Gamma(1 - a - b) / Gamma(1 - a)
        // `sinpi_sub` avoids the rounding error of a + b
        let (y, sign) = lgam_sgn(b);
        let s = sinpi_sub(a, -b) / cephes_sinpi(a);
        let y = y + s.abs().ln() - lgam_ratio(1.0 - a - b, b);
        Some((y, f64::from(if s < 0.0 { -sign } else { sign })))
    } else if a > 0.0 && b <= -0.5 * a {
        // B(a, b) = pi / sin(pi b) Gamma(a) / Gamma(1 - b) / Gamma(a + b), with exact a + b
        let y = (PI / cephes_sinpi(b).abs()).ln() - lgam_ratio(a, 1.0 - (a + b)) - lgam(a + b);
        Some((y, gammasgn(b)))
    } else if a > 0.0 && b > 0.0 {
        // Stirling's formula for all three Gamma functions, with b >= a / 8 > 0 and a + b > MAXGAM
        // here, so b > MAXGAM / 9. Its large terms are combined into ln(q) and ln(1 - q), with
        // q = b / (a + b) from b / a, as a + b may overflow.
        let r = b / a;
        let q = r / (1.0 + r);
        let y =
            (b - 0.5) * q.ln() + a * (-q).ln_1p() - 0.5 * a.ln() + LN_SQRT_2PI - incbet_corr(a, b);
        Some((y, 1.0))
    } else {
        None
    }
}

/// Whether `beta(a, b)` and `lbeta(a, b)` reflect a negative `a` in `lbeta_ratio`
pub(crate) fn beta_reflects_a(a: f64, b: f64) -> bool {
    !is_nonpos_int(a) && !is_nonpos_int(b) && lbeta_reflects_a(a, b) && beta_uses_lgam(a, b)
}

/// Whether `lbeta_ratio` reflects a negative `a`
fn lbeta_reflects_a(a: f64, b: f64) -> bool {
    a < 0.0 && b.abs() < -0.125 * a
}

/// Whether `beta` and `lbeta` use `lbeta_ratio` or `lgam` instead of `gamma`
fn beta_uses_lgam(a: f64, b: f64) -> bool {
    let y = a + b;
    y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM || gamma_out_of_range(a, b, y)
}

/// Whether `gamma(a)`, `gamma(b)` or `rgamma(a + b)` under- or overflows, for `|a|`, `|b|` and
/// `|a + b|` at most `MAXGAM`
///
/// Unlike in xsf (scipy/xsf#357), `lgam` is then used instead.
fn gamma_out_of_range(a: f64, b: f64, y: f64) -> bool {
    (a < 1.0 - MAXGAM || b < 1.0 - MAXGAM || y < 1.0 - MAXGAM)
        && !(gamma(a).is_normal() && gamma(b).is_normal() && rgamma(y).is_normal())
}

/// Whether `|b| < 1` and `4 < a < inf`, where `Gamma(a) / Gamma(a + b)` is computed as
/// `exp(-lgam_ratio(a, b))`
///
/// Unlike in xsf, the rounding error of `a + b` is then not amplified ~`a ln(a)` times by
/// `Gamma(a + b)`, e.g. in `beta(1e-10, 171)`. For `a <= 4` that rounding error is (measured)
/// smaller than the error of `lgam_ratio`.
fn small_b(a: f64, b: f64) -> bool {
    b.abs() < 1.0 && a > 4.0 && a.is_finite()
}

/// Special case for a negative integer argument
fn beta_negint(a: i32, b: f64) -> f64 {
    // Like the C++ `static_cast<int>(b)` wherever that is defined; the saturating cast makes the
    // comparison false otherwise.
    if b == f64::from(b as i32) && 1.0 - f64::from(a) - b > 0.0 {
        let sgn = if (b as i32) % 2 == 0 { 1 } else { -1 };
        f64::from(sgn) * beta(1.0 - f64::from(a) - b, b)
    } else {
        f64::INFINITY
    }
}

fn lbeta_negint(a: i32, b: f64) -> f64 {
    if b == f64::from(b as i32) && 1.0 - f64::from(a) - b > 0.0 {
        lbeta(1.0 - f64::from(a) - b, b)
    } else {
        f64::INFINITY
    }
}

/// Beta function
///
/// ```text
///                   -     -
///                  | (a) | (b)
/// beta( a, b )  =  -----------.
///                     -
///                    | (a+b)
/// ```
///
/// For large arguments the logarithm of the function is
/// evaluated using lgam(), then exponentiated.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE       0,30       30000       1.6e-14     6.3e-15
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition          value returned
/// beta overflow    log(beta) > MAXLOG       +-INFINITY
///                  a or b <=0 integer       INFINITY
/// ```
///
/// If `a <= 0` and `b` are integers that fit in an `int` (32 bits) and `1 - a - b > 0`, the limit
/// `(-1)^b * beta(1 - a - b, b)` is returned instead (and likewise with `a` and `b` swapped).
pub(crate) fn beta(mut a: f64, mut b: f64) -> f64 {
    let mut y;
    let mut sign = 1;

    if is_nonpos_int(a) {
        if a == f64::from(a as i32) {
            return beta_negint(a as i32, b);
        }
        // overflow
        return f64::from(sign) * f64::INFINITY;
    }

    if is_nonpos_int(b) {
        if b == f64::from(b as i32) {
            return beta_negint(b as i32, a);
        }
        // overflow
        return f64::from(sign) * f64::INFINITY;
    }

    if a.abs() < b.abs() {
        y = a;
        a = b;
        b = y;
    }

    if small_b(a, b) {
        return gamma(b) * (-lgam_ratio(a, b)).exp();
    }

    if a.abs() > BETA_ASYMP_FACTOR * b.abs() && a > BETA_ASYMP_FACTOR {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        (y, sign) = lbeta_asymp(a, b);
        return f64::from(sign) * y.exp();
    }

    y = a + b;
    if beta_uses_lgam(a, b) {
        if let Some((y, sign)) = lbeta_ratio(a, b) {
            return sign * y.exp();
        }
        let mut sgngam;
        (y, sgngam) = lgam_sgn(y);
        sign *= sgngam; /* keep track of the sign */
        let lgam_b;
        (lgam_b, sgngam) = lgam_sgn(b);
        y = lgam_b - y;
        sign *= sgngam;
        let lgam_a;
        (lgam_a, sgngam) = lgam_sgn(a);
        y += lgam_a;
        sign *= sgngam;
        if y > MAXLOG {
            // overflow
            return f64::from(sign) * f64::INFINITY;
        }
        return f64::from(sign) * y.exp();
    }

    y = rgamma(y);
    a = gamma(a);
    b = gamma(b);
    if y.is_infinite() {
        // overflow
        return f64::from(sign) * f64::INFINITY;
    }

    if ((a * y).abs() - 1.0).abs() > ((b * y).abs() - 1.0).abs() {
        y *= b;
        y *= a;
    } else {
        y *= a;
        y *= b;
    }

    y
}

/// Natural log of |beta|.
pub(crate) fn lbeta(mut a: f64, mut b: f64) -> f64 {
    let mut y;
    let sign = 1;

    if is_nonpos_int(a) {
        if a == f64::from(a as i32) {
            return lbeta_negint(a as i32, b);
        }
        // over
        return f64::from(sign) * f64::INFINITY;
    }

    if is_nonpos_int(b) {
        if b == f64::from(b as i32) {
            return lbeta_negint(b as i32, a);
        }
        // over
        return f64::from(sign) * f64::INFINITY;
    }

    if a.abs() < b.abs() {
        y = a;
        a = b;
        b = y;
    }

    if a.abs() > BETA_ASYMP_FACTOR * b.abs() && a > BETA_ASYMP_FACTOR {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        (y, _) = lbeta_asymp(a, b);
        return y;
    }

    if small_b(a, b) {
        return lgam(b) - lgam_ratio(a, b);
    }

    y = a + b;
    if beta_uses_lgam(a, b) {
        if let Some((y, _)) = lbeta_ratio(a, b) {
            return y;
        }
        (y, _) = lgam_sgn(y);
        y = lgam_sgn(b).0 - y;
        y += lgam_sgn(a).0;
        return y;
    }

    y = rgamma(y);
    a = gamma(a);
    b = gamma(b);
    if y.is_infinite() {
        // over
        return f64::from(sign) * f64::INFINITY;
    }

    if ((a * y).abs() - 1.0).abs() > ((b * y).abs() - 1.0).abs() {
        y *= b;
        y *= a;
    } else {
        y *= a;
        y *= b;
    }

    if y < 0.0 {
        y = -y;
    }

    y.ln()
}
