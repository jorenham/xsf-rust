//! Translated into pure Rust from `xsf/cephes/beta.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `1 - a` in the negative integer helpers is computed in floating point instead of
//! in `int`, where it overflows (UB) for `a = -2^31` and `a = -2^31 + 1`.

use crate::xsf::cephes::{gamma, lgam_sgn, rgamma};

// cephes/const.h

/// log(DBL_MAX)
const MAXLOG: f64 = 709.782_712_893_384;
/// Largest x such that Gamma(x) is finite
const MAXGAM: f64 = 171.624_376_956_302_7;

const BETA_ASYMP_FACTOR: f64 = 1e6;

/// Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
///
/// Returns the sign as well, instead of through the `int *sgn` out-parameter.
fn lbeta_asymp(a: f64, b: f64) -> (f64, i32) {
    let (mut r, sgn) = lgam_sgn(b);
    r -= b * a.ln();

    r += b * (1.0 - b) / (2.0 * a);
    r += b * (1.0 - b) * (1.0 - 2.0 * b) / (12.0 * a * a);
    r += -b * b * (1.0 - b) * (1.0 - b) / (12.0 * a * a * a);

    (r, sgn)
}

/// Special case for a negative integer argument
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
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

#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
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
/// (Unlike in xsf, where the value returned is documented as 0.0.) If `a <= 0` and `b` are
/// integers that fit in an `int` (32 bits) and `1 - a - b > 0`, the limit
/// `(-1)^b * beta(1 - a - b, b)` is returned instead (and likewise with `a` and `b` swapped).
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn beta(mut a: f64, mut b: f64) -> f64 {
    let mut y;
    let mut sign = 1;

    if a <= 0.0 && a == a.floor() {
        if a == f64::from(a as i32) {
            return beta_negint(a as i32, b);
        }
        // overflow
        return f64::from(sign) * f64::INFINITY;
    }

    if b <= 0.0 && b == b.floor() {
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

    if a.abs() > BETA_ASYMP_FACTOR * b.abs() && a > BETA_ASYMP_FACTOR {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        (y, sign) = lbeta_asymp(a, b);
        return f64::from(sign) * y.exp();
    }

    y = a + b;
    if y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM {
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
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn lbeta(mut a: f64, mut b: f64) -> f64 {
    let mut y;
    let sign = 1;

    if a <= 0.0 && a == a.floor() {
        if a == f64::from(a as i32) {
            return lbeta_negint(a as i32, b);
        }
        // over
        return f64::from(sign) * f64::INFINITY;
    }

    if b <= 0.0 && b == b.floor() {
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

    y = a + b;
    if y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM {
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
