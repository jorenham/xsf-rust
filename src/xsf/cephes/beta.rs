//! Translated into pure Rust from `xsf/cephes/beta.h` (scipy/xsf@8e32110).
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf (scipy/xsf#338), the negative integer helpers take any non-positive integer `a`
//! and integer `b`, not only those that fit in an `int`.

use crate::utils::{LN_SQRT_TAU, is_nonpos_int, neg_one_pow};
use crate::xsf::cephes::consts::{MAXGAM, MAXLOG};
use crate::xsf::cephes::dd_real::two_sum;
use crate::xsf::cephes::incbet::incbet_corr;
use crate::xsf::cephes::{gamma, gammasgn, lgam, lgam_ratio, lgam_sgn, lgam1p, psi, rgamma};
use crate::xsf::trig::{cephes_sinpi, sinpi_sub};

use core::f64::consts::{PI, TAU};

const BETA_ASYMP_FACTOR: f64 = 1e6;
/// Below this `|b / a|`, `lbeta_ratio` uses `lgam_ratio`, and at or above it, for `a, b > 0`,
/// Stirling's formula
const STIRLING_RATIO: f64 = 0.125;
/// Below this `|a|` and `|b|`, `B(a, b) = (a + b) / (a b)` to within `~1.6 |a b|` (DLMF 5.7.3)
const BETA_TINY: f64 = 1e-9;

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
    if a > 0.0 && b.abs() < STIRLING_RATIO * a {
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
        // here, so b > MAXGAM / 9. Its large terms are combined into ln(q) and ln(1 - q).
        let q = stirling_q(a, b);
        let y =
            (b - 0.5) * q.ln() + a * (-q).ln_1p() - 0.5 * a.ln() + LN_SQRT_TAU - incbet_corr(a, b);
        Some((y, 1.0))
    } else {
        None
    }
}

/// `q = b / (a + b)` for Stirling's formula, from `b / a`, as `a + b` may overflow
fn stirling_q(a: f64, b: f64) -> f64 {
    let r = b / a;
    r / (1.0 + r)
}

/// `B(a, b)` by Stirling's formula, for `a >= b >= 1` and `a + b > MAXGAM`
///
/// Not as `exp(ln(B))`, whose rounding error `exp` amplifies `|ln(B)|` (up to ~745) times. With
/// `q = b / (a + b)`, as in `lbeta_ratio`:
///
/// ```text
/// B(a, b) = sqrt(2 pi / a) q^(b - 1/2) (1 - q)^a exp(-incbet_corr(a, b)).
/// ```
///
/// The rounding error `d` of `q` contributes only `~d / (2 q)`, as `q^b (1 - q)^a` is stationary
/// at `q = b / (a + b)`.
fn beta_stirling(a: f64, b: f64) -> f64 {
    let q = stirling_q(a, b);
    let p = 1.0 - q;
    // `1 - q = p + e` exactly (Sterbenz), so `(1 - q)^a = p^a (1 + a e / p)` to first order
    let e = (1.0 - p) - q;
    let pow = q.powf(b - 0.5) * p.powf(a);
    if pow == 0.0 {
        // also where `1 + a e / p < 0`, which needs `a > ~1e15`
        0.0
    } else {
        pow * (1.0 + a * e / p) * (-incbet_corr(a, b)).exp() * (TAU / a).sqrt()
    }
}

/// Whether `lbeta_ratio` reflects a negative `a`
fn lbeta_reflects_a(a: f64, b: f64) -> bool {
    a < 0.0 && b.abs() < -STIRLING_RATIO * a
}

/// Whether `beta` and `lbeta` use `lbeta_ratio` or `lgam` instead of `gamma`
fn beta_uses_lgam(a: f64, b: f64) -> bool {
    let y = a + b;
    y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM || gamma_out_of_range(a, b, y)
}

/// Whether `gamma(a)`, `gamma(b)` or `rgamma(a + b)` under- or overflows, for `|a|`, `|b|` and
/// `|a + b|` at most `MAXGAM`
///
/// Unlike in xsf (scipy/xsf#357, scipy/xsf#540), `lgam` is then used instead, also if `a` or `b`
/// is subnormal.
fn gamma_out_of_range(a: f64, b: f64, y: f64) -> bool {
    (a < 1.0 - MAXGAM
        || b < 1.0 - MAXGAM
        || y < 1.0 - MAXGAM
        || a.abs().min(b.abs()) < f64::MIN_POSITIVE)
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

/// `1 / Gamma(a + b)` for `|a + b| <= MAXGAM`
///
/// Unlike in xsf, without the rounding error `e` of `y = a + b`: to first order in `e`,
/// `1 / Gamma(y + e) = (1 - e psi(y)) / Gamma(y)`, with `psi(y) ~ ln(y - 1/2)` for `y > 1`, and
/// `y < -1/2` is reflected (DLMF 5.5.3), where `sinpi_sub` gives the exact sine, and
/// `Gamma(1 - y - e) = -y Gamma(-y) (1 - e psi(1 - y))`.
fn rgamma_sum(a: f64, b: f64) -> f64 {
    let (y, e) = two_sum(a, b);
    if y < -0.5 && y > -MAXGAM {
        sinpi_sub(a, -b) * (-y / PI) * gamma(-y) * (1.0 - e * (0.5 - y).ln())
    } else if y > 1.0 {
        rgamma(y) * (1.0 - e * (y - 0.5).ln())
    } else {
        rgamma(y)
    }
}

/// Special case for a negative integer argument
fn beta_negint(a: f64, b: f64) -> f64 {
    // `b <= -a` is `1 - a - b > 0` for integers, without its rounding error
    if b.fract() == 0.0 && b <= -a {
        let sgn = neg_one_pow(b);
        // B(1 - a - b, b) is a pole for b <= 0, which the rounding of 1 - a - b could hide
        sgn * if b > 0.0 {
            beta(1.0 + (-a - b), b)
        } else {
            f64::INFINITY
        }
    } else {
        f64::INFINITY
    }
}

fn lbeta_negint(a: f64, b: f64) -> f64 {
    if b.fract() == 0.0 && 0.0 < b && b <= -a {
        lbeta(1.0 + (-a - b), b)
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
/// If `a <= 0` and `b` are integers and `1 - a - b > 0`, the limit `(-1)^b * beta(1 - a - b, b)`
/// is returned instead (and likewise with `a` and `b` swapped).
pub(crate) fn beta(mut a: f64, mut b: f64) -> f64 {
    let mut y;
    let mut sign = 1;

    if is_nonpos_int(a) {
        if a.is_finite() {
            return beta_negint(a, b);
        }
        // overflow
        return f64::from(sign) * f64::INFINITY;
    }

    if is_nonpos_int(b) {
        if b.is_finite() {
            return beta_negint(b, a);
        }
        // overflow
        return f64::from(sign) * f64::INFINITY;
    }

    if a.abs() < b.abs() {
        y = a;
        a = b;
        b = y;
    }

    // Unlike in xsf, `(a + b) / (a b)` for a tiny `a` and `b`, where `gamma` can overflow
    if a.abs() < BETA_TINY {
        return (a + b) / a / b;
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
        // Unlike in xsf, not as `exp(ln(B))`, which amplifies the rounding error of `ln(B)`. A
        // positive `b < 1` is subnormal here, as `small_b` takes the others.
        if a > 0.0 && b >= 1.0 {
            return beta_stirling(a, b);
        }
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

    y = rgamma_sum(a, b);
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
        if a.is_finite() {
            return lbeta_negint(a, b);
        }
        // over
        return f64::from(sign) * f64::INFINITY;
    }

    if is_nonpos_int(b) {
        if b.is_finite() {
            return lbeta_negint(b, a);
        }
        // over
        return f64::from(sign) * f64::INFINITY;
    }

    if a.abs() < b.abs() {
        y = a;
        a = b;
        b = y;
    }

    // Unlike in xsf, `ln|(a + b) / (a b)|` for a tiny `a` and `b`, where `gamma` can overflow
    if a.abs() < BETA_TINY {
        // `b / a` rounds before `ln_1p` if they nearly cancel, where `a + b` is exact (Sterbenz)
        let r = b / a;
        let ln_y = if r < -0.5 {
            ((a + b) / a).abs().ln()
        } else {
            r.ln_1p()
        };
        return ln_y - b.abs().ln();
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

    y = rgamma_sum(a, b);
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

/// Below this `a`, [`ln_a_beta`] uses Simpson's rule, which errs by `~a^5 / 100`, instead of
/// `lbeta`, which loses `~|ln(a)| eps` against `ln(a)`
const LN_A_BETA_TINY: f64 = 1e-3;

/// `ln(a B(a, b))`, also where it is `O(a)`, for `a, b > 0`
///
/// For a small `a`, `a B(a, b) = (1 + a / b) Gamma(1 + a) Gamma(1 + b) / Gamma(1 + b + a)`, where
/// `ln(Gamma(1 + b + a) / Gamma(1 + b))` is the integral of `psi` over `[1 + b, 1 + b + a]`.
pub(crate) fn ln_a_beta(a: f64, b: f64) -> f64 {
    if a < LN_A_BETA_TINY {
        let c = 1.0 + b;
        let int_psi = a / 6.0 * (psi(c) + 4.0 * psi(c + 0.5 * a) + psi(c + a));
        (a / b).ln_1p() + lgam1p(a) - int_psi
    } else {
        lbeta(a, b) + a.ln()
    }
}
