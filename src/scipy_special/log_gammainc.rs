//! The edge cases follow `lgamma_p_wrap` and `lgamma_q_wrap` in
//! `scipy/special/boost_special_functions.h`. Unlike SciPy, which uses Boost's `lgamma_p` and
//! `lgamma_q`, these use the Cephes `igam` and `igamc`, with the same split as Boost: the log of
//! the complement where it is small, and a log-space series where the function underflows.

use crate::xsf::cephes::{igam, igam_exp_scaled, igamc, igamc_exp_scaled};

/// Log of the regularized lower incomplete gamma function, $\ln P(a,x)$
///
/// Accurate also where $P(a,x)$ underflows, or is close to 1.
///
/// Corresponds to [`scipy.special.log_gammainc`][scipy].
///
/// [scipy]: https://scipy.github.io/devdocs/reference/generated/scipy.special.log_gammainc.html
///
/// # See also
/// - [`gammainc`](crate::gammainc): Regularized lower incomplete gamma function, $P(a,x)$
/// - [`log_gammaincc`]: Log of the regularized upper incomplete gamma function, $\ln Q(a,x)$
#[must_use]
#[inline]
pub fn log_gammainc(a: f64, x: f64) -> f64 {
    ln_igam(a, x, igam, igamc, igam_exp_scaled)
}

/// Log of the regularized upper incomplete gamma function, $\ln Q(a,x)$
///
/// Accurate also where $Q(a,x)$ underflows, or is close to 1.
///
/// Corresponds to [`scipy.special.log_gammaincc`][scipy].
///
/// [scipy]: https://scipy.github.io/devdocs/reference/generated/scipy.special.log_gammaincc.html
///
/// # See also
/// - [`gammaincc`](crate::gammaincc): Regularized upper incomplete gamma function, $Q(a,x)$
/// - [`log_gammainc`]: Log of the regularized lower incomplete gamma function, $\ln P(a,x)$
#[must_use]
#[inline]
pub fn log_gammaincc(a: f64, x: f64) -> f64 {
    ln_igam(a, x, igamc, igam, |a, x| Some(igamc_exp_scaled(a, x)))
}

/// `ln(f(a, x))`, with `f` either `igam` or `igamc`, `complement` the other one, and `exp_scaled`
/// as `igam_exp_scaled` for `f`
fn ln_igam(
    a: f64,
    x: f64,
    f: fn(f64, f64) -> f64,
    complement: fn(f64, f64) -> f64,
    exp_scaled: fn(f64, f64) -> Option<(f64, f64)>,
) -> f64 {
    // unlike `igam` and `igamc`, NaN for `a = 0`
    if !(a > 0.0) {
        return f64::NAN;
    }
    let v = f(a, x);
    // exact for `x = 0`, `x = inf` or `a = inf`, and NaN for `x < 0`
    if !(x > 0.0) || x.is_infinite() || a.is_infinite() {
        return v.ln();
    }
    // `ln(v)` loses `~eps / (1 - v)` relative to the complement
    if v > 0.5 {
        return (-complement(a, x)).ln_1p();
    }
    if v >= f64::MIN_POSITIVE {
        return v.ln();
    }
    exp_scaled(a, x).map_or(v.ln(), |(m, ln_f)| m.ln() + ln_f)
}
