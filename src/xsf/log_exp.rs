/// Expit (a.k.a. logistic sigmoid) function, $1 / (1 + e^{-x})$
///
/// Corresponds to [`scipy.special.expit`][expit]. Translated into pure Rust from xsf.
///
/// [expit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expit.html
///
/// # See also
/// - [`logit`]
#[must_use]
#[inline]
pub fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Relative error exponential, $(e^x - 1) / x$
///
/// Corresponds to [`scipy.special.exprel`][exprel]. Translated into pure Rust from xsf.
///
/// [exprel]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exprel.html
///
/// # See also
/// - [`expm1`](crate::expm1)
#[must_use]
#[inline]
pub fn exprel(x: f64) -> f64 {
    if x.abs() < f64::EPSILON {
        1.0
    } else if x > 717.0 {
        // near f64::MAX.log()
        f64::INFINITY
    } else {
        x.exp_m1() / x
    }
}

/// Logit function, $\ln(x / (1 - x))$
///
/// Corresponds to [`scipy.special.logit`][logit]. Translated into pure Rust from xsf.
///
/// [logit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logit.html
///
/// # See also
/// - [`expit`]
#[must_use]
#[inline]
pub fn logit(x: f64) -> f64 {
    // The standard formula is log(x/(1 - x)), but this expression
    // loses precision near x=0.5, as does log(x) - log1p(-x).
    // We use the standard formula away from p=0.5, and use
    // log1p(2*(x - 0.5)) - log1p(-2*(x - 0.5)) around p=0.5, which
    // provides very good precision in this interval.
    if (0.3..=0.65).contains(&x) {
        let s = 2.0 * (x - 0.5);
        s.ln_1p() - (-s).ln_1p()
    } else {
        (x / (1.0 - x)).ln()
    }
}

/// Natural logarithm of [`expit`]
///
/// Corresponds to [`scipy.special.log_expit`][log_expit]. Translated into pure Rust from xsf.
///
/// [log_expit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_expit.html
#[must_use]
#[inline]
pub fn log_expit(x: f64) -> f64 {
    if x < 0.0 {
        x - x.exp().ln_1p()
    } else {
        -(-x).exp().ln_1p()
    }
}

/// Compute $\ln(1 - e^x)$
///
/// Translated into pure Rust from xsf.
/// This function has no analogue in `scipy.special`.
///
/// # See also
/// - [`expm1`](crate::expm1)
#[must_use]
#[inline]
pub fn log1mexp(x: f64) -> f64 {
    if x > 0.0 {
        f64::NAN
    } else if x == 0.0 {
        f64::NEG_INFINITY
    } else if x < -1.0 {
        (-x.exp()).ln_1p()
    } else {
        (-x.exp_m1()).ln()
    }
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    #[test]
    fn test_expit() {
        xsref::test("expit", "d-d", |x| crate::expit(x[0]));
    }

    #[test]
    fn test_exprel() {
        xsref::test("exprel", "d-d", |x| crate::exprel(x[0]));
    }

    #[test]
    fn test_logit() {
        xsref::test("logit", "d-d", |x| crate::logit(x[0]));
    }

    #[test]
    fn test_log_expit() {
        xsref::test("log_expit", "d-d", |x| crate::log_expit(x[0]));
    }

    #[test]
    fn test_log1mexp() {
        let xs: [f64; 5] = [-20.0, -5.0, -1.5, -1.0, -0.75];
        let expected = xs.map(|x| (-x.exp_m1()).ln());
        np_assert_allclose!(xs.map(crate::log1mexp), expected, atol = 1e-15);

        assert!(crate::log1mexp(0.5).is_nan());
        let zero_limit = crate::log1mexp(0.0);
        assert!(zero_limit.is_infinite() && zero_limit.is_sign_negative());
    }
}
