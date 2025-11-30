use num_traits::Float;

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
pub fn expit<T: Float>(x: T) -> T {
    (T::one() + (-x).exp()).recip()
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
pub fn exprel<T: Float>(x: T) -> T {
    if x.abs() < T::epsilon() {
        T::one()
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
#[allow(clippy::missing_panics_doc)]
pub fn logit<T: Float>(x: T) -> T {
    // The standard formula is log(x/(1 - x)), but this expression
    // loses precision near x=0.5, as does log(x) - log1p(-x).
    // We use the standard formula away from p=0.5, and use
    // log1p(2*(x - 0.5)) - log1p(-2*(x - 0.5)) around p=0.5, which
    // provides very good precision in this interval.
    let a = T::from(0.3).unwrap();
    let b = T::from(0.65).unwrap();
    if a <= x && x <= b {
        let s = x + x - T::one();
        s.ln_1p() - (-s).ln_1p()
    } else {
        (x / (T::one() - x)).ln()
    }
}

/// Natural logarithm of [`expit`]
///
/// Corresponds to [`scipy.special.log_expit`][log_expit]. Translated into pure Rust from xsf.
///
/// [log_expit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_expit.html
#[must_use]
#[inline]
pub fn log_expit<T: Float>(x: T) -> T {
    if x.is_sign_negative() {
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
pub fn log1mexp<T: Float>(x: T) -> T {
    if x > T::zero() {
        T::nan()
    } else if x.is_zero() {
        T::neg_infinity()
    } else if x < -T::one() {
        (-x.exp()).ln_1p()
    } else {
        (-x.exp_m1()).ln()
    }
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation)]
mod tests {
    use crate::np_assert_allclose;

    #[test]
    fn test_expit_f32() {
        xsref::test("expit", "f-f", |x| crate::expit(x[0] as f32));
    }

    #[test]
    fn test_expit_f64() {
        xsref::test("expit", "d-d", |x| crate::expit(x[0]));
    }

    #[test]
    fn test_logit_f32() {
        xsref::test("logit", "f-f", |x| crate::logit(x[0] as f32));
    }

    #[test]
    fn test_logit_f64() {
        xsref::test("logit", "d-d", |x| crate::logit(x[0]));
    }

    #[test]
    fn test_log_expit_f32() {
        xsref::test("log_expit", "f-f", |x| crate::log_expit(x[0] as f32));
    }

    #[test]
    fn test_log_expit_f64() {
        xsref::test("log_expit", "d-d", |x| crate::log_expit(x[0]));
    }

    #[test]
    fn test_exprel() {
        xsref::test("exprel", "d-d", |x| crate::exprel(x[0]));
    }

    #[test]
    fn test_log1mexp() {
        let xs: [f64; 5] = [-20.0, -5.0, -1.5, -1.0, -0.75];
        let expected = xs.map(|x| (-x.exp_m1()).ln());
        np_assert_allclose!(xs.map(crate::log1mexp), expected, atol = 1e-15);

        assert!(crate::log1mexp(0.5_f64).is_nan());

        let zero_limit = crate::log1mexp(0.0_f64);
        assert!(zero_limit.is_infinite() && zero_limit.is_sign_negative());
    }
}
