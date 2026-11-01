// ported from `scipy/special/_boxcox.pxd`

const MIN_ABS_LAMBDA: f64 = 1e-19;
const MAX_ABS_LAMBDA: f64 = 1e+273;
const MIN_ABS_LOG_X: f64 = 1e-289;
const LOG_F64_MAX: f64 = 709.78;
const SQRT_F64_MIN_POSITIVE: f64 = 1.49e-154;

/// Box-Cox transformation
///
/// The Box-Cox transformation is:
///
/// - `x.ln()` if `lambda == 0.0`
/// - `(x.powf(lambda) - 1.0) / lambda` if `lambda != 0.0`
///
/// Returns [`f64::NAN`] if `x < 0.0`.
/// Returns [`f64::NEG_INFINITY`] if `x == 0.0 && lambda < 0.0`.
///
/// Corresponds to [`scipy.special.boxcox`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.boxcox.html
///
/// # See also
/// - [`inv_boxcox`]: Inverse of the Box-Cox transformation
/// - [`boxcox1p`]: Box-Cox transformation of 1 + `x`
#[must_use]
#[inline]
pub fn boxcox(x: f64, lambda: f64) -> f64 {
    // if lmbda << 1 and x.ln() < 1.0, the lmbda*x.ln() product can lose
    // precision, furthermore, expm1(x) == x for x < eps.
    // For f64, the range of log is -744.44 to +709.78, with eps being
    // the smallest value produced. This range means that we will have
    // abs(lmbda)*x.ln() < eps whenever abs(lmbda) <= eps/-log(min double)
    // which is ~2.98e-19.
    let abs_lambda = lambda.abs();
    if abs_lambda < MIN_ABS_LAMBDA {
        return x.ln();
    }

    let lambda_log_x = lambda * x.ln();
    if lambda_log_x < LOG_F64_MAX {
        lambda_log_x.exp_m1() / lambda
    } else {
        1.0_f64.copysign(lambda) * (lambda_log_x - abs_lambda.ln()).exp() - 1.0 / lambda
    }
}

/// Box-Cox transformation of 1 + `x`
///
/// For small `x` and `lambda`, this transformation is more numerically stable than [`boxcox`].
///
/// Returns [`f64::NAN`] if `x < -1.0`.
/// Returns [`f64::NEG_INFINITY`] if `x == -1.0 && lambda < 0.0`.
///
/// Corresponds to [`scipy.special.boxcox1p`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.boxcox1p.html
///
/// # See also
/// - [`inv_boxcox1p`]: Inverse of the Box-Cox transformation of 1 + `x`
/// - [`boxcox`]: Box-Cox transformation of `x`
#[doc(alias = "boxcox_1p")]
#[must_use]
#[inline]
pub fn boxcox1p(x: f64, lambda: f64) -> f64 {
    // The argument given above in boxcox applies here with the modification
    // that the smallest value produced by log1p is the minimum representable
    // value, rather than eps.  The second condition here prevents underflow
    // when log1p(x) is < eps.
    let log_x = x.ln_1p();
    let abs_lambda = lambda.abs();
    if abs_lambda < MIN_ABS_LAMBDA || (abs_lambda < MAX_ABS_LAMBDA && log_x.abs() < MIN_ABS_LOG_X) {
        return log_x;
    }

    let lambda_log_x = lambda * log_x;
    if lambda_log_x < LOG_F64_MAX {
        lambda_log_x.exp_m1() / lambda
    } else {
        1.0_f64.copysign(lambda) * (lambda_log_x - abs_lambda.ln()).exp() - 1.0 / lambda
    }
}

/// Inverse of the Box-Cox transformation
///
/// Compute the compositional inverse of `y = boxcox(x, lambda)` for `x`.
///
/// Corresponds to [`scipy.special.inv_boxcox`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.inv_boxcox.html
///
/// # See also
/// - [`boxcox`]: Box-Cox transformation of `x`
/// - [`inv_boxcox1p`]: Inverse of the Box-Cox transformation of `1 + x`
#[doc(alias = "boxcox_inv")]
#[must_use]
#[inline]
pub fn inv_boxcox(y: f64, lambda: f64) -> f64 {
    if lambda == 0.0 {
        return y.exp();
    }

    let lambda_y = lambda * y;
    let log_term = if lambda_y < f64::MAX {
        lambda_y.ln_1p()
    } else {
        (1.0_f64.copysign(lambda) * (y + 1.0 / lambda)).ln() + lambda.abs().ln()
    };
    (log_term / lambda).exp()
}

/// Inverse of the Box-Cox transformation of 1 + `x`
///
/// Compute the compositional inverse of `y = boxcox1p(x, lambda)` for `x`.
///
/// Corresponds to [`scipy.special.inv_boxcox1p`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.inv_boxcox1p.html
///
/// # See also
/// - [`boxcox1p`]: Box-Cox transformation of `1 + x`
/// - [`inv_boxcox`]: Inverse of the Box-Cox transformation of `x`
#[doc(alias = "boxcox_1p_inv", alias = "inv_boxcox_1p")]
#[must_use]
#[inline]
pub fn inv_boxcox1p(y: f64, lambda: f64) -> f64 {
    if lambda == 0.0 {
        return y.exp_m1();
    }

    let lambda_y = lambda * y;
    if lambda_y.abs() < SQRT_F64_MIN_POSITIVE {
        return y.exp_m1();
    }

    let log_term = if lambda_y < f64::MAX {
        lambda_y.ln_1p()
    } else {
        (1.0_f64.copysign(lambda) * (y + 1.0 / lambda)).ln() + lambda.abs().ln()
    };
    (log_term / lambda).exp_m1()
}
