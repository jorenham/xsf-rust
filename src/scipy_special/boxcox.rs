// ported from `scipy/special/_boxcox.pxd`

use crate::xsf::exp::{exp, expm1};
use crate::xsf::log::{log, log1p};

const MIN_ABS_LAMBDA: f64 = 1e-19;
const MAX_ABS_LAMBDA: f64 = 1e+273;
const MIN_ABS_LOG_X: f64 = 1e-289;
const LOG_F64_MAX: f64 = 709.78;
const SQRT_F64_MIN_POSITIVE: f64 = 1.49e-154;

/// Box-Cox transformation
///
/// The Box-Cox transformation is:
///
/// - `log(x)` if `lambda == 0.0`
/// - `(x.powf(lambda) - 1.0) / lambda` if `lambda != 0.0`
///
/// Returns [`f64::NAN`] if `x < 0.0`.
/// Returns [`f64::NEG_INFINITY`] if `x == 0.0 && lambda < 0.0`.
///
/// ## See also:
/// - [`inv_boxcox`]: Inverse of the Box-Cox transformation
/// - [`boxcox1p`]: Box-Cox transformation of 1 + `x`
#[inline]
pub fn boxcox(x: f64, lambda: f64) -> f64 {
    // if lmbda << 1 and log(x) < 1.0, the lmbda*log(x) product can lose
    // precision, furthermore, expm1(x) == x for x < eps.
    // For f64, the range of log is -744.44 to +709.78, with eps being
    // the smallest value produced. This range means that we will have
    // abs(lmbda)*log(x) < eps whenever abs(lmbda) <= eps/-log(min double)
    // which is ~2.98e-19.
    let abs_lambda = lambda.abs();
    if abs_lambda < MIN_ABS_LAMBDA {
        return log(x);
    }

    let lambda_log_x = lambda * log(x);
    if lambda_log_x < LOG_F64_MAX {
        expm1(lambda_log_x) / lambda
    } else {
        1.0_f64.copysign(lambda) * exp(lambda_log_x - log(abs_lambda)) - 1.0 / lambda
    }
}

/// Box-Cox transformation of 1 + `x`
///
/// For small `x` and `lambda`, this transformation is more numerically stable than [`boxcox`].
///
/// Returns [`f64::NAN`] if `x < -1.0`.
/// Returns [`f64::NEG_INFINITY`] if `x == -1.0 && lambda < 0.0`.
///
/// ## See also:
/// - [`inv_boxcox1p`]: Inverse of the Box-Cox transformation of 1 + `x`
/// - [`boxcox`]: Box-Cox transformation of `x`
#[doc(alias = "boxcox_1p")]
#[inline]
pub fn boxcox1p(x: f64, lambda: f64) -> f64 {
    // The argument given above in boxcox applies here with the modification
    // that the smallest value produced by log1p is the minimum representable
    // value, rather than eps.  The second condition here prevents underflow
    // when log1p(x) is < eps.
    let log_x = log1p(x);
    let abs_lambda = lambda.abs();
    if abs_lambda < MIN_ABS_LAMBDA || (abs_lambda < MAX_ABS_LAMBDA && log_x.abs() < MIN_ABS_LOG_X) {
        return log_x;
    }

    let lambda_log_x = lambda * log_x;
    if lambda_log_x < LOG_F64_MAX {
        expm1(lambda_log_x) / lambda
    } else {
        1.0_f64.copysign(lambda) * exp(lambda_log_x - log(abs_lambda)) - 1.0 / lambda
    }
}

/// Inverse of the Box-Cox transformation
///
/// Compute the compositional inverse of `y = boxcox(x, lambda)` for `x`.
///
/// ## See also:
/// - [`boxcox`]: Box-Cox transformation of `x`
/// - [`inv_boxcox1p`]: Inverse of the Box-Cox transformation of `1 + x`
#[doc(alias = "boxcox_inv")]
#[inline]
pub fn inv_boxcox(y: f64, lambda: f64) -> f64 {
    if lambda == 0.0 {
        return exp(y);
    }

    let lambda_y = lambda * y;
    let log_term = if lambda_y < f64::MAX {
        log1p(lambda_y)
    } else {
        log(1.0_f64.copysign(lambda) * (y + 1.0 / lambda)) + log(lambda.abs())
    };
    exp(log_term / lambda)
}

/// Inverse of the Box-Cox transformation of 1 + `x`
///
/// Compute the compositional inverse of `y = boxcox1p(x, lambda)` for `x`.
///
/// ## See also:
/// - [`boxcox1p`]: Box-Cox transformation of `1 + x`
/// - [`inv_boxcox`]: Inverse of the Box-Cox transformation of `x`
#[doc(alias = "boxcox_1p_inv")]
#[doc(alias = "inv_boxcox_1p")]
#[inline]
pub fn inv_boxcox1p(y: f64, lambda: f64) -> f64 {
    if lambda == 0.0 {
        return expm1(y);
    }

    let lambda_y = lambda * y;
    if lambda_y.abs() < SQRT_F64_MIN_POSITIVE {
        return y;
    }

    let log_term = if lambda_y < f64::MAX {
        log1p(lambda_y)
    } else {
        log(1.0_f64.copysign(lambda) * (y + 1.0 / lambda)) + log(lambda.abs())
    };
    expm1(log_term / lambda)
}

#[cfg(test)]
mod tests {
    use crate::testing::{np_assert_allclose, np_assert_equal};

    fn map2<X1, X2, Y, const N: usize>(f: fn(X1, X2) -> Y, x1: &[X1; N], x2: &[X2; N]) -> [Y; N]
    where
        X1: Copy,
        X2: Copy,
    {
        core::array::from_fn(|i| f(x1[i], x2[i]))
    }

    // Ported from `scipy.special.tests.test_boxcox`

    #[test]
    fn test_boxcox_basic() {
        let x = [0.5, 1.0, 2.0, 4.0];

        // lambda = 0  =>  y = log(x)
        let y = map2(crate::boxcox, &x, &[0.0; 4]);
        let expected = x.map(|xi| xi.ln());
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);

        // lambda = 1  =>  y = x - 1
        let y = map2(crate::boxcox, &x, &[1.0; 4]);
        let expected = x.map(|xi| xi - 1.0);
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);

        // lambda = 2  =>  y = 0.5*(x**2 - 1)
        let y = map2(crate::boxcox, &x, &[2.0; 4]);
        let expected = x.map(|xi| 0.5 * (xi * xi - 1.0));
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);

        // x = 0 and lambda > 0  =>  y = -1 / lambda
        let lam = [0.5, 1.0, 2.0];
        let y = map2(crate::boxcox, &[0.0; 3], &lam);
        let expected = lam.map(|l| -1.0 / l);
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);
    }

    #[test]
    fn test_boxcox_underflow() {
        let x = 1.0 + 1e-15;
        let lmbda = 1e-306;
        let y = crate::boxcox(x, lmbda);
        np_assert_allclose(&[y], &[x.ln()], 1e-14, 0.0);
    }

    #[test]
    fn test_boxcox_nonfinite() {
        // x < 0  =>  y = nan
        let x = [-1.0, -1.0, -0.5];
        let y = map2(crate::boxcox, &x, &[0.5, 2.0, -1.5]);
        np_assert_equal(&y, &[f64::NAN; 3]);

        // x = 0 and lambda <= 0  =>  y = -inf
        let x = 0.0;
        let y = map2(crate::boxcox, &[x; 2], &[-2.5, 0.0]);
        np_assert_equal(&y, &[f64::NEG_INFINITY; 2]);
    }

    #[test]
    fn test_boxcox1p_basic() {
        let x = [-0.25, -1e-20, 0.0, 1e-20, 0.25, 1.0, 3.0];

        // lambda = 0  =>  y = log(1+x)
        let y = map2(crate::boxcox1p, &x, &[0.0; 7]);
        let expected = x.map(|xi| (1.0 + xi).ln());
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);

        // lambda = 1  =>  y = x
        let y = map2(crate::boxcox1p, &x, &[1.0; 7]);
        np_assert_allclose(&y, &x, 0.0, 1.5e-7);

        // lambda = 2  =>  y = 0.5*((1+x)**2 - 1) = 0.5*x*(2 + x)
        let y = map2(crate::boxcox1p, &x, &[2.0; 7]);
        let expected = x.map(|xi| 0.5 * (xi * xi + 2.0 * xi));
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);

        // x = -1 and lambda > 0  =>  y = -1 / lambda
        let lam = [0.5, 1.0, 2.0];
        let y = map2(crate::boxcox1p, &[-1.0; 3], &lam);
        let expected = lam.map(|l| -1.0 / l);
        np_assert_allclose(&y, &expected, 0.0, 1.5e-7);
    }

    #[test]
    fn test_boxcox1p_underflow() {
        let x = [1e-15, 1e-306];
        let lmbda = [1e-306, 1e-18];
        let y = map2(crate::boxcox1p, &x, &lmbda);
        np_assert_allclose(&y, &x, 1e-14, 0.0);
    }

    #[test]
    fn test_boxcox1p_nonfinite() {
        // x < -1  =>  y = nan
        let x = [-2.0, -2.0, -1.5];
        let y = map2(crate::boxcox1p, &x, &[0.5, 2.0, -1.5]);
        np_assert_equal(&y, &[f64::NAN; 3]);

        // x = -1 and lambda <= 0  =>  y = -inf
        let x = -1.0;
        let y = map2(crate::boxcox1p, &[x; 2], &[-2.5, 0.0]);
        np_assert_equal(&y, &[f64::NEG_INFINITY; 2]);
    }

    #[test]
    fn test_inv_boxcox() {
        let x = [0.0, 1.0, 2.0];
        let lam = [0.0, 1.0, 2.0];
        let y = map2(crate::boxcox, &x, &lam);
        let x2 = map2(crate::inv_boxcox, &y, &lam);
        np_assert_allclose(&x, &x2, 0.0, 1.5e-7);

        let x = [0.0, 1.0, 2.0];
        let lam = [0.0, 1.0, 2.0];
        let y = map2(crate::boxcox1p, &x, &lam);
        let x2 = map2(crate::inv_boxcox1p, &y, &lam);
        np_assert_allclose(&x, &x2, 0.0, 1.5e-7);
    }

    #[test]
    fn test_inv_boxcox1p_underflow() {
        let x = 1e-15;
        let lam = 1e-306;
        let y = crate::inv_boxcox1p(x, lam);
        np_assert_allclose(&[y], &[x], 1e-14, 0.0);
    }

    #[test]
    fn test_boxcox_premature_overflow() {
        let x = [100.0, 155.0];
        let lmb = [0.01, -155.0];

        // test boxcox & inv_boxcox
        let y = map2(crate::boxcox, &x, &lmb);
        assert!(y.iter().all(|&v| v.is_finite()));
        let x_inv = map2(crate::inv_boxcox, &y, &lmb);
        np_assert_allclose(&x, &x_inv, 1e-7, f64::EPSILON);

        // test boxcox1p & inv_boxcox1p
        let x1m = x.map(|v| v - 1.0);
        let y1p = map2(crate::boxcox1p, &x1m, &lmb);
        assert!(y1p.iter().all(|&v| v.is_finite()));
        let x1p_inv = map2(crate::inv_boxcox1p, &y1p, &lmb);
        np_assert_allclose(&x1m, &x1p_inv, 1e-7, f64::EPSILON);
    }
}
