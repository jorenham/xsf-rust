//! Translated from Boost.Math `boost/math/special_functions/powm1.hpp` (boostorg/math@51fde38),
//! with the edge cases of `powm1_wrap` in `scipy/special/boost_special_functions.h`
//! (scipy/scipy@4fb704a).

/// $x^y - 1$
///
/// Accurate when $x^y$ is close to 1, e.g. for $x$ near 1 or small $y$. Returns NaN if $x < 0$
/// and $y$ is not a finite integer.
///
/// Corresponds to [`scipy.special.powm1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.powm1.html
///
/// # See also
/// - [`expm1`](crate::expm1) for $e^x - 1$
/// - [`log1p`](crate::log1p) for $\ln(1 + x)$
#[must_use]
pub fn powm1(x: f64, y: f64) -> f64 {
    // SciPy's `powm1_wrap`: x^0 = 1^y = 1, also for NaN or infinite x or y
    if y == 0.0 || x == 1.0 {
        return 0.0;
    }
    if x == 0.0 {
        if y < 0.0 {
            return f64::INFINITY;
        } else if y > 0.0 {
            return -1.0;
        }
    }

    // Boost's `powm1_imp_dispatch`
    if x < 0.0 {
        // y must be an integer; `fract` is NaN for infinite y, which Boost's `trunc` rejects
        if y.fract() != 0.0 {
            return f64::NAN;
        }
        if y % 2.0 == 0.0 {
            return powm1_imp(-x, y);
        }
    }
    powm1_imp(x, y)
}

/// Boost's `powm1_imp`
fn powm1_imp(x: f64, y: f64) -> f64 {
    if (y * (x - 1.0)).abs() < 0.5 || y.abs() < 0.2 {
        // We don't have any good/quick approximation for log(x) * y
        // so just try it and see:
        let l = y * x.ln();
        if l < 0.5 {
            return l.exp_m1();
        }
        // fall through....
    }

    // Unlike in SciPy, which returns 0 if x^y overflows for y < 0, the result is infinite.
    x.powf(y) - 1.0
}
