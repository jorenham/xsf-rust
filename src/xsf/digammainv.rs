//! Translated into pure Rust from `xsf/digammainv.h` (xsf main).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::utils::{MAXLOG, SCIPY_EULER};
use crate::xsf::cephes::{psi, zeta};
use crate::xsf::tools::{NewtonRootFinderStatus, find_root_newton};

/// Initial guess for `digammainv` using Minka (2000):
///   x0 = exp(y) + 0.5      if y >= -2.22
///   x0 = -1 / (y + gamma)  if y < -2.22
/// where gamma = 0.5772... is the Euler-Mascheroni constant.
fn digammainv_initial_guess(y: f64) -> f64 {
    if y >= -2.22 {
        return y.exp() + 0.5;
    }
    -1.0 / (y + SCIPY_EULER)
}

/// Refine an initial guess x using Newton-Raphson iteration:
///   x_new = x_old - (psi(x_old) - y) / psi'(x_old)
/// where psi'(x) = zeta(2, x) (the trigamma function).
fn digammainv_newton_raphson(x: f64, y: f64) -> f64 {
    let (res, status) = find_root_newton(
        |x| (psi(x) - y, zeta(2.0, x)),
        x,
        4.0 * f64::EPSILON,
        0.0,
        100,
    );
    match status {
        // Unlike in xsf (scipy/xsf#361), `MaxIterationsExceeded` returns the last iterate instead
        // of NaN, e.g. `digammainv(9.2)` was NaN instead of 9898.
        NewtonRootFinderStatus::DerivativeZero | NewtonRootFinderStatus::ObjectiveReturnedNan => {
            f64::NAN
        },
        _ => res,
    }
}

/// Inverse of the digamma function
///
/// Computes $x > 0$ such that $\psi(x) = y$, where $\psi$ is the [`digamma`](crate::digamma)
/// function, using Newton-Raphson iteration with an initial guess from Minka (2000) \[1\].
///
/// Corresponds to `xsf::digammainv` in xsf main.
///
/// \[1\] T. Minka, "Estimating a Dirichlet distribution", 2000. Appendix C
///     <https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf>
///
/// # See also
/// - [`digamma`](crate::digamma): The digamma function, $\psi(z)$
#[must_use]
#[inline]
pub fn digammainv(y: f64) -> f64 {
    if y.is_nan() {
        return f64::NAN;
    }
    // digamma(x) app. log(x) for large x, so digammainv(y) must be inf for y > log(DBL_MAX)
    if y > MAXLOG {
        return f64::INFINITY;
    }
    // digamma(x) -> - inf for x -> 0+, so digammainv must be zero for y = - inf
    if y.is_infinite() && y < 0.0 {
        return 0.0;
    }
    let x = digammainv_initial_guess(y);
    digammainv_newton_raphson(x, y)
}
