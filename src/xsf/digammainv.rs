//! Translated into pure Rust from `xsf/digammainv.h` (xsf main).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::cephes::{psi, zeta};
use crate::xsf::tools::{NewtonRootFinderStatus, find_root_newton};

const SCIPY_EULER: f64 = 0.577_215_664_901_532_9;
const MAXLOG: f64 = 709.782_712_893_384; // log(DBL_MAX)

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
        // of NaN. Near the root, a 1 ulp change of `psi(x) - y` moves `x` by up to ~`y` eps
        // relative (because `x psi'(x) > 1`), so for `y > 4` Newton can end up cycling between two
        // values that are both accurate, but too far apart for `rtol = 4 eps`, e.g. for
        // `digammainv(9.220027678716601)`.
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

#[cfg(test)]
mod tests {
    #[test]
    fn test_digammainv() {
        // mpmath: x such that digamma(x) == y
        let ys = [0.0, 1.0, -1.0, -10.0, 5.0, 700.0];
        let xs = [
            1.461_632_144_968_362_2,
            3.203_171_468_376_931_2,
            0.785_003_325_374_520_7,
            0.104_357_198_770_116_51,
            148.912_878_356_218_87,
            1.014_232_054_735_004_5e304,
        ];
        xsref::np_assert_allclose!(ys.map(crate::digammainv), xs, rtol = 1e-14);

        assert!(crate::digammainv(f64::NAN).is_nan());
        assert_eq!(crate::digammainv(f64::NEG_INFINITY), 0.0);
        assert_eq!(crate::digammainv(710.0), f64::INFINITY);
    }

    #[test]
    fn test_digammainv_newton_cycle() {
        // unlike in xsf (scipy/xsf#361): Newton cycles here, which used to return NaN
        let x = crate::digammainv(9.220_027_678_716_601);
        xsref::np_assert_allclose!([x], [10_097.843_801_671_605], rtol = 1e-14);
    }
}
