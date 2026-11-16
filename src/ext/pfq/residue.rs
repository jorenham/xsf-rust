//! The expansions of pFq in `1 / z` for `a_i` that differ by integers
//!
//! These expansions (DLMF 16.8.8 and 16.11.2) sum the residues of
//!
//! ```text
//! Gamma(-s) prod_i Gamma(a_i + s) / prod_j Gamma(b_j + s) (-z)^s
//! ```
//!
//! (DLMF 16.5.1) at its poles `s = -a_i - n`. The `a_i` that differ by integers, a class, share
//! poles of a higher order, where the terms of DLMF 16.8.8 and 16.11.2 are infinite. At
//! `s = -a_m - n + eps`, for the smallest `a_m` of a class, the integrand is `eps^-order`
//! times a power series in `eps`, whose coefficient of `eps^(order - 1)` is the residue.
//!
//! That series is `t_n u(eps)` with `u(0) = 1`, where `t_n` is the term of the expansion
//! without its zero factors. Each factor `x + n - eps` of the ratio of the integrands at `n + 1`
//! and at `n` is `(x + n) (1 - eps / (x + n))`, or `-eps` for `x + n = 0`, which changes the
//! order. So `u` follows from the one for `n = 0`, which is the exponential of the Taylor series
//! of the logarithms of the Gamma functions.

use core::f64::consts::PI;

use super::{Dd, is_pole};
use crate::utils::neg_one_pow;
use crate::xsf::cephes::psi;
use crate::zeta;

/// `u (1 + r eps)`, for a power series `u` in `eps`
fn mul_linear(u: &mut [f64], r: f64) {
    for i in (1..u.len()).rev() {
        u[i] += r * u[i - 1];
    }
}

/// `u / (1 + r eps)`, for a power series `u` in `eps`
fn div_linear(u: &mut [f64], r: f64) {
    for i in 1..u.len() {
        u[i] -= r * u[i - 1];
    }
}

/// `exp(c_1 eps + c_2 eps^2 + ...)` as a power series in `eps` up to the degree of the sum
fn exp_series(c: &[f64]) -> Vec<f64> {
    let mut u = vec![1.0];
    for k in 1..=c.len() {
        // (by u' = u c')
        let sum: f64 = (1..=k).map(|j| j as f64 * c[j - 1] * u[k - j]).sum();
        u.push(sum / k as f64);
    }
    u
}

/// The coefficient of `eps^j` of `ln Gamma(x + eps)`, or of `ln(eps Gamma(x + eps))` if `x` is a
/// pole, for `j >= 1`
fn ln_gamma_coef(j: u32, x: Dd) -> f64 {
    let (sign, power) = (neg_one_pow(j), f64::from(j));
    // psi^(j-1)(x) / j!, which is (-1)^j zeta(j, x) / j for j > 1 (DLMF 25.11.12)
    // (by `psi`, as `digamma` is much slower near its negative root)
    let taylor = |x: f64| match j {
        1 => psi(x),
        _ => sign * zeta(power, x) / power,
    };
    if x.hi > 0.0 {
        return taylor(x.hi);
    }
    // By the reflection formula (DLMF 5.5.3), Gamma(x + eps) is pi / sin(pi (x + eps)) /
    // Gamma(1 - x - eps), and the coefficient of ln sin(pi (x + eps)) is -(-1)^j / j times the
    // sum of the (d + n)^-j over the integers n, where d is the distance of x to an integer:
    // pi cot(pi d) for j = 1. A pole is without its infinite term, that of n = 0.
    //
    // That d is exact, unlike the one of x.hi, which `psi` and `zeta` would use.
    let d = f64::from(x - x.hi.round());
    let sum = match (is_pole(&x), j) {
        (true, 1) => 0.0,
        (true, _) => (1.0 + sign) * zeta(power, 1.0),
        (false, 1) => PI / (PI * d).tan(),
        (false, _) => {
            let odd = if d < 0.0 { sign } else { 1.0 };
            odd * (zeta(power, d.abs()) + sign * zeta(power, 1.0 - d.abs()))
        },
    };
    sign * (sum / power - taylor(1.0 - x.hi))
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{LN_2, PI};

    use super::*;
    use crate::utils::EGAMMA;

    const ZETA_3: f64 = 1.202_056_903_159_594_2;

    #[test]
    fn test_linear() {
        // (1 + 2 eps + 3 eps^2) (1 - 4 eps), and back
        let mut u = [1.0, 2.0, 3.0];
        mul_linear(&mut u, -4.0);
        assert_eq!(u, [1.0, -2.0, -5.0]);
        div_linear(&mut u, -4.0);
        assert_eq!(u, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_exp_series() {
        assert_eq!(exp_series(&[]), [1.0]);
        // exp(2 eps) and exp(eps^2 - eps)
        assert_eq!(exp_series(&[2.0, 0.0, 0.0]), [1.0, 2.0, 2.0, 4.0 / 3.0]);
        assert_eq!(exp_series(&[-1.0, 1.0, 0.0]), [1.0, -1.0, 1.5, -7.0 / 6.0]);
    }

    #[test]
    fn test_ln_gamma_coef() {
        let zeta_2 = PI * PI / 6.0;
        for (j, x, expected) in [
            // psi(x), psi'(x) / 2, and psi''(1) / 6
            (1, 1.0, -EGAMMA),
            (1, -0.5, 2.0 - EGAMMA - 2.0 * LN_2),
            (2, 1.0, zeta_2 / 2.0),
            (2, 3.0, (zeta_2 - 1.25) / 2.0),
            (3, 1.0, -ZETA_3 / 3.0),
            // at the poles 0 and -2: of Gamma(1 + eps), over (eps - 1) (eps - 2) for the latter
            (1, 0.0, -EGAMMA),
            (2, 0.0, zeta_2 / 2.0),
            (1, -2.0, 1.5 - EGAMMA),
            (2, -2.0, (zeta_2 + 1.25) / 2.0),
            (3, -2.0, (1.125 - ZETA_3) / 3.0),
        ] {
            let c = ln_gamma_coef(j, Dd::new(x));
            assert!((c - expected).abs() <= 4.0 * f64::EPSILON, "{j} {x}: {c}");
        }
    }

    #[test]
    fn test_ln_gamma_coef_low() {
        // where the low part of x is most of its distance to the nearest pole, or a large part
        // of it
        for (j, [x, k], expected) in [
            (1, [1e-10, 1.0], -9_999_999_999.577_215),
            (2, [1e-10, 1.0], 5e19),
            (3, [1e-10, 1.0], -3.333_333_333_333_333e29),
            (1, [-1e-10, 3.0], 10_000_000_001.256_117),
            (3, [-1e-10, 3.0], 3.333_333_333_333_333e29),
            // (the high part of x is a pole, unlike x)
            (1, [1e-20, 1.0], -1e20),
            (1, [-0.3, 40.0], 5.991_207_777_718_763),
            (2, [0.3, 40.0], 7.527_269_681_689_948),
            (3, [0.3, 40.0], -11.473_027_069_234_465),
            (1, [0.1, 1e9], 11.054_437_846_882_387),
            (4, [0.1, 1e9], 2_500.593_828_786_834_4),
        ] {
            let c = ln_gamma_coef(j, Dd::new(x) - k);
            assert!(
                (c / expected - 1.0).abs() <= 4.0 * f64::EPSILON,
                "{j} {x} {k}: {c}"
            );
        }
    }
}
