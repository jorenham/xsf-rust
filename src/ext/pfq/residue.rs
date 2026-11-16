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
use core::iter::{once, repeat};

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

/// The residues of the class of `a_m` at `s = -a_m - n` for `n = 0, 1, ...`, without the factor
/// of the one at `n = 0` if that were a simple pole
///
/// Each comes with a bound on its size: `|t_n|` times the largest of the bounds on the
/// coefficients of `u` up to its own. They end before a `t_n` that is not normal.
pub(super) struct Residues {
    /// The parameters `x` of the factors `x + n - eps` of the numerator and of the denominator
    /// of the ratio of two integrands, in double-double as `x + n` is exact where it vanishes
    up: Vec<Dd>,
    low: Vec<Dd>,
    /// The other factor of that ratio: `+-1 / z`
    w: f64,
    n: f64,
    /// `t_n`
    term: f64,
    order: i32,
    /// `u`, and bounds on the sizes of its coefficients
    u: Vec<f64>,
    size: Vec<f64>,
}

impl Residues {
    /// For the class of `a_m`, where `b` is without the 1 of `n!`, and no `a_i` is a non-positive
    /// integer
    ///
    /// For `z > 0`, where `ln(-z)` is `ln|z| -+ pi i`, only the bounds are valid.
    pub(super) fn of_class(a: &[f64], b: &[f64], m: usize, z: f64) -> Self {
        let am = Dd::new(a[m]);
        let diffs = a.iter().map(|&ai| Dd::new(ai) - am);
        let degree = diffs.clone().filter(|d| d.is_integer()).count() - 1;
        // The logarithm of the integrand at n = 0 without its poles and zeros, as its
        // coefficients of eps, eps^2, ..., and bounds on the sizes of their parts, from
        // (-z)^eps on.
        let ln_z = z.abs().ln();
        let pi = if z > 0.0 { PI } else { 0.0 };
        let first = |v: f64| once(v).chain(repeat(0.0)).take(degree).collect::<Vec<_>>();
        let (mut coef, mut size, mut order) = (first(ln_z), first(ln_z.abs() + pi), 0);
        // Gamma(-s), the Gamma(a_i + s), and the 1 / Gamma(b_j + s), as the x, power and sign of
        // Gamma(x + sign eps)^power
        let gammas = once((am, 1, -1.0))
            .chain(diffs.map(|d| (d, 1, 1.0)))
            .chain(b.iter().map(|&bj| (Dd::new(bj) - am, -1, 1.0)));
        for (x, power, sign) in gammas {
            order += if is_pole(&x) { power } else { 0 };
            let mut signs = f64::from(power);
            for (j, (coef, size)) in (1..).zip(coef.iter_mut().zip(&mut size)) {
                signs *= sign;
                let c = signs * ln_gamma_coef(j, x);
                (*coef, *size) = (*coef + c, *size + c.abs());
            }
        }
        // a_m and the 1 + a_m - b_j, over the 1 + a_m - a_i
        let up = once(am).chain(b.iter().map(|&bj| am + 1.0 - bj)).collect();
        let low = a.iter().map(|&ai| am + 1.0 - ai).collect();
        // (-1)^(p+q+1) / z
        let w = neg_one_pow(a.len() + b.len() + 1) / z;
        let (u, size) = (exp_series(&coef), exp_series(&size));
        let (n, term) = (0.0, 1.0);
        Self {
            up,
            low,
            w,
            n,
            term,
            order,
            u,
            size,
        }
    }

    /// The parameters `up` and `low` as `f64`
    pub(super) fn params(&self) -> [Vec<f64>; 2] {
        [&self.up, &self.low].map(|v| v.iter().copied().map(f64::from).collect())
    }

    /// Bound on the ratio of the bounds on two successive residues from term `n` on, relative
    /// to that on the ratio of their `t_n`, for the parameters `up` and `low` as `f64`, which
    /// plus `n` must exceed 1
    ///
    /// It is the product of the sums of the sizes of the coefficients of `1 - eps / (x + n)`
    /// and of `1 / (1 - eps / (x + n))`.
    pub(super) fn growth_bound(up: &[f64], low: &[f64], n: f64) -> f64 {
        let up = up.iter().map(|x| 1.0 + 1.0 / (x + n));
        let low = low.iter().map(|x| (x + n) / (x + n - 1.0));
        up.chain(low).product()
    }
}

impl Iterator for Residues {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<(f64, f64)> {
        if !self.term.is_normal() {
            return None;
        }
        // (there is no residue for an order below 1)
        let k = usize::try_from(self.order).unwrap_or(0);
        let size = self.size[..k].iter().copied().fold(0.0, f64::max);
        let weight = if k > 0 { self.u[k - 1] } else { 0.0 };
        let item = (self.term * weight, self.term.abs() * size);
        for &x in &self.up {
            let c = f64::from(x + self.n);
            if c == 0.0 {
                (self.term, self.order) = (-self.term, self.order - 1);
            } else {
                self.term *= c;
                mul_linear(&mut self.u, -1.0 / c);
                mul_linear(&mut self.size, 1.0 / c.abs());
            }
        }
        for &x in &self.low {
            let c = f64::from(x + self.n);
            if c == 0.0 {
                (self.term, self.order) = (-self.term, self.order + 1);
            } else {
                self.term /= c;
                div_linear(&mut self.u, -1.0 / c);
                div_linear(&mut self.size, -1.0 / c.abs());
            }
        }
        self.term *= self.w;
        self.n += 1.0;
        Some(item)
    }
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

    /// Whether `y` is within a few ulp of an `expected` value of up to about 10
    fn close(y: f64, expected: f64) -> bool {
        (y - expected).abs() <= 1e-14
    }

    #[test]
    fn test_of_class() {
        // 2F1(1, 1; 2; z) = -ln(1 - z) / z = (ln(-z) + ln(1 - 1 / z)) / -z: a pole of order 2
        let residues = Residues::of_class(&[1.0; 2], &[2.0], 0, -3.0);
        assert_eq!(residues.order, 2);
        assert!(close(residues.u[1], 3.0_f64.ln()));
        // 3F2(1, 1, 1; 2, 2; z) = Li_2(z) / z, which is
        // (ln(-z)^2 / 2 + pi^2 / 6 + Li_2(1 / z)) / -z: a pole of order 3
        let (l, zeta_2) = (40.0_f64.ln(), PI * PI / 6.0);
        let residues = Residues::of_class(&[1.0; 3], &[2.0; 2], 0, -40.0);
        assert_eq!(residues.order, 3);
        assert!(close(residues.u[2], l * l / 2.0 + zeta_2));
        // for z > 0, with ln(z) -+ pi i, the bound holds for (ln(z)^2 - pi^2) / 2 + pi^2 / 6
        let residues = Residues::of_class(&[1.0; 3], &[2.0; 2], 0, 40.0);
        assert!(residues.size[2] >= (l * l + PI * PI) / 2.0 + zeta_2);
    }

    #[test]
    fn test_residues() {
        // no other a_i differs from a_m by an integer: the term of a_1 = 1 of
        // 2F1(1, 1.5; 2; z), a series with the parameters 1 and 0 over 0.5
        let mut residues = Residues::of_class(&[1.0, 1.5], &[2.0], 0, -3.0);
        assert_eq!(residues.next(), Some((1.0, 1.0)));
        assert_eq!(residues.next(), Some((0.0, 0.0)));
        // 2F1(1, 1; 2; z): ln(-z), and the terms of ln(1 - 1 / z) from the simple poles after
        // the zero factor of b_1
        let mut residues = Residues::of_class(&[1.0; 2], &[2.0], 0, -3.0);
        for expected in [3.0_f64.ln(), 1.0 / 3.0, -1.0 / 18.0] {
            let (t, size) = residues.next().unwrap();
            assert!(close(t, expected) && size >= t.abs(), "{t} {size}");
        }
    }
}
