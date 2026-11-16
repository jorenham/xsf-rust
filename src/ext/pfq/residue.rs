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

use core::iter::{once, repeat, zip};

use super::Dd;
use crate::utils::neg_one_pow;
use crate::xsf::cephes::psi;
use crate::zeta;

/// `u (1 + r eps)^power` for `power = +-1`, where `u` is a power series in `eps`
fn mul_linear(u: &mut [f64], r: f64, power: i32) {
    if power > 0 {
        for i in (1..u.len()).rev() {
            u[i] += r * u[i - 1];
        }
    } else {
        for i in 1..u.len() {
            u[i] -= r * u[i - 1];
        }
    }
}

/// `mul_linear` for `u` and for the bounds on the sizes of the parts of its coefficients
fn mul_linear_sized(u: &mut [Vec<f64>; 2], r: f64, power: i32) {
    mul_linear(&mut u[0], r, power);
    mul_linear(&mut u[1], r.abs() * f64::from(power), power);
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

/// The distance `d` of `x` to the pole of the Gamma function nearest to it, for `x < 1/2`
fn pole_distance(x: Dd) -> Option<f64> {
    (x.hi < 0.5).then(|| f64::from(x - x.hi.round()))
}

/// The coefficient of `eps^j` of `ln Gamma(x + eps)`, or of `ln((d + eps) Gamma(x + eps))` if `x`
/// has a `pole_distance`, for `j >= 1`
fn ln_gamma_coef(j: u32, x: Dd) -> f64 {
    let (sign, power) = (neg_one_pow(j), f64::from(j));
    // psi^(j-1)(x) / j!, which is (-1)^j zeta(j, x) / j for j > 1 (DLMF 25.11.12)
    // (by `psi`, as `digamma` is much slower near its negative root)
    let taylor = |x: f64| match j {
        1 => psi(x),
        _ => sign * zeta(power, x) / power,
    };
    let Some(d) = pole_distance(x) else {
        return taylor(x.hi);
    };
    // by the reflection formula (DLMF 5.5.3), for x and for d, (d + eps) Gamma(x + eps) is
    // +-Gamma(1 + d + eps) Gamma(1 - d - eps) / Gamma(1 - x - eps)
    taylor(1.0 + d) + sign * (taylor(1.0 - d) - taylor(1.0 - x.hi))
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
    /// `u`, and bounds on the sizes of the parts of its coefficients
    u: [Vec<f64>; 2],
}

impl Residues {
    /// For the class of `a_m` and `z < 0`, where `b` is without the 1 of `n!`, and no `a_i` is
    /// a non-positive integer
    pub(super) fn of_class(a: &[f64], b: &[f64], m: usize, z: f64) -> Self {
        let am = Dd::new(a[m]);
        let diffs = a.iter().map(|&ai| Dd::new(ai) - am);
        let degree = diffs.clone().filter(|d| d.is_integer()).count() - 1;
        // Gamma(-s), the Gamma(a_i + s), and the 1 / Gamma(b_j + s), as the x, power and sign of
        // Gamma(x + sign eps)^power
        let gammas = once((am, 1, -1.0))
            .chain(diffs.map(|d| (d, 1, 1.0)))
            .chain(b.iter().map(|&bj| (Dd::new(bj) - am, -1, 1.0)));
        // The logarithm of the integrand at n = 0 without the factors of the nearest poles, as
        // its coefficients of eps, eps^2, ..., and the sizes of their parts, from (-z)^eps on.
        let first = |l: f64| once(l).chain(repeat(0.0)).take(degree).collect::<Vec<_>>();
        let mut coef = [(-z).ln(), (-z).ln().abs()].map(first);
        for (x, power, sign) in gammas.clone() {
            let mut signs = f64::from(power);
            for j in 0..degree {
                signs *= sign;
                let c = signs * ln_gamma_coef(j as u32 + 1, x);
                (coef[0][j], coef[1][j]) = (coef[0][j] + c, coef[1][j] + c.abs());
            }
        }
        // those factors (d + sign eps)^-power: eps^-power for d = 0, else times d^-power
        let (mut order, mut u) = (0, coef.map(|c| exp_series(&c)));
        for (x, power, sign) in gammas {
            match pole_distance(x) {
                Some(0.0) => order += power,
                Some(d) => mul_linear_sized(&mut u, sign / d, -power),
                None => {},
            }
        }
        // a_m and the 1 + a_m - b_j, over the 1 + a_m - a_i
        let up = once(am).chain(b.iter().map(|&bj| am + 1.0 - bj)).collect();
        let low = a.iter().map(|&ai| am + 1.0 - ai).collect();
        // (-1)^(p+q+1) / z
        let w = neg_one_pow(a.len() + b.len() + 1) / z;
        let (n, term) = (0.0, 1.0);
        Self {
            up,
            low,
            w,
            n,
            term,
            order,
            u,
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
        let weight = k.checked_sub(1).map_or(0.0, |i| self.u[0][i]);
        let size = self.u[1][..k].iter().copied().fold(0.0, f64::max);
        let item = (self.term * weight, self.term.abs() * size);
        let factors = zip(&self.up, repeat(1)).chain(zip(&self.low, repeat(-1)));
        for (&x, power) in factors {
            let c = f64::from(x + self.n);
            if c == 0.0 {
                (self.term, self.order) = (-self.term, self.order - power);
            } else {
                self.term = if power > 0 {
                    self.term * c
                } else {
                    self.term / c
                };
                mul_linear_sized(&mut self.u, -1.0 / c, power);
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
        mul_linear(&mut u, -4.0, 1);
        assert_eq!(u, [1.0, -2.0, -5.0]);
        mul_linear(&mut u, -4.0, -1);
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
    fn test_pole_distance() {
        assert_eq!(pole_distance(Dd::new(0.5)), None);
        assert_eq!(pole_distance(Dd::new(0.25)), Some(0.25));
        assert_eq!(pole_distance(Dd::new(-2.0)), Some(0.0));
        // exact, unlike that of the rounded -0.9999999999
        assert_eq!(pole_distance(Dd::new(1e-10) - 1.0), Some(1e-10));
    }

    #[test]
    fn test_ln_gamma_coef() {
        let check = |j, x: Dd, expected: f64| {
            let c = ln_gamma_coef(j, x);
            let tol = 4.0 * f64::EPSILON * expected.abs().max(1.0);
            assert!((c - expected).abs() <= tol, "{j} {}: {c}", x.hi);
        };
        let zeta_2 = PI * PI / 6.0;
        for (j, x, expected) in [
            // psi(x), psi'(x) / 2, and psi''(1) / 6
            (1, 1.0, -EGAMMA),
            (2, 1.0, zeta_2 / 2.0),
            (2, 3.0, (zeta_2 - 1.25) / 2.0),
            (3, 1.0, -ZETA_3 / 3.0),
            // at the poles 0 and -2: of Gamma(1 + eps), over (eps - 1) (eps - 2) for the latter
            (1, 0.0, -EGAMMA),
            (2, 0.0, zeta_2 / 2.0),
            (1, -2.0, 1.5 - EGAMMA),
            (2, -2.0, (zeta_2 + 1.25) / 2.0),
            (3, -2.0, (1.125 - ZETA_3) / 3.0),
            // at a distance 1/2 from the pole -1: psi(x) + 2
            (1, -0.5, 4.0 - EGAMMA - 2.0 * LN_2),
        ] {
            check(j, Dd::new(x), expected);
        }
        // x - k, where the low part of x is most of its distance to the nearest pole, or
        // a large part of it (Arb references)
        for (j, [x, k], expected) in [
            (1, [1e-10, 1.0], 0.422_784_335_362_960_56),
            (2, [1e-10, 1.0], 1.322_467_033_403_907_5),
            (3, [1e-10, 1.0], -0.067_352_300_844_965_77),
            (1, [-1e-10, 3.0], 1.256_117_668_131_195_9),
            (3, [-1e-10, 3.0], -0.013_339_955_589_902_641),
            // (the high part of x is a pole, unlike x)
            (1, [1e-20, 1.0], 1.0 - EGAMMA),
            (1, [-0.3, 40.0], 2.657_874_444_385_43),
            (2, [0.3, 40.0], 1.971_714_126_134_391_6),
            (3, [0.3, 40.0], 0.872_651_943_111_215_1),
            (1, [0.1, 1e9], 21.054_437_846_882_387),
            (4, [0.1, 1e9], 0.593_828_786_835_143_2),
        ] {
            check(j, Dd::new(x) - k, expected);
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
        assert!(close(residues.u[0][1], 3.0_f64.ln()));
        // 3F2(1, 1, 1; 2, 2; z) = Li_2(z) / z, which is
        // (ln(-z)^2 / 2 + pi^2 / 6 + Li_2(1 / z)) / -z: a pole of order 3
        let (l, zeta_2) = (40.0_f64.ln(), PI * PI / 6.0);
        let residues = Residues::of_class(&[1.0; 3], &[2.0; 2], 0, -40.0);
        assert_eq!(residues.order, 3);
        assert!(close(residues.u[0][2], l * l / 2.0 + zeta_2));
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
