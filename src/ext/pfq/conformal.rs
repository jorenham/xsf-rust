//! `_{q+1}F_q` by its series in the conformal variable `w`, where
//! `z = 1 - ((1 - w) / (1 + w))^m`
//!
//! Ported from `pFqconformalrational` of
//! [HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl) v0.3.31:
//!
//! > MIT License
//! >
//! > Copyright (c) 2018-2023 Richard Mikael Slevinsky and other contributors:
//! >
//! > <https://github.com/JuliaMath/HypergeometricFunctions.jl/graphs/contributors>
//! >
//! > Permission is hereby granted, free of charge, to any person obtaining a copy
//! > of this software and associated documentation files (the "Software"), to deal
//! > in the Software without restriction, including without limitation the rights
//! > to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! > copies of the Software, and to permit persons to whom the Software is
//! > furnished to do so, subject to the following conditions:
//! >
//! > The above copyright notice and this permission notice shall be included in all
//! > copies or substantial portions of the Software.
//! >
//! > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//! > IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! > FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//! > AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! > LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//! > OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//! > SOFTWARE.
//!
//! The coefficients of the series follow from a recurrence, which equates two products of
//! `q + 1` factors of order `m + 1`. The functions are generic over `M = m + 2`.

use core::array::from_fn;
use core::f64::consts::PI;
use core::iter::zip;

use super::{Dd, PFQ_DD_MAXITER, PFQ_TOL, Params, SeriesReal, pfq_settled, root_dd, skip};
use crate::utils::pow2;

/// Largest size of a negligible term, relative to the sum of the series
///
/// On 400 points for which the other methods of `pfq_accurate` are not accurate, as many sums
/// are accurate to `PFQ_TOL` as with `2^-105`, in 0.8 times the time.
const CONFORMAL_EPS: f64 = pow2(-70);

/// The number of successive negligible terms after which the series ends
///
/// On 5525 points, this alone prevents 1 of the 8 wrong results of ending at the first
/// negligible term, and `min_terms` alone all of them.
const CONFORMAL_NEGLIGIBLE: u32 = 3;

/// Smallest `-ln |w^n|` at the term `n = PFQ_DD_MAXITER` for which the series are summed
///
/// It is half of `-ln CONFORMAL_EPS`: on 2320 points, each series took at least 0.55 times
/// `ln CONFORMAL_EPS / ln |w|` terms, so that below it the series would not end. Only a
/// series with tiny parameters can end earlier.
const CONFORMAL_MIN_DECAY: f64 = 24.0;

/// The factor from the difference of the sums for `m = 3` and `m = 4` to the estimated error
///
/// On 5525 points, the error was at most 22 times that difference where it is below `PFQ_TOL`.
/// As the difference is no bound, this keeps a wide margin: 100 would return a value for only
/// 26 more of the 681 inputs in `docs/accuracy.md` that give `None`.
const CONFORMAL_SAFETY: f64 = 1e3;

/// The binomial coefficient `C(n, k)`, which is zero for `k > n`
const fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    // (each division is exact)
    let (mut c, mut i) = (1_u64, 0);
    while i < k {
        c = c * (n - i) as u64 / (i + 1) as u64;
        i += 1;
    }
    c as f64
}

/// The integers of a kind of factor: for each `k`, the factors of its parameter and of its
/// shift in the coefficient `c_k` for `n = 0`, and the step of `c_k` in `n`
type Rows<const M: usize> = [[f64; M]; 3];

/// The index in `Rows` of the factors of the parameter
const BY_PARAM: usize = 0;
/// The index in `Rows` of the factors of the shift
const BY_SHIFT: usize = 1;
/// The index in `Rows` of the steps
const STEP: usize = 2;

/// `Factor::INNER`, `Factor::LAST` and `Factor::TAIL`
const fn rows<const M: usize>() -> [Rows<M>; 3] {
    let m = M - 2;
    let (mut inner, mut last, mut tail) = ([[0.0; M]; 3], [[0.0; M]; 3], [[0.0; M]; 3]);
    let mut k = 0;
    while k < M {
        let even = k % 2 == 0;
        let sign = if even { 1.0 } else { -1.0 };
        // (`k - 1` wraps for `k = 0`, to a `k > m`)
        let before = k.wrapping_sub(1);
        // the odd numbers next to `k`
        let (below, above) = if even { (before, k + 1) } else { (k, k) };
        let diff = binom(m, k + 1) - binom(m, before);
        inner[BY_PARAM][k] = sign * binom(m, k);
        inner[BY_SHIFT][k] = 2.0 * binom(m, below);
        inner[STEP][k] = if even { 2.0 * diff } else { 0.0 };
        last[BY_PARAM][k] = sign * binom(m - 1, k);
        last[STEP][k] = 2.0 * binom(m, above);
        tail[BY_SHIFT][k] = binom(m + 1, k);
        tail[STEP][k] = diff;
        k += 1;
    }
    [inner, last, tail]
}

/// A factor of the recurrence: its coefficients `c_k` for term `n`, and their steps in `n`
///
/// It takes a sequence `x` to the sequence `sum_{k=0}^{m+1} c_k x_{n-k}`.
pub(super) struct Factor<const M: usize> {
    coef: [Dd; M],
    step: [f64; M],
}

impl<const M: usize> Factor<M> {
    /// For a parameter but the last `a_i` and `b_j`: `L[:, j]` and `R[:, j]` in the Julia code,
    /// and their steps
    pub(super) const INNER: Rows<M> = rows::<M>()[0];
    /// For the last `a_i` and `b_j`: `L[:, q]` and `R[:, q + 1]`, and their steps
    pub(super) const LAST: Rows<M> = rows::<M>()[1];
    /// For no parameter: `Lt0`, and its steps
    pub(super) const TAIL: Rows<M> = rows::<M>()[2];

    /// The factor of a kind for `n = 0`, where `c_k` is `-k` times its step for a zero
    /// parameter and shift
    pub(super) fn new(rows: &Rows<M>, param: Dd, shift: f64) -> Self {
        let [by_param, by_shift, step] = *rows;
        let coef = from_fn(|k| param * by_param[k] + (shift * by_shift[k] - k as f64 * step[k]));
        Self { coef, step }
    }

    /// To the next `n`
    pub(super) fn advance(&mut self) {
        for (c, step) in zip(&mut self.coef, self.step) {
            *c = *c + step;
        }
    }

    /// `sum_{k >= skip} c_k x_{n+skip-k}` for `x = [x_{n-m-1}, ..., x_n]`
    pub(super) fn dot(&self, x: &[Dd; M], skip: usize) -> Dd {
        Dd::dot(zip(&self.coef[skip..], x.iter().rev()).map(|(&c, &x)| (c, x)))
    }
}

/// `w = (1 - t) / (1 + t)` for `t = (1 - z)^(1 / m)` and `m >= 2`, without the cancellation
/// for a small `z`
fn conformal_map(z: f64, m: usize) -> Dd {
    let t = root_dd(1.0 - Dd::new(z), m);
    // z = 1 - t^m = (1 - t) (1 + t + ... + t^(m-1))
    let u = t + 1.0;
    let s = (2..m).fold(u, |s, _| t * s + 1.0);
    z / (u * s)
}

/// Appends `new` to the values `x_{n-k} w^k`
fn push<const M: usize>(x: &mut [Dd; M], w: Dd, new: Dd) {
    for k in 1..M {
        x[k - 1] = w * x[k];
    }
    x[M - 1] = new;
}

/// The number of terms after which the series in `w` for `m` may end, or infinity if the
/// product of the `2 m a_i` in its first term underflows, so that its terms are lost
fn min_terms(alpha: &[f64], beta: &[f64], w: f64, m: usize) -> f64 {
    let first: f64 = alpha
        .iter()
        .map(|a| ((2 * m) as f64 * a).abs().min(1.0))
        .product();
    if first < Dd::MIN_POSITIVE {
        return f64::INFINITY;
    }
    // The coefficients of `w^n` grow at most about like `n^(g - 1) / r^n`, where `-g` is the
    // exponent of a singularity at `|w| = r`: at `w = 1` (`z = 1`), at `w = -1` (`z = inf`),
    // and at `w = +-i tan(pi / m)` (`z = 0` on the other sheets). The terms can grow up to
    // the `n` with `(n + g) |w| = (n + 1) r`, also if the ones before it are negligible.
    let last_growing = |g: f64, r: f64| (g * w.abs() - r) / (r - w.abs());
    let excess = alpha.iter().sum::<f64>() - beta.iter().sum::<f64>();
    let min_alpha = alpha.iter().copied().fold(f64::INFINITY, f64::min);
    let max_beta = beta.iter().copied().fold(0.0, f64::max);
    // past the last term with a divisor `b_j + n - 1` that may be near zero, by the `m + 2`
    // terms in `w` that its jump takes to show
    let last_pole = pfq_settled(alpha, beta).floor() + 2.0 + (m + 2) as f64;
    last_pole
        .max(last_growing(m as f64 * excess.max(-min_alpha), 1.0))
        .max(last_growing(max_beta - 1.0, (PI / m as f64).tan()))
}

/// The sum of the series in `w` for `m = M - 2`, and its largest partial sum
///
/// `alpha` are the `q + 1` upper parameters and `beta` the `q >= 1` lower ones, without the 1
/// of `n!`. Unlike in HypergeometricFunctions.jl, the sum is NaN if it is not finite or if the
/// series did not converge within `PFQ_DD_MAXITER` terms, and it ends at a term that is
/// negligible in `f64` rather than in the arithmetic of the sum. It also ends only after
/// `CONFORMAL_NEGLIGIBLE` successive negligible terms that do not grow, and after
/// `min_terms` terms, as larger terms can follow a negligible one.
pub(super) fn conformal_sum<const M: usize>(alpha: &[f64], beta: &[f64], z: f64) -> (Dd, f64) {
    let (m, q) = (M - 2, beta.len());
    let scale = (2 * m) as f64;
    // the factors are numbered from the last one, for the last parameter
    let shift = |j: usize| (j * m) as f64 - 1.0;
    let factor = |(j, param): (usize, Dd)| {
        let rows = if j == 0 {
            &Factor::<M>::LAST
        } else {
            &Factor::INNER
        };
        Factor::new(rows, param, shift(j))
    };
    // 2 m (b_j - 1), in the lower factors `L`, and 2 m a_i, in the upper factors `R`
    let lower: Vec<Dd> = beta
        .iter()
        .rev()
        .map(|&b| (Dd::new(b) - 1.0) * scale)
        .collect();
    let upper = alpha.iter().rev().map(|&a| Dd::new(a) * scale);
    let mut l: Vec<_> = lower.iter().copied().enumerate().map(factor).collect();
    let mut r: Vec<_> = upper.enumerate().map(factor).collect();
    let mut tail = Factor::new(&Factor::<M>::TAIL, Dd::ZERO, shift(q));
    let w = conformal_map(z, m);
    let n_min = min_terms(alpha, beta, w.hi, m);
    if n_min >= f64::from(PFQ_DD_MAXITER) {
        return (Dd::NAN, 1.0);
    }
    // The values `x_{n-k} w^k` of the sequences between the factors: `L[j]` takes `zs[j]` to
    // `zs[j + 1]`, and `R[j]` takes `zt[j - 1]`, or `zs[0]`, to `zt[j]`, or to `zs[q]` with
    // `tail`. `zs[0]` are the terms of the series.
    let zero = [Dd::ZERO; M];
    let (mut zs, mut zt) = (vec![zero; q + 1], vec![zero; q]);
    zs[0][M - 1] = Dd::ONE;
    for j in 0..q {
        zs[j + 1][M - 1] = l[j].dot(&zs[j], 0) / w;
    }
    let (mut sum, mut top, mut n) = (Dd::ONE, 1.0_f64, 0);
    // the number of the last terms that are negligible, and the sizes of the `M` terms before
    // (which can alternate in size, with a period up to the order of the recurrence)
    let (mut negligible, mut recent) = (0, [1.0_f64; M]);
    let ended = |negligible, n| negligible >= CONFORMAL_NEGLIGIBLE && f64::from(n) >= n_min;
    while n < PFQ_DD_MAXITER && sum.hi.is_finite() && !ended(negligible, n) {
        // `zt` for term n
        let new = r[0].dot(&zs[0], 0);
        push(&mut zt[0], w, new);
        for j in 1..q {
            let new = r[j].dot(&zt[j - 1], 0) / w;
            push(&mut zt[j], w, new);
        }
        // `zs` for term n + 1
        n += 1;
        let new = (r[q].dot(&zt[q - 1], 0) - w * tail.dot(&zs[q], 0)) / f64::from(n);
        push(&mut zs[q], w, new);
        l.iter_mut().for_each(Factor::advance);
        for j in (0..q).rev() {
            let divisor = lower[j] + scale * f64::from(n);
            let new = w * (zs[j + 1][M - 1] - l[j].dot(&zs[j], 1)) / divisor;
            push(&mut zs[j], w, new);
        }
        r.iter_mut().for_each(Factor::advance);
        tail.advance();
        let term = zs[0][M - 1];
        sum = sum + term;
        top = top.max(sum.hi.abs());
        let size = term.hi.abs();
        let small = size <= CONFORMAL_EPS * sum.hi.abs().max((sum - term).hi.abs())
            && size <= recent.into_iter().fold(0.0, f64::max);
        recent.rotate_left(1);
        recent[M - 1] = size;
        negligible = if small { negligible + 1 } else { 0 };
    }
    let converged = ended(negligible, n) && sum.hi.is_finite();
    (if converged { sum } else { Dd::NAN }, top)
}

/// `_{q+1}F_q` for `z < 1` by its series in `w` for `m = 4`, and an estimate of its absolute
/// error, or `None` if the series converge too slowly or the rounding errors are too large
///
/// `b` contains the 1 of `n!`. Unlike in HypergeometricFunctions.jl, which has no estimate of
/// the error, it is `CONFORMAL_SAFETY` times the difference to the sum for `m = 3`.
pub(super) fn pfq_conformal<P: Params>(a: &P, b: &P, z: f64) -> Option<(f64, f64)> {
    let one = b.as_ref().iter().position(|&v| v == 1.0)?;
    let beta: Vec<f64> = skip(b, one).collect();
    // the terms decrease at most about like `w^n`, the slowest for `m = 3`
    let decay = -conformal_map(z, 3).hi.abs().ln() * f64::from(PFQ_DD_MAXITER);
    if beta.is_empty() || !(decay >= CONFORMAL_MIN_DECAY) {
        return None;
    }
    let (y4, top) = conformal_sum::<6>(a.as_ref(), &beta, z);
    // The rounding errors of the partial sums, and those of the first terms, which the
    // divisions by the `b_j + n` near zero amplify alike for both `m`. On those 5525 points,
    // the difference accepts none of the sums that this rejects, which saves the other sum.
    let distance = |b: &f64| (b - b.round().min(0.0)).abs().min(1.0);
    let near_pole: f64 = beta.iter().map(distance).product();
    let noise = Dd::EPSILON * top.max(1.0 / near_pole);
    (noise <= PFQ_TOL * y4.hi.abs()).then(|| {
        let (y3, _) = conformal_sum::<5>(a.as_ref(), &beta, z);
        (y4.hi, CONFORMAL_SAFETY * (y4 - y3).hi.abs())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    type F3 = Factor<5>;
    type F4 = Factor<6>;

    /// The coefficients of a factor for term `n`
    fn coef<const M: usize>(rows: &Rows<M>, param: f64, shift: f64, n: usize) -> [f64; M] {
        let mut factor = Factor::new(rows, Dd::new(param), shift);
        (0..n).for_each(|_| factor.advance());
        factor.coef.map(f64::from)
    }

    #[test]
    fn test_factor() {
        // by the formulas of the Julia code for q = 3, with 2 m a_j = 2.5 and
        // 2 m (b_j - 1) = -1.5: `R[:, 1]`, `R[:, 4]` and `Lt0` for n = 2, and `L[:, 1]`
        assert_eq!(coef(&F3::INNER, 2.5, 8.0, 2), [14.5, 40.5, 55.5, 13.5, 20.0]);
        assert_eq!(coef(&F3::LAST, 2.5, 0.0, 2), [14.5, 1.0, 2.5, -2.0, 0.0]);
        assert_eq!(coef(&F3::TAIL, 0.0, 8.0, 2), [14.0, 34.0, 48.0, 35.0, 10.0]);
        assert_eq!(coef(&F3::INNER, -1.5, 5.0, 0), [-1.5, 34.5, 33.5, 11.5, 18.0]);
        let expected = [18.5, 78.0, 103.0, 78.0, 106.5, 0.0];
        assert_eq!(coef(&F4::INNER, 2.5, 11.0, 2), expected);
        assert_eq!(coef(&F4::LAST, 2.5, 0.0, 2), [18.5, 0.5, 7.5, -10.5, 0.0, 0.0]);
        let expected = [19.0, 60.0, 110.0, 115.0, 63.0, 14.0];
        assert_eq!(coef(&F4::TAIL, 0.0, 11.0, 2), expected);
        let expected = [-1.5, 62.0, 47.0, 62.0, 86.5, 0.0];
        assert_eq!(coef(&F4::INNER, -1.5, 7.0, 0), expected);
    }

    /// `conformal_sum` for m = 3 and m = 4
    fn sums(alpha: &[f64], beta: &[f64], z: f64) -> [Dd; 2] {
        let (y3, y4) = (conformal_sum::<5>(alpha, beta, z), conformal_sum::<6>(alpha, beta, z));
        [y3.0, y4.0]
    }

    #[test]
    fn test_conformal_map() {
        for m in [2, 3, 4] {
            for z in [-1e6, -3.0, 0.5, 0.999_999] {
                // z = 1 - ((1 - w) / (1 + w))^m
                let w = conformal_map(z, m);
                let t = (1.0 - w) / (w + 1.0);
                let pow = (0..m).fold(Dd::ONE, |p, _| p * t);
                let err = (pow - (1.0 - Dd::new(z))).hi.abs();
                assert!(err <= pow2(-90) * (1.0 - z), "m = {m}, z = {z}");
            }
            // w = z / (2 m) + O(z^2)
            let w = conformal_map(1e-20, m);
            assert!((w * (2 * m) as f64 - 1e-20).hi.abs() <= 1e-40);
        }
    }

    #[test]
    fn test_conformal_sum() {
        // the closed forms -ln(1 - z) / z, Li_2(z) / z, Li_3(z) / z and sum_n (n + 1)^4 z^n,
        // and a Wolfram reference near z = 1
        let cases: [(&[f64], &[f64], f64, f64); 5] = [
            (&[1.0; 2], &[2.0], -3.0, 4.0_f64.ln() / 3.0),
            (&[1.0; 3], &[2.0; 2], -50.0, 0.185_539_903_706_652_44),
            (&[1.0; 4], &[2.0; 3], -40.0, 0.361_479_025_355_404_44),
            (&[2.0; 4], &[1.0; 3], -3.0, 0.039_062_5),
            (&[0.5; 3], &[1.5; 2], 0.999_999, 1.088_792_804_483_117),
        ];
        for (alpha, beta, z, expected) in cases {
            for y in sums(alpha, beta, z) {
                assert!((y.hi - expected).abs() <= 2.0 * f64::EPSILON * expected, "z = {z}");
            }
        }
        // beyond `f64` precision, with its largest partial sum
        let (y, top) = conformal_sum::<6>(&[2.0; 4], &[1.0; 3], -3.0);
        assert!((y - 0.039_062_5).hi.abs() <= CONFORMAL_EPS && top > 1e3);
        // no convergence within `PFQ_DD_MAXITER` terms
        assert!(
            sums(&[1.0; 3], &[2.0; 2], -1e12)
                .iter()
                .all(|y| y.hi.is_nan())
        );
    }

    #[test]
    fn test_conformal_sum_end() {
        // negligible first terms, and large later ones; Arb reference
        let beta = [1.0; 2];
        for y in sums(&[60.0, 1e-12, 1e-12], &beta, 0.9) {
            assert!((y.hi - 3.628_990_674_015_697_6e30).abs() <= PFQ_TOL * y.hi);
        }
        // negligible terms that grow, past the near-pole `b_2 + 79`; Arb reference
        let (a4, b3) = ([-1.5, 0.75, -0.125, -77.75], [-3.75, -78.875, 21.5]);
        for y in sums(&a4, &b3, -2.0) {
            assert!((y.hi - 1.441_597_903_532_934).abs() <= 1e-5);
        }
        // negligible terms up to the near-pole `b_1 + 49`, and large ones after it, which only
        // the rounding errors before it make less accurate; Arb reference
        let (a5, b4) = (
            [0.17549, -33.943, -0.1108, 0.011587, 0.04128],
            [-48.912, 0.0012253, 29.137, 0.0093308],
        );
        for y in sums(&a5, &b4, -6.5618) {
            assert!((y.hi - 1.118_370_664_887_746).abs() <= 4e-8);
        }
        // a negligible term before the near-pole `b_1 + 3`; Arb reference
        for y in sums(&[1.5, 1e-23], &[-3.0 + pow2(-51)], 0.5) {
            assert!((y.hi - 0.999_999_925_733_457).abs() <= PFQ_TOL, "{y:?}");
        }
        // the later terms overflow, and the first term underflows
        for a in [1e-30, 1e-200] {
            let ys = sums(&[500.0, a, a], &beta, 0.9);
            assert!(ys.iter().all(|y| y.hi.is_nan()), "a = {a}");
        }
    }

    #[test]
    fn test_pfq_conformal() {
        let (y, err) = pfq_conformal(&[1.0; 3], &[2.0, 2.0, 1.0], -50.0).unwrap();
        assert!((y - 0.185_539_903_706_652_44).abs() <= f64::EPSILON && err <= PFQ_TOL * y);
        // the series converge too slowly
        assert_eq!(pfq_conformal(&[1.0; 3], &[2.0, 2.0, 1.0], -1e7), None);
        // the partial sums are too large
        assert_eq!(pfq_conformal(&[70.5, 60.5, 50.5], &[2.5, 1.5, 1.0], -0.75), None);
        // 1F0
        assert_eq!(pfq_conformal(&[0.5], &[1.0], -0.75), None);
    }

    #[test]
    fn test_factor_dot() {
        // the coefficients 8, 30, 52, 41, 12
        let factor = Factor::new(&F3::TAIL, Dd::ZERO, 8.0);
        let x = [1.0, 2.0, 3.0, 4.0, 5.0].map(Dd::new);
        assert_eq!(f64::from(factor.dot(&x, 0)), 410.0);
        assert_eq!(f64::from(factor.dot(&x, 1)), 505.0);
    }
}
