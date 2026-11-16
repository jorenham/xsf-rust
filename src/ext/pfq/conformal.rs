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
use core::iter::zip;

use super::Dd;

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
        assert_eq!(
            coef(&F3::INNER, 2.5, 8.0, 2),
            [14.5, 40.5, 55.5, 13.5, 20.0]
        );
        assert_eq!(coef(&F3::LAST, 2.5, 0.0, 2), [14.5, 1.0, 2.5, -2.0, 0.0]);
        assert_eq!(coef(&F3::TAIL, 0.0, 8.0, 2), [14.0, 34.0, 48.0, 35.0, 10.0]);
        assert_eq!(
            coef(&F3::INNER, -1.5, 5.0, 0),
            [-1.5, 34.5, 33.5, 11.5, 18.0]
        );
        let expected = [18.5, 78.0, 103.0, 78.0, 106.5, 0.0];
        assert_eq!(coef(&F4::INNER, 2.5, 11.0, 2), expected);
        assert_eq!(
            coef(&F4::LAST, 2.5, 0.0, 2),
            [18.5, 0.5, 7.5, -10.5, 0.0, 0.0]
        );
        let expected = [19.0, 60.0, 110.0, 115.0, 63.0, 14.0];
        assert_eq!(coef(&F4::TAIL, 0.0, 11.0, 2), expected);
        let expected = [-1.5, 62.0, 47.0, 62.0, 86.5, 0.0];
        assert_eq!(coef(&F4::INNER, -1.5, 7.0, 0), expected);
    }

    #[test]
    fn test_factor_dot() {
        // the coefficients 8, 30, 52, 41, 12
        let factor = Factor::new(&F3::TAIL, Dd::new(0.0), 8.0);
        let x = [1.0, 2.0, 3.0, 4.0, 5.0].map(Dd::new);
        assert_eq!(f64::from(factor.dot(&x, 0)), 410.0);
        assert_eq!(f64::from(factor.dot(&x, 1)), 505.0);
    }
}
