//! Translated from `cephes/polevl.h`

/// Evaluate the polynomial with coefficients `coef`, stored from the highest degree down
#[inline]
pub(crate) fn polevl(x: f64, coef: &[f64]) -> f64 {
    coef.iter().copied().reduce(|acc, c| acc * x + c).unwrap()
}

/// Like [`polevl`], with an implicit leading coefficient of 1
#[inline]
pub(crate) fn p1evl(x: f64, coef: &[f64]) -> f64 {
    coef.iter().fold(1.0, |acc, &c| acc * x + c)
}

/// Evaluate a rational function. See \[1\].
///
/// \[1\] Holin et. al., "Polynomial and Rational Function Evaluation",
///     <https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots/rational.html>
#[inline]
pub(crate) fn ratevl(x: f64, num: &[f64], denom: &[f64]) -> f64 {
    ratevl_with(x, num, denom, |acc, y, c| acc * y + c)
}

/// Like [`ratevl`], but with fused multiply-adds.
///
/// Unlike the C++ implementation, this reduces the max error of `lanczos_sum_expg_scaled` on
/// `(0, 172)` from 8.5 to 6.8 ulp (rms 0.85 to 0.71 ulp).
#[inline]
pub(crate) fn ratevl_fma(x: f64, num: &[f64], denom: &[f64]) -> f64 {
    ratevl_with(x, num, denom, f64::mul_add)
}

/// [`ratevl`] with `mul_add(acc, y, c)` as the Horner step
#[inline]
fn ratevl_with(x: f64, num: &[f64], denom: &[f64], mul_add: impl Fn(f64, f64, f64) -> f64) -> f64 {
    fn horner(
        y: f64,
        coef: impl Iterator<Item = f64>,
        mul_add: impl Fn(f64, f64, f64) -> f64,
    ) -> f64 {
        coef.reduce(|acc, c| mul_add(acc, y, c)).unwrap()
    }

    let absx = x.abs();

    if absx > 1.0 {
        /* Evaluate as a polynomial in 1/x. */
        let y = 1.0 / x;
        let num_ans = horner(y, num.iter().rev().copied(), &mul_add);
        let denom_ans = horner(y, denom.iter().rev().copied(), &mul_add);

        let i = f64::from(i32::try_from(num.len()).unwrap() - i32::try_from(denom.len()).unwrap());
        x.powf(i) * num_ans / denom_ans
    } else {
        horner(x, num.iter().copied(), &mul_add) / horner(x, denom.iter().copied(), &mul_add)
    }
}
