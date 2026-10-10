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
/// The function ratevl is only used once in cephes/lanczos.h.
///
/// Unlike the C++ implementation, this uses fused multiply-adds, which reduces the max error of
/// `lanczos_sum_expg_scaled` on `(0, 172)` from 6.7 to 5.1 ulp.
///
/// \[1\] Holin et. al., "Polynomial and Rational Function Evaluation",
///     <https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots/rational.html>
#[inline]
pub(crate) fn ratevl(x: f64, num: &[f64], denom: &[f64]) -> f64 {
    fn horner(y: f64, coef: impl Iterator<Item = f64>) -> f64 {
        coef.reduce(|acc, c| acc.mul_add(y, c)).unwrap()
    }

    let absx = x.abs();

    if absx > 1.0 {
        /* Evaluate as a polynomial in 1/x. */
        let y = 1.0 / x;
        let num_ans = horner(y, num.iter().rev().copied());
        let denom_ans = horner(y, denom.iter().rev().copied());

        let i = f64::from(i32::try_from(num.len()).unwrap() - i32::try_from(denom.len()).unwrap());
        x.powf(i) * num_ans / denom_ans
    } else {
        horner(x, num.iter().copied()) / horner(x, denom.iter().copied())
    }
}
