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
