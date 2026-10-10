/// Round to nearest or even integer-valued float
///
/// Returns the nearest integer to x as a f64 precision floating point result.
/// If x ends in 0.5 exactly, the nearest even integer is chosen.
///
/// This uses [`f64::round_ties_even`], which, unlike the Cephes implementation in xsf, returns
/// `-0.0` for `-0.5 <= x < 0`.
///
/// Corresponds to [`scipy.special.round`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.round.html
#[doc(alias = "round_even")]
#[must_use]
#[inline]
pub fn round(x: f64) -> f64 {
    x.round_ties_even()
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_round() {
        xsref::test("round", "d-d", |x| crate::round(x[0]));
    }
}
