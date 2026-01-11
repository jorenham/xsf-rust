/// Round to nearest or even integer-valued float
///
/// Returns the nearest integer to x as a f64 precision floating point result.
/// If x ends in 0.5 exactly, the nearest even integer is chosen.
///
/// Corresponds to [`scipy.special.round`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.round.html
#[doc(alias = "round_even")]
#[must_use]
#[inline]
pub fn round(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::round(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_round() {
        xsref::test("round", "d-d", |x| crate::round(x[0]));
    }
}
