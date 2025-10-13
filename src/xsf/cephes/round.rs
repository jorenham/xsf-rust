/// Round to nearest or even integer-valued float
///
/// Returns the nearest integer to x as a f64 precision floating point result.
/// If x ends in 0.5 exactly, the nearest even integer is chosen.
#[doc(alias = "round_even")]
pub fn round(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::round(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_round() {
        crate::xsref::test("round", "d-d", |x| crate::round(x[0]));
    }
}
