use crate::ffi;

/// Rising factorial
///
/// It is defined as `gamma(x + m) / gamma(x)`.
///
/// Note that in the Cephes library and `scipy.special` this function is called `poch`.
///
/// See [`pow_falling`] for the falling factorial.
pub fn pow_rising(x: f64, m: f64) -> f64 {
    unsafe { ffi::poch(x, m) }
}

/// Falling factorial
///
/// It is defined as `gamma(x + 1) / gamma(x - m + 1)`.
///
/// Note that there is no `scipy.special` analogue for this function, but it can be expressed in
/// terms of the rising factorial as `pow_rising(x - m + 1, m)`. See [`pow_rising`] for the details.
pub fn pow_falling(x: f64, m: f64) -> f64 {
    unsafe { ffi::poch(x - m + 1.0, m) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_pow_rising() {
        testing::test::<f64, _>("poch", "d_d-d", |x: &[f64]| pow_rising(x[0], x[1]));
    }
}
