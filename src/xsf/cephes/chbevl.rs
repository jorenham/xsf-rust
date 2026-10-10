//! Translated into pure Rust from `xsf/cephes/chbevl.h` (xsf v0.2.2).

/// Evaluate a Chebyshev series, with the coefficients in reverse order (`n = array.len()`)
pub(crate) fn chbevl(x: f64, array: &[f64]) -> f64 {
    let mut b0 = array[0];
    let mut b1 = 0.0;
    let mut b2 = 0.0;

    for &c in &array[1..] {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + c;
    }

    0.5 * (b0 - b2)
}
