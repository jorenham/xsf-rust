/// Cube root of $x$, $\sqrt\[3\]{x}$
///
/// This corresponds to [`scipy.special.cbrt`][cbrt] in SciPy.
///
/// [cbrt]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cbrt.html
#[must_use]
#[inline]
pub fn cbrt(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cbrt(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_cbrt() {
        xsref::test("cbrt", "d-d", |x| crate::cbrt(x[0]));
    }
}
