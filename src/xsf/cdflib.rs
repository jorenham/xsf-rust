/// Inverse of [`gdtr`](xsf::gdtr) with respect to `b`
///
/// Corresponds to [`scipy.special.gdtrib`][gdtrib].
///
/// [gdtrib]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtrib.html
///
/// # See also
/// - [`gdtr`]: Gamma CDF
/// - [`gdtrc`]: Gamma survival function
#[must_use]
#[inline]
pub fn gdtrib(a: f64, p: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gdtrib(a, p, x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_gdtrib() {
        xsref::test("gdtrib", "d_d_d-d", |x| crate::gdtrib(x[0], x[1], x[2]));
    }
}
