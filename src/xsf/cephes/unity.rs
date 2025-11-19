/// see `scipy.special._ufuncs._lgam1p`
#[doc(hidden)]
#[inline(always)]
pub fn lgam1p(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::lgam1p(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_lgam1p() {
        xsref::test("lgam1p", "d-d", |x| crate::lgam1p(x[0]));
    }
}
