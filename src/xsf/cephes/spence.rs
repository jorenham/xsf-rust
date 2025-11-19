//! NOTE: The Cephes library only provides a real implementation of the Spence function.
//! For complex argument, `scipy.special.spence` uses the Cython implementation at
//! `scipy/special/_spence.pxd`.

/// Spence's function, also known as the dilogarithm.
///
/// It is defined to be
///
/// ∫₁ˣ log(t)/(1 - t) dt
///
/// It currently only supports real arguments.
#[doc(alias = "dilogarithm")]
pub fn spence(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::spence(x) }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_spence() {
        xsref::test("spence", "d-d", |x| crate::spence(x[0]));
    }
}
