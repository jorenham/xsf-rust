//! NOTE: The Cephes library only provides a real implementation of the Spence function.
//! For complex argument, `scipy.special.spence` uses the Cython implementation at
//! `scipy/special/_spence.pxd`.

use crate::ffi;

/// Spence's function, also known as the dilogarithm.
///
/// It is defined to be
///
/// ∫₁ˣ log(t)/(1 - t) dt
///
/// It currently only supports real arguments.
pub fn spence(x: f64) -> f64 {
    unsafe { ffi::xsf::spence(x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_spence() {
        xsref::test::<f64, _>("spence", "d-d", |x: &[f64]| spence(x[0]));
    }
}
