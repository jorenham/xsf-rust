use crate::bindings;
use crate::bindings::xsf_impl;
use num_complex::Complex;

xsf_impl!(ber, (x: f64), "Kelvin function `ber`");
xsf_impl!(bei, (x: f64), "Kelvin function `bei`");
xsf_impl!(ker, (x: f64), "Kelvin function `ker`");
xsf_impl!(kei, (x: f64), "Kelvin function `kei`");

xsf_impl!(berp, (x: f64), "Derivative of the Kelvin function `ber`");
xsf_impl!(beip, (x: f64), "Derivative of the Kelvin function `bei`");
xsf_impl!(kerp, (x: f64), "Derivative of the Kelvin function `ker`");
xsf_impl!(keip, (x: f64), "Derivative of the Kelvin function `kei`");

/// Kelvin functions as complex numbers
///
/// # Arguments
/// - `x` - Real argument
///
/// # Returns
/// - *Be*: Value of the Kelvin function [`ber`] + *i* [`bei`]
/// - *Ke*: Value of the Kelvin function [`ker`] + *i* [`kei`]
/// - *Be*': Derivative of the Kelvin function [`berp`] + *i* [`beip`]
/// - *Ke*': Derivative of the Kelvin function [`kerp`] + *i* [`keip`]
pub fn kelvin(x: f64) -> (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>) {
    let mut be = bindings::complex_nan();
    let mut ke = bindings::complex_nan();
    let mut bep = bindings::complex_nan();
    let mut kep = bindings::complex_nan();
    unsafe {
        bindings::kelvin(x, &mut be, &mut ke, &mut bep, &mut kep);
    }
    (be.into(), ke.into(), bep.into(), kep.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::Complex;

    #[test]
    fn test_ber_f64() {
        xsref::test::<f64, _>("ber", "d-d", |x: &[f64]| ber(x[0]));
    }

    #[test]
    fn test_bei_f64() {
        xsref::test::<f64, _>("bei", "d-d", |x: &[f64]| bei(x[0]));
    }

    #[test]
    fn test_ker_f64() {
        xsref::test::<f64, _>("ker", "d-d", |x: &[f64]| ker(x[0]));
    }

    #[test]
    fn test_kei_f64() {
        xsref::test::<f64, _>("kei", "d-d", |x: &[f64]| kei(x[0]));
    }

    #[test]
    fn test_berp_f64() {
        xsref::test::<f64, _>("berp", "d-d", |x: &[f64]| berp(x[0]));
    }

    #[test]
    fn test_beip_f64() {
        xsref::test::<f64, _>("beip", "d-d", |x: &[f64]| beip(x[0]));
    }

    #[test]
    fn test_kerp_f64() {
        xsref::test::<f64, _>("kerp", "d-d", |x: &[f64]| kerp(x[0]));
    }

    #[test]
    fn test_keip_f64() {
        xsref::test::<f64, _>("keip", "d-d", |x: &[f64]| keip(x[0]));
    }

    #[test]
    fn test_kelvin_c64() {
        xsref::test::<(Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>), _>(
            "kelvin",
            "d-cd_cd_cd_cd",
            |x: &[f64]| kelvin(x[0]),
        );
    }
}
