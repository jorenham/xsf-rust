use alloc::vec::Vec;
use core::ffi::c_int;
use num_complex::Complex;

/// Kelvin function `ber`
pub fn ber(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::ber(x) }
}

/// Kelvin function `bei`
pub fn bei(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::bei(x) }
}

/// Kelvin function `ker`
pub fn ker(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::ker(x) }
}

/// Kelvin function `kei`
pub fn kei(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kei(x) }
}

/// Derivative of the Kelvin function `ber`
pub fn berp(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::berp(x) }
}

/// Derivative of the Kelvin function `bei`
pub fn beip(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::beip(x) }
}

/// Derivative of the Kelvin function `ker`
pub fn kerp(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kerp(x) }
}

/// Derivative of the Kelvin function `kei`
pub fn keip(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::keip(x) }
}

/// Kelvin functions as complex numbers
///
/// Returned in a length-4 array of complex numbers:
/// - *Be*: Value of [`ber`] + *i* [`bei`]
/// - *Ke*: Value of [`ker`] + *i* [`kei`]
/// - *Be*': Derivative of [`berp`] + *i* [`beip`]
/// - *Ke*': Derivative of [`kerp`] + *i* [`keip`]
pub fn kelvin(x: f64) -> [Complex<f64>; 4] {
    let mut be = f64::NAN.into();
    let mut ke = f64::NAN.into();
    let mut bep = f64::NAN.into();
    let mut kep = f64::NAN.into();
    unsafe {
        crate::ffi::xsf::kelvin(x, &mut be, &mut ke, &mut bep, &mut kep);
    }
    [be.into(), ke.into(), bep.into(), kep.into()]
}

#[repr(C)]
enum KelvinFunction {
    Ber = 1,
    Bei = 2,
    Ker = 3,
    Kei = 4,
    Berp = 5,
    Beip = 6,
    Kerp = 7,
    Keip = 8,
}

#[inline(always)]
fn klvnzo(nt: usize, kd: KelvinFunction) -> Vec<f64> {
    assert!(nt <= c_int::MAX as usize);

    let mut zs = alloc::vec![f64::NAN; nt];
    unsafe {
        crate::ffi::xsf::klvnzo(nt as c_int, kd as c_int, zs.as_mut_ptr());
    }
    zs
}

/// First `nt` zeros of Kelvin function [`ber`]
pub fn ber_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Ber)
}

/// First `nt` zeros of Kelvin function [`bei`]
pub fn bei_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Bei)
}

/// First `nt` zeros of Kelvin function [`ker`]
pub fn ker_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Ker)
}

/// First `nt` zeros of Kelvin function [`kei`]
pub fn kei_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Kei)
}

/// First `nt` zeros of Kelvin function derivative [`berp`]
pub fn berp_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Berp)
}

/// First `nt` zeros of Kelvin function derivative [`beip`]
pub fn beip_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Beip)
}

/// First `nt` zeros of Kelvin function derivative [`kerp`]
pub fn kerp_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Kerp)
}

/// First `nt` zeros of Kelvin function derivative [`keip`]
pub fn keip_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Keip)
}

/// First `nt` zeros of all Kelvin functions and their derivatives
///
/// Returned in a length-8 array of vectors of length nt. The array contains the vectors of zeros of
/// - [`ber`]
/// - [`bei`]
/// - [`ker`]
/// - [`kei`]
/// - [`berp`]
/// - [`beip`]
/// - [`kerp`]
/// - [`keip`]
pub fn kelvin_zeros(nt: usize) -> [Vec<f64>; 8] {
    [
        klvnzo(nt, KelvinFunction::Ber),
        klvnzo(nt, KelvinFunction::Bei),
        klvnzo(nt, KelvinFunction::Ker),
        klvnzo(nt, KelvinFunction::Kei),
        klvnzo(nt, KelvinFunction::Berp),
        klvnzo(nt, KelvinFunction::Beip),
        klvnzo(nt, KelvinFunction::Kerp),
        klvnzo(nt, KelvinFunction::Keip),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::Complex;

    #[test]
    fn test_ber_f64() {
        testing::test::<f64, _>("ber", "d-d", |x: &[f64]| ber(x[0]));
    }

    #[test]
    fn test_bei_f64() {
        testing::test::<f64, _>("bei", "d-d", |x: &[f64]| bei(x[0]));
    }

    #[test]
    fn test_ker_f64() {
        testing::test::<f64, _>("ker", "d-d", |x: &[f64]| ker(x[0]));
    }

    #[test]
    fn test_kei_f64() {
        testing::test::<f64, _>("kei", "d-d", |x: &[f64]| kei(x[0]));
    }

    #[test]
    fn test_berp_f64() {
        testing::test::<f64, _>("berp", "d-d", |x: &[f64]| berp(x[0]));
    }

    #[test]
    fn test_beip_f64() {
        testing::test::<f64, _>("beip", "d-d", |x: &[f64]| beip(x[0]));
    }

    #[test]
    fn test_kerp_f64() {
        testing::test::<f64, _>("kerp", "d-d", |x: &[f64]| kerp(x[0]));
    }

    #[test]
    fn test_keip_f64() {
        testing::test::<f64, _>("keip", "d-d", |x: &[f64]| keip(x[0]));
    }

    #[test]
    fn test_kelvin_c64() {
        testing::test::<(Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>), _>(
            "kelvin",
            "d-cd_cd_cd_cd",
            |x: &[f64]| kelvin(x[0]).into(),
        );
    }

    // Kelvin function zeros tests

    fn assert_allclose(actual: &[f64], expected: &[f64], atol: f64) {
        assert_eq!(actual.len(), expected.len());
        for (&a, &e) in actual.iter().zip(expected.iter()) {
            assert!(
                (a - e).abs() <= atol,
                "actual: {}, expected: {}, diff: {}",
                a,
                e,
                (a - e).abs()
            );
        }
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_ber_zeros`
    #[test]
    fn test_ber_zeros() {
        let ber = ber_zeros(5);
        assert_allclose(
            &ber,
            &[2.84892, 7.23883, 11.67396, 16.11356, 20.55463],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_bei_zeros`
    #[test]
    fn test_bei_zeros() {
        // Abramowitz & Stegun, Table 9.12
        let bi = bei_zeros(5);
        assert_allclose(
            &bi,
            &[5.02622, 9.45541, 13.89349, 18.33398, 22.77544],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_ker_zeros`
    #[test]
    fn test_ker_zeros() {
        let ker = ker_zeros(5);
        assert_allclose(
            &ker,
            &[1.71854, 6.12728, 10.56294, 15.00269, 19.44381],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_kei_zeros`
    #[test]
    fn test_kei_zeros() {
        let kei = kei_zeros(5);
        assert_allclose(
            &kei,
            &[3.91467, 8.34422, 12.78256, 17.22314, 21.66464],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_berp_zeros`
    #[test]
    fn test_berp_zeros() {
        let brp = berp_zeros(5);
        assert_allclose(
            &brp,
            &[6.03871, 10.51364, 14.96844, 19.41758, 23.86430],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_beip_zeros`
    #[test]
    fn test_beip_zeros() {
        let bip = beip_zeros(5);
        assert_allclose(
            &bip,
            &[
                3.772_673_304_934_953,
                8.280_987_849_760_042,
                12.742_147_523_633_703,
                17.193_431_752_512_54,
                21.641_143_941_167_325,
            ],
            1.5e-8,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_kerp_zeros`
    #[test]
    fn test_kerp_zeros() {
        let kerp = kerp_zeros(5);
        assert_allclose(
            &kerp,
            &[2.66584, 7.17212, 11.63218, 16.08312, 20.53068],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_keip_zeros`
    #[test]
    fn test_keip_zeros() {
        let keip = keip_zeros(5);
        assert_allclose(
            &keip,
            &[4.93181, 9.40405, 13.85827, 18.30717, 22.75379],
            1.5e-4,
        );
    }

    /// Ported from `scipy.special.tests.test_basic.TestKelvin.test_kelvin_zeros`
    #[test]
    fn test_kelvin_zeros() {
        // numbers come from 9.9 of A&S pg. 381
        let tmp = kelvin_zeros(5);
        let [berz, beiz, kerz, keiz, berpz, beipz, kerpz, keipz] = tmp;

        assert_allclose(
            &berz,
            &[2.84892, 7.23883, 11.67396, 16.11356, 20.55463],
            1.5e-4,
        );
        assert_allclose(
            &beiz,
            &[5.02622, 9.45541, 13.89349, 18.33398, 22.77544],
            1.5e-4,
        );
        assert_allclose(
            &kerz,
            &[1.71854, 6.12728, 10.56294, 15.00269, 19.44382],
            1.5e-4,
        );
        assert_allclose(
            &keiz,
            &[3.91467, 8.34422, 12.78256, 17.22314, 21.66464],
            1.5e-4,
        );
        assert_allclose(
            &berpz,
            &[6.03871, 10.51364, 14.96844, 19.41758, 23.86430],
            1.5e-4,
        );
        assert_allclose(
            &beipz,
            // table from 1927 had 3.77320 but this is more accurate
            &[3.77267, 8.28099, 12.74215, 17.19343, 21.64114],
            1.5e-4,
        );
        assert_allclose(
            &kerpz,
            &[2.66584, 7.17212, 11.63218, 16.08312, 20.53068],
            1.5e-4,
        );
        assert_allclose(
            &keipz,
            &[4.93181, 9.40405, 13.85827, 18.30717, 22.75379],
            1.5e-4,
        );
    }
}
