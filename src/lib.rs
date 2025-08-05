mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

pub mod xsf64 {
    use crate::ffi;

    pub fn gamma(x: f64) -> f64 {
        unsafe { ffi::xsf_gamma(x) }
    }

    pub fn gammaln(x: f64) -> f64 {
        unsafe { ffi::xsf_gammaln(x) }
    }

    pub fn gammasgn(x: f64) -> f64 {
        unsafe { ffi::xsf_gammasgn(x) }
    }

    pub fn gammainc(a: f64, x: f64) -> f64 {
        unsafe { ffi::xsf_gammainc(a, x) }
    }

    pub fn gammaincinv(a: f64, p: f64) -> f64 {
        unsafe { ffi::xsf_gammaincinv(a, p) }
    }

    pub fn gammaincc(a: f64, x: f64) -> f64 {
        unsafe { ffi::xsf_gammaincc(a, x) }
    }

    pub fn gammainccinv(a: f64, p: f64) -> f64 {
        unsafe { ffi::xsf_gammainccinv(a, p) }
    }

    pub fn gamma_ratio(a: f64, b: f64) -> f64 {
        unsafe { ffi::xsf_gamma_ratio(a, b) }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::LN_2;

    use super::xsf64::*;
    use float_cmp::assert_approx_eq;

    const LN_3: f64 = 1.098_612_288_668_109_7;
    const SQRT_PI: f64 = 1.772_453_850_905_516;
    const LN_PI_HALF: f64 = 0.572_364_942_924_700;

    #[test]
    fn test_gamma() {
        assert_eq!(gamma(1.0), 1.0);
        assert_eq!(gamma(2.0), 1.0);
        assert_eq!(gamma(3.0), 2.0);
        assert_eq!(gamma(5.0), 24.0);
        assert_eq!(gamma(8.0), 5_040.0);

        assert_approx_eq!(f64, gamma(0.5), SQRT_PI, ulps = 1);
        assert_approx_eq!(f64, gamma(-0.5), -2.0 * SQRT_PI, ulps = 1);

        assert_eq!(gamma(0.0), f64::INFINITY);
        assert_eq!(gamma(f64::INFINITY), f64::INFINITY);
        assert!(gamma(f64::NEG_INFINITY).is_nan());
        assert!(gamma(f64::NAN).is_nan());
    }

    #[test]
    fn test_gammaln() {
        assert_eq!(gammaln(1.0), 0.0);
        assert_eq!(gammaln(2.0), 0.0);
        assert_eq!(gammaln(3.0), LN_2);
        assert_eq!(gammaln(4.0), LN_2 + LN_3);

        assert_approx_eq!(f64, gammaln(0.5), LN_PI_HALF, ulps = 1);
        assert_approx_eq!(f64, gammaln(-0.5), LN_2 + LN_PI_HALF, ulps = 1);

        assert_eq!(gammaln(0.0), f64::INFINITY);
        assert_eq!(gammaln(f64::INFINITY), f64::INFINITY);
        assert_eq!(gammaln(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert!(gammaln(f64::NAN).is_nan());
    }

    #[test]
    fn test_gammasgn() {
        // Test positive gamma values (sign should be 1.0)
        assert_eq!(gammasgn(1.0), 1.0);
        assert_eq!(gammasgn(2.0), 1.0);
        assert_eq!(gammasgn(0.5), 1.0);

        // Test negative gamma values
        // gamma(-0.5) = -2*sqrt(pi), so sign should be -1.0
        assert_eq!(gammasgn(-0.5), -1.0);
        // gamma(-1.5) = 4*sqrt(pi)/3, so sign should be 1.0
        assert_eq!(gammasgn(-1.5), 1.0);

        // Special values
        assert_eq!(gammasgn(f64::INFINITY), 1.0);
        assert!(gammasgn(f64::NAN).is_nan());
    }

    #[test]
    fn test_gammainc() {
        // Lower incomplete gamma function: gammainc(a, x) = γ(a, x) / Γ(a)

        // Test basic values
        assert_approx_eq!(f64, gammainc(1.0, 0.0), 0.0, ulps = 1);
        assert_approx_eq!(f64, gammainc(1.0, 1.0), 1.0 - (-1.0f64).exp(), ulps = 2);

        // Test with different a values
        assert_approx_eq!(
            f64,
            gammainc(2.0, 1.0),
            1.0 - 2.0 * (-1.0f64).exp(),
            ulps = 2
        );

        // Test edge cases
        assert_approx_eq!(f64, gammainc(1.0, f64::INFINITY), 1.0, ulps = 1);
        assert!(gammainc(f64::NAN, 1.0).is_nan());
        assert!(gammainc(1.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_gammaincc() {
        // Upper incomplete gamma function: gammaincc(a, x) = Γ(a, x) / Γ(a)

        // Test basic values
        assert_approx_eq!(f64, gammaincc(1.0, 0.0), 1.0, ulps = 1);
        assert_approx_eq!(f64, gammaincc(1.0, 1.0), (-1.0f64).exp(), ulps = 2);

        // Test relationship: gammainc + gammaincc = 1
        let a = 2.5;
        let x = 1.5;
        assert_approx_eq!(f64, gammainc(a, x) + gammaincc(a, x), 1.0, ulps = 2);

        // Test edge cases
        assert_approx_eq!(f64, gammaincc(1.0, f64::INFINITY), 0.0, ulps = 1);
        assert!(gammaincc(f64::NAN, 1.0).is_nan());
        assert!(gammaincc(1.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_gammaincinv() {
        // Inverse of lower incomplete gamma function

        // Test that gammaincinv is inverse of gammainc
        let a = 2.0;
        let p = 0.5;
        let x = gammaincinv(a, p);
        assert_approx_eq!(f64, gammainc(a, x), p, ulps = 3);

        // Test boundary values
        assert_eq!(gammaincinv(1.0, 0.0), 0.0);
        assert_eq!(gammaincinv(1.0, 1.0), f64::INFINITY);

        // Test edge cases
        assert!(gammaincinv(f64::NAN, 0.5).is_nan());
        assert!(gammaincinv(1.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_gammainccinv() {
        // Inverse of upper incomplete gamma function

        // Test that gammainccinv is inverse of gammaincc
        let a = 2.0;
        let q = 0.3;
        let x = gammainccinv(a, q);
        assert_approx_eq!(f64, gammaincc(a, x), q, ulps = 3);

        // Test boundary values
        assert_eq!(gammainccinv(1.0, 1.0), 0.0);
        assert_eq!(gammainccinv(1.0, 0.0), f64::INFINITY);

        // Test edge cases
        assert!(gammainccinv(f64::NAN, 0.5).is_nan());
        assert!(gammainccinv(1.0, f64::NAN).is_nan());
    }

    #[test]
    fn test_gamma_ratio() {
        // Test gamma_ratio(a, b) = Γ(a) / Γ(b)

        // Test basic ratios
        assert_approx_eq!(f64, gamma_ratio(2.0, 1.0), 1.0, ulps = 1); // Γ(2)/Γ(1) = 1/1 = 1
        assert_approx_eq!(f64, gamma_ratio(3.0, 2.0), 2.0, ulps = 1); // Γ(3)/Γ(2) = 2/1 = 2
        assert_approx_eq!(f64, gamma_ratio(4.0, 3.0), 3.0, ulps = 1); // Γ(4)/Γ(3) = 6/2 = 3

        // Test identity
        assert_approx_eq!(f64, gamma_ratio(5.0, 5.0), 1.0, ulps = 1);

        // Test with fractional values
        // Γ(1.5) = sqrt(π)/2, Γ(0.5) = sqrt(π), so Γ(1.5)/Γ(0.5) = 0.5
        assert_approx_eq!(f64, gamma_ratio(1.5, 0.5), 0.5, ulps = 2);

        // Test edge cases
        assert!(gamma_ratio(f64::NAN, 1.0).is_nan());
        assert!(gamma_ratio(1.0, f64::NAN).is_nan());
        assert_eq!(gamma_ratio(1.0, f64::INFINITY), 0.0);
        assert_eq!(gamma_ratio(f64::INFINITY, 1.0), f64::INFINITY);
    }
}
