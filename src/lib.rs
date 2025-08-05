use num_traits::Float;

mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

macro_rules! xsf_function {
    ($fn_name:ident, 1, $(($rust_type:ty, $ffi_func:ident)),+) => {
        paste::paste! {
            pub trait [<$fn_name:camel>]<T: Float> {
                fn [<xsf_ $fn_name>](x: T) -> T;
            }

            $(
                impl [<$fn_name:camel>]<$rust_type> for $rust_type {
                    fn [<xsf_ $fn_name>](x: $rust_type) -> $rust_type {
                        unsafe { ffi::$ffi_func(x) }
                    }
                }
            )+

            pub fn $fn_name<T>(x: T) -> T
            where
                T: [<$fn_name:camel>]<T> + Float,
            {
                T::[<xsf_ $fn_name>](x)
            }
        }
    };
    ($fn_name:ident, 2, $(($rust_type:ty, $ffi_func:ident)),+) => {
        paste::paste! {
            pub trait [<$fn_name:camel>]<T: Float> {
                fn [<xsf_ $fn_name>](a: T, x: T) -> T;
            }

            $(
                impl [<$fn_name:camel>]<$rust_type> for $rust_type {
                    fn [<xsf_ $fn_name>](a: $rust_type, x: $rust_type) -> $rust_type {
                        unsafe { ffi::$ffi_func(a, x) }
                    }
                }
            )+

            pub fn $fn_name<T>(a: T, x: T) -> T
            where
                T: [<$fn_name:camel>]<T> + Float,
            {
                T::[<xsf_ $fn_name>](a, x)
            }
        }
    };
}

xsf_function!(gamma, 1, (f32, gamma_f), (f64, gamma_d));
xsf_function!(gammaln, 1, (f32, gammaln_f), (f64, gammaln_d));
xsf_function!(gammasgn, 1, (f32, gammasgn_f), (f64, gammasgn_d));
xsf_function!(gammainc, 2, (f32, gammainc_f), (f64, gammainc_d));
xsf_function!(gammaincinv, 2, (f32, gammaincinv_f), (f64, gammaincinv_d));
xsf_function!(gammaincc, 2, (f32, gammaincc_f), (f64, gammaincc_d));
xsf_function!(
    gammainccinv,
    2,
    (f32, gammainccinv_f),
    (f64, gammainccinv_d)
);
xsf_function!(gamma_ratio, 2, (f32, gamma_ratio_f), (f64, gamma_ratio_d));

#[cfg(test)]
mod tests {
    use std::f64::consts::LN_2;

    use super::*;
    use float_cmp::assert_approx_eq;

    const LN_3: f64 = 1.098_612_288_668_109_7;
    const SQRT_PI: f64 = 1.772_453_850_905_516;
    const LN_PI_HALF: f64 = 0.572_364_942_924_700;

    #[test]
    fn test_gamma() {
        assert_eq!(gamma(1.0f32), 1.0);
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
        assert_eq!(gammaln(1.0f32), 0.0);
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
        assert_eq!(gammasgn(1.0f32), 1.0);
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

        // Test f32 version
        assert_approx_eq!(f32, gammainc(1.0f32, 0.0f32), 0.0, ulps = 1);

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

        // Test f32 version
        assert_approx_eq!(f32, gammaincc(1.0f32, 0.0f32), 1.0, ulps = 1);

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

        // Test f32 version
        let x_f32 = gammaincinv(2.0f32, 0.5f32);
        assert_approx_eq!(f32, gammainc(2.0f32, x_f32), 0.5, ulps = 3);

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

        // Test f32 version
        let x_f32 = gammainccinv(2.0f32, 0.3f32);
        assert_approx_eq!(f32, gammaincc(2.0f32, x_f32), 0.3, ulps = 3);

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

        // Test f32 version
        assert_approx_eq!(f32, gamma_ratio(3.0f32, 2.0f32), 2.0, ulps = 1);

        // Test edge cases
        assert!(gamma_ratio(f64::NAN, 1.0).is_nan());
        assert!(gamma_ratio(1.0, f64::NAN).is_nan());
        assert_eq!(gamma_ratio(1.0, f64::INFINITY), 0.0);
        assert_eq!(gamma_ratio(f64::INFINITY, 1.0), f64::INFINITY);
    }
}
