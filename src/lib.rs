mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

pub mod xsf64 {
    use crate::ffi;
    use std::ffi::c_int;

    macro_rules! xsf_impl {
        ($name:ident, $($param:ident: $type:ty),*) => {
            paste::paste! {
                pub fn $name($($param: $type),*) -> f64 {
                    unsafe { ffi::[<xsf_ $name>]($($param),*) }
                }
            }
        };
    }

    // alg.h
    xsf_impl!(cbrt, x: f64);
    // beta.h
    xsf_impl!(beta, a: f64, b: f64);
    xsf_impl!(betaln, a: f64, b: f64);
    // binom.h
    xsf_impl!(binom, n: f64, k: f64);
    // digamma.h
    xsf_impl!(digamma, z: f64);
    // erf.h
    xsf_impl!(erf, x: f64);
    xsf_impl!(erfc, x: f64);
    xsf_impl!(erfcx, x: f64);
    xsf_impl!(erfi, x: f64);
    xsf_impl!(voigt_profile, x: f64, sigma: f64, gamma: f64);
    xsf_impl!(dawsn, x: f64);
    // exp.h
    xsf_impl!(expm1, x: f64);
    xsf_impl!(exp2, x: f64);
    xsf_impl!(exp10, x: f64);
    // expint.h
    xsf_impl!(exp1, x: f64);
    xsf_impl!(expi, x: f64);
    xsf_impl!(scaled_exp1, x: f64);
    // fp_error_metrics.h
    xsf_impl!(extended_absolute_error, actual: f64, desired: f64);
    xsf_impl!(extended_relative_error, actual: f64, desired: f64);
    // gamma.h
    xsf_impl!(gamma, x: f64);
    xsf_impl!(gammaln, x: f64);
    xsf_impl!(gammasgn, x: f64);
    xsf_impl!(gammainc, a: f64, x: f64);
    xsf_impl!(gammaincinv, a: f64, p: f64);
    xsf_impl!(gammaincc, a: f64, x: f64);
    xsf_impl!(gammainccinv, a: f64, p: f64);
    xsf_impl!(gamma_ratio, a: f64, b: f64);
    // hyp2f1.h
    xsf_impl!(hyp2f1, a: f64, b: f64, c: f64, x: f64);
    // iv_ratio.h
    xsf_impl!(iv_ratio, v: f64, x: f64);
    xsf_impl!(iv_ratio_c, v: f64, x: f64);
    // kelvin.h
    xsf_impl!(ber, x: f64);
    xsf_impl!(bei, x: f64);
    xsf_impl!(ker, x: f64);
    xsf_impl!(kei, x: f64);
    xsf_impl!(berp, x: f64);
    xsf_impl!(beip, x: f64);
    xsf_impl!(kerp, x: f64);
    xsf_impl!(keip, x: f64);
    // legendre.h
    xsf_impl!(legendre_p, n: c_int, z: f64);
    xsf_impl!(sph_legendre_p, n: c_int, m: c_int, theta: f64);
    // log_exp.h
    xsf_impl!(expit, x: f64);
    xsf_impl!(exprel, x: f64);
    xsf_impl!(logit, x: f64);
    xsf_impl!(log_expit, x: f64);
    xsf_impl!(log1mexp, x: f64);
    // log.h
    xsf_impl!(log1p, x: f64);
    xsf_impl!(log1pmx, x: f64);
    xsf_impl!(xlogy, x: f64, y: f64);
    xsf_impl!(xlog1py, x: f64, y: f64);
    // loggamma.h
    xsf_impl!(loggamma, x: f64);
    xsf_impl!(rgamma, z: f64);
    // mathieu.h
    xsf_impl!(cem_cva, m: f64, q: f64);
    xsf_impl!(sem_cva, m: f64, q: f64);
    // specfun.h
    xsf_impl!(hypu, a: f64, b: f64, x: f64);
    xsf_impl!(hyp1f1, a: f64, b: f64, x: f64);
    xsf_impl!(pmv, m: f64, v: f64, x: f64);
    // sphd_wave.h
    xsf_impl!(prolate_segv, m: f64, n: f64, c: f64);
    xsf_impl!(oblate_segv, m: f64, n: f64, c: f64);
    // stats.h
    xsf_impl!(bdtr, k: f64, n: c_int, p: f64);
    xsf_impl!(bdtrc, k: f64, n: c_int, p: f64);
    xsf_impl!(bdtri, k: f64, n: c_int, y: f64);
    xsf_impl!(chdtr, df: f64, x: f64);
    xsf_impl!(chdtrc, df: f64, x: f64);
    xsf_impl!(chdtri, df: f64, y: f64);
    xsf_impl!(fdtr, a: f64, b: f64, x: f64);
    xsf_impl!(fdtrc, a: f64, b: f64, x: f64);
    xsf_impl!(fdtri, a: f64, b: f64, y: f64);
    xsf_impl!(gdtr, a: f64, b: f64, x: f64);
    xsf_impl!(gdtrc, a: f64, b: f64, x: f64);
    xsf_impl!(kolmogorov, x: f64);
    xsf_impl!(kolmogc, x: f64);
    xsf_impl!(kolmogi, x: f64);
    xsf_impl!(kolmogp, x: f64);
    xsf_impl!(ndtr, x: f64);
    xsf_impl!(ndtri, x: f64);
    xsf_impl!(log_ndtr, x: f64);
    xsf_impl!(nbdtr, k: c_int, n: c_int, p: f64);
    xsf_impl!(nbdtrc, k: c_int, n: c_int, p: f64);
    xsf_impl!(nbdtri, k: c_int, n: c_int, p: f64);
    xsf_impl!(owens_t, h: f64, a: f64);
    xsf_impl!(pdtr, k: f64, m: f64);
    xsf_impl!(pdtrc, k: f64, m: f64);
    xsf_impl!(pdtri, k: c_int, y: f64);
    xsf_impl!(smirnov, n: c_int, x: f64);
    xsf_impl!(smirnovc, n: c_int, x: f64);
    xsf_impl!(smirnovi, n: c_int, x: f64);
    xsf_impl!(smirnovp, n: c_int, x: f64);
    xsf_impl!(tukeylambdacdf, x: f64, lmbda: f64);
    // struve.h
    xsf_impl!(itstruve0, x: f64);
    xsf_impl!(it2struve0, x: f64);
    xsf_impl!(itmodstruve0, x: f64);
    xsf_impl!(struve_h, v: f64, z: f64);
    xsf_impl!(struve_l, v: f64, z: f64);
    // trig.h
    xsf_impl!(sinpi, x: f64);
    xsf_impl!(cospi, x: f64);
    xsf_impl!(sindg, x: f64);
    xsf_impl!(cosdg, x: f64);
    xsf_impl!(tandg, x: f64);
    xsf_impl!(cotdg, x: f64);
    xsf_impl!(radian, d: f64, m: f64, s: f64);
    xsf_impl!(cosm1, x: f64);
    // wright_bessel.h
    // xsf_impl!(wright_bessel_t, a: f64, b: f64, x: f64);
    xsf_impl!(wright_bessel, a: f64, b: f64, x: f64);
    xsf_impl!(log_wright_bessel, a: f64, b: f64, x: f64);
    // zeta.h
    xsf_impl!(riemann_zeta, s: f64);
    xsf_impl!(zeta, s: f64, q: f64);
    xsf_impl!(zetac, s: f64);
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
