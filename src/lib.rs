use std::ffi::c_int;

mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

    unsafe extern "C" {
        pub fn cbrt(x: f64) -> f64;
        pub fn erf(x: f64) -> f64;
        pub fn erfc(x: f64) -> f64;
        pub fn tgamma(x0: f64) -> f64;
        pub fn expm1(x: f64) -> f64;
        pub fn exp2(x: f64) -> f64;
        pub fn exp10(x: f64) -> f64;
    }
}

macro_rules! xsf_impl {
    ($name:ident, ($($param:ident: $type:ty),*), $docs:expr) => {
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { ffi::$name($($param),*) }
        }
    };
}

// alg.h
/// Cube root
pub fn cbrt(x: f64) -> f64 {
    unsafe { ffi::cbrt(x) }
}

// bessel.h
xsf_impl!(cyl_bessel_j, (v: f64, x: f64), "Bessel function, 1st kind");
xsf_impl!(cyl_bessel_je, (v: f64, x: f64), "Exponentially scaled Bessel function, 1st kind");
xsf_impl!(cyl_bessel_y, (v: f64, x: f64), "Bessel function, 2nd kind");
xsf_impl!(cyl_bessel_ye, (v: f64, x: f64), "Exponentially scaled Bessel function, 2nd kind");
xsf_impl!(cyl_bessel_i, (v: f64, x: f64), "Modified Bessel function, 1st kind");
xsf_impl!(
    cyl_bessel_ie,
    (v: f64, x: f64),
    "Exponentially scaled modified Bessel function, 1st kind"
);
xsf_impl!(cyl_bessel_k, (v: f64, x: f64), "Modified Bessel function, 2nd kind");
xsf_impl!(
    cyl_bessel_ke,
    (v: f64, x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind"
);
xsf_impl!(cyl_bessel_j0, (x: f64), "Bessel function, 1st kind, order 0");
xsf_impl!(cyl_bessel_j1, (x: f64), "Bessel function, 1st kind, order 1");
xsf_impl!(cyl_bessel_y0, (x: f64), "Bessel function, 2nd kind, order 0");
xsf_impl!(cyl_bessel_y1, (x: f64), "Bessel function, 2nd kind, order 1");
xsf_impl!(cyl_bessel_i0, (x: f64), "Modified Bessel function, 1st kind, order 0");
xsf_impl!(
    cyl_bessel_i0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 0"
);
xsf_impl!(cyl_bessel_i1, (x: f64), "Modified Bessel function, 1st kind, order 1");
xsf_impl!(
    cyl_bessel_i1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 1"
);
xsf_impl!(cyl_bessel_k0, (x: f64), "Modified Bessel function, 2nd kind, order 0");
xsf_impl!(
    cyl_bessel_k0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 0"
);
xsf_impl!(cyl_bessel_k1, (x: f64), "Modified Bessel function, 2nd kind, order 1");
xsf_impl!(
    cyl_bessel_k1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 1"
);
xsf_impl!(
    besselpoly,
    (a: f64, lambda: f64, nu: f64),
    "Weighted integral of the Bessel function of the first kind"
);

// beta.h
xsf_impl!(beta, (a: f64, b: f64), "Beta function");
xsf_impl!(betaln, (a: f64, b: f64), "Natural log of `|beta|`");

// binom.h
xsf_impl!(binom, (n: f64, k: f64), "Binomial coefficient");

// digamma.h
xsf_impl!(digamma, (x: f64), "Digamma function");

// erf.h
/// Error function
pub fn erf(x: f64) -> f64 {
    unsafe { ffi::erf(x) }
}
/// Complementary error function `1 - erf(x)`
pub fn erfc(x: f64) -> f64 {
    unsafe { ffi::erfc(x) }
}
xsf_impl!(erfcx, (x: f64), "Scaled complementary error function `exp(x^2) * erfc(x)`");
xsf_impl!(erfi, (x: f64), "Imaginary error function `-i erf(ix)`");
xsf_impl!(voigt_profile, (x: f64, sigma: f64, gamma: f64), "Voigt profile");
xsf_impl!(dawsn, (x: f64), "Dawson function `sqrt(pi)/2 * exp(-x^2) * erfi(x)`");

// exp.h
/// `exp(x) - 1`
pub fn expm1(x: f64) -> f64 {
    unsafe { ffi::expm1(x) }
}
/// `2^x`
pub fn exp2(x: f64) -> f64 {
    unsafe { ffi::exp2(x) }
}
/// `10^x`
pub fn exp10(x: f64) -> f64 {
    unsafe { ffi::exp10(x) }
}

// expint.h
xsf_impl!(exp1, (x: f64), "Exponential integral `E_1(x)`");
xsf_impl!(expi, (x: f64), "Exponential integral `E_i(x)`");
xsf_impl!(scaled_exp1, (x: f64), "Scaled version of the exponential integral `E_1(x)`");

// gamma.h
/// Gamma function
pub fn gamma(x: f64) -> f64 {
    unsafe { ffi::tgamma(x) }
}
xsf_impl!(gamma_ratio, (a: f64, b: f64), "`gamma(a) / gamma(b)`");
xsf_impl!(gammainc, (a: f64, x: f64), "Incomplete Gamma integral");
xsf_impl!(gammaincc, (a: f64, x: f64), "Complemented incomplete Gamma integral");
xsf_impl!(gammaincinv, (a: f64, p: f64), "Inverse of `gammainc`");
xsf_impl!(gammainccinv, (a: f64, p: f64), "Inverse of `gammaincc`");
xsf_impl!(gammaln, (x: f64), "Natural logarithm of Gamma function");
xsf_impl!(gammasgn, (x: f64), "Sign of the Gamma function");

// hyp2f1.h
xsf_impl!(hyp2f1, (a: f64, b: f64, c: f64, x: f64), "Gauss hypergeometric function `2F1`");

// iv_ratio.h
xsf_impl!(
    iv_ratio,
    (v: f64, x: f64),
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function, 1st kind"
);
xsf_impl!(
    iv_ratio_c,
    (v: f64, x: f64),
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function, 1st kind"
);

// kelvin.h
xsf_impl!(ber, (x: f64), "Kelvin function `ber`");
xsf_impl!(bei, (x: f64), "Kelvin function `bei`");
xsf_impl!(ker, (x: f64), "Kelvin function `ker`");
xsf_impl!(kei, (x: f64), "Kelvin function `kei`");
xsf_impl!(berp, (x: f64), "Derivative of the Kelvin function `ber`");
xsf_impl!(beip, (x: f64), "Derivative of the Kelvin function `bei`");
xsf_impl!(kerp, (x: f64), "Derivative of the Kelvin function `ker`");
xsf_impl!(keip, (x: f64), "Derivative of the Kelvin function `kei`");

// legendre.h
xsf_impl!(legendre_p, (n: c_int, x: f64), "Legendre polynomial of degree n");
xsf_impl!(
    sph_legendre_p,
    (n: c_int, m: c_int, theta: f64),
    "Spherical Legendre polynomial of degree n and order m"
);

// log_exp.h
xsf_impl!(expit, (x: f64), "Expit function, `1/(1 + exp(-x))`");
xsf_impl!(exprel, (x: f64), "Relative error exponential, `(exp(x) - 1)/x`");
xsf_impl!(logit, (x: f64), "Logit function, `log(x / (1 - x))`");
xsf_impl!(log_expit, (x: f64), "Log of the expit function, `log(expit(x))`");
xsf_impl!(log1mexp, (x: f64), "Compute `log(1 - exp(x))`");

// log.h
xsf_impl!(log1pmx, (x: f64), "Compute `log(1 + x) - x`");
xsf_impl!(xlogy, (x: f64, y: f64), "Compute `x * log(y)`");
xsf_impl!(xlog1py, (x: f64, y: f64), "Compute `x * log(1 + y)`");

// loggamma.h
xsf_impl!(loggamma, (x: f64), "Principal branch of the logarithm of `gamma(x)`");
xsf_impl!(rgamma, (x: f64), "Reciprocal Gamma function");

// mathieu.h
xsf_impl!(cem_cva, (m: f64, q: f64), "Characteristic value of even Mathieu functions");
xsf_impl!(sem_cva, (m: f64, q: f64), "Characteristic value of odd Mathieu functions");

// specfun.h
xsf_impl!(hypu, (a: f64, b: f64, x: f64), "Confluent hypergeometric function `U`");
xsf_impl!(hyp1f1, (a: f64, b: f64, x: f64), "Confluent hypergeometric function `1F1`");
xsf_impl!(pmv, (m: f64, v: f64, x: f64), "Associated Legendre function");

// sphd_wave.h
xsf_impl!(
    prolate_segv,
    (m: f64, n: f64, c: f64),
    "Characteristic value of prolate spheroidal function"
);
xsf_impl!(
    oblate_segv,
    (m: f64, n: f64, c: f64),
    "Characteristic value of oblate spheroidal function"
);

// stats.h
xsf_impl!(bdtr, (k: f64, n: c_int, p: f64), "Binomial distribution function");
xsf_impl!(bdtrc, (k: f64, n: c_int, p: f64), "Binomial survival function");
xsf_impl!(bdtri, (k: f64, n: c_int, y: f64), "Binomial quantile function");
xsf_impl!(chdtr, (df: f64, x: f64), "Chi-squared distribution function");
xsf_impl!(chdtrc, (df: f64, x: f64), "Chi-squared survival function");
xsf_impl!(chdtri, (df: f64, y: f64), "Chi-squared quantile function");
xsf_impl!(fdtr, (a: f64, b: f64, x: f64), "F distribution function");
xsf_impl!(fdtrc, (a: f64, b: f64, x: f64), "F survival function");
xsf_impl!(fdtri, (a: f64, b: f64, y: f64), "F quantile function");
xsf_impl!(gdtr, (a: f64, b: f64, x: f64), "Gamma distribution function");
xsf_impl!(gdtrc, (a: f64, b: f64, x: f64), "Gamma survival function");
xsf_impl!(kolmogorov, (x: f64), "Kolmogorov survival function");
xsf_impl!(kolmogc, (x: f64), "Kolmogorov distribution function");
xsf_impl!(kolmogi, (x: f64), "Inverse of `kolmogorov`");
xsf_impl!(kolmogp, (x: f64), "Derivative of `kolmogorov`");
xsf_impl!(ndtr, (x: f64), "Normal distribution function");
xsf_impl!(ndtri, (x: f64), "Normal quantile function");
xsf_impl!(log_ndtr, (x: f64), "Log of the normal distribution function");
xsf_impl!(nbdtr, (k: c_int, n: c_int, p: f64), "Negative binomial distribution function");
xsf_impl!(nbdtrc, (k: c_int, n: c_int, p: f64), "Negative binomial survival function");
xsf_impl!(nbdtri, (k: c_int, n: c_int, p: f64), "Negative binomial quantile function");
xsf_impl!(owens_t, (h: f64, a: f64), "Owen's T function");
xsf_impl!(pdtr, (k: f64, m: f64), "Poisson distribution function");
xsf_impl!(pdtrc, (k: f64, m: f64), "Poisson survival function");
xsf_impl!(pdtri, (k: c_int, y: f64), "Poisson quantile function");
xsf_impl!(smirnov, (n: c_int, x: f64), "Kolmogorov-Smirnov survival function");
xsf_impl!(smirnovc, (n: c_int, x: f64), "Kolmogorov-Smirnov distribution function");
xsf_impl!(smirnovi, (n: c_int, x: f64), "Inverse of `smirnov`");
xsf_impl!(smirnovp, (n: c_int, x: f64), "Derivative of `smirnov`");
xsf_impl!(tukeylambdacdf, (x: f64, lmbda: f64), "Tukey-Lambda distribution function");

// struve.h
xsf_impl!(itstruve0, (x: f64), "Integral of the Struve function of order 0");
xsf_impl!(it2struve0, (x: f64), "Integral related to the Struve function of order 0");
xsf_impl!(itmodstruve0, (x: f64), "Integral of the modified Struve function of order 0");
xsf_impl!(struve_h, (v: f64, x: f64), "Struve `H` function");
xsf_impl!(struve_l, (v: f64, x: f64), "Struve `L` function");

// trig.h
xsf_impl!(sinpi, (x: f64), "Compute `sin(pi * x)`");
xsf_impl!(cospi, (x: f64), "Compute `cos(pi * x)`");
xsf_impl!(sindg, (x: f64), "Circular sine of angle in degrees");
xsf_impl!(cosdg, (x: f64), "Circular cosine of angle in degrees");
xsf_impl!(tandg, (x: f64), "Circular tangent of argument in degrees");
xsf_impl!(cotdg, (x: f64), "Circular cotangent of argument in degrees");
xsf_impl!(radian, (d: f64, m: f64, s: f64), "Degrees, minutes, seconds to radians");
xsf_impl!(cosm1, (x: f64), "Compute `cos(x) - 1`");

// wright_bessel.h
xsf_impl!(wright_bessel, (a: f64, b: f64, x: f64), "Wright's generalized Bessel function");
xsf_impl!(log_wright_bessel, (a: f64, b: f64, x: f64), "Logarithm of `wright_bessel`");

// zeta.h
xsf_impl!(riemann_zeta, (x: f64), "Riemann zeta function");
xsf_impl!(zeta, (x: f64, q: f64), "Riemann zeta function of two arguments");
xsf_impl!(zetac, (x: f64), "Riemann zeta function, minus one");

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
        assert_eq!(gamma(1.0), 1.0);
        assert_eq!(gamma(2.0), 1.0);
        assert_eq!(gamma(3.0), 2.0);
        assert_eq!(gamma(5.0), 24.0);
        assert_approx_eq!(f64, gamma(8.0), 5_040.0, ulps = 1);

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
