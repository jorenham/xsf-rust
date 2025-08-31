use std::ffi::c_int;

mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

macro_rules! xsf_impl {
    ($name:ident, $xsf_name:ident, $($param:ident: $type:ty),*) => {
        xsf_impl!($name, $xsf_name, "", $($param: $type),*);
    };
    ($name:ident, $xsf_name:ident, $docs:expr, $($param:ident: $type:ty),*) => {
        #[allow(clippy::empty_docs)]
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { ffi::$xsf_name($($param),*) }
        }
    };
}

// alg.h
xsf_impl!(cbrt, xsf_cbrt, "Cube root", x: f64);

// bessel.h
xsf_impl!(cyl_bessel_j, xsf_cyl_bessel_j, "Bessel function, 1st kind", v: f64, x: f64);
xsf_impl!(
    cyl_bessel_je,
    xsf_cyl_bessel_je,
    "Exponentially scaled Bessel function, 1st kind",
    v: f64,
    x: f64
);
xsf_impl!(cyl_bessel_y, xsf_cyl_bessel_y, "Bessel function, 2nd kind", v: f64, x: f64);
xsf_impl!(
    cyl_bessel_ye,
    xsf_cyl_bessel_ye,
    "Exponentially scaled Bessel function, 2nd kind",
    v: f64,
    x: f64
);
xsf_impl!(cyl_bessel_i, xsf_cyl_bessel_i, "Modified Bessel function, 1st kind", v: f64, x: f64);
xsf_impl!(
    cyl_bessel_ie,
    xsf_cyl_bessel_ie,
    "Exponentially scaled modified Bessel function, 1st kind",
    v: f64,
    x: f64
);
xsf_impl!(cyl_bessel_k, xsf_cyl_bessel_k, "Modified Bessel function, 2nd kind", v: f64, x: f64);
xsf_impl!(
    cyl_bessel_ke,
    xsf_cyl_bessel_ke,
    "Exponentially scaled modified Bessel function, 2nd kind",
    v: f64,
    x: f64
);

xsf_impl!(cyl_bessel_j0, xsf_cyl_bessel_j0, "Bessel function, 1st kind, order 0", x: f64);
xsf_impl!(cyl_bessel_j1, xsf_cyl_bessel_j1, "Bessel function, 1st kind, order 1", x: f64);
xsf_impl!(cyl_bessel_y0, xsf_cyl_bessel_y0, "Bessel function, 2nd kind, order 0", x: f64);
xsf_impl!(cyl_bessel_y1, xsf_cyl_bessel_y1, "Bessel function, 2nd kind, order 1", x: f64);
xsf_impl!(cyl_bessel_i0, xsf_cyl_bessel_i0, "Modified Bessel function, 1st kind, order 0", x: f64);
xsf_impl!(
    cyl_bessel_i0e,
    xsf_cyl_bessel_i0e,
    "Exponentially scaled modified Bessel function, 1st kind, order 0",
    x: f64
);
xsf_impl!(cyl_bessel_i1, xsf_cyl_bessel_i1, "Modified Bessel function, 1st kind, order 1", x: f64);
xsf_impl!(
    cyl_bessel_i1e,
    xsf_cyl_bessel_i1e,
    "Exponentially scaled modified Bessel function, 1st kind, order 1",
    x: f64
);
xsf_impl!(cyl_bessel_k0, xsf_cyl_bessel_k0, "Modified Bessel function, 2nd kind, order 0", x: f64);
xsf_impl!(
    cyl_bessel_k0e,
    xsf_cyl_bessel_k0e,
    "Exponentially scaled modified Bessel function, 2nd kind, order 0",
    x: f64
);
xsf_impl!(cyl_bessel_k1, xsf_cyl_bessel_k1, "Modified Bessel function, 2nd kind, order 1", x: f64);
xsf_impl!(
    cyl_bessel_k1e,
    xsf_cyl_bessel_k1e,
    "Exponentially scaled modified Bessel function, 2nd kind, order 1",
    x: f64
);

xsf_impl!(
    besselpoly,
    xsf_besselpoly,
    "Weighted integral of the Bessel function of the first kind",
    a: f64,
    lambda: f64,
    nu: f64
);

// beta.h
xsf_impl!(beta, xsf_beta, "Beta function", a: f64, b: f64);
xsf_impl!(betaln, xsf_betaln, "Natural log of `|beta|`", a: f64, b: f64);

// binom.h
xsf_impl!(binom, xsf_binom, "Binomial coefficient", n: f64, k: f64);

// digamma.h
xsf_impl!(digamma, xsf_digamma, "Digamma function", x: f64);

// erf.h
xsf_impl!(erf, xsf_erf, "Error function", x: f64);
xsf_impl!(erfc, xsf_erfc, "Complementary error function `1 - erf(x)`", x: f64);
xsf_impl!(erfcx, xsf_erfcx, "Scaled complementary error function `exp(x^2) * erfc(x)`", x: f64);
xsf_impl!(erfi, xsf_erfi, "Imaginary error function `-i erf(ix)`", x: f64);
xsf_impl!(voigt_profile, xsf_voigt_profile, "Voigt profile", x: f64, sigma: f64, gamma: f64);
xsf_impl!(dawsn, xsf_dawsn, "Dawson function `sqrt(pi)/2 * exp(-x^2) * erfi(x)`", x: f64);

// exp.h
xsf_impl!(expm1, xsf_expm1, "`exp(x) - 1`", x: f64);
xsf_impl!(exp2, xsf_exp2, "`2^x`", x: f64);
xsf_impl!(exp10, xsf_exp10, "`10^x`", x: f64);

// expint.h
xsf_impl!(exp1, xsf_exp1, "Exponential integral `E_1(x)`", x: f64);
xsf_impl!(expi, xsf_expi, "Exponential integral `E_i(x)`", x: f64);
xsf_impl!(
    scaled_exp1,
    xsf_scaled_exp1,
    "Scaled version of the exponential integral `E_1(x)`",
    x: f64
);

// fp_error_metrics.h
// xsf_impl!(extended_absolute_error, xsf_extended_absolute_error, actual: f64, desired: f64);
// xsf_impl!(extended_relative_error, xsf_extended_relative_error, actual: f64, desired: f64);

// gamma.h
xsf_impl!(gamma, xsf_gamma, "Gamma function", x: f64);
xsf_impl!(gamma_ratio, xsf_gamma_ratio, "`gamma(a) / gamma(b)`", a: f64, b: f64);
xsf_impl!(gammainc, xsf_gammainc, "Incomplete Gamma integral", a: f64, x: f64);
xsf_impl!(gammaincc, xsf_gammaincc, "Complemented incomplete Gamma integral", a: f64, x: f64);
xsf_impl!(gammaincinv, xsf_gammaincinv, "Inverse of `gammainc`", a: f64, p: f64);
xsf_impl!(gammainccinv, xsf_gammainccinv, "Inverse of `gammaincc`", a: f64, p: f64);
xsf_impl!(gammaln, xsf_gammaln, "Natural logarithm of Gamma function", x: f64);
xsf_impl!(gammasgn, xsf_gammasgn, "Sign of the Gamma function", x: f64);

// hyp2f1.h
xsf_impl!(hyp2f1, xsf_hyp2f1, "Gauss hypergeometric function `2F1`", a: f64, b: f64, c: f64, x: f64);

// iv_ratio.h
xsf_impl!(
    iv_ratio,
    xsf_iv_ratio,
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function, 1st kind",
    v: f64,
    x: f64
);
xsf_impl!(
    iv_ratio_c,
    xsf_iv_ratio_c,
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function, 1st kind",
    v: f64,
    x: f64
);

// kelvin.h
xsf_impl!(ber, xsf_ber, "Kelvin function `ber`", x: f64);
xsf_impl!(bei, xsf_bei,  "Kelvin function `bei`", x: f64);
xsf_impl!(ker, xsf_ker, "Kelvin function `ker`", x: f64);
xsf_impl!(kei, xsf_kei, "Kelvin function `kei`", x: f64);
xsf_impl!(berp, xsf_berp, "Derivative of the Kelvin function `ber`", x: f64);
xsf_impl!(beip, xsf_beip, "Derivative of the Kelvin function `bei`", x: f64);
xsf_impl!(kerp, xsf_kerp, "Derivative of the Kelvin function `ker`", x: f64);
xsf_impl!(keip, xsf_keip, "Derivative of the Kelvin function `kei`", x: f64);

// legendre.h
xsf_impl!(legendre_p, xsf_legendre_p, "Legendre polynomial of degree n", n: c_int, x: f64);
xsf_impl!(
    sph_legendre_p,
    xsf_sph_legendre_p,
    "Spherical Legendre polynomial of degree n and order m",
    n: c_int,
    m: c_int,
    theta: f64
);

// log_exp.h
xsf_impl!(expit, xsf_expit, "Expit function, `1/(1 + exp(-x))`", x: f64);
xsf_impl!(exprel, xsf_exprel, "Relative error exponential, `(exp(x) - 1)/x`", x: f64);
xsf_impl!(logit, xsf_logit, "Logit function, `log(x / (1 - x))`", x: f64);
xsf_impl!(log_expit, xsf_log_expit, "Log of the expit function, `log(expit(x))`", x: f64);
xsf_impl!(log1mexp, xsf_log1mexp, "Compute `log(1 - exp(x))`", x: f64);

// log.h
xsf_impl!(log1p, xsf_log1p, "Compute `log(1 + x)`", x: f64);
xsf_impl!(log1pmx, xsf_log1pmx, "Compute `log(1 + x) - x`", x: f64);
xsf_impl!(xlogy, xsf_xlogy, "Compute `x * log(y)`", x: f64, y: f64);
xsf_impl!(xlog1py, xsf_xlog1py, "Compute `x * log(1 + y)`", x: f64, y: f64);

// loggamma.h
xsf_impl!(loggamma, xsf_loggamma, "Principal branch of the logarithm of `gamma(x)`", x: f64);
xsf_impl!(rgamma, xsf_rgamma, "Reciprocal Gamma function", x: f64);

// mathieu.h
xsf_impl!(cem_cva, xsf_cem_cva, "Characteristic value of even Mathieu functions", m: f64, q: f64);
xsf_impl!(sem_cva, xsf_sem_cva, "Characteristic value of odd Mathieu functions", m: f64, q: f64);

// specfun.h
xsf_impl!(hypu, xsf_hypu, "Confluent hypergeometric function `U`", a: f64, b: f64, x: f64);
xsf_impl!(hyp1f1, xsf_hyp1f1, "Confluent hypergeometric function `1F1`", a: f64, b: f64, x: f64);
xsf_impl!(
    pmv,
    xsf_pmv,
    "Associated Legendre function of integer order and real degree",
    m: f64,
    v: f64,
    x: f64
);

// sphd_wave.h
xsf_impl!(
    prolate_segv,
    xsf_prolate_segv,
    "Characteristic value of prolate spheroidal function",
    m: f64,
    n: f64,
    c: f64
);
xsf_impl!(
    oblate_segv,
    xsf_oblate_segv,
    "Characteristic value of oblate spheroidal function",
    m: f64,
    n: f64,
    c: f64
);

// stats.h
xsf_impl!(bdtr, xsf_bdtr, "Binomial distribution function", k: f64, n: c_int, p: f64);
xsf_impl!(bdtrc, xsf_bdtrc, "Binomial survival function", k: f64, n: c_int, p: f64);
xsf_impl!(bdtri, xsf_bdtri, "Binomial quantile function", k: f64, n: c_int, y: f64);
xsf_impl!(chdtr, xsf_chdtr, "Chi-squared distribution function", df: f64, x: f64);
xsf_impl!(chdtrc, xsf_chdtrc, "Chi-squared survival function", df: f64, x: f64);
xsf_impl!(chdtri, xsf_chdtri, "Chi-squared quantile function", df: f64, y: f64);
xsf_impl!(fdtr, xsf_fdtr, "F distribution function", a: f64, b: f64, x: f64);
xsf_impl!(fdtrc, xsf_fdtrc, "F survival function", a: f64, b: f64, x: f64);
xsf_impl!(fdtri, xsf_fdtri, "F quantile function", a: f64, b: f64, y: f64);
xsf_impl!(gdtr, xsf_gdtr, "Gamma distribution function", a: f64, b: f64, x: f64);
xsf_impl!(gdtrc, xsf_gdtrc, "Gamma survival function", a: f64, b: f64, x: f64);
xsf_impl!(kolmogorov, xsf_kolmogorov, "Kolmogorov survival function", x: f64);
xsf_impl!(kolmogc, xsf_kolmogc, "Kolmogorov distribution function", x: f64);
xsf_impl!(kolmogi, xsf_kolmogi, "Inverse of `kolmogorov`", x: f64);
xsf_impl!(kolmogp, xsf_kolmogp, "Derivative of `kolmogorov`", x: f64);
xsf_impl!(ndtr, xsf_ndtr, "Normal distribution function", x: f64);
xsf_impl!(ndtri, xsf_ndtri, "Normal quantile function", x: f64);
xsf_impl!(log_ndtr, xsf_log_ndtr, "Log of the normal distribution function", x: f64);
xsf_impl!(nbdtr, xsf_nbdtr, "Negative binomial distribution function", k: c_int, n: c_int, p: f64);
xsf_impl!(nbdtrc, xsf_nbdtrc, "Negative binomial survival function", k: c_int, n: c_int, p: f64);
xsf_impl!(nbdtri, xsf_nbdtri, "Negative binomial quantile function", k: c_int, n: c_int, p: f64);
xsf_impl!(owens_t, xsf_owens_t, "Owen's T function", h: f64, a: f64);
xsf_impl!(pdtr, xsf_pdtr, "Poisson distribution function", k: f64, m: f64);
xsf_impl!(pdtrc, xsf_pdtrc, "Poisson survival function", k: f64, m: f64);
xsf_impl!(pdtri, xsf_pdtri, "Poisson quantile function", k: c_int, y: f64);
xsf_impl!(smirnov, xsf_smirnov, "Kolmogorov-Smirnov survival function", n: c_int, x: f64);
xsf_impl!(smirnovc, xsf_smirnovc, "Kolmogorov-Smirnov distribution function", n: c_int, x: f64);
xsf_impl!(smirnovi, xsf_smirnovi, "Inverse of `smirnov`", n: c_int, x: f64);
xsf_impl!(smirnovp, xsf_smirnovp, "Derivative of `smirnov`", n: c_int, x: f64);
xsf_impl!(
    tukeylambdacdf,
    xsf_tukeylambdacdf,
    "Tukey-Lambda distribution function",
    x: f64,
    lmbda: f64
);

// struve.h
xsf_impl!(itstruve0, xsf_itstruve0, "Integral of the Struve function of order 0", x: f64);
xsf_impl!(it2struve0, xsf_it2struve0, "Integral related to the Struve function of order 0", x: f64);
xsf_impl!(
    itmodstruve0,
    xsf_itmodstruve0,
    "Integral of the modified Struve function of order 0",
    x: f64
);
xsf_impl!(struve_h, xsf_struve_h, "Struve `H` function", v: f64, x: f64);
xsf_impl!(struve_l, xsf_struve_l, "Struve `L` function", v: f64, x: f64);

// trig.h
xsf_impl!(sinpi, xsf_sinpi, "Compute `sin(pi * x)`", x: f64);
xsf_impl!(cospi, xsf_cospi, "Compute `cos(pi * x)`", x: f64);
xsf_impl!(sindg, xsf_sindg, "Circular sine of angle in degrees", x: f64);
xsf_impl!(cosdg, xsf_cosdg, "Circular cosine of angle in degrees", x: f64);
xsf_impl!(tandg, xsf_tandg, "Circular tangent of argument in degrees", x: f64);
xsf_impl!(cotdg, xsf_cotdg, "Circular cotangent of argument in degrees", x: f64);
xsf_impl!(radian, xsf_radian, "Degrees, minutes, seconds to radians", d: f64, m: f64, s: f64);
xsf_impl!(cosm1, xsf_cosm1, "Compute `cos(x) - 1`", x: f64);

// wright_bessel.h
xsf_impl!(
    wright_bessel,
    xsf_wright_bessel,
    "Wright's generalized Bessel function for scalar arguments",
    a: f64,
    b: f64,
    x: f64
);
xsf_impl!(
    log_wright_bessel,
    xsf_log_wright_bessel,
    "Logarithm of `wright_bessel`",
    a: f64,
    b: f64,
    x: f64
);

// zeta.h
xsf_impl!(riemann_zeta, xsf_riemann_zeta, "Riemann zeta function", x: f64);
xsf_impl!(zeta, xsf_zeta, "Riemann zeta function of two arguments", x: f64, q: f64);
xsf_impl!(zetac, xsf_zetac, "Riemann zeta function, minus one", x: f64);

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
