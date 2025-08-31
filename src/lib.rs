use std::ffi::c_int;

mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

macro_rules! xsf_impl {
    ($xsf_name:ident, $name:ident, ($($param:ident: $type:ty),*)) => {
        xsf_impl!($xsf_name, $name, ($($param: $type),*), "");
    };
    ($xsf_name:ident, $name:ident, ($($param:ident: $type:ty),*), $docs:expr) => {
        #[allow(clippy::empty_docs)]
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { ffi::$xsf_name($($param),*) }
        }
    };
}

// alg.h
xsf_impl!(xsf_cbrt, cbrt, (x: f64), "Cube root");

// bessel.h
xsf_impl!(xsf_cyl_bessel_j, cyl_bessel_j, (v: f64, x: f64), "Bessel function, 1st kind");
xsf_impl!(
    xsf_cyl_bessel_je,
    cyl_bessel_je,
    (v: f64, x: f64),
    "Exponentially scaled Bessel function, 1st kind"
);
xsf_impl!(xsf_cyl_bessel_y, cyl_bessel_y, (v: f64, x: f64), "Bessel function, 2nd kind");
xsf_impl!(
    xsf_cyl_bessel_ye,
    cyl_bessel_ye,
    (v: f64, x: f64),
    "Exponentially scaled Bessel function, 2nd kind"
);
xsf_impl!(xsf_cyl_bessel_i, cyl_bessel_i, (v: f64, x: f64), "Modified Bessel function, 1st kind");
xsf_impl!(
    xsf_cyl_bessel_ie,
    cyl_bessel_ie,
    (v: f64, x: f64),
    "Exponentially scaled modified Bessel function, 1st kind"
);
xsf_impl!(xsf_cyl_bessel_k, cyl_bessel_k, (v: f64, x: f64), "Modified Bessel function, 2nd kind");
xsf_impl!(
    xsf_cyl_bessel_ke,
    cyl_bessel_ke,
    (v: f64, x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind"
);

xsf_impl!(xsf_cyl_bessel_j0, cyl_bessel_j0, (x: f64), "Bessel function, 1st kind, order 0");
xsf_impl!(xsf_cyl_bessel_j1, cyl_bessel_j1, (x: f64), "Bessel function, 1st kind, order 1");
xsf_impl!(xsf_cyl_bessel_y0, cyl_bessel_y0, (x: f64), "Bessel function, 2nd kind, order 0");
xsf_impl!(xsf_cyl_bessel_y1, cyl_bessel_y1, (x: f64), "Bessel function, 2nd kind, order 1");
xsf_impl!(
    xsf_cyl_bessel_i0,
    cyl_bessel_i0,
    (x: f64),
    "Modified Bessel function, 1st kind, order 0"
);
xsf_impl!(
    xsf_cyl_bessel_i0e,
    cyl_bessel_i0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 0"
);
xsf_impl!(
    xsf_cyl_bessel_i1,
    cyl_bessel_i1,
    (x: f64),
    "Modified Bessel function, 1st kind, order 1"
);
xsf_impl!(
    xsf_cyl_bessel_i1e,
    cyl_bessel_i1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 1st kind, order 1"
);
xsf_impl!(
    xsf_cyl_bessel_k0,
    cyl_bessel_k0,
    (x: f64),
    "Modified Bessel function, 2nd kind, order 0"
);
xsf_impl!(
    xsf_cyl_bessel_k0e,
    cyl_bessel_k0e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 0"
);
xsf_impl!(
    xsf_cyl_bessel_k1,
    cyl_bessel_k1,
    (x: f64),
    "Modified Bessel function, 2nd kind, order 1"
);
xsf_impl!(
    xsf_cyl_bessel_k1e,
    cyl_bessel_k1e,
    (x: f64),
    "Exponentially scaled modified Bessel function, 2nd kind, order 1"
);

xsf_impl!(
    xsf_besselpoly,
    besselpoly,
    (a: f64, lambda: f64, nu: f64),
    "Weighted integral of the Bessel function of the first kind"
);

// beta.h
xsf_impl!(xsf_beta, beta, (a: f64, b: f64), "Beta function");
xsf_impl!(xsf_betaln, betaln, (a: f64, b: f64), "Natural log of `|beta|`");

// binom.h
xsf_impl!(xsf_binom, binom, (n: f64, k: f64), "Binomial coefficient");

// digamma.h
xsf_impl!(xsf_digamma, digamma, (x: f64), "Digamma function");

// erf.h
xsf_impl!(xsf_erf, erf, (x: f64), "Error function");
xsf_impl!(xsf_erfc, erfc, (x: f64), "Complementary error function `1 - erf(x)`");
xsf_impl!(xsf_erfcx, erfcx, (x: f64), "Scaled complementary error function `exp(x^2) * erfc(x)`");
xsf_impl!(xsf_erfi, erfi, (x: f64), "Imaginary error function `-i erf(ix)`");
xsf_impl!(xsf_voigt_profile, voigt_profile, (x: f64, sigma: f64, gamma: f64), "Voigt profile");
xsf_impl!(xsf_dawsn, dawsn, (x: f64), "Dawson function `sqrt(pi)/2 * exp(-x^2) * erfi(x)`");

// exp.h
xsf_impl!(xsf_expm1, expm1, (x: f64), "`exp(x) - 1`");
xsf_impl!(xsf_exp2, exp2, (x: f64), "`2^x`");
xsf_impl!(xsf_exp10, exp10, (x: f64), "`10^x`");

// expint.h
xsf_impl!(xsf_exp1, exp1, (x: f64), "Exponential integral `E_1(x)`");
xsf_impl!(xsf_expi, expi, (x: f64), "Exponential integral `E_i(x)`");
xsf_impl!(
    xsf_scaled_exp1,
    scaled_exp1,
    (x: f64),
    "Scaled version of the exponential integral `E_1(x)`"
);

// fp_error_metrics.h
// xsf_impl!(extended_absolute_error, xsf_extended_absolute_error, actual: f64, desired: f64);
// xsf_impl!(extended_relative_error, xsf_extended_relative_error, actual: f64, desired: f64);

// gamma.h
xsf_impl!(xsf_gamma, gamma, (x: f64), "Gamma function");
xsf_impl!(xsf_gamma_ratio, gamma_ratio, (a: f64, b: f64), "`gamma(a) / gamma(b)`");
xsf_impl!(xsf_gammainc, gammainc, (a: f64, x: f64), "Incomplete Gamma integral");
xsf_impl!(xsf_gammaincc, gammaincc, (a: f64, x: f64), "Complemented incomplete Gamma integral");
xsf_impl!(xsf_gammaincinv, gammaincinv, (a: f64, p: f64), "Inverse of `gammainc`");
xsf_impl!(xsf_gammainccinv, gammainccinv, (a: f64, p: f64), "Inverse of `gammaincc`");
xsf_impl!(xsf_gammaln, gammaln, (x: f64), "Natural logarithm of Gamma function");
xsf_impl!(xsf_gammasgn, gammasgn, (x: f64), "Sign of the Gamma function");

// hyp2f1.h
xsf_impl!(
    xsf_hyp2f1,
    hyp2f1,
    (a: f64, b: f64, c: f64, x: f64),
    "Gauss hypergeometric function `2F1`"
);

// iv_ratio.h
xsf_impl!(
    xsf_iv_ratio,
    iv_ratio,
    (v: f64, x: f64),
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function, 1st kind"
);
xsf_impl!(
    xsf_iv_ratio_c,
    iv_ratio_c,
    (v: f64, x: f64),
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function, 1st kind"
);

// kelvin.h
xsf_impl!(xsf_ber, ber, (x: f64), "Kelvin function `ber`");
xsf_impl!(xsf_bei, bei, (x: f64), "Kelvin function `bei`");
xsf_impl!(xsf_ker, ker, (x: f64), "Kelvin function `ker`");
xsf_impl!(xsf_kei, kei, (x: f64), "Kelvin function `kei`");
xsf_impl!(xsf_berp, berp, (x: f64), "Derivative of the Kelvin function `ber`");
xsf_impl!(xsf_beip, beip, (x: f64), "Derivative of the Kelvin function `bei`");
xsf_impl!(xsf_kerp, kerp, (x: f64), "Derivative of the Kelvin function `ker`");
xsf_impl!(xsf_keip, keip, (x: f64), "Derivative of the Kelvin function `kei`");

// legendre.h
xsf_impl!(xsf_legendre_p, legendre_p, (n: c_int, x: f64), "Legendre polynomial of degree n");
xsf_impl!(
    xsf_sph_legendre_p,
    sph_legendre_p,
    (n: c_int, m: c_int, theta: f64),
    "Spherical Legendre polynomial of degree n and order m"
);

// log_exp.h
xsf_impl!(xsf_expit, expit, (x: f64), "Expit function, `1/(1 + exp(-x))`");
xsf_impl!(xsf_exprel, exprel, (x: f64), "Relative error exponential, `(exp(x) - 1)/x`");
xsf_impl!(xsf_logit, logit, (x: f64), "Logit function, `log(x / (1 - x))`");
xsf_impl!(xsf_log_expit, log_expit, (x: f64), "Log of the expit function, `log(expit(x))`");
xsf_impl!(xsf_log1mexp, log1mexp, (x: f64), "Compute `log(1 - exp(x))`");

// log.h
xsf_impl!(xsf_log1p, log1p, (x: f64), "Compute `log(1 + x)`");
xsf_impl!(xsf_log1pmx, log1pmx, (x: f64), "Compute `log(1 + x) - x`");
xsf_impl!(xsf_xlogy, xlogy, (x: f64, y: f64), "Compute `x * log(y)`");
xsf_impl!(xsf_xlog1py, xlog1py, (x: f64, y: f64), "Compute `x * log(1 + y)`");

// loggamma.h
xsf_impl!(xsf_loggamma, loggamma, (x: f64), "Principal branch of the logarithm of `gamma(x)`");
xsf_impl!(xsf_rgamma, rgamma, (x: f64), "Reciprocal Gamma function");

// mathieu.h
xsf_impl!(xsf_cem_cva, cem_cva, (m: f64, q: f64), "Characteristic value of even Mathieu functions");
xsf_impl!(xsf_sem_cva, sem_cva, (m: f64, q: f64), "Characteristic value of odd Mathieu functions");

// specfun.h
xsf_impl!(xsf_hypu, hypu, (a: f64, b: f64, x: f64), "Confluent hypergeometric function `U`");
xsf_impl!(xsf_hyp1f1, hyp1f1, (a: f64, b: f64, x: f64), "Confluent hypergeometric function `1F1`");
xsf_impl!(
    xsf_pmv,
    pmv,
    (m: f64, v: f64, x: f64),
    "Associated Legendre function of integer order and real degree"
);

// sphd_wave.h
xsf_impl!(
    xsf_prolate_segv,
    prolate_segv,
    (m: f64, n: f64, c: f64),
    "Characteristic value of prolate spheroidal function"
);
xsf_impl!(
    xsf_oblate_segv,
    oblate_segv,
    (m: f64, n: f64, c: f64),
    "Characteristic value of oblate spheroidal function"
);

// stats.h
xsf_impl!(xsf_bdtr, bdtr, (k: f64, n: c_int, p: f64), "Binomial distribution function");
xsf_impl!(xsf_bdtrc, bdtrc, (k: f64, n: c_int, p: f64), "Binomial survival function");
xsf_impl!(xsf_bdtri, bdtri, (k: f64, n: c_int, y: f64), "Binomial quantile function");
xsf_impl!(xsf_chdtr, chdtr, (df: f64, x: f64), "Chi-squared distribution function");
xsf_impl!(xsf_chdtrc, chdtrc, (df: f64, x: f64), "Chi-squared survival function");
xsf_impl!(xsf_chdtri, chdtri, (df: f64, y: f64), "Chi-squared quantile function");
xsf_impl!(xsf_fdtr, fdtr, (a: f64, b: f64, x: f64), "F distribution function");
xsf_impl!(xsf_fdtrc, fdtrc, (a: f64, b: f64, x: f64), "F survival function");
xsf_impl!(xsf_fdtri, fdtri, (a: f64, b: f64, y: f64), "F quantile function");
xsf_impl!(xsf_gdtr, gdtr, (a: f64, b: f64, x: f64), "Gamma distribution function");
xsf_impl!(xsf_gdtrc, gdtrc, (a: f64, b: f64, x: f64), "Gamma survival function");
xsf_impl!(xsf_kolmogorov, kolmogorov, (x: f64), "Kolmogorov survival function");
xsf_impl!(xsf_kolmogc, kolmogc, (x: f64), "Kolmogorov distribution function");
xsf_impl!(xsf_kolmogi, kolmogi, (x: f64), "Inverse of `kolmogorov`");
xsf_impl!(xsf_kolmogp, kolmogp, (x: f64), "Derivative of `kolmogorov`");
xsf_impl!(xsf_ndtr, ndtr, (x: f64), "Normal distribution function");
xsf_impl!(xsf_ndtri, ndtri, (x: f64), "Normal quantile function");
xsf_impl!(xsf_log_ndtr, log_ndtr, (x: f64), "Log of the normal distribution function");
xsf_impl!(
    xsf_nbdtr,
    nbdtr,
    (k: c_int, n: c_int, p: f64),
    "Negative binomial distribution function"
);
xsf_impl!(xsf_nbdtrc, nbdtrc, (k: c_int, n: c_int, p: f64), "Negative binomial survival function");
xsf_impl!(xsf_nbdtri, nbdtri, (k: c_int, n: c_int, p: f64), "Negative binomial quantile function");
xsf_impl!(xsf_owens_t, owens_t, (h: f64, a: f64), "Owen's T function");
xsf_impl!(xsf_pdtr, pdtr, (k: f64, m: f64), "Poisson distribution function");
xsf_impl!(xsf_pdtrc, pdtrc, (k: f64, m: f64), "Poisson survival function");
xsf_impl!(xsf_pdtri, pdtri, (k: c_int, y: f64), "Poisson quantile function");
xsf_impl!(xsf_smirnov, smirnov, (n: c_int, x: f64), "Kolmogorov-Smirnov survival function");
xsf_impl!(xsf_smirnovc, smirnovc, (n: c_int, x: f64), "Kolmogorov-Smirnov distribution function");
xsf_impl!(xsf_smirnovi, smirnovi, (n: c_int, x: f64), "Inverse of `smirnov`");
xsf_impl!(xsf_smirnovp, smirnovp, (n: c_int, x: f64), "Derivative of `smirnov`");
xsf_impl!(
    xsf_tukeylambdacdf,
    tukeylambdacdf,
    (x: f64, lmbda: f64),
    "Tukey-Lambda distribution function"
);

// struve.h
xsf_impl!(xsf_itstruve0, itstruve0, (x: f64), "Integral of the Struve function of order 0");
xsf_impl!(
    xsf_it2struve0,
    it2struve0,
    (x: f64),
    "Integral related to the Struve function of order 0"
);
xsf_impl!(
    xsf_itmodstruve0,
    itmodstruve0,
    (x: f64),
    "Integral of the modified Struve function of order 0"
);
xsf_impl!(xsf_struve_h, struve_h, (v: f64, x: f64), "Struve `H` function");
xsf_impl!(xsf_struve_l, struve_l, (v: f64, x: f64), "Struve `L` function");

// trig.h
xsf_impl!(xsf_sinpi, sinpi, (x: f64), "Compute `sin(pi * x)`");
xsf_impl!(xsf_cospi, cospi, (x: f64), "Compute `cos(pi * x)`");
xsf_impl!(xsf_sindg, sindg, (x: f64), "Circular sine of angle in degrees");
xsf_impl!(xsf_cosdg, cosdg, (x: f64), "Circular cosine of angle in degrees");
xsf_impl!(xsf_tandg, tandg, (x: f64), "Circular tangent of argument in degrees");
xsf_impl!(xsf_cotdg, cotdg, (x: f64), "Circular cotangent of argument in degrees");
xsf_impl!(xsf_radian, radian, (d: f64, m: f64, s: f64), "Degrees, minutes, seconds to radians");
xsf_impl!(xsf_cosm1, cosm1, (x: f64), "Compute `cos(x) - 1`");

// wright_bessel.h
xsf_impl!(
    xsf_wright_bessel,
    wright_bessel,
    (a: f64, b: f64, x: f64),
    "Wright's generalized Bessel function for scalar arguments"
);
xsf_impl!(
    xsf_log_wright_bessel,
    log_wright_bessel,
    (a: f64, b: f64, x: f64),
    "Logarithm of `wright_bessel`"
);

// zeta.h
xsf_impl!(xsf_riemann_zeta, riemann_zeta, (x: f64), "Riemann zeta function");
xsf_impl!(xsf_zeta, zeta, (x: f64, q: f64), "Riemann zeta function of two arguments");
xsf_impl!(xsf_zetac, zetac, (x: f64), "Riemann zeta function, minus one");

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
