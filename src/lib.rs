#![allow(non_snake_case)]

use num_complex::Complex;
use std::os::raw::c_int;

mod xsf {
    use num_complex::Complex;

    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

    impl<T> std_complex<T> {
        fn new(re: T, im: T) -> Self {
            Self {
                _phantom_0: std::marker::PhantomData,
                _M_real: re,
                _M_imag: im,
            }
        }
    }

    impl<T> From<std_complex<T>> for Complex<T> {
        fn from(val: std_complex<T>) -> Self {
            Complex::new(val._M_real, val._M_imag)
        }
    }

    impl<T> From<Complex<T>> for std_complex<T> {
        fn from(c: Complex<T>) -> Self {
            Self::new(c.re, c.im)
        }
    }
}

macro_rules! xsf_impl {
    ($name:ident, ($($param:ident: $type:ty),*), $docs:expr) => {
        #[doc = $docs]
        pub fn $name($($param: $type),*) -> f64 {
            unsafe { xsf::$name($($param),*) }
        }
    };
}

// alg.h
/// Cube root
pub fn cbrt(x: f64) -> f64 {
    unsafe { xsf::cbrt_(x) }
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
/// Hankel function, 1st kind
pub fn cyl_hankel_1(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_1(v, z.into()) }.into()
}
/// Exponentially scaled Hankel function, 1st kind
pub fn cyl_hankel_1e(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_1e(v, z.into()) }.into()
}
/// Hankel function, 2nd kind
pub fn cyl_hankel_2(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_2(v, z.into()) }.into()
}
/// Exponentially scaled Hankel function, 2nd kind
pub fn cyl_hankel_2e(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_2e(v, z.into()) }.into()
}
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
    unsafe { xsf::erf_(x) }
}
/// Complementary error function `1 - erf(x)`
pub fn erfc(x: f64) -> f64 {
    unsafe { xsf::erfc_(x) }
}
xsf_impl!(erfcx, (x: f64), "Scaled complementary error function `exp(x^2) * erfc(x)`");
xsf_impl!(erfi, (x: f64), "Imaginary error function `-i erf(ix)`");
xsf_impl!(voigt_profile, (x: f64, sigma: f64, gamma: f64), "Voigt profile");
xsf_impl!(dawsn, (x: f64), "Dawson function `sqrt(pi)/2 * exp(-x^2) * erfi(x)`");

// exp.h
/// `exp(x) - 1`
pub fn expm1(x: f64) -> f64 {
    unsafe { xsf::expm1_(x) }
}
/// `2^x`
pub fn exp2(x: f64) -> f64 {
    unsafe { xsf::exp2_(x) }
}
/// `10^x`
pub fn exp10(x: f64) -> f64 {
    unsafe { xsf::exp10_(x) }
}

// expint.h
xsf_impl!(exp1, (x: f64), "Exponential integral `E_1(x)`");
xsf_impl!(expi, (x: f64), "Exponential integral `E_i(x)`");
xsf_impl!(scaled_exp1, (x: f64), "Scaled version of the exponential integral `E_1(x)`");

// gamma.h
/// Gamma function
pub fn gamma(x: f64) -> f64 {
    unsafe { xsf::gamma_(x) }
}
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
