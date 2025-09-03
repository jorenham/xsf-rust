use std::os::raw::c_int;

pub(crate) mod bindings;
use bindings::xsf_impl;

// alg.h
mod alg;
pub use alg::cbrt;

// bessel.h
mod bessel;
pub use bessel::{
    besselpoly, cyl_bessel_i, cyl_bessel_i0, cyl_bessel_i0e, cyl_bessel_i1, cyl_bessel_i1e,
    cyl_bessel_ie, cyl_bessel_j, cyl_bessel_j0, cyl_bessel_j1, cyl_bessel_je, cyl_bessel_k,
    cyl_bessel_k0, cyl_bessel_k0e, cyl_bessel_k1, cyl_bessel_k1e, cyl_bessel_ke, cyl_bessel_y,
    cyl_bessel_y0, cyl_bessel_y1, cyl_bessel_ye, cyl_hankel_1, cyl_hankel_1e, cyl_hankel_2,
    cyl_hankel_2e,
};

// beta.h
mod beta;
pub use beta::{beta, betaln};

// binom.h
mod binom;
pub use binom::binom;

// cdflib.h
mod cdflib;
pub use cdflib::gdtrib;

// digamma.h
mod digamma;
pub use digamma::digamma;

// ellip.h
mod ellip;
pub use ellip::{ellipe, ellipeinc, ellipk, ellipkinc, ellipkm1};

// erf.h
mod erf;
pub use erf::{dawsn, erf, erfc, erfcx, erfi, voigt_profile, wofz};

// exp.h
mod exp;
pub use exp::{exp2, exp10, expm1};

// expint.h
mod expint;
pub use expint::{exp1, expi, scaled_exp1};

// gamma.h
mod gamma;
pub use gamma::{gamma, gammainc, gammaincc, gammainccinv, gammaincinv, gammaln, gammasgn};

// hyp2f1.h
mod hyp2f1;
pub use hyp2f1::hyp2f1;

// iv_ratio.h
mod iv_ratio;
pub use iv_ratio::{iv_ratio, iv_ratio_c};

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
