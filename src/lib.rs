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
mod kelvin;
pub use kelvin::{bei, beip, ber, berp, kei, keip, ker, kerp};

// lambertw.h
mod lambertw;
pub use lambertw::lambertw;

// legendre.h
mod legendre;
pub use legendre::{assoc_legendre_p, assoc_legendre_p_norm, legendre_p, sph_legendre_p};

// log_exp.h
mod log_exp;
pub use log_exp::{expit, exprel, log_expit, log1mexp, logit};

// log.h
mod log;
pub use log::{log1p, log1pmx, xlog1py, xlogy};

// loggamma.h
mod loggamma;
pub use loggamma::{loggamma, rgamma};

// mathieu.h
mod mathieu;
pub use mathieu::{cem_cva, sem_cva};

// specfun.h
mod specfun;
pub use specfun::{hyp1f1, hyperu, hypu, pmv};

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
