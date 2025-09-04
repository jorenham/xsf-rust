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

// sph_bessel.h
mod sph_bessel;
pub use sph_bessel::{
    sph_bessel_i, sph_bessel_i_jac, sph_bessel_j, sph_bessel_j_jac, sph_bessel_k, sph_bessel_k_jac,
    sph_bessel_y, sph_bessel_y_jac,
};

// sph_harm.h
mod sph_harm;
pub use sph_harm::sph_harm_y;

// sphd_wave.h
mod sphd_wave;
pub use sphd_wave::{oblate_segv, prolate_segv};

// stats.h
mod stats;
pub use stats::{
    bdtr, bdtrc, bdtri, chdtr, chdtrc, chdtri, fdtr, fdtrc, fdtri, gdtr, gdtrc, kolmogc, kolmogci,
    kolmogi, kolmogorov, kolmogp, log_ndtr, nbdtr, nbdtrc, nbdtri, ndtr, ndtri, owens_t, pdtr,
    pdtrc, pdtri, smirnov, smirnovc, smirnovci, smirnovi, smirnovp, tukeylambdacdf,
};

// struve.h
mod struve;
pub use struve::{it2struve0, itmodstruve0, itstruve0, struve_h, struve_l};

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
