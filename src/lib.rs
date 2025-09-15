#![cfg_attr(not(test), no_std)]

extern crate alloc;

mod bindings;
mod utils;

#[cfg(test)]
mod testing;

mod scipy_special;
pub use scipy_special::boxcox::{boxcox, boxcox1p, inv_boxcox, inv_boxcox1p};

// airy.h
mod airy;
pub use airy::{airy, airyb, airye, airyzo, itairy};

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
    cyl_hankel_2e, it1i0k0, it1j0y0, it2i0k0, it2j0y0, rctj, rcty,
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
pub use ellip::{ellipe, ellipeinc, ellipj, ellipk, ellipkinc, ellipkm1};

// erf.h
mod erf;
pub use erf::{dawsn, erf, erfc, erfcx, erfi, voigt_profile, wofz};

// evalpoly.h
mod evalpoly;
pub use evalpoly::cevalpoly;

// exp.h
mod exp;
pub use exp::{exp2, exp10, expm1};

// expint.h
mod expint;
pub use expint::{exp1, expi, scaled_exp1};

// fp_error_metrics.h
mod fp_error_metrics;
pub use fp_error_metrics::{extended_absolute_error, extended_relative_error};

// fresnel.h
mod fresnel;
pub use fresnel::{fresnel, fresnel_zeros, modified_fresnel_minus, modified_fresnel_plus};

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
pub use kelvin::{
    bei, bei_zeros, beip, beip_zeros, ber, ber_zeros, berp, berp_zeros, kei, kei_zeros, keip,
    keip_zeros, kelvin, kelvin_zeros, ker, ker_zeros, kerp, kerp_zeros,
};

// lambertw.h
mod lambertw;
pub use lambertw::lambertw;

// legendre.h
mod legendre;
pub use legendre::{
    assoc_legendre_p, assoc_legendre_p_all, assoc_legendre_p_norm, assoc_legendre_p_norm_all,
    assoc_legendre_q_all, legendre_p, legendre_p_all, legendre_q_all, sph_legendre_p,
    sph_legendre_p_all,
};

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
pub use mathieu::{cem, cem_cva, mcm1, mcm2, msm1, msm2, sem, sem_cva};

// par_cyl.h
mod par_cyl;
pub use par_cyl::{pbdv, pbvv, pbwa};

// sici.h
mod sici;
pub use sici::{shichi, sici};

// specfun.h
mod specfun;
pub use specfun::{hyp1f1, hypu, pmv};

// sph_bessel.h
mod sph_bessel;
pub use sph_bessel::{
    sph_bessel_i, sph_bessel_i_jac, sph_bessel_j, sph_bessel_j_jac, sph_bessel_k, sph_bessel_k_jac,
    sph_bessel_y, sph_bessel_y_jac,
};

// sph_harm.h
mod sph_harm;
pub use sph_harm::{sph_harm_y, sph_harm_y_all};

// sphd_wave.h
mod sphd_wave;
pub use sphd_wave::{
    oblate_aswfa, oblate_aswfa_nocv, oblate_radial1, oblate_radial1_nocv, oblate_radial2,
    oblate_radial2_nocv, oblate_segv, prolate_aswfa, prolate_aswfa_nocv, prolate_radial1,
    prolate_radial1_nocv, prolate_radial2, prolate_radial2_nocv, prolate_segv,
};

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
mod trig;
pub use trig::{cosdg, cosm1, cospi, cotdg, radian, sindg, sinpi, tandg};

// wright_bessel.h
mod wright_bessel;
pub use wright_bessel::{log_wright_bessel, wright_bessel};

// zeta.h
mod zeta;
pub use zeta::{riemann_zeta, zeta, zetac};
