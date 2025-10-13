//! Bindings to the scipy XSF C++ library of `scipy.special` functions.
//!
//! See the [`scipy.special` documentation](https://docs.scipy.org/doc/scipy/reference/special.html)
//! for additional information.
//!
//! # Airy functions
//!
//! | Function      | Description                                                 |
//! | ------------- | ----------------------------------------------------------- |
//! | [`airy`]      | Airy functions and their derivatives                        |
//! | [`airye`]     | Exponentially scaled Airy functions and their derivatives   |
//! | [`ai_zeros`]  | Zeros and values of the Airy function Ai and its derivative |
//! | [`bi_zeros`]  | Zeros and values of the Airy function Bi and its derivative |
//! | [`itairy`]    | Integrals of Airy functions                                 |
//!
//! # Elliptic functions and integrals
//!
//! | Function         | Description                                                 |
//! | ---------------- | ----------------------------------------------------------- |
//! | [`ellipj`]       | Jacobian elliptic functions                                 |
//! | [`ellipk`]       | Complete elliptic integral of the first kind                |
//! | [`ellipkm1`]     | Complete elliptic integral of the first kind around `m = 1` |
//! | [`ellipkinc`]    | Incomplete elliptic integral of the first kind              |
//! | [`ellipe`]       | Complete elliptic integral of the second kind               |
//! | [`ellipeinc`]    | Incomplete elliptic integral of the second kind             |
//!
//! # Bessel functions
//!
//! | Function                | Description                                                       |
//! | ----------------------- | ----------------------------------------------------------------- |
//! | [`bessel_j`]            | Bessel function of the first kind                                 |
//! | [`bessel_je`]           | Exponentially scaled Bessel function of the first kind            |
//! | [`bessel_y`]            | Bessel function of the second kind                                |
//! | [`bessel_ye`]           | Exponentially scaled Bessel function of the second kind           |
//! | [`bessel_i`]            | Modified Bessel function of the first kind                        |
//! | [`bessel_ie`]           | Exponentially scaled modified Bessel function of the first kind   |
//! | [`bessel_k`]            | Modified Bessel function of the second kind                       |
//! | [`bessel_ke`]           | Exponentially scaled modified Bessel function of the second kind  |
//! | [`hankel_1`]            | Hankel function of the first kind                                 |
//! | [`hankel_1e`]           | Exponentially scaled Hankel function of the first kind            |
//! | [`hankel_2`]            | Hankel function of the second kind                                |
//! | [`hankel_2e`]           | Exponentially scaled Hankel function of the second kind           |
//! | [`wright_bessel`]       | Wright's generalized Bessel function                              |
//! | [`log_wright_bessel`]   | Natural logarithm of Wright's generalized Bessel function         |
//! | [`jahnke_emden_lambda`] | Jahnke-Emden Lambda function Λ<sub>v</sub>(x) and its derivatives |
//!
//! ## Zeros of Bessel functions
//!
//! | Function         | Description                                                             |
//! | ---------------- | ----------------------------------------------------------------------- |
//! | [`bessel_zeros`] | Compute `NT` zeros of Bessel functions Jn(x), Jn'(x), Yn(x), and Yn'(x) |
//!
//! ## Faster versions of common Bessel functions
//!
//! | Function       | Description                                                                 |
//! | -------------- | --------------------------------------------------------------------------- |
//! | [`bessel_j0`]  | Bessel function of the first kind of order 0                                |
//! | [`bessel_j1`]  | Bessel function of the first kind of order 1                                |
//! | [`bessel_y0`]  | Bessel function of the second kind of order 0                               |
//! | [`bessel_y1`]  | Bessel function of the second kind of order 1                               |
//! | [`bessel_i0`]  | Modified Bessel function of the first kind of order 0                       |
//! | [`bessel_i0e`] | Exponentially scaled modified Bessel function of the first kind of order 0  |
//! | [`bessel_i1`]  | Modified Bessel function of the first kind of order 1                       |
//! | [`bessel_i1e`] | Exponentially scaled modified Bessel function of the first kind of order 1  |
//! | [`bessel_k0`]  | Modified Bessel function of the second kind of order 0                      |
//! | [`bessel_k0e`] | Exponentially scaled modified Bessel function of the second kind of order 0 |
//! | [`bessel_k1`]  | Modified Bessel function of the second kind of order 1                      |
//! | [`bessel_k1e`] | Exponentially scaled modified Bessel function of the second kind of order 1 |
//!
//! ## Integrals of Bessel functions
//!
//! | Function       | Description                                                                 |
//! | -------------- | --------------------------------------------------------------------------- |
//! | [`it1j0y0`]    | Integral of Bessel functions of the first kind of order 0                   |
//! | [`it2j0y0`]    | Integral related to Bessel functions of the first kind of order 0           |
//! | [`it1i0k0`]    | Integral of modified Bessel functions of the second kind of order 0         |
//! | [`it2i0k0`]    | Integral related to modified Bessel functions of the second kind of order 0 |
//! | [`besselpoly`] | Weighted integral of the Bessel function of the first kind                  |
//!
//! ## Derivatives of Bessel functions
//!
//! | Function           | Description                                 |
//! | ------------------ | ------------------------------------------- |
//! | [`bessel_j_prime`] | *n*<sup>th</sup> derivative of [`bessel_j`] |
//! | [`bessel_y_prime`] | *n*<sup>th</sup> derivative of [`bessel_y`] |
//! | [`bessel_i_prime`] | *n*<sup>th</sup> derivative of [`bessel_i`] |
//! | [`bessel_k_prime`] | *n*<sup>th</sup> derivative of [`bessel_k`] |
//! | [`hankel_1_prime`] | *n*<sup>th</sup> derivative of [`hankel_1`] |
//! | [`hankel_2_prime`] | *n*<sup>th</sup> derivative of [`hankel_2`] |
//!
//! ## Spherical Bessel functions
//!
//! | Function               | Description                                           |
//! | ---------------------- | ----------------------------------------------------- |
//! | [`sph_bessel_j`]       | Spherical Bessel function of the first kind           |
//! | [`sph_bessel_j_prime`] | Derivative of [`sph_bessel_j`]                        |
//! | [`sph_bessel_y`]       | Spherical Bessel function of the second kind          |
//! | [`sph_bessel_y_prime`] | Derivative of [`sph_bessel_y`]                        |
//! | [`sph_bessel_i`]       | Modified Spherical Bessel function of the first kind  |
//! | [`sph_bessel_i_prime`] | Derivative of [`sph_bessel_i`]                        |
//! | [`sph_bessel_k`]       | Modified Spherical Bessel function of the second kind |
//! | [`sph_bessel_k_prime`] | Derivative of [`sph_bessel_k`]                        |
//!
//! ## Riccati-Bessel functions
//!
//! | Function      | Description                                                   |
//! | ------------- | ------------------------------------------------------------- |
//! | [`riccati_j`] | Riccati-Bessel function of the first kind and its derivative  |
//! | [`riccati_y`] | Riccati-Bessel function of the second kind and its derivative |
//!
//! # Struve functions
//!
//! | Function         | Description                                                             |
//! | ---------------- | ----------------------------------------------------------------------- |
//! | [`struve_h`]     | Struve function *H<sub>v</sub>(x)*                                      |
//! | [`struve_l`]     | Modified Struve function *L<sub>v</sub>(x)*                             |
//! | [`itstruve0`]    | Integral of the Struve function of order 0, *H<sub>0</sub>(x)*          |
//! | [`it2struve0`]   | Integral related to the Struve function of order 0                      |
//! | [`itmodstruve0`] | Integral of the modified Struve function of order 0, *L<sub>0</sub>(x)* |
//!
//! # Raw statistical functions
//!
//! ## Binomial distribution
//!
//! | Function  | Description                      |
//! | --------- | -------------------------------- |
//! | [`bdtr`]  | Cumulative distribution function |
//! | [`bdtrc`] | Survival function                |
//! | [`bdtri`] | Quantile function                |
//!
//! ## F distribution
//!
//! | Function  | Description                      |
//! | --------- | -------------------------------- |
//! | [`fdtr`]  | Cumulative distribution function |
//! | [`fdtrc`] | Survival function                |
//! | [`fdtri`] | Quantile function                |
//!
//! ## Gamma distribution
//!
//! | Function   | Description                                            |
//! | ---------- | ------------------------------------------------------ |
//! | [`gdtr`]   | Cumulative distribution function                       |
//! | [`gdtrc`]  | Survival function                                      |
//! | [`gdtrib`] | Inverse of [`gdtr(a, b, x)`](gdtr) with respect to `b` |
//!
//! ## Negative binomial distribution
//!
//! | Function   | Description                      |
//! | ---------- | -------------------------------- |
//! | [`nbdtr`]  | Cumulative distribution function |
//! | [`nbdtrc`] | Survival function                |
//! | [`nbdtri`] | Quantile function                |
//!
//! ## Normal distribution
//!
//! | Function     | Description                      |
//! | ------------ | -------------------------------- |
//! | [`ndtr`]     | Cumulative distribution function |
//! | [`log_ndtr`] | Logarithm of [`ndtr`]            |
//! | [`ndtri`]    | Quantile function                |
//!
//! ## Poisson distribution
//!
//! | Function  | Description                      |
//! | --------- | -------------------------------- |
//! | [`pdtr`]  | Cumulative distribution function |
//! | [`pdtrc`] | Survival function                |
//! | [`pdtri`] | Quantile function                |
//!
//! ## Chi square distribution
//!
//! | Function   | Description                      |
//! | ---------- | -------------------------------- |
//! | [`chdtr`]  | Cumulative distribution function |
//! | [`chdtrc`] | Survival function                |
//! | [`chdtri`] | Quantile function                |
//!
//! ## Kolmogorov distribution
//!
//! | Function       | Description                      |
//! | -------------- | -------------------------------- |
//! | [`kolmogorov`] | Survival function                |
//! | [`kolmogp`]    | Derivative of [`kolmogorov`]     |
//! | [`kolmogi`]    | Inverse of [`kolmogorov`]        |
//! | [`kolmogc`]    | Cumulative distribution function |
//! | [`kolmogci`]   | Inverse of [`kolmogc`]           |
//!
//! ## Kolmogorov-Smirnov distribution
//!
//! | Function      | Description                      |
//! | ------------- | -------------------------------- |
//! | [`smirnov`]   | Survival function                |
//! | [`smirnovp`]  | Derivative of [`smirnov`]        |
//! | [`smirnovi`]  | Inverse of [`smirnov`]           |
//! | [`smirnovc`]  | Cumulative distribution function |
//! | [`smirnovci`] | Inverse of [`smirnovc`]          |
//!
//! ## Box-Cox transformation
//!
//! | Function         | Description                       |
//! | ---------------- | --------------------------------- |
//! | [`boxcox`]       | Box-Cox transformation of *x*     |
//! | [`boxcox1p`]     | Box-Cox transformation of *1 + x* |
//! | [`inv_boxcox`]   | Inverse of [`boxcox`]             |
//! | [`inv_boxcox1p`] | Inverse of [`boxcox1p`]           |
//!
//! ## Sigmoidal functions
//!
//! | Function       | Description                              |
//! | -------------- | ---------------------------------------- |
//! | [`logit`]      | Logit function, *log(x / (1 - x))*       |
//! | [`expit`]      | Expit function, *1/(1 + e<sup>-x</sup>)* |
//! | [`log_expit`]  | Logarithm of [`expit`]                   |
//!
//! ## Miscellaneous
//!
//! | Function           | Description                                   |
//! | ------------------ | --------------------------------------------- |
//! | [`tukeylambdacdf`] | Tukey-Lambda cumulative distribution function |
//! | [`owens_t`]        | Owen's T function                             |
//!
//! # Gamma and related functions
//!
//! | Function         | Description                                                |
//! | ---------------- | ---------------------------------------------------------- |
//! | [`gamma`]        | Gamma function                                             |
//! | [`gammaln`]      | Logarithm of the absolute value of the gamma function      |
//! | [`loggamma`]     | Principal branch of the logarithm of the gamma function    |
//! | [`gammasgn`]     | Sign of the gamma function                                 |
//! | [`gammainc`]     | Regularized lower incomplete gamma function                |
//! | [`gammaincinv`]  | Inverse to the regularized lower incomplete gamma function |
//! | [`gammaincc`]    | Regularized upper incomplete gamma function                |
//! | [`gammainccinv`] | Inverse to the regularized upper incomplete gamma function |
//! | [`beta`]         | Beta function                                              |
//! | [`betaln`]       | Logarithm of the absolute value of the beta function       |
//! | [`betainc`]      | Regularized incomplete beta function                       |
//! | [`betaincinv`]   | Inverse to the regularized incomplete beta function        |
//! | [`digamma`]      | The digamma function                                       |
//! | [`rgamma`]       | Reciprocal of the gamma function                           |
//! | [`pow_rising`]   | Rising factorial  (Pochhammer symbol)                      |
//! | [`pow_falling`]  | Falling factorial                                          |
//!
//! # Error function and Fresnel integrals
//!
//! | Function                   | Description                                  |
//! | -------------------------- | -------------------------------------------- |
//! | [`erf`]                    | Error function                               |
//! | [`erfc`]                   | Complementary error function, *1 - erf(z)*   |
//! | [`erfcx`]                  | Scaled complementary error function          |
//! | [`erfi`]                   | Imaginary error function *-i erf(i z)*       |
//! | [`erfinv`]                 | Inverse of [`erf`]                           |
//! | [`erfcinv`]                | Inverse of [`erfc`]                          |
//! | [`erf_zeros`]              | Zeros of [`erf`]                             |
//! | [`wofz`]                   | Faddeeva function                            |
//! | [`dawsn`]                  | Dawson's integral                            |
//! | [`fresnel`]                | Fresnel integrals *S(z)* and *C(z)*          |
//! | [`fresnel_zeros`]          | Zeros of Fresnel integrals *S(z)* and *C(z)* |
//! | [`modified_fresnel_plus`]  | Modified Fresnel positive integrals          |
//! | [`modified_fresnel_minus`] | Modified Fresnel negative integrals          |
//! | [`voigt_profile`]          | Voigt profile                                |
//!
//! # Legendre functions
//!
//! | Function                      | Description                                                  |
//! | ----------------------------- | ------------------------------------------------------------ |
//! | [`legendre_p`]                | Legendre polynomial of the first kind                        |
//! | [`legendre_p_all`]            | All Legendre polynomials of the first kind                   |
//! | [`assoc_legendre_p`]          | Associated Legendre polynomial of the 1st kind               |
//! | [`assoc_legendre_p_all`]      | All associated Legendre polynomials of the 1st kind          |
//! | [`assoc_legendre_p_norm`]     | Normalized associated Legendre polynomial                    |
//! | [`assoc_legendre_p_norm_all`] | All normalized associated Legendre polynomials               |
//! | [`sph_legendre_p`]            | Spherical Legendre polynomial of the first kind              |
//! | [`sph_legendre_p_all`]        | All spherical Legendre polynomials of the first kind         |
//! | [`sph_harm_y`]                | Spherical harmonics                                          |
//! | [`sph_harm_y_all`]            | All spherical harmonics                                      |
//! | [`legendre_q_all`]       | All Legendre functions of the 2nd kind and derivatives            |
//! | [`assoc_legendre_q_all`] | All associated Legendre functions of the 2nd kind and derivatives |
//!
//! # Hypergeometric functions
//!
//! | Function   | Description                                                   |
//! | ---------- | ------------------------------------------------------------- |
//! | [`hyp2f1`] | Gauss hypergeometric function *<sub>2</sub>F<sub>1</sub>*     |
//! | [`hyp1f1`] | Confluent hypergeometric function *<sub>2</sub>F<sub>1</sub>* |
//! | [`hypu`]   | Confluent hypergeometric function *U*                         |
//!
//! # Parabolic cylinder functions
//!
//! | Function | Description                                    |
//! | -------- | ---------------------------------------------- |
//! | [`pbdv`] | Parabolic cylinder function *D* and derivative |
//! | [`pbvv`] | Parabolic cylinder function *V* and derivative |
//! | [`pbwa`] | Parabolic cylinder function *W* and derivative |
//!
//! # Mathieu and related functions
//!
//! | Function    | Description                                                          |
//! | ----------- | -------------------------------------------------------------------- |
//! | [`cem_cva`] | Even Mathieu function characteristic value                           |
//! | [`sem_cva`] | Odd Mathieu function characteristic value                            |
//! | [`cem`]     | Even Mathieu function and its derivative                             |
//! | [`sem`]     | Odd Mathieu function and its derivative                              |
//! | [`mcm1`]    | Even modified Mathieu function of the first kind and its derivative  |
//! | [`msm1`]    | Odd modified Mathieu function of the first kind and its derivative   |
//! | [`mcm2`]    | Even modified Mathieu function of the second kind and its derivative |
//! | [`msm2`]    | Odd modified Mathieu function of the second kind and its derivative  |
//!
//! # Spherical wave functions
//!
//! | Function                 | Description                                           |
//! | ------------------------ | ----------------------------------------------------- |
//! | [`prolate_aswfa_nocv`]   | Prolate spheroidal angular function of the first kind |
//! | [`prolate_radial1_nocv`] | Prolate spheroidal radial function of the first kind  |
//! | [`prolate_radial2_nocv`] | Prolate spheroidal radial function of the second kind |
//! | [`oblate_aswfa_nocv`]    | Oblate spheroidal angular function of the first kind  |
//! | [`oblate_radial1_nocv`]  | Oblate spheroidal radial function of the first kind   |
//! | [`oblate_radial2_nocv`]  | Oblate spheroidal radial function of the second kind  |
//! | [`prolate_segv`]         | Characteristic value of prolate spheroidal function   |
//! | [`oblate_segv`]          | Characteristic value of oblate spheroidal function    |
//!
//! The following functions require pre-computed characteristic value:
//!
//! | Function            | Description                                           |
//! | ------------------- | ----------------------------------------------------- |
//! | [`prolate_aswfa`]   | Prolate spheroidal angular function of the first kind |
//! | [`prolate_radial1`] | Prolate spheroidal radial function of the first kind  |
//! | [`prolate_radial2`] | Prolate spheroidal radial function of the second kind |
//! | [`oblate_aswfa`]    | Oblate spheroidal angular function of the first kind  |
//! | [`oblate_radial1`]  | Oblate spheroidal radial function of the first kind   |
//! | [`oblate_radial2`]  | Oblate spheroidal radial function of the second kind  |
//!
//! # Kelvin functions
//!
//! | Function   | Zeros            | Description                               |
//! | ---------- | ---------------- | ----------------------------------------- |
//! | [`kelvin`] | [`kelvin_zeros`] | Kelvin functions as complex numbers       |
//! | [`ber`]    | [`ber_zeros`]    | Kelvin function *ber*                     |
//! | [`bei`]    | [`bei_zeros`]    | Kelvin function *bei*                     |
//! | [`berp`]   | [`berp_zeros`]   | Derivative of the Kelvin function [`ber`] |
//! | [`beip`]   | [`beip_zeros`]   | Derivative of the Kelvin function [`bei`] |
//! | [`ker`]    | [`ker_zeros`]    | Kelvin function *ker*                     |
//! | [`kei`]    | [`kei_zeros`]    | Kelvin function *kei*                     |
//! | [`kerp`]   | [`kerp_zeros`]   | Derivative of the Kelvin function [`ker`] |
//! | [`keip`]   | [`keip_zeros`]   | Derivative of the Kelvin function [`kei`] |
//!
//! # Exponential integrals
//!
//! | Function        | Description                                                    |
//! | --------------- | -------------------------------------------------------------- |
//! | [`expn`]        | Generalized exponential integral *E<sub>n</sub>(x)*            |
//! | [`expi`]        | Exponential integral *Ei(x)*                                   |
//! | [`exp1`]        | Exponential integral *E<sub>1</sub>(x)*                        |
//! | [`scaled_exp1`] | Scaled exponential integral *x e<sup>x</sup> E<sub>1</sub>(x)* |
//!
//! # Zeta functions
//!
//! | Function         | Description                                                    |
//! | ---------------- | -------------------------------------------------------------- |
//! | [`zeta`]         | Hurwitz zeta function of two arguments for real or complex `z` |
//! | [`riemann_zeta`] | Riemann zeta function for real or complex input                |
//! | [`zetac`]        | Riemann zeta function minus 1                                  |
//!
//! # Other special functions
//!
//! | Function     | Description                                                |
//! | ------------ | ---------------------------------------------------------- |
//! | [`binom`]    | Binomial coefficient                                       |
//! | [`lambertw`] | Lambert W function                                         |
//! | [`sici`]     | Sine and cosine integrals *Si(z)* and *Ci(z)*              |
//! | [`shichi`]   | Hyperbolic sine and cosine integrals *Shi(z)* and *Chi(z)* |
//! | [`spence`]   | Spence's function, also known as the dilogarithm           |
//!
//! # Convenience functions
//!
//! | Function    | Description                                         |
//! | ----------- | --------------------------------------------------- |
//! | [`cbrt`]    | *∛x*                                                |
//! | [`exp10`]   | *10<sup>x</sup>*                                    |
//! | [`exp2`]    | *2<sup>x</sup>*                                     |
//! | [`radian`]  | Convert from degrees to radians                     |
//! | [`cosdg`]   | Cosine of an angle in degrees                       |
//! | [`sindg`]   | Sine of an angle in degrees                         |
//! | [`tandg`]   | Tangent of an angle in degrees                      |
//! | [`cotdg`]   | Cotangent of an angle in degrees                    |
//! | [`log1p`]   | *log(1+x)*                                          |
//! | [`expm1`]   | *e<sup>x</sup> - 1*                                 |
//! | [`cosm1`]   | *cos(x) - 1*                                        |
//! | [`round`]   | Round to nearest or even integer-valued float       |
//! | [`xlogy`]   | *x log(y)* or *0* if *x = 0*                        |
//! | [`xlog1py`] | *x log(1+y)* or *0* if *x = 0*                      |
//! | [`exprel`]  | Relative error exponential, *(e<sup>x</sup> - 1)/x* |
//!

#![cfg_attr(not(test), no_std)]
#![warn(
    missing_debug_implementations,
    missing_docs,
    rust_2018_idioms,
    unreachable_pub
)]

extern crate alloc;

#[cfg(test)]
mod macros;
#[cfg(test)]
mod xsref;

mod ffi;
mod scipy_special;
mod utils;
mod xsf;

pub use scipy_special::*;
pub use xsf::cephes::*;
pub use xsf::*;
