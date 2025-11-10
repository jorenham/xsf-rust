//! Bindings to the [scipy/xsf](https://github.com/scipy/xsf) C++ library that powers
//! `scipy.special`.
//! See the [`scipy.special` documentation](https://docs.scipy.org/doc/scipy/reference/special.html)
//! for additional information.
//!
//! # Airy functions
//!
//! | Function           | Description                                                      |
//! | ------------------ | ---------------------------------------------------------------- |
//! | [`airy`]           | Airy functions and derivatives                                   |
//! | [`airy_scaled`]    | Exponentially scaled Airy functions and derivatives              |
//! | [`airy_ai_zeros`]  | Zeros and values of the Airy function Ai and its derivative      |
//! | [`airy_bi_zeros`]  | Zeros and values of the Airy function Bi and its derivative      |
//! | [`airy_integrals`] | Integrals of Airy functions                                      |
//!
//! # Elliptic functions and integrals
//!
//! | Function         | Description                                                 |
//! | ---------------- | ----------------------------------------------------------- |
//! | [`ellipj`]       | Jacobian elliptic functions                                 |
//! | [`ellipk`]       | Complete elliptic integral of the first kind                |
//! | [`ellipkm1`]     | Complete elliptic integral of the first kind around $m = 1$ |
//! | [`ellipkinc`]    | Incomplete elliptic integral of the first kind              |
//! | [`ellipe`]       | Complete elliptic integral of the second kind               |
//! | [`ellipeinc`]    | Incomplete elliptic integral of the second kind             |
//!
//! # Bessel functions
//!
//! | Function                | Description                                                       |
//! | ----------------------- | ----------------------------------------------------------------- |
//! | [`bessel_j`]            | Bessel function of the first kind, $J_v(z)$                       |
//! | [`bessel_je`]           | Exponentially scaled Bessel function of the first kind            |
//! | [`bessel_y`]            | Bessel function of the second kind, $Y_v(z)$                      |
//! | [`bessel_ye`]           | Exponentially scaled Bessel function of the second kind           |
//! | [`bessel_i`]            | Modified Bessel function of the first kind, $I_v(z)$              |
//! | [`bessel_ie`]           | Exponentially scaled modified Bessel function of the first kind   |
//! | [`bessel_k`]            | Modified Bessel function of the second kind, $K_v(z)$             |
//! | [`bessel_ke`]           | Exponentially scaled modified Bessel function of the second kind  |
//! | [`hankel_1`]            | Hankel function of the first kind, $H_v^{(1)}(z)$                 |
//! | [`hankel_1e`]           | Exponentially scaled Hankel function of the first kind            |
//! | [`hankel_2`]            | Hankel function of the second kind, $H_v^{(2)}(z)$                |
//! | [`hankel_2e`]           | Exponentially scaled Hankel function of the second kind           |
//! | [`wright_bessel`]       | Wright's generalized Bessel function                              |
//! | [`log_wright_bessel`]   | Natural logarithm of Wright's generalized Bessel function         |
//! | [`jahnke_emden_lambda`] | Jahnke-Emden Lambda function $\Lambda_{\nu}(x)$ and derivatives   |
//!
//! ## Zeros of Bessel functions
//!
//! | Function         | Description                                                              |
//! | ---------------- | ------------------------------------------------------------------------ |
//! | [`bessel_zeros`] | Zeros of Bessel functions $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$ |
//!
//! ## Faster versions of common Bessel functions
//!
//! | Function       | Description                                                                 |
//! | -------------- | --------------------------------------------------------------------------- |
//! | [`bessel_j0`]  | Bessel function of the first kind of order 0, $J_0(x)$                      |
//! | [`bessel_j1`]  | Bessel function of the first kind of order 1, $J_1(x)$                      |
//! | [`bessel_y0`]  | Bessel function of the second kind of order 0, $Y_0(x)$                     |
//! | [`bessel_y1`]  | Bessel function of the second kind of order 1, $Y_1(x)$                     |
//! | [`bessel_i0`]  | Modified Bessel function of the first kind of order 0, $I_0(x)$             |
//! | [`bessel_i0e`] | Exponentially scaled modified Bessel function of the first kind of order 0  |
//! | [`bessel_i1`]  | Modified Bessel function of the first kind of order 1, $I_1(x)$             |
//! | [`bessel_i1e`] | Exponentially scaled modified Bessel function of the first kind of order 1  |
//! | [`bessel_k0`]  | Modified Bessel function of the second kind of order 0, $K_0(x)$            |
//! | [`bessel_k0e`] | Exponentially scaled modified Bessel function of the second kind of order 0 |
//! | [`bessel_k1`]  | Modified Bessel function of the second kind of order 1, $K_1(x)$            |
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
//! | Function           | Description                       |
//! | ------------------ | --------------------------------- |
//! | [`bessel_j_prime`] | $n$-th derivative of [`bessel_j`] |
//! | [`bessel_y_prime`] | $n$-th derivative of [`bessel_y`] |
//! | [`bessel_i_prime`] | $n$-th derivative of [`bessel_i`] |
//! | [`bessel_k_prime`] | $n$-th derivative of [`bessel_k`] |
//! | [`hankel_1_prime`] | $n$-th derivative of [`hankel_1`] |
//! | [`hankel_2_prime`] | $n$-th derivative of [`hankel_2`] |
//!
//! ## Spherical Bessel functions
//!
//! | Function               | Description                                                     |
//! | ---------------------- | --------------------------------------------------------------- |
//! | [`sph_bessel_j`]       | Spherical Bessel function of the first kind, $j_n(z)$           |
//! | [`sph_bessel_j_prime`] | Derivative of [`sph_bessel_j`], $j_n\'(z)$                      |
//! | [`sph_bessel_y`]       | Spherical Bessel function of the second kind, $y_n(z)$          |
//! | [`sph_bessel_y_prime`] | Derivative of [`sph_bessel_y`], $y_n\'(z)$                      |
//! | [`sph_bessel_i`]       | Modified Spherical Bessel function of the first kind, $i_n(z)$  |
//! | [`sph_bessel_i_prime`] | Derivative of [`sph_bessel_i`], $i_n\'(z)$                      |
//! | [`sph_bessel_k`]       | Modified Spherical Bessel function of the second kind, $k_n(z)$ |
//! | [`sph_bessel_k_prime`] | Derivative of [`sph_bessel_k`], $k_n\'(z)$                      |
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
//! | Function         | Description                                                   |
//! | ---------------- | ------------------------------------------------------------- |
//! | [`struve_h`]     | Struve function $H_{\nu}(x)$                                  |
//! | [`struve_l`]     | Modified Struve function $L_{\nu}(x)$                         |
//! | [`itstruve0`]    | Integral of the Struve function of order 0, $H_0(x)$          |
//! | [`it2struve0`]   | Integral related to the Struve function of order 0            |
//! | [`itmodstruve0`] | Integral of the modified Struve function of order 0, $L_0(x)$ |
//!
//! # Raw statistical functions
//!
//! ## Binomial distribution
//!
//! | Function  | Description                      |
//! | --------- | -------------------------------- |
//! | [`bdtr`]  | Cumulative distribution function |
//! | [`bdtrc`] | Complement of [`bdtr`]           |
//! | [`bdtri`] | Inverse of [`bdtr`]              |
//!
//! ## F distribution
//!
//! | Function  | Description                      |
//! | --------- | -------------------------------- |
//! | [`fdtr`]  | Cumulative distribution function |
//! | [`fdtrc`] | Complement of [`fdtr`]           |
//! | [`fdtri`] | Inverse of [`fdtr`]              |
//!
//! ## Gamma distribution
//!
//! | Function   | Description                                            |
//! | ---------- | ------------------------------------------------------ |
//! | [`gdtr`]   | Cumulative distribution function                       |
//! | [`gdtrc`]  | Complement of [`gdtr`]                                 |
//! | [`gdtrib`] | Inverse of [`gdtr(a, b, x)`](gdtr) with respect to `b` |
//!
//! ## Negative binomial distribution
//!
//! | Function   | Description                      |
//! | ---------- | -------------------------------- |
//! | [`nbdtr`]  | Cumulative distribution function |
//! | [`nbdtrc`] | Complement of [`nbdtr`]          |
//! | [`nbdtri`] | Inverse of [`nbdtr`]             |
//!
//! ## Normal distribution
//!
//! | Function     | Description                      |
//! | ------------ | -------------------------------- |
//! | [`ndtr`]     | Cumulative distribution function |
//! | [`log_ndtr`] | Logarithm of [`ndtr`]            |
//! | [`ndtri`]    | Inverse of [`ndtr`]              |
//!
//! ## Poisson distribution
//!
//! | Function  | Description                      |
//! | --------- | -------------------------------- |
//! | [`pdtr`]  | Cumulative distribution function |
//! | [`pdtrc`] | Complement of [`pdtr`]           |
//! | [`pdtri`] | Inverse of [`pdtr`]              |
//!
//! ## Student's t distribution
//!
//! | Function   | Description                      |
//! | ---------- | -------------------------------- |
//! | [`stdtr`]  | Cumulative distribution function |
//! | [`stdtri`] | Inverse of [`stdtr`]             |
//!
//! ## Chi square distribution
//!
//! | Function   | Description                      |
//! | ---------- | -------------------------------- |
//! | [`chdtr`]  | Cumulative distribution function |
//! | [`chdtrc`] | Complement of [`chdtr`]          |
//! | [`chdtri`] | Inverse of [`chdtr`]             |
//!
//! ## Kolmogorov distribution
//!
//! | Function       | Description                      |
//! | -------------- | -------------------------------- |
//! | [`kolmogorov`] | Survival function                |
//! | [`kolmogp`]    | Derivative of [`kolmogorov`]     |
//! | [`kolmogi`]    | Inverse of [`kolmogorov`]        |
//! | [`kolmogc`]    | Complement of [`kolmogorov`]     |
//! | [`kolmogci`]   | Inverse of [`kolmogc`]           |
//!
//! ## Kolmogorov-Smirnov distribution
//!
//! | Function      | Description                      |
//! | ------------- | -------------------------------- |
//! | [`smirnov`]   | Survival function                |
//! | [`smirnovp`]  | Derivative of [`smirnov`]        |
//! | [`smirnovi`]  | Inverse of [`smirnov`]           |
//! | [`smirnovc`]  | Complement of [`smirnov`]        |
//! | [`smirnovci`] | Inverse of [`smirnovc`]          |
//!
//! ## Box-Cox transformation
//!
//! | Function         | Description                       |
//! | ---------------- | --------------------------------- |
//! | [`boxcox`]       | Box-Cox transformation of $x$     |
//! | [`boxcox1p`]     | Box-Cox transformation of $1 + x$ |
//! | [`inv_boxcox`]   | Inverse of [`boxcox`]             |
//! | [`inv_boxcox1p`] | Inverse of [`boxcox1p`]           |
//!
//! ## Sigmoidal functions
//!
//! | Function       | Description                              |
//! | -------------- | ---------------------------------------- |
//! | [`logit`]      | Logit function, $\ln \( \frac{x}{1-x} \) $ |
//! | [`expit`]      | Expit function, $\frac{1}{1 + \exp(-x)}$   |
//! | [`log_expit`]  | Logarithm of [`expit`]                   |
//!
//! ## Miscellaneous
//!
//! | Function           | Description                                   |
//! | ------------------ | --------------------------------------------- |
//! | [`tukeylambdacdf`] | Tukey-Lambda cumulative distribution function |
//! | [`owens_t`]        | Owen's T function                             |
//!
//! # Information Theory functions
//!
//! | Function         | Description                                                            |
//! | ---------------- | ---------------------------------------------------------------------- |
//! | [`entr`]         | Elementwise function for computing entropy, $H\[X\]$                   |
//! | [`rel_entr`]     | Elementwise function for computing relative entropy, $H\[X \rvert Y\]$ |
//! | [`kl_div`]       | Elementwise function for computing Kullback-Leibler divergence         |
//! | [`huber`]        | Huber loss function, $L_\delta(r)$                                     |
//! | [`pseudo_huber`] | Pseudo-Huber loss function, $\widetilde{L}_\delta(r)$                  |
//!
//! # Gamma and related functions
//!
//! | Function         | Description                                                       |
//! | ---------------- | ----------------------------------------------------------------- |
//! | [`gamma`]        | Gamma function, $\Gamma(z)$                                       |
//! | [`gammaln`]      | Log-gamma function, $\ln\abs{\Gamma(z)}$                          |
//! | [`loggamma`]     | Principal branch of $\ln \Gamma(z)$                               |
//! | [`gammasgn`]     | Sign of [`gamma`], $\sgn(\Gamma(z))$                              |
//! | [`gammainc`]     | Regularized lower incomplete gamma function $P(a,x) = 1 - Q(a,x)$ |
//! | [`gammaincinv`]  | Inverse of [`gammainc`], $P^{-1}(a,y)$                            |
//! | [`gammaincc`]    | Regularized upper incomplete gamma function $Q(a,x) = 1 - P(a,x)$ |
//! | [`gammainccinv`] | Inverse of [`gammaincc`], $Q^{-1}(a,y)$                           |
//! | [`beta`]         | Beta function, $\B(a,b) = {\Gamma(a)\Gamma(b) \over \Gamma(a+b)}$ |
//! | [`betaln`]       | Log-Beta function, $\ln\abs{\B(a,b)}$                             |
//! | [`betainc`]      | Regularized incomplete beta function, $\I_x(a,b)$                 |
//! | [`betaincinv`]   | Inverse of [`betainc`], $\I_y^{-1}(a,b)$                          |
//! | [`digamma`]      | The digamma function, $\psi(z)$                                   |
//! | [`polygamma`]    | The polygamma function, $\psi^{(n)}(x)$                           |
//! | [`rgamma`]       | Reciprocal of the gamma function, $\frac{1}{\Gamma(z)}$           |
//! | [`pow_rising`]   | Rising factorial $\rpow x m = {\Gamma(x+m) \over \Gamma(x)}$      |
//! | [`pow_falling`]  | Falling factorial $\fpow x m = {\Gamma(x+1) \over \Gamma(x+1-m)}$ |
//!
//! # Error function and Fresnel integrals
//!
//! | Function                   | Description                                                     |
//! | -------------------------- | --------------------------------------------------------------- |
//! | [`erf`]                    | Error function, $\erf(z)$                                       |
//! | [`erfc`]                   | Complementary error function, $\erfc(z) = 1 - \erf(z)$          |
//! | [`erfcx`]                  | Scaled complementary error function, $e^{z^2} \erfc(z)$         |
//! | [`erfi`]                   | Imaginary error function $\erfi(z) = -i \erf(i z)$              |
//! | [`erfinv`]                 | Inverse of [`erf`], $\erf^{-1}(z)$                              |
//! | [`erfcinv`]                | Inverse of [`erfc`], $\erfc^{-1}(z) = \erf^{-1}(1 - z)$         |
//! | [`erf_zeros`]              | Zeros (roots) of [`erf`]                                        |
//! | [`wofz`]                   | Faddeeva function, $w(z) = e^{-z^2} \erfc(-iz)$                 |
//! | [`dawsn`]                  | Dawson function $D(z) = \frac{\sqrt{\pi}}{2} e^{-z^2} \erfi(z)$ |
//! | [`fresnel`]                | Fresnel integrals $S(z)$ and $C(z)$                             |
//! | [`fresnel_zeros`]          | Zeros (roots) of Fresnel integrals $S(z)$ and $C(z)$            |
//! | [`modified_fresnel_plus`]  | Modified Fresnel positive integrals                             |
//! | [`modified_fresnel_minus`] | Modified Fresnel negative integrals                             |
//! | [`voigt_profile`]          | Voigt profile                                                   |
//!
//! # Legendre functions
//!
//! | Function                      | Description                                                  |
//! | ----------------------------- | ------------------------------------------------------------ |
//! | [`legendre_p`]                | Legendre polynomial of the first kind, $P_n(z)$              |
//! | [`legendre_p_all`]            | All Legendre polynomials of the first kind                   |
//! | [`assoc_legendre_p`]          | Associated Legendre polynomial of the 1st kind, $P_n^m(z)$   |
//! | [`assoc_legendre_p_all`]      | All associated Legendre polynomials of the 1st kind          |
//! | [`assoc_legendre_p_norm`]     | Normalized associated Legendre polynomial                    |
//! | [`assoc_legendre_p_norm_all`] | All normalized associated Legendre polynomials               |
//! | [`sph_legendre_p`]            | Spherical Legendre polynomial of the first kind              |
//! | [`sph_legendre_p_all`]        | All spherical Legendre polynomials of the first kind         |
//! | [`sph_harm_y`]                | Spherical harmonics, $Y_n^m(\theta,\phi)$                    |
//! | [`sph_harm_y_all`]            | All spherical harmonics                                      |
//! | [`legendre_q_all`]       | All Legendre functions of the 2nd kind and derivatives            |
//! | [`assoc_legendre_q_all`] | All associated Legendre functions of the 2nd kind and derivatives |
//!
//! # Orthogonal polynomials
//!
//! The following functions evaluate values of orthogonal polynomials:
//!
//! | Function          | Name     | Notation                   |
//! | ----------------- | -------- | -------------------------- |
//! | [`eval_legendre`] | Legendre | $P_n(z)$                   |
//! | [`eval_jacobi`]   | Jacobi   | $P_n^{(\alpha, \beta)}(z)$ |
//!
//! # Hypergeometric functions
//!
//! | Function   | Description                             | Notation                         |
//! | ---------- | --------------------------------------- | -------------------------------- |
//! | [`hyp0f0`] | Generalized hypergeometric function     | $_0F_0\left[ \middle\| z\right]$ |
//! | [`hyp1f0`] | Generalized hypergeometric function     | $_1F_0\left[a\middle\| z\right]$ |
//! | [`hyp0f1`] | Confluent hypergeometric limit function | $_0F_1\left[b\middle\| z\right]$ |
//! | [`hypu`]   | Confluent hypergeometric function       | $U(a,b,x)$                       |
//! | [`hyp1f1`] | Confluent hypergeometric function       | $\hyp 1 1 a b z$                 |
//! | [`hyp2f1`] | Gauss' hypergeometric function          | $\hyp 2 1 {a_1\enspace a_2} b z$ |
//!
//! # Parabolic cylinder functions
//!
//! | Function | Description                                    |
//! | -------- | ---------------------------------------------- |
//! | [`pbdv`] | Parabolic cylinder function $D_v(x)$ and its derivative $D_v\'(x)$ |
//! | [`pbvv`] | Parabolic cylinder function $V_v(x)$ and its derivative $V_v\'(x)$ |
//! | [`pbwa`] | Parabolic cylinder function $W_a(x)$ and its derivative $W_a\'(x)$ |
//!
//! # Mathieu and related functions
//!
//! | Even                  | Odd                  | Description                                   |
//! | --------------------- | -------------------- | --------------------------------------------- |
//! | [`mathieu_a`]         | [`mathieu_b`]        | Characteristic value of the Mathieu functions |
//! | [`mathieu_cem`]       | [`mathieu_sem`]      | Mathieu functions                             |
//! | [`mathieu_modcem1`]   | [`mathieu_modsem1`]  | Modified Mathieu functions of the first kind  |
//! | [`mathieu_modcem2`]   | [`mathieu_modsem2`]  | Modified Mathieu functions of the second kind |
//! | [`mathieu_even_coef`] | [`mathieu_odd_coef`] | Fourier coefficients for Mathieu functions    |
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
//! | Function   | Zeros            | Description                         |
//! | ---------- | ---------------- | ----------------------------------- |
//! | [`kelvin`] | [`kelvin_zeros`] | Kelvin functions as complex numbers |
//! | [`ber`]    | [`ber_zeros`]    | Kelvin function $\ber(x)$           |
//! | [`berp`]   | [`berp_zeros`]   | Derivative of [`ber`], $\ber\'(x)$   |
//! | [`bei`]    | [`bei_zeros`]    | Kelvin function $\bei(x)$           |
//! | [`beip`]   | [`beip_zeros`]   | Derivative of [`bei`], $\bei\'(x)$   |
//! | [`ker`]    | [`ker_zeros`]    | Kelvin function $\ker(x)$           |
//! | [`kerp`]   | [`kerp_zeros`]   | Derivative of [`ker`], $\ker\'(x)$   |
//! | [`kei`]    | [`kei_zeros`]    | Kelvin function $\kei(x)$           |
//! | [`keip`]   | [`keip_zeros`]   | Derivative of [`kei`], $\kei\'(x)$   |
//!
//! # Combinatorics
//!
//! | Function       | Description                                                              |
//! | -------------- | ------------------------------------------------------------------------ |
//! | [`comb`]       | $k$-combinations of $n$ things, $_nC_k = {n \choose k}$                  |
//! | [`comb_rep`]   | $k$-combinations with replacement, $\big(\\!\\!{n \choose k}\\!\\!\big)$ |
//! | [`perm`]       | $k$-permutations of $n$ things, $_nP_k = {n! \over (n-k)!}$              |
//! | [`stirling2`]  | Stirling number of the second kind $S(n,k)$                              |
//!
//! # Factorials
//!
//! | Function                   | Description                               |
//! | -------------------------- | ----------------------------------------- |
//! | [`factorial`]              | Factorial $n!$                            |
//! | [`factorial_checked`]      | [`factorial`] with overflow checking      |
//! | [`multifactorial`]         | Multifactorial $n!_{(k)}$                 |
//! | [`multifactorial_checked`] | [`multifactorial`] with overflow checking |
//!
//! # Exponential integrals
//!
//! | Function        | Description                                |
//! | --------------- | ------------------------------------------ |
//! | [`expn`]        | Generalized exponential integral $E_n(x)$  |
//! | [`expi`]        | Exponential integral $Ei(x)$               |
//! | [`exp1`]        | Exponential integral $E_1(x)$              |
//! | [`scaled_exp1`] | Scaled exponential integral $x e^x E_1(x)$ |
//!
//! # Zeta functions
//!
//! | Function         | Description                                                |
//! | ---------------- | ---------------------------------------------------------- |
//! | [`zeta`]         | Hurwitz zeta function $\zeta(z,q)$ for real or complex $z$ |
//! | [`riemann_zeta`] | Riemann zeta function $\zeta(z)$ for real or complex $z$   |
//! | [`zetac`]        | $\zeta(x) - 1$ for real $x$                                |
//!
//! # Other special functions
//!
//! | Function      | Description                                                  |
//! | ------------- | ------------------------------------------------------------ |
//! | [`bernoulli`] | Bernoulli numbers $B_0,\dotsc,B_{N-1}$                       |
//! | [`binom`]     | Binomial coefficient $\binom{n}{k}$ for real input           |
//! | [`diric`]     | Periodic sinc function, also called the Dirichlet kernel     |
//! | [`euler`]     | Euler numbers $E_0,\dotsc,E_{N-1}$                           |
//! | [`lambertw`]  | Lambert W function, $W(z)$                                   |
//! | [`sici`]      | Sine and cosine integrals $\Si(z)$ and $\Ci(z)$              |
//! | [`shichi`]    | Hyperbolic sine and cosine integrals $\Shi(z)$ and $\Chi(z)$ |
//! | [`spence`]    | Spence's function, also known as the dilogarithm             |
//! | [`softplus`]  | $\ln(1 + e^x)$                                               |
//!
//! # Convenience functions
//!
//! | Function       | Description                                           |
//! | -------------- | ----------------------------------------------------- |
//! | [`cbrt`]       | $\sqrt\[3\]{x}$                                       |
//! | [`exp10`]      | $10^x$                                                |
//! | [`exp2`]       | $2^x$                                                 |
//! | [`radian`]     | Convert from degrees to radians                       |
//! | [`cosdg`]      | Cosine of an angle in degrees                         |
//! | [`sindg`]      | Sine of an angle in degrees                           |
//! | [`tandg`]      | Tangent of an angle in degrees                        |
//! | [`cotdg`]      | Cotangent of an angle in degrees                      |
//! | [`log1p`]      | $\ln(1+x)$                                            |
//! | [`expm1`]      | $e^x - 1$                                             |
//! | [`cosm1`]      | $\cos(x) - 1$                                         |
//! | [`round`]      | Round to nearest or even integer-valued float         |
//! | [`xlogy`]      | $x \ln(y)$ or $0$ if $x = 0$                          |
//! | [`xlog1py`]    | $x \ln(1+y)$ or $0$ if $x = 0$                        |
//! | [`logaddexp`]  | $\ln(e^x + e^y)$                                      |
//! | [`logaddexp2`] | $\log_2(2^x + 2^y)$                                   |
//! | [`exprel`]     | Relative error exponential, $e^x - 1 \over x$         |
//! | [`sinc`]       | Normalized sinc function, $\sin(\pi x) \over \pi x$   |
//!

#![warn(
    missing_debug_implementations,
    missing_docs,
    rust_2018_idioms,
    unreachable_pub
)]

#[cfg(test)]
mod macros;
#[cfg(test)]
mod xsref;

mod ffi;
mod numpy;
mod scipy_special;
mod utils;
mod xsf;

pub use numpy::*;
pub use scipy_special::*;
pub use xsf::cephes::*;
pub use xsf::*;
