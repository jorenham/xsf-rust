# xsf-rust

[![GitHub License](https://img.shields.io/github/license/jorenham/xsf-rust?style=flat-square&color=333)](https://github.com/jorenham/xsf-rust/blob/master/LICENSE)
[![Crates.io Version](https://img.shields.io/crates/v/xsf?style=flat-square&color=333333)](https://crates.io/crates/xsf)
[![docs.rs](https://img.shields.io/docsrs/xsf?style=flat-square&color=333333)](https://docs.rs/xsf/)

- docs (release): [docs.rs/xsf](https://docs.rs/xsf/)
- docs (dev): [jorenham.github.io/xsf-rust](https://jorenham.github.io/xsf-rust/xsf/)
- crate: [crates.io/crates/xsf](https://crates.io/crates/xsf)
- xsf (C++): [scipy/xsf](https://github.com/scipy/xsf/)

This crate provides pure Rust implementations of `scipy.special` functions, translated from the
[xsf](https://github.com/scipy/xsf/) C++ library that powers `scipy.special`.

All `scipy.special` functions are available, except for the orthogonal polynomial objects
(e.g. `legendre`) and the error handling (`errstate`, `geterr`, `seterr`). See
<https://docs.rs/xsf/> and the
[`scipy.special` docs](https://docs.scipy.org/doc/scipy/reference/special.html) for additional
information. [`docs/scipy.md`](https://github.com/jorenham/xsf-rust/blob/master/docs/scipy.md)
lists the functions that SciPy computes with another library than xsf.

## Functions

| `scipy.special`                                                                                                        | `xsf`                                               |
| ---------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------- |
| [`agm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.agm.html)                                   | `agm`                                               |
| [`ai_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ai_zeros.html)                         | `airy_ai_zeros`                                     |
| [`airy`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airy.html)                                 | `airy`                                              |
| [`airye`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airye.html)                               | `airy_scaled`                                       |
| [`assoc_laguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_laguerre.html)             | `eval_genlaguerre`                                  |
| [`assoc_legendre_p`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_legendre_p.html)         | `assoc_legendre_p`, `assoc_legendre_p_norm`         |
| [`assoc_legendre_p_all`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_legendre_p_all.html) | `assoc_legendre_p_all`, `assoc_legendre_p_norm_all` |
| [`bdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtr.html)                                 | `bdtr`                                              |
| [`bdtrc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrc.html)                               | `bdtrc`                                             |
| [`bdtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtri.html)                               | `bdtri`                                             |
| [`bdtrik`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrik.html)                             | `bdtrik`                                            |
| [`bdtrin`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrin.html)                             | `bdtrin`                                            |
| [`bei`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bei.html)                                   | `bei`                                               |
| [`bei_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bei_zeros.html)                       | `bei_zeros`                                         |
| [`beip`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.beip.html)                                 | `beip`                                              |
| [`beip_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.beip_zeros.html)                     | `beip_zeros`                                        |
| [`ber`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ber.html)                                   | `ber`                                               |
| [`ber_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ber_zeros.html)                       | `ber_zeros`                                         |
| [`bernoulli`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bernoulli.html)                       | `bernoulli`                                         |
| [`berp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.berp.html)                                 | `berp`                                              |
| [`berp_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.berp_zeros.html)                     | `berp_zeros`                                        |
| [`besselpoly`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.besselpoly.html)                     | `besselpoly`                                        |
| [`beta`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.beta.html)                                 | `beta`                                              |
| [`betainc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html)                           | `betainc`                                           |
| [`betaincc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincc.html)                         | `betaincc`                                          |
| [`betainccinv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainccinv.html)                   | `betainccinv`                                       |
| [`betaincinv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincinv.html)                     | `betaincinv`                                        |
| [`betaln`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaln.html)                             | `betaln`                                            |
| [`bi_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bi_zeros.html)                         | `airy_bi_zeros`                                     |
| [`binom`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.binom.html)                               | `binom`                                             |
|                                                                                                                        | `bivariate_normal_cdf`                              |
| [`boxcox`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.boxcox.html)                             | `boxcox`                                            |
| [`boxcox1p`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.boxcox1p.html)                         | `boxcox1p`                                          |
| [`btdtria`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.btdtria.html)                           | `btdtria`                                           |
| [`btdtrib`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.btdtrib.html)                           | `btdtrib`                                           |
| [`cbrt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cbrt.html)                                 | `cbrt`                                              |
|                                                                                                                        | `cdf_cvm`                                           |
|                                                                                                                        | `cdf_cvm_inf`                                       |
| [`chdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtr.html)                               | `chdtr`                                             |
| [`chdtrc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtrc.html)                             | `chdtrc`                                            |
| [`chdtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtri.html)                             | `chdtri`                                            |
| [`chdtriv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtriv.html)                           | `chdtriv`                                           |
| [`chebyc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chebyc.html)                             |                                                     |
| [`chebys`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chebys.html)                             |                                                     |
| [`chebyt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chebyt.html)                             |                                                     |
| [`chebyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chebyu.html)                             |                                                     |
| [`chndtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtr.html)                             | `chndtr`                                            |
| [`chndtridf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtridf.html)                       | `chndtridf`                                         |
| [`chndtrinc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtrinc.html)                       | `chndtrinc`                                         |
| [`chndtrix`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chndtrix.html)                         | `chndtri`                                           |
| [`comb`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.comb.html)                                 | `comb`, `comb_rep`                                  |
| [`cosdg`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cosdg.html)                               | `cosdg`                                             |
|                                                                                                                        | `cosine_cdf`                                        |
|                                                                                                                        | `cosine_invcdf`                                     |
| [`cosm1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cosm1.html)                               | `cosm1`                                             |
|                                                                                                                        | `cospi`                                             |
| [`cotdg`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cotdg.html)                               | `cotdg`                                             |
| [`dawsn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.dawsn.html)                               | `dawsn`                                             |
| [`digamma`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.digamma.html)                           | `digamma`                                           |
| [`digammainv`](https://scipy.github.io/devdocs/reference/generated/scipy.special.digammainv.html)                      | `digammainv`                                        |
| [`diric`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.diric.html)                               | `diric`                                             |
| [`ellip_harm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellip_harm.html)                     | `ellip_harm`                                        |
| [`ellip_harm_2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellip_harm_2.html)                 | `ellip_harm_2`                                      |
| [`ellip_normal`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellip_normal.html)                 | `ellip_normal`                                      |
| [`ellipe`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html)                             | `ellipe`                                            |
| [`ellipeinc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipeinc.html)                       | `ellipeinc`                                         |
| [`ellipj`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipj.html)                             | `ellipj`                                            |
| [`ellipk`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html)                             | `ellipk`                                            |
| [`ellipkinc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkinc.html)                       | `ellipkinc`                                         |
| [`ellipkm1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkm1.html)                         | `ellipkm1`                                          |
| [`elliprc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprc.html)                           | `elliprc`                                           |
| [`elliprd`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprd.html)                           | `elliprd`                                           |
| [`elliprf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprf.html)                           | `elliprf`                                           |
| [`elliprg`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprg.html)                           | `elliprg`                                           |
| [`elliprj`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprj.html)                           | `elliprj`                                           |
| [`entr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.entr.html)                                 | `entr`                                              |
| [`erf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf.html)                                   | `erf`                                               |
| [`erf_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf_zeros.html)                       | `erf_zeros`                                         |
| [`erfc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfc.html)                                 | `erfc`                                              |
| [`erfcinv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcinv.html)                           | `erfcinv`                                           |
| [`erfcx`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcx.html)                               | `erfcx`                                             |
| [`erfi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfi.html)                                 | `erfi`                                              |
| [`erfinv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfinv.html)                             | `erfinv`                                            |
| [`euler`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.euler.html)                               | `euler`                                             |
| [`eval_chebyc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebyc.html)                   | `eval_chebyshev_c`                                  |
| [`eval_chebys`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebys.html)                   | `eval_chebyshev_s`                                  |
| [`eval_chebyt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebyt.html)                   | `eval_chebyshev_t`                                  |
| [`eval_chebyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_chebyu.html)                   | `eval_chebyshev_u`                                  |
| [`eval_gegenbauer`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_gegenbauer.html)           | `eval_gegenbauer`                                   |
| [`eval_genlaguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_genlaguerre.html)         | `eval_genlaguerre`                                  |
| [`eval_hermite`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermite.html)                 | `eval_hermite_h`                                    |
| [`eval_hermitenorm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermitenorm.html)         | `eval_hermite_he`                                   |
| [`eval_jacobi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_jacobi.html)                   | `eval_jacobi`                                       |
| [`eval_laguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_laguerre.html)               | `eval_laguerre`                                     |
| [`eval_legendre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_legendre.html)               | `eval_legendre`, `legendre_p`                       |
| [`eval_sh_chebyt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_chebyt.html)             | `eval_sh_chebyshev_t`                               |
| [`eval_sh_chebyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_chebyu.html)             | `eval_sh_chebyshev_u`                               |
| [`eval_sh_jacobi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_jacobi.html)             | `eval_gengegenbauer`                                |
|                                                                                                                        | `eval_sh_jacobi`                                    |
| [`eval_sh_legendre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_sh_legendre.html)         | `eval_sh_legendre`                                  |
| [`exp1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp1.html)                                 | `exp1`                                              |
| [`exp10`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp10.html)                               | `exp10`                                             |
| [`exp2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp2.html)                                 | `exp2`                                              |
| [`expi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expi.html)                                 | `expi`                                              |
| [`expit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expit.html)                               | `expit`                                             |
| [`expm1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expm1.html)                               | `expm1`                                             |
| [`expn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expn.html)                                 | `expn`                                              |
| [`exprel`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exprel.html)                             | `exprel`                                            |
| [`factorial`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.factorial.html)                       | `factorial`, `factorial_checked`                    |
| [`factorial2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.factorial2.html)                     | `multifactorial`, `multifactorial_checked`          |
| [`factorialk`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.factorialk.html)                     | `multifactorial`, `multifactorial_checked`          |
| [`fdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtr.html)                                 | `fdtr`                                              |
| [`fdtrc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtrc.html)                               | `fdtrc`                                             |
| [`fdtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtri.html)                               | `fdtri`                                             |
| [`fdtridfd`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtridfd.html)                         | `fdtridfd`                                          |
| [`fresnel`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fresnel.html)                           | `fresnel`                                           |
| [`fresnel_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fresnel_zeros.html)               | `fresnel_zeros`                                     |
| [`fresnelc_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fresnelc_zeros.html)             | `fresnel_zeros`                                     |
| [`fresnels_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fresnels_zeros.html)             | `fresnel_zeros`                                     |
| [`gamma`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gamma.html)                               | `gamma`                                             |
| [`gammainc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainc.html)                         | `gammainc`                                          |
| [`gammaincc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaincc.html)                       | `gammaincc`                                         |
| [`gammainccinv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainccinv.html)                 | `gammainccinv`                                      |
| [`gammaincinv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaincinv.html)                   | `gammaincinv`                                       |
| [`gammaln`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaln.html)                           | `gammaln`                                           |
| [`gammasgn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammasgn.html)                         | `gammasgn`                                          |
| [`gdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtr.html)                                 | `gdtr`                                              |
| [`gdtrc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtrc.html)                               | `gdtrc`                                             |
| [`gdtria`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtria.html)                             | `gdtria`                                            |
| [`gdtrib`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtrib.html)                             | `gdtrib`                                            |
| [`gdtrix`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtrix.html)                             | `gdtrix`                                            |
| [`gegenbauer`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gegenbauer.html)                     |                                                     |
|                                                                                                                        | `gen_harmonic`                                      |
|                                                                                                                        | `genhyperbolic_logpdf`                              |
|                                                                                                                        | `genhyperbolic_pdf`                                 |
|                                                                                                                        | `geninvgauss_logpdf`                                |
|                                                                                                                        | `geninvgauss_pdf`                                   |
| [`genlaguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.genlaguerre.html)                   |                                                     |
| [`geterr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.geterr.html)                             |                                                     |
| [`h1vp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.h1vp.html)                                 | `hankel_1_prime`                                    |
| [`h2vp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.h2vp.html)                                 | `hankel_2_prime`                                    |
|                                                                                                                        | `hahn_q`                                            |
| [`hankel1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1.html)                           | `hankel_1`                                          |
| [`hankel1e`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1e.html)                         | `hankel_1e`                                         |
| [`hankel2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2.html)                           | `hankel_2`                                          |
| [`hankel2e`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2e.html)                         | `hankel_2e`                                         |
| [`hermite`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hermite.html)                           |                                                     |
| [`hermitenorm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hermitenorm.html)                   |                                                     |
| [`huber`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.huber.html)                               | `huber`                                             |
|                                                                                                                        | `hyp0f0`                                            |
| [`hyp0f1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp0f1.html)                             | `hyp0f1`                                            |
|                                                                                                                        | `hyp1f0`                                            |
| [`hyp1f1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp1f1.html)                             | `hyp1f1`                                            |
|                                                                                                                        | `hyp1f2`                                            |
| [`hyp2f1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp2f1.html)                             | `hyp2f1`                                            |
|                                                                                                                        | `hyp2f2`                                            |
|                                                                                                                        | `hyp3f2`                                            |
| [`hyperu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyperu.html)                             | `hypu`                                              |
|                                                                                                                        | `hyppfq`                                            |
| [`i0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0.html)                                     | `bessel_i0`                                         |
| [`i0e`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0e.html)                                   | `bessel_i0e`                                        |
| [`i1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1.html)                                     | `bessel_i1`                                         |
| [`i1e`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1e.html)                                   | `bessel_i1e`                                        |
| [`inv_boxcox`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.inv_boxcox.html)                     | `inv_boxcox`                                        |
| [`inv_boxcox1p`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.inv_boxcox1p.html)                 | `inv_boxcox1p`                                      |
| [`it2i0k0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2i0k0.html)                           | `it2i0k0`                                           |
| [`it2j0y0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2j0y0.html)                           | `it2j0y0`                                           |
| [`it2struve0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2struve0.html)                     | `it2struve0`                                        |
| [`itairy`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itairy.html)                             | `airy_integrals`                                    |
| [`iti0k0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iti0k0.html)                             | `it1i0k0`                                           |
| [`itj0y0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itj0y0.html)                             | `it1j0y0`                                           |
| [`itmodstruve0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itmodstruve0.html)                 | `itmodstruve0`                                      |
| [`itstruve0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itstruve0.html)                       | `itstruve0`                                         |
| [`iv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iv.html)                                     | `bessel_i`                                          |
|                                                                                                                        | `iv_ratio`                                          |
|                                                                                                                        | `iv_ratio_c`                                        |
|                                                                                                                        | `iv_ratioinv`                                       |
| [`ive`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ive.html)                                   | `bessel_ie`                                         |
| [`ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ivp.html)                                   | `bessel_i_prime`                                    |
| [`j0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j0.html)                                     | `bessel_j0`                                         |
| [`j1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j1.html)                                     | `bessel_j1`                                         |
| [`jacobi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jacobi.html)                             |                                                     |
| [`jn_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jn_zeros.html)                         | `bessel_zeros`                                      |
| [`jnjnp_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnjnp_zeros.html)                   | `bessel_j_jprime_zeros`                             |
| [`jnp_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnp_zeros.html)                       | `bessel_zeros`                                      |
| [`jnyn_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnyn_zeros.html)                     | `bessel_zeros`                                      |
| [`jv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html)                                     | `bessel_j`                                          |
| [`jve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jve.html)                                   | `bessel_je`                                         |
| [`jvp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jvp.html)                                   | `bessel_j_prime`                                    |
| [`k0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0.html)                                     | `bessel_k0`                                         |
| [`k0e`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0e.html)                                   | `bessel_k0e`                                        |
| [`k1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1.html)                                     | `bessel_k1`                                         |
| [`k1e`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1e.html)                                   | `bessel_k1e`                                        |
| [`kei`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kei.html)                                   | `kei`                                               |
| [`kei_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kei_zeros.html)                       | `kei_zeros`                                         |
| [`keip`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.keip.html)                                 | `keip`                                              |
| [`keip_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.keip_zeros.html)                     | `keip_zeros`                                        |
| [`kelvin`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin.html)                             | `kelvin`                                            |
| [`kelvin_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin_zeros.html)                 | `kelvin_zeros`                                      |
| [`ker`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ker.html)                                   | `ker`                                               |
| [`ker_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ker_zeros.html)                       | `ker_zeros`                                         |
| [`kerp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kerp.html)                                 | `kerp`                                              |
| [`kerp_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kerp_zeros.html)                     | `kerp_zeros`                                        |
| [`kl_div`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kl_div.html)                             | `kl_div`                                            |
| [`kn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kn.html)                                     | `bessel_k`                                          |
|                                                                                                                        | `kolmogc`                                           |
|                                                                                                                        | `kolmogci`                                          |
| [`kolmogi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kolmogi.html)                           | `kolmogi`                                           |
| [`kolmogorov`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kolmogorov.html)                     | `kolmogorov`                                        |
|                                                                                                                        | `kolmogp`                                           |
| [`kv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kv.html)                                     | `bessel_k`                                          |
| [`kve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kve.html)                                   | `bessel_ke`                                         |
| [`kvp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kvp.html)                                   | `bessel_k_prime`                                    |
| [`laguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.laguerre.html)                         |                                                     |
| [`lambertw`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lambertw.html)                         | `lambertw`                                          |
| [`legendre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre.html)                         |                                                     |
| [`legendre_p`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre_p.html)                     | `legendre_p`                                        |
| [`legendre_p_all`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre_p_all.html)             | `legendre_p_all`                                    |
| [`lmbda`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lmbda.html)                               | `jahnke_emden_lambda`                               |
| [`log1p`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log1p.html)                               | `log1p`                                             |
|                                                                                                                        | `log1pmx`                                           |
| [`log_expit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_expit.html)                       | `log_expit`                                         |
| [`log_ndtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_ndtr.html)                         | `log_ndtr`                                          |
| [`log_softmax`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_softmax.html)                   | `log_softmax`                                       |
| [`log_wright_bessel`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_wright_bessel.html)       | `log_wright_bessel`                                 |
|                                                                                                                        | `logaddexp2`                                        |
| [`loggamma`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.loggamma.html)                         | `loggamma`                                          |
| [`logit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logit.html)                               | `logit`                                             |
| [`logsumexp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html)                       | `log1mexp`, `logaddexp`, `logsumexp`                |
| [`lpmv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lpmv.html)                                 | `assoc_legendre_p`                                  |
| [`lqmn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lqmn.html)                                 | `assoc_legendre_q_all`                              |
| [`lqn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lqn.html)                                   | `legendre_q_all`                                    |
| [`mathieu_a`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_a.html)                       | `mathieu_a`                                         |
| [`mathieu_b`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_b.html)                       | `mathieu_b`                                         |
| [`mathieu_cem`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_cem.html)                   | `mathieu_cem`                                       |
| [`mathieu_even_coef`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_even_coef.html)       | `mathieu_even_coef`                                 |
| [`mathieu_modcem1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modcem1.html)           | `mathieu_modcem1`                                   |
| [`mathieu_modcem2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modcem2.html)           | `mathieu_modcem2`                                   |
| [`mathieu_modsem1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modsem1.html)           | `mathieu_modsem1`                                   |
| [`mathieu_modsem2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_modsem2.html)           | `mathieu_modsem2`                                   |
| [`mathieu_odd_coef`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_odd_coef.html)         | `mathieu_odd_coef`                                  |
| [`mathieu_sem`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.mathieu_sem.html)                   | `mathieu_sem`                                       |
| [`modfresnelm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modfresnelm.html)                   | `modified_fresnel_minus`                            |
| [`modfresnelp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modfresnelp.html)                   | `modified_fresnel_plus`                             |
| [`modstruve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modstruve.html)                       | `struve_l`                                          |
| [`multigammaln`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.multigammaln.html)                 | `multigammaln`                                      |
| [`nbdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtr.html)                               | `nbdtr`                                             |
| [`nbdtrc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtrc.html)                             | `nbdtrc`                                            |
| [`nbdtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtri.html)                             | `nbdtri`                                            |
| [`nbdtrik`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtrik.html)                           | `nbdtrik`                                           |
| [`nbdtrin`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtrin.html)                           | `nbdtrin`                                           |
| [`ncfdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtr.html)                             | `ncfdtr`                                            |
| [`ncfdtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtri.html)                           | `ncfdtri`                                           |
| [`ncfdtridfd`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtridfd.html)                     | `ncfdtridfd`                                        |
| [`ncfdtridfn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtridfn.html)                     | `ncfdtridfn`                                        |
| [`ncfdtrinc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ncfdtrinc.html)                       | `ncfdtrinc`                                         |
| [`nctdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtr.html)                             | `nctdtr`                                            |
| [`nctdtridf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtridf.html)                       | `nctdtridf`                                         |
| [`nctdtrinc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtrinc.html)                       | `nctdtrinc`                                         |
| [`nctdtrit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nctdtrit.html)                         | `nctdtrit`                                          |
| [`ndtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtr.html)                                 | `ndtr`                                              |
| [`ndtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri.html)                               | `ndtri`                                             |
| [`ndtri_exp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri_exp.html)                       | `ndtri_exp`                                         |
|                                                                                                                        | `normalized_gen_harmonic`                           |
| [`nrdtrimn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nrdtrimn.html)                         | `nrdtrimn`                                          |
| [`nrdtrisd`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nrdtrisd.html)                         | `nrdtrisd`                                          |
| [`obl_ang1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_ang1.html)                         | `oblate_aswfa_nocv`                                 |
| [`obl_ang1_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_ang1_cv.html)                   | `oblate_aswfa`                                      |
| [`obl_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_cv.html)                             | `oblate_segv`                                       |
| [`obl_cv_seq`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_cv_seq.html)                     | `oblate_segv_seq`                                   |
| [`obl_rad1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad1.html)                         | `oblate_radial1_nocv`                               |
| [`obl_rad1_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad1_cv.html)                   | `oblate_radial1`                                    |
| [`obl_rad2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad2.html)                         | `oblate_radial2_nocv`                               |
| [`obl_rad2_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad2_cv.html)                   | `oblate_radial2`                                    |
| [`owens_t`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.owens_t.html)                           | `owens_t`                                           |
| [`pbdn_seq`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbdn_seq.html)                         | `pbdn_seq`                                          |
| [`pbdv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbdv.html)                                 | `pbdv`                                              |
| [`pbdv_seq`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbdv_seq.html)                         | `pbdv_seq`                                          |
| [`pbvv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbvv.html)                                 | `pbvv`                                              |
| [`pbvv_seq`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbvv_seq.html)                         | `pbvv_seq`                                          |
| [`pbwa`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbwa.html)                                 | `pbwa`                                              |
| [`pdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtr.html)                                 | `pdtr`                                              |
| [`pdtrc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtrc.html)                               | `pdtrc`                                             |
| [`pdtri`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtri.html)                               | `pdtri`                                             |
| [`pdtrik`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtrik.html)                             | `pdtrik`                                            |
| [`perm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.perm.html)                                 | `perm`, `pow_falling`                               |
| [`poch`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.poch.html)                                 | `pow_rising`                                        |
| [`poisson_binom_cdf`](https://scipy.github.io/devdocs/reference/generated/scipy.special.poisson_binom_cdf.html)        | `poisson_binom_cdf`                                 |
|                                                                                                                        | `poisson_binom_cdf_all`                             |
|                                                                                                                        | `poisson_binom_pmf_all`                             |
| [`polygamma`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.polygamma.html)                       | `polygamma`, `trigamma`                             |
| [`powm1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.powm1.html)                               | `powm1`                                             |
| [`pro_ang1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_ang1.html)                         | `prolate_aswfa_nocv`                                |
| [`pro_ang1_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_ang1_cv.html)                   | `prolate_aswfa`                                     |
| [`pro_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_cv.html)                             | `prolate_segv`                                      |
| [`pro_cv_seq`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_cv_seq.html)                     | `prolate_segv_seq`                                  |
| [`pro_rad1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad1.html)                         | `prolate_radial1_nocv`                              |
| [`pro_rad1_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad1_cv.html)                   | `prolate_radial1`                                   |
| [`pro_rad2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad2.html)                         | `prolate_radial2_nocv`                              |
| [`pro_rad2_cv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad2_cv.html)                   | `prolate_radial2`                                   |
| [`pseudo_huber`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pseudo_huber.html)                 | `pseudo_huber`                                      |
| [`psi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.psi.html)                                   | `digamma`                                           |
| [`radian`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.radian.html)                             | `radian`                                            |
| [`rel_entr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html)                         | `rel_entr`                                          |
| [`rgamma`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rgamma.html)                             | `rgamma`                                            |
| [`riccati_jn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_jn.html)                     | `riccati_j`                                         |
| [`riccati_yn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_yn.html)                     | `riccati_y`                                         |
| [`roots_chebyc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebyc.html)                 | `roots_chebyshev_c`                                 |
| [`roots_chebys`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebys.html)                 | `roots_chebyshev_s`                                 |
| [`roots_chebyt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebyt.html)                 | `roots_chebyshev_t`                                 |
| [`roots_chebyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_chebyu.html)                 | `roots_chebyshev_u`                                 |
| [`roots_gegenbauer`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_gegenbauer.html)         | `roots_gegenbauer`                                  |
| [`roots_genlaguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_genlaguerre.html)       | `roots_genlaguerre`                                 |
| [`roots_hermite`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_hermite.html)               | `roots_hermite_h`                                   |
| [`roots_hermitenorm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_hermitenorm.html)       | `roots_hermite_he`                                  |
| [`roots_jacobi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_jacobi.html)                 | `roots_jacobi`                                      |
| [`roots_laguerre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_laguerre.html)             | `roots_laguerre`                                    |
| [`roots_legendre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_legendre.html)             | `roots_legendre`                                    |
| [`roots_sh_chebyt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_chebyt.html)           | `roots_sh_chebyshev_t`                              |
| [`roots_sh_chebyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_chebyu.html)           | `roots_sh_chebyshev_u`                              |
| [`roots_sh_jacobi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_jacobi.html)           | `roots_gengegenbauer`                               |
|                                                                                                                        | `roots_sh_jacobi`                                   |
| [`roots_sh_legendre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.roots_sh_legendre.html)       | `roots_sh_legendre`                                 |
| [`round`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.round.html)                               | `round`                                             |
|                                                                                                                        | `scaled_exp1`                                       |
| [`seterr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.seterr.html)                             |                                                     |
| [`sh_chebyt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sh_chebyt.html)                       |                                                     |
| [`sh_chebyu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sh_chebyu.html)                       |                                                     |
| [`sh_jacobi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sh_jacobi.html)                       |                                                     |
| [`sh_legendre`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sh_legendre.html)                   |                                                     |
| [`shichi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.shichi.html)                             | `shichi`                                            |
| [`sici`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sici.html)                                 | `sici`                                              |
| [`sinc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sinc.html)                                 | `sinc`                                              |
| [`sindg`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sindg.html)                               | `sindg`                                             |
|                                                                                                                        | `sinpi`                                             |
| [`smirnov`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.smirnov.html)                           | `smirnov`                                           |
|                                                                                                                        | `smirnovc`                                          |
|                                                                                                                        | `smirnovci`                                         |
| [`smirnovi`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.smirnovi.html)                         | `smirnovi`                                          |
|                                                                                                                        | `smirnovp`                                          |
| [`softmax`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.softmax.html)                           | `softmax`                                           |
| [`softplus`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.softplus.html)                         | `softplus`                                          |
| [`spence`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spence.html)                             | `spence`                                            |
| [`sph_harm_y`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html)                     | `sph_harm_y`                                        |
| [`sph_harm_y_all`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y_all.html)             | `sph_harm_y_all`                                    |
| [`sph_legendre_p`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_legendre_p.html)             | `sph_legendre_p`                                    |
| [`sph_legendre_p_all`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_legendre_p_all.html)     | `sph_legendre_p_all`                                |
| [`spherical_in`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_in.html)                 | `sph_bessel_i`, `sph_bessel_i_prime`                |
| [`spherical_jn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_jn.html)                 | `sph_bessel_j`, `sph_bessel_j_prime`                |
| [`spherical_kn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_kn.html)                 | `sph_bessel_k`, `sph_bessel_k_prime`                |
| [`spherical_yn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_yn.html)                 | `sph_bessel_y`, `sph_bessel_y_prime`                |
| [`stdtr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtr.html)                               | `stdtr`                                             |
| [`stdtridf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtridf.html)                         | `stdtridf`                                          |
| [`stdtrit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtrit.html)                           | `stdtri`                                            |
| [`stirling2`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stirling2.html)                       | `stirling2`                                         |
| [`struve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.struve.html)                             | `struve_h`                                          |
| [`tandg`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.tandg.html)                               | `tandg`                                             |
| [`tklmbda`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.tklmbda.html)                           | `tukeylambdacdf`                                    |
|                                                                                                                        | `trivariate_normal_cdf`                             |
| [`voigt_profile`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.voigt_profile.html)               | `voigt_profile`                                     |
|                                                                                                                        | `von_mises_cdf`                                     |
|                                                                                                                        | `wilcoxon_cdf`                                      |
|                                                                                                                        | `wilcoxon_cdf_all`                                  |
|                                                                                                                        | `wilcoxon_pmf`                                      |
|                                                                                                                        | `wilcoxon_pmf_all`                                  |
|                                                                                                                        | `wilcoxon_sf`                                       |
| [`wofz`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wofz.html)                                 | `wofz`                                              |
| [`wright_bessel`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wright_bessel.html)               | `wright_bessel`                                     |
| [`wrightomega`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wrightomega.html)                   | `wrightomega`                                       |
| [`xlog1py`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.xlog1py.html)                           | `xlog1py`                                           |
| [`xlogy`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.xlogy.html)                               | `xlogy`                                             |
| [`y0`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y0.html)                                     | `bessel_y0`                                         |
| [`y0_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y0_zeros.html)                         | `bessel_y0_zeros`                                   |
| [`y1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1.html)                                     | `bessel_y1`                                         |
| [`y1_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1_zeros.html)                         | `bessel_y1_zeros`                                   |
| [`y1p_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1p_zeros.html)                       | `bessel_y1_prime_zeros`                             |
| [`yn`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yn.html)                                     | `bessel_y`                                          |
| [`yn_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yn_zeros.html)                         | `bessel_zeros`                                      |
| [`ynp_zeros`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ynp_zeros.html)                       | `bessel_zeros`                                      |
| [`yv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yv.html)                                     | `bessel_y`                                          |
| [`yve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yve.html)                                   | `bessel_ye`                                         |
| [`yvp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yvp.html)                                   | `bessel_y_prime`                                    |
| [`zeta`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html)                                 | `riemann_zeta`, `zeta`                              |
| [`zetac`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zetac.html)                               | `zetac`                                             |
