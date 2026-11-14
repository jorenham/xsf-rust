# `scipy.special` without xsf

SciPy [v1.18.1](https://github.com/scipy/scipy/tree/v1.18.1/scipy/special) with xsf
[`f7b85f5`](https://github.com/scipy/xsf/tree/f7b85f505a94fd398024478d2e747b01dd2d7a6d): the 138
of 353 functions whose algorithm isn't in xsf.

| `scipy.special`                         | SciPy                  | xsf-rust                   |
| --------------------------------------- | ---------------------- | -------------------------- |
| `betainc`                               | Boost                  | same                       |
| `betaincc`, `betainccinv`               | Boost                  | same                       |
| `betaincinv`                            | Boost                  | same                       |
| `erfinv`                                | Boost                  | same                       |
| `fdtr`                                  | Boost                  | same                       |
| `fdtrc`                                 | Boost                  | same                       |
| `fdtri`                                 | Boost                  | same                       |
| `hyp1f1`                                | Boost; complex: xsf    | same                       |
| `powm1`                                 | Boost                  | same                       |
| `stdtr`                                 | Boost                  | same                       |
| `stdtrit`                               | Boost                  | `stdtri`                   |
| `chndtr`                                | Boost                  | same                       |
| `bdtrik`, `bdtrin`                      | Boost                  |                            |
| `btdtria`, `btdtrib`                    | Boost                  |                            |
| `chdtriv`                               | Boost                  |                            |
| `chndtrix`                              | Boost                  |                            |
| `chndtridf`, `chndtrinc`                | Boost                  |                            |
| `nbdtrik`, `nbdtrin`                    | Boost                  |                            |
| `ncfdtr`, `ncfdtri`                     | Boost                  |                            |
| `nctdtr`, `nctdtrit`                    | Boost                  |                            |
| `pdtrik`                                | Boost                  |                            |
| `fdtridfd`                              | cdflib                 | same                       |
| `ncfdtridfd`, `ncfdtridfn`, `ncfdtrinc` | cdflib                 | same                       |
| `nctdtridf`, `nctdtrinc`                | cdflib                 | same                       |
| `stdtridf`                              | cdflib                 | same                       |
| `ellipr{c,d,f,g,j}`                     | C++                    | same                       |
| `stirling2`                             | C++; Python if `exact` | same                       |
| `wrightomega`                           | C++                    | same                       |
| `agm`                                   | Cython                 | same                       |
| `boxcox`, `boxcox1p`                    | Cython                 | same                       |
| `inv_boxcox`, `inv_boxcox1p`            | Cython                 | same                       |
| `comb`                                  | Cython if `exact`; xsf | same                       |
| `entr`, `rel_entr`, `kl_div`            | Cython                 | same                       |
| other `eval_*`                          | Cython                 | same                       |
| `eval_cheby{c,s,t,u}`                   | Cython                 | `eval_chebyshev_{c,s,t,u}` |
| `eval_sh_cheby{t,u}`                    | Cython                 | `eval_sh_chebyshev_{t,u}`  |
| `eval_hermite`                          | Cython                 | `eval_hermite_h`           |
| `eval_hermitenorm`                      | Cython                 | `eval_hermite_he`          |
| `eval_sh_jacobi`                        | Cython                 | `eval_gengegenbauer`       |
| `huber`, `pseudo_huber`                 | Cython                 | same                       |
| `hyp0f1`                                | Cython                 | same                       |
| `ndtri_exp`                             | Cython                 | same                       |
| `ellip_harm`, `ellip_harm_2`            | Cython                 | same                       |
| `ellip_normal`                          | Cython                 | same                       |
| `spence`                                | xsf; complex: Cython   | same                       |
| `assoc_laguerre`                        | Python                 | `eval_genlaguerre`         |
| `diric`                                 | Python                 | same                       |
| `factorial`                             | Python                 | same                       |
| `factorial2`, `factorialk`              | Python                 | `multifactorial`           |
| `jvp`, `yvp`, `ivp`, `kvp`              | Python                 | `bessel_{j,y,i,k}_prime`   |
| `h1vp`, `h2vp`                          | Python                 | `hankel_{1,2}_prime`       |
| `logsumexp`, `softplus`                 | Python                 | same                       |
| `softmax`, `log_softmax`                | Python                 | same                       |
| `multigammaln`                          | Python                 | same                       |
| `perm`                                  | Python                 | same                       |
| `polygamma`                             | Python                 | same                       |
| `roots_{jacobi,gegenbauer,legendre}`    | Python                 | same                       |
| `roots_{laguerre,genlaguerre}`          | Python                 | same                       |
| `roots_sh_legendre`                     | Python                 | same                       |
| `roots_cheby{c,s,t,u}`                  | Python                 | `roots_chebyshev_*`        |
| `roots_sh_cheby{t,u}`                   | Python                 | `roots_sh_chebyshev_*`     |
| `roots_hermite`                         | Python                 | `roots_hermite_h`          |
| `roots_hermitenorm`                     | Python                 | `roots_hermite_he`         |
| `roots_sh_jacobi`                       | Python                 | `roots_gengegenbauer`      |
| aliases of `roots_*`, e.g. `p_roots`    | Python                 |                            |
| orthogonal polynomial classes           | Python                 |                            |
| `sinc`                                  | NumPy                  | same                       |
