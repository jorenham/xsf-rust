// Evaluates xsf (C++) functions for the accuracy table (see accuracy.py).
//
// Usage: harness < INPUT > OUTPUT
//
// Each input line is `<function> <arg>...`; each output line is the double result (%.17g), or the
// real and imaginary parts of a complex result.
// The function names are those of the xsf-rust functions in specs.py, and must match
// src/main.rs.

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "xsf/agm.h"
#include "xsf/alg.h"
#include "xsf/beta.h"
#include "xsf/binom.h"
#include "xsf/boxcox.h"
#include "xsf/cdflib.h"
#include "xsf/cephes/erfinv.h"
#include "xsf/cephes/expn.h"
#include "xsf/cephes/incbet.h"
#include "xsf/cephes/incbi.h"
#include "xsf/cephes/poch.h"
#include "xsf/cephes/unity.h"
#include "xsf/convex_analysis.h"
#include "xsf/cosine.h"
#include "xsf/cpu/airy.h"
#include "xsf/cpu/bessel.h"
#include "xsf/cpu/ellint_carlson.h"
#include "xsf/cpu/exp.h"
#include "xsf/cpu/fresnel.h"
#include "xsf/cpu/hyp0f1.h"
#include "xsf/cpu/hyperu.h"
#include "xsf/cpu/kelvin.h"
#include "xsf/cpu/legendre.h"
#include "xsf/cpu/orthogonal_eval.h"
#include "xsf/cpu/par_cyl.h"
#include "xsf/cpu/specfun.h"
#include "xsf/cpu/sph_bessel.h"
#include "xsf/cpu/sph_harm.h"
#include "xsf/cpu/stats.h"
#include "xsf/cpu/struve.h"
#include "xsf/digamma.h"
#include "xsf/digammainv.h"
#include "xsf/ellip.h"
#include "xsf/erf.h"
#include "xsf/expint.h"
#include "xsf/gamma.h"
#include "xsf/gen_harmonic.h"
#include "xsf/hyp2f1.h"
#include "xsf/iv_ratio.h"
#include "xsf/lambertw.h"
#include "xsf/log.h"
#include "xsf/log_exp.h"
#include "xsf/loggamma.h"
#include "xsf/multivariate_normal.h"
#include "xsf/ndtri_exp.h"
#include "xsf/orthogonal_eval.h"
#include "xsf/sici.h"
#include "xsf/spence.h"
#include "xsf/stats.h"
#include "xsf/stirling2.h"
#include "xsf/trig.h"
#include "xsf/wright.h"
#include "xsf/wright_bessel.h"
#include "xsf/zeta.h"

using cd = std::complex<double>;

// the (index)th output of a function `f(args..., out0, out1, ...)` with output arguments
template <int I, typename F> double nth(F f) {
    double o[4];
    f(o[0], o[1], o[2], o[3]);
    return o[I];
}

template <int I, typename F> cd nthc(F f) {
    cd o[2];
    f(o[0], o[1]);
    return o[I];
}

static bool eval(const std::string &f, const std::vector<double> &a, double &r, cd &rc, bool &cplx) {
#define F1(name, expr)                                                                             \
    if (f == name && a.size() == 1) {                                                              \
        const double x = a[0];                                                                     \
        r = (expr);                                                                                \
        return true;                                                                               \
    }
#define F2(name, expr)                                                                             \
    if (f == name && a.size() == 2) {                                                              \
        const double x = a[0], y = a[1];                                                           \
        r = (expr);                                                                                \
        return true;                                                                               \
    }
#define F3(name, expr)                                                                             \
    if (f == name && a.size() == 3) {                                                              \
        const double x = a[0], y = a[1], z = a[2];                                                 \
        r = (expr);                                                                                \
        return true;                                                                               \
    }
#define F4(name, expr)                                                                             \
    if (f == name && a.size() == 4) {                                                              \
        const double x = a[0], y = a[1], z = a[2], w = a[3];                                       \
        r = (expr);                                                                                \
        return true;                                                                               \
    }
// complex results; the complex arguments are passed as (real, imaginary) pairs
#define C(name, n, expr)                                                                           \
    if (f == name && a.size() == n) {                                                              \
        rc = (expr);                                                                               \
        cplx = true;                                                                               \
        return true;                                                                               \
    }
#define OUT2(I, call) nth<I>([&](double &o0, double &o1, double &, double &) { call; })
#define OUT4(I, call) nth<I>([&](double &o0, double &o1, double &o2, double &o3) { call; })
    const auto l = [](double x) { return static_cast<long>(x); };
    const auto i = [](double x) { return static_cast<int>(x); };
    // gamma and related
    F1("gamma", xsf::gamma(x))
    F1("gammaln", xsf::gammaln(x))
    F1("rgamma", xsf::rgamma(x))
    F1("digamma", xsf::digamma(x))
    F1("digammainv", xsf::digammainv(x))
    F2("beta", xsf::beta(x, y))
    F2("betaln", xsf::betaln(x, y))
    F2("binom", xsf::binom(x, y))
    F2("pow_rising", xsf::cephes::poch(x, y))
    // error function and normal distribution
    F1("erf", xsf::erf(x))
    F1("erfc", xsf::erfc(x))
    F1("erfcx", xsf::erfcx(x))
    F1("erfinv", xsf::cephes::erfinv(x))
    F1("erfcinv", xsf::cephes::erfcinv(x))
    F1("ndtr", xsf::ndtr(x))
    F1("ndtri", xsf::ndtri(x))
    F1("log_ndtr", xsf::log_ndtr(x))
    // elementary functions
    F1("expm1", xsf::expm1(x))
    F1("log1p", xsf::log1p(x))
    F1("cbrt", xsf::cbrt(x))
    F1("exp2", xsf::exp2(x))
    F1("exp10", xsf::exp10(x))
    F1("sinpi", xsf::sinpi(x))
    F1("cospi", xsf::cospi(x))
    F1("cosine_cdf", xsf::cosine_cdf(x))
    F1("cosine_invcdf", xsf::cosine_invcdf(x))
    // Bessel and Struve functions
    F1("bessel_j0", xsf::cyl_bessel_j0(x))
    F1("bessel_j1", xsf::cyl_bessel_j1(x))
    F1("bessel_y0", xsf::cyl_bessel_y0(x))
    F1("bessel_y1", xsf::cyl_bessel_y1(x))
    F1("bessel_i0", xsf::cyl_bessel_i0(x))
    F1("bessel_i1", xsf::cyl_bessel_i1(x))
    F1("bessel_k0", xsf::cyl_bessel_k0(x))
    F1("bessel_k1", xsf::cyl_bessel_k1(x))
    F2("bessel_j", xsf::cyl_bessel_j(x, y))
    F2("bessel_y", xsf::cyl_bessel_y(x, y))
    F2("bessel_i", xsf::cyl_bessel_i(x, y))
    F2("bessel_k", xsf::cyl_bessel_k(x, y))
    F2("struve_h", xsf::struve_h(x, y))
    F2("struve_l", xsf::struve_l(x, y))
    // zeta functions
    F2("zeta", xsf::zeta(x, y))
    F1("riemann_zeta", xsf::riemann_zeta(x))
    F1("zetac", xsf::zetac(x))
    // exponential integrals
    F1("expi", xsf::expi(x))
    F2("expn", xsf::cephes::expn(static_cast<int>(x), y))
    // incomplete gamma and beta functions
    F2("gammainc", xsf::gammainc(x, y))
    F2("gammaincc", xsf::gammaincc(x, y))
    F3("betainc", xsf::cephes::incbet(x, y, z))
    F3("fdtr", xsf::fdtr(x, y, z))
    F3("fdtrc", xsf::fdtrc(x, y, z))
    // elliptic integrals
    F1("ellipk", xsf::ellipk(x))
    F1("ellipe", xsf::ellipe(x))
    F2("ellipkinc", xsf::ellipkinc(x, y))
    F2("ellipeinc", xsf::ellipeinc(x, y))
    // other
    F1("lambertw", xsf::lambertw(std::complex<double>(x, 0.0), 0, 1e-8).real())
    F1("spence", xsf::spence(x))
    F1("wrightomega", xsf::wrightomega(x))
    F2("stirling2", xsf::stirling2(x, y))
    // gamma and related
    F1("loggamma", xsf::loggamma(x))
    F1("gammasgn", xsf::gammasgn(x))
    F1("lgam1p", xsf::cephes::lgam1p(x))
    // error function and normal distribution
    F1("erfi", xsf::erfi(x))
    F1("dawsn", xsf::dawsn(x))
    F3("voigt_profile", xsf::voigt_profile(x, y, z))
    F1("ndtri_exp", xsf::ndtri_exp(x))
    F2("owens_t", xsf::owens_t(x, y))
    F3("bivariate_normal_cdf", xsf::bivariate_normal_cdf(x, y, z))
    if (f == "trivariate_normal_cdf" && a.size() == 6) {
        r = xsf::trivariate_normal_cdf(a[0], a[1], a[2], a[3], a[4], a[5], 1e-14);
        return true;
    }
    F3("nrdtrimn", xsf::nrdtrimn(x, y, z))
    F3("nrdtrisd", xsf::nrdtrisd(x, y, z))
    // elementary functions
    F1("cosm1", xsf::cosm1(x))
    F1("sindg", xsf::sindg(x))
    F1("cosdg", xsf::cosdg(x))
    F1("tandg", xsf::tandg(x))
    F1("cotdg", xsf::cotdg(x))
    F3("radian", xsf::radian(x, y, z))
    F2("agm", xsf::agm(x, y))
    // exponential and logarithmic functions
    F1("expit", xsf::expit(x))
    F1("logit", xsf::logit(x))
    F1("log_expit", xsf::log_expit(x))
    F1("exprel", xsf::exprel(x))
    F1("log1mexp", xsf::log1mexp(x))
    F1("log1pmx", xsf::log1pmx(x))
    F2("xlogy", xsf::xlogy(x, y))
    F2("xlog1py", xsf::xlog1py(x, y))
    // Box-Cox transformations and convex analysis
    F2("boxcox", xsf::boxcox(x, y))
    F2("boxcox1p", xsf::boxcox1p(x, y))
    F2("inv_boxcox", xsf::inv_boxcox(x, y))
    F2("inv_boxcox1p", xsf::inv_boxcox1p(x, y))
    F1("entr", xsf::entr(x))
    F2("rel_entr", xsf::rel_entr(x, y))
    F2("kl_div", xsf::kl_div(x, y))
    F2("huber", xsf::huber(x, y))
    F2("pseudo_huber", xsf::pseudo_huber(x, y))
    // scaled Bessel functions, derivatives, and ratios
    F1("bessel_i0e", xsf::cyl_bessel_i0e(x))
    F1("bessel_i1e", xsf::cyl_bessel_i1e(x))
    F1("bessel_k0e", xsf::cyl_bessel_k0e(x))
    F1("bessel_k1e", xsf::cyl_bessel_k1e(x))
    F2("bessel_je", xsf::cyl_bessel_je(x, y))
    F2("bessel_ye", xsf::cyl_bessel_ye(x, y))
    F2("bessel_ie", xsf::cyl_bessel_ie(x, y))
    F2("bessel_ke", xsf::cyl_bessel_ke(x, y))
    F2("iv_ratio", xsf::iv_ratio(x, y))
    F2("iv_ratio_c", xsf::iv_ratio_c(x, y))
    F2("iv_ratioinv", xsf::iv_ratioinv(x, y))
    // spherical Bessel functions
    F2("sph_bessel_j", xsf::sph_bessel_j(l(x), y))
    F2("sph_bessel_y", xsf::sph_bessel_y(l(x), y))
    F2("sph_bessel_i", xsf::sph_bessel_i(l(x), y))
    F2("sph_bessel_k", xsf::sph_bessel_k(l(x), y))
    F2("sph_bessel_j_prime", xsf::sph_bessel_j_jac(l(x), y))
    F2("sph_bessel_y_prime", xsf::sph_bessel_y_jac(l(x), y))
    F2("sph_bessel_i_prime", xsf::sph_bessel_i_jac(l(x), y))
    F2("sph_bessel_k_prime", xsf::sph_bessel_k_jac(l(x), y))
    // integrals of Bessel and Struve functions
    F1("it1j0y0.0", OUT2(0, xsf::it1j0y0(x, o0, o1)))
    F1("it1j0y0.1", OUT2(1, xsf::it1j0y0(x, o0, o1)))
    F1("it2j0y0.0", OUT2(0, xsf::it2j0y0(x, o0, o1)))
    F1("it2j0y0.1", OUT2(1, xsf::it2j0y0(x, o0, o1)))
    F1("it1i0k0.0", OUT2(0, xsf::it1i0k0(x, o0, o1)))
    F1("it1i0k0.1", OUT2(1, xsf::it1i0k0(x, o0, o1)))
    F1("it2i0k0.0", OUT2(0, xsf::it2i0k0(x, o0, o1)))
    F1("it2i0k0.1", OUT2(1, xsf::it2i0k0(x, o0, o1)))
    F3("besselpoly", xsf::besselpoly(x, y, z))
    F1("itstruve0", xsf::itstruve0(x))
    F1("it2struve0", xsf::it2struve0(x))
    F1("itmodstruve0", xsf::itmodstruve0(x))
    // Kelvin functions
    F1("ber", xsf::ber(x))
    F1("bei", xsf::bei(x))
    F1("ker", xsf::ker(x))
    F1("kei", xsf::kei(x))
    F1("berp", xsf::berp(x))
    F1("beip", xsf::beip(x))
    F1("kerp", xsf::kerp(x))
    F1("keip", xsf::keip(x))
    // Airy functions
    F1("airy[0]", OUT4(0, xsf::airy(x, o0, o1, o2, o3)))
    F1("airy[1]", OUT4(1, xsf::airy(x, o0, o1, o2, o3)))
    F1("airy[2]", OUT4(2, xsf::airy(x, o0, o1, o2, o3)))
    F1("airy[3]", OUT4(3, xsf::airy(x, o0, o1, o2, o3)))
    F1("airy_scaled[0]", OUT4(0, xsf::airye(x, o0, o1, o2, o3)))
    F1("airy_scaled[1]", OUT4(1, xsf::airye(x, o0, o1, o2, o3)))
    F1("airy_scaled[2]", OUT4(2, xsf::airye(x, o0, o1, o2, o3)))
    F1("airy_scaled[3]", OUT4(3, xsf::airye(x, o0, o1, o2, o3)))
    F1("airy_integrals[0]", OUT4(0, xsf::itairy(x, o0, o1, o2, o3)))
    F1("airy_integrals[1]", OUT4(1, xsf::itairy(x, o0, o1, o2, o3)))
    F1("airy_integrals[2]", OUT4(2, xsf::itairy(x, o0, o1, o2, o3)))
    F1("airy_integrals[3]", OUT4(3, xsf::itairy(x, o0, o1, o2, o3)))
    // parabolic cylinder functions
    F2("pbdv.0", OUT2(0, xsf::pbdv(x, y, o0, o1)))
    F2("pbdv.1", OUT2(1, xsf::pbdv(x, y, o0, o1)))
    F2("pbvv.0", OUT2(0, xsf::pbvv(x, y, o0, o1)))
    F2("pbvv.1", OUT2(1, xsf::pbvv(x, y, o0, o1)))
    F2("pbwa.0", OUT2(0, xsf::pbwa(x, y, o0, o1)))
    F2("pbwa.1", OUT2(1, xsf::pbwa(x, y, o0, o1)))
    // zeta functions
    F2("gen_harmonic", xsf::gen_harmonic(l(x), y))
    F4("normalized_gen_harmonic", xsf::normalized_gen_harmonic(l(x), l(y), l(z), w))
    // exponential, sine, cosine, and Fresnel integrals
    F1("exp1", xsf::exp1(x))
    F1("scaled_exp1", xsf::scaled_exp1(x))
    F1("sici.0", OUT2(0, xsf::sici(x, o0, o1)))
    F1("sici.1", OUT2(1, xsf::sici(x, o0, o1)))
    F1("shichi.0", OUT2(0, xsf::shichi(x, o0, o1)))
    F1("shichi.1", OUT2(1, xsf::shichi(x, o0, o1)))
    F1("fresnel.0", OUT2(0, xsf::fresnel(x, o0, o1)))
    F1("fresnel.1", OUT2(1, xsf::fresnel(x, o0, o1)))
    // inverse incomplete gamma and beta functions
    F2("gammaincinv", xsf::gammaincinv(x, y))
    F2("gammainccinv", xsf::gammainccinv(x, y))
    F3("betaincinv", xsf::cephes::incbi(x, y, z))
    // elliptic functions and integrals
    F1("ellipkm1", xsf::ellipkm1(x))
    F2("elliprc", xsf::cpu::elliprc(x, y))
    F3("elliprd", xsf::cpu::elliprd(x, y, z))
    F3("elliprf", xsf::cpu::elliprf(x, y, z))
    F3("elliprg", xsf::cpu::elliprg(x, y, z))
    F4("elliprj", xsf::cpu::elliprj(x, y, z, w))
    F2("ellipj.0", OUT4(0, xsf::ellipj(x, y, o0, o1, o2, o3)))
    F2("ellipj.1", OUT4(1, xsf::ellipj(x, y, o0, o1, o2, o3)))
    F2("ellipj.2", OUT4(2, xsf::ellipj(x, y, o0, o1, o2, o3)))
    F2("ellipj.3", OUT4(3, xsf::ellipj(x, y, o0, o1, o2, o3)))
    // hypergeometric functions
    F2("hyp0f1", xsf::hyp0f1(x, y))
    F3("hyp1f1", xsf::hyp1f1(x, y, z))
    F3("hypu", xsf::hyperu(x, y, z))
    F4("hyp2f1", xsf::hyp2f1(x, y, z, w))
    // orthogonal polynomials; the `_l` variants have an integer degree
    F4("eval_jacobi", xsf::eval_jacobi(x, y, z, w))
    F4("eval_jacobi_l", xsf::eval_jacobi(l(x), y, z, w))
    F4("eval_gengegenbauer", xsf::eval_sh_jacobi(x, y, z, w))
    F4("eval_gengegenbauer_l", xsf::eval_sh_jacobi(l(x), y, z, w))
    F3("eval_gegenbauer", xsf::eval_gegenbauer(x, y, z))
    F3("eval_gegenbauer_l", xsf::eval_gegenbauer(l(x), y, z))
    F2("eval_chebyshev_t", xsf::eval_chebyt(x, y))
    F2("eval_chebyshev_t_l", xsf::eval_chebyt(l(x), y))
    F2("eval_chebyshev_u", xsf::eval_chebyu(x, y))
    F2("eval_chebyshev_u_l", xsf::eval_chebyu(l(x), y))
    F2("eval_chebyshev_c", xsf::eval_chebyc(x, y))
    F2("eval_chebyshev_c_l", xsf::eval_chebyc(l(x), y))
    F2("eval_chebyshev_s", xsf::eval_chebys(x, y))
    F2("eval_chebyshev_s_l", xsf::eval_chebys(l(x), y))
    F2("eval_sh_chebyshev_t", xsf::eval_sh_chebyt(x, y))
    F2("eval_sh_chebyshev_t_l", xsf::eval_sh_chebyt(l(x), y))
    F2("eval_sh_chebyshev_u", xsf::eval_sh_chebyu(x, y))
    F2("eval_sh_chebyshev_u_l", xsf::eval_sh_chebyu(l(x), y))
    F2("eval_legendre", xsf::eval_legendre(x, y))
    F2("eval_legendre_l", xsf::eval_legendre(l(x), y))
    F2("eval_sh_legendre", xsf::eval_sh_legendre(x, y))
    F2("eval_sh_legendre_l", xsf::eval_sh_legendre(l(x), y))
    F3("eval_genlaguerre", xsf::eval_genlaguerre(x, y, z))
    F3("eval_genlaguerre_l", xsf::eval_genlaguerre(l(x), y, z))
    F2("eval_laguerre", xsf::eval_laguerre(x, y))
    F2("eval_laguerre_l", xsf::eval_laguerre(l(x), y))
    F2("eval_hermite_h_l", xsf::eval_hermite(l(x), y))
    F2("eval_hermite_he_l", xsf::eval_hermitenorm(l(x), y))
    // Legendre functions
    F2("legendre_p", xsf::legendre_p(i(x), y))
    F3("assoc_legendre_p", xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, i(x), i(y), z, 2))
    F3("assoc_legendre_p_norm", xsf::assoc_legendre_p(xsf::assoc_legendre_norm, i(x), i(y), z, 2))
    F3("sph_legendre_p", xsf::sph_legendre_p(i(x), i(y), z))
    F3("pmv", xsf::pmv(x, y, z))
    // discrete distributions
    F3("bdtr", xsf::bdtr(x, i(y), z))
    F3("bdtrc", xsf::bdtrc(x, i(y), z))
    F3("bdtri", xsf::bdtri(x, i(y), z))
    F3("nbdtr", xsf::nbdtr(i(x), i(y), z))
    F3("nbdtrc", xsf::nbdtrc(i(x), i(y), z))
    F3("nbdtri", xsf::nbdtri(i(x), i(y), z))
    F2("pdtr", xsf::pdtr(x, y))
    F2("pdtrc", xsf::pdtrc(x, y))
    F2("pdtri", xsf::pdtri(i(x), y))
    // continuous distributions
    F3("gdtr", xsf::gdtr(x, y, z))
    F3("gdtrc", xsf::gdtrc(x, y, z))
    F3("gdtria", xsf::gdtria(x, y, z))
    F3("gdtrib", xsf::gdtrib(x, y, z))
    F3("gdtrix", xsf::gdtrix(x, y, z))
    F2("chdtr", xsf::chdtr(x, y))
    F2("chdtrc", xsf::chdtrc(x, y))
    F2("chdtri", xsf::chdtri(x, y))
    F3("fdtri", xsf::fdtri(x, y, z))
    F2("tukeylambdacdf", xsf::tukeylambdacdf(x, y))
    F2("von_mises_cdf", xsf::von_mises_cdf(x, y))
    F3("geninvgauss_pdf", xsf::cpu::geninvgauss_pdf(x, y, z))
    F3("geninvgauss_logpdf", xsf::cpu::geninvgauss_logpdf(x, y, z))
    F4("genhyperbolic_pdf", xsf::cpu::genhyperbolic_pdf(x, y, z, w))
    F4("genhyperbolic_logpdf", xsf::cpu::genhyperbolic_logpdf(x, y, z, w))
    F1("cdf_cvm_inf", xsf::cpu::cdf_cvm_inf(x))
    // Kolmogorov-Smirnov distributions
    F2("smirnov", xsf::cpu::smirnov(i(x), y))
    F2("smirnovc", xsf::cpu::smirnovc(i(x), y))
    F2("smirnovp", xsf::cpu::smirnovp(i(x), y))
    F2("smirnovi", xsf::cpu::smirnovi(i(x), y))
    F2("smirnovci", xsf::cpu::smirnovci(i(x), y))
    F1("kolmogorov", xsf::cpu::kolmogorov(x))
    F1("kolmogc", xsf::cpu::kolmogc(x))
    F1("kolmogp", xsf::cpu::kolmogp(x))
    F1("kolmogi", xsf::cpu::kolmogi(x))
    F1("kolmogci", xsf::cpu::kolmogci(x))
    // Wright's generalized Bessel function
    F3("wright_bessel", xsf::wright_bessel(x, y, z))
    F3("log_wright_bessel", xsf::log_wright_bessel(x, y, z))
    // complex functions
    C("gamma", 2, xsf::gamma(cd(a[0], a[1])))
    C("loggamma", 2, xsf::loggamma(cd(a[0], a[1])))
    C("rgamma", 2, xsf::rgamma(cd(a[0], a[1])))
    C("digamma", 2, xsf::digamma(cd(a[0], a[1])))
    C("erf", 2, xsf::erf(cd(a[0], a[1])))
    C("erfc", 2, xsf::erfc(cd(a[0], a[1])))
    C("wofz", 2, xsf::wofz(cd(a[0], a[1])))
    C("exp1", 2, xsf::exp1(cd(a[0], a[1])))
    C("expi", 2, xsf::expi(cd(a[0], a[1])))
    C("lambertw", 2, xsf::lambertw(cd(a[0], a[1]), 0, 1e-8))
    C("hyp1f1", 4, xsf::hyp1f1(a[0], a[1], cd(a[2], a[3])))
    C("hyp2f1", 5, xsf::hyp2f1(a[0], a[1], a[2], cd(a[3], a[4])))
    C("bessel_j", 3, xsf::cyl_bessel_j(a[0], cd(a[1], a[2])))
    C("bessel_y", 3, xsf::cyl_bessel_y(a[0], cd(a[1], a[2])))
    C("bessel_i", 3, xsf::cyl_bessel_i(a[0], cd(a[1], a[2])))
    C("bessel_k", 3, xsf::cyl_bessel_k(a[0], cd(a[1], a[2])))
    C("hankel_1", 3, xsf::cyl_hankel_1(a[0], cd(a[1], a[2])))
    C("hankel_2", 3, xsf::cyl_hankel_2(a[0], cd(a[1], a[2])))
    C("hankel_1e", 3, xsf::cyl_hankel_1e(a[0], cd(a[1], a[2])))
    C("hankel_2e", 3, xsf::cyl_hankel_2e(a[0], cd(a[1], a[2])))
    C("modified_fresnel_plus.0", 1, nthc<0>([&](cd &o0, cd &o1) { xsf::modified_fresnel_plus(a[0], o0, o1); }))
    C("modified_fresnel_plus.1", 1, nthc<1>([&](cd &o0, cd &o1) { xsf::modified_fresnel_plus(a[0], o0, o1); }))
    C("modified_fresnel_minus.0", 1, nthc<0>([&](cd &o0, cd &o1) { xsf::modified_fresnel_minus(a[0], o0, o1); }))
    C("modified_fresnel_minus.1", 1, nthc<1>([&](cd &o0, cd &o1) { xsf::modified_fresnel_minus(a[0], o0, o1); }))
    C("sph_harm_y", 4, xsf::sph_harm_y(i(a[0]), i(a[1]), a[2], a[3]))
#undef F1
#undef F2
#undef F3
#undef F4
#undef C
#undef OUT2
#undef OUT4
    return false;
}

int main() {
    std::string line, f;
    std::vector<double> a;
    while (std::getline(std::cin, line)) {
        std::istringstream ss(line);
        ss >> f;
        a.clear();
        std::string tok;
        while (ss >> tok) {
            a.push_back(std::strtod(tok.c_str(), nullptr));
        }
        double r;
        cd rc;
        bool cplx = false;
        if (!eval(f, a, r, rc, cplx)) {
            std::fprintf(stderr, "unknown function: %s\n", line.c_str());
            return 1;
        }
        if (cplx) {
            std::printf("%.17g %.17g\n", rc.real(), rc.imag());
        } else {
            std::printf("%.17g\n", r);
        }
    }
    return 0;
}
