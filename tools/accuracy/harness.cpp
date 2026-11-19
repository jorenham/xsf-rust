// Evaluates xsf (C++) functions for the accuracy table and the test tables (see accuracy.py and
// tools/tables/generate.py).
//
// Usage: harness < INPUT > OUTPUT
//
// Each input line is `<function> <arg>...`; each output line is the double result (%.17g), or the
// real and imaginary parts of a complex result.
// The function names are those of the xsf-rust functions in specs.py, and must match
// src/main.rs.

#include <algorithm>
#include <cmath>
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
#include "xsf/cephes/round.h"
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
#include "xsf/cpu/mathieu_legacy.h"
#include "xsf/cpu/orthogonal_eval.h"
#include "xsf/cpu/par_cyl.h"
#include "xsf/cpu/specfun.h"
#include "xsf/cpu/sph_bessel.h"
#include "xsf/cpu/sph_harm.h"
#include "xsf/cpu/sphd_wave.h"
#include "xsf/cpu/stats.h"
#include "xsf/cpu/struve.h"
#include "xsf/digamma.h"
#include "xsf/digammainv.h"
#include "xsf/ellip.h"
#include "xsf/erf.h"
#include "xsf/evalpoly.h"
#include "xsf/expint.h"
#include "xsf/gamma.h"
#include "xsf/gen_harmonic.h"
#include "xsf/hyp2f1.h"
#include "xsf/iv_ratio.h"
#include "xsf/lambertw.h"
#include "xsf/log.h"
#include "xsf/log_exp.h"
#include "xsf/loggamma.h"
#include "xsf/mdspan.h"
#include "xsf/multivariate_normal.h"
#include "xsf/ndtri_exp.h"
#include "xsf/orthogonal_eval.h"
#include "xsf/sici.h"
#include "xsf/spence.h"
#include "xsf/stats.h"
#include "xsf/stirling2.h"
#include "xsf/trig.h"
#include "xsf/trigamma.h"
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

template <typename T> using vec = xsf::cxx::mdspan<T, xsf::cxx::dextents<int, 1>>;
template <typename T> using mat = xsf::cxx::mdspan<T, xsf::cxx::dextents<int, 2>>;

// the number of orders of xsf::rctj and xsf::rcty; must match src/main.rs
constexpr int RICCATI_N = 101;

// order `n` of the values (I = 0) or derivatives (I = 1) of `f(x, values, derivatives)`, with the
// orders 0, ..., RICCATI_N - 1
template <int I, typename F> double riccati(F f, int n, double x) {
    std::vector<double> o(2 * RICCATI_N, NAN);
    f(x, vec<double>(o.data(), RICCATI_N), vec<double>(o.data() + RICCATI_N, RICCATI_N));
    return o[I * RICCATI_N + n];
}

// the element of degree `n` and order `m` of the output of `f(res)` with the degrees 0, ..., n and
// the orders 0, ..., |m|, followed by -|m|, ..., -1
template <typename T, typename F> T all_nm(F f, int n, int m) {
    const int cols = 2 * std::abs(m) + 1;
    std::vector<T> o(std::size_t(n + 1) * cols, NAN);
    f(mat<T>(o.data(), n + 1, cols));
    return o[n * cols + (m < 0 ? m + cols : m)];
}

// degree `n` of xsf::legendre_p_all with the degrees 0, ..., n
template <typename T> T legendre_p_all(int n, T x) {
    std::vector<T> o(n + 1, NAN);
    xsf::legendre_p_all(x, vec<T>(o.data(), n + 1));
    return o[n];
}

// degree `n` of the values (I = 0) or derivatives (I = 1) of xsf::lqn with the degrees 0, ..., n
// (at least 1, like in SciPy)
template <int I, typename T> T lqn(int n, T x) {
    const int len = std::max(n, 1) + 1;
    std::vector<T> o(2 * len, NAN);
    xsf::lqn(x, vec<T>(o.data(), len), vec<T>(o.data() + len, len));
    return o[I * len + n];
}

// the element of degree `n` and order `m` of the values (I = 0) or derivatives (I = 1) of xsf::lqmn
// with the degrees 0, ..., n and orders 0, ..., m (at least 1, like in SciPy)
template <int I, typename T> T lqmn(int n, int m, T x) {
    const int rows = std::max(m, 1) + 1, cols = std::max(n, 1) + 1;
    std::vector<T> o(std::size_t(2) * rows * cols, NAN);
    xsf::lqmn(x, mat<T>(o.data(), rows, cols), mat<T>(o.data() + rows * cols, rows, cols));
    return o[I * rows * cols + m * cols + n];
}

// element `k` of the output of xsf::poisson_binom_pmf_all (I = 0) or xsf::poisson_binom_cdf_all
// (I = 1) for the probabilities `p`
template <int I> double poisson_binom_all(vec<const double> p, int k) {
    std::vector<double> o(p.extent(0) + 1, NAN);
    const vec<double> r(o.data(), p.extent(0) + 1);
    if constexpr (I == 0) {
        xsf::poisson_binom_pmf_all(p, r);
    } else {
        xsf::poisson_binom_cdf_all(p, r);
    }
    return o[k];
}

// `f` of the output of xsf::wilcoxon_pmf_all (I = 0) or xsf::wilcoxon_cdf_all (I = 1) for `n`
template <int I, typename F> double wilcoxon_all(int n, F f) {
    std::vector<double> o(n * (n + 1) / 2 + 1, NAN);
    const vec<double> r(o.data(), n * (n + 1) / 2 + 1);
    if constexpr (I == 0) {
        xsf::wilcoxon_pmf_all(n, r);
    } else {
        xsf::wilcoxon_cdf_all(n, r);
    }
    return f(r);
}

// order `v` of the values (I = 0) or derivatives (I = 1) of xsf::detail::pbdv or pbvv as `f`, with
// the orders v - floor(v), ..., v (and v + 1 if v < 1)
template <int I, typename F> double pb_seq(F f, double v, double x) {
    const int num = std::abs(static_cast<int>(v)) + 3; // like xsf::pbdv and xsf::pbvv
    std::vector<double> o(2 * num, NAN);
    double pf, pd;
    f(x, v, o.data(), o.data() + num, &pf, &pd);
    return o[I * num + static_cast<int>(v)];
}

// order `n` of the values (I = 0) or derivatives (I = 1) of xsf::specfun::cpbdn, with the orders
// 0, ..., n (at least 1 if n >= 0, like in SciPy)
template <int I> cd pbdn_seq(int n, cd z) {
    const int m = n < 0 ? n : std::max(n, 1), len = std::abs(m) + 2;
    std::vector<cd> o(2 * len, NAN);
    xsf::specfun::cpbdn(m, z, o.data(), o.data() + len);
    return o[I * len + std::abs(n)];
}

// coefficient `k` of scipy.special.mathieu_even_coef (or mathieu_odd_coef, if `odd`), which is
// Python code around xsf::specfun::fcoef
static double mathieu_coef(bool odd, int m, double q, int k) {
    // one more than in SciPy, as fcoef reads fc[km] (scipy/xsf#312)
    double fc[252];
    const int kd = odd ? (m % 2 ? 3 : 4) : (m % 2 ? 2 : 1);
    const double a = odd ? xsf::sem_cva<double>(m, q) : xsf::cem_cva<double>(m, q);
    xsf::specfun::fcoef(kd, m, q, a, fc);
    return fc[k];
}

// the last mode of the characteristic values of `segv_seq`; must match src/main.rs
constexpr int SEGV_N = 20;

// mode `n` of the characteristic values of the modes m, ..., SEGV_N, like `pro_cv_seq` (kd = 1)
// and `obl_cv_seq` (kd = -1) of SciPy
static double segv_seq(int m, int n, double c, int kd) {
    double cv;
    std::vector<double> eg(SEGV_N - m + 2, NAN);
    xsf::specfun::segv(m, SEGV_N, c, kd, &cv, eg.data());
    return eg[n - m];
}

// the number of zeros or numbers of the functions with a const `N` in xsf-rust; must match
// src/main.rs
constexpr int ZEROS_N = 500;

// output I at zero `k` (from 1) of xsf::airyzo for Ai (kf = 1) or Bi (kf = 2)
template <int I> double airyzo(int kf, int k) {
    std::vector<double> o(4 * ZEROS_N, NAN);
    double *p = o.data();
    xsf::airyzo(ZEROS_N, kf, p, p + ZEROS_N, p + 2 * ZEROS_N, p + 3 * ZEROS_N);
    return o[I * ZEROS_N + k - 1];
}

// output I at zero `k` (from 1) of xsf::specfun::jyzo for the order `v`
template <int I> double jyzo(int v, int k) {
    std::vector<double> o(4 * ZEROS_N, NAN);
    double *p = o.data();
    xsf::specfun::jyzo(v, ZEROS_N, p, p + ZEROS_N, p + 2 * ZEROS_N, p + 3 * ZEROS_N);
    return o[I * ZEROS_N + k - 1];
}

// output I (x, n, m, or p) of zero `k` (from 1) of xsf::specfun::jdzo for `nt` zeros, with the
// array sizes of SciPy
template <int I> double jdzo(int nt, int k) {
    std::vector<double> zo(1401);
    std::vector<int> n(1400), m(1400), p(1400);
    xsf::specfun::jdzo(nt, zo.data(), n.data(), m.data(), p.data());
    // J_0' = -J_1, so their shared zeros come in either order (by rounding); as in the
    // reference, the zero of J_1 is first
    auto tie = [&](int a, int b) {
        return b < nt && n[a] == 0 && p[a] == 1 && n[b] == 1 && p[b] == 0 && m[a] == m[b];
    };
    int j = k - 1, nj = n[j], pj = p[j];
    if (tie(j, j + 1)) {
        nj = 1, pj = 0;
    } else if (j > 0 && tie(j - 1, j)) {
        nj = 0, pj = 1;
    }
    const double o[] = {zo[k], double(nj), double(m[j]), double(pj)};
    return o[I];
}

// zero `k` (from 1) of the first `k` zeros of the Kelvin function `kd` (1 to 8)
static double klvnzo(int kd, int k) {
    std::vector<double> o(k, NAN);
    xsf::klvnzo(k, kd, o.data());
    return o[k - 1];
}

// output I at zero `k` (from 1) of xsf::specfun::cyzo, for the complex (`c` = 1) or real zeros
template <int I> cd cyzo(int kf, int c, int k) {
    std::vector<cd> o(2 * ZEROS_N);
    xsf::specfun::cyzo(ZEROS_N, kf, c ? 0 : 1, o.data(), o.data() + ZEROS_N);
    return o[I * ZEROS_N + k - 1];
}

// zero `k` (from 1) of xsf::specfun::cerzo
static cd cerzo(int k) {
    std::vector<cd> o(ZEROS_N);
    xsf::specfun::cerzo(ZEROS_N, o.data());
    return o[k - 1];
}

// zero `k` (from 1) of the first `k` zeros of C (kf = 1) or S (kf = 2) of xsf::fcszo
static cd fcszo(int kf, int k) {
    std::vector<cd> o(k);
    xsf::fcszo(kf, k, o.data());
    return o[k - 1];
}

// number `n` (from 0) of `f(ZEROS_N - 1, numbers)`, i.e. xsf::specfun::bernob or eulerb
template <typename F> double numbers(F f, int n) {
    std::vector<double> o(ZEROS_N, 0.0);
    f(ZEROS_N - 1, o.data());
    return o[n];
}

// Be, Ke, Be', or Ke' (index 0 to 3) of xsf::kelvin
static cd kelvin(double x, int index) {
    cd o[4];
    xsf::kelvin(x, o[0], o[1], o[2], o[3]);
    return o[index];
}

// Ai, Ai', Bi, or Bi' (index 0 to 3) of xsf::airy, or of xsf::airye if `scaled`, at x + iy
static cd airyc(double x, double y, int index, bool scaled) {
    cd o[4];
    if (scaled) {
        xsf::airye(cd(x, y), o[0], o[1], o[2], o[3]);
    } else {
        xsf::airy(cd(x, y), o[0], o[1], o[2], o[3]);
    }
    return o[index];
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
// functions of an array: the K scalar arguments, followed by the elements `xs` of the array
#define A(name, K, expr)                                                                           \
    if (f == name && a.size() >= K) {                                                              \
        const vec<const double> xs(a.data() + K, static_cast<int>(a.size()) - K);                  \
        r = (expr);                                                                                \
        return true;                                                                               \
    }
#define F4(name, expr)                                                                             \
    if (f == name && a.size() == 4) {                                                              \
        const double x = a[0], y = a[1], z = a[2], w = a[3];                                       \
        r = (expr);                                                                                \
        return true;                                                                               \
    }
#define F5(name, expr)                                                                             \
    if (f == name && a.size() == 5) {                                                              \
        const double x = a[0], y = a[1], z = a[2], w = a[3], v = a[4];                             \
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
// the value (.0) and the derivative (.1) of a spheroidal wave function, without and with a given
// characteristic value
#define SW(name)                                                                                   \
    F4(#name "_nocv.0", OUT2(0, xsf::name##_nocv(x, y, z, w, o0, o1)))                             \
    F4(#name "_nocv.1", OUT2(1, xsf::name##_nocv(x, y, z, w, o0, o1)))                             \
    F5(#name ".0", OUT2(0, xsf::name(x, y, z, w, v, o0, o1)))                                      \
    F5(#name ".1", OUT2(1, xsf::name(x, y, z, w, v, o0, o1)))
    const auto l = [](double x) { return static_cast<long>(x); };
    const auto i = [](double x) { return static_cast<int>(x); };
    // gamma and related
    F1("gamma", xsf::gamma(x))
    F1("gammaln", xsf::gammaln(x))
    F1("rgamma", xsf::rgamma(x))
    F1("digamma", xsf::digamma(x))
    F1("digammainv", xsf::digammainv(x))
    F1("trigamma", xsf::trigamma(x))
    F2("beta", xsf::beta(x, y))
    F2("betaln", xsf::betaln(x, y))
    F2("binom", xsf::binom(x, y))
    F2("pow_rising", xsf::cephes::poch(x, y))
    F1("lanczos_sum_expg_scaled", xsf::cephes::lanczos_sum_expg_scaled(x))
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
    F1("round", xsf::cephes::round(x))
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
    // Riccati-Bessel functions
    const auto rctj = [](double x, vec<double> rj, vec<double> dj) { xsf::rctj(x, rj, dj); };
    const auto rcty = [](double x, vec<double> ry, vec<double> dy) { xsf::rcty(x, ry, dy); };
    F2("riccati_j.0", riccati<0>(rctj, i(x), y))
    F2("riccati_j.1", riccati<1>(rctj, i(x), y))
    F2("riccati_y.0", riccati<0>(rcty, i(x), y))
    F2("riccati_y.1", riccati<1>(rcty, i(x), y))
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
    F2("pbdv_seq.0", pb_seq<0>(xsf::detail::pbdv<double>, x, y))
    F2("pbdv_seq.1", pb_seq<1>(xsf::detail::pbdv<double>, x, y))
    F2("pbvv_seq.0", pb_seq<0>(xsf::detail::pbvv<double>, x, y))
    F2("pbvv_seq.1", pb_seq<1>(xsf::detail::pbvv<double>, x, y))
    F2("pbwa.0", OUT2(0, xsf::pbwa(x, y, o0, o1)))
    F2("pbwa.1", OUT2(1, xsf::pbwa(x, y, o0, o1)))
    // spheroidal wave functions
    F3("prolate_segv", xsf::prolate_segv(x, y, z))
    F3("oblate_segv", xsf::oblate_segv(x, y, z))
    F3("prolate_segv_seq", segv_seq(i(x), i(y), z, 1))
    F3("oblate_segv_seq", segv_seq(i(x), i(y), z, -1))
    SW(prolate_aswfa)
    SW(prolate_radial1)
    SW(prolate_radial2)
    SW(oblate_aswfa)
    SW(oblate_radial1)
    SW(oblate_radial2)
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
    // the last degree n (and order m) of the outputs up to n (and |m|)
    F2("legendre_p_all", legendre_p_all(i(x), y))
    F3("assoc_legendre_p_all", all_nm<double>([&](mat<double> r) {
           xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, 2, r);
       }, i(x), i(y)))
    F3("assoc_legendre_p_norm_all", all_nm<double>([&](mat<double> r) {
           xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, 2, r);
       }, i(x), i(y)))
    F3("sph_legendre_p_all", all_nm<double>([&](mat<double> r) { xsf::sph_legendre_p_all(z, r); }, i(x), i(y)))
    F2("legendre_q_all.0", lqn<0>(i(x), y))
    F2("legendre_q_all.1", lqn<1>(i(x), y))
    F3("assoc_legendre_q_all.0", lqmn<0>(i(x), i(y), z))
    F3("assoc_legendre_q_all.1", lqmn<1>(i(x), i(y), z))
    // Mathieu functions
    F2("mathieu_a", xsf::cem_cva(x, y))
    F2("mathieu_b", xsf::sem_cva(x, y))
    F3("mathieu_cem.0", OUT2(0, xsf::cem(x, y, z, o0, o1)))
    F3("mathieu_cem.1", OUT2(1, xsf::cem(x, y, z, o0, o1)))
    F3("mathieu_sem.0", OUT2(0, xsf::sem(x, y, z, o0, o1)))
    F3("mathieu_sem.1", OUT2(1, xsf::sem(x, y, z, o0, o1)))
    F3("mathieu_modcem1.0", OUT2(0, xsf::mcm1(x, y, z, o0, o1)))
    F3("mathieu_modcem1.1", OUT2(1, xsf::mcm1(x, y, z, o0, o1)))
    F3("mathieu_modcem2.0", OUT2(0, xsf::mcm2(x, y, z, o0, o1)))
    F3("mathieu_modcem2.1", OUT2(1, xsf::mcm2(x, y, z, o0, o1)))
    F3("mathieu_modsem1.0", OUT2(0, xsf::msm1(x, y, z, o0, o1)))
    F3("mathieu_modsem1.1", OUT2(1, xsf::msm1(x, y, z, o0, o1)))
    F3("mathieu_modsem2.0", OUT2(0, xsf::msm2(x, y, z, o0, o1)))
    F3("mathieu_modsem2.1", OUT2(1, xsf::msm2(x, y, z, o0, o1)))
    F3("mathieu_even_coef", mathieu_coef(false, i(x), y, i(z)))
    F3("mathieu_odd_coef", mathieu_coef(true, i(x), y, i(z)))
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
    // Poisson binomial and Wilcoxon signed-rank distributions: element `k` of the `_all` outputs
    A("poisson_binom_pmf_all", 1, poisson_binom_all<0>(xs, i(a[0])))
    A("poisson_binom_cdf_all", 1, poisson_binom_all<1>(xs, i(a[0])))
    F2("wilcoxon_pmf_all", wilcoxon_all<0>(i(x), [&](auto t) { return t(i(y)); }))
    F2("wilcoxon_cdf_all", wilcoxon_all<1>(i(x), [&](auto t) { return t(i(y)); }))
    F2("wilcoxon_pmf", wilcoxon_all<0>(i(x), [&](auto t) { return xsf::wilcoxon_pmf(y, t); }))
    F2("wilcoxon_cdf", wilcoxon_all<1>(i(x), [&](auto t) { return xsf::wilcoxon_cdf(y, t); }))
    F2("wilcoxon_sf", wilcoxon_all<1>(i(x), [&](auto t) { return xsf::wilcoxon_sf(y, t); }))
    // zeros and numbers: element `k` (from 1), or `n` (from 0)
    F1("airy_ai_zeros.0", airyzo<0>(1, i(x)))
    F1("airy_ai_zeros.1", airyzo<1>(1, i(x)))
    F1("airy_ai_zeros.2", airyzo<2>(1, i(x)))
    F1("airy_ai_zeros.3", airyzo<3>(1, i(x)))
    F1("airy_bi_zeros.0", airyzo<0>(2, i(x)))
    F1("airy_bi_zeros.1", airyzo<1>(2, i(x)))
    F1("airy_bi_zeros.2", airyzo<2>(2, i(x)))
    F1("airy_bi_zeros.3", airyzo<3>(2, i(x)))
    F2("bessel_zeros.0", jyzo<0>(i(x), i(y)))
    F2("bessel_zeros.1", jyzo<1>(i(x), i(y)))
    F2("bessel_zeros.2", jyzo<2>(i(x), i(y)))
    F2("bessel_zeros.3", jyzo<3>(i(x), i(y)))
    F2("bessel_j_jprime_zeros.0", jdzo<0>(i(x), i(y)))
    F2("bessel_j_jprime_zeros.1", jdzo<1>(i(x), i(y)))
    F2("bessel_j_jprime_zeros.2", jdzo<2>(i(x), i(y)))
    F2("bessel_j_jprime_zeros.3", jdzo<3>(i(x), i(y)))
    F1("ber_zeros", klvnzo(1, i(x)))
    F1("bei_zeros", klvnzo(2, i(x)))
    F1("ker_zeros", klvnzo(3, i(x)))
    F1("kei_zeros", klvnzo(4, i(x)))
    F1("berp_zeros", klvnzo(5, i(x)))
    F1("beip_zeros", klvnzo(6, i(x)))
    F1("kerp_zeros", klvnzo(7, i(x)))
    F1("keip_zeros", klvnzo(8, i(x)))
    F1("bernoulli", numbers(xsf::specfun::bernob, i(x)))
    F1("euler", numbers(xsf::specfun::eulerb, i(x)))
    C("bessel_y0_zeros.0", 2, cyzo<0>(0, i(a[0]), i(a[1])))
    C("bessel_y0_zeros.1", 2, cyzo<1>(0, i(a[0]), i(a[1])))
    C("bessel_y1_zeros.0", 2, cyzo<0>(1, i(a[0]), i(a[1])))
    C("bessel_y1_zeros.1", 2, cyzo<1>(1, i(a[0]), i(a[1])))
    C("bessel_y1_prime_zeros.0", 2, cyzo<0>(2, i(a[0]), i(a[1])))
    C("bessel_y1_prime_zeros.1", 2, cyzo<1>(2, i(a[0]), i(a[1])))
    C("erf_zeros", 1, cerzo(i(a[0])))
    C("fresnel_zeros.0", 1, fcszo(2, i(a[0])))
    C("fresnel_zeros.1", 1, fcszo(1, i(a[0])))
    // complex functions
    C("gamma", 2, xsf::gamma(cd(a[0], a[1])))
    C("loggamma", 2, xsf::loggamma(cd(a[0], a[1])))
    C("rgamma", 2, xsf::rgamma(cd(a[0], a[1])))
    C("digamma", 2, xsf::digamma(cd(a[0], a[1])))
    C("trigamma", 2, xsf::trigamma(cd(a[0], a[1])))
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
    C("bessel_je", 3, xsf::cyl_bessel_je(a[0], cd(a[1], a[2])))
    C("bessel_ye", 3, xsf::cyl_bessel_ye(a[0], cd(a[1], a[2])))
    C("bessel_ie", 3, xsf::cyl_bessel_ie(a[0], cd(a[1], a[2])))
    C("bessel_ke", 3, xsf::cyl_bessel_ke(a[0], cd(a[1], a[2])))
    C("airy[0]", 2, airyc(a[0], a[1], 0, false))
    C("airy[1]", 2, airyc(a[0], a[1], 1, false))
    C("airy[2]", 2, airyc(a[0], a[1], 2, false))
    C("airy[3]", 2, airyc(a[0], a[1], 3, false))
    C("airy_scaled[0]", 2, airyc(a[0], a[1], 0, true))
    C("airy_scaled[1]", 2, airyc(a[0], a[1], 1, true))
    C("airy_scaled[2]", 2, airyc(a[0], a[1], 2, true))
    C("airy_scaled[3]", 2, airyc(a[0], a[1], 3, true))
    C("erfcx", 2, xsf::erfcx(cd(a[0], a[1])))
    C("erfi", 2, xsf::erfi(cd(a[0], a[1])))
    C("dawsn", 2, xsf::dawsn(cd(a[0], a[1])))
    C("ndtr", 2, xsf::ndtr(cd(a[0], a[1])))
    C("log_ndtr", 2, xsf::log_ndtr(cd(a[0], a[1])))
    C("fresnel.0", 2, nthc<0>([&](cd &o0, cd &o1) { xsf::fresnel(cd(a[0], a[1]), o0, o1); }))
    C("fresnel.1", 2, nthc<1>([&](cd &o0, cd &o1) { xsf::fresnel(cd(a[0], a[1]), o0, o1); }))
    C("hyp0f1", 3, xsf::hyp0f1(a[0], cd(a[1], a[2])))
    C("riemann_zeta", 2, xsf::riemann_zeta(cd(a[0], a[1])))
    C("sph_bessel_j", 3, xsf::sph_bessel_j(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_y", 3, xsf::sph_bessel_y(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_i", 3, xsf::sph_bessel_i(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_k", 3, xsf::sph_bessel_k(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_j_prime", 3, xsf::sph_bessel_j_jac(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_y_prime", 3, xsf::sph_bessel_y_jac(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_i_prime", 3, xsf::sph_bessel_i_jac(l(a[0]), cd(a[1], a[2])))
    C("sph_bessel_k_prime", 3, xsf::sph_bessel_k_jac(l(a[0]), cd(a[1], a[2])))
    C("eval_jacobi", 5, xsf::eval_jacobi(a[0], a[1], a[2], cd(a[3], a[4])))
    C("eval_gengegenbauer", 5, xsf::eval_sh_jacobi(a[0], a[1], a[2], cd(a[3], a[4])))
    C("eval_gegenbauer", 4, xsf::eval_gegenbauer(a[0], a[1], cd(a[2], a[3])))
    C("eval_chebyshev_t", 3, xsf::eval_chebyt(a[0], cd(a[1], a[2])))
    C("eval_chebyshev_u", 3, xsf::eval_chebyu(a[0], cd(a[1], a[2])))
    C("eval_chebyshev_c", 3, xsf::eval_chebyc(a[0], cd(a[1], a[2])))
    C("eval_chebyshev_s", 3, xsf::eval_chebys(a[0], cd(a[1], a[2])))
    C("eval_sh_chebyshev_t", 3, xsf::eval_sh_chebyt(a[0], cd(a[1], a[2])))
    C("eval_sh_chebyshev_u", 3, xsf::eval_sh_chebyu(a[0], cd(a[1], a[2])))
    C("eval_legendre", 3, xsf::eval_legendre(a[0], cd(a[1], a[2])))
    C("eval_sh_legendre", 3, xsf::eval_sh_legendre(a[0], cd(a[1], a[2])))
    C("eval_genlaguerre", 4, xsf::eval_genlaguerre(a[0], a[1], cd(a[2], a[3])))
    C("eval_laguerre", 3, xsf::eval_laguerre(a[0], cd(a[1], a[2])))
    C("legendre_p", 3, xsf::legendre_p(i(a[0]), cd(a[1], a[2])))
    C("assoc_legendre_p", 4, xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, i(a[0]), i(a[1]), cd(a[2], a[3]), 2))
    C("assoc_legendre_p_norm", 4, xsf::assoc_legendre_p(xsf::assoc_legendre_norm, i(a[0]), i(a[1]), cd(a[2], a[3]), 2))
    C("sph_legendre_p", 4, xsf::sph_legendre_p(i(a[0]), i(a[1]), cd(a[2], a[3])))
    C("legendre_p_all", 3, legendre_p_all(i(a[0]), cd(a[1], a[2])))
    C("assoc_legendre_p_all", 4, all_nm<cd>([&](mat<cd> r) {
          xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, cd(a[2], a[3]), 2, r);
      }, i(a[0]), i(a[1])))
    C("assoc_legendre_p_norm_all", 4, all_nm<cd>([&](mat<cd> r) {
          xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, cd(a[2], a[3]), 2, r);
      }, i(a[0]), i(a[1])))
    C("sph_legendre_p_all", 4, all_nm<cd>([&](mat<cd> r) { xsf::sph_legendre_p_all(cd(a[2], a[3]), r); }, i(a[0]), i(a[1])))
    C("legendre_q_all.0", 3, lqn<0>(i(a[0]), cd(a[1], a[2])))
    C("legendre_q_all.1", 3, lqn<1>(i(a[0]), cd(a[1], a[2])))
    C("assoc_legendre_q_all.0", 4, lqmn<0>(i(a[0]), i(a[1]), cd(a[2], a[3])))
    C("assoc_legendre_q_all.1", 4, lqmn<1>(i(a[0]), i(a[1]), cd(a[2], a[3])))
    C("sici.0", 2, nthc<0>([&](cd &o0, cd &o1) { xsf::sici(cd(a[0], a[1]), o0, o1); }))
    C("sici.1", 2, nthc<1>([&](cd &o0, cd &o1) { xsf::sici(cd(a[0], a[1]), o0, o1); }))
    C("shichi.0", 2, nthc<0>([&](cd &o0, cd &o1) { xsf::shichi(cd(a[0], a[1]), o0, o1); }))
    C("shichi.1", 2, nthc<1>([&](cd &o0, cd &o1) { xsf::shichi(cd(a[0], a[1]), o0, o1); }))
    C("spence", 2, xsf::spence(cd(a[0], a[1])))
    C("wrightomega", 2, xsf::wrightomega(cd(a[0], a[1])))
    C("elliprc", 4, xsf::cpu::elliprc(cd(a[0], a[1]), cd(a[2], a[3])))
    C("elliprd", 6, xsf::cpu::elliprd(cd(a[0], a[1]), cd(a[2], a[3]), cd(a[4], a[5])))
    C("elliprf", 6, xsf::cpu::elliprf(cd(a[0], a[1]), cd(a[2], a[3]), cd(a[4], a[5])))
    C("elliprg", 6, xsf::cpu::elliprg(cd(a[0], a[1]), cd(a[2], a[3]), cd(a[4], a[5])))
    C("elliprj", 8, xsf::cpu::elliprj(cd(a[0], a[1]), cd(a[2], a[3]), cd(a[4], a[5]), cd(a[6], a[7])))
    C("sinpi", 2, xsf::sinpi(cd(a[0], a[1])))
    C("cospi", 2, xsf::cospi(cd(a[0], a[1])))
    C("log1p", 2, xsf::log1p(cd(a[0], a[1])))
    C("expm1", 2, xsf::expm1(cd(a[0], a[1])))
    C("xlogy", 4, xsf::xlogy(cd(a[0], a[1]), cd(a[2], a[3])))
    C("xlog1py", 4, xsf::xlog1py(cd(a[0], a[1]), cd(a[2], a[3])))
    C("modified_fresnel_plus.0", 1, nthc<0>([&](cd &o0, cd &o1) { xsf::modified_fresnel_plus(a[0], o0, o1); }))
    C("modified_fresnel_plus.1", 1, nthc<1>([&](cd &o0, cd &o1) { xsf::modified_fresnel_plus(a[0], o0, o1); }))
    C("modified_fresnel_minus.0", 1, nthc<0>([&](cd &o0, cd &o1) { xsf::modified_fresnel_minus(a[0], o0, o1); }))
    C("modified_fresnel_minus.1", 1, nthc<1>([&](cd &o0, cd &o1) { xsf::modified_fresnel_minus(a[0], o0, o1); }))
    C("sph_harm_y", 4, xsf::sph_harm_y(i(a[0]), i(a[1]), a[2], a[3]))
    C("pbdn_seq.0", 3, pbdn_seq<0>(i(a[0]), cd(a[1], a[2])))
    C("pbdn_seq.1", 3, pbdn_seq<1>(i(a[0]), cd(a[1], a[2])))
    C("sph_harm_y_all", 4, all_nm<cd>([&](mat<cd> r) { xsf::sph_harm_y_all(a[2], a[3], r); }, i(a[0]), i(a[1])))
    // the coefficients follow z
    if (f == "cevalpoly" && a.size() >= 3) {
        rc = xsf::evalpoly(a.data() + 2, static_cast<int>(a.size()) - 3, cd(a[0], a[1]));
        cplx = true;
        return true;
    }
    // for the reference tables (tools/tables/mapping.py); the parts are above
    C("kelvin.0", 1, kelvin(a[0], 0))
    C("kelvin.1", 1, kelvin(a[0], 1))
    C("kelvin.2", 1, kelvin(a[0], 2))
    C("kelvin.3", 1, kelvin(a[0], 3))
#undef F1
#undef F2
#undef F3
#undef F4
#undef A
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
