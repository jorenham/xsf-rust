// Evaluates xsf (C++) functions for the accuracy table (see accuracy.py).
//
// Usage: harness < INPUT > OUTPUT
//
// Each input line is `<function> <arg>...`; each output line is the double result (%.17g).
// The function names are those of the xsf-rust functions in specs.py, and must match
// src/main.rs.

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "xsf/alg.h"
#include "xsf/beta.h"
#include "xsf/binom.h"
#include "xsf/cephes/erfinv.h"
#include "xsf/cephes/expn.h"
#include "xsf/cephes/incbet.h"
#include "xsf/cephes/poch.h"
#include "xsf/cosine.h"
#include "xsf/cpu/bessel.h"
#include "xsf/cpu/exp.h"
#include "xsf/cpu/hyp0f1.h"
#include "xsf/cpu/struve.h"
#include "xsf/digamma.h"
#include "xsf/digammainv.h"
#include "xsf/ellip.h"
#include "xsf/erf.h"
#include "xsf/expint.h"
#include "xsf/gamma.h"
#include "xsf/lambertw.h"
#include "xsf/log.h"
#include "xsf/loggamma.h"
#include "xsf/spence.h"
#include "xsf/stats.h"
#include "xsf/trig.h"
#include "xsf/zeta.h"

static bool eval(const std::string &f, const std::vector<double> &a, double &r) {
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
    F2("hyp0f1", xsf::hyp0f1(x, y))
    F1("lambertw", xsf::lambertw(std::complex<double>(x, 0.0), 0, 1e-8).real())
    F1("spence", xsf::spence(x))
#undef F1
#undef F2
#undef F3
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
        if (!eval(f, a, r)) {
            std::fprintf(stderr, "unknown function: %s\n", line.c_str());
            return 1;
        }
        std::printf("%.17g\n", r);
    }
    return 0;
}
