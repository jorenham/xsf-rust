#include "wrapper.h"
#include "xsf/gamma.h"

extern "C" {

double xsf_gamma(double x) { return xsf::gamma(x); }
double xsf_gammaln(double x) { return xsf::gammaln(x); }
double xsf_gammasgn(double x) { return xsf::gammasgn(x); }
double xsf_gammainc(double a, double x) { return xsf::gammainc(a, x); }
double xsf_gammaincinv(double a, double p) { return xsf::gammaincinv(a, p); }
double xsf_gammaincc(double a, double x) { return xsf::gammaincc(a, x); }
double xsf_gammainccinv(double a, double p) { return xsf::gammainccinv(a, p); }
double xsf_gamma_ratio(double a, double b) { return xsf::gamma_ratio(a, b); }

} // extern "C"
