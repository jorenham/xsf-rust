#include "wrapper.h"
#include "xsf/gamma.h"

extern "C" {

float gamma_f(float x) { return xsf::gamma(x); }
double gamma_d(double x) { return xsf::gamma(x); }

float gammaln_f(float x) { return xsf::gammaln(x); }
double gammaln_d(double x) { return xsf::gammaln(x); }

float gammasgn_f(float x) { return xsf::gammasgn(x); }
double gammasgn_d(double x) { return xsf::gammasgn(x); }

float gammainc_f(float a, float x) { return xsf::gammainc(a, x); }
double gammainc_d(double a, double x) { return xsf::gammainc(a, x); }

float gammaincinv_f(float a, float p) { return xsf::gammaincinv(a, p); }
double gammaincinv_d(double a, double p) { return xsf::gammaincinv(a, p); }

float gammaincc_f(float a, float x) { return xsf::gammaincc(a, x); }
double gammaincc_d(double a, double x) { return xsf::gammaincc(a, x); }

float gammainccinv_f(float a, float p) { return xsf::gammainccinv(a, p); }
double gammainccinv_d(double a, double p) { return xsf::gammainccinv(a, p); }

float gamma_ratio_f(float a, float b) { return xsf::gamma_ratio(a, b); }
double gamma_ratio_d(double a, double b) { return xsf::gamma_ratio(a, b); }

} // extern "C"
