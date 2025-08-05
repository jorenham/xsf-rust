#pragma once

#ifdef __cplusplus
extern "C" {
#endif

double xsf_gamma(double x);
double xsf_gammaln(double x);
double xsf_gammasgn(double x);
double xsf_gammainc(double a, double x);
double xsf_gammaincinv(double a, double p);
double xsf_gammaincc(double a, double x);
double xsf_gammainccinv(double a, double p);
double xsf_gamma_ratio(double a, double b);

#ifdef __cplusplus
} // extern "C"
#endif
