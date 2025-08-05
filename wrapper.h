#pragma once

#ifdef __cplusplus
extern "C" {
#endif

float gamma_f(float x);
double gamma_d(double x);

float gammaln_f(float x);
double gammaln_d(double x);

float gammasgn_f(float x);
double gammasgn_d(double x);

float gammainc_f(float a, float x);
double gammainc_d(double a, double x);

float gammaincinv_f(float a, float p);
double gammaincinv_d(double a, double p);

float gammaincc_f(float a, float x);
double gammaincc_d(double a, double x);

float gammainccinv_f(float a, float p);
double gammainccinv_d(double a, double p);

float gamma_ratio_f(float a, float b);
double gamma_ratio_d(double a, double b);

#ifdef __cplusplus
} // extern "C"
#endif
