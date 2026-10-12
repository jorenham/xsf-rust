//! Translated into pure Rust from `xsf/cephes/kolmogorov.h` (xsf v0.2.2), which was translated
//! into C++ by SciPy developers in 2024. Original header comment appears below.
//!
//! File altered for inclusion in cephes module for Python:
//! Main loop commented out....
//! Travis Oliphant Nov. 1998
//!
//! Re Kolmogorov statistics, here is Birnbaum and Tingey's (actually it was already present
//! in Smirnov's paper) formula for the
//! distribution of D+, the maximum of all positive deviations between a
//! theoretical distribution function P(x) and an empirical one Sn(x)
//! from n samples.
//!
//! ```text
//!     +
//!    D  =         sup     [P(x) - S (x)]
//!     n     -inf < x < inf         n
//!
//!
//!                  [n(1-d)]
//!        +            -                    v-1              n-v
//!    Pr{D   > d} =    >    C    d (d + v/n)    (1 - d - v/n)
//!        n            -   n v
//!                    v=0
//! ```
//!
//! (also equals the following sum, but note the terms may be large and alternating in sign)
//! See Smirnov 1944, Dwass 1959
//!
//! ```text
//!                         n
//!                         -                         v-1              n-v
//!                =  1 -   >         C    d (d + v/n)    (1 - d - v/n)
//!                         -        n v
//!                       v=[n(1-d)]+1
//! ```
//!
//! \[n(1-d)\] is the largest integer not exceeding n(1-d).
//! nCv is the number of combinations of n things taken v at a time.
//!
//! Sources:
//! \[1\] Smirnov, N.V. "Approximate laws of distribution of random variables from empirical data"
//!     Usp. Mat. Nauk, 1944. <http://mi.mathnet.ru/umn8798>
//! \[2\] Birnbaum, Z. W. and Tingey, Fred H.
//!     "One-Sided Confidence Contours for Probability Distribution Functions",
//!     Ann. Math. Statist. 1951.  <https://doi.org/10.1214/aoms/1177729550>
//! \[3\] Dwass, Meyer, "The Distribution of a Generalized $\mathrm{D}^+_n$ Statistic",
//!     Ann. Math. Statist., 1959.  <https://doi.org/10.1214/aoms/1177706085>
//! \[4\] van Mulbregt, Paul, "Computing the Cumulative Distribution Function and Quantiles of the
//!     One-sided Kolmogorov-Smirnov Statistic"
//!     <http://arxiv.org/abs/1802.06966>
//! \[5\] van Mulbregt, Paul,  "Computing the Cumulative Distribution Function and Quantiles of the
//!     limit of the Two-sided Kolmogorov-Smirnov Statistic"
//!     <https://arxiv.org/abs/1803.00426>
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use core::f64::consts::{PI, SQRT_2, TAU};

// cephes/const.h
const LOGSQRT2PI: f64 = 0.918_938_533_204_672_8; // log(sqrt(2*pi))

/* ************************************************************************ */
/* Algorithm Configuration */

const KOLMOG_MAXITER: i32 = 500;

/*
 * Kolmogorov Two-sided:
 * Switchover between the two series to compute K(x)
 *  0 <=  x <= KOLMOG_CUTOVER and
 *  KOLMOG_CUTOVER < x < infty
 */
const KOLMOG_CUTOVER: f64 = 0.82;

/* ************************************************************************ */
/* ************************************************************************ */

/* exp() of anything below this returns 0 */
const MIN_EXPABLE: i32 = -708 - 38;

/* Class to hold the CDF, SF and PDF, which are computed simultaneously */
struct ThreeProbs {
    sf: f64,
    cdf: f64,
    pdf: f64,
}

const XTOL: f64 = f64::EPSILON;
const RTOL: f64 = 2.0 * XTOL;

fn within_tol(x: f64, y: f64, atol: f64, rtol: f64) -> bool {
    let diff = (x - y).abs();
    diff <= (atol + rtol * y.abs())
}

/* ************************************************************************ */
/* Kolmogorov : Two-sided                      **************************** */
/* ************************************************************************ */

#[allow(clippy::float_cmp)]
fn kolmogorov_(x: f64) -> ThreeProbs {
    let mut p = 1.0;
    let mut d = 0.0;
    let (sf, cdf, pdf);

    if x.is_nan() {
        return ThreeProbs {
            sf: f64::NAN,
            cdf: f64::NAN,
            pdf: f64::NAN,
        };
    }
    if x <= 0.0 {
        return ThreeProbs {
            sf: 1.0,
            cdf: 0.0,
            pdf: 0.0,
        };
    }
    /* x <= 0.040611972203751713 */
    if x <= PI / f64::from(-MIN_EXPABLE * 8).sqrt() {
        return ThreeProbs {
            sf: 1.0,
            cdf: 0.0,
            pdf: 0.0,
        };
    }

    if x <= KOLMOG_CUTOVER {
        /*
         *  u = e^(-pi^2/(8x^2))
         *  w = sqrt(2pi)/x
         *  P = w*u * (1 + u^8 + u^24 + u^48 + ...)
         */
        let w = TAU.sqrt() / x;
        let logu8 = -PI * PI / (x * x); /* log(u^8) */
        let u = (logu8 / 8.0).exp();
        if u == 0.0 {
            /*
             * P = w*u, but u < 1e-308, and w > 1,
             * so compute as logs, then exponentiate
             */
            let log_p = logu8 / 8.0 + w.ln();
            p = log_p.exp();
        } else {
            /* Just unroll the loop, 3 iterations */
            let u8 = logu8.exp();
            let u8cub = u8.powf(3.0);
            p = 1.0 + u8cub * p;
            d = 5.0 * 5.0 + u8cub * d;
            p = 1.0 + u8 * u8 * p;
            d = 3.0 * 3.0 + u8 * u8 * d;
            p = 1.0 + u8 * p;
            d = 1.0 * 1.0 + u8 * d;

            d = PI * PI / 4.0 / (x * x) * d - p;
            d *= w * u / x;
            p *= w * u;
        }
        cdf = p;
        sf = 1.0 - p;
        pdf = d;
    } else {
        /*
         *  v = e^(-2x^2)
         *  P = 2 (v - v^4 + v^9 - v^16 + ...)
         *    = 2v(1 - v^3*(1 - v^5*(1 - v^7*(1 - ...)))
         */
        let logv = -2.0 * x * x;
        let v = logv.exp();
        /*
         * Want q^((2k-1)^2)(1-q^(4k-1)) / q(1-q^3) < epsilon to break out of loop.
         * With KOLMOG_CUTOVER ~ 0.82, k <= 4.  Just unroll the loop, 4 iterations
         */
        let vsq = v * v;
        let v3 = v.powf(3.0);

        let mut vpwr = v3 * v3 * v; /* v**7 */
        p = 1.0 - vpwr * p; /* P <- 1 - (1-v**(2k-1)) * P */
        d = 3.0 * 3.0 - vpwr * d;

        vpwr = v3 * vsq;
        p = 1.0 - vpwr * p;
        d = 2.0 * 2.0 - vpwr * d;

        vpwr = v3;
        p = 1.0 - vpwr * p;
        d = 1.0 * 1.0 - vpwr * d;

        p *= 2.0 * v;
        d *= 8.0 * v * x;
        sf = p;
        cdf = 1.0 - sf;
        pdf = d;
    }
    ThreeProbs {
        sf: sf.clamp(0.0, 1.0),
        cdf: cdf.clamp(0.0, 1.0),
        pdf: 0.0_f64.max(pdf),
    }
}

/* Find x such kolmogorov(x)=psf, kolmogc(x)=pcdf */
#[allow(clippy::float_cmp, clippy::manual_range_contains)]
fn kolmogi_(psf: f64, pcdf: f64) -> f64 {
    let mut iterations = 0;
    let (mut a, mut b);
    let mut x;

    if !(psf >= 0.0 && pcdf >= 0.0 && pcdf <= 1.0 && psf <= 1.0) {
        return f64::NAN;
    }
    if (1.0 - pcdf - psf).abs() > 4.0 * f64::EPSILON {
        return f64::NAN;
    }
    if pcdf == 0.0 {
        return 0.0;
    }
    if psf == 0.0 {
        return f64::INFINITY;
    }

    if pcdf <= 0.5 {
        /* p ~ (sqrt(2pi)/x) *exp(-pi^2/8x^2).  Generate lower and upper bounds  */
        let logpcdf = pcdf.ln();
        /* Now that 1 >= x >= sqrt(p) */
        /* Iterate twice: x <- pi/(sqrt(8) sqrt(log(sqrt(2pi)) - log(x) - log(pdf))) */
        a = PI / (2.0 * SQRT_2 * (-(logpcdf + logpcdf / 2.0 - LOGSQRT2PI)).sqrt());
        b = PI / (2.0 * SQRT_2 * (-(logpcdf + 0.0 - LOGSQRT2PI)).sqrt());
        a = PI / (2.0 * SQRT_2 * (-(logpcdf + a.ln() - LOGSQRT2PI)).sqrt());
        b = PI / (2.0 * SQRT_2 * (-(logpcdf + b.ln() - LOGSQRT2PI)).sqrt());
        x = f64::midpoint(a, b);
    } else {
        /*
         * Based on the approximation p ~ 2 exp(-2x^2)
         * Found that needed to replace psf with a slightly smaller number in the second element
         *  as otherwise _kolmogorov(b) came back as a very small number but with
         *  the same sign as _kolmogorov(a)
         *  kolmogi(0.5) = 0.82757355518990772
         *  so (1-q^(-(4-1)*2*x^2)) = (1-exp(-6*0.8275^2) ~ (1-exp(-4.1)
         */
        const JIGGERB: f64 = 256.0 * f64::EPSILON;
        let pba = psf / (1.0 - (-4.0_f64).exp()) / 2.0;
        let pbb = psf * (1.0 - JIGGERB) / 2.0;
        a = (-0.5 * pba.ln()).sqrt();
        b = (-0.5 * pbb.ln()).sqrt();
        /*
         * Use inversion of
         *   p = q - q^4 + q^9 - q^16 + ...:
         *   q = p + p^4 + 4p^7 - p^9 + 22p^10  - 13p^12 + 140*p^13 ...
         */
        let q0 = {
            let p = psf / 2.0;
            let p2 = p * p;
            let p3 = p * p * p;
            let q0 =
                1.0 + p3 * (1.0 + p3 * (4.0 + p2 * (-1.0 + p * (22.0 + p2 * (-13.0 + 140.0 * p)))));
            q0 * p
        };
        x = (-q0.ln() / 2.0).sqrt();
        if x < a || x > b {
            x = f64::midpoint(a, b);
        }
    }
    debug_assert!(a <= b);

    loop {
        let x0 = x;
        let probs = kolmogorov_(x0);
        let df = if pcdf < 0.5 {
            pcdf - probs.cdf
        } else {
            probs.sf - psf
        };

        if df.abs() == 0.0 {
            break;
        }
        /* Update the bracketing interval */
        if df > 0.0 && x > a {
            a = x;
        } else if df < 0.0 && x < b {
            b = x;
        }

        let dfdx = -probs.pdf;
        if dfdx.abs() <= 0.0 {
            x = f64::midpoint(a, b);
        } else {
            let t = df / dfdx;
            x = x0 - t;
        }

        /*
         * Check out-of-bounds.
         * Not expecting this to happen often --- kolmogorov is convex near x=infinity and
         * concave near x=0, and we should be approaching from the correct side.
         * If out-of-bounds, replace x with a midpoint of the bracket.
         */
        if x >= a && x <= b {
            if within_tol(x, x0, XTOL, RTOL) {
                break;
            }
            if (x == a) || (x == b) {
                x = f64::midpoint(a, b);
                /* If the bracket is already so small ... */
                if x == a || x == b {
                    break;
                }
            }
        } else {
            x = f64::midpoint(a, b);
            if within_tol(x, x0, XTOL, RTOL) {
                break;
            }
        }

        iterations += 1;
        if iterations > KOLMOG_MAXITER {
            break;
        }
    }
    x
}

pub(crate) fn kolmogorov(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    kolmogorov_(x).sf
}

pub(crate) fn kolmogc(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    kolmogorov_(x).cdf
}

pub(crate) fn kolmogp(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return -0.0;
    }
    -kolmogorov_(x).pdf
}

/// Functional inverse of Kolmogorov survival statistic for two-sided test.
/// Finds x such that kolmogorov(x) = p.
pub(crate) fn kolmogi(p: f64) -> f64 {
    if p.is_nan() {
        return f64::NAN;
    }
    kolmogi_(p, 1.0 - p)
}

/// Functional inverse of Kolmogorov cumulative statistic for two-sided test.
/// Finds x such that kolmogc(x) = p = (or kolmogorov(x) = 1-p).
pub(crate) fn kolmogci(p: f64) -> f64 {
    if p.is_nan() {
        return f64::NAN;
    }
    kolmogi_(1.0 - p, p)
}
