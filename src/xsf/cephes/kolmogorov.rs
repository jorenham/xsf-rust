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
//!
//! `std::pow(6.0 * n * x + 1, 2)` is written as `t * t`, which is what GCC compiles it to, and
//! what LLVM only does when optimizing.
//!
//! Unlike in xsf (scipy/xsf#315), the `int` products that can overflow are computed in `f64`,
//! e.g. `smirnovp(50000, 0.002)` was `-0.0` instead of `-268.21...`.
//!
//! The `long double` arithmetic in `_smirnovi` (with `SCIPY_El`, `xmax` and `xmax6`) is done in
//! `f64`. This only affects the initial bracket and estimate of the root, and occasionally the
//! last bit of the result.

use core::f64::consts::{E, LN_2, PI, SQRT_2};

use crate::utils::{LN_SQRT_2PI, MINLOG, SQRT_TAU};
use crate::xsf::cephes::dd_real::DoubleDouble;

// cephes/const.h
const SCIPY_EL: f64 = E; // (a `long double` in xsf)

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

/*
 * Smirnov One-sided:
 * n larger than SMIRNOV_MAX_COMPUTE_N will result in an approximation
 */
const SMIRNOV_MAX_COMPUTE_N: i32 = 1_000_000;

/*
 * Use the upper sum formula, if the number of terms is at most SM_UPPER_MAX_TERMS,
 * and n is at least SM_UPPERSUM_MIN_N
 * Don't use the upper sum if lots of terms are involved as the series alternates
 *  sign and the terms get much bigger than 1.
 */
const SM_UPPER_MAX_TERMS: i32 = 3;
const SM_UPPERSUM_MIN_N: i32 = 10;

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
        let w = SQRT_TAU / x;
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
        a = PI / (2.0 * SQRT_2 * (-(logpcdf + logpcdf / 2.0 - LN_SQRT_2PI)).sqrt());
        b = PI / (2.0 * SQRT_2 * (-(logpcdf + 0.0 - LN_SQRT_2PI)).sqrt());
        a = PI / (2.0 * SQRT_2 * (-(logpcdf + a.ln() - LN_SQRT_2PI)).sqrt());
        b = PI / (2.0 * SQRT_2 * (-(logpcdf + b.ln() - LN_SQRT_2PI)).sqrt());
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

/* ************************************************************************ */
/* ********** Smirnov : One-sided ***************************************** */
/* ************************************************************************ */

#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
fn next_power_of2(x: f64) -> f64 {
    let q = x * f64::EPSILON; // std::ldexp(x, 1 - std::numeric_limits<double>::digits)
    let mut l = (q + x).abs();
    if l == 0.0 {
        l = x.abs();
    } else {
        let lint = l as i32;
        if f64::from(lint) == l {
            l = f64::from(lint);
        }
    }
    l
}

/// Returns `(alpha, nxfloor, nx)`
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
fn mod_nx(n: i32, x: f64) -> (f64, i32, f64) {
    /*
     * Compute floor(n*x) and remainder *exactly*.
     * If remainder is too close to 1 (E.g. (1, -std::numeric_limits<double>::epsilon()/2))
     *  round up and adjust   */
    let nx_d = f64::from(n) * DoubleDouble::new(x);
    let nxfloor_d = nx_d.floor();
    let alpha_d = nx_d - nxfloor_d;
    let mut alpha = alpha_d.hi;
    let mut nxfloor = nxfloor_d.hi as i32;
    debug_assert!(alpha >= 0.0);
    debug_assert!(alpha <= 1.0);
    if alpha == 1.0 {
        nxfloor += 1;
        alpha = 0.0;
    }
    debug_assert!(alpha < 1.0);
    (alpha, nxfloor, f64::from(nx_d))
}

/*
 * The binomial coefficient C overflows a 64 bit double, as the 11-bit
 * exponent is too small.
 * Store C as (Cman:double_double, Cexpt:int).
 *  I.e a Mantissa/significand, and an exponent.
 *  Cman lies between 0.5 and 1, and the exponent has >=32-bit.
 */
fn update_binomial(cman: &mut DoubleDouble, cexpt: &mut i32, n: i32, j: i32) {
    let rat = DoubleDouble::new(f64::from(n - j)) / (f64::from(j) + 1.0);
    let man2 = *cman * rat;
    let (man2, expt) = man2.frexp();
    debug_assert!(man2 != 0.0);
    *cexpt += expt;
    *cman = man2;
}

fn pow_d(a: DoubleDouble, m: i32) -> DoubleDouble {
    /*
     * Using dd_npwr() here would be quite time-consuming.
     * Tradeoff accuracy-time by using pow().
     */

    if m <= 0 {
        if m == 0 {
            return DoubleDouble::new(1.0);
        }
        return 1.0 / pow_d(a, -m);
    }
    if a == 0.0 {
        return DoubleDouble::new(0.0);
    }
    let ans = a.hi.powf(f64::from(m));
    let r = a.lo / a.hi;
    let mut adj = f64::from(m) * r;
    if adj.abs() > 1e-8 {
        if adj.abs() < 1e-4 {
            /* Take 1st two terms of Taylor Series for (1+r)^m */
            adj += (f64::from(m) * r) * (f64::from(m - 1) / 2.0 * r);
        } else {
            /* Take exp of scaled log */
            adj = (f64::from(m) * r.ln_1p()).exp_m1();
        }
    }
    DoubleDouble::new(ans) + ans * adj
}

fn pow2(a: f64, b: f64, m: i32) -> f64 {
    f64::from(pow_d(DoubleDouble::new(a) + b, m))
}

/*
 * Not 1024 as too big.  Want _MAX_EXPONENT < 1023-52 so as to keep both
 * elements of the double_double normalized
 */
const SM_MAX_EXPONENT: i32 = 960;

#[allow(clippy::cast_possible_truncation)]
fn pow2_scaled_d(a: DoubleDouble, m: i32) -> (DoubleDouble, i32) {
    /* Compute a^m = significand*2^expt and return as (significand, expt) */
    let mut max_expt = SM_MAX_EXPONENT;

    if m <= 0 {
        if m == 0 {
            return (DoubleDouble::new(1.0), 0);
        }
        let (ans, a_e1) = pow2_scaled_d(a, -m);
        let (ans, a_e2) = (1.0 / ans).frexp();
        let ans_e = -a_e1 + a_e2;
        return (ans, ans_e);
    }
    let (y, y_e) = a.frexp();
    if m == 1 {
        return (y, y_e);
    }
    /*
     *  y ^ maxExpt >= 2^{-960}
     *  =>  maxExpt = 960 / log2(y.x[0]) = 708 / log(y.x[0])
     *              = 665/((1-y.x[0] + y.x[0]^2/2 - ...)
     *              <= 665/(1-y.x[0])
     * Quick check to see if we might need to break up the exponentiation
     */
    if f64::from(m) * (y.hi - 1.0) / y.hi < f64::from(-SM_MAX_EXPONENT) * LN_2 {
        /* Now do it carefully, calling log() */
        let lg2y = y.hi.ln() / LN_2;
        let lg_ans = f64::from(m) * lg2y;
        if lg_ans <= f64::from(-SM_MAX_EXPONENT) {
            max_expt = (next_power_of2(f64::from(-SM_MAX_EXPONENT) / lg2y + 1.0) / 2.0) as i32;
        }
    }
    if m <= max_expt {
        let ans1 = pow_d(y, m);
        let (ans, mut ans_e) = ans1.frexp();
        ans_e += m * y_e;
        return (ans, ans_e);
    }

    let q = m / max_expt;
    let r = m % max_expt;
    /* y^m = (y^maxExpt)^q * y^r */
    let (y2r, y2r_e) = pow2_scaled_d(y, r);
    let (y2m, y2m_e) = pow2_scaled_d(y, max_expt);
    let (y2mq, mut y2mq_e) = pow2_scaled_d(y2m, q);
    let (ans, mut ans_e) = (y2r * y2mq).frexp();
    y2mq_e += y2m_e * q;
    ans_e += y2mq_e + y2r_e;
    ans_e += m * y_e;
    (ans, ans_e)
}

#[allow(non_snake_case)]
fn pow4_d(a: f64, b: f64, c: f64, d: f64, m: i32) -> DoubleDouble {
    /* Compute ((a+b)/(c+d)) ^ m */
    if m <= 0 {
        if m == 0 {
            return DoubleDouble::new(1.0);
        }
        return pow4_d(c, d, a, b, -m);
    }
    let A = DoubleDouble::new(a) + b;
    let C = DoubleDouble::new(c) + d;
    if A == 0.0 {
        return if C == 0.0 {
            DoubleDouble::NAN
        } else {
            DoubleDouble::new(0.0)
        };
    }
    if C == 0.0 {
        return if A < 0.0 {
            -DoubleDouble::INFINITY
        } else {
            DoubleDouble::INFINITY
        };
    }
    let X = A / C;
    pow_d(X, m)
}

fn pow4(a: f64, b: f64, c: f64, d: f64, m: i32) -> f64 {
    let ret = pow4_d(a, b, c, d, m);
    f64::from(ret)
}

#[allow(non_snake_case)]
fn logpow4_d(a: f64, b: f64, c: f64, d: f64, m: i32) -> DoubleDouble {
    /*
     * Compute log(((a+b)/(c+d)) ^ m)
     *    == m * log((a+b)/(c+d))
     *    == m * log( 1 + (a+b-c-d)/(c+d))
     */
    if m == 0 {
        return DoubleDouble::new(0.0);
    }
    let A = DoubleDouble::new(a) + b;
    let C = DoubleDouble::new(c) + d;
    if A == 0.0 {
        return if C == 0.0 {
            DoubleDouble::new(0.0)
        } else {
            -DoubleDouble::INFINITY
        };
    }
    if C == 0.0 {
        return DoubleDouble::INFINITY;
    }
    let X = A / C;
    debug_assert!(X.hi >= 0.0);
    let ans = if 0.5 <= X.hi && X.hi <= 1.5 {
        let A1 = A - C;
        let X1 = A1 / C;
        X1.log1p()
    } else {
        X.log()
    };
    f64::from(m) * ans
}

fn logpow4(a: f64, b: f64, c: f64, d: f64, m: i32) -> f64 {
    let ans = logpow4_d(a, b, c, d, m);
    f64::from(ans)
}

/*
 *  Compute a single term in the summation, A_v(n, x):
 *  A_v(n, x) =  Binomial(n,v) * (1-x-v/n)^(n-v) * (x+v/n)^(v-1)
 */
fn compute_av(n: i32, x: f64, v: i32, cman: DoubleDouble, cexpt: i32) -> DoubleDouble {
    let t2x = DoubleDouble::new(f64::from(n - v)) / f64::from(n) - x; /*  1 - x - v/n */
    let (t2, t2_e) = pow2_scaled_d(t2x, n - v);
    let t1x = DoubleDouble::new(f64::from(v)) / f64::from(n) + x; /* x + v/n */
    let (t1, t1_e) = pow2_scaled_d(t1x, v - 1);
    let mut ans = t1 * t2;
    ans = ans * cman;
    let ans_e = cexpt + t1_e + t2_e;
    ans.ldexp(ans_e)
}

#[allow(
    clippy::bool_to_int_with_if,
    clippy::float_cmp,
    clippy::manual_range_contains,
    clippy::too_many_lines
)]
fn smirnov_(n: i32, x: f64) -> ThreeProbs {
    let mut aj_sum = DoubleDouble::new(0.0);
    let mut daj_sum = DoubleDouble::new(0.0);
    let (cdf, sf, mut pdf);

    let mut use_upper_sum;

    if !(n > 0 && x >= 0.0 && x <= 1.0) {
        return ThreeProbs {
            sf: f64::NAN,
            cdf: f64::NAN,
            pdf: f64::NAN,
        };
    }
    if n == 1 {
        return ThreeProbs {
            sf: 1.0 - x,
            cdf: x,
            pdf: 1.0,
        };
    }
    if x == 0.0 {
        return ThreeProbs {
            sf: 1.0,
            cdf: 0.0,
            pdf: 1.0,
        };
    }
    if x == 1.0 {
        return ThreeProbs {
            sf: 0.0,
            cdf: 1.0,
            pdf: 0.0,
        };
    }

    let (alpha, nxfl, nx) = mod_nx(n, x);
    let mut n1mxfl = n - nxfl - if alpha == 0.0 { 0 } else { 1 };
    let mut n1mxceil = n - nxfl;
    /*
     * If alpha is 0, don't actually want to include the last term
     * in either the lower or upper summations.
     */
    if alpha == 0.0 {
        n1mxfl -= 1;
        n1mxceil += 1;
    }

    /* Special case:  x <= 1/n  */
    if nxfl == 0 || (nxfl == 1 && alpha == 0.0) {
        let t = pow2(1.0, x, n - 1);
        pdf = (nx + 1.0) * t / (1.0 + x);
        cdf = x * t;
        sf = 1.0 - cdf;
        /* Adjust if x=1/n *exactly* */
        if nxfl == 1 {
            debug_assert_eq!(alpha, 0.0);
            pdf -= 0.5;
        }
        return ThreeProbs { sf, cdf, pdf };
    }
    /* Special case:  x is so big, the sf underflows double64 */
    if -2.0 * f64::from(n) * x * x < MINLOG {
        return ThreeProbs {
            sf: 0.0,
            cdf: 1.0,
            pdf: 0.0,
        };
    }
    /* Special case:  x >= 1 - 1/n */
    if nxfl >= n - 1 {
        sf = pow2(1.0, -x, n);
        cdf = 1.0 - sf;
        pdf = f64::from(n) * sf / (1.0 - x);
        return ThreeProbs { sf, cdf, pdf };
    }
    /* Special case:  n is so big, take too long to compute */
    if n > SMIRNOV_MAX_COMPUTE_N {
        /* p ~ e^(-(6nx+1)^2 / 18n) */
        let t = 6.0 * f64::from(n) * x + 1.0;
        let logp = -(t * t) / 18.0 / f64::from(n);
        /* Maximise precision for small p-value. */
        if logp < -LN_2 {
            sf = logp.exp();
            cdf = 1.0 - sf;
        } else {
            cdf = -logp.exp_m1();
            sf = 1.0 - cdf;
        }
        pdf = (6.0 * f64::from(n) * x + 1.0) * 2.0 * sf / 3.0;
        return ThreeProbs { sf, cdf, pdf };
    }
    {
        /*
         * Use the upper sum if n is large enough, and x is small enough and
         * the number of terms is going to be small enough.
         * Otherwise it just drops accuracy, about 1.6bits * nUpperTerms
         */
        let n_upper_terms = n - n1mxceil + 1;
        use_upper_sum = n_upper_terms <= 1 && x < 0.5;
        use_upper_sum = use_upper_sum
            || ((n >= SM_UPPERSUM_MIN_N)
                && (n_upper_terms <= SM_UPPER_MAX_TERMS)
                && (x <= 0.5 / f64::from(n).sqrt()));
    }
    {
        let (mut start, mut step, mut n_terms) = (0, 1, n1mxfl + 1);
        let mut first_j = 0;
        let vmid = n / 2;
        let mut cman = DoubleDouble::new(1.0);
        let mut cexpt = 0;
        let (mut aj, mut daj_coeff);
        let one_over_x = DoubleDouble::new(1.0) / x;

        if use_upper_sum {
            start = n;
            step = -1;
            n_terms = n - n1mxceil + 1;

            let t1 = pow4_d(1.0, x, 1.0, 0.0, n - 1);
            aj = t1;

            daj_coeff = f64::from(n - 1) / (DoubleDouble::new(1.0) + x);
            daj_coeff = daj_coeff + one_over_x;
        } else {
            let t2 = pow4_d(1.0, -x, 1.0, 0.0, n);
            aj = t2 / x;

            daj_coeff =
                (-1.0 - DoubleDouble::new(f64::from(n - 1)) * x) / (DoubleDouble::new(1.0) - x);
            daj_coeff = daj_coeff / x;
            daj_coeff = daj_coeff + one_over_x;
        }

        let mut daj = aj * daj_coeff;
        aj_sum = aj_sum + aj;
        daj_sum = daj_sum + daj;

        update_binomial(&mut cman, &mut cexpt, n, 0);
        first_j += 1;

        for j in first_j..n_terms {
            let v = start + j * step;

            aj = compute_av(n, x, v, cman, cexpt);

            if aj.is_finite() && (aj != 0.0) {
                /* coeff = 1/x + (j-1)/(x+j/n) - (n-j)/(1-x-j/n) */
                daj_coeff = (f64::from(n) * f64::from(v - 1))
                    / (DoubleDouble::new(f64::from(nxfl + v)) + alpha)
                    - (f64::from(n - v) * f64::from(n))
                        / (DoubleDouble::new(f64::from(n - nxfl - v)) - alpha);
                daj_coeff = daj_coeff + one_over_x;
                daj = aj * daj_coeff;

                debug_assert!(aj.is_finite());
                aj_sum = aj_sum + aj;
                daj_sum = daj_sum + daj;
            }
            /* Safe to terminate early? */
            if aj != 0.0 {
                if ((f64::from(4 * (n_terms - j)) * f64::from(aj).abs())
                    < (f64::EPSILON * f64::from(aj_sum)))
                    && (j != n_terms - 1)
                {
                    break;
                }
            } else if j > vmid {
                debug_assert_eq!(aj, 0.0);
                break;
            }
            update_binomial(&mut cman, &mut cexpt, n, j);
        }
        debug_assert!(aj_sum.is_finite());
        debug_assert!(daj_sum.is_finite());
        {
            let deriv_d = x * daj_sum;
            let prob_d = x * aj_sum;
            let deriv = f64::from(deriv_d);
            let prob = f64::from(prob_d);

            debug_assert!(nx != 1.0 || alpha > 0.0);
            if step < 0 {
                cdf = prob;
                sf = 1.0 - prob;
                pdf = deriv;
            } else {
                cdf = 1.0 - prob;
                sf = prob;
                pdf = -deriv;
            }
        }
    }
    ThreeProbs {
        sf: sf.clamp(0.0, 1.0),
        cdf: cdf.clamp(0.0, 1.0),
        pdf: 0.0_f64.max(pdf),
    }
}

/*
 * Functional inverse of Smirnov distribution
 * finds x such that smirnov(n, x) = psf; smirnovc(n, x) = pcdf).
 */
#[allow(
    clippy::float_cmp,
    clippy::manual_range_contains,
    clippy::too_many_lines
)]
fn smirnovi_(n: i32, psf: f64, pcdf: f64) -> f64 {
    /*
     * Need to use a bracketing NR algorithm here and be very careful
     *  about the starting point.
     */
    let mut x;
    let mut iterations = 0;
    let (mut a, mut b);
    let (mut dx, mut dxold);

    if !(n > 0 && psf >= 0.0 && pcdf >= 0.0 && pcdf <= 1.0 && psf <= 1.0) {
        return f64::NAN;
    }
    if (1.0 - pcdf - psf).abs() > 4.0 * f64::EPSILON {
        return f64::NAN;
    }
    /* STEP 1: Handle psf==0, or pcdf == 0 */
    if pcdf == 0.0 {
        return 0.0;
    }
    if psf == 0.0 {
        return 1.0;
    }
    /* STEP 2: Handle n=1 */
    if n == 1 {
        return pcdf;
    }

    /* STEP 3 Handle psf *very* close to 0.  Correspond to (n-1)/n < x < 1  */
    let psfrootn = psf.powf(1.0 / f64::from(n));
    /* xmin > 1 - 1.0 / n */
    if n < 150 && f64::from(n) * psfrootn <= 1.0 {
        /* Solve exactly. */
        x = 1.0 - psfrootn;
        return x;
    }

    let logpcdf = if pcdf < 0.5 {
        pcdf.ln()
    } else {
        (-psf).ln_1p()
    };

    /*
     * STEP 4 Find bracket and initial estimate for use in N-R
     * 4(a)  Handle 0 < x <= 1/n:   pcdf = x * (1+x)^*(n-1)
     */
    let maxlogpcdf = logpow4(1.0, 0.0, f64::from(n), 0.0, 1)
        + logpow4(f64::from(n), 1.0, f64::from(n), 0.0, n - 1);
    if logpcdf <= maxlogpcdf {
        let xmin = pcdf / SCIPY_EL;
        let xmax = pcdf;
        let p1 = pow4(f64::from(n), 1.0, f64::from(n), 0.0, n - 1) / f64::from(n);
        let mut r = pcdf / p1;
        let mut z0 = r;
        /*
         * Do one iteration of N-R solving: z*e^(z-1) = R, with z0=pcdf/P1
         *  z  <-  z - (z exp(z-1) - pcdf)/((z+1)exp(z-1))
         *  If z_0 = R, z_1 = R(1-exp(1-R))/(R+1)
         */
        if r >= 1.0 {
            /*
             * R=1 is OK;
             * R>1 can happen due to truncation error for x = (1-1/n)+-eps
             */
            r = 1.0;
            x = r / f64::from(n);
            return x;
        }
        z0 = (z0 * z0 + r * (1.0 - z0).exp()) / (1.0 + z0);
        x = z0 / f64::from(n);
        a = xmin * (1.0 - 4.0 * f64::EPSILON);
        a = a.max(0.0);
        b = xmax * (1.0 + 4.0 * f64::EPSILON);
        b = b.min(1.0 / f64::from(n));
        // std::clamp(x, a, b)
        x = if x < a {
            a
        } else if b < x {
            b
        } else {
            x
        };
    } else {
        /* 4(b) : 1/n < x < (n-1)/n */
        let xmin = 1.0 - psfrootn;
        let logpsf = if psf < 0.5 { psf.ln() } else { (-pcdf).ln_1p() };
        let xmax = (-logpsf / (2.0 * f64::from(n))).sqrt();
        let xmax6 = xmax - 1.0 / (6.0 * f64::from(n));
        a = xmin.max(1.0 / f64::from(n));
        b = xmax.min(1.0 - 1.0 / f64::from(n));
        x = xmax6;
    }
    if x < a || x > b {
        x = f64::midpoint(a, b);
    }
    debug_assert!(x < 1.0);

    /*
     * Skip computing fa, fb as that takes cycles and the exact values
     * are not needed.
     */

    /* STEP 5 Run N-R.
     * smirnov should be well-enough behaved for NR starting at this location.
     * Use smirnov(n, x)-psf, or pcdf - smirnovc(n, x), whichever has smaller p.
     */
    dxold = b - a;
    dx = dxold;
    loop {
        let x0 = x;
        let (dfdx, deltax, df);
        debug_assert!(x < 1.0);
        debug_assert!(x > 0.0);
        {
            let probs = smirnov_(n, x0);
            df = if pcdf < 0.5 {
                pcdf - probs.cdf
            } else {
                probs.sf - psf
            };
            dfdx = -probs.pdf;
        }
        if df == 0.0 {
            return x;
        }
        /* Update the bracketing interval */
        if df > 0.0 && x > a {
            a = x;
        } else if df < 0.0 && x < b {
            b = x;
        }

        if dfdx == 0.0 {
            /*
             * x was not within tolerance, but now we hit a 0 derivative.
             * This implies that x >> 1/sqrt(n), and even then |smirnovp| >= |smirnov|
             * so this condition is unexpected.  Do a bisection step.
             */
            x = f64::midpoint(a, b);
            deltax = x0 - x;
        } else {
            deltax = df / dfdx;
            x = x0 - deltax;
        }
        /*
         * Check out-of-bounds.
         * Not expecting this to happen ofen --- smirnov is convex near x=1 and
         * concave near x=0, and we should be approaching from the correct side.
         * If out-of-bounds, replace x with a midpoint of the bracket.
         * Also check fast enough convergence.
         */
        if (a <= x)
            && (x <= b)
            && ((2.0 * deltax).abs() <= dxold.abs() || dxold.abs() < 256.0 * f64::EPSILON)
        {
            dxold = dx;
            dx = deltax;
        } else {
            dxold = dx;
            dx /= 2.0;
            x = f64::midpoint(a, b);
        }
        /*
         * Note that if psf is close to 1, f(x) -> 1, f'(x) -> -1.
         *  => abs difference |x-x0| is approx |f(x)-p| >= std::numeric_limits<double>::epsilon(),
         *  => |x-x0|/x >= std::numeric_limits<double>::epsilon()/x.
         *  => cannot use a purely relative criteria as it will fail for x close to 0.
         */
        if within_tol(x, x0, if psf < 0.5 { 0.0 } else { XTOL }, RTOL) {
            break;
        }
        iterations += 1;
        if iterations > KOLMOG_MAXITER {
            return x;
        }
    }
    x
}

pub(crate) fn smirnov(n: i32, d: f64) -> f64 {
    if d.is_nan() {
        return f64::NAN;
    }
    smirnov_(n, d).sf
}

pub(crate) fn smirnovc(n: i32, d: f64) -> f64 {
    if d.is_nan() {
        return f64::NAN;
    }
    smirnov_(n, d).cdf
}

/// Derivative of smirnov(n, d)
///  One interior point of discontinuity at d=1/n.
#[allow(clippy::float_cmp, clippy::manual_range_contains)]
pub(crate) fn smirnovp(n: i32, d: f64) -> f64 {
    if !(n > 0 && d >= 0.0 && d <= 1.0) {
        return f64::NAN;
    }
    if n == 1 {
        /* Slope is always -1 for n=1, even at d = 1.0 */
        return -1.0;
    }
    if d == 1.0 {
        return -0.0;
    }
    /*
     * If d is 0, the derivative is discontinuous, but approaching
     * from the right the limit is -1
     */
    if d == 0.0 {
        return -1.0;
    }
    -smirnov_(n, d).pdf
}

pub(crate) fn smirnovi(n: i32, p: f64) -> f64 {
    if p.is_nan() {
        return f64::NAN;
    }
    smirnovi_(n, p, 1.0 - p)
}

pub(crate) fn smirnovci(n: i32, p: f64) -> f64 {
    if p.is_nan() {
        return f64::NAN;
    }
    smirnovi_(n, 1.0 - p, p)
}
