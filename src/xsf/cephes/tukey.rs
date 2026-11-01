//! Translated into pure Rust from `xsf/cephes/tukey.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024. Original header comment appears below.
//!
//! Compute the CDF of the Tukey-Lambda distribution
//! using a bracketing search with special checks
//!
//! The PPF of the Tukey-lambda distribution is
//! G(p) = (p**lam + (1-p)**lam) / lam
//!
//! Author:  Travis Oliphant

const TUKEY_SMALLVAL: f64 = 1e-4;
const TUKEY_EPS: f64 = 1.0e-18;
const TUKEY_MAXCOUNT: i32 = 60;

#[allow(clippy::float_cmp)]
pub(crate) fn tukeylambdacdf(x: f64, lmbda: f64) -> f64 {
    if x.is_nan() || lmbda.is_nan() {
        return f64::NAN;
    }

    let mut xeval = 1.0 / lmbda;
    if lmbda > 0.0 {
        if x <= (-xeval) {
            return 0.0;
        }
        if x >= xeval {
            return 1.0;
        }
    }

    // Unlike in xsf, the logistic CDF only for (subnormal) lambda = 0, instead of |lambda| < 1e-4,
    // as its relative error is ~lambda x^2 / 2, e.g. `tukeylambdacdf(-20.0, 9e-5)` was 2.06e-9,
    // not 2.02e-9.
    if lmbda.abs() < f64::MIN_POSITIVE {
        if x >= 0.0 {
            return 1.0 / (1.0 + (-x).exp());
        }
        return x.exp() / (1.0 + x.exp());
    }

    let pmin: f64 = 0.0;
    let mut pmid: f64 = 0.5;
    let pmax = 1.0;
    let mut plow = pmin;
    let mut phigh = pmax;
    let mut count = 0;

    // Unlike in xsf, a relative tolerance (and more iterations for small p), since e.g.
    // `tukeylambdacdf(-4.6870665239279194, 0.21296947909944208)` was 1.287382e-13, not 1.287389e-13.
    while (f64::from(count) < f64::from(TUKEY_MAXCOUNT) - pmid.log2())
        && ((pmid - plow).abs() > TUKEY_EPS * pmid)
    {
        xeval = if lmbda.abs() < TUKEY_SMALLVAL {
            // Unlike in xsf, without the cancellation of p^lambda - (1-p)^lambda for small lambda
            ((lmbda * pmid.ln()).exp_m1() - (lmbda * (-pmid).ln_1p()).exp_m1()) / lmbda
        } else {
            (pmid.powf(lmbda) - (1.0 - pmid).powf(lmbda)) / lmbda
        };
        if xeval == x {
            return pmid;
        }
        if xeval > x {
            phigh = pmid;
            pmid = f64::midpoint(pmid, plow);
        } else {
            plow = pmid;
            pmid = f64::midpoint(pmid, phigh);
        }
        count += 1;
    }
    pmid
}
