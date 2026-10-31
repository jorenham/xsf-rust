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

    if (-TUKEY_SMALLVAL < lmbda) && (lmbda < TUKEY_SMALLVAL) {
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

    while (count < TUKEY_MAXCOUNT) && ((pmid - plow).abs() > TUKEY_EPS) {
        xeval = (pmid.powf(lmbda) - (1.0 - pmid).powf(lmbda)) / lmbda;
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
