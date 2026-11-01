//! Translated into pure Rust from `xsf/cephes/expn.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 1.1:  March, 1985
//! Copyright 1985 by Stephen L. Moshier
//!
//! Sources
//! \[1\] NIST, "The Digital Library of Mathematical Functions", dlmf.nist.gov
//!
//! Scipy changes:
//! - 09-10-2016: improved asymptotic expansion for large n
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::utils::{MACHEP, MAXLOG, SCIPY_EULER};
use crate::xsf::cephes::polevl::polevl;
use crate::xsf::cephes::rgamma;

const EXPN_NA: usize = 13;
const EXPN_A0: [f64; 1] = [1.0];
const EXPN_A1: [f64; 1] = [1.0];
const EXPN_A2: [f64; 2] = [-2.0, 1.0];
const EXPN_A3: [f64; 3] = [6.0, -8.0, 1.0];
const EXPN_A4: [f64; 4] = [-24.0, 58.0, -22.0, 1.0];
const EXPN_A5: [f64; 5] = [120.0, -444.0, 328.0, -52.0, 1.0];
const EXPN_A6: [f64; 6] = [-720.0, 3_708.0, -4_400.0, 1_452.0, -114.0, 1.0];
const EXPN_A7: [f64; 7] = [
    5_040.0, -33_984.0, 58_140.0, -32_120.0, 5_610.0, -240.0, 1.0,
];
const EXPN_A8: [f64; 8] = [
    -40_320.0, 341_136.0, -785_304.0, 644_020.0, -195_800.0, 19_950.0, -494.0, 1.0,
];
const EXPN_A9: [f64; 9] = [
    362_880.0,
    -3_733_920.0,
    11_026_296.0,
    -12_440_064.0,
    5_765_500.0,
    -1_062_500.0,
    67_260.0,
    -1_004.0,
    1.0,
];
const EXPN_A10: [f64; 10] = [
    -3_628_800.0,
    44_339_040.0,
    -162_186_912.0,
    238_904_904.0,
    -155_357_384.0,
    44_765_000.0,
    -5_326_160.0,
    218_848.0,
    -2_026.0,
    1.0,
];
const EXPN_A11: [f64; 11] = [
    39_916_800.0,
    -568_356_480.0,
    2_507_481_216.0,
    -4_642_163_952.0,
    4_002_695_088.0,
    -1_648_384_304.0,
    314_369_720.0,
    -25_243_904.0,
    695_038.0,
    -4_072.0,
    1.0,
];
const EXPN_A12: [f64; 12] = [
    -479_001_600.0,
    7_827_719_040.0,
    -40_788_301_824.0,
    92_199_790_224.0,
    -101_180_433_024.0,
    56_041_398_784.0,
    -15_548_960_784.0,
    2_051_482_776.0,
    -114_876_376.0,
    2_170_626.0,
    -8_166.0,
    1.0,
];
const EXPN_A: [&[f64]; EXPN_NA] = [
    &EXPN_A0, &EXPN_A1, &EXPN_A2, &EXPN_A3, &EXPN_A4, &EXPN_A5, &EXPN_A6, &EXPN_A7, &EXPN_A8,
    &EXPN_A9, &EXPN_A10, &EXPN_A11, &EXPN_A12,
];

/* Asymptotic expansion for large n, DLMF 8.20(ii) */
#[allow(clippy::float_cmp)]
fn expn_large_n(n: u32, x: f64) -> f64 {
    let p = f64::from(n);
    let lambda = x / p;
    let multiplier = 1.0 / p / (lambda + 1.0) / (lambda + 1.0);
    let mut fac = 1.0;
    let mut res = 1.0; /* A[0] = 1 */

    let expfac = (-lambda * p).exp() / (lambda + 1.0) / p;
    if expfac == 0.0 {
        return 0.0;
    }

    /* Do the k = 1 term outside the loop since A[1] = 1 */
    fac *= multiplier;
    res += fac;

    for a in &EXPN_A[2..] {
        fac *= multiplier;
        let term = fac * polevl(lambda, a);
        res += term;
        if term.abs() < MACHEP * res.abs() {
            break;
        }
    }

    expfac * res
}

/// Generalized exponential integral $E_n(x)$
///
/// Corresponds to [`scipy.special.expn`][expn]
///
/// [expn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expn.html
///
/// # Definition
/// For integer $n \geq 0$ and real $x \geq 0$ the
/// generalized exponential integral is defined as:
///
/// $$ E_n(x) = x^{n-1} \int_x^\infty \frac{e^{-t}}{t^n} \dd t $$
///
/// # Accuracy
///
/// The routine employs either a power series, a continued
/// fraction, or an asymptotic formula depending on the
/// relative values of n and x.
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       10000       1.6e-15     3.6e-16
/// ```
///
/// for uniformly distributed integer $1 \le n \le 50$. For $50 < n \le 200$ (the asymptotic
/// expansion) and $0 < x < 30$, the peak and rms relative errors are 4.1e-15 and 9.0e-16.
///
/// # See also
/// - [`exp1`](fn.exp1.html): Special case of $E_n$ for $n = 1$
/// - [`expi`](fn.expi.html): Related to $E_n$ when $n = 1$
#[must_use]
#[allow(
    clippy::float_cmp,
    clippy::if_not_else,
    clippy::neg_cmp_op_on_partial_ord
)]
pub fn expn(n: u32, x: f64) -> f64 {
    const BIG: f64 = 1.441_151_880_758_558_7e17;

    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }

    if x > MAXLOG {
        return 0.0;
    }

    if x == 0.0 {
        if n < 2 {
            return f64::INFINITY;
        }
        return 1.0 / (f64::from(n) - 1.0);
    }

    if n == 0 {
        return (-x).exp() / x;
    }

    /* Asymptotic expansion for large n, DLMF 8.20(ii) */
    if n > 50 {
        return expn_large_n(n, x);
    }

    /* Continued fraction, DLMF 8.19.17 */
    if x > 1.0 {
        let mut k: u32 = 1;
        let mut pkm2 = 1.0;
        let mut qkm2 = x;
        let mut pkm1 = 1.0;
        let mut qkm1 = x + f64::from(n);
        let mut ans = pkm1 / qkm1;

        loop {
            k += 1;
            let (yk, xk) = if k & 1 != 0 {
                (1.0, f64::from(n + (k - 1) / 2))
            } else {
                (x, f64::from(k / 2))
            };
            let pk = pkm1 * yk + pkm2 * xk;
            let qk = qkm1 * yk + qkm2 * xk;
            let t;
            if qk != 0.0 {
                let r = pk / qk;
                t = ((ans - r) / r).abs();
                ans = r;
            } else {
                t = 1.0;
            }
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            if pk.abs() > BIG {
                pkm2 /= BIG;
                pkm1 /= BIG;
                qkm2 /= BIG;
                qkm1 /= BIG;
            }
            if !(t > MACHEP) {
                break;
            }
        }

        ans *= (-x).exp();
        return ans;
    }

    /* Power series expansion, DLMF 8.19.8 */
    let mut psi = -SCIPY_EULER - x.ln();
    for i in 1..n {
        psi += 1.0 / f64::from(i);
    }

    let z = -x;
    let mut xk = 0.0;
    let mut yk = 1.0;
    let mut pk = 1.0 - f64::from(n);
    let mut ans = if n == 1 { 0.0 } else { 1.0 / pk };
    loop {
        xk += 1.0;
        yk *= z / xk;
        pk += 1.0;
        if pk != 0.0 {
            ans += yk / pk;
        }
        let t = if ans != 0.0 { (yk / ans).abs() } else { 1.0 };
        if !(t > MACHEP) {
            break;
        }
    }
    let t = f64::from(n);
    let r = f64::from(n - 1);
    (z.powf(r) * psi * rgamma(t)) - ans
}
