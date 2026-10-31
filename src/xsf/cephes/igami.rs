//! Translated into pure Rust from `xsf/cephes/igami.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024. Original header comment appears below.
//!
//! (C) Copyright John Maddock 2006.
//! Use, modification and distribution are subject to the
//! Boost Software License, Version 1.0. (See accompanying file
//!  LICENSE_1_0.txt or copy at <https://www.boost.org/LICENSE_1_0.txt>)
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::utils::SCIPY_EULER;
use crate::xsf::cephes::igam::igam_fac;
use crate::xsf::cephes::polevl::polevl;
use crate::xsf::cephes::{gamma, igam, igamc, lgam};

fn find_inverse_s(p: f64, q: f64) -> f64 {
    /*
     * Computation of the Incomplete Gamma Function Ratios and their Inverse
     * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
     * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
     * December 1986, Pages 377-393.
     *
     * See equation 32.
     */
    const A: [f64; 4] = [
        0.213_623_493_715_853,
        4.283_421_559_671_04,
        11.661_672_028_896_8,
        3.311_259_221_087_41,
    ];
    const B: [f64; 5] = [
        0.036_117_081_018_842_03,
        1.273_644_897_822_23,
        6.406_915_977_600_39,
        6.610_537_656_254_62,
        1.0,
    ];

    let t = if p < 0.5 {
        (-2.0 * p.ln()).sqrt()
    } else {
        (-2.0 * q.ln()).sqrt()
    };
    let s = t - polevl(t, &A) / polevl(t, &B);
    if p < 0.5 { -s } else { s }
}

fn didonato_sn(a: f64, x: f64, n: u32, tolerance: f64) -> f64 {
    /*
     * Computation of the Incomplete Gamma Function Ratios and their Inverse
     * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
     * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
     * December 1986, Pages 377-393.
     *
     * See equation 34.
     */
    let mut sum = 1.0;

    if n >= 1 {
        let mut partial = x / (a + 1.0);

        sum += partial;
        for i in 2..=n {
            partial *= x / (a + f64::from(i));
            sum += partial;
            if partial < tolerance {
                break;
            }
        }
    }
    sum
}

/// DiDonato and Morris Eq 25
fn didonato_eq25(a: f64, y: f64) -> f64 {
    let c1 = (a - 1.0) * y.ln();
    let c1_2 = c1 * c1;
    let c1_3 = c1_2 * c1;
    let c1_4 = c1_2 * c1_2;
    let a_2 = a * a;
    let a_3 = a_2 * a;

    let c2 = (a - 1.0) * (1.0 + c1);
    let c3 = (a - 1.0) * (-(c1_2 / 2.0) + (a - 2.0) * c1 + (3.0 * a - 5.0) / 2.0);
    let c4 = (a - 1.0)
        * ((c1_3 / 3.0) - (3.0 * a - 5.0) * c1_2 / 2.0
            + (a_2 - 6.0 * a + 7.0) * c1
            + (11.0 * a_2 - 46.0 * a + 47.0) / 6.0);
    let c5 = (a - 1.0)
        * (-(c1_4 / 4.0)
            + (11.0 * a - 17.0) * c1_3 / 6.0
            + (-3.0 * a_2 + 13.0 * a - 13.0) * c1_2
            + (2.0 * a_3 - 25.0 * a_2 + 72.0 * a - 61.0) * c1 / 2.0
            + (25.0 * a_3 - 195.0 * a_2 + 477.0 * a - 379.0) / 12.0);

    let y_2 = y * y;
    let y_3 = y_2 * y;
    let y_4 = y_2 * y_2;
    y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4)
}

#[allow(clippy::float_cmp)]
fn find_inverse_gamma(a: f64, p: f64, q: f64) -> f64 {
    /*
     * In order to understand what's going on here, you will
     * need to refer to:
     *
     * Computation of the Incomplete Gamma Function Ratios and their Inverse
     * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
     * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
     * December 1986, Pages 377-393.
     */
    if a == 1.0 {
        if q > 0.9 { -(-p).ln_1p() } else { -q.ln() }
    } else if a < 1.0 {
        let g = gamma(a);
        let b = q * g;

        if (b > 0.6) || ((b >= 0.45) && (a >= 0.3)) {
            /* DiDonato & Morris Eq 21:
             *
             * There is a slight variation from DiDonato and Morris here:
             * the first form given here is unstable when p is close to 1,
             * making it impossible to compute the inverse of Q(a,x) for small
             * q. Fortunately the second form works perfectly well in this case.
             */
            let u = if (b * q > 1e-8) && (q > 1e-5) {
                (p * g * a).powf(1.0 / a)
            } else {
                ((-q / a) - SCIPY_EULER).exp()
            };
            u / (1.0 - (u / (a + 1.0)))
        } else if (a < 0.3) && (b >= 0.35) {
            /* DiDonato & Morris Eq 22: */
            let t = (-SCIPY_EULER - b).exp();
            let u = t * t.exp();
            t * u.exp()
        } else if (b > 0.15) || (a >= 0.3) {
            /* DiDonato & Morris Eq 23: */
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            y - (1.0 - a) * u.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln()
        } else if b > 0.1 {
            /* DiDonato & Morris Eq 24: */
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            y - (1.0 - a) * u.ln()
                - ((u * u + 2.0 * (3.0 - a) * u + (2.0 - a) * (3.0 - a))
                    / (u * u + (5.0 - a) * u + 2.0))
                    .ln()
        } else {
            /* DiDonato & Morris Eq 25: */
            didonato_eq25(a, -b.ln())
        }
    } else {
        /* DiDonato and Morris Eq 31: */
        let s = find_inverse_s(p, q);

        let s_2 = s * s;
        let s_3 = s_2 * s;
        let s_4 = s_2 * s_2;
        let s_5 = s_4 * s;
        let ra = a.sqrt();

        let mut w = a + s * ra + (s_2 - 1.0) / 3.0;
        w += (s_3 - 7.0 * s) / (36.0 * ra);
        w -= (3.0 * s_4 + 7.0 * s_2 - 16.0) / (810.0 * a);
        w += (9.0 * s_5 + 256.0 * s_3 - 433.0 * s) / (38880.0 * a * ra);

        if (a >= 500.0) && ((1.0 - w / a).abs() < 1e-6) {
            w
        } else if p > 0.5 {
            if w < 3.0 * a {
                w
            } else {
                let d = 2.0_f64.max(a * (a - 1.0));
                let lg = lgam(a);
                let lb = q.ln() + lg;
                if lb < -d * 2.3 {
                    /* DiDonato and Morris Eq 25: */
                    didonato_eq25(a, -lb)
                } else {
                    /* DiDonato and Morris Eq 33: */
                    let u = -lb + (a - 1.0) * w.ln() - (1.0 + (1.0 - a) / (1.0 + w)).ln();
                    -lb + (a - 1.0) * u.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln()
                }
            }
        } else {
            let mut z = w;
            let ap1 = a + 1.0;
            let ap2 = a + 2.0;
            if w < 0.15 * ap1 {
                /* DiDonato and Morris Eq 35: */
                let v = p.ln() + lgam(ap1);
                z = ((v + w) / a).exp();
                let mut s = (z / ap1 * (1.0 + z / ap2)).ln_1p();
                z = ((v + z - s) / a).exp();
                s = (z / ap1 * (1.0 + z / ap2)).ln_1p();
                z = ((v + z - s) / a).exp();
                s = (z / ap1 * (1.0 + z / ap2 * (1.0 + z / (a + 3.0)))).ln_1p();
                z = ((v + z - s) / a).exp();
            }

            if (z <= 0.01 * ap1) || (z > 0.7 * ap1) {
                z
            } else {
                /* DiDonato and Morris Eq 36: */
                let ls = didonato_sn(a, z, 100, 1e-4).ln();
                let v = p.ln() + lgam(ap1);
                z = ((v + z - ls) / a).exp();
                z * (1.0 - (a * z.ln() - z - v + ls) / (a - z))
            }
        }
    }
}

/// Inverse of the regularized lower incomplete gamma function
#[allow(clippy::float_cmp)]
pub(crate) fn igami(a: f64, p: f64) -> f64 {
    if a.is_nan() || p.is_nan() {
        return f64::NAN;
    } else if (a < 0.0) || (p < 0.0) || (p > 1.0) {
        // domain error (no early return, as in xsf)
    } else if p == 0.0 {
        return 0.0;
    } else if p == 1.0 {
        return f64::INFINITY;
    } else if p > 0.9 {
        return igamci(a, 1.0 - p);
    }

    let mut x = find_inverse_gamma(a, p, 1.0 - p);
    /* Halley's method */
    for _ in 0..3 {
        let fac = igam_fac(a, x);
        if fac == 0.0 {
            return x;
        }
        let f_fp = (igam(a, x) - p) * x / fac;
        /* The ratio of the first and second derivatives simplifies */
        let fpp_fp = -1.0 + (a - 1.0) / x;
        if fpp_fp.is_infinite() {
            /* Resort to Newton's method in the case of overflow */
            x -= f_fp;
        } else {
            x -= f_fp / (1.0 - 0.5 * f_fp * fpp_fp);
        }
    }

    x
}

/// Inverse of the regularized upper incomplete gamma function
#[allow(clippy::float_cmp)]
pub(crate) fn igamci(a: f64, q: f64) -> f64 {
    if a.is_nan() || q.is_nan() {
        return f64::NAN;
    } else if (a < 0.0) || (q < 0.0) || (q > 1.0) {
        // domain error (no early return, as in xsf)
    } else if q == 0.0 {
        return f64::INFINITY;
    } else if q == 1.0 {
        return 0.0;
    } else if q > 0.9 {
        return igami(a, 1.0 - q);
    }

    let mut x = find_inverse_gamma(a, 1.0 - q, q);
    for _ in 0..3 {
        let fac = igam_fac(a, x);
        if fac == 0.0 {
            return x;
        }
        let f_fp = (igamc(a, x) - q) * x / (-fac);
        let fpp_fp = -1.0 + (a - 1.0) / x;
        if fpp_fp.is_infinite() {
            x -= f_fp;
        } else {
            x -= f_fp / (1.0 - 0.5 * f_fp * fpp_fp);
        }
    }

    x
}
