//! Mathieu functions, translated into pure Rust from `xsf/specfun/specfun.h` (xsf v0.2.2).
//!
//! The `cv0` polynomials in `q` are stored as a coefficient table (highest degree first) and
//! evaluated with `polevl`, instead of the unrolled `if` branches in the C++ code, which evaluates
//! them in the same order.
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.

use super::{Status, jynb, try_zeroed};
use crate::xsf::cephes::polevl::polevl;

/// The `cv0` polynomials in `q` (highest degree first), as `(m, kd, q_max, coefficients)`
#[allow(clippy::unreadable_literal)]
#[rustfmt::skip]
const CV0_POLY: [(f64, f64, f64, &[f64]); 25] = [
    (0.0, 1.0, 10.0, &[3.999267e-3, -9.638957e-2, -0.88297, 0.5542818]),
    (1.0, 2.0, 10.0, &[-4.94603e-4, 1.92917e-2, -0.3089229, 1.33372, 0.811752]),
    (1.0, 3.0, 10.0, &[1.971096e-3, -5.482465e-2, -1.152218, 1.10427]),
    (2.0, 1.0, 15.0, &[3.200972e-4, -8.667445e-3, -1.829032e-4, 0.9919999, 3.3290504]),
    (2.0, 4.0, 10.0, &[2.38446e-3, -0.08725329, -4.732542e-3, 4.00909]),
    (3.0, 2.0, 20.0, &[3.035731e-4, -1.453021e-2, 0.19069602, -0.1039356, 8.9449274]),
    (3.0, 3.0, 15.0, &[9.369364e-5, -0.03569325, 0.2689874, 8.771735]),
    (4.0, 1.0, 25.0, &[1.076676e-4, -7.9684875e-3, 0.17344854, -0.5924058, 16.620847]),
    (4.0, 4.0, 20.0, &[-7.08719e-4, 3.8216144e-3, 0.1907493, 15.744]),
    (5.0, 2.0, 35.0, &[2.238231e-5, -2.983416e-3, 0.10706975, -0.600205, 25.93515]),
    (5.0, 3.0, 25.0, &[-7.425364e-4, 2.18225e-2, 4.16399e-2, 24.897]),
    (6.0, 1.0, 40.0, &[-1.66846e-5, 4.80263e-4, 2.53998e-2, -0.181233, 36.423]),
    (6.0, 4.0, 35.0, &[-4.57146e-4, 2.16609e-2, -2.349616e-2, 35.99251]),
    (7.0, 2.0, 50.0, &[-1.411114e-5, 9.730514e-4, -3.097887e-3, 3.533597e-2, 49.0547]),
    (7.0, 3.0, 40.0, &[-3.043872e-4, 2.05511e-2, -9.16292e-2, 49.19035]),
    // for m >= 8, q_max is m*m
    (8.0, 1.0, 64.0, &[8.634308e-6, -2.100289e-3, 0.169072, -4.64336, 109.4211]),
    (8.0, 4.0, 64.0, &[-6.7842e-5, 2.2057e-3, 0.48296, 56.59]),
    (9.0, 2.0, 81.0, &[2.906435e-6, -1.019893e-3, 0.1101965, -3.821851, 127.6098]),
    (9.0, 3.0, 81.0, &[-9.577289e-5, 0.01043839, 0.06588934, 78.0198]),
    (10.0, 1.0, 100.0, &[5.44927e-7, -3.926119e-4, 0.0612099, -2.600805, 138.1923]),
    (10.0, 4.0, 100.0, &[-7.660143e-5, 0.01132506, -0.09746023, 99.29494]),
    (11.0, 2.0, 121.0, &[-5.67615e-7, 7.152722e-6, 0.01920291, -1.081583, 140.88]),
    (11.0, 3.0, 121.0, &[-6.310551e-5, 0.0119247, -0.2681195, 123.667]),
    (12.0, 1.0, 144.0, &[-2.38351e-7, -2.90139e-5, 0.02023088, -1.289, 171.2723]),
    (12.0, 4.0, 144.0, &[3.08902e-7, -1.577869e-4, 0.0247911, -1.05454, 161.471]),
];

#[allow(
    clippy::float_cmp,
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::cast_possible_truncation
)]
pub(crate) fn cv0(kd: f64, m: f64, q: f64) -> f64 {
    // =====================================================
    // Purpose: Compute the initial characteristic value of
    //          Mathieu functions for m ≤ 12  or q ≤ 300 or
    //          q ≥ m*m
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    // Output:  A0 --- Characteristic value
    // Routines called:
    //       (1) CVQM for computing initial characteristic
    //           value for q ≤ 3*m
    //       (2) CVQL for computing initial characteristic
    //           value for q ≥ m*m
    // ====================================================

    // NOTE: `kd` and `m` are `double` here (the only caller, `cva2`, passes `int`s), so the calls to `cvql` and
    // `cvqm` (which take `int`) implicitly convert them back with `(int)`; this is exact for those values.
    let q2 = q * q;
    let poly = CV0_POLY.iter().find(|&&(pm, pkd, ..)| pm == m && pkd == kd);

    if m < 8.0 {
        let a0 = match (m as i32, kd as i32) {
            (0, _) if q <= 1.0 => (((0.0036392 * q2 - 0.0125868) * q2 + 0.0546875) * q2 - 0.5) * q2,
            (1, 2) if q <= 1.0 => (((-6.51e-4 * q - 0.015625) * q - 0.125) * q + 1.0) * q + 1.0,
            (1, 3) if q <= 1.0 => (((-6.51e-4 * q + 0.015625) * q - 0.125) * q - 1.0) * q + 1.0,
            (2, 1) if q <= 1.0 => {
                (((-0.0036391 * q2 + 0.0125888) * q2 - 0.0551939) * q2 + 0.416667) * q2 + 4.0
            },
            (2, 4) if q <= 1.0 => (0.0003617 * q2 - 0.0833333) * q2 + 4.0,
            (3, 2) if q <= 1.0 => ((6.348e-4 * q + 0.015625) * q + 0.0625) * q2 + 9.0,
            (3, 3) if q <= 1.0 => ((6.348e-4 * q - 0.015625) * q + 0.0625) * q2 + 9.0,
            (4, 1) if q <= 1.0 => ((-2.1e-6 * q2 + 5.012e-4) * q2 + 0.0333333) * q2 + 16.0,
            (4, 4) if q <= 1.0 => ((3.7e-6 * q2 - 3.669e-4) * q2 + 0.0333333) * q2 + 16.0,
            (5, 2) if q <= 1.0 => ((6.8e-6 * q + 1.42e-5) * q2 + 0.0208333) * q2 + 25.0,
            (5, 3) if q <= 1.0 => ((-6.8e-6 * q + 1.42e-5) * q2 + 0.0208333) * q2 + 25.0,
            (6, _) if q <= 1.0 => (4e-6 * q2 + 0.0142857) * q2 + 36.0,
            (7, _) if q <= 10.0 => cvqm(m as i32, q),
            _ => match poly {
                Some(&(.., q_max, c)) if q <= q_max => polevl(q, c),
                _ => cvql(kd as i32, m as i32, q),
            },
        };
        return a0;
    }
    if q <= 3.0 * m {
        cvqm(m as i32, q)
    } else if q > m * m {
        cvql(kd as i32, m as i32, q)
    } else if let Some(&(.., c)) = poly {
        polevl(q, c)
    } else {
        // NOTE: for `m > 12` (unreachable from `cva2`), or a `kd` that doesn't match the parity of `m`,
        // `a0` stays `0.0`
        0.0
    }
}

#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::assign_op_pattern
)]
pub(crate) fn cva2(kd: i32, m: i32, q: f64) -> f64 {
    // ======================================================
    // Purpose: Calculate a specific characteristic value of
    //          Mathieu functions
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    //          KD --- Case code
    //                 KD=1 for cem(x,q)  ( m = 0,2,4,...)
    //                 KD=2 for cem(x,q)  ( m = 1,3,5,...)
    //                 KD=3 for sem(x,q)  ( m = 1,3,5,...)
    //                 KD=4 for sem(x,q)  ( m = 2,4,6,...)
    // Output:  A  --- Characteristic value
    // Routines called:
    //       (1) REFINE for finding accurate characteristic
    //           value using an iteration method
    //       (2) CV0 for finding initial characteristic
    //           values using polynomial approximation
    //       (3) CVQM for computing initial characteristic
    //           values for q ≤ 3*m
    //       (3) CVQL for computing initial characteristic
    //           values for q ≥ m*m
    // ======================================================

    let mut a = 0.0;

    // NOTE: the `int` product `m * m` overflows for `m > 46340`, which is UB in C++ (wrapping on x86-64)
    if (m <= 12) || (q <= 3.0 * f64::from(m)) || (q > f64::from(m.wrapping_mul(m))) {
        a = cv0(f64::from(kd), f64::from(m), q);
        if (q != 0.0) && (m != 2) {
            a = refine(kd, m, q, a);
        }
        if (q > 2.0e-3) && (m == 2) {
            a = refine(kd, m, q, a);
        }
    } else {
        let ndiv = 10;
        let mut delta = (f64::from(m) - 3.0) * f64::from(m) / f64::from(ndiv);

        if (q - 3.0 * f64::from(m)) <= (f64::from(m.wrapping_mul(m)) - q) {
            let nn = ((q - 3.0 * f64::from(m)) / delta) as i32 + 1;
            delta = (q - 3.0 * f64::from(m)) / f64::from(nn);
            let mut q1 = 2.0 * f64::from(m);
            let mut a1 = cvqm(m, q1);
            let mut q2 = 3.0 * f64::from(m);
            let mut a2 = cvqm(m, q2);
            let mut qq = 3.0 * f64::from(m);
            for _i in 1..=nn {
                qq = qq + delta;
                a = (a1 * q2 - a2 * q1 + (a2 - a1) * qq) / (q2 - q1);
                a = refine(kd, m, qq, a);
                q1 = q2;
                q2 = qq;
                a1 = a2;
                a2 = a;
            }
        } else {
            // Unlike in xsf, a NaN `q` (with `m > 12`) results in NaN instead of 0: In C++ `(int)NaN` is UB
            // (INT_MIN on x86-64, so the loop is skipped), whereas `NaN as i32` is 0, so the loop runs once.
            let nn = ((f64::from(m.wrapping_mul(m)) - q) / delta) as i32 + 1;
            delta = (f64::from(m.wrapping_mul(m)) - q) / f64::from(nn);
            let mut q1 = f64::from(m) * (f64::from(m) - 1.0);
            let mut a1 = cvql(kd, m, q1);
            let mut q2 = f64::from(m.wrapping_mul(m));
            let mut a2 = cvql(kd, m, q2);
            let mut qq = f64::from(m.wrapping_mul(m));
            for _i in 1..=nn {
                qq = qq - delta;
                a = (a1 * q2 - a2 * q1 + (a2 - a1) * qq) / (q2 - q1);
                a = refine(kd, m, qq, a);
                q1 = q2;
                q2 = qq;
                a1 = a2;
                a2 = a;
            }
        }
    }
    a
}

#[allow(clippy::let_and_return, clippy::many_single_char_names)]
pub(crate) fn cvf(kd: i32, m: i32, q: f64, a: f64, mj: i32) -> f64 {
    // ======================================================
    // Purpose: Compute the value of F for characteristic
    //          equation of Mathieu functions
    // Input :  m --- Order of Mathieu functions
    //          q --- Parameter of Mathieu functions
    //          A --- Characteristic value
    // Output:  F --- Value of F for characteristic equation
    // ======================================================

    let ic = m / 2;
    let mut l = 0;
    let mut l0 = 0;
    let mut j0 = 2;
    let mut jf = ic;
    let mut t0 = 0.0;
    let mut t1 = 0.0;
    let mut t2 = 0.0;
    let b = a;

    if kd == 1 {
        j0 = 3;
        l0 = 2;
    }
    if (kd == 2) || (kd == 3) {
        l = 1;
    }
    if kd == 4 {
        jf = ic - 1;
    }

    for j in (ic + 1..=mj).rev() {
        // `pow(2.0 * j + l, 2)`
        let p = 2.0 * f64::from(j) + f64::from(l);
        t1 = -q * q / (p * p - b + t1);
    }

    if m <= 2 {
        if (kd == 1) && (m == 0) {
            t1 += t1;
        }
        if (kd == 1) && (m == 2) {
            t1 = -2.0 * q * q / (4.0 - b + t1) - 4.0;
        }
        if (kd == 2) && (m == 1) {
            t1 += q;
        }
        if (kd == 3) && (m == 1) {
            t1 -= q;
        }
    } else {
        if kd == 1 {
            t0 = 4.0 - b + 2.0 * q * q / b;
        }
        if kd == 2 {
            t0 = 1.0 - b + q;
        }
        if kd == 3 {
            t0 = 1.0 - b - q;
        }
        if kd == 4 {
            t0 = 4.0 - b;
        }
        t2 = -q * q / t0;
        for j in j0..=jf {
            // `pow(2.0 * j - l - l0, 2.0)`
            let p = 2.0 * f64::from(j) - f64::from(l) - f64::from(l0);
            t2 = -q * q / (p * p - b + t2);
        }
    }
    // `pow(2.0 * ic + l, 2)`
    let p = 2.0 * f64::from(ic) + f64::from(l);
    let f = p * p + t1 + t2 - b;
    f
}

#[allow(
    clippy::let_and_return,
    clippy::similar_names,
    clippy::assign_op_pattern
)]
pub(crate) fn cvql(kd: i32, m: i32, q: f64) -> f64 {
    // ========================================================
    // Purpose: Compute the characteristic value of Mathieu
    //          functions  for q ≥ 3m
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    // Output:  A0 --- Initial characteristic value
    // ========================================================

    let mut w = 0.0;
    if (kd == 1) || (kd == 2) {
        w = 2.0 * f64::from(m) + 1.0;
    }
    if (kd == 3) || (kd == 4) {
        w = 2.0 * f64::from(m) - 1.0;
    }
    let w2 = w * w;
    let w3 = w * w2;
    let w4 = w2 * w2;
    let w6 = w2 * w4;
    let d1 = 5.0 + 34.0 / w2 + 9.0 / w4;
    let d2 = (33.0 + 410.0 / w2 + 405.0 / w4) / w;
    let d3 = (63.0 + 1260.0 / w2 + 2943.0 / w4 + 486.0 / w6) / w2;
    let d4 = (527.0 + 15617.0 / w2 + 69001.0 / w4 + 41607.0 / w6) / w3;
    let c1 = 128.0;
    let p2 = q / w4;
    let p1 = p2.sqrt();
    let cv1 = -2.0 * q + 2.0 * w * q.sqrt() - (w2 + 1.0) / 8.0;
    let mut cv2 = (w + 3.0 / w) + d1 / (32.0 * p1) + d2 / (8.0 * c1 * p2);
    cv2 = cv2 + d3 / (64.0 * c1 * p1 * p2) + d4 / (16.0 * c1 * c1 * p2 * p2);
    let a0 = cv1 - cv2 / (c1 * p1);
    a0
}

#[allow(clippy::let_and_return)]
pub(crate) fn cvqm(m: i32, q: f64) -> f64 {
    // =====================================================
    // Purpose: Compute the characteristic value of Mathieu
    //          functions for q ≤ m*m
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    // Output:  A0 --- Initial characteristic value
    // =====================================================

    // NOTE: `m * m` is an `int` product here, which overflows for `m > 46340` (UB in C++, wrapping on x86-64)
    let hm1 = 0.5 * q / (f64::from(m.wrapping_mul(m)) - 1.0);
    let hm3 = 0.25 * hm1.powf(3.0) / (f64::from(m.wrapping_mul(m)) - 4.0);
    let hm5 = hm1 * hm3 * q
        / ((f64::from(m.wrapping_mul(m)) - 1.0) * (f64::from(m.wrapping_mul(m)) - 9.0));
    let a0 = f64::from(m.wrapping_mul(m))
        + q * (hm1
            + (5.0 * f64::from(m) * f64::from(m) + 7.0) * hm3
            + (9.0 * f64::from(m).powf(4.0) + 58.0 * f64::from(m) * f64::from(m) + 29.0) * hm5);
    a0
}

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::assign_op_pattern,
    clippy::needless_range_loop,
    clippy::needless_late_init,
    clippy::redundant_else,
    clippy::range_minus_one,
    clippy::if_not_else,
    unused_assignments
)]
pub(crate) fn fcoef(kd: i32, m: i32, q: f64, a: f64, fc: &mut [f64]) {
    // =====================================================
    // Purpose: Compute expansion coefficients for Mathieu
    //          functions and modified Mathieu functions
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    //          KD --- Case code
    //                 KD=1 for cem(x,q)  ( m = 0,2,4,...)
    //                 KD=2 for cem(x,q)  ( m = 1,3,5,...)
    //                 KD=3 for sem(x,q)  ( m = 1,3,5,...)
    //                 KD=4 for sem(x,q)  ( m = 2,4,6,...)
    //          A  --- Characteristic value of Mathieu
    //                 functions for given m and q
    // Output:  FC(k) --- Expansion coefficients of Mathieu
    //                 functions ( k= 1,2,...,KM )
    //                 FC(1),FC(2),FC(3),... correspond to
    //                 A0,A2,A4,... for KD=1 case, A1,A3,
    //                 A5,... for KD=2 case, B1,B3,B5,...
    //                 for KD=3 case and B2,B4,B6,... for
    //                 KD=4 case
    // =====================================================

    let mut jm = 0;
    let mut kb: i32;
    let mut f1: f64;
    let qm: f64;
    let mut s: f64;
    let mut f: f64;
    let mut u: f64;
    let mut v: f64;
    let mut f2: f64;
    let f3: f64;
    let mut sp: f64;
    let ss: f64;
    let s0: f64;

    for i in 0..251 {
        fc[i] = 0.0;
    }

    if q.abs() <= 1.0e-7 {
        // Expansion up to order Q^1 (Abramowitz & Stegun 20.2.27-28)
        if kd == 1 {
            jm = m / 2 + 1;
        } else if (kd == 2) || (kd == 3) {
            jm = (m - 1) / 2 + 1;
        } else if kd == 4 {
            jm = m / 2;
        }

        if jm + 1 > 251 {
            let fnan = f64::NAN;
            for i in 0..251 {
                fc[i] = fnan;
            }
            return;
        }
        // Proceed using the simplest expansion
        // NOTE: the `fc[jm - 2]` writes below index out of bounds (C++ UB; panics here) for
        // invalid (kd, m) combinations, e.g. kd == 4 with m == 3, or a negative m.
        if kd == 1 || kd == 2 {
            if m == 0 {
                fc[0] = 1.0 / 2.0_f64.sqrt();
                fc[1] = -q / (2.0 * 2.0_f64.sqrt());
            } else if m == 1 {
                fc[0] = 1.0;
                fc[1] = -q / 8.0;
            } else if m == 2 {
                fc[0] = q / 4.0;
                fc[1] = 1.0;
                fc[2] = -q / 12.0;
            } else {
                fc[(jm - 1) as usize] = 1.0;
                fc[jm as usize] = -q / (4.0 * f64::from(m + 1));
                fc[(jm - 2) as usize] = q / (4.0 * f64::from(m - 1));
            }
        } else if kd == 3 || kd == 4 {
            if m == 1 {
                fc[0] = 1.0;
                fc[1] = -q / 8.0;
            } else if m == 2 {
                fc[0] = 1.0;
                fc[1] = -q / 12.0;
            } else {
                fc[(jm - 1) as usize] = 1.0;
                fc[jm as usize] = -q / (4.0 * f64::from(m + 1));
                fc[(jm - 2) as usize] = q / (4.0 * f64::from(m - 1));
            }
        }
        return;
    } else if q <= 1.0 {
        qm = 7.5 + 56.1 * q.sqrt() - 134.7 * q + 90.7 * q.sqrt() * q;
    } else {
        qm = 17.0 + 3.1 * q.sqrt() - 0.126 * q + 0.0037 * q.sqrt() * q;
    }

    // NOTE: if `qm` is NaN (q < -1e-7 or q NaN) or +inf (q == +inf), the C++ `(int)` cast is UB
    // (x86-64 yields INT_MIN, after which `fc[km - 1]` is out of bounds); `as i32` gives 0
    // (-> out-of-bounds panic below) resp. i32::MAX (-> the NaN branch below). For finite
    // q > 1e-7, qm >= 7.5.
    let km = (qm + 0.5 * f64::from(m)) as i32;

    if km > 251 {
        // Overflow, generate NaNs
        for i in 0..251 {
            fc[i] = f64::NAN;
        }
        return;
    }

    kb = 0;
    s = 0.0;
    f = 1.0e-100;
    u = 0.0;
    fc[(km - 1) as usize] = 0.0;
    f2 = 0.0;

    if kd == 1 {
        for k in (3..=km).rev() {
            v = u;
            u = f;
            f = (a - 4.0 * f64::from(k) * f64::from(k)) * u / q - v;

            // NOTE: reads `fc[251]` when km == 251 (the first iteration): out of bounds for the
            // 251-element buffers that the C++ callers (`mtu0`, `mtu12`) pass (C++ UB; panics here
            // if `fc.len() == 251`). For a 252-element `fc`, `fc[251]` is not zeroed above.
            if f.abs() < fc[k as usize].abs() {
                kb = k;
                fc[0] = 1.0e-100;
                sp = 0.0;
                f3 = fc[k as usize];
                fc[1] = a / q * fc[0];
                fc[2] = (a - 4.0) * fc[1] / q - 2.0 * fc[0];
                u = fc[1];
                f1 = fc[2];

                for i in 3..=kb {
                    v = u;
                    u = f1;
                    f1 = (a - 4.0 * (f64::from(i) - 1.0) * (f64::from(i) - 1.0)) * u / q - v;
                    fc[i as usize] = f1;

                    if i == kb {
                        f2 = f1;
                    }
                    if i != kb {
                        sp += f1 * f1;
                    }
                }

                sp += 2.0 * fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2];
                // (unlike in xsf) see `sum_sq_rel`
                ss = if sp.is_finite() {
                    s + sp * (f3 / f2) * (f3 / f2)
                } else {
                    let sp_rel = sum_sq_rel(&fc[..kb as usize], f2) + (fc[0] / f2) * (fc[0] / f2);
                    s + sp_rel * f3 * f3
                };
                s0 = (1.0 / ss).sqrt();

                for j in 1..=km {
                    if j <= kb + 1 {
                        fc[(j - 1) as usize] = s0 * fc[(j - 1) as usize] * f3 / f2;
                    } else {
                        fc[(j - 1) as usize] *= s0;
                    }
                }
                if fc[0] < 0.0 {
                    for j in 0..km {
                        fc[j as usize] = -fc[j as usize];
                    }
                }
                return;
            } else {
                fc[(k - 1) as usize] = f;
                s += f * f;
            }
        }

        fc[1] = q * fc[2] / (a - 4.0 - 2.0 * q * q / a);
        fc[0] = q / a * fc[1];
        s += 2.0 * fc[0] * fc[0] + fc[1] * fc[1];
        s0 = (1.0 / s).sqrt();

        for k in 1..=km {
            fc[(k - 1) as usize] *= s0;
        }
    } else if (kd == 2) || (kd == 3) {
        'l45: {
            for k in (3..=km).rev() {
                v = u;
                u = f;
                f = (a - (2.0 * f64::from(k) - 1.0) * (2.0 * f64::from(k) - 1.0)) * u / q - v;

                if f.abs() >= fc[(k - 1) as usize].abs() {
                    fc[(k - 2) as usize] = f;
                    s += f * f;
                } else {
                    kb = k;
                    f3 = fc[(k - 1) as usize];
                    break 'l45;
                }
            }

            fc[0] = q / (a - 1.0 - (-1.0_f64).powf(f64::from(kd)) * q) * fc[1];
            s += fc[0] * fc[0];
            s0 = (1.0 / s).sqrt();

            for k in 1..=km {
                fc[(k - 1) as usize] *= s0;
            }
            if fc[0] < 0.0 {
                for j in 0..km {
                    fc[j as usize] = -fc[j as usize];
                }
            }
            return;
        }
        // L45:
        fc[0] = 1.0e-100;
        fc[1] = (a - 1.0 - (-1.0_f64).powf(f64::from(kd)) * q) / q * fc[0];
        sp = 0.0;
        u = fc[0];
        f1 = fc[1];

        for i in 2..=kb - 1 {
            v = u;
            u = f1;
            f1 = (a - (2.0 * f64::from(i) - 1.0) * (2.0 * f64::from(i) - 1.0)) * u / q - v;

            if i != kb - 1 {
                fc[i as usize] = f1;
                sp += f1 * f1;
            } else {
                f2 = f1;
            }
        }

        sp += fc[0] * fc[0] + fc[1] * fc[1];
        // (unlike in xsf) see `sum_sq_rel`
        ss = if sp.is_finite() {
            s + sp * (f3 / f2) * (f3 / f2)
        } else {
            s + sum_sq_rel(&fc[..(kb - 1) as usize], f2) * f3 * f3
        };
        s0 = (1.0 / ss).sqrt();

        for j in 1..=km {
            if j < kb {
                fc[(j - 1) as usize] *= s0 * f3 / f2;
            }

            if j >= kb {
                fc[(j - 1) as usize] *= s0;
            }
        }
    } else if kd == 4 {
        'l70: {
            for k in (3..=km).rev() {
                v = u;
                u = f;
                f = (a - 4.0 * f64::from(k) * f64::from(k)) * u / q - v;

                // Unlike in xsf, this compares against `fc[k - 1]` instead of `fc[k]`, like the
                // Fortran original (FC(K)) and the KD=2,3 case do (scipy/xsf#343)
                if f.abs() >= fc[(k - 1) as usize].abs() {
                    fc[(k - 2) as usize] = f;
                    s += f * f;
                } else {
                    kb = k;
                    f3 = fc[(k - 1) as usize];
                    break 'l70;
                }
            }

            fc[0] = q / (a - 4.0) * fc[1];
            s += fc[0] * fc[0];
            s0 = (1.0 / s).sqrt();

            for k in 1..=km {
                fc[(k - 1) as usize] *= s0;
            }
            if fc[0] < 0.0 {
                for j in 0..km {
                    fc[j as usize] = -fc[j as usize];
                }
            }
            return;
        }
        // L70:
        fc[0] = 1.0e-100;
        fc[1] = (a - 4.0) / q * fc[0];
        sp = 0.0;
        u = fc[0];
        f1 = fc[1];

        for i in 2..=kb - 1 {
            v = u;
            u = f1;
            f1 = (a - 4.0 * f64::from(i) * f64::from(i)) * u / q - v;

            if i != kb - 1 {
                fc[i as usize] = f1;
                sp = sp + f1 * f1;
            } else {
                f2 = f1;
            }
        }

        sp += fc[0] * fc[0] + fc[1] * fc[1];
        // (unlike in xsf) see `sum_sq_rel`
        ss = if sp.is_finite() {
            s + sp * (f3 / f2) * (f3 / f2)
        } else {
            s + sum_sq_rel(&fc[..(kb - 1) as usize], f2) * f3 * f3
        };
        s0 = (1.0 / ss).sqrt();

        for j in 1..=km {
            if j < kb {
                fc[(j - 1) as usize] *= s0 * f3 / f2;
            } else {
                fc[(j - 1) as usize] *= s0;
            }
        }
    }
    if fc[0] < 0.0 {
        for j in 0..km {
            fc[j as usize] = -fc[j as usize];
        }
    }
}

/// `sum((fc / f2)^2)`, for computing `sp * (f3 / f2)^2` in `fcoef` without overflowing `sp`.
///
/// Unlike in xsf, `fcoef` uses this if the sum of squares `sp` of the forward recurrence overflows
/// (e.g. for m >= 92 and q = 0.01), instead of returning zeros (scipy/xsf#344).
fn sum_sq_rel(fc: &[f64], f2: f64) -> f64 {
    fc.iter().map(|c| (c / f2) * (c / f2)).sum()
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::too_many_lines,
    clippy::unreadable_literal
)]
pub(crate) fn mtu0(kf: i32, m: i32, q: f64, x: f64, csf: &mut f64, csd: &mut f64) -> Status {
    // ===============================================================
    // Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
    //          and their derivatives ( q ≥ 0 )
    // Input :  KF  --- Function code
    //                  KF=1 for computing cem(x,q) and cem'(x,q)
    //                  KF=2 for computing sem(x,q) and sem'(x,q)
    //          m   --- Order of Mathieu functions
    //          q   --- Parameter of Mathieu functions
    //          x   --- Argument of Mathieu functions (in degrees)
    // Output:  CSF --- cem(x,q) or sem(x,q)
    //          CSD --- cem'x,q) or sem'x,q)
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed. The output
    //              values will be set to nan.
    //          Status::Other
    //              An internal check failed. For mtu0, this occurs
    //              when km (see the code) is too big.  km is a
    //              function of q and m.  The output values will be
    //              set to nan.
    //
    // Routines called:
    //      (1) CVA2 for computing the characteristic values
    //      (2) FCOEF for computing the expansion coefficients
    // ===============================================================

    let mut kd = 0;
    let eps = 1.0e-14;
    let rd = 1.74532925199433e-2;

    if kf == 1 && m == 2 * (m / 2) {
        kd = 1;
    }
    if kf == 1 && m != 2 * (m / 2) {
        kd = 2;
    }
    if kf == 2 && m != 2 * (m / 2) {
        kd = 3;
    }
    if kf == 2 && m == 2 * (m / 2) {
        kd = 4;
    }

    let a = cva2(kd, m, q);

    let qm = if q <= 1.0 {
        7.5 + 56.1 * q.sqrt() - 134.7 * q + 90.7 * q.sqrt() * q
    } else {
        17.0 + 3.1 * q.sqrt() - 0.126 * q + 0.0037 * q.sqrt() * q
    };

    // NOTE: `(int)` is UB in C++ if `qm + 0.5 * m` is out of the int range or NaN (e.g. q >~ 7e7, or
    // q NaN). Rust's `as i32` saturates (NaN -> 0), so a huge qm gives km = i32::MAX > 251 (returns
    // Status::Other with NaNs), whereas g++ on x86-64 gives km = INT_MIN (no terms summed, so
    // csf = csd = 0 with Status::OK).
    let km = (qm + 0.5 * f64::from(m)) as i32;

    if km > 251 {
        *csf = f64::NAN;
        *csd = f64::NAN;
        return Status::Other;
    }

    // (unlike in xsf) 252, since fcoef reads fc[km] and km can be 251
    let Some(mut fg) = try_zeroed(252) else {
        *csf = f64::NAN;
        *csd = f64::NAN;
        return Status::NoMemory;
    };
    fcoef(kd, m, q, a, &mut fg);

    let ic = m / 2 + 1;
    let xr = x * rd;
    *csf = 0.0;
    for k in 1..=km {
        if kd == 1 {
            *csf += fg[(k - 1) as usize] * (f64::from(2 * k - 2) * xr).cos();
        } else if kd == 2 {
            *csf += fg[(k - 1) as usize] * (f64::from(2 * k - 1) * xr).cos();
        } else if kd == 3 {
            *csf += fg[(k - 1) as usize] * (f64::from(2 * k - 1) * xr).sin();
        } else if kd == 4 {
            *csf += fg[(k - 1) as usize] * (f64::from(2 * k) * xr).sin();
        }
        // NOTE: this reads fg[k] (not fg[k - 1] as in the csd loop below); for k == km == 251 that
        // is fg[251], out of bounds of the 251-element buffer in C++.
        if (k >= ic) && (fg[k as usize].abs() < csf.abs() * eps) {
            break;
        }
    }

    *csd = 0.0;
    for k in 1..=km {
        if kd == 1 {
            *csd -= f64::from(2 * k - 2) * fg[(k - 1) as usize] * (f64::from(2 * k - 2) * xr).sin();
        } else if kd == 2 {
            *csd -= f64::from(2 * k - 1) * fg[(k - 1) as usize] * (f64::from(2 * k - 1) * xr).sin();
        } else if kd == 3 {
            *csd += f64::from(2 * k - 1) * fg[(k - 1) as usize] * (f64::from(2 * k - 1) * xr).cos();
        } else if kd == 4 {
            *csd += 2.0 * f64::from(k) * fg[(k - 1) as usize] * (f64::from(2 * k) * xr).cos();
        }
        if (k >= ic) && (fg[(k - 1) as usize].abs() < csd.abs() * eps) {
            break;
        }
    }
    Status::Ok
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
pub(crate) fn mtu12(
    kf: i32,
    kc: i32,
    m: i32,
    q: f64,
    x: f64,
    f1r: &mut f64,
    d1r: &mut f64,
    f2r: &mut f64,
    d2r: &mut f64,
) -> Status {
    // ==============================================================
    // Purpose: Compute modified Mathieu functions of the first and
    //          second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
    //          and their derivatives
    // Input:   KF --- Function code
    //                 KF=1 for computing Mcm(x,q)
    //                 KF=2 for computing Msm(x,q)
    //          KC --- Function Code
    //                 KC=1 for computing the first kind
    //                 KC=2 for computing the second kind
    //                      or Msm(2)(x,q) and Msm(2)'(x,q)
    //                 KC=3 for computing both the first
    //                      and second kinds
    //          m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions ( q ≥ 0 )
    //          x  --- Argument of Mathieu functions
    // Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
    //          D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
    //          F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
    //          D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed. The output
    //              values will be set to nan.
    //          Status::Other
    //              An internal check failed. For mtu12, this occurs
    //              when km (see the code) is too big.  km is a
    //              function of q and m.  The output values will be
    //              set to nan.
    // Routines called:
    //      (1) CVA2 for computing the characteristic values
    //      (2) FCOEF for computing expansion coefficients
    //      (3) JYNB for computing Jn(x), Yn(x) and their
    //          derivatives
    // ==============================================================

    let eps = 1.0e-14;
    // NOTE: in C++ `kd` is left uninitialized (UB to read) if kf is neither 1 nor 2; here it's 0.
    let mut kd = 0;
    let mut nm = 0;

    if (kf == 1) && (m % 2 == 0) {
        kd = 1;
    }
    if (kf == 1) && (m % 2 != 0) {
        kd = 2;
    }
    if (kf == 2) && (m % 2 != 0) {
        kd = 3;
    }
    if (kf == 2) && (m % 2 == 0) {
        kd = 4;
    }

    let a = cva2(kd, m, q);

    let qm = if q <= 1.0 {
        7.5 + 56.1 * q.sqrt() - 134.7 * q + 90.7 * q.sqrt() * q
    } else {
        17.0 + 3.1 * q.sqrt() - 0.126 * q + 0.0037 * q.sqrt() * q
    };

    // NOTE: `(int)` is UB in C++ if `qm + 0.5 * m` is out of the int range or NaN (e.g. q >~ 7e7, or
    // q NaN). Rust's `as i32` saturates (NaN -> 0), so a huge qm gives km = i32::MAX >= 251 (returns
    // Status::Other with NaNs), whereas g++ on x86-64 gives km = INT_MIN (then jynb(INT_MIN + 1, ...)).
    let km = (qm + 0.5 * f64::from(m)) as i32;
    if km >= 251 {
        *f1r = f64::NAN;
        *d1r = f64::NAN;
        *f2r = f64::NAN;
        *d2r = f64::NAN;
        return Status::Other;
    }

    // allocate memory after a possible NAN return
    let fg = try_zeroed(251);
    let bj1 = try_zeroed(252);
    let dj1 = try_zeroed(252);
    let bj2 = try_zeroed(252);
    let dj2 = try_zeroed(252);
    let by1 = try_zeroed(252);
    let dy1 = try_zeroed(252);
    let by2 = try_zeroed(252);
    let dy2 = try_zeroed(252);

    let (
        Some(mut fg),
        Some(mut bj1),
        Some(mut dj1),
        Some(mut bj2),
        Some(mut dj2),
        Some(mut by1),
        Some(mut dy1),
        Some(mut by2),
        Some(mut dy2),
    ) = (fg, bj1, dj1, bj2, dj2, by1, dy1, by2, dy2)
    else {
        *f1r = f64::NAN;
        *d1r = f64::NAN;
        *f2r = f64::NAN;
        *d2r = f64::NAN;
        return Status::NoMemory;
    };

    fcoef(kd, m, q, a, &mut fg);
    let mut ic = m / 2 + 1;
    if kd == 4 {
        ic = m / 2;
    }

    let c1 = (-x).exp();
    let c2 = x.exp();
    let u1 = q.sqrt() * c1;
    let u2 = q.sqrt() * c2;
    jynb(km + 1, u1, &mut nm, &mut bj1, &mut dj1, &mut by1, &mut dy1);
    jynb(km + 1, u2, &mut nm, &mut bj2, &mut dj2, &mut by2, &mut dy2);
    let mut w1 = 0.0;
    let mut w2 = 0.0;

    if kc != 2 {
        *f1r = 0.0;
        for k in 1..=km {
            let k_ = k as usize;
            if kd == 1 {
                *f1r += (-1.0_f64).powf(f64::from(ic + k)) * fg[k_ - 1] * bj1[k_ - 1] * bj2[k_ - 1];
            } else if kd == 2 || kd == 3 {
                *f1r += (-1.0_f64).powf(f64::from(ic + k))
                    * fg[k_ - 1]
                    * (bj1[k_ - 1] * bj2[k_]
                        + (-1.0_f64).powf(f64::from(kd)) * bj1[k_] * bj2[k_ - 1]);
            } else {
                *f1r += (-1.0_f64).powf(f64::from(ic + k))
                    * fg[k_ - 1]
                    * (bj1[k_ - 1] * bj2[k_ + 1] - bj1[k_ + 1] * bj2[k_ - 1]);
            }

            if k >= 5 && (*f1r - w1).abs() < f1r.abs() * eps {
                break;
            }
            w1 = *f1r;
        }

        *f1r /= fg[0];

        *d1r = 0.0;
        for k in 1..=km {
            let k_ = k as usize;
            if kd == 1 {
                *d1r += (-1.0_f64).powf(f64::from(ic + k))
                    * fg[k_ - 1]
                    * (c2 * bj1[k_ - 1] * dj2[k_ - 1] - c1 * dj1[k_ - 1] * bj2[k_ - 1]);
            } else if kd == 2 || kd == 3 {
                *d1r += (-1.0_f64).powf(f64::from(ic + k))
                    * fg[k_ - 1]
                    * (c2
                        * (bj1[k_ - 1] * dj2[k_]
                            + (-1.0_f64).powf(f64::from(kd)) * bj1[k_] * dj2[k_ - 1])
                        - c1 * (dj1[k_ - 1] * bj2[k_]
                            + (-1.0_f64).powf(f64::from(kd)) * dj1[k_] * bj2[k_ - 1]));
            } else {
                *d1r += (-1.0_f64).powf(f64::from(ic + k))
                    * fg[k_ - 1]
                    * (c2 * (bj1[k_ - 1] * dj2[k_ + 1] - bj1[k_ + 1] * dj2[k_ - 1])
                        - c1 * (dj1[k_ - 1] * bj2[k_ + 1] - dj1[k_ + 1] * bj2[k_ - 1]));
            }

            if k >= 5 && (*d1r - w2).abs() < d1r.abs() * eps {
                break;
            }
            w2 = *d1r;
        }
        *d1r *= q.sqrt() / fg[0];
        if kc == 1 {
            return Status::Ok;
        }
    }

    *f2r = 0.0;
    for k in 1..=km {
        let k_ = k as usize;
        if kd == 1 {
            *f2r += (-1.0_f64).powf(f64::from(ic + k)) * fg[k_ - 1] * bj1[k_ - 1] * by2[k_ - 1];
        } else if kd == 2 || kd == 3 {
            *f2r += (-1.0_f64).powf(f64::from(ic + k))
                * fg[k_ - 1]
                * (bj1[k_ - 1] * by2[k_] + (-1.0_f64).powf(f64::from(kd)) * bj1[k_] * by2[k_ - 1]);
        } else {
            *f2r += (-1.0_f64).powf(f64::from(ic + k))
                * fg[k_ - 1]
                * (bj1[k_ - 1] * by2[k_ + 1] - bj1[k_ + 1] * by2[k_ - 1]);
        }

        if k >= 5 && (*f2r - w1).abs() < f2r.abs() * eps {
            break;
        }
        w1 = *f2r;
    }
    *f2r /= fg[0];

    *d2r = 0.0;
    for k in 1..=km {
        let k_ = k as usize;
        if kd == 1 {
            *d2r += (-1.0_f64).powf(f64::from(ic + k))
                * fg[k_ - 1]
                * (c2 * bj1[k_ - 1] * dy2[k_ - 1] - c1 * dj1[k_ - 1] * by2[k_ - 1]);
        } else if kd == 2 || kd == 3 {
            *d2r += (-1.0_f64).powf(f64::from(ic + k))
                * fg[k_ - 1]
                * (c2
                    * (bj1[k_ - 1] * dy2[k_]
                        + (-1.0_f64).powf(f64::from(kd)) * bj1[k_] * dy2[k_ - 1])
                    - c1 * (dj1[k_ - 1] * by2[k_]
                        + (-1.0_f64).powf(f64::from(kd)) * dj1[k_] * by2[k_ - 1]));
        } else {
            *d2r += (-1.0_f64).powf(f64::from(ic + k))
                * fg[k_ - 1]
                * (c2 * (bj1[k_ - 1] * dy2[k_ + 1] - bj1[k_ + 1] * dy2[k_ - 1])
                    - c1 * (dj1[k_ - 1] * by2[k_ + 1] - dj1[k_ + 1] * by2[k_ - 1]));
        }

        if k >= 5 && (*d2r - w2).abs() < d2r.abs() * eps {
            break;
        }
        w2 = *d2r;
    }
    *d2r = *d2r * q.sqrt() / fg[0];
    Status::Ok
}

#[allow(clippy::float_cmp, clippy::many_single_char_names)]
pub(crate) fn refine(kd: i32, m: i32, q: f64, a: f64) -> f64 {
    // =====================================================
    // Purpose: calculate the accurate characteristic value
    //          by the secant method
    // Input :  m --- Order of Mathieu functions
    //          q --- Parameter of Mathieu functions
    //          A --- Initial characteristic value
    // Output:  A --- Refineed characteristic value
    // Routine called:  CVF for computing the value of F for
    //                  characteristic equation
    // ========================================================

    const EPS: f64 = 1e-14;

    // NOTE: the `int` sums `10 + m` and `mj++` overflow for `m > INT_MAX - 110`, which is UB in C++
    // (wrapping on x86-64)
    let mut mj = m.wrapping_add(10);
    let mut x0 = a;
    let mut f0 = cvf(kd, m, q, x0, mj);
    let mut x1 = 1.002 * a;
    let mut f1 = cvf(kd, m, q, x1, mj);
    // `x` is uninitialized in C++, but the loop body always runs, so it is always assigned before it's returned
    let mut x = 0.0;
    for _it in 1..=100 {
        mj = mj.wrapping_add(1);
        x = x1 - (x1 - x0) / (1.0 - f0 / f1);
        let f = cvf(kd, m, q, x, mj);
        if ((1.0 - x1 / x).abs() < EPS) || (f == 0.0) {
            break;
        }
        x0 = x1;
        f0 = f1;
        x1 = x;
        f1 = f;
    }
    x
}
