//! Bessel functions of integer order, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.

use core::f64::consts::PI;

#[allow(
    clippy::approx_constant,
    clippy::excessive_precision,
    clippy::unreadable_literal
)]
pub(crate) fn gam0(x: f64) -> f64 {
    // ================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 )
    // Output:  GA --- Г(x)
    // ================================================
    const G: [f64; 25] = [
        1.0e0,
        0.5772156649015329e0,
        -0.6558780715202538e0,
        -0.420026350340952e-1,
        0.1665386113822915e0,
        -0.421977345555443e-1,
        -0.96219715278770e-2,
        0.72189432466630e-2,
        -0.11651675918591e-2,
        -0.2152416741149e-3,
        0.1280502823882e-3,
        -0.201348547807e-4,
        -0.12504934821e-5,
        0.11330272320e-5,
        -0.2056338417e-6,
        0.61160950e-8,
        0.50020075e-8,
        -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
        -0.36968e-11,
        0.51e-12,
        -0.206e-13,
        -0.54e-14,
        0.14e-14,
    ];
    let mut gr = G[24];
    for k in (0..=23).rev() {
        gr = gr * x + G[k];
    }
    1.0 / (gr * x)
}

#[allow(clippy::cast_sign_loss)]
pub(crate) fn jynb(
    n: i32,
    x: f64,
    nm: &mut i32,
    bj: &mut [f64],
    dj: &mut [f64],
    by: &mut [f64],
    dy: &mut [f64],
) {
    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x) and
    //          their derivatives
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Order of Jn(x) and Yn(x)
    // Output:  BJ(n) --- Jn(x)
    //          DJ(n) --- Jn'(x)
    //          BY(n) --- Yn(x)
    //          DY(n) --- Yn'(x)
    //          NM --- Highest order computed
    // Routines called:
    //          JYNBH to calculate the Jn and Yn
    // =====================================================

    jynbh(n, 0, x, nm, bj, by);
    // Compute derivatives by differentiation formulas
    if x < 1.0e-100 {
        for k in 0..=n {
            dj[k as usize] = 0.0;
            dy[k as usize] = 1.0e+300;
        }
        dj[1] = 0.5;
    } else {
        dj[0] = -bj[1];
        for k in 1..=*nm {
            dj[k as usize] = bj[(k - 1) as usize] - f64::from(k) / x * bj[k as usize];
        }

        dy[0] = -by[1];
        for k in 1..=*nm {
            dy[k as usize] = by[(k - 1) as usize] - f64::from(k) * by[k as usize] / x;
        }
    }
}

#[allow(
    clippy::approx_constant,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::excessive_precision,
    clippy::many_single_char_names,
    clippy::needless_late_init,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::unreadable_literal
)]
pub(crate) fn jynbh(n: i32, nmin: i32, x: f64, nm: &mut i32, bj: &mut [f64], by: &mut [f64]) {
    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x)
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
    //          nmin -- Lowest order computed  ( nmin ≥ 0 )
    // Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
    //          BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 to calculate the starting
    //          point for backward recurrence
    // =====================================================

    let ky: i32;
    let r2p = 0.63661977236758;
    let mut bj0: f64;
    let mut bj1: f64;
    let mut by0: f64;
    let mut by1: f64;

    let a: [f64; 4] = [
        -0.0703125,
        0.112152099609375,
        -0.5725014209747314,
        0.6074042001273483e+01,
    ];
    let b: [f64; 4] = [
        0.0732421875,
        -0.2271080017089844,
        0.1727727502584457e+01,
        -0.2438052969955606e+02,
    ];
    let a1: [f64; 4] = [
        0.1171875,
        -0.144195556640625,
        0.6765925884246826,
        -0.6883914268109947e+01,
    ];
    let b1: [f64; 4] = [
        -0.1025390625,
        0.2775764465332031,
        -0.1993531733751297e+01,
        0.2724882731126854e+02,
    ];

    *nm = n;
    if x < 1.0e-100 {
        for k in nmin..=n {
            bj[(k - nmin) as usize] = 0.0;
            by[(k - nmin) as usize] = -1.0e+300;
        }

        if nmin == 0 {
            bj[0] = 1.0;
        }
        return;
    }

    // Unlike in xsf (scipy/xsf#351), the Hankel expansion is used for x >= 2^31 / 0.9 (~2.39e9),
    // since `as i32` saturates.
    if (x <= 300.0) || (n > (0.9 * x) as i32) {
        // Backward recurrence for Jn
        if n == 0 {
            *nm = 1;
        }
        let mut m = msta1(x, 200);
        if m < *nm {
            *nm = m;
        } else {
            m = msta2(x, *nm, 15);
        }
        let mut bs = 0.0;
        let mut su = 0.0;
        let mut sv = 0.0;
        let mut f2 = 0.0;
        let mut f1 = 1.0e-100;
        let mut f = 0.0;

        for k in (0..=m).rev() {
            f = 2.0 * (f64::from(k) + 1.0) / x * f1 - f2;
            if (k <= *nm) && (k >= nmin) {
                bj[(k - nmin) as usize] = f;
            }
            if k == 2 * (k / 2) && k != 0 {
                bs += 2.0 * f;
                su += (-1.0_f64).powf(f64::from(k / 2)) * f / f64::from(k);
            } else if k > 1 {
                // Unlike in xsf, `k * k` is computed in `double` instead of `int`, which overflows for
                // k > 46340 (e.g. for jyzo with n >~ 4e4); this is exact for k < 2^26
                sv += (-1.0_f64).powf(f64::from(k / 2)) * f64::from(k)
                    / (f64::from(k) * f64::from(k) - 1.0)
                    * f;
            }
            f2 = f1;
            f1 = f;
        }
        let s0 = bs + f;

        for k in nmin..=*nm {
            bj[(k - nmin) as usize] /= s0;
        }
        // Estimates for Yn at start of recurrence
        bj0 = f1 / s0;
        bj1 = f2 / s0;
        let ec = (x / 2.0).ln() + 0.5772156649015329;
        by0 = r2p * (ec * bj0 - 4.0 * su / s0);
        by1 = r2p * ((ec - 1.0) * bj1 - bj0 / x - 4.0 * sv / s0);

        if 0 >= nmin {
            by[(0 - nmin) as usize] = by0;
        }
        if 1 >= nmin {
            by[(1 - nmin) as usize] = by1;
        }
        ky = 2;
    } else {
        // Hankel expansion
        let t1 = x - 0.25 * PI;
        let mut p0 = 1.0;
        let mut q0 = -0.125 / x;

        for k in 1..=4 {
            p0 += a[(k - 1) as usize] * x.powf(f64::from(-2 * k));
            q0 += b[(k - 1) as usize] * x.powf(f64::from(-2 * k - 1));
        }

        let cu = (r2p / x).sqrt();
        bj0 = cu * (p0 * t1.cos() - q0 * t1.sin());
        by0 = cu * (p0 * t1.sin() + q0 * t1.cos());

        if 0 >= nmin {
            bj[(0 - nmin) as usize] = bj0;
            by[(0 - nmin) as usize] = by0;
        }

        let t2 = x - 0.75 * PI;
        let mut p1 = 1.0;
        let mut q1 = 0.375 / x;

        for k in 1..=4 {
            p1 += a1[(k - 1) as usize] * x.powf(f64::from(-2 * k));
            q1 += b1[(k - 1) as usize] * x.powf(f64::from(-2 * k - 1));
        }

        bj1 = cu * (p1 * t2.cos() - q1 * t2.sin());
        by1 = cu * (p1 * t2.sin() + q1 * t2.cos());

        if 1 >= nmin {
            bj[(1 - nmin) as usize] = bj1;
            by[(1 - nmin) as usize] = by1;
        }

        for k in 2..=*nm {
            let bjk = 2.0 * (f64::from(k) - 1.0) / x * bj1 - bj0;
            if k >= nmin {
                bj[(k - nmin) as usize] = bjk;
            }
            bj0 = bj1;
            bj1 = bjk;
        }
        ky = 2;
    }
    // Forward recurrence for Yn
    for k in ky..=*nm {
        let byk = 2.0 * (f64::from(k) - 1.0) * by1 / x - by0;

        if k >= nmin {
            by[(k - nmin) as usize] = byk;
        }

        by0 = by1;
        by1 = byk;
    }
}

#[allow(clippy::similar_names, clippy::too_many_arguments)]
pub(crate) fn jyndd(
    n: i32,
    x: f64,
    bjn: &mut f64,
    djn: &mut f64,
    fjn: &mut f64,
    byn: &mut f64,
    dyn_: &mut f64,
    fyn: &mut f64,
) {
    // ===========================================================
    // purpose: compute bessel functions jn(x) and yn(x), and
    //          their first and second derivatives
    // input:   x   ---  argument of jn(x) and yn(x) ( x > 0 )
    //          n   ---  order of jn(x) and yn(x)
    // output:  bjn ---  jn(x)
    //          djn ---  jn'(x)
    //          fjn ---  jn"(x)
    //          byn ---  yn(x)
    //          dyn ---  yn'(x)
    //          fyn ---  yn"(x)
    // routines called:
    //          jynbh to compute jn and yn
    // ===========================================================

    // NOTE: in C++ `bj` and `by` are uninitialized, and partly read uninitialized if
    // msta1(x, 200) <= n (x << n). jyzo only calls this near the zeros (x > n), so that's fine.
    let mut nm = 0;
    let mut bj = [0.0; 2];
    let mut by = [0.0; 2];

    jynbh(n + 1, n, x, &mut nm, &mut bj, &mut by);
    // compute derivatives by differentiation formulas
    *bjn = bj[0];
    *byn = by[0];
    *djn = -bj[1] + f64::from(n) * bj[0] / x;
    *dyn_ = -by[1] + f64::from(n) * by[0] / x;
    // Unlike in xsf, `n * n` is computed in `double` instead of `int`, which overflows for n > 46340
    // (e.g. for jyzo); this is exact for n < 2^26
    *fjn = (f64::from(n) * f64::from(n) / (x * x) - 1.0) * (*bjn) - (*djn) / x;
    *fyn = (f64::from(n) * f64::from(n) / (x * x) - 1.0) * (*byn) - (*dyn_) / x;
}

#[allow(
    non_snake_case,
    clippy::cast_sign_loss,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::unreadable_literal
)]
pub(crate) fn jyzo(
    n: i32,
    nt: i32,
    rj0: &mut [f64],
    rj1: &mut [f64],
    ry0: &mut [f64],
    ry1: &mut [f64],
) {
    // ======================================================
    // Purpose: Compute the zeros of Bessel functions Jn(x),
    //          Yn(x), and their derivatives
    // Input :  n  --- Order of Bessel functions  (n >= 0)
    //          NT --- Number of zeros (roots)
    // Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
    //          RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
    //          RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
    //          RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
    // Routine called: JYNDD for computing Jn(x), Yn(x), and
    //                 their first and second derivatives
    // ======================================================

    /*
     * SciPy Note:
     * See GH-18859 for additional changes done by SciPy for
     * better initial condition selection in Newton iteration
     */

    // NOTE: each of rj0, rj1, ry0, ry1 gets max(nt, 1) elements written (index 0 is written even
    // if nt < 1), and only rj0[0..nt] / ry0[0..nt] are read, so buffers of length max(nt, 1) are
    // in bounds. The caller (bessel_zeros) returns early for N == 0.
    // NOTE: like in C++, the Newton iterations loop forever if they don't converge.
    let mut L: i32;
    let mut x: f64;
    let mut bjn = 0.0;
    let mut djn = 0.0;
    let mut fjn = 0.0;
    let mut byn = 0.0;
    let mut dyn_ = 0.0;
    let mut fyn = 0.0;
    // -- Newton method for j_{N,L}
    // initial guess for j_{N,1}
    if n == 0 {
        x = 2.4;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 1.85576 * f64::from(n).powf(0.33333)
            + 1.03315 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l10: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= bjn / djn;
        if (x - x0).abs() > 1e-11 {
            continue 'l10;
        }

        L += 1;
        rj0[(L - 1) as usize] = x;
        // initial guess for j_{N,L+1}
        if L == 1 {
            if n == 0 {
                x = 5.52;
            } else {
                // Expansion from https://dlmf.nist.gov/10.21#E32 and
                // coefficients from Olver 1951
                x = f64::from(n)
                    + 3.24460 * f64::from(n).powf(0.33333)
                    + 3.15824 / f64::from(n).powf(0.33333);
            }
        } else {
            // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
            x = rj0[(L - 1) as usize] + (rj0[(L - 1) as usize] - rj0[(L - 2) as usize]);
        }
        // NOTE: `n + 10` overflows i32 for n > i32::MAX - 10 (C++ UB; Rust panics in debug builds).
        if L <= (n + 10) {
            jyndd(
                n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
            );
            let h = djn.abs().atan2((fjn * bjn).abs().sqrt());
            // NOTE: `atan(h)` here, but `tan(h)` in the corresponding step for y_{N,L} below;
            // reproduced as is.
            let b = -djn / (bjn * h.atan());
            x -= (h - PI / 2.0) / b;
        }

        if L < nt {
            continue 'l10;
        }
        break;
    }

    // -- Newton method for j_{N,L+1}'
    if n == 0 {
        x = 3.8317;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 0.80861 * f64::from(n).powf(0.33333)
            + 0.07249 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l15: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= djn / fjn;
        if (x - x0).abs() > 1e-11 {
            continue 'l15;
        }
        L += 1;
        rj1[(L - 1) as usize] = x;
        if L < nt {
            // https://dlmf.nist.gov/10.21#E20
            x = rj1[(L - 1) as usize] + (rj0[L as usize] - rj0[(L - 1) as usize]);
            continue 'l15;
        }
        break;
    }

    // -- Newton method for y_{N,L}
    // initial guess for y_{N,1}
    if n == 0 {
        x = 0.89357697;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 0.93158 * f64::from(n).powf(0.33333)
            + 0.26035 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l20: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= byn / dyn_;
        if (x - x0).abs() > 1.0e-11 {
            continue 'l20;
        }
        L += 1;
        ry0[(L - 1) as usize] = x;
        // initial guess for y_{N,L+1}
        if L == 1 {
            if n == 0 {
                x = 3.957678419314858;
            } else {
                // Expansion from https://dlmf.nist.gov/10.21#E33 and
                // coefficients from Olver 1951
                x = f64::from(n)
                    + 2.59626 * f64::from(n).powf(0.33333)
                    + 2.022183 / f64::from(n).powf(0.33333);
            }
        } else {
            // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
            x = ry0[(L - 1) as usize] + (ry0[(L - 1) as usize] - ry0[(L - 2) as usize]);
        }
        if L <= n + 10 {
            jyndd(
                n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
            );
            let h = dyn_.abs().atan2((fyn * byn).abs().sqrt());
            let b = -dyn_ / (byn * h.tan());
            x -= (h - PI / 2.0) / b;
        }

        if L < nt {
            continue 'l20;
        }
        break;
    }

    // -- Newton method for y_{N,L+1}'
    if n == 0 {
        x = 2.67257;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 1.8211 * f64::from(n).powf(0.33333)
            + 0.94001 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l25: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= dyn_ / fyn;
        if (x - x0).abs() > 1.0e-11 {
            continue 'l25;
        }
        L += 1;
        ry1[(L - 1) as usize] = x;
        if L < nt {
            // https://dlmf.nist.gov/10.21#E20
            x = ry1[(L - 1) as usize] + (ry0[L as usize] - ry0[(L - 1) as usize]);
            continue 'l25;
        }
        break;
    }
}

#[allow(
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::similar_names
)]
pub(crate) fn lamn(n: i32, x: f64, nm: &mut i32, bl: &mut [f64], dl: &mut [f64]) {
    // =========================================================
    // Purpose: Compute lambda functions and their derivatives
    // Input:   x --- Argument of lambda function
    //          n --- Order of lambda function
    // Output:  BL(n) --- Lambda function of order n
    //          DL(n) --- Derivative of lambda function
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the start
    //          point for backward recurrence
    // =========================================================

    // NOTE: bl[0..=n] and dl[0..=n] can be written, and also dl[1] (for |x| < 1e-100) or bl[1] /
    // dl[1] (for x > 12, where *nm is 1 if n == 0), which is out of bounds for n == 0 with buffers
    // of length n + 1. The caller (jahnke_emden_lambda) passes n >= 1, so that's in bounds.
    *nm = n;
    if x.abs() < 1e-100 {
        for k in 0..=n {
            bl[k as usize] = 0.0;
            dl[k as usize] = 0.0;
        }
        bl[0] = 1.0;
        dl[1] = 0.5;
        return;
    }
    if x <= 12.0 {
        let x2 = x * x;
        for k in 0..=n {
            let mut bk = 1.0;
            let mut r = 1.0;
            for i in 1..=50 {
                // Unlike in xsf, `i * (i + k)` is computed in `double` instead of `int`, which
                // overflows for k >~ 4e7; this is exact, since `i <= 50`
                r = -0.25 * r * x2 / (f64::from(i) * f64::from(i + k));
                bk += r;

                if r.abs() < bk.abs() * 1.0e-15 {
                    break;
                }
            }
            bl[k as usize] = bk;
            if k >= 1 {
                dl[(k - 1) as usize] = -0.5 * x / f64::from(k) * bk;
            }
        }
        let mut uk = 1.0;
        let mut r = 1.0;
        for i in 1..=50 {
            r = -0.25 * r * x2 / (f64::from(i) * (f64::from(i + n) + 1.0));
            uk += r;

            if r.abs() < uk.abs() * 1.0e-15 {
                break;
            }
        }
        dl[n as usize] = -0.5 * x / (f64::from(n) + 1.0) * uk;
        return;
    }
    if n == 0 {
        *nm = 1;
    }
    let mut m = msta1(x, 200);
    if m < *nm {
        *nm = m;
    } else {
        m = msta2(x, *nm, 15);
    }
    let mut bs = 0.0;
    let mut f = 0.0;
    let mut f0 = 0.0;
    let mut f1 = 1e-100;
    for k in (0..=m).rev() {
        f = 2.0 * (f64::from(k) + 1.0) * f1 / x - f0;
        if k <= *nm {
            bl[k as usize] = f;
        }
        if k % 2 == 0 {
            bs += 2.0 * f;
        }
        f0 = f1;
        f1 = f;
    }
    let bg = bs - f;
    for k in 0..=*nm {
        bl[k as usize] /= bg;
    }
    let mut r0 = 1.0;
    for k in 1..=*nm {
        r0 = 2.0 * r0 * f64::from(k) / x;
        bl[k as usize] *= r0;
    }
    dl[0] = -0.5 * x * bl[1];
    for k in 1..=*nm {
        dl[k as usize] = 2.0 * f64::from(k) / x * (bl[(k - 1) as usize] - bl[k as usize]);
    }
}

#[allow(
    clippy::approx_constant,
    clippy::assign_op_pattern,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::manual_midpoint,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::unreadable_literal
)]
pub(crate) fn lamv(v: f64, x: f64, vm: &mut f64, vl: &mut [f64], dl: &mut [f64]) {
    // =========================================================
    // Purpose: Compute lambda function with arbitrary order v,
    //          and their derivative
    // Input :  x --- Argument of lambda function
    //          v --- Order of lambda function
    // Output:  VL(n) --- Lambda function of order n+v0
    //          DL(n) --- Derivative of lambda function
    //          VM --- Highest order computed
    // Routines called:
    //      (1) MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    //      (2) GAM0 for computing gamma function (|x| ≤ 1)
    // =========================================================

    // NOTE: vl[0..=n] and dl[0..=n] are written with n = (int)v, and for x > 12 also vl[1] / dl[1]
    // (out of bounds for v < 1 with buffers of length n + 1). The caller
    // (jahnke_emden_lambda) passes 1 <= v < 2^31 with buffers of length n + 1, so that's in bounds,
    // and `(int)v` is not UB.
    let rp2 = 0.63661977236758;

    let x = x.abs();
    let x2 = x * x;
    let mut n = v as i32;
    let v0 = v - f64::from(n);
    *vm = v;

    if x <= 12.0 {
        for k in 0..=n {
            let vk = v0 + f64::from(k);
            let mut bk = 1.0;
            let mut r = 1.0;

            for i in 1..=50 {
                r = -0.25 * r * x2 / (f64::from(i) * (f64::from(i) + vk));
                bk = bk + r;

                if r.abs() < bk.abs() * 1.0e-15 {
                    break;
                }
            }
            vl[k as usize] = bk;

            let mut uk = 1.0;
            r = 1.0;

            for i in 1..=50 {
                r = -0.25 * r * x2 / (f64::from(i) * (f64::from(i) + vk + 1.0));
                uk = uk + r;

                if r.abs() < uk.abs() * 1.0e-15 {
                    break;
                }
            }
            dl[k as usize] = -0.5 * x / (vk + 1.0) * uk;
        }
        return;
    }

    let k0 = if x >= 50.0 {
        8
    } else if x >= 35.0 {
        10
    } else {
        11
    };
    let mut bjv0 = 0.0;
    let mut bjv1 = 0.0;

    for j in 0..=1 {
        let vv = 4.0 * (f64::from(j) + v0) * (f64::from(j) + v0);
        let mut px = 1.0;
        let mut rp = 1.0;
        for k in 1..=k0 {
            // NOTE: `pow(y, 2.0)` is written as `y * y`, which is what GCC compiles it to.
            rp = -0.78125e-2
                * rp
                * (vv - (4.0 * f64::from(k) - 3.0) * (4.0 * f64::from(k) - 3.0))
                * (vv - (4.0 * f64::from(k) - 1.0) * (4.0 * f64::from(k) - 1.0))
                / (f64::from(k) * (2.0 * f64::from(k) - 1.0) * x2);
            px += rp;
        }

        let mut qx = 1.0;
        rp = 1.0;
        for k in 1..=k0 {
            rp = -0.78125e-2
                * rp
                * (vv - (4.0 * f64::from(k) - 1.0) * (4.0 * f64::from(k) - 1.0))
                * (vv - (4.0 * f64::from(k) + 1.0) * (4.0 * f64::from(k) + 1.0))
                / (f64::from(k) * (2.0 * f64::from(k) + 1.0) * x2);
            qx += rp;
        }

        qx = 0.125 * (vv - 1.0) * qx / x;
        let xk = x - (0.5 * (f64::from(j) + v0) + 0.25) * PI;
        let a0 = (rp2 / x).sqrt();
        let ck = xk.cos();
        let sk = xk.sin();

        if j == 0 {
            bjv0 = a0 * (px * ck - qx * sk);
        }
        if j == 1 {
            bjv1 = a0 * (px * ck - qx * sk);
        }
    }

    let mut ga;
    if v0 == 0.0 {
        ga = 1.0;
    } else {
        ga = gam0(v0);
        ga *= v0;
    }

    let fac = (2.0 / x).powf(v0) * ga;
    vl[0] = bjv0;
    dl[0] = -bjv1 + v0 / x * bjv0;
    vl[1] = bjv1;
    dl[1] = bjv0 - (1.0 + v0) / x * bjv1;
    let mut r0 = 2.0 * (1.0 + v0) / x;

    if n <= 1 {
        vl[0] *= fac;
        dl[0] = fac * dl[0] - v0 / x * vl[0];
        vl[1] *= fac * r0;
        dl[1] = fac * r0 * dl[1] - (1.0 + v0) / x * vl[1];
        return;
    }

    // NOTE: `(int)(0.9 * x)` is UB in C++ for x >= 2^31 / 0.9 (Rust saturates); the caller only
    // passes |x| <= 1e9.
    if n >= 2 && n <= (0.9 * x) as i32 {
        let mut f0 = bjv0;
        let mut f1 = bjv1;

        for k in 2..=n {
            let f = 2.0 * (f64::from(k) + v0 - 1.0) / x * f1 - f0;
            f0 = f1;
            f1 = f;
            vl[k as usize] = f;
        }
    } else if n >= 2 {
        let mut m = msta1(x, 200);
        if m < n {
            n = m;
        } else {
            m = msta2(x, n, 15);
        }

        let mut f = 0.0;
        let mut f2 = 0.0;
        let mut f1 = 1.0e-100;

        for k in (0..=m).rev() {
            f = 2.0 * (v0 + f64::from(k) + 1.0) / x * f1 - f2;
            if k <= n {
                vl[k as usize] = f;
            }
            f2 = f1;
            f1 = f;
        }

        let cs = if bjv0.abs() > bjv1.abs() {
            bjv0 / f
        } else {
            bjv1 / f2
        };

        for k in 0..=n {
            vl[k as usize] *= cs;
        }
    }

    vl[0] *= fac;
    for j in 1..=n {
        vl[j as usize] *= fac * r0;
        dl[(j - 1) as usize] = -0.5 * x / (f64::from(j) + v0) * vl[j as usize];
        r0 = 2.0 * (f64::from(j) + v0 + 1.0) / x * r0;
    }

    dl[n as usize] = 2.0 * (v0 + f64::from(n)) * (vl[(n - 1) as usize] - vl[n as usize]) / x;
    *vm = f64::from(n) + v0;
}

#[allow(clippy::approx_constant, clippy::cast_possible_truncation)]
pub(crate) fn msta1(x: f64, mp: i32) -> i32 {
    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that the magnitude of
    //          Jn(x) at that point is about 10^(-MP)
    // Input :  x     --- Argument of Jn(x)
    //          MP    --- Value of magnitude
    // Output:  MSTA1 --- Starting point
    // ===================================================

    // NOTE: the `(int)` casts below are UB in C++ for |x| >= 2^31 / 1.1 (~1.95e9) or NaN, and if
    // the secant step is NaN, +-inf (e.g. for x == 0) or huge. Rust's `as i32` saturates
    // (NaN -> 0), whereas g++ on x86-64 yields INT_MIN. The `+ 1` / `+ 5` / `nn - n1` below can
    // then overflow i32 (panics in debug builds). jynbh only calls this for x >= 1e-100 with
    // n > (int)(0.9 * x) or x <= 300.
    let a0 = x.abs();
    let mut n0 = (1.1 * a0) as i32 + 1;
    let mut f0 = 0.5 * (6.28 * f64::from(n0)).log10()
        - f64::from(n0) * (1.36 * a0 / f64::from(n0)).log10()
        - f64::from(mp);
    let mut n1 = n0 + 5;
    let mut f1 = 0.5 * (6.28 * f64::from(n1)).log10()
        - f64::from(n1) * (1.36 * a0 / f64::from(n1)).log10()
        - f64::from(mp);
    let mut nn = 0;
    for _it in 1..=20 {
        nn = (f64::from(n1) - f64::from(n1 - n0) / (1.0 - f0 / f1)) as i32;
        let f = 0.5 * (6.28 * f64::from(nn)).log10()
            - f64::from(nn) * (1.36 * a0 / f64::from(nn)).log10()
            - f64::from(mp);
        if (nn - n1).abs() < 1 {
            break;
        }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    nn
}

#[allow(clippy::approx_constant, clippy::cast_possible_truncation)]
pub(crate) fn msta2(x: f64, n: i32, mp: i32) -> i32 {
    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that all Jn(x) has MP
    //          significant digits
    // Input :  x  --- Argument of Jn(x)
    //          n  --- Order of Jn(x)
    //          MP --- Significant digit
    // Output:  MSTA2 --- Starting point
    // ===================================================

    // NOTE: the `(int)` casts below have the same UB as in msta1 (see there), and `nn + 10` can
    // overflow as well.
    let mut n0: i32;
    let obj: f64;

    let a0 = x.abs();
    let hmp = 0.5 * f64::from(mp);
    let ejn =
        0.5 * (6.28 * f64::from(n)).log10() - f64::from(n) * (1.36 * a0 / f64::from(n)).log10();
    if ejn <= hmp {
        obj = f64::from(mp);
        n0 = (1.1 * a0) as i32 + 1;
    } else {
        obj = hmp + ejn;
        n0 = n;
    }
    let mut f0 = 0.5 * (6.28 * f64::from(n0)).log10()
        - f64::from(n0) * (1.36 * a0 / f64::from(n0)).log10()
        - obj;
    let mut n1 = n0 + 5;
    let mut f1 = 0.5 * (6.28 * f64::from(n1)).log10()
        - f64::from(n1) * (1.36 * a0 / f64::from(n1)).log10()
        - obj;
    let mut nn = 0;
    for _it in 1..=20 {
        nn = (f64::from(n1) - f64::from(n1 - n0) / (1.0 - f0 / f1)) as i32;
        let f = 0.5 * (6.28 * f64::from(nn)).log10()
            - f64::from(nn) * (1.36 * a0 / f64::from(nn)).log10()
            - obj;
        if (nn - n1).abs() < 1 {
            break;
        }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    nn + 10
}

#[allow(clippy::cast_sign_loss)]
pub(crate) fn sphj(x: f64, n: i32, nm: &mut i32, sj: &mut [f64], dj: &mut [f64]) {
    //  MODIFIED to ALLOW N=0 CASE (ALSO IN SPHY)
    //
    //  =======================================================
    //  Purpose: Compute spherical Bessel functions jn(x) and
    //           their derivatives
    //  Input :  x --- Argument of jn(x)
    //           n --- Order of jn(x)  ( n = 0,1,… )
    //  Output:  SJ(n) --- jn(x)
    //           DJ(n) --- jn'(x)
    //           NM --- Highest order computed
    //  Routines called:
    //           MSTA1 and MSTA2 for computing the starting
    //           point for backward recurrence
    //  =======================================================

    // NOTE: `sj[0..=*nm]` / `dj[0..=*nm]` are written, with *nm = min(n, msta1(x, 200)) for
    // |x| >= 1e-100. The only caller, rmn1, uses 252-element buffers, so this overflows if both n
    // and msta1(x, 200) exceed 251 (c * x >~ 30), e.g. prolate_radial1_nocv(2, 200, 1, 300)
    // (confirmed with ASan).
    // NOTE: for x NaN, msta1/msta2 hit their `(int)NaN` UB (see there), so `*nm` (and therefore
    // which elements of `sj` / `dj` are written) may differ from C++; the values are all NaN.
    *nm = n;
    if x.abs() < 1e-100 {
        for k in 0..=n {
            sj[k as usize] = 0.0;
            dj[k as usize] = 0.0;
        }
        sj[0] = 1.0;
        if n > 0 {
            dj[1] = 1.0 / 3.0;
        }
        return;
    }
    sj[0] = x.sin() / x;
    dj[0] = (x.cos() - x.sin() / x) / x;
    if n < 1 {
        return;
    }
    sj[1] = (sj[0] - x.cos()) / x;
    if n >= 2 {
        let sa = sj[0];
        let sb = sj[1];
        let mut m = msta1(x, 200);
        if m < n {
            *nm = m;
        } else {
            m = msta2(x, n, 15);
        }
        let mut f = 0.0;
        let mut f0 = 0.0;
        let mut f1 = 1e-100;
        for k in (0..=m).rev() {
            f = (2.0 * f64::from(k) + 3.0) * f1 / x - f0;
            if k <= *nm {
                sj[k as usize] = f;
            }
            f0 = f1;
            f1 = f;
        }
        let cs = if sa.abs() > sb.abs() { sa / f } else { sb / f0 };
        for k in 0..=*nm {
            sj[k as usize] *= cs;
        }
    }
    for k in 1..=*nm {
        dj[k as usize] = sj[(k - 1) as usize] - (f64::from(k) + 1.0) * sj[k as usize] / x;
    }
}

#[allow(clippy::cast_sign_loss)]
pub(crate) fn sphy(x: f64, n: i32, nm: &mut i32, sy: &mut [f64], dy: &mut [f64]) {
    // ======================================================
    // Purpose: Compute spherical Bessel functions yn(x) and
    //          their derivatives
    // Input :  x --- Argument of yn(x) ( x ≥ 0 )
    //          n --- Order of yn(x) ( n = 0,1,… )
    // Output:  SY(n) --- yn(x)
    //          DY(n) --- yn'(x)
    //          NM --- Highest order computed
    // ======================================================

    // NOTE: `sy[0..=n]` / `dy[0..=n]` are written, but the loop stops after the first
    // |sy[k]| >= 1e300. The only caller, rmn2l, uses 252-element buffers, so this overflows if
    // n > 251 and sy[251] is finite (c * x >~ 11.75), e.g. prolate_radial2_nocv(2, 200, 1, 300)
    // (confirmed with ASan).
    //
    // NOTE: unlike the Fortran original, the C++ computes no dy[k], k >= 1, if |f| >= 1e300, and
    // otherwise sets *nm = n - 1, so dy[n] is never computed. Reproduced as is.
    if x < 1.0e-60 {
        for k in 0..=n {
            sy[k as usize] = -1.0e300;
            dy[k as usize] = 1.0e300;
        }
        *nm = n;
        return;
    }
    sy[0] = -x.cos() / x;
    let mut f0 = sy[0];
    dy[0] = (x.sin() + x.cos() / x) / x;

    if n < 1 {
        *nm = n;
        return;
    }

    sy[1] = (sy[0] - x.sin()) / x;
    let mut f1 = sy[1];

    for k in 2..=n {
        let f = ((2.0 * f64::from(k) - 1.0) * f1 / x) - f0;
        sy[k as usize] = f;
        if f.abs() >= 1.0e300 {
            *nm = k - 1;
            return;
        }
        f0 = f1;
        f1 = f;
    }
    *nm = n - 1;
    for k in 1..=*nm {
        dy[k as usize] = sy[(k - 1) as usize] - (f64::from(k) + 1.0) * sy[k as usize] / x;
    }
}
