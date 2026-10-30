//! Bessel functions of integer order, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.

use core::f64::consts::PI;

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

    // Unlike in xsf, the Hankel expansion is used for x >= 2^31 / 0.9 (~2.39e9), because `as i32`
    // saturates. In C++ `(int)(0.9 * x)` is UB there, which on x86-64 yields INT_MIN, so that the
    // backward recurrence is taken, resulting in 0 or NaN.
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
                // NOTE: `k * k` is `int` arithmetic, as in C++; it overflows i32 for k > 46340 (C++ UB, Rust
                // panics in debug), i.e. only if msta1/msta2 return such an m (x > 300 with n > 0.9 * x >~ 4e4).
                sv +=
                    (-1.0_f64).powf(f64::from(k / 2)) * f64::from(k) / (f64::from(k * k) - 1.0) * f;
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

    // NOTE: the `(int)` casts below are UB in C++ when the double is out of range or NaN: for
    // |x| >= 2^31 / 1.1 (~1.95e9) or NaN in `(int)(1.1 * a0)`, and in the `nn = ...` assignment
    // whenever the secant step is NaN/+-inf (e.g. x == 0, where f0 = f1 = +inf) or huge. Rust's
    // `as i32` saturates (NaN -> 0), whereas g++ on x86-64 yields INT_MIN. Moreover, with
    // saturation, the `+ 1` / `+ 5` and `nn - n1` below can overflow i32 (panics in debug builds)
    // for such inputs. jynbh only calls this for x >= 1e-100 with n > (int)(0.9 * x) or x <= 300.
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

    // NOTE: as in msta1, the `(int)` casts below are UB in C++ when the double is out of range or
    // NaN (|x| >= 2^31 / 1.1, x NaN, or a NaN/+-inf/huge secant step); Rust's `as i32` saturates
    // (NaN -> 0), whereas g++ on x86-64 yields INT_MIN, and the subsequent `+ 1` / `+ 5` /
    // `nn - n1` / `nn + 10` can then overflow i32 (panics in debug builds).
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

    // NOTE: `sj[0..=n]` / `dj[0..=n]` are written (and `sj[1]` for n >= 1, `dj[1]` for n > 0).
    // The only caller, rmn1, passes n = 2 * (25 + (n - m) / 2 + (int)c) + m into 252-element
    // buffers. Since only sj[0..=*nm] / dj[0..=*nm] are written for |x| >= 1e-100 (with
    // *nm = min(n, msta1(x, 200))), this overflows when that n exceeds 251, i.e.
    // n_caller - (n_caller - m) % 2 + 2 * (int)c >= 202, and also msta1(c * x, 200) >= 252, i.e.
    // c * x >~ 30 (e.g. prolate_radial1_nocv(2, 200, 1, 300), confirmed with ASan), which the
    // public API allows.
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

    // NOTE: `sy[0..=n]` / `dy[0..=n]` are written. The only caller, rmn2l, passes
    // n = 2 * (25 + (n - m) / 2 + (int)c) + m into 252-element buffers. Since the loop stops after
    // writing the first |sy[k]| >= 1e300, this overflows when that n exceeds 251 (same condition as
    // for sphj) and sy[251] is still finite, i.e. c * x >~ 11.75 (e.g.
    // prolate_radial2_nocv(2, 200, 1, 300), confirmed with ASan).
    // NOTE: unlike the Fortran original (where NM = K - 1 is N after a completed loop, and DY is
    // computed up to NM in both cases), the C++ returns early without computing any dy[k], k >= 1,
    // when |f| >= 1e300, and otherwise sets *nm = n - 1, so dy[n] is never computed (for n == 1,
    // dy[1] isn't either). Reproduced as is.
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
