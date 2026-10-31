//! Spheroidal wave functions, translated into pure Rust from `xsf/specfun/specfun.h` (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.

use core::f64::consts::PI;

use super::{Status, lpmns, lqmns, sphj, sphy, try_zeroed};
use crate::xsf::tools::int_x86;

#[allow(
    clippy::too_many_arguments,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::bool_to_int_with_if
)]
pub(crate) fn aswfa(
    x: f64,
    m: i32,
    n: i32,
    c: f64,
    kd: i32,
    cv: f64,
    s1f: &mut f64,
    s1d: &mut f64,
) -> Status {
    // ===========================================================
    // Purpose: Compute the prolate and oblate spheroidal angular
    //          functions of the first kind and their derivatives
    // Input :  m  --- Mode parameter,  m = 0,1,2,...
    //          n  --- Mode parameter,  n = m,m+1,...
    //          c  --- Spheroidal parameter
    //          x  --- Argument of angular function, |x| < 1.0
    //          KD --- Function code
    //                 KD=1 for prolate;  KD=-1 for oblate
    //          cv --- Characteristic value
    // Output:  S1F --- Angular function of the first kind
    //          S1D --- Derivative of the angular function of
    //                  the first kind
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //
    // Routine called:
    //          SDMN for computing expansion coefficients df
    //          SCKB for computing expansion coefficients ck
    // ===========================================================

    let ck = try_zeroed(200);
    let df = try_zeroed(200);
    let (Some(mut ck), Some(mut df)) = (ck, df) else {
        return Status::NoMemory;
    };
    let eps = 1e-14;
    let x0 = x;
    let x = x.abs();
    let ip = if (n - m) % 2 == 0 { 0 } else { 1 };
    // NOTE: `(int)` is UB in C++ if `(n - m) / 2 + c` is NaN or out of the int range; `as i32`
    // saturates (NaN -> 0).
    let nm = 40 + (f64::from((n - m) / 2) + c) as i32;
    let nm2 = nm / 2 - 2;
    if sdmn(m, n, c, cv, kd, &mut df) == Status::NoMemory {
        return Status::NoMemory;
    }
    sckb(m, n, c, &df, &mut ck);
    let x1 = 1.0 - x * x;
    let a0 = if (m == 0) && (x1 == 0.0) {
        1.0
    } else {
        x1.powf(0.5 * f64::from(m))
    };
    let mut su1 = ck[0];
    // NOTE: ck has 200 elements, so ck[k] (here and in the su2 loop below) is out of bounds (UB in
    // C++) if nm2 >= 200, i.e. if (n - m) / 2 + c >= 364 (only if the loop doesn't break before
    // that); but sdmn and sckb already go out of bounds for much smaller c (see there).
    for k in 1..=nm2 {
        let r = ck[k as usize] * x1.powf(f64::from(k));
        su1 += r;
        if (k >= 10) && ((r / su1).abs() < eps) {
            break;
        }
    }
    *s1f = a0 * x.powf(f64::from(ip)) * su1;
    if x == 1.0 {
        if m == 0 {
            *s1d = f64::from(ip) * ck[0] - 2.0 * ck[1];
        } else if m == 1 {
            *s1d = -1e100;
        } else if m == 2 {
            *s1d = -2.0 * ck[0];
        } else if m >= 3 {
            *s1d = 0.0;
        }
    } else {
        let d0 = f64::from(ip) - f64::from(m) / x1 * x.powf(f64::from(ip) + 1.0);
        let d1 = -2.0 * a0 * x.powf(f64::from(ip) + 1.0);
        let mut su2 = ck[1];
        for k in 2..=nm2 {
            let r = f64::from(k) * ck[k as usize] * x1.powf(f64::from(k) - 1.0);
            su2 += r;
            if (k >= 10) && ((r / su2).abs() < eps) {
                break;
            }
        }
        *s1d = d0 * a0 * su1 + d1 * su2;
    }
    if (x0 < 0.0) && (ip == 0) {
        *s1d = -*s1d;
    }
    if (x0 < 0.0) && (ip == 1) {
        *s1f = -*s1f;
    }
    // x = x0;  (dead store in C++)
    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::bool_to_int_with_if
)]
pub(crate) fn cbk(m: i32, n: i32, c: f64, cv: f64, qt: f64, ck: &[f64], bk: &mut [f64]) -> Status {
    // ==========================================================
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    // ==========================================================

    let eps = 1.0e-14;

    let ip = if (n - m) % 2 == 0 { 0 } else { 1 };
    // NOTE: `(int)` is UB in C++ if `0.5 * (n - m) + c` is NaN or out of the int range; `as i32`
    // saturates (NaN -> 0).
    let nm = 25 + (0.5 * f64::from(n - m) + c) as i32;

    let u = try_zeroed(200);
    let v = try_zeroed(200);
    let w = try_zeroed(200);
    let (Some(mut u), Some(mut v), Some(mut w)) = (u, v, w) else {
        return Status::NoMemory;
    };

    // NOTE: u, v, w (and, in the caller rmn2so, ck and bk) have 200 elements. With
    // nm = 25 + (int)(0.5 * (n - m) + c), the accesses below are out of bounds (UB in C++) if
    // nm >= 200 (ck[i] with i <= nm; only if the inner loop doesn't break before), nm >= 202
    // (u[j] with j <= nm - 2, w[j - 1] with j <= nm - 1), or nm >= 203 (v[j - 1], bk[k] with
    // k <= nm - 3), i.e. for 0.5 * (n - m) + c >= 175.
    u[0] = 0.0;
    let n2 = nm - 2;

    for j in 1..=n2 {
        u[j as usize] = c * c;
    }
    for j in 1..=n2 {
        v[(j - 1) as usize] = (2.0 * f64::from(j) - 1.0 - f64::from(ip))
            * (2.0 * f64::from(j - m) - f64::from(ip))
            + f64::from(m) * (f64::from(m) - 1.0)
            - cv;
    }
    for j in 1..nm {
        w[(j - 1) as usize] =
            (2.0 * f64::from(j) - f64::from(ip)) * (2.0 * f64::from(j) + 1.0 - f64::from(ip));
    }

    if ip == 0 {
        let mut sw = 0.0;
        for k in 0..n2 {
            let mut s1 = 0.0;
            let i1 = k - m + 1;

            for i in i1..=nm {
                if i < 0 {
                    continue;
                }
                let mut r1 = 1.0;
                for j in 1..=k {
                    r1 = r1 * f64::from(i + m - j) / (1.0 * f64::from(j));
                }
                s1 += ck[i as usize] * (2.0 * f64::from(i) + f64::from(m)) * r1;
                if (s1 - sw).abs() < s1.abs() * eps {
                    break;
                }
                sw = s1;
            }
            bk[k as usize] = qt * s1;
        }
    } else if ip == 1 {
        let mut sw = 0.0;
        for k in 0..n2 {
            let mut s1 = 0.0;
            let i1 = k - m + 1;

            for i in i1..=nm {
                if i < 0 {
                    continue;
                }
                let mut r1 = 1.0;
                for j in 1..=k {
                    r1 = r1 * f64::from(i + m - j) / (1.0 * f64::from(j));
                }
                if i > 0 {
                    s1 += ck[(i - 1) as usize] * (2.0 * f64::from(i) + f64::from(m) - 1.0) * r1;
                }
                s1 -= ck[i as usize] * (2.0 * f64::from(i) + f64::from(m)) * r1;

                if (s1 - sw).abs() < s1.abs() * eps {
                    break;
                }
                sw = s1;
            }
            bk[k as usize] = qt * s1;
        }
    }

    w[0] /= v[0];
    bk[0] /= v[0];

    for k in 2..=n2 {
        let t = v[(k - 1) as usize] - w[(k - 2) as usize] * u[(k - 1) as usize];
        w[(k - 1) as usize] /= t;
        bk[(k - 1) as usize] =
            (bk[(k - 1) as usize] - bk[(k - 2) as usize] * u[(k - 1) as usize]) / t;
    }

    for k in (1..n2).rev() {
        bk[(k - 1) as usize] -= w[(k - 1) as usize] * bk[k as usize];
    }
    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::bool_to_int_with_if
)]
pub(crate) fn gmn(m: i32, n: i32, c: f64, x: f64, bk: &[f64], gf: &mut f64, gd: &mut f64) {
    // ===========================================================
    // Purpose: Compute gmn(-ic,ix) and its derivative for oblate
    //          radial functions with a small argument
    // ===========================================================

    let eps = 1.0e-14;
    let ip = if (n - m) % 2 == 0 { 0 } else { 1 };
    // NOTE: `(int)` is UB in C++ if `0.5 * (n - m) + c` is NaN or out of the int range; `as i32`
    // saturates (NaN -> 0).
    let nm = 25 + (0.5 * f64::from(n - m) + c) as i32;
    let xm = (1.0 + x * x).powf(-0.5 * f64::from(m));
    let mut gf0 = 0.0;
    let mut gw = 0.0;

    // NOTE: bk has 200 elements in the caller (rmn2so), so bk[k - 1] is out of bounds (UB in C++)
    // if nm > 200, i.e. 0.5 * (n - m) + c >= 176 (only if the loops don't break before that).
    for k in 1..=nm {
        gf0 += bk[(k - 1) as usize] * x.powf(2.0 * f64::from(k) - 2.0);
        if (((gf0 - gw) / gf0).abs() < eps) && (k >= 10) {
            break;
        }
        gw = gf0;
    }

    *gf = xm * gf0 * x.powf(f64::from(1 - ip));
    let gd1 = f64::from(-m) * x / (1.0 + x * x) * (*gf);
    let mut gd0 = 0.0;

    // NOTE: gw is not reset here, so the first convergence check compares against gf0
    for k in 1..nm {
        if ip == 0 {
            gd0 += (2.0 * f64::from(k) - 1.0)
                * bk[(k - 1) as usize]
                * x.powf(2.0 * f64::from(k) - 2.0);
        } else {
            // Unlike in xsf, `bk[k]` (as `BK(K+1)` in the Fortran original), instead of `bk[k - 1]`
            gd0 += 2.0 * f64::from(k) * bk[k as usize] * x.powf(2.0 * f64::from(k) - 1.0);
        }
        if (((gd0 - gw) / gd0).abs() < eps) && (k >= 10) {
            break;
        }
        gw = gd0;
    }
    *gd = gd1 + xm * gd0;
}

#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::assign_op_pattern,
    clippy::bool_to_int_with_if,
    clippy::range_minus_one
)]
pub(crate) fn kmn(
    m: i32,
    n: i32,
    c: f64,
    cv: f64,
    kd: i32,
    df: &[f64],
    dn: &mut [f64],
    ck1: &mut f64,
    ck2: &mut f64,
) -> Status {
    // ===================================================
    // Purpose: Compute the expansion coefficients of the
    //          prolate and oblate spheroidal functions
    //          and joining factors
    //
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    // ===================================================

    // NOTE: out-of-bounds indexing (C++ UB; panics in Rust), where
    // nm = 25 + (int)(0.5 * (n - m) + c) and nn = nm + m:
    // - `dn[0..nn]` is written, but the callers allocate `dn` with 200 (rmn2so) or 201 (rmn2sp)
    //   elements, so this overflows when nn > 200 (rmn2so) / nn > 201 (rmn2sp), e.g. m = n = 180,
    //   c = 0.5 (nn = 205); this is reachable through the public API since only n - m <= 198 is
    //   checked (e.g. prolate_radial2_nocv(180, 180, 0.5, x = 1.0001) and
    //   oblate_radial2_nocv(180, 180, 0.5, x = 0.1), confirmed with ASan).
    // - `df[k - 1]` is read for k up to nm (unless the sum converges earlier); the callers'
    //   `df` has 200 elements, so this overflows when nm > 200, i.e. 0.5 * (n - m) + c >= 176
    //   (e.g. n - m = 198 and c >= 77); sdmn has already overflowed `df` in that case.
    // NOTE: `rk[nn - 1]` is never written (the second loop only writes rk[m + 1..=nn - 2]) but is
    // read in the `dn` loop, so `dn[nn - 1]` is computed from uninitialized memory in C++ (here it
    // is 0 * dn[nn - 2], i.e. +-0, or NaN if dn[nn - 2] is +-inf/NaN). kmn itself never reads
    // `dn[nn - 1]` (only dn[m - 1]), but rmn2sp may read it in its `su2` / `sd2` loops. The
    // Fortran original has the same bug (RK(NN) is read but never set).
    let nm: i32 = 25 + (0.5 * f64::from(n - m) + c) as i32;
    let nn: i32 = nm + m;
    let u = try_zeroed((nn + 4) as usize);
    let v = try_zeroed((nn + 4) as usize);
    let w = try_zeroed((nn + 4) as usize);
    let tp = try_zeroed((nn + 4) as usize);
    let rk = try_zeroed((nn + 4) as usize);
    let (Some(mut u), Some(mut v), Some(mut w), Some(mut tp), Some(mut rk)) = (u, v, w, tp, rk)
    else {
        return Status::NoMemory;
    };

    let eps = 1.0e-14;

    let cs = c * c * f64::from(kd);
    *ck1 = 0.0;
    *ck2 = 0.0;

    let ip: i32 = if (n - m) % 2 == 0 { 0 } else { 1 };

    for i in 1..=nn + 3 {
        let k: i32 = if ip == 0 { -2 * (i - 1) } else { -(2 * i - 3) };
        let gk0 = 2.0 * f64::from(m) + f64::from(k);
        let gk1 = f64::from(m + k) * (f64::from(m + k) + 1.0);
        let gk2 = 2.0 * f64::from(m + k) - 1.0;
        let gk3 = 2.0 * f64::from(m + k) + 3.0;

        u[(i - 1) as usize] = gk0 * (gk0 - 1.0) * cs / (gk2 * (gk2 + 2.0));
        v[(i - 1) as usize] = gk1 - cv + (2.0 * (gk1 - f64::from(m * m)) - 1.0) * cs / (gk2 * gk3);
        w[(i - 1) as usize] =
            (f64::from(k) + 1.0) * (f64::from(k) + 2.0) * cs / ((gk2 + 2.0) * gk3);
    }

    for k in 1..=m {
        let mut t = v[m as usize];
        for l in 0..=m - k - 1 {
            t = v[(m - l - 1) as usize] - w[(m - l) as usize] * u[(m - l - 1) as usize] / t;
        }

        rk[(k - 1) as usize] = -u[(k - 1) as usize] / t;
    }

    let mut r = 1.0;
    for k in 1..=m {
        r = r * rk[(k - 1) as usize];
        dn[(k - 1) as usize] = df[0] * r;
    }

    tp[(nn - 1) as usize] = v[nn as usize];

    for k in (m + 1..=nn - 1).rev() {
        tp[(k - 1) as usize] = v[k as usize] - w[(k + 1) as usize] * u[k as usize] / tp[k as usize];

        if k > m + 1 {
            rk[(k - 1) as usize] = -u[(k - 1) as usize] / tp[(k - 1) as usize];
        }
    }

    let dnp = if m == 0 { df[0] } else { dn[(m - 1) as usize] };
    dn[m as usize] = (-1.0_f64).powf(f64::from(ip)) * dnp * cs
        / ((2.0 * f64::from(m) - 1.0)
            * (2.0 * f64::from(m) + 1.0 - 4.0 * f64::from(ip))
            * tp[m as usize]);

    for k in m + 2..=nn {
        dn[(k - 1) as usize] = rk[(k - 1) as usize] * dn[(k - 2) as usize];
    }

    let mut r1 = 1.0;
    for j in 1..=(n + m + ip) / 2 {
        r1 = r1 * (f64::from(j) + 0.5 * f64::from(n + m + ip));
    }
    r = 1.0;
    for j in 1..=2 * m + ip {
        r *= f64::from(j);
    }
    let mut su0 = r * df[0];
    let mut sw = 0.0;

    for k in 2..=nm {
        r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
            / (f64::from(k) - 1.0)
            / (f64::from(k + ip) - 1.5);
        su0 = su0 + r * df[(k - 1) as usize];
        if k > (n - m) / 2 && ((su0 - sw) / su0).abs() < eps {
            break;
        }
        sw = su0;
    }

    if kd != 1 {
        let mut r2 = 1.0;

        for j in 1..=m {
            r2 = 2.0 * c * r2 * f64::from(j);
        }

        let mut r3 = 1.0;

        for j in 1..=(n - m - ip) / 2 {
            r3 = r3 * f64::from(j);
        }

        let sa0 = (2.0 * f64::from(m + ip) + 1.0) * r1
            / (f64::from(n).exp2() * c.powf(f64::from(ip)) * r2 * r3 * df[0]);
        *ck1 = sa0 * su0;

        if kd == -1 {
            return Status::Ok;
        }
    }

    let mut r4 = 1.0;
    for j in 1..=(n - m - ip) / 2 {
        r4 *= 4.0 * f64::from(j);
    }
    let mut r5 = 1.0;
    for j in 1..=m {
        r5 = r5 * f64::from(j + m) / c;
    }

    let g0 = if m == 0 { df[0] } else { dn[(m - 1) as usize] };

    let sb0 = (f64::from(ip) + 1.0) * c.powf(f64::from(ip + 1))
        / (2.0 * f64::from(ip) * (f64::from(m) - 2.0) + 1.0)
        / (2.0 * f64::from(m) - 1.0);
    *ck2 = (-1.0_f64).powf(f64::from(ip)) * sb0 * r4 * r5 * g0 / r1 * su0;

    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::bool_to_int_with_if,
    clippy::assign_op_pattern
)]
pub(crate) fn qstar(
    m: i32,
    n: i32,
    c: f64,
    ck1: f64,
    ck: &[f64],
    qs: &mut f64,
    qt: &mut f64,
) -> Status {
    // ==========================================================
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    // ==========================================================

    // NOTE: uninitialized (`new T[200]`) in C++; the only element that can be read before being
    // written is ap[-1] (for m == 0, see below).
    let Some(mut ap) = try_zeroed(200) else {
        return Status::NoMemory;
    };
    let ip = if (n - m) == 2 * ((n - m) / 2) { 0 } else { 1 };
    let mut r = 1.0 / (ck[0] * ck[0]);
    ap[0] = r;

    // NOTE: ap and ck (in the caller rmn2so) have 200 elements, so ap[i] and ck[k] (k <= m) are
    // out of bounds (UB in C++) if m >= 200.
    for i in 1..=m {
        let mut s = 0.0;
        for l in 1..=i {
            let mut sk = 0.0;
            for k in 0..=l {
                sk += ck[k as usize] * ck[(l - k) as usize];
            }
            s += sk * ap[(i - l) as usize];
        }
        ap[i as usize] = -r * s;
    }
    // Unlike in xsf, `ap[m]` and `l <= m` (as `QS0=AP(M+1)` and `DO 30 L=1,M` in the Fortran
    // original), instead of `ap[m - 1]` and `l < m`, which reads out of bounds for `m == 0`
    let mut qs0 = ap[m as usize];

    for l in 1..=m {
        r = 1.0;
        for k in 1..=l {
            r = r
                * (2.0 * f64::from(k) + f64::from(ip))
                * (2.0 * f64::from(k) - 1.0 + f64::from(ip))
                / ((2.0 * f64::from(k)) * (2.0 * f64::from(k)));
        }
        qs0 += ap[(m - l) as usize] * r;
    }
    *qs = (-1.0_f64).powf(f64::from(ip)) * ck1 * (ck1 * qs0) / c;
    *qt = -2.0 / ck1 * (*qs);
    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::bool_to_int_with_if,
    clippy::assign_op_pattern
)]
pub(crate) fn rmn1(
    m: i32,
    n: i32,
    c: f64,
    x: f64,
    kd: i32,
    df: &[f64],
    r1f: &mut f64,
    r1d: &mut f64,
) -> Status {
    // =======================================================
    // Purpose: Compute prolate and oblate spheroidal radial
    //          functions of the first kind for given m, n,
    //          c and x
    //
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //
    // Routines called:
    //      (1) SCKB for computing expansion coefficients c2k
    //      (2) SPHJ for computing the spherical Bessel
    //          functions of the first kind
    // =======================================================

    let (Some(mut ck), Some(mut dj), Some(mut sj)) =
        (try_zeroed(200), try_zeroed(252), try_zeroed(252))
    else {
        return Status::NoMemory;
    };
    let eps = 1.0e-14;

    let nm1 = (n - m) / 2;
    let ip = if n - m == 2 * nm1 { 0 } else { 1 };
    // NOTE: `(int)c` is UB in C++ for NaN or |c| >= 2^31 (g++ on x86-64 gives INT_MIN, so nm < 0 and
    // all sums below are empty); `as i32` saturates (NaN -> 0, +inf -> i32::MAX, after which the
    // addition overflows).
    let nm = 25 + nm1 + c as i32;
    let reg = if m + nm > 80 { 1.0e-200 } else { 1.0 };
    let mut r0 = reg;

    for j in 1..=2 * m + ip {
        r0 *= f64::from(j);
    }
    let mut r = r0;
    // NOTE: here and below, `df[k - 1]` for k <= nm is out of bounds if nm > df.len() (200 in
    // rswfp/rswfo) and the sum hasn't converged yet; from rswfp/rswfo, sdmn already indexes out of
    // bounds (df[nm] with an nm that is >= this nm) in that case.
    let mut suc = r * df[0];
    let mut sw = 0.0;

    for k in 2..=nm {
        r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
            / (f64::from(k) - 1.0)
            / (f64::from(k + ip) - 1.5);
        suc += r * df[(k - 1) as usize];
        if (k > nm1) && ((suc - sw).abs() < suc.abs() * eps) {
            break;
        }
        sw = suc;
    }

    if x == 0.0 {
        sckb(m, n, c, df, &mut ck);

        let mut sum = 0.0;
        let mut sw1 = 0.0;
        // NOTE: `ck[j - 1]` is out of bounds if nm > 200 (sckb, which writes ck[0..nm'] with nm' >= nm
        // for c >= 0, already does that first).
        for j in 1..=nm {
            sum += ck[(j - 1) as usize];
            if (sum - sw1).abs() < sum.abs() * eps {
                break;
            }
            sw1 = sum;
        }

        let mut r1 = 1.0;
        for j in 1..=(n + m + ip) / 2 {
            r1 = r1 * (f64::from(j) + 0.5 * f64::from(n + m + ip));
        }

        let mut r2 = 1.0;
        for j in 1..=m {
            r2 *= 2.0 * c * f64::from(j);
        }

        let mut r3 = 1.0;
        for j in 1..=(n - m - ip) / 2 {
            r3 *= f64::from(j);
        }
        let sa0 = (2.0 * f64::from(m + ip) + 1.0) * r1
            / (f64::from(n).exp2() * c.powf(f64::from(ip)) * r2 * r3);

        if ip == 0 {
            *r1f = sum / (sa0 * suc) * df[0] * reg;
            *r1d = 0.0;
        } else if ip == 1 {
            *r1f = 0.0;
            *r1d = sum / (sa0 * suc) * df[0] * reg;
        }
        return Status::Ok;
    }

    let cx = c * x;
    let mut nm2 = 2 * nm + m;
    // NOTE: sphj writes sj/dj up to index nm2 (or less if msta1(cx, 200) < nm2), which is out of
    // bounds of the 252-element buffers if nm2 = 2 * (25 + (n - m) / 2 + (int)c) + m >= 252.
    sphj(cx, nm2, &mut nm2, &mut sj, &mut dj);

    let a0 = (1.0 - f64::from(kd) / (x * x)).powf(0.5 * f64::from(m)) / suc;
    *r1f = 0.0;
    sw = 0.0;
    // lg = 0;  (dead store in the C++)

    for k in 1..=nm {
        let l = 2 * k + m - n - 2 + ip;
        let lg = if l % 4 == 0 { 1 } else { -1 };

        if k == 1 {
            r = r0;
        } else {
            r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
                / (f64::from(k) - 1.0)
                / (f64::from(k + ip) - 1.5);
        }

        let np = m + 2 * k - 2 + ip;
        // NOTE: `sj[np]` (and `dj[np]` below) is out of bounds (C++ UB, panics here) once np >= 252
        // before convergence; since convergence is only checked for k > nm1 and np == n at
        // k == nm1 + 1, this always happens for n >= 252 (if sphj didn't overflow first).
        *r1f += f64::from(lg) * r * df[(k - 1) as usize] * sj[np as usize];

        if (k > nm1) && ((*r1f - sw).abs() < r1f.abs() * eps) {
            break;
        }
        sw = *r1f;
    }

    *r1f *= a0;
    let b0 = f64::from(kd * m) / x.powf(3.0) / (1.0 - f64::from(kd) / (x * x)) * (*r1f);

    let mut sud = 0.0;
    sw = 0.0;

    for k in 1..=nm {
        let l = 2 * k + m - n - 2 + ip;
        let lg = if l % 4 == 0 { 1 } else { -1 };

        if k == 1 {
            r = r0;
        } else {
            r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
                / (f64::from(k) - 1.0)
                / (f64::from(k + ip) - 1.5);
        }

        let np = m + 2 * k - 2 + ip;
        sud = sud + f64::from(lg) * r * df[(k - 1) as usize] * dj[np as usize];

        if (k > nm1) && ((sud - sw).abs() < sud.abs() * eps) {
            break;
        }
        sw = sud;
    }
    *r1d = b0 + a0 * c * sud;
    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
pub(crate) fn rmn2l(
    m: i32,
    n: i32,
    c: f64,
    x: f64,
    kd: i32,
    df: &[f64],
    r2f: &mut f64,
    r2d: &mut f64,
    id: &mut i32,
) -> Status {
    // ========================================================
    // Purpose: Compute prolate and oblate spheroidal radial
    //          functions of the second kind for given m, n,
    //          c and a large cx
    //
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //          Status::Other
    //              An internal convergence check failed.  When
    //              this happens, *Id is set to 10 on return.
    //
    // Routine called:
    //          SPHY for computing the spherical Bessel
    //          functions of the second kind
    // ========================================================

    let eps = 1.0e-14;

    let (Some(mut sy), Some(mut dy)) = (try_zeroed(252), try_zeroed(252)) else {
        return Status::NoMemory;
    };

    let mut ip = 1;
    let nm1 = (n - m) / 2;
    if n - m == 2 * nm1 {
        ip = 0;
    }
    // NOTE: `(int)c` is UB in C++ for NaN or |c| >= 2^31 (g++ on x86-64 gives INT_MIN); `as i32`
    // saturates (NaN -> 0).
    let nm = 25 + nm1 + c as i32;
    let mut reg = 1.0;
    if m + nm > 80 {
        reg = 1.0e-200;
    }
    let mut nm2 = 2 * nm + m;
    let cx = c * x;
    // NOTE: sphy writes sy up to index nm2 (unless |sy[k]| >= 1e300 for some smaller k) and dy up to
    // nm2 - 1, which is out of bounds of the 252-element buffers if
    // nm2 = 2 * (25 + (n - m) / 2 + (int)c) + m >= 252.
    sphy(cx, nm2, &mut nm2, &mut sy, &mut dy);
    let mut r0 = reg;

    for j in 1..=2 * m + ip {
        r0 *= f64::from(j);
    }

    let mut r = r0;
    // NOTE: `df[k - 1]` for k <= nm is out of bounds if nm > df.len() (200 in rswfp/rswfo) before
    // convergence (sdmn already indexes out of bounds in that case).
    let mut suc = r * df[0];
    let mut sw = 0.0;

    for k in 2..=nm {
        r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
            / (f64::from(k) - 1.0)
            / (f64::from(k + ip) - 1.5);
        suc += r * df[(k - 1) as usize];
        if (k > nm1) && ((suc - sw).abs() < suc.abs() * eps) {
            break;
        }
        sw = suc;
    }

    let a0 = (1.0 - f64::from(kd) / (x * x)).powf(0.5 * f64::from(m)) / suc;
    *r2f = 0.0;
    let mut eps1 = 0.0;
    let mut np = 0;

    // NOTE: `sw` is not reset here (nor before the `sud` loop below), unlike in rmn1.
    for k in 1..=nm {
        let l = 2 * k + m - n - 2 + ip;
        let lg = if l % 4 == 0 { 1 } else { -1 };

        if k == 1 {
            r = r0;
        } else {
            r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
                / (f64::from(k) - 1.0)
                / (f64::from(k + ip) - 1.5);
        }

        np = m + 2 * k - 2 + ip;
        // NOTE: `sy[np]` (and `dy[np]` below) is out of bounds (C++ UB, panics here) once np >= 252
        // before convergence; since convergence is only checked for k > nm1 and np == n at
        // k == nm1 + 1, this always happens for n >= 252 (if sphy didn't overflow first).
        *r2f += f64::from(lg) * r * (df[(k - 1) as usize] * sy[np as usize]);
        eps1 = (*r2f - sw).abs();
        if k > nm1 && eps1 < r2f.abs() * eps {
            break;
        }
        sw = *r2f;
    }

    // NOTE: `(int)` of NaN or +-inf is UB in C++, which is reachable here, e.g. when `*R2f` overflows,
    // as for `c == 0`. Like g++ on x86-64, `int_x86` gives INT_MIN, so that rswfp/rswfo skip
    // rmn2sp/rmn2so, which would otherwise return NaN or 1e300.
    let id1 = int_x86((eps1 / r2f.abs() + eps).log10());
    *r2f *= a0;

    if np >= nm2 {
        // On page 584 of "Computation of Special functions" by Shanjie Zhang
        // and Jian-Ming Jin, there is a comment next to this code that says
        // this condition indicates that convergence is not achieved, so we
        // return Status::Other in addition to setting *Id = 10.  But the
        // functions that call rmn2l (namely rswfp and rswfo) will actually
        // look at the value of Id to check for convergence; the returned
        // Status value is only checked for the occurrence of Status::NoMemory.
        *id = 10;
        return Status::Other;
    }

    let b0 = f64::from(kd * m) / x.powf(3.0) / (1.0 - f64::from(kd) / (x * x)) * (*r2f);
    let mut sud = 0.0;
    let mut eps2 = 0.0;

    for k in 1..=nm {
        let l = 2 * k + m - n - 2 + ip;
        let lg = if l % 4 == 0 { 1 } else { -1 };

        if k == 1 {
            r = r0;
        } else {
            r = r * (f64::from(m + k) - 1.0) * (f64::from(m + k + ip) - 1.5)
                / (f64::from(k) - 1.0)
                / (f64::from(k + ip) - 1.5);
        }

        np = m + 2 * k - 2 + ip;
        sud += f64::from(lg) * r * (df[(k - 1) as usize] * dy[np as usize]);
        eps2 = (sud - sw).abs();
        if (k > nm1) && (eps2 < sud.abs() * eps) {
            break;
        }
        sw = sud;
    }

    *r2d = b0 + a0 * c * sud;
    // NOTE: same `(int)` of NaN or +-inf as for id1 above
    let id2 = int_x86((eps2 / sud.abs() + eps).log10());
    *id = if id1 > id2 { id1 } else { id2 };
    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp
)]
pub(crate) fn rmn2so(
    m: i32,
    n: i32,
    c: f64,
    x: f64,
    cv: f64,
    kd: i32,
    df: &[f64],
    r2f: &mut f64,
    r2d: &mut f64,
) -> Status {
    // =============================================================
    // Purpose: Compute oblate radial functions of the second kind
    //          with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix)
    //
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //
    // Routines called:
    //      (1) SCKB for computing the expansion coefficients c2k
    //      (2) KMN for computing the joining factors
    //      (3) QSTAR for computing the factor defined in (15.7.3)
    //      (4) CBK for computing the expansion coefficient
    //          defined in (15.7.6)
    //      (5) GMN for computing the function defined in (15.7.4)
    //      (6) RMN1 for computing the radial function of the first
    //          kind
    // =============================================================

    let mut ck1 = 0.0;
    let mut ck2 = 0.0;
    let mut r1f = 0.0;
    let mut r1d = 0.0;
    let mut qs = 0.0;
    let mut qt = 0.0;
    let mut gf = 0.0;
    let mut gd = 0.0;
    let eps = 1.0e-14;

    if df[0].abs() <= 1.0e-280 {
        *r2f = 1.0e+300;
        *r2d = 1.0e+300;
        return Status::Ok;
    }

    let (Some(mut bk), Some(mut ck), Some(mut dn)) =
        (try_zeroed(200), try_zeroed(200), try_zeroed(200))
    else {
        return Status::NoMemory;
    };

    // NOTE: `(int)` is UB in C++ for NaN or out-of-range `(n - m) / 2 + c` (g++ on x86-64 gives
    // INT_MIN); `as i32` saturates (NaN -> 0).
    let nm = 25 + (f64::from((n - m) / 2) + c) as i32;
    let ip = (n - m) % 2;
    // NOTE: out of bounds (C++ UB, panics here) for the 200-element buffers: sckb writes ck[0..nm']
    // with nm' = 25 + (int)(0.5 * (n - m) + max(c, 1e-10)) (OOB if nm' > 200); kmn writes
    // dn[0..nn] with nn = 25 + (int)(0.5 * (n - m) + c) + m (OOB if nn > 200, e.g. m >= 176 for
    // n == m, c < 1); qstar reads ck[0..=m] (OOB if m >= 200; ck[nm'..] read as 0); cbk reads ck up
    // to index nm' (OOB if nm' >= 200; ck[nm'] is the zero from the zero-initialization) and writes
    // bk[0..nm' - 2]; gmn reads bk[nm' - 2] and bk[nm' - 1], which cbk never writes (zeros here as
    // in the C++).
    sckb(m, n, c, df, &mut ck);
    if kmn(m, n, c, cv, kd, df, &mut dn, &mut ck1, &mut ck2) == Status::NoMemory {
        return Status::NoMemory;
    }
    if qstar(m, n, c, ck1, &ck, &mut qs, &mut qt) == Status::NoMemory {
        return Status::NoMemory;
    }
    if cbk(m, n, c, cv, qt, &ck, &mut bk) == Status::NoMemory {
        return Status::NoMemory;
    }

    if x == 0.0 {
        let mut sum = 0.0;
        let mut sw = 0.0;

        for j in 0..nm {
            sum += ck[j as usize];
            if (sum - sw).abs() < sum.abs() * eps {
                break;
            }
            sw = sum;
        }

        if ip == 0 {
            r1f = sum / ck1;
            *r2f = -0.5 * PI * qs * r1f;
            *r2d = qs * r1f + bk[0];
        } else {
            r1d = sum / ck1;
            *r2f = bk[0];
            *r2d = -0.5 * PI * qs * r1d;
        }
    } else {
        gmn(m, n, c, x, &bk, &mut gf, &mut gd);
        if rmn1(m, n, c, x, kd, df, &mut r1f, &mut r1d) == Status::NoMemory {
            return Status::NoMemory;
        }
        let h0 = x.atan() - 0.5 * PI;
        *r2f = qs * r1f * h0 + gf;
        *r2d = qs * (r1d * h0 + r1f / (1.0 + x * x)) + gd;
    }
    Status::Ok
}

#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::manual_midpoint
)]
pub(crate) fn rmn2sp(
    m: i32,
    n: i32,
    c: f64,
    x: f64,
    cv: f64,
    kd: i32,
    df: &[f64],
    r2f: &mut f64,
    r2d: &mut f64,
) -> Status {
    // ======================================================
    // Purpose: Compute prolate spheroidal radial function
    //          of the second kind with a small argument
    //
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //
    // Routines called:
    //      (1) LPMNS for computing the associated Legendre
    //          functions of the first kind
    //      (2) LQMNS for computing the associated Legendre
    //          functions of the second kind
    //      (3) KMN for computing expansion coefficients
    //          and joining factors
    // ======================================================

    // NOTE: in the C++, `ip`, `nm1`, `nm` and `nm2` are `T` (double), not `int`.
    let mut ck1 = 0.0;
    let mut ck2 = 0.0;
    // NOTE: `gb` is uninitialized in the C++, but always assigned before it's read.
    let mut gb = 0.0;

    // NOTE: these buffers are uninitialized in the C++ (`new T[n]` without `()`); see the NOTEs
    // at the `su2` loop below for the elements that can be read before being written.
    // fortran index start from 0
    let (Some(mut pm), Some(mut pd), Some(mut qm), Some(mut qd), Some(mut dn)) = (
        try_zeroed(252),
        try_zeroed(252),
        try_zeroed(252),
        try_zeroed(252),
        // fortran index start from 1
        try_zeroed(201),
    ) else {
        return Status::NoMemory;
    };

    let eps = 1.0e-14;

    let nm1 = f64::from((n - m) / 2);
    let nm = 25.0 + nm1 + c;
    let nm2 = 2.0 * nm + f64::from(m);
    let ip = f64::from((n - m) % 2);

    // NOTE: kmn writes dn[0..nn] with nn = 25 + (int)(0.5 * (n - m) + c) + m, which is out of bounds
    // of the 201-element `dn` (C++ UB, panics here) if nn > 201.
    if kmn(m, n, c, cv, kd, df, &mut dn, &mut ck1, &mut ck2) == Status::NoMemory {
        return Status::NoMemory;
    }
    // NOTE: implicit double -> int conversion of nm2 (UB in C++ if NaN/out of range; g++ on x86-64
    // gives INT_MIN, `as i32` gives 0 for NaN). lpmns/lqmns write pm/pd/qm/qd[0..=nm2i] with
    // nm2i = (int)nm2 = 50 + 2 * ((n - m) / 2) + (int)(2c) + m (c >= 0), out of bounds of the
    // 252-element buffers if nm2i >= 252 (e.g. n - m = 198 with m + (int)(2c) >= 4).
    lpmns(m, nm2 as i32, x, &mut pm, &mut pd);
    lqmns(m, nm2 as i32, x, &mut qm, &mut qd);

    let mut su0 = 0.0;
    let mut sw = 0.0;
    let mut k = 1;
    while f64::from(k) <= nm {
        let j = (f64::from(2 * k - 2 + m) + ip) as i32;
        su0 += df[(k - 1) as usize] * qm[j as usize];
        if (f64::from(k) > nm1) && ((su0 - sw).abs() < su0.abs() * eps) {
            break;
        }
        sw = su0;
        k += 1;
    }

    // NOTE: `sw` is not reset here (nor before the su2/sd2 loops below).
    let mut sd0 = 0.0;
    k = 1;
    while f64::from(k) <= nm {
        let j = (f64::from(2 * k - 2 + m) + ip) as i32;
        sd0 += df[(k - 1) as usize] * qd[j as usize];
        if f64::from(k) > nm1 && (sd0 - sw).abs() < sd0.abs() * eps {
            break;
        }
        sw = sd0;
        k += 1;
    }

    let mut su1 = 0.0;
    let mut sd1 = 0.0;
    for k in 1..=m {
        let mut j = (f64::from(m - 2 * k) + ip) as i32;
        if j < 0 {
            j = -j - 1;
        }
        su1 += dn[(k - 1) as usize] * qm[j as usize];
        sd1 += dn[(k - 1) as usize] * qd[j as usize];
    }

    let ga = ((x - 1.0) / (x + 1.0)).powf(0.5 * f64::from(m));
    for k in 1..=m {
        let mut j = (f64::from(m - 2 * k) + ip) as i32;
        if j >= 0 {
            continue;
        }
        if j < 0 {
            j = -j - 1;
        }
        let mut r1 = 1.0;
        for j1 in 0..j {
            r1 *= f64::from(m + j1);
        }
        let mut r2 = 1.0;
        for j2 in 1..=(m - j - 2) {
            r2 *= f64::from(j2);
        }
        let mut r3 = 1.0;
        let mut sf = 1.0;
        for l1 in 1..=j {
            r3 = 0.5 * r3 * (f64::from(-j + l1) - 1.0) * f64::from(j + l1)
                / f64::from((m + l1) * l1)
                * (1.0 - x);
            sf += r3;
        }
        if m - j >= 2 {
            gb = (f64::from(m - j) - 1.0) * r2;
        }
        if m - j <= 1 {
            gb = 1.0;
        }
        let spl = r1 * ga * gb * sf;
        su1 += (-1.0_f64).powf(f64::from(j + m)) * dn[(k - 1) as usize] * spl;
        let spd1 = f64::from(m) / (x * x - 1.0) * spl;
        let gc = 0.5 * f64::from(j) * (f64::from(j) + 1.0) / (f64::from(m) + 1.0);
        let mut sd = 1.0;
        let mut r4 = 1.0;
        // for (l1 = 1; l1 <= j - 1; l1++)
        for l1 in 1..j {
            r4 = 0.5 * r4 * f64::from(-j + l1) * (f64::from(j + l1) + 1.0)
                / ((f64::from(m + l1) + 1.0) * f64::from(l1))
                * (1.0 - x);
            sd += r4;
        }
        let spd2 = r1 * ga * gb * gc * sd;
        sd1 += (-1.0_f64).powf(f64::from(j + m)) * dn[(k - 1) as usize] * (spd1 + spd2);
    }
    let mut su2 = 0.0;
    let mut ki = ((f64::from(2 * m + 1) + ip) / 2.0) as i32;
    ki = ki.max(1);
    assert!((ki - 1) >= 0);
    // NOTE: implicit double -> int conversion (same NaN/out-of-range UB as for nm2 above).
    let nm3 = (nm + f64::from(ki)) as i32;
    // NOTE: if this loop doesn't break early (e.g. when all terms vanish, as for c == 0), it reads
    // dn[nm3 - 1] == dn[nn] (kmn writes only dn[0..nn]) for (ip == 0, m == 0) and for (ip == 1,
    // frac(c) < 0.5), and pm[nm2i + 1] (lpmns writes only pm[0..=nm2i]) for (ip == 0, m == 0,
    // frac(c) < 0.5). These are uninitialized reads in the C++ (here they read 0). They are out of
    // bounds if nm3 > 201 resp. nm2i + 1 >= 252.
    for k in ki..=nm3 {
        let j = (f64::from(2 * k - 1 - m) - ip) as i32;
        su2 += dn[(k - 1) as usize] * pm[j as usize];
        if (j > m) && ((su2 - sw).abs() < su2.abs() * eps) {
            break;
        }
        sw = su2;
    }
    let mut sd2 = 0.0;
    for k in ki..nm3 {
        let j = (f64::from(2 * k - 1 - m) - ip) as i32;
        sd2 += dn[(k - 1) as usize] * pd[j as usize];
        if j > m && (sd2 - sw).abs() < sd2.abs() * eps {
            break;
        }
        sw = sd2;
    }
    let sum = su0 + su1 + su2;
    let sdm = sd0 + sd1 + sd2;
    *r2f = sum / ck2;
    *r2d = sdm / ck2;
    Status::Ok
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn rswfp(
    m: i32,
    n: i32,
    c: f64,
    x: f64,
    cv: f64,
    kf: i32,
    r1f: &mut f64,
    r1d: &mut f64,
    r2f: &mut f64,
    r2d: &mut f64,
) -> Status {
    // ==============================================================
    // Purpose: Compute prolate spheriodal radial functions of the
    //          first and second kinds, and their derivatives
    // Input :  m  --- Mode parameter, m = 0,1,2,...
    //          n  --- Mode parameter, n = m,m+1,m+2,...
    //          c  --- Spheroidal parameter
    //          x  --- Argument of radial function ( x > 1.0 )
    //          cv --- Characteristic value
    //          KF --- Function code
    //                 KF=1 for the first kind
    //                 KF=2 for the second kind
    //                 KF=3 for both the first and second kinds
    // Output:  R1F --- Radial function of the first kind
    //          R1D --- Derivative of the radial function of
    //                  the first kind
    //          R2F --- Radial function of the second kind
    //          R2D --- Derivative of the radial function of
    //                  the second kind
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //
    // Routines called:
    //      (1) SDMN for computing expansion coefficients dk
    //      (2) RMN1 for computing prolate and oblate radial
    //          functions of the first kind
    //      (3) RMN2L for computing prolate and oblate radial
    //          functions of the second kind for a large argument
    //      (4) RMN2SP for computing the prolate radial function
    //          of the second kind for a small argument
    // ==============================================================

    // NOTE: `df` is uninitialized in the C++ (`new T[200]` without `()`). sdmn writes df[0..=nm]
    // (nm = 25 + (int)(0.5 * (n - m) + c); out of bounds if nm >= 200, e.g. n - m = 198 with
    // c >= 76), but for c < 1e-10 only df[0..nm]; then df[nm] can be read uninitialized by sckb (only
    // reached via rmn1 with x == 0, i.e. not from here) and, for negative c, by rmn1/rmn2l/rmn2sp.
    let Some(mut df) = try_zeroed(200) else {
        return Status::NoMemory;
    };
    // NOTE: `id` is uninitialized in the C++, but always written by rmn2l before it's read.
    let mut id = 0;
    let kd = 1;

    if sdmn(m, n, c, cv, kd, &mut df) == Status::NoMemory {
        return Status::NoMemory;
    }

    if kf != 2 && rmn1(m, n, c, x, kd, &df, r1f, r1d) == Status::NoMemory {
        return Status::NoMemory;
    }
    if kf > 1 {
        if rmn2l(m, n, c, x, kd, &df, r2f, r2d, &mut id) == Status::NoMemory {
            return Status::NoMemory;
        }
        if id > -8 && rmn2sp(m, n, c, x, cv, kd, &df, r2f, r2d) == Status::NoMemory {
            return Status::NoMemory;
        }
    }
    Status::Ok
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn rswfo(
    m: i32,
    n: i32,
    c: f64,
    x: f64,
    cv: f64,
    kf: i32,
    r1f: &mut f64,
    r1d: &mut f64,
    r2f: &mut f64,
    r2d: &mut f64,
) -> Status {
    // ==========================================================
    // Purpose: Compute oblate radial functions of the first
    //          and second kinds, and their derivatives
    // Input :  m  --- Mode parameter,  m = 0,1,2,...
    //          n  --- Mode parameter,  n = m,m+1,m+2,...
    //          c  --- Spheroidal parameter
    //          x  --- Argument (x ≥ 0)
    //          cv --- Characteristic value
    //          KF --- Function code
    //                 KF=1 for the first kind
    //                 KF=2 for the second kind
    //                 KF=3 for both the first and second kinds
    // Output:  R1F --- Radial function of the first kind
    //          R1D --- Derivative of the radial function of
    //                  the first kind
    //          R2F --- Radial function of the second kind
    //          R2D --- Derivative of the radial function of
    //                  the second kind
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    //
    // Routines called:
    //      (1) SDMN for computing expansion coefficients dk
    //      (2) RMN1 for computing prolate or oblate radial
    //          function of the first kind
    //      (3) RMN2L for computing prolate or oblate radial
    //          function of the second kind for a large argument
    //      (4) RMN2SO for computing oblate radial functions of
    //          the second kind for a small argument
    // ==========================================================

    // NOTE: `df` is uninitialized in the C++ (`new T[200]` without `()`). sdmn writes df[0..=nm]
    // (nm = 25 + (int)(0.5 * (n - m) + c); out of bounds if nm >= 200, e.g. n - m = 198 with
    // c >= 76), but for c < 1e-10 only df[0..nm]. In that case sckb (called by rmn2so, and by rmn1
    // when x == 0) reads df[nm] uninitialized whenever its inner sum doesn't converge, which happens
    // for c == 0 (e.g. m = n, c = 0, kf >= 2: for k >= 1 all df[i] are 0); here it reads 0.
    let Some(mut df) = try_zeroed(200) else {
        return Status::NoMemory;
    };
    let mut id;
    let kd = -1;

    if sdmn(m, n, c, cv, kd, &mut df) == Status::NoMemory {
        return Status::NoMemory;
    }

    if kf != 2 && rmn1(m, n, c, x, kd, &df, r1f, r1d) == Status::NoMemory {
        return Status::NoMemory;
    }
    if kf > 1 {
        id = 10;
        if x > 1e-8 && rmn2l(m, n, c, x, kd, &df, r2f, r2d, &mut id) == Status::NoMemory {
            return Status::NoMemory;
        }
        if id > -1 && rmn2so(m, n, c, x, cv, kd, &df, r2f, r2d) == Status::NoMemory {
            return Status::NoMemory;
        }
    }
    Status::Ok
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::similar_names
)]
pub(crate) fn sckb(m: i32, n: i32, mut c: f64, df: &[f64], ck: &mut [f64]) {
    // ======================================================
    // Purpose: Compute the expansion coefficients of the
    //          prolate and oblate spheroidal functions
    // Input :  m  --- Mode parameter
    //          n  --- Mode parameter
    //          c  --- Spheroidal parameter
    //          DF(k) --- Expansion coefficients dk
    // Output:  CK(k) --- Expansion coefficients ck;
    //                    CK(1), CK(2), ... correspond to
    //                    c0, c2, ...
    // ======================================================

    if c <= 1.0e-10 {
        c = 1.0e-10;
    }
    // NOTE: `(int)` is UB in C++ if `0.5 * (n - m) + c` is NaN or out of the int range; `as i32`
    // saturates (NaN -> 0).
    let nm = 25 + (0.5 * f64::from(n - m) + c) as i32;
    let ip = (n - m) % 2;
    let reg = if (m + nm) > 80 { 1.0e-200 } else { 1.0 };
    let mut fac = -(0.5_f64.powf(f64::from(m)));
    let mut sw = 0.0;

    // NOTE: df and ck have 200 elements in all callers, so ck[k] (k <= nm - 1) is out of bounds
    // (UB in C++) if nm > 200, and df[i] (i <= nm, only if the inner loop doesn't break before) if
    // nm >= 200, i.e. for 0.5 * (n - m) + c >= 175.
    for k in 0..nm {
        fac = -fac;
        let i1 = 2 * k + ip + 1;
        let mut r = reg;

        for i in i1..i1 + 2 * m {
            r *= f64::from(i);
        }
        let i2 = k + m + ip;
        for i in i2..i2 + k {
            r *= f64::from(i) + 0.5;
        }
        let mut sum = r * df[k as usize];
        for i in k + 1..=nm {
            let d1 = 2.0 * f64::from(i) + f64::from(ip);
            let d2 = 2.0 * f64::from(m) + d1;
            let d3 = f64::from(i + m + ip) - 0.5;
            r = r * d2 * (d2 - 1.0) * f64::from(i) * (d3 + f64::from(k))
                / (d1 * (d1 - 1.0) * f64::from(i - k) * d3);
            sum += r * df[i as usize];
            if (sw - sum).abs() < sum.abs() * 1.0e-14 {
                break;
            }
            sw = sum;
        }
        let mut r1 = reg;
        for i in 2..=m + k {
            r1 *= f64::from(i);
        }
        ck[k as usize] = fac * sum / r1;
    }
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines
)]
pub(crate) fn sdmn(m: i32, n: i32, c: f64, cv: f64, kd: i32, df: &mut [f64]) -> Status {
    // =====================================================
    // Purpose: Compute the expansion coefficients of the
    //          prolate and oblate spheroidal functions, dk
    // Input :  m  --- Mode parameter
    //          n  --- Mode parameter
    //          c  --- Spheroidal parameter
    //          cv --- Characteristic value
    //          KD --- Function code
    //                 KD=1 for prolate; KD=-1 for oblate
    // Output:  DF(k) --- Expansion coefficients dk;
    //                    DF(1), DF(2), ... correspond to
    //                    d0, d2, ... for even n-m and d1,
    //                    d3, ... for odd n-m
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    // =====================================================

    // NOTE: `(int)` is UB in C++ if `0.5 * (n - m) + c` is NaN or out of the int range (g++ on
    // x86-64 gives INT_MIN, so nm + 2 < 0 and the allocations below fail with Status::NoMemory).
    // Here `as i32` saturates: NaN -> 0 (nm = 25, so it computes NaNs with Status::Ok), and
    // c >= 2^31 -> i32::MAX, so `25 + ...` overflows (panics in debug; wraps in release, after
    // which nm + 2 < 0 and the allocations fail with Status::NoMemory, like in C++).
    let nm = 25 + (0.5 * f64::from(n - m) + c) as i32;

    if c < 1e-10 {
        // NOTE: df[nm] is not written in this branch (it is in the other one), and sckb can read
        // it; the callers rswfp and rswfo pass an uninitialized `new T[200]` buffer.
        for i in 1..=nm {
            df[(i - 1) as usize] = 0.0;
        }
        df[((n - m) / 2) as usize] = 1.0;
        return Status::Ok;
    }

    let a = try_zeroed((nm + 2) as usize);
    let d = try_zeroed((nm + 2) as usize);
    let g = try_zeroed((nm + 2) as usize);
    let (Some(mut a), Some(mut d), Some(mut g)) = (a, d, g) else {
        return Status::NoMemory;
    };
    let cs = c * c * f64::from(kd);
    let ip = (n - m) % 2;

    for i in 1..=nm + 2 {
        let k = if ip == 0 { 2 * (i - 1) } else { 2 * i - 1 };

        let dk0 = f64::from(m + k);
        let dk1 = f64::from(m + k + 1);
        let dk2 = f64::from(2 * (m + k));
        let d2k = f64::from(2 * m + k);

        a[(i - 1) as usize] = (d2k + 2.0) * (d2k + 1.0) / ((dk2 + 3.0) * (dk2 + 5.0)) * cs;
        d[(i - 1) as usize] = dk0 * dk1
            + (2.0 * dk0 * dk1 - 2.0 * f64::from(m) * f64::from(m) - 1.0)
                / ((dk2 - 1.0) * (dk2 + 3.0))
                * cs;
        g[(i - 1) as usize] =
            f64::from(k) * (f64::from(k) - 1.0) / ((dk2 - 3.0) * (dk2 - 1.0)) * cs;
    }

    let mut fs = 1.0;
    let mut f1 = 0.0;
    let mut f0 = 1.0e-100;
    let mut kb = 0;
    // NOTE: df has 200 elements in all callers, so df[nm] is out of bounds (UB in C++) if
    // nm >= 200, i.e. for 0.5 * (n - m) + c >= 175. (Can't happen for c < 1e-10 above, as long as
    // n - m <= 198; but the `*_cv` public functions don't check n - m.)
    df[nm as usize] = 0.0;
    let mut fl = 0.0;

    for k in (1..=nm).rev() {
        let mut f = -((d[k as usize] - cv) * f0 + a[k as usize] * f1) / g[k as usize];

        if f.abs() > df[k as usize].abs() {
            df[(k - 1) as usize] = f;
            f1 = f0;
            f0 = f;

            if f.abs() > 1.0e+100 {
                for k1 in k..=nm {
                    df[(k1 - 1) as usize] *= 1.0e-100;
                }
                f1 *= 1.0e-100;
                f0 *= 1.0e-100;
            }
        } else {
            kb = k;
            fl = df[k as usize];
            f1 = 1.0e-100;
            let mut f2 = -((d[0] - cv) / a[0]) * f1;
            df[0] = f1;

            if kb == 1 {
                fs = f2;
            } else if kb == 2 {
                df[1] = f2;
                fs = -((d[1] - cv) * f2 + g[1] * f1) / a[1];
            } else {
                df[1] = f2;

                for j in 3..=kb + 1 {
                    f = -((d[(j - 2) as usize] - cv) * f2 + g[(j - 2) as usize] * f1)
                        / a[(j - 2) as usize];
                    if j <= kb {
                        df[(j - 1) as usize] = f;
                    }
                    if f.abs() > 1.0e+100 {
                        for k1 in 1..=j {
                            df[(k1 - 1) as usize] *= 1.0e-100;
                        }
                        f *= 1.0e-100;
                        f2 *= 1.0e-100;
                    }
                    f1 = f2;
                    f2 = f;
                }
                fs = f;
            }
            break;
        }
    }

    // su1 = 0.0;  (dead store in C++)
    let mut r1 = 1.0;

    for j in m + ip + 1..=2 * (m + ip) {
        r1 *= f64::from(j);
    }
    let mut su1 = df[0] * r1;

    for k in 2..=kb {
        r1 = -r1 * (f64::from(k + m + ip) - 1.5) / (f64::from(k) - 1.0);
        su1 += r1 * df[(k - 1) as usize];
    }

    let mut su2 = 0.0;
    let mut sw = 0.0;

    for k in kb + 1..=nm {
        if k != 1 {
            r1 = -r1 * (f64::from(k + m + ip) - 1.5) / (f64::from(k) - 1.0);
        }
        su2 += r1 * df[(k - 1) as usize];

        if (sw - su2).abs() < su2.abs() * 1.0e-14 {
            break;
        }
        sw = su2;
    }
    let mut r3 = 1.0;

    for j in 1..=(m + n + ip) / 2 {
        r3 *= f64::from(j) + 0.5 * f64::from(n + m + ip);
    }
    let mut r4 = 1.0;

    for j in 1..=(n - m - ip) / 2 {
        r4 *= -4.0 * f64::from(j);
    }
    let s0 = r3 / (fl * (su1 / fs) + su2) / r4;

    for k in 1..=kb {
        df[(k - 1) as usize] *= fl / fs * s0;
    }
    for k in kb + 1..=nm {
        df[(k - 1) as usize] *= s0;
    }
    Status::Ok
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::float_cmp,
    clippy::collapsible_if,
    clippy::int_plus_one,
    clippy::manual_midpoint
)]
pub(crate) fn segv(m: i32, n: i32, c: f64, kd: i32, cv: &mut f64, eg: &mut [f64]) -> Status {
    // =========================================================
    // Purpose: Compute the characteristic values of spheroidal
    //          wave functions
    // Input :  m  --- Mode parameter
    //          n  --- Mode parameter
    //          c  --- Spheroidal parameter
    //          KD --- Function code
    //                 KD=1 for Prolate; KD=-1 for Oblate
    // Output:  CV --- Characteristic value for given m, n and c
    //          EG(L) --- Characteristic value for mode m and n'
    //                    ( L = n' - m + 1 )
    // Return value:
    //          Status::OK
    //              Normal return.
    //          Status::NoMemory
    //              An internal memory allocation failed.
    // =========================================================

    // eg[<=200] is supplied by the caller

    if c < 1e-10 {
        for i in 1..=(n - m + 1) {
            // Unlike in xsf, the product is computed in `double` instead of `int`, which overflows for
            // `n > 46340`; this is exact, since both factors are less than 2^31
            eg[(i - 1) as usize] = f64::from(i + m) * f64::from(i + m - 1);
        }
        *cv = eg[(n - m) as usize];
        return Status::Ok;
    }

    // TODO: Following array sizes should be decided dynamically
    let a = try_zeroed(300);
    let b = try_zeroed(100);
    let cv0 = try_zeroed(100);
    let d = try_zeroed(300);
    let e = try_zeroed(300);
    let f = try_zeroed(300);
    let g = try_zeroed(300);
    let h = try_zeroed(100);

    let (
        Some(mut a),
        Some(mut b),
        Some(mut cv0),
        Some(mut d),
        Some(mut e),
        Some(mut f),
        Some(mut g),
        Some(mut h),
    ) = (a, b, cv0, d, e, f, g, h)
    else {
        return Status::NoMemory;
    };
    // NOTE: b, cv0 and h have 100 elements, so b[icm - 1] etc. are out of bounds (UB in C++) if
    // icm > 100, i.e. n - m >= 200 (the public callers require n - m <= 198).
    let icm = (n - m + 2) / 2;
    // NOTE: `(int)` is UB in C++ if `0.5 * (n - m) + c` is NaN or out of the int range; `as i32`
    // saturates (NaN -> 0). For NaN c, g++ on x86-64 gives nm = 10 + INT_MIN, so d[nm - 1] below is
    // wildly out of bounds (in practice a segfault); here nm = 10 and everything is NaN, so the
    // bisection loop below never terminates (`fabs(NaN) < 1e-14` is never true).
    let nm = 10 + (0.5 * f64::from(n - m) + c) as i32;
    let cs = c * c * f64::from(kd);
    // k = 0;  (dead store in C++)
    // NOTE: a, d, e, f and g have 300 elements, so the accesses with index nm - 1 are out of bounds
    // (UB in C++) if nm > 300, i.e. for 0.5 * (n - m) + c >= 291.
    for l in 0..=1 {
        for i in 1..=nm {
            let k = if l == 0 { 2 * (i - 1) } else { 2 * i - 1 };
            let dk0 = f64::from(m + k);
            let dk1 = f64::from(m + k + 1);
            let dk2 = f64::from(2 * (m + k));
            let d2k = f64::from(2 * m + k);
            a[(i - 1) as usize] = (d2k + 2.0) * (d2k + 1.0) / ((dk2 + 3.0) * (dk2 + 5.0)) * cs;
            d[(i - 1) as usize] = dk0 * dk1
                + (2.0 * dk0 * dk1 - 2.0 * f64::from(m) * f64::from(m) - 1.0)
                    / ((dk2 - 1.0) * (dk2 + 3.0))
                    * cs;
            g[(i - 1) as usize] =
                f64::from(k) * (f64::from(k) - 1.0) / ((dk2 - 3.0) * (dk2 - 1.0)) * cs;
        }
        for k in 2..=nm {
            e[(k - 1) as usize] = (a[(k - 2) as usize] * g[(k - 1) as usize]).sqrt();
            f[(k - 1) as usize] = e[(k - 1) as usize] * e[(k - 1) as usize];
        }
        f[0] = 0.0;
        e[0] = 0.0;
        let mut xa = d[(nm - 1) as usize] + e[(nm - 1) as usize].abs();
        let mut xb = d[(nm - 1) as usize] - e[(nm - 1) as usize].abs();
        let nm1 = nm - 1;
        for i in 1..=nm1 {
            let t = e[(i - 1) as usize].abs() + e[i as usize].abs();
            let mut t1 = d[(i - 1) as usize] + t;
            if xa < t1 {
                xa = t1;
            }
            t1 = d[(i - 1) as usize] - t;
            if t1 < xb {
                xb = t1;
            }
        }
        for i in 1..=icm {
            b[(i - 1) as usize] = xa;
            h[(i - 1) as usize] = xb;
        }
        for k in 1..=icm {
            for k1 in k..=icm {
                if b[(k1 - 1) as usize] < b[(k - 1) as usize] {
                    b[(k - 1) as usize] = b[(k1 - 1) as usize];
                    break;
                }
            }
            if k != 1 {
                if h[(k - 1) as usize] < h[(k - 2) as usize] {
                    h[(k - 1) as usize] = h[(k - 2) as usize];
                }
            }
            let mut x1;
            loop {
                x1 = (b[(k - 1) as usize] + h[(k - 1) as usize]) / 2.0;
                cv0[(k - 1) as usize] = x1;
                if ((b[(k - 1) as usize] - h[(k - 1) as usize]) / x1).abs() < 1e-14 {
                    break;
                }
                let mut j = 0;
                let mut s = 1.0;
                for i in 1..=nm {
                    if s == 0.0 {
                        s += 1e-30;
                    }
                    let t = f[(i - 1) as usize] / s;
                    s = d[(i - 1) as usize] - t - x1;
                    if s < 0.0 {
                        j += 1;
                    }
                }
                if j < k {
                    h[(k - 1) as usize] = x1;
                } else {
                    b[(k - 1) as usize] = x1;
                    if j >= icm {
                        b[(icm - 1) as usize] = x1;
                    } else {
                        if h[j as usize] < x1 {
                            h[j as usize] = x1;
                        }
                        if x1 < b[(j - 1) as usize] {
                            b[(j - 1) as usize] = x1;
                        }
                    }
                }
            }
            cv0[(k - 1) as usize] = x1;
            let mut index = 0;
            if l == 0 {
                index = 2 * k - 2;
            }
            if l == 1 {
                index = 2 * k - 1;
            }
            // boundary check
            if index >= n - m + 1 {
                break;
            }
            eg[index as usize] = cv0[(k - 1) as usize];
        }
    }
    *cv = eg[(n - m) as usize];
    Status::Ok
}
