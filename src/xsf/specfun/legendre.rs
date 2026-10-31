//! Associated Legendre functions, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.

use core::f64::consts::PI;

use super::{gamma2, psi_spec};
use crate::xsf::tools::int_x86;

#[allow(clippy::float_cmp, clippy::cast_sign_loss)]
pub(crate) fn lpmns(m: i32, n: i32, x: f64, pm: &mut [f64], pd: &mut [f64]) {
    // ========================================================
    // Purpose: Compute associated Legendre functions Pmn(x)
    //          and Pmn'(x) for a given order
    // Input :  x --- Argument of Pmn(x)
    //          m --- Order of Pmn(x),  m = 0,1,2,...,n
    //          n --- Degree of Pmn(x), n = 0,1,2,...,N
    // Output:  PM(n) --- Pmn(x)
    //          PD(n) --- Pmn'(x)
    // ========================================================

    // NOTE: only `pm[0..=n]` / `pd[0..=n]` belong to the output, but (for |x| != 1) `pm[m]` and
    // `pm[m + 1]` are always written and `pm[1]` is always read, so this indexes out of bounds
    // when m + 1 > n (with an `n + 1`-element buffer). The only caller, rmn2sp, passes
    // n = (int)(2 * (25.0 + (n - m) / 2 + c) + m) >= m + 50 into 252-element buffers, which
    // overflow when that n exceeds 251, i.e. when
    // n_caller - (n_caller - m) % 2 + (int)(2 * c) >= 202, e.g. prolate_radial2_nocv(m, n, c, x)
    // with (0, 198, 2, 1.0001) or (2, 200, 1, 1.5) (confirmed with ASan), which the public API
    // allows.
    for k in 0..=n {
        pm[k as usize] = 0.0;
        pd[k as usize] = 0.0;
    }

    if x.abs() == 1.0 {
        for k in 0..=n {
            if m == 0 {
                pm[k as usize] = 1.0;
                pd[k as usize] = 0.5 * f64::from(k) * (f64::from(k) + 1.0);
                if x < 0.0 {
                    pm[k as usize] *= f64::from(if (k % 2) == 0 { 1 } else { -1 });
                    pd[k as usize] *= f64::from(if ((k + 1) % 2) == 0 { 1 } else { -1 });
                }
            } else if m == 1 {
                pd[k as usize] = 1e300;
            } else if m == 2 {
                pd[k as usize] = -0.25
                    * (f64::from(k) + 2.0)
                    * (f64::from(k) + 1.0)
                    * f64::from(k)
                    * (f64::from(k) - 1.0);
                if x < 0.0 {
                    pd[k as usize] *= f64::from(if ((k + 1) % 2) == 0 { 1 } else { -1 });
                }
            }
        }
        return;
    }

    let x0 = (1.0 - x * x).abs();
    let mut pm0 = 1.0;
    let mut pmk = pm0;
    for k in 1..=m {
        pmk = (2.0 * f64::from(k) - 1.0) * x0.sqrt() * pm0;
        pm0 = pmk;
    }
    let mut pm1 = (2.0 * f64::from(m) + 1.0) * x * pm0;
    pm[m as usize] = pmk;
    pm[(m + 1) as usize] = pm1;
    for k in m + 2..=n {
        let pm2 = ((2.0 * f64::from(k) - 1.0) * x * pm1 - (f64::from(k + m) - 1.0) * pmk)
            / f64::from(k - m);
        pm[k as usize] = pm2;
        pmk = pm1;
        pm1 = pm2;
    }

    pd[0] = ((1.0 - f64::from(m)) * pm[1] - x * pm[0]) / (x * x - 1.0);
    for k in 1..=n {
        pd[k as usize] = (f64::from(k) * x * pm[k as usize]
            - f64::from(k + m) * pm[(k - 1) as usize])
            / (x * x - 1.0);
    }
    let coef = f64::from(if (m % 2) == 0 { 1 } else { -1 });
    for k in 1..=n {
        pm[k as usize] *= coef;
        pd[k as usize] *= coef;
    }
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
pub(crate) fn lpmv(x: f64, m: i32, v: f64) -> f64 {
    // =======================================================
    // Purpose: Compute the associated Legendre function
    //          Pmv(x) with an integer order and an arbitrary
    //          degree v, using recursion for large degrees
    // Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
    //          m   --- Order of Pmv(x)
    //          v   --- Degree of Pmv(x)
    // Output:  PMV --- Pmv(x)
    // Routine called:  LPMV0
    // =======================================================

    let mut pmv;
    if (x == -1.0) && (v != f64::from(v as i32)) {
        if m == 0 {
            let vs = (v * PI).sin() / (v * PI);
            pmv = if vs > 0.0 { -1e300 } else { 1e300 };
        } else {
            pmv = 1e300;
        }
        return pmv;
    }
    let mut vx = v;
    let mut mx = m;
    // DLMF 14.9.5
    if v < 0.0 {
        vx = -vx - 1.0;
    }
    let mut neg_m = 0;
    if m < 0 {
        if ((vx + f64::from(m) + 1.0) > 0.0) || (vx != f64::from(vx as i32)) {
            neg_m = 1;
            // Unlike in xsf, NaN for m = INT_MIN, where `-m` overflows
            if m == i32::MIN {
                return f64::NAN;
            }
            mx = -m;
        } else {
            // We don't handle cases where DLMF 14.9.3 doesn't help
            return f64::NAN;
        }
    }
    // NOTE: `(int)vx` is UB for vx >= 2^31; like g++ on x86-64, `int_x86` gives INT_MIN, so that
    // lpmv0 is called directly (instead of up-recursion with ~2^31 steps).
    let nv = int_x86(vx);
    let v0 = vx - f64::from(nv);
    if (nv > 2) && (nv > mx) {
        // Up-recursion on degree, AMS 8.5.3 / DLMF 14.10.3
        let mut p0 = lpmv0(v0 + f64::from(mx), mx, x);
        let mut p1 = lpmv0(v0 + f64::from(mx) + 1.0, mx, x);
        pmv = p1;
        // NOTE: `mx + 2` overflows for mx = INT_MAX - 1 (with nv = INT_MAX), where the loop is
        // empty here. In C++, `j <= nv` with nv = INT_MAX can never be false (UB on `j++`), but
        // g++ -O2 terminates the loop like the Rust range does.
        for j in mx.saturating_add(2)..=nv {
            pmv = ((2.0 * (v0 + f64::from(j)) - 1.0) * x * p1
                - (v0 + f64::from(j) - 1.0 + f64::from(mx)) * p0)
                / (v0 + f64::from(j) - f64::from(mx));
            p0 = p1;
            p1 = pmv;
        }
    } else {
        pmv = lpmv0(vx, mx, x);
    }
    if (neg_m != 0) && (pmv.abs() < 1.0e300) {
        // DLMF 14.9.3
        let g1 = gamma2(vx - f64::from(mx) + 1.0);
        let g2 = gamma2(vx + f64::from(mx) + 1.0);
        pmv = pmv * g1 / g2 * (-1.0_f64).powf(f64::from(mx));
    }
    pmv
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::approx_constant,
    clippy::assign_op_pattern,
    clippy::manual_midpoint,
    clippy::unreadable_literal
)]
pub(crate) fn lpmv0(v: f64, m: i32, x: f64) -> f64 {
    // =======================================================
    // Purpose: Compute the associated Legendre function
    //          Pmv(x) with an integer order and an arbitrary
    //          nonnegative degree v
    // Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
    //          m   --- Order of Pmv(x)
    //          v   --- Degree of Pmv(x)
    // Output:  PMV --- Pmv(x)
    // Routine called:  PSI_SPEC for computing Psi function
    // =======================================================

    // Unlike in xsf, the `int` products and sums that would overflow for large orders `m` (e.g.
    // `j * j` for m > 46340) are computed in `double`, which is exact where they don't overflow.
    // NOTE: `(int)v` is UB for v >= 2^31 (x86: INT_MIN; Rust saturates), see lpmv.

    let el = 0.5772156649015329;
    let eps = 1e-14;

    let nv = v as i32;
    let v0 = v - f64::from(nv);

    if x == -1.0 && v != f64::from(nv) {
        if m == 0 {
            return -1.0e+300;
        }
        if m != 0 {
            return 1.0e+300;
        }
    }

    let mut c0 = 1.0;
    if m != 0 {
        let mut rg = v * (v + f64::from(m));
        for j in 1..m {
            rg *= v * v - f64::from(j) * f64::from(j);
        }
        let xq = (1.0 - x * x).sqrt();
        let mut r0 = 1.0;
        for j in 1..=m {
            r0 = 0.5 * r0 * xq / f64::from(j);
        }
        c0 = r0 * rg;
    }

    if v0 == 0.0 {
        // DLMF 14.3.4, 14.7.17, 15.2.4
        let mut pmv = 1.0;
        let mut r = 1.0;
        for k in 1..=(nv - m) {
            r = 0.5
                * r
                * (-f64::from(nv) + f64::from(m) + f64::from(k) - 1.0)
                * (f64::from(nv) + f64::from(m) + f64::from(k))
                / (f64::from(k) * (f64::from(k) + f64::from(m)))
                * (1.0 + x);
            pmv += r;
        }
        (-1.0_f64).powf(f64::from(nv)) * c0 * pmv
    } else if x >= -0.35 {
        // DLMF 14.3.4, 15.2.1
        let mut pmv = 1.0;
        let mut r = 1.0;
        for k in 1..=100 {
            r = 0.5
                * r
                * (-v + f64::from(m) + f64::from(k) - 1.0)
                * (v + f64::from(m) + f64::from(k))
                / (f64::from(k) * (f64::from(m) + f64::from(k)))
                * (1.0 - x);
            pmv += r;
            if k > 12 && (r / pmv).abs() < eps {
                break;
            }
        }
        (-1.0_f64).powf(f64::from(m)) * c0 * pmv
    } else {
        // DLMF 14.3.5, 15.8.10
        let vs = (v * PI).sin() / PI;
        let mut pv0 = 0.0;
        if m != 0 {
            let qr = ((1.0 - x) / (1.0 + x)).sqrt();
            let mut r2 = 1.0;
            for j in 1..=m {
                r2 *= qr * f64::from(j);
            }
            let mut s0 = 1.0;
            let mut r1 = 1.0;
            for k in 1..m {
                r1 = 0.5 * r1 * (-v + f64::from(k) - 1.0) * (v + f64::from(k))
                    / (f64::from(k) * (f64::from(k) - f64::from(m)))
                    * (1.0 + x);
                s0 += r1;
            }
            pv0 = -vs * r2 / f64::from(m) * s0;
        }

        let pa = 2.0 * (psi_spec(v) + el) + PI / (PI * v).tan() + 1.0 / v;
        let mut s1 = 0.0;
        for j in 1..=m {
            s1 += (f64::from(j) * f64::from(j) + v * v)
                / (f64::from(j) * (f64::from(j) * f64::from(j) - v * v));
        }
        let mut pmv = pa + s1 - 1.0 / (f64::from(m) - v) + (0.5 * (1.0 + x)).ln();
        let mut r = 1.0;
        for k in 1..=100 {
            r = 0.5
                * r
                * (-v + f64::from(m) + f64::from(k) - 1.0)
                * (v + f64::from(m) + f64::from(k))
                / (f64::from(k) * (f64::from(k) + f64::from(m)))
                * (1.0 + x);
            let mut s = 0.0;
            for j in 1..=m {
                let kj = f64::from(k) + f64::from(j);
                s += (kj * kj + v * v) / (kj * (kj * kj - v * v));
            }

            let mut s2 = 0.0;
            for j in 1..=k {
                s2 = s2 + 1.0 / (f64::from(j) * (f64::from(j) * f64::from(j) - v * v));
            }

            let pss = pa + s + 2.0 * v * v * s2 - 1.0 / (f64::from(m) + f64::from(k) - v)
                + (0.5 * (1.0 + x)).ln();
            let r2 = pss * r;
            pmv += r2;
            if (r2 / pmv).abs() < eps {
                break;
            }
        }
        pv0 + pmv * vs * c0
    }
}

#[allow(
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::assign_op_pattern
)]
pub(crate) fn lqmns(m: i32, n: i32, x: f64, qm: &mut [f64], qd: &mut [f64]) {
    // ========================================================
    // Purpose: Compute associated Legendre functions Qmn(x)
    //          and Qmn'(x) for a given order
    // Input :  x --- Argument of Qmn(x)
    //          m --- Order of Qmn(x),  m = 0,1,2,...
    //          n --- Degree of Qmn(x), n = 0,1,2,...
    // Output:  QM(n) --- Qmn(x)
    //          QD(n) --- Qmn'(x)
    // ========================================================

    // NOTE: only `qm[0..=n]` / `qd[0..=n]` belong to the output, but `qm[1]` is always read (and,
    // for |x| < 1.0001, written), so this indexes out of bounds when n == 0 (with a 1-element
    // buffer); for 1.0001 <= |x| and n == 0, `qm[1]` is read without having been written. The
    // only caller, rmn2sp, has the same 252-element overflow as for lpmns (n > 251).
    let mut val = 0.0;
    if x.abs() == 1.0 {
        val = 1e300;
    }
    for k in 0..=n {
        qm[k as usize] = val;
        qd[k as usize] = val;
    }

    if x.abs() == 1.0 {
        return;
    }
    let ls: i32 = if x.abs() > 1.0 { -1 } else { 1 };

    let xq = (f64::from(ls) * (1.0 - x * x)).sqrt();
    let q0 = 0.5 * ((x + 1.0) / (x - 1.0)).abs().ln();
    let q00 = q0;
    let q10 = -1.0 / xq;
    let q01 = x * q0 - 1.0;
    let q11 = f64::from(-ls) * xq * (q0 + x / (1.0 - x * x));
    let mut qf0 = q00;
    let mut qf1 = q10;
    let mut qm0 = 0.0;
    let mut qm1 = 0.0;

    for k in 2..=m {
        qm0 = -2.0 * (f64::from(k) - 1.0) / xq * x * qf1
            - f64::from(ls) * (f64::from(k) - 1.0) * (2.0 - f64::from(k)) * qf0;
        qf0 = qf1;
        qf1 = qm0;
    }

    if m == 0 {
        qm0 = q00;
    }
    if m == 1 {
        qm0 = q10;
    }

    qm[0] = qm0;

    if x.abs() < 1.0001 {
        if (m == 0) && (n > 0) {
            qf0 = q00;
            qf1 = q01;
            for k in 2..=n {
                let qf2 = ((2.0 * f64::from(k) - 1.0) * x * qf1 - (f64::from(k) - 1.0) * qf0)
                    / f64::from(k);
                qm[k as usize] = qf2;
                qf0 = qf1;
                qf1 = qf2;
            }
        }

        let mut qg0 = q01;
        let mut qg1 = q11;
        for k in 2..=m {
            qm1 = -2.0 * (f64::from(k) - 1.0) / xq * x * qg1
                - f64::from(ls * k) * (3.0 - f64::from(k)) * qg0;
            qg0 = qg1;
            qg1 = qm1;
        }

        if m == 0 {
            qm1 = q01;
        }
        if m == 1 {
            qm1 = q11;
        }

        qm[1] = qm1;

        if (m == 1) && (n > 1) {
            let mut qh0 = q10;
            let mut qh1 = q11;
            for k in 2..=n {
                let qh2 = ((2.0 * f64::from(k) - 1.0) * x * qh1 - f64::from(k) * qh0)
                    / (f64::from(k) - 1.0);
                qm[k as usize] = qh2;
                qh0 = qh1;
                qh1 = qh2;
            }
        } else if m >= 2 {
            qg0 = q00;
            qg1 = q01;
            let mut qh0 = q10;
            let mut qh1 = q11;
            let mut qmk = 0.0;
            for l in 2..=n {
                let q0l = ((2.0 * f64::from(l) - 1.0) * x * qg1 - (f64::from(l) - 1.0) * qg0)
                    / f64::from(l);
                let q1l = ((2.0 * f64::from(l) - 1.0) * x * qh1 - f64::from(l) * qh0)
                    / (f64::from(l) - 1.0);
                qf0 = q0l;
                qf1 = q1l;
                for k in 2..=m {
                    qmk = -2.0 * (f64::from(k) - 1.0) / xq * x * qf1
                        - f64::from(ls)
                            * (f64::from(k + l) - 1.0)
                            * (f64::from(l) + 2.0 - f64::from(k))
                            * qf0;
                    qf0 = qf1;
                    qf1 = qmk;
                }
                qm[l as usize] = qmk;
                qg0 = qg1;
                qg1 = q0l;
                qh0 = qh1;
                qh1 = q1l;
            }
        }
    } else {
        // NOTE: for -1.1 <= x <= -1.0001 (or x NaN), `log(x - 1.0)` is NaN and `(int)NaN` is UB in
        // C++ (g++ on x86-64 yields INT_MIN, after which `(40 + m + n) * INT_MIN` is signed
        // overflow UB as well; in practice km is 0 or INT_MIN, so the loop below runs once or not
        // at all). Rust's `as i32` gives 0, so the loop runs once (k = 0). Not reachable with a
        // finite x from rmn2sp (x > 1), but a NaN x passes the public `x <= 1.0` check.
        let km = if x.abs() > 1.1 {
            40 + m + n
        } else {
            (40 + m + n) * (-1.0 - 1.8 * (x - 1.0).ln()) as i32
        };
        let mut qf2 = 0.0;
        qf1 = 1.0;
        for k in (0..=km).rev() {
            qf0 = ((2.0 * f64::from(k) + 3.0) * x * qf1
                - (f64::from(k) + 2.0 - f64::from(m)) * qf2)
                / (f64::from(k + m) + 1.0);
            if k <= n {
                qm[k as usize] = qf0;
            }
            qf2 = qf1;
            qf1 = qf0;
        }
        for k in 0..=n {
            qm[k as usize] = qm[k as usize] * qm0 / qf0;
        }
    }

    if x.abs() < 1.0 {
        for k in 0..=n {
            qm[k as usize] = (-1.0_f64).powf(f64::from(m)) * qm[k as usize];
        }
    }

    qd[0] = ((1.0 - f64::from(m)) * qm[1] - x * qm[0]) / (x * x - 1.0);
    for k in 1..=n {
        qd[k as usize] = (f64::from(k) * x * qm[k as usize]
            - f64::from(k + m) * qm[(k - 1) as usize])
            / (x * x - 1.0);
    }
}
