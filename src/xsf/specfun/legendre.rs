//! Associated Legendre functions, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).
//!
//! Unlike in xsf, the `double` -> `int` conversions saturate instead of being undefined behavior.

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
