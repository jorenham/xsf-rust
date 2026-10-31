//! Bernoulli and Euler numbers, translated into pure Rust from `xsf/specfun/specfun.h`
//! (xsf v0.2.2).

use core::f64::consts::PI;

#[allow(
    clippy::approx_constant,
    clippy::cast_sign_loss,
    clippy::excessive_precision,
    clippy::unreadable_literal
)]
pub(crate) fn bernob(n: i32, bn: &mut [f64]) {
    // ======================================
    // Purpose: Compute Bernoulli number Bn
    // Input :  n >= 3 --- Serial number
    // Output:  BN(n)  --- Bn
    // ======================================

    // NOTE: writes bn[0..=2] and bn[m] for even 4 <= m <= n, so it needs
    // bn.len() >= max(3, n + 1); with n < 2 and an (n + 1)-element buffer, the write to bn[2]
    // (and bn[1] for n = 0) is out of bounds. The only caller, `bernoulli::<N>`, passes n = N - 1
    // with N >= 3 and an N-element buffer, so that's fine. (`m < (n + 1)` is written as
    // `4..=n`, which is the same except that `n + 1` would overflow for n = i32::MAX.)

    let tpi = 6.283185307179586;

    bn[0] = 1.0;
    bn[1] = -0.5;
    bn[2] = 1.0 / 6.0;
    // pow(2.0 / tpi, 2)
    let mut r1 = (2.0 / tpi) * (2.0 / tpi);
    // C++: for (m = 4; m < (n + 1); m += 2)
    for m in (4..=n).step_by(2) {
        r1 = -r1 * f64::from(m - 1) * f64::from(m) / (tpi * tpi);
        let mut r2 = 1.0;
        for k in 2..10001 {
            let s = (1.0 / f64::from(k)).powf(f64::from(m));
            r2 += s;
            if s < 1e-15 {
                break;
            }
        }
        bn[m as usize] = r1 * r2;
    }
}

#[allow(clippy::cast_sign_loss)]
pub(crate) fn eulerb(n: i32, en: &mut [f64]) {
    // ======================================
    // Purpose: Compute Euler number En
    // Input :  n --- Serial number
    // Output:  EN(n) --- En
    // ======================================

    // NOTE: writes en[0], en[2] and en[m] for even 4 <= m <= n (en[1] and the odd entries are
    // never written), so it needs en.len() >= max(3, n + 1); with n < 2 and an (n + 1)-element
    // buffer, the write to en[2] is out of bounds. The only caller, `euler::<N>`, passes
    // n = N - 1 with N >= 3 and a zeroed N-element buffer, so that's fine.

    let hpi = 2.0 / PI;
    en[0] = 1.0;
    en[2] = -1.0;
    let mut r1 = -4.0 * hpi.powf(3.0);
    for m in (4..=n).step_by(2) {
        r1 = -r1 * f64::from(m - 1) * f64::from(m) * hpi * hpi;
        let mut r2 = 1.0;
        let mut isgn = 1;
        for k in (3..=1000).step_by(2) {
            isgn = -isgn;
            let s = (1.0 / f64::from(k)).powf(f64::from(m + 1));
            r2 += f64::from(isgn) * s;
            if s < 1e-15 {
                break;
            }
        }
        en[m as usize] = r1 * r2;
    }
}
