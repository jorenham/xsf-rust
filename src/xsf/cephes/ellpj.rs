//! Translated into pure Rust from `xsf/cephes/ellpj.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Scipy changes:
//! - 07-18-2016: improve evaluation of dn near quarter periods
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `(sn, cn, dn, phi)` is returned instead of written to output pointers, and the
//! (unused) status code is dropped.

use crate::utils::MACHEP;
use core::f64::consts::{FRAC_PI_2, PI};

/// Jacobian Elliptic Functions
///
/// Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
/// and dn(u|m) of parameter m less than or equal to 1, and real
/// argument u.
///
/// These functions are periodic, with quarter-period on the
/// real axis equal to the complete elliptic integral
/// ellpk(m).
///
/// Relation to incomplete elliptic integral:
/// If u = ellik(phi,m), then sn(u|m) = sin(phi),
/// and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
///
/// Computation is by means of the arithmetic-geometric mean
/// algorithm, except when m is within 1e-9 of 0 or 1.  In the
/// latter case with m close to 1, the approximation applies
/// only for phi < pi/2.
///
/// ACCURACY:
///
/// Tested at random points with u between 0 and 10, m between
/// 0 and 1.
///
/// ```text
///            Absolute error (* = relative error):
/// arithmetic   function   # trials      peak         rms
///    IEEE      phi         20000       6.8e-16*    1.3e-16*
///    IEEE      sn          20000       3.4e-15     4.3e-16
///    IEEE      cn          20000       2.9e-15     4.3e-16
///    IEEE      dn          20000       8.9e-15     9.2e-16
/// ```
///
///  Peak error observed in consistency check using addition
/// theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
/// the above relation to the incomplete elliptic integral.
/// Accuracy deteriorates when u is large.
#[allow(clippy::many_single_char_names, clippy::manual_midpoint)] // `phi` can overflow
pub(crate) fn ellpj(u: f64, m: f64) -> (f64, f64, f64, f64) {
    let mut a = [0.0; 9];
    let mut c = [0.0; 9];

    /* Check for special cases */
    if m > 1.0 || m.is_nan() {
        return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    }
    /* Handle negative m by transformation (DLMF 22.17) */
    if m < 0.0 {
        let m1 = -m / (1.0 - m);
        let sqrt_1_minus_m = (1.0 - m).sqrt();
        let u1 = sqrt_1_minus_m * u;

        // Unlike in xsf, the outputs are NaN (instead of unset) if this fails, i.e. if m = -inf
        let (sn_p, cn_p, dn_p, ph_p) = ellpj(u1, m1);

        let sn = sn_p / (sqrt_1_minus_m * dn_p);
        let cn = cn_p / dn_p;
        let dn = 1.0 / dn_p;
        let ph0 = sn_p.atan2(sqrt_1_minus_m * cn_p);
        let ph = ph0 + PI * ((ph_p - ph0) / PI).round();
        return (sn, cn, dn, ph);
    }
    if m < 1.0e-9 {
        let t = u.sin();
        let b = u.cos();
        let ai = 0.25 * m * (u - t * b);
        let sn = t - ai * b;
        let cn = b + ai * t;
        let ph = u - ai;
        let dn = 1.0 - 0.5 * m * t * t;
        return (sn, cn, dn, ph);
    }
    // Unlike in xsf, m = 1 is handled separately, and the first-order expansion in 1 - m below is
    // only used if x = (1 - m) cosh(u) is small. Otherwise the AGM is used.
    //
    // - Its error was observed to grow like x^2, e.g. 3e-8 at x = 1e-3, where the AGM's is ~1e-9.
    // - For |u| > 355, `twon` overflows, so the results were NaN or inf, even for m = 1 (0 * inf).
    //
    // E.g. ellpj(-607.2815430467333, 0.9999999999994756) was (NaN, -inf, inf, -inf).
    #[allow(clippy::float_cmp)]
    if m == 1.0 {
        // sech(u), also for |u| > 710 where cosh(u) overflows
        let b = u.cosh();
        let phi = if b.is_finite() {
            1.0 / b
        } else {
            let t = (-0.5 * u.abs()).exp();
            2.0 * t * t
        };
        return (u.tanh(), phi, phi, 2.0 * u.exp().atan() - FRAC_PI_2);
    }
    if m >= 0.999_999_999_9 && (1.0 - m) * u.cosh() < 1e-3 {
        let mut ai = 0.25 * (1.0 - m);
        let b = u.cosh();
        let t = u.tanh();
        let phi = 1.0 / b;
        let twon = b * u.sinh();
        let sn = t + ai * (twon - u) / (b * b);
        let ph = 2.0 * u.exp().atan() - FRAC_PI_2 + ai * (twon - u) / b;
        ai *= t * phi;
        let cn = phi - ai * (twon - u);
        let dn = phi + ai * (twon + u);
        return (sn, cn, dn, ph);
    }

    /* A. G. M. scale. See DLMF 22.20(ii) */
    a[0] = 1.0;
    let mut b = (1.0 - m).sqrt();
    c[0] = m.sqrt();
    let mut twon = 1.0;
    let mut i = 0;

    while (c[i] / a[i]).abs() > MACHEP {
        if i > 7 {
            break;
        }
        let ai = a[i];
        i += 1;
        c[i] = (ai - b) / 2.0;
        let t = (ai * b).sqrt();
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }

    // done:
    /* backward recurrence */
    let mut phi = twon * a[i] * u;
    loop {
        let t = c[i] * phi.sin() / a[i];
        b = phi;
        phi = (t.asin() + phi) / 2.0;
        i -= 1;
        if i == 0 {
            break;
        }
    }

    let sn = phi.sin();
    let t = phi.cos();
    let cn = t;
    let dnfac = (phi - b).cos();
    /* See discussion after DLMF 22.20.5 */
    let dn = if dnfac.abs() < 0.1 {
        (1.0 - m * sn * sn).sqrt()
    } else {
        t / dnfac
    };
    let ph = phi;
    (sn, cn, dn, ph)
}
