//! Translated into pure Rust from `xsf/cephes/ellie.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987, 1993 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Copyright 2014, Eric W. Moore

use crate::cosm1;
use crate::xsf::cephes::ellpe::ellpe;
use crate::xsf::cephes::ellpk::ellpk;
use core::f64::consts::{FRAC_PI_2, PI};

const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53

/// To calculate legendre's incomplete elliptical integral of the second kind for
/// negative m, we use a power series in phi for small m*phi*phi, an asymptotic
/// series in m for large m*phi*phi* and the relation to Carlson's symmetric
/// integrals, R_F(x,y,z) and R_D(x,y,z).
///
/// ```text
/// E(phi, m) = sin(phi) * R_F(cos(phi)^2, 1 - m * sin(phi)^2, 1.0)
///             - m * sin(phi)^3 * R_D(cos(phi)^2, 1 - m * sin(phi)^2, 1.0) / 3
///
///           = R_F(c-1, c-m, c) - m * R_D(c-1, c-m, c) / 3
/// ```
///
/// where c = csc(phi)^2. We use the second form of this for (approximately)
/// phi > 1/(sqrt(DBL_MAX) ~ 1e-154, where csc(phi)^2 overflows. Elsewhere we
/// use the first form, accounting for the smallness of phi.
///
/// The algorithm used is described in Carlson, B. C. Numerical computation of
/// real or complex elliptic integrals. (1994) <https://arxiv.org/abs/math/9409227>
/// Most variable names reflect Carlson's usage.
///
/// In this routine, we assume m < 0 and  0 > phi > pi/2.
#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names
)]
fn ellie_neg_m(phi: f64, m: f64) -> f64 {
    let (x, y, z, scalef, scaled);
    let mut n = 0;
    let mpp = (m * phi) * phi;

    if -mpp < 1e-6 && phi < -m {
        return phi + (mpp * phi * phi / 30.0 - mpp * mpp / 40.0 - mpp / 6.0) * phi;
    }

    if -mpp > 1e6 {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = -cosm1(phi);
        let b1 = (4.0 * sp * sm / (1.0 + cp)).ln();
        let b = -(0.5 + b1) / 2.0 / m;
        let c = (0.75 + cp / sp / sp - b1) / 16.0 / m / m;
        return (a + b + c) * sm;
    }

    if phi > 1e-153 && m > -1e200 {
        let s = phi.sin();
        let csc2 = 1.0 / s / s;
        scalef = 1.0;
        scaled = m / 3.0;
        x = 1.0 / phi.tan() / phi.tan();
        y = csc2 - m;
        z = csc2;
    } else {
        scalef = phi;
        scaled = mpp * phi / 3.0;
        x = 1.0;
        y = 1.0 - mpp;
        z = 1.0;
    }

    if x == y && x == z {
        return (scalef + scaled / x) / x.sqrt();
    }

    let a0f = (x + y + z) / 3.0;
    let mut af = a0f;
    let a0d = (x + y + 3.0 * z) / 5.0;
    let mut ad = a0d;
    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    let mut seriesd = 0.0;
    let mut seriesn = 1.0;
    /* Carlson gives 1/pow(3*r, 1.0/6.0) for this constant. if r == eps,
     * it is ~338.38. */

    let mut q = 400.0 * (a0f - x).abs().max((a0f - y).abs().max((a0f - z).abs()));

    while q > af.abs() && q > ad.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        seriesd += seriesn / (sz * (z1 + lam));
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        af = (x1 + y1 + z1) / 3.0;
        ad = (ad + lam) / 4.0;
        n += 1;
        q /= 4.0;
        seriesn /= 4.0;
    }

    let xf = (a0f - x) / af / f64::from(1 << (2 * n));
    let yf = (a0f - y) / af / f64::from(1 << (2 * n));
    let zf = -(xf + yf);

    let e2f = xf * yf - zf * zf;
    let e3f = xf * yf * zf;

    let mut ret = scalef
        * (1.0 - e2f / 10.0 + e3f / 14.0 + e2f * e2f / 24.0 - 3.0 * e2f * e3f / 44.0)
        / af.sqrt();

    let xd = (a0d - x) / ad / f64::from(1 << (2 * n));
    let yd = (a0d - y) / ad / f64::from(1 << (2 * n));
    let zd = -(xd + yd) / 3.0;

    let e2d = xd * yd - 6.0 * zd * zd;
    let e3d = (3.0 * xd * yd - 8.0 * zd * zd) * zd;
    let e4d = 3.0 * (xd * yd - zd * zd) * zd * zd;
    let e5d = xd * yd * zd * zd * zd;

    ret -= scaled
        * (1.0 - 3.0 * e2d / 14.0 + e3d / 6.0 + 9.0 * e2d * e2d / 88.0
            - 3.0 * e4d / 22.0
            - 9.0 * e2d * e3d / 52.0
            + 3.0 * e5d / 26.0)
        / f64::from(1 << (2 * n))
        / ad
        / ad.sqrt();
    ret -= 3.0 * scaled * seriesd;
    ret
}

/// Incomplete elliptic integral of the second kind
///
/// Approximates the integral
///
/// ```text
///                 phi
///                  -
///                 | |
///                 |                   2
/// E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
///                 |
///               | |
///                -
///                 0
/// ```
///
/// of amplitude phi and modulus m, using the arithmetic -
/// geometric mean algorithm.
///
/// ACCURACY:
///
/// Tested at random arguments with phi in [-10, 10] and m in
/// [0, 1].
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE     -10,10       20000       1.0e-15     1.2e-16
/// ```
#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names
)]
pub(crate) fn ellie(phi: f64, m: f64) -> f64 {
    if phi.is_nan() || m.is_nan() {
        return f64::NAN;
    }
    if m > 1.0 {
        return f64::NAN;
    }
    if phi.is_infinite() {
        return phi;
    }
    if m.is_infinite() {
        return -m;
    }
    if m == 0.0 {
        return phi;
    }
    let mut lphi = phi;
    let mut npio2 = (lphi / FRAC_PI_2).floor();
    if npio2.abs() % 2.0 == 1.0 {
        npio2 += 1.0;
    }
    lphi -= npio2 * FRAC_PI_2;
    let sign = if lphi < 0.0 {
        lphi = -lphi;
        -1
    } else {
        1
    };
    let mut a = 1.0 - m;
    let e_ = ellpe(m);
    let mut temp = 'done: {
        if a == 0.0 {
            break 'done lphi.sin();
        }
        if a > 1.0 {
            let temp = ellie_neg_m(lphi, m);
            // Unlike in xsf, a NaN from a residual lphi > pi is retried with lphi reduced again.
            // For |phi| >~ 1e16, the rounding error of `npio2 * FRAC_PI_2` exceeds pi/2, but
            // `ellie_neg_m` assumes 0 <= phi <= pi/2, e.g. `ellipeinc(-1.79e19, -1.24e35)` was
            // NaN (the log of a negative number) instead of -4.02e36.
            if temp.is_nan() && lphi > PI {
                break 'done ellie(lphi, m);
            }
            break 'done temp;
        }

        if lphi < 0.135 {
            let m11 = (((((-7.0 / 2816.0) * m + (5.0 / 1056.0)) * m - (7.0 / 2640.0)) * m
                + (17.0 / 41580.0))
                * m
                - (1.0 / 155_925.0))
                * m;
            let m9 = ((((-5.0 / 1152.0) * m + (1.0 / 144.0)) * m - (1.0 / 360.0)) * m
                + (1.0 / 5670.0))
                * m;
            let m7 = ((-m / 112.0 + (1.0 / 84.0)) * m - (1.0 / 315.0)) * m;
            let m5 = (-m / 40.0 + (1.0 / 30.0)) * m;
            let m3 = -m / 6.0;
            let p2 = lphi * lphi;

            break 'done ((((m11 * p2 + m9) * p2 + m7) * p2 + m5) * p2 + m3) * p2 * lphi + lphi;
        }
        let mut t = lphi.tan();
        let mut b = a.sqrt();
        /* Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
         * for pointing out an instability near odd multiples of pi/2.  */
        if t.abs() > 10.0 {
            /* Transform the amplitude */
            let e = 1.0 / (b * t);
            /* ... but avoid multiple recursions.  */
            if e.abs() < 10.0 {
                let e = e.atan();
                break 'done e_ + m * lphi.sin() * e.sin() - ellie(e, m);
            }
        }
        let mut c = m.sqrt();
        a = 1.0;
        let mut d = 1;
        let mut e = 0.0;
        let mut r#mod = 0;

        while (c / a).abs() > MACHEP {
            let temp = b / a;
            lphi = lphi + (t * temp).atan() + f64::from(r#mod) * PI;
            let denom = 1.0 - temp * t * t;
            if denom.abs() > 10.0 * MACHEP {
                t = t * (1.0 + temp) / denom;
                r#mod = ((lphi + FRAC_PI_2) / PI) as i32;
            } else {
                t = lphi.tan();
                r#mod = ((lphi - t.atan()) / PI).floor() as i32;
            }
            c = (a - b) / 2.0;
            let temp = (a * b).sqrt();
            a = f64::midpoint(a, b);
            b = temp;
            d += d;
            e += c * lphi.sin();
        }

        let mut temp = e_ / ellpk(1.0 - m);
        temp *= (t.atan() + f64::from(r#mod) * PI) / (f64::from(d) * a);
        temp += e;
        temp
    };

    // done:

    if sign < 0 {
        temp = -temp;
    }
    temp += npio2 * e_;
    temp
}
