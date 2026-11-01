//! Translated into pure Rust from `xsf/cephes/ellik.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Copyright 2014, Eric W. Moore
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::utils::MACHEP;
use crate::xsf::cephes::ellpk::ellpk;
use core::f64::consts::{FRAC_PI_2, PI};

/// To calculate legendre's incomplete elliptical integral of the first kind for
/// negative m, we use a power series in phi for small m*phi*phi, an asymptotic
/// series in m for large m*phi*phi* and the relation to Carlson's symmetric
/// integral of the first kind.
///
/// ```text
/// F(phi, m) = sin(phi) * R_F(cos(phi)^2, 1 - m * sin(phi)^2, 1.0)
///           = R_F(c-1, c-m, c)
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
#[allow(clippy::float_cmp, clippy::many_single_char_names)]
fn ellik_neg_m(phi: f64, m: f64) -> f64 {
    let (x, y, z, scale);
    let mut n = 0;
    let mpp = (m * phi) * phi;

    if -mpp < 1e-6 && phi < -m {
        return phi + (-mpp * phi * phi / 30.0 + 3.0 * mpp * mpp / 40.0 + mpp / 6.0) * phi;
    }

    if -mpp > 4e7 {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = (4.0 * sp * sm / (1.0 + cp)).ln();
        let b = -(1.0 + cp / sp / sp - a) / 4.0 / m;
        return (a + b) / sm;
    }

    if phi > 1e-153 && m > -1e305 {
        let s = phi.sin();
        let csc2 = 1.0 / (s * s);
        scale = 1.0;
        x = 1.0 / (phi.tan() * phi.tan());
        y = csc2 - m;
        z = csc2;
    } else {
        scale = phi;
        x = 1.0;
        y = 1.0 - m * scale * scale;
        z = 1.0;
    }

    if x == y && x == z {
        return scale / x.sqrt();
    }

    let a0 = (x + y + z) / 3.0;
    let mut a = a0;
    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    /* Carlson gives 1/pow(3*r, 1.0/6.0) for this constant. if r == eps,
     * it is ~338.38. */
    let mut q = 400.0 * (a0 - x).abs().max((a0 - y).abs().max((a0 - z).abs()));

    while q > a.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        a = (x1 + y1 + z1) / 3.0;
        n += 1;
        q /= 4.0;
    }
    let x = (a0 - x) / a / f64::from(1 << (2 * n));
    let y = (a0 - y) / a / f64::from(1 << (2 * n));
    let z = -(x + y);

    let e2 = x * y - z * z;
    let e3 = x * y * z;

    scale * (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0) / a.sqrt()
}

/// Incomplete elliptic integral of the first kind
///
/// Approximates the integral
///
/// ```text
///                phi
///                 -
///                | |
///                |           dt
/// F(phi | m) =   |    ------------------
///                |                   2
///              | |    sqrt( 1 - m sin t )
///               -
///                0
/// ```
///
/// of amplitude phi and modulus m, using the arithmetic -
/// geometric mean algorithm.
///
/// ACCURACY:
///
/// Tested at random points with m in [0, 1] and phi as indicated.
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE     -10,10        20000      7.0e-16     1.0e-16
/// ```
#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::if_not_else,
    clippy::many_single_char_names
)]
pub(crate) fn ellik(mut phi: f64, m: f64) -> f64 {
    if phi.is_nan() || m.is_nan() {
        return f64::NAN;
    }
    if m > 1.0 {
        return f64::NAN;
    }
    if phi.is_infinite() || m.is_infinite() {
        if m.is_infinite() && phi.is_finite() {
            return 0.0;
        } else if phi.is_infinite() && m.is_finite() {
            return phi;
        }
        return f64::NAN;
    }
    if m == 0.0 {
        return phi;
    }
    let mut a = 1.0 - m;
    if a == 0.0 {
        if phi.abs() >= FRAC_PI_2 {
            return f64::INFINITY;
        }
        /* DLMF 19.6.8, and 4.23.42 */
        return phi.tan().asinh();
    }
    let mut npio2 = (phi / FRAC_PI_2).floor();
    if npio2.abs() % 2.0 == 1.0 {
        npio2 += 1.0;
    }
    let mut k = if npio2 != 0.0 {
        phi -= npio2 * FRAC_PI_2;
        // Unlike in xsf (scipy/xsf#359), the residual is reduced again if |phi| >~ 1e16 left it
        // outside [-pi, pi], e.g. `ellipkinc(-1.79e19, -1.24e35)` was NaN instead of -1348.95.
        if phi.abs() > PI {
            return ellik(phi, m) + npio2 * ellpk(a);
        }
        ellpk(a)
    } else {
        0.0
    };
    let sign = if phi < 0.0 {
        phi = -phi;
        -1
    } else {
        0
    };
    let mut temp = 'done: {
        if a > 1.0 {
            break 'done ellik_neg_m(phi, m);
        }
        let mut b = a.sqrt();
        let mut t = phi.tan();
        if t.abs() > 10.0 {
            /* Transform the amplitude */
            let e = 1.0 / (b * t);
            /* ... but avoid multiple recursions.  */
            if e.abs() < 10.0 {
                let e = e.atan();
                if npio2 == 0.0 {
                    k = ellpk(a);
                }
                break 'done k - ellik(e, m);
            }
        }
        a = 1.0;
        let mut c = m.sqrt();
        let mut d = 1;
        let mut r#mod = 0;

        while (c / a).abs() > MACHEP {
            let temp = b / a;
            phi = phi + (t * temp).atan() + f64::from(r#mod) * PI;
            let denom = 1.0 - temp * t * t;
            if denom.abs() > 10.0 * MACHEP {
                t = t * (1.0 + temp) / denom;
                r#mod = ((phi + FRAC_PI_2) / PI) as i32;
            } else {
                t = phi.tan();
                r#mod = ((phi - t.atan()) / PI).floor() as i32;
            }
            c = (a - b) / 2.0;
            let temp = (a * b).sqrt();
            a = f64::midpoint(a, b);
            b = temp;
            d += d;
        }

        (t.atan() + f64::from(r#mod) * PI) / (f64::from(d) * a)
    };

    // done:
    if sign < 0 {
        temp = -temp;
    }
    temp += npio2 * k;
    temp
}
