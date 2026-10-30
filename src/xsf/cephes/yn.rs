//! Translated into pure Rust from `xsf/cephes/yn.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `-n` wraps for `n = INT_MIN` (undefined behavior in C++).

use crate::xsf::cephes::{y0, y1};

/*                                                     yn.c
 *
 *     Bessel function of second kind of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, yn();
 * int n;
 *
 * y = yn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The function is evaluated by forward recurrence on
 * n, starting with values computed by the routines
 * y0() and y1().
 *
 * If n = 0 or 1 the routine for y0 or y1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *
 *                      Absolute error, except relative
 *                      when y > 1:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       4.1e-15     6.1e-16
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * yn singularity   x = 0              INFINITY
 * yn overflow                         INFINITY
 *
 * Spot checked against tables for x, n between 0 and 100.
 *
 */
pub(crate) fn yn(mut n: i32, x: f64) -> f64 {
    let sign;

    if n < 0 {
        n = n.wrapping_neg();
        if (n & 1) == 0 {
            /* -1**n */
            sign = 1.0;
        } else {
            sign = -1.0;
        }
    } else {
        sign = 1.0;
    }

    if n == 0 {
        return sign * y0(x);
    }
    if n == 1 {
        return sign * y1(x);
    }

    /* test for overflow */
    if x == 0.0 {
        return f64::NEG_INFINITY * sign;
    } else if x < 0.0 {
        return f64::NAN;
    }

    /* forward recurrence on n */

    let mut anm2 = y0(x);
    let mut anm1 = y1(x);
    let mut k = 1;
    let mut r = f64::from(2 * k);
    let mut an;
    loop {
        an = r * anm1 / x - anm2;
        anm2 = anm1;
        anm1 = an;
        r += 2.0;
        k += 1;
        if !(k < n && an.is_finite()) {
            break;
        }
    }

    sign * an
}
