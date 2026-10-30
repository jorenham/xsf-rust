//! Translated into pure Rust from `xsf/cephes/yv.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `(int) v` saturates for `v` outside the `int` range (undefined behavior in C++).

use core::f64::consts::PI;

use crate::xsf::cephes::yn;

fn jv(v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::cephes_jv(v, x) }
}

/*
 * Bessel function of noninteger order
 */
#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
pub(crate) fn yv(v: f64, x: f64) -> f64 {
    let n = v as i32;
    if f64::from(n) == v {
        let y = yn(n, x);
        return y;
    } else if v == v.floor() {
        /* Zero in denominator. */
        return f64::NAN;
    }

    let t = PI * v;
    let y = (t.cos() * jv(v, x) - jv(-v, x)) / t.sin();

    if y.is_infinite() {
        if v > 0.0 {
            return f64::NEG_INFINITY;
        } else if v < -1e10 {
            /* Whether it's +inf or -inf is numerically ill-defined. */
            return f64::NAN;
        }
    }

    y
}
