//! Translated into pure Rust from `xsf/cephes/struve.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024. Original header comment appears below.
//!
//! Compute the Struve function.
//!
//! # Notes
//!
//! We use three expansions for the Struve function discussed in \[1\]:
//!
//! - power series
//! - expansion in Bessel functions
//! - asymptotic large-z expansion
//!
//! Rounding errors are estimated based on the largest terms in the sums.
//!
//! ``struve_convergence.py`` plots the convergence regions of the different
//! expansions.
//!
//! (i)
//!
//! Looking at the error in the asymptotic expansion, one finds that
//! it's not worth trying if z ~> 0.7 * v + 12 for v > 0.
//!
//! (ii)
//!
//! The Bessel function expansion tends to fail for |z| >~ |v| and is not tried
//! there.
//!
//! For Struve H it covers the quadrant v > z where the power series may fail to
//! produce reasonable results.
//!
//! (iii)
//!
//! The three expansions together cover for Struve H the region z > 0, v real.
//!
//! They also cover Struve L, except that some loss of precision may occur around
//! the transition region z ~ 0.7 |v|, v < 0, |v| >> 1 where the function changes
//! rapidly.
//!
//! (iv)
//!
//! The power series is evaluated in double-double precision. This fixes accuracy
//! issues in Struve H for |v| << |z| before the asymptotic expansion kicks in.
//! Moreover, it improves the Struve L behavior for negative v.
//!
//! # References
//!
//! \[1\] NIST Digital Library of Mathematical Functions
//!     <https://dlmf.nist.gov/11>
//!
//! ```text
//! Copyright (C) 2013  Pauli Virtanen
//!
//! Redistribution and use in source and binary forms, with or without
//! modification, are permitted provided that the following conditions are met:
//!
//! a. Redistributions of source code must retain the above copyright notice,
//!    this list of conditions and the following disclaimer.
//! b. Redistributions in binary form must reproduce the above copyright
//!    notice, this list of conditions and the following disclaimer in the
//!    documentation and/or other materials provided with the distribution.
//! c. Neither the name of Enthought nor the names of the SciPy Developers
//!    may be used to endorse or promote products derived from this software
//!    without specific prior written permission.
//!
//! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
//! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
//! OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
//! THE POSSIBILITY OF SUCH DAMAGE.
//! ```
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `(int) v` and `(int) (-v - 0.5)` saturate for `v` outside the `int` range
//! (undefined behavior in C++).

use core::f64::consts::{FRAC_2_SQRT_PI, PI};

use crate::utils::SQRT_PI;
use crate::xsf::bessel::{cyl_bessel_j_f64, cyl_bessel_y_f64};
use crate::xsf::cephes::dd_real::DoubleDouble;
use crate::xsf::cephes::{gammasgn, iv, lgam, rgamma};

const STRUVE_MAXITER: i32 = 10000;
const STRUVE_SUM_EPS: f64 = 1e-16; /* be sure we are in the tail of the sum */
const STRUVE_SUM_TINY: f64 = 1e-100;
const STRUVE_GOOD_EPS: f64 = 1e-12;
const STRUVE_ACCEPTABLE_EPS: f64 = 1e-7;
const STRUVE_ACCEPTABLE_ATOL: f64 = 1e-300;

/*
 * Large-z expansion for Struve H and L
 * https://dlmf.nist.gov/11.6.1
 */
#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn struve_asymp_large_z(v: f64, z: f64, is_h: bool, err: &mut f64) -> f64 {
    let sgn = if is_h { -1 } else { 1 };

    /* Asymptotic expansion divergenge point */
    let m = z / 2.0;
    let maxiter = if m <= 0.0 {
        0
    } else if m > f64::from(STRUVE_MAXITER) {
        STRUVE_MAXITER
    } else {
        m as i32
    };
    if maxiter == 0 {
        *err = f64::INFINITY;
        return f64::NAN;
    }

    if z < v {
        /* Exclude regions where our error estimation fails */
        *err = f64::INFINITY;
        return f64::NAN;
    }

    /* Evaluate sum */
    let mut term = f64::from(-sgn) / SQRT_PI
        * (-lgam(v + 0.5) + (v - 1.0) * (z / 2.0).ln()).exp()
        * gammasgn(v + 0.5);
    let mut sum = term;
    let mut maxterm: f64 = 0.0;

    for n in 0..maxiter {
        term *= f64::from(sgn * (1 + 2 * n)) * (f64::from(1 + 2 * n) - 2.0 * v) / (z * z);
        sum += term;
        if term.abs() > maxterm {
            maxterm = term.abs();
        }
        if term.abs() < STRUVE_SUM_EPS * sum.abs() || term == 0.0 || !sum.is_finite() {
            break;
        }
    }

    if is_h {
        sum += cyl_bessel_y_f64(v, z);
    } else {
        // Unlike in xsf (scipy/xsf#359), this uses I_{-v} if it differs from I_v by more than
        // 1e-12, since the expansion also holds for L_v - I_{-v} (DLMF 10.27.4), and only then is
        // exact at half-integers (DLMF 11.4.7), e.g. `struve_l(-15.5000001, 15)` was 162.7429
        // instead of 162.7430.
        let (i_v, i_mv) = (iv(v, z), iv(-v, z));
        sum += if (i_mv - i_v).abs() > STRUVE_GOOD_EPS * i_v.abs() {
            i_mv
        } else {
            i_v
        };
    }

    /*
     * This error estimate is strictly speaking valid only for
     * n > v - 0.5, but numerical results indicate that it works
     * reasonably.
     */
    *err = term.abs() + maxterm.abs() * STRUVE_SUM_EPS;

    sum
}

/*
 * Power series for Struve H and L
 * https://dlmf.nist.gov/11.2.1
 *
 * Starts to converge roughly at |n| > |z|
 */
#[allow(clippy::float_cmp, clippy::manual_range_contains)]
fn struve_power_series(v: f64, z: f64, is_h: bool, err: &mut f64) -> f64 {
    let mut term;
    let mut sum;
    let scaleexp;
    let mut cdiv;
    let mut ctmp;

    let sgn = if is_h { -1 } else { 1 };

    let mut tmp = -lgam(v + 1.5) + (v + 1.0) * (z / 2.0).ln();
    if tmp < -600.0 || tmp > 600.0 {
        /* Scale exponent to postpone underflow/overflow */
        scaleexp = tmp / 2.0;
        tmp -= scaleexp;
    } else {
        scaleexp = 0.0;
    }

    term = FRAC_2_SQRT_PI * tmp.exp() * gammasgn(v + 1.5);
    sum = term;
    let mut maxterm: f64 = 0.0;

    let mut cterm = DoubleDouble::new(term);
    let mut csum = DoubleDouble::new(sum);
    let z2 = DoubleDouble::new(f64::from(sgn) * z * z);
    let c2v = DoubleDouble::new(2.0 * v);

    for n in 0..STRUVE_MAXITER {
        /* cdiv = (3 + 2*n) * (3 + 2*n + 2*v)) */
        cdiv = DoubleDouble::new(f64::from(3 + 2 * n));
        ctmp = DoubleDouble::new(f64::from(3 + 2 * n));
        ctmp = ctmp + c2v;
        cdiv = cdiv * ctmp;

        /* cterm *= z2 / cdiv */
        cterm = cterm * z2;
        cterm = cterm / cdiv;

        csum = csum + cterm;

        term = f64::from(cterm);
        sum = f64::from(csum);

        if term.abs() > maxterm {
            maxterm = term.abs();
        }
        if term.abs() < STRUVE_SUM_TINY * sum.abs() || term == 0.0 || !sum.is_finite() {
            break;
        }
    }

    *err = term.abs() + maxterm.abs() * 1e-22;

    if scaleexp != 0.0 {
        sum *= scaleexp.exp();
        *err *= scaleexp.exp();
    }

    if sum == 0.0 && term == 0.0 && v < 0.0 && !is_h {
        /* Spurious underflow */
        *err = f64::INFINITY;
        return f64::NAN;
    }

    sum
}

/*
 * Bessel series
 * https://dlmf.nist.gov/11.4.19
 */
#[allow(clippy::float_cmp)]
fn struve_bessel_series(v: f64, z: f64, is_h: bool, err: &mut f64) -> f64 {
    let mut term = 0.0;

    if is_h && v < 0.0 {
        /* Works less reliably in this region */
        *err = f64::INFINITY;
        return f64::NAN;
    }

    let mut sum = 0.0;
    let mut maxterm: f64 = 0.0;

    let mut cterm = (z / (2.0 * PI)).sqrt();

    for n in 0..STRUVE_MAXITER {
        if is_h {
            term = cterm * cyl_bessel_j_f64(f64::from(n) + v + 0.5, z) / (f64::from(n) + 0.5);
            cterm *= z / 2.0 / f64::from(n + 1);
        } else {
            term = cterm * iv(f64::from(n) + v + 0.5, z) / (f64::from(n) + 0.5);
            cterm *= -z / 2.0 / f64::from(n + 1);
        }
        sum += term;
        if term.abs() > maxterm {
            maxterm = term.abs();
        }
        if term.abs() < STRUVE_SUM_EPS * sum.abs() || term == 0.0 || !sum.is_finite() {
            break;
        }
    }

    *err = term.abs() + maxterm.abs() * 1e-16;

    /* Account for potential underflow of the Bessel functions */
    *err += 1e-300 * cterm.abs();

    sum
}

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn struve_hl(v: f64, z: f64, is_h: bool) -> f64 {
    let mut value = [0.0; 4];
    let mut err = [0.0; 4];
    let mut tmp;

    if z < 0.0 {
        let n = v as i32;
        if v == f64::from(n) {
            tmp = if n % 2 == 0 { -1.0 } else { 1.0 };
            return tmp * struve_hl(v, -z, is_h);
        }
        return f64::NAN;
    } else if z == 0.0 {
        if v < -1.0 {
            return gammasgn(v + 1.5) * f64::INFINITY;
        } else if v == -1.0 {
            return FRAC_2_SQRT_PI * rgamma(0.5);
        }
        return 0.0;
    }

    let n = (-v - 0.5) as i32;
    if f64::from(n) == -v - 0.5 && n > 0 {
        if is_h {
            return f64::from(if n % 2 == 0 { 1 } else { -1 })
                * cyl_bessel_j_f64(f64::from(n) + 0.5, z);
        }
        return iv(f64::from(n) + 0.5, z);
    }

    /* Try the asymptotic expansion */
    if z >= 0.7 * v + 12.0 {
        value[0] = struve_asymp_large_z(v, z, is_h, &mut err[0]);
        if err[0] < STRUVE_GOOD_EPS * value[0].abs() {
            return value[0];
        }
    } else {
        err[0] = f64::INFINITY;
    }

    /* Try power series */
    value[1] = struve_power_series(v, z, is_h, &mut err[1]);
    if err[1] < STRUVE_GOOD_EPS * value[1].abs() {
        return value[1];
    }

    /* Try bessel series */
    if z.abs() < v.abs() + 20.0 {
        value[2] = struve_bessel_series(v, z, is_h, &mut err[2]);
        if err[2] < STRUVE_GOOD_EPS * value[2].abs() {
            return value[2];
        }
    } else {
        err[2] = f64::INFINITY;
    }

    /* Return the best of the three, if it is acceptable */
    let mut n = 0;
    if err[1] < err[n] {
        n = 1;
    }
    if err[2] < err[n] {
        n = 2;
    }
    if err[n] < STRUVE_ACCEPTABLE_EPS * value[n].abs() || err[n] < STRUVE_ACCEPTABLE_ATOL {
        return value[n];
    }

    /* Maybe it really is an overflow? */
    tmp = -lgam(v + 1.5) + (v + 1.0) * (z / 2.0).ln();
    if !is_h {
        tmp = tmp.abs();
    }
    if tmp > 700.0 {
        return f64::INFINITY * gammasgn(v + 1.5);
    }

    /* Failure */
    f64::NAN
}

pub(crate) fn struve_h(v: f64, z: f64) -> f64 {
    struve_hl(v, z, true)
}

pub(crate) fn struve_l(v: f64, z: f64) -> f64 {
    struve_hl(v, z, false)
}
