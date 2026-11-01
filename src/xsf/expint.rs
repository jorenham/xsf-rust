//! Translated into pure Rust from `xsf/expint.h` (xsf v0.2.2).
//!
//! The functions exp1, expi below are based on translations of the Fortran code
//! by Shanjie Zhang and Jianming Jin from the book
//!
//!  Shanjie Zhang, Jianming Jin,
//!  Computation of Special Functions,
//!  Wiley, 1996,
//!  ISBN: 0-471-11963-6,
//!  LC: QA351.C45.

use core::f64::consts::PI;

use num_complex::{Complex, c64};

use crate::utils::SCIPY_EULER;
use crate::xsf::complex::{cdiv, cmul, rsub};
use crate::xsf::exp::cexp;
use crate::xsf::log::clog;

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn exp1_f64(x: f64) -> f64 {
    // ============================================
    // Purpose: Compute exponential integral E1(x)
    // Input :  x  --- Argument of E1(x)
    // Output:  E1 --- E1(x)  ( x > 0 )
    // ============================================
    const GA: f64 = SCIPY_EULER;

    if x == 0.0 {
        return f64::INFINITY;
    }
    if x <= 1.0 {
        let mut e1 = 1.0;
        let mut r = 1.0;
        for k in 1..26 {
            let k = f64::from(k);
            r = -r * k * x / (k + 1.0).powi(2);
            e1 += r;
            if r.abs() <= e1.abs() * 1e-15 {
                break;
            }
        }
        return -GA - x.ln() + x * e1;
    }
    let m = 20 + (80.0 / x) as i32;
    let mut t0 = 0.0;
    for k in (1..=m).rev() {
        let k = f64::from(k);
        t0 = k / (1.0 + k / (x + t0));
    }
    let t = 1.0 / (x + t0);
    (-x).exp() * t
}

#[allow(clippy::float_cmp)]
fn exp1_c64(z: Complex<f64>) -> Complex<f64> {
    // ====================================================
    // Purpose: Compute complex exponential integral E1(z)
    // Input :  z   --- Argument of E1(z)
    // Output:  CE1 --- E1(z)
    // ====================================================
    const EL: f64 = SCIPY_EULER;

    let x = z.re;
    let a0 = z.norm();
    // Continued fraction converges slowly near negative real axis,
    // so use power series in a wedge around it until radius 40.0
    let xt = -2.0 * z.im.abs();

    if a0 == 0.0 {
        return c64(f64::INFINITY, 0.0);
    }
    let mut ce1;
    if (a0 < 5.0) || ((x < xt) && (a0 < 40.0)) {
        // Power series
        ce1 = c64(1.0, 0.0);
        let mut cr = c64(1.0, 0.0);
        for k in 1..501 {
            cr = cmul(-cr, z) * (f64::from(k) / f64::from(k + 1).powi(2));
            ce1 += cr;
            if cr.norm() < ce1.norm() * 1e-15 {
                break;
            }
        }
        if (x <= 0.0) && (z.im == 0.0) {
            // Careful on the branch cut -- use the sign of the imaginary part
            //  to get the right sign on the factor if pi.
            ce1 = rsub(-EL, clog(-z)) + cmul(z, ce1) - PI.copysign(z.im) * c64(0.0, 1.0);
        } else {
            ce1 = rsub(-EL, clog(z)) + cmul(z, ce1);
        }
    } else {
        // Continued fraction https://dlmf.nist.gov/6.9
        //                  1     1     1     2     2     3     3
        // E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
        //                Z +   1 +   Z +   1 +   Z +   1 +   Z +
        let mut zc = c64(0.0, 0.0);
        let mut zd = cdiv(c64(1.0, 0.0), z);
        let mut zdc = zd;
        zc += zdc;
        for k in 1..501 {
            let fk = f64::from(k);
            zd = cdiv(c64(1.0, 0.0), zd * fk + 1.0);
            zdc = cmul(zdc, zd - 1.0);
            zc += zdc;

            zd = cdiv(c64(1.0, 0.0), zd * fk + z);
            zdc = cmul(zdc, cmul(z, zd) - 1.0);
            zc += zdc;
            if (zdc.norm() <= zc.norm() * 1e-15) && (k > 20) {
                break;
            }
        }
        ce1 = cmul(cexp(-z), zc);
        if (x <= 0.0) && (z.im == 0.0) {
            ce1 -= PI * c64(0.0, 1.0);
        }
    }
    ce1
}

#[allow(clippy::float_cmp)]
fn expi_f64(x: f64) -> f64 {
    // ============================================
    // Purpose: Compute exponential integral Ei(x)
    // Input :  x  --- Argument of Ei(x)
    // Output:  EI --- Ei(x)
    // ============================================
    const GA: f64 = SCIPY_EULER;

    if x == 0.0 {
        f64::NEG_INFINITY
    } else if x < 0.0 {
        -exp1_f64(-x)
    } else if x.abs() <= 40.0 {
        // Power series around x=0
        let mut ei = 1.0;
        let mut r = 1.0;

        for k in 1..=100 {
            let k = f64::from(k);
            r = r * k * x / ((k + 1.0) * (k + 1.0));
            ei += r;
            if (r / ei).abs() <= 1.0e-15 {
                break;
            }
        }
        GA + x.ln() + x * ei
    } else if x.is_infinite() {
        // Special use-case needed because exp(inf) / inf is undefined/NaN
        f64::INFINITY
    } else {
        // Asymptotic expansion (the series is not convergent)
        let mut ei = 1.0;
        let mut r = 1.0;
        for k in 1..=20 {
            r = r * f64::from(k) / x;
            ei += r;
        }
        x.exp() / x * ei
    }
}

fn expi_c64(z: Complex<f64>) -> Complex<f64> {
    // ============================================
    // Purpose: Compute exponential integral Ei(x)
    // Input :  x  --- Complex argument of Ei(x)
    // Output:  EI --- Ei(x)
    // ============================================
    let mut cei = -exp1_c64(-z);
    if z.im > 0.0 {
        cei += c64(0.0, PI);
    } else if z.im < 0.0 {
        cei -= c64(0.0, PI);
    } else if z.re > 0.0 {
        cei += c64(0.0, PI.copysign(z.im));
    }
    cei
}

//
// Compute a factor of the exponential integral E1.
// This is used in scaled_exp1(x) for moderate values of x.
//
// The function uses the continued fraction expansion given in equation 5.1.22
// of Abramowitz & Stegun, "Handbook of Mathematical Functions".
// For n=1, this is
//
//    E1(x) = exp(-x)*C(x)
//
// where C(x), expressed in the notation used in A&S, is the continued fraction
//
//            1    1    1    2    2    3    3
//    C(x) = ---  ---  ---  ---  ---  ---  ---  ...
//           x +  1 +  x +  1 +  x +  1 +  x +
//
// Here, we pull a factor of 1/z out of C(x), so
//
//    E1(x) = (exp(-x)/x)*F(x)
//
// and a bit of algebra gives the continued fraction expansion of F(x) to be
//
//            1    1    1    2    2    3    3
//    F(x) = ---  ---  ---  ---  ---  ---  ---  ...
//           1 +  x +  1 +  x +  1 +  x +  1 +
//
#[allow(clippy::cast_possible_truncation)]
fn expint1_factor_cont_frac(x: f64) -> f64 {
    // The number of terms to use in the truncated continued fraction
    // depends on x.  Larger values of x require fewer terms.
    let m = 20 + (80.0 / x) as i32;
    let mut t0 = 0.0;
    for k in (1..=m).rev() {
        let k = f64::from(k);
        t0 = k / (x + k / (1.0 + t0));
    }
    1.0 / (1.0 + t0)
}

pub trait ExpIntArg: crate::sealed::Sealed {
    fn expi(self) -> Self;
    fn exp1(self) -> Self;
}

impl ExpIntArg for f64 {
    #[inline]
    fn expi(self) -> Self {
        expi_f64(self)
    }

    #[inline]
    fn exp1(self) -> Self {
        exp1_f64(self)
    }
}

impl ExpIntArg for Complex<f64> {
    #[inline]
    fn expi(self) -> Self {
        expi_c64(self)
    }

    #[inline]
    fn exp1(self) -> Self {
        exp1_c64(self)
    }
}

/// Exponential integral $E_1$ for real or complex input
///
/// Corresponds to [`scipy.special.exp1`][exp1]
///
/// [exp1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp1.html
///
/// # Definition
/// For complex $z \neq 0$ the exponential integral can be defined as
///
/// $$ E_1(z) = \int_z^\infty \frac{e^{-t}}{t} \dd t $$
///
/// where the path of the integral does not cross the negative real axis or pass through the origin.
///
/// # See also
/// - [`scaled_exp1`]: scaled version of $E_1$
/// - [`expn`](crate::expn): generalization of $E_1$
/// - [`expi`]: exponential integral $Ei$
#[must_use]
#[inline]
pub fn exp1<T: ExpIntArg>(z: T) -> T {
    z.exp1()
}

/// Exponential integral $Ei$ for real or complex input
///
/// # Definition
/// For real $x$, the exponential integral is defined as
///
/// $$ Ei(x) = \int_{-\infty}^x \frac{e^t}{t} \dd t $$
///
/// For x > 0 the integral is understood as a Cauchy principal value.
///
/// It is extended to the complex plane by analytic continuation of the function on the interval
/// $(0, \infty)$. The complex variant has a branch cut on the negative real axis.
///
/// # See also
/// - [`exp1`]: exponential integral $E_1$
/// - [`expn`](crate::expn): generalization of $E_1$
#[must_use]
#[inline]
pub fn expi<T: ExpIntArg>(z: T) -> T {
    z.expi()
}

/// Scaled version of the exponential integral $E_1$ for real input
///
/// There exists no corresponding function in `scipy.special`.
///
/// # Definition
/// This function computes $F(x)$, where $F$ is the factor remaining in $E_1(x)$
/// when $e^{-x}/x$ is factored out. That is:
///
/// $$ E_1(x) = x^{-1} e^{-x} F(x) $$
///
/// or
///
/// $$ F(x) = x e^{x} E_1(x) $$
///
/// The function is defined for real $x \geq 0$. For $x < 0$, NaN is returned.
///
/// $F$ has the properties:
///
/// - $F(0) = 0$
/// - $F(x)$ is increasing on $[0, \infty)$.
/// - $F(x) = 1$ in the limit as $x \to \infty$.
///
/// # See also
/// - [`exp1`]: exponential integral $E_1$
/// - [`expi`]: exponential integral $Ei$
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
pub fn scaled_exp1(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return 0.0;
    }

    if x <= 1.0 {
        // For small x, the naive implementation is sufficiently accurate.
        return x * x.exp() * exp1_f64(x);
    }

    if x <= 1250.0 {
        // For moderate x, use the continued fraction expansion.
        return expint1_factor_cont_frac(x);
    }

    // For large x, use the asymptotic expansion.  This is equation 5.1.51
    // from Abramowitz & Stegun, "Handbook of Mathematical Functions".
    1.0 + (-1.0 + (2.0 + (-6.0 + (24.0 - 120.0 / x) / x) / x) / x) / x
}
