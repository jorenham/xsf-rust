//! Translated into pure Rust from `xsf/cephes/incbi.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.4:  March,1996
//! Copyright 1984, 1996 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use super::incbet::incbet;
use crate::xsf::cephes::{lgam, ndtri};

// cephes/const.h
const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53
const MAXLOG: f64 = 709.782_712_893_384; // log(DBL_MAX)
const MINLOG: f64 = -745.133_219_101_941_2; // log(2**(-1075))

/// The `goto` labels in `incbi`
#[derive(Clone, Copy)]
enum Label {
    Ihalve,
    Under,
    Newt,
    Done,
}

/// Inverse of incomplete beta integral
///
/// Given y, the function finds x such that
///
/// ```text
///  incbet( a, b, x ) = y .
/// ```
///
/// The routine performs interval halving or Newton iterations to find the
/// root of incbet(a,b,x) - y = 0.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
///                x     a,b
/// arithmetic   domain  domain  # trials    peak       rms
///    IEEE      0,1    .5,10000   50000    1.9e-11   2.3e-13
///    IEEE      0,1   .25,100    100000    1.4e-13   4.1e-15
///    IEEE      0,1     0,5       50000    1.8e-12   9.0e-15
///    VAX       0,1    .5,100     25000    3.5e-14   1.1e-15
/// With a and b constrained to half-integer or integer values:
///    IEEE      0,1    .5,10000   50000    4.9e-11   3.0e-13
///    IEEE      0,1    .5,100    100000    2.4e-14   1.4e-15
/// With a = .5, b constrained to half-integer or integer values:
///    IEEE      0,1    .5,10000   10000    1.1e-10   1.5e-11
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
#[allow(clippy::float_cmp, clippy::too_many_lines)]
fn incbi(aa: f64, bb: f64, yy0: f64) -> f64 {
    let mut a;
    let mut b;
    let mut y0;
    let mut d;
    let mut y = 0.0;
    let mut x;
    let mut lgm;
    let mut yp;
    let mut di;
    let mut dithresh;
    let mut xt;
    let mut rflg;
    let mut dir;
    let mut nflg;

    if yy0 <= 0.0 {
        return 0.0;
    }
    if yy0 >= 1.0 {
        return 1.0;
    }
    let mut x0 = 0.0;
    let mut yl = 0.0;
    let mut x1 = 1.0;
    let mut yh = 1.0;
    nflg = false;

    let mut label = 'start: {
        if aa <= 1.0 || bb <= 1.0 {
            dithresh = 1.0e-6;
            rflg = false;
            a = aa;
            b = bb;
            y0 = yy0;
            x = a / (a + b);
            y = incbet(a, b, x);
            break 'start Label::Ihalve;
        }
        dithresh = 1.0e-4;

        /* approximation to inverse function */

        yp = -ndtri(yy0);

        if yy0 > 0.5 {
            rflg = true;
            a = bb;
            b = aa;
            y0 = 1.0 - yy0;
            yp = -yp;
        } else {
            rflg = false;
            a = aa;
            b = bb;
            y0 = yy0;
        }

        lgm = (yp * yp - 3.0) / 6.0;
        x = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
        d = yp * (x + lgm).sqrt() / x
            - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * x));
        d *= 2.0;
        if d < MINLOG {
            x = 1.0;
            break 'start Label::Under;
        }
        x = a / (a + b * d.exp());
        y = incbet(a, b, x);
        yp = (y - y0) / y0;
        if yp.abs() < 0.2 {
            break 'start Label::Newt;
        }

        /* Resort to interval halving if not close enough. */
        Label::Ihalve
    };

    loop {
        label = match label {
            Label::Ihalve => 'ihalve: {
                dir = 0;
                di = 0.5;
                for i in 0..100 {
                    if i != 0 {
                        x = x0 + di * (x1 - x0);
                        if x == 1.0 {
                            x = 1.0 - MACHEP;
                        }
                        if x == 0.0 {
                            di = 0.5;
                            x = x0 + di * (x1 - x0);
                            if x == 0.0 {
                                break 'ihalve Label::Under;
                            }
                        }
                        y = incbet(a, b, x);
                        yp = (x1 - x0) / (x1 + x0);
                        if yp.abs() < dithresh {
                            break 'ihalve Label::Newt;
                        }
                        yp = (y - y0) / y0;
                        if yp.abs() < dithresh {
                            break 'ihalve Label::Newt;
                        }
                    }
                    if y < y0 {
                        x0 = x;
                        yl = y;
                        if dir < 0 {
                            dir = 0;
                            di = 0.5;
                        } else if dir > 3 {
                            di = 1.0 - (1.0 - di) * (1.0 - di);
                        } else if dir > 1 {
                            di = 0.5 * di + 0.5;
                        } else {
                            di = (y0 - y) / (yh - yl);
                        }
                        dir += 1;
                        if x0 > 0.75 {
                            if rflg {
                                rflg = false;
                                a = aa;
                                b = bb;
                                y0 = yy0;
                            } else {
                                rflg = true;
                                a = bb;
                                b = aa;
                                y0 = 1.0 - yy0;
                            }
                            x = 1.0 - x;
                            y = incbet(a, b, x);
                            x0 = 0.0;
                            yl = 0.0;
                            x1 = 1.0;
                            yh = 1.0;
                            break 'ihalve Label::Ihalve;
                        }
                    } else {
                        x1 = x;
                        if rflg && x1 < MACHEP {
                            x = 0.0;
                            break 'ihalve Label::Done;
                        }
                        yh = y;
                        if dir > 0 {
                            dir = 0;
                            di = 0.5;
                        } else if dir < -3 {
                            di = di * di;
                        } else if dir < -1 {
                            di *= 0.5;
                        } else {
                            di = (y - y0) / (yh - yl);
                        }
                        dir -= 1;
                    }
                }
                if x0 >= 1.0 {
                    x = 1.0 - MACHEP;
                    break 'ihalve Label::Done;
                }
                if x <= 0.0 {
                    break 'ihalve Label::Under;
                }
                Label::Newt
            },
            Label::Under => {
                x = 0.0;
                Label::Done
            },
            Label::Newt => 'newt: {
                if nflg {
                    break 'newt Label::Done;
                }
                nflg = true;
                lgm = lgam(a + b) - lgam(a) - lgam(b);

                for i in 0..8 {
                    /* Compute the function at this point. */
                    if i != 0 {
                        y = incbet(a, b, x);
                    }
                    if y < yl {
                        x = x0;
                        y = yl;
                    } else if y > yh {
                        x = x1;
                        y = yh;
                    } else if y < y0 {
                        x0 = x;
                        yl = y;
                    } else {
                        x1 = x;
                        yh = y;
                    }
                    if x == 1.0 || x == 0.0 {
                        break;
                    }
                    /* Compute the derivative of the function at this point. */
                    d = (a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() + lgm;
                    if d < MINLOG {
                        break 'newt Label::Done;
                    }
                    if d > MAXLOG {
                        break;
                    }
                    d = d.exp();
                    /* Compute the step to the next approximation of x. */
                    d = (y - y0) / d;
                    xt = x - d;
                    if xt <= x0 {
                        y = (x - x0) / (x1 - x0);
                        xt = x0 + 0.5 * y * (x - x0);
                        if xt <= 0.0 {
                            break;
                        }
                    }
                    if xt >= x1 {
                        y = (x1 - x) / (x1 - x0);
                        xt = x1 - 0.5 * y * (x1 - x);
                        if xt >= 1.0 {
                            break;
                        }
                    }
                    x = xt;
                    if (d / x).abs() < 128.0 * MACHEP {
                        break 'newt Label::Done;
                    }
                }
                /* Did not converge.  */
                dithresh = 256.0 * MACHEP;
                Label::Ihalve
            },
            Label::Done => break,
        };
    }

    // done:
    if rflg {
        if x <= MACHEP { 1.0 - MACHEP } else { 1.0 - x }
    } else {
        x
    }
}

/// Inverse of the regularized incomplete Beta function, $y = \I_x(a, b)$
///
/// Computes $x \in \[0, 1\]$ such that [^DLMF]:
///
/// $$
/// y = \I_x(a, b) = {1 \over \B(a, b)} \int_0^x t^{a-1} (1-t)^{b-1} \dd t ,
/// $$
///
/// for $0 \le y \le 1$ and $a,b \ge 0$, with $\B(a, b)$ the Beta function.
///
/// This function is the quantile function (inverse CDF) of the Beta distribution.
///
/// # Notes
/// This function uses the `incbi` Cephes routine [^CEPHES], making it less accurate than
/// [`scipy.special.betaincinv`][scipy], which wraps the Boost `ibeta_inv` routine [^BOOST].
/// Especially for small `y` and large `a - b` the accuracy may be poor. For example, with
/// $10^3 < a < 10^5$ and $1/2 < b < 10$, the relative error is up to 6e-2 for
/// $10^{-16} \le y < 10^{-5}$, and the result is meaningless for $y < 10^{-16}$, where $1 - y$
/// (almost) rounds to $1$.
///
/// # See also
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function $\I_x(a, b)$
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`scipy.special.betaincinv`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincinv.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17>.
/// [^CEPHES]: Cephes Math Library Release 2.4. Translated into C++ by SciPy developers in 2024.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(
    alias = "inc_beta_inv",
    alias = "beta_inc_inv",
    alias = "ibeta_inv",
    alias = "incbi"
)]
#[must_use]
#[inline]
pub fn betaincinv(a: f64, b: f64, y: f64) -> f64 {
    if a.is_nan() || b.is_nan() || y.is_nan() || a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&y) {
        f64::NAN
    } else {
        incbi(a, b, y)
    }
}

#[cfg(test)]
mod tests {
    //! based on scipy.special.tests.test_basic.TestBetaInc

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_betaincinv_a1_b1() {
        // betaincinv(1, 1, x) is x.
        let x = [0.0, 0.25, 1.0];
        assert_eq!(x.map(|xv| crate::betaincinv(1.0, 1.0, xv)), x);
    }

    #[test]
    fn test_betaincinv_nontrivial() {
        for &(a, b, x, p) in &[
            (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
            (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
            // scipy/scipy#3761:
            (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
            // scipy/scipy#4244 (NOTE: this results in a relative error of 1.872e-6)
            // (
            //     0.0002742794749792665,
            //     289206.03125,
            //     1.639984034231756e-56,
            //     0.9688708782196045,
            // ),
            // scipy/scipy#12796 (NOTE: this results in a relative error of 1.787e2)
            // (4.0, 99997.0, 0.0001947841578892121, 0.999995),
        ] {
            let x1 = crate::betaincinv(a, b, p);
            crate::np_assert_allclose!(&[x1], &[x], rtol = 5e-13);
        }
    }

    #[test]
    fn test_betaincinv_tiny_y() {
        // Test with extremely small y values.  This test includes
        // a regression test for an issue in the boost code;
        // see https://github.com/boostorg/math/issues/961
        //
        // The reference values were computed with mpmath. For example,
        //
        //   from mpmath import mp
        //   mp.dps = 1000
        //   a = 14.208308325339239
        //   p = 7.703145458496392e-307
        //   x = mp.findroot(lambda t: mp.betainc(a, a, 0, t,
        //                                        regularized=True) - p,
        //                   x0=8.566e-23)
        //   print(float(x))

        for &(a, b, y, ref_) in &[
            (
                14.208_308_325_339_239,
                14.208_308_325_339_239,
                7.703_145_458_496_392e-307,
                8.566_004_561_846_704e-023,
            ),
            (14.0, 14.5, 1e-280, 2.934_391_500_664_242_4e-21),
            (3.5, 15.0, 4e-95, 1.329_075_142_928_922_7e-28),
            (10.0, 1.25, 2e-234, 3.982_659_092_143_654e-24),
            // NOTE: This results in a relative error of 3.806e-11
            // (4.0, 99997.0, 5e-88, 3.309800566862242e-27),
        ] {
            let x = crate::betaincinv(a, b, y);
            crate::np_assert_allclose!(&[x], &[ref_], rtol = 1e-14);
        }
    }

    #[test]
    fn test_betaincinv_gh21426() {
        // Test for scipy/scipy#21426: betaincinv must not return NaN
        let a = 5.0;
        let x = 0.5;
        let result = crate::betaincinv(a, a, x);
        crate::np_assert_allclose!(&[result], &[0.5], rtol = 10.0 * f64::EPSILON);
    }
}
