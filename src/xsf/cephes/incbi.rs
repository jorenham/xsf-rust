//! Translated into pure Rust from `xsf/cephes/incbi.h` (scipy/xsf@09c93c7).
//!
//! Cephes Math Library Release 2.4:  March,1996
//! Copyright 1984, 1996 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use super::incbet::{INCBET_TINY, betainc_xy, incbet};
use crate::expit;
use crate::utils::{MACHEP, MAXLOG, MINLOG};
use crate::xsf::cephes::{lbeta, lgam, ndtri};

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
///    IEEE      0,1    .5,10000     350    3.2e-15   3.9e-16
///    IEEE      0,1   .25,100       350    6.8e-14   4.8e-15
///    IEEE      0,1     0,5       50000    1.8e-12   9.0e-15
///    VAX       0,1    .5,100     25000    3.5e-14   1.1e-15
/// With a and b constrained to half-integer or integer values:
///    IEEE      0,1    .5,10000     350    3.0e-15   2.6e-16
///    IEEE      0,1    .5,100       350    2.9e-15   4.4e-16
/// With a = .5, b constrained to half-integer or integer values:
///    IEEE      0,1    .5,10000     350    3.7e-15   1.2e-15
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
///
/// Unlike in xsf, the complement of the root is returned as well, `(x, 1 - x)`, which is not
/// rounded if `incbi` ends in the reflected frame (`rflg`), where it is the root.
pub(super) fn incbi(aa: f64, bb: f64, yy0: f64) -> (f64, f64) {
    let mut d;
    let mut y = 0.0;
    let mut x;
    let mut lgm;
    let mut yp;
    let mut di;
    let mut dithresh;
    let mut xt;
    let mut dir;
    let mut nflg;

    if yy0 <= 0.0 {
        return (0.0, 1.0);
    }
    if yy0 >= 1.0 {
        return (1.0, 0.0);
    }
    // `incbet` is NaN for `a = 0` or `b = 0`
    if aa == 0.0 || bb == 0.0 {
        return (f64::NAN, f64::NAN);
    }
    // Unlike in xsf, for `a, b < eps / 1000`, where `I_x(a, b)` is nearly flat, the root of its
    // first order `c + a c ln(x / (1 - x))` with `c = b / (a + b)` is used.
    if aa.min(bb) > 0.0 && aa.max(bb) < INCBET_TINY {
        let c = bb / (aa + bb);
        let t = (yy0 - c) / aa / c;
        return (expit(t), expit(-t));
    }
    let mut x0 = 0.0;
    let mut yl = 0.0;
    let mut x1 = 1.0;
    let mut yh = 1.0;
    nflg = false;

    // `I_x(aa, bb) = yy0`, or if `rflg`, `I_(1-x)(bb, aa) = 1 - yy0`
    let frame = |rflg| {
        if rflg {
            (bb, aa, 1.0 - yy0)
        } else {
            (aa, bb, yy0)
        }
    };
    // Unlike in xsf, `a <= 1` or `b <= 1` with `y0 > 1/2` also starts with `rflg`, since
    // `I_x(a, b)` near 1 has lost the low digits of its complement.
    let mut rflg = yy0 > 0.5;
    let (mut a, mut b, mut y0) = frame(rflg);

    let mut label = 'start: {
        if aa <= 1.0 || bb <= 1.0 {
            dithresh = 1.0e-6;
            x = a / (a + b);
            y = incbet(a, b, x);
            break 'start Label::Ihalve;
        }
        dithresh = 1.0e-4;

        /* approximation to inverse function */

        yp = if rflg { ndtri(yy0) } else { -ndtri(yy0) };

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
                        // Unlike in xsf (scipy/xsf#474), the interval width and the error of y
                        // are relative to the distance to 0 or 1, whichever is closer, since
                        // otherwise, near 1 they can stop the halving too far from the root for
                        // Newton to converge.
                        yp = (x1 - x0) / (x1 + x0).min(2.0 - x1 - x0);
                        if yp.abs() < dithresh {
                            break 'ihalve Label::Newt;
                        }
                        yp = (y - y0) / y0.min(1.0 - y0);
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
                        // Unlike in xsf (scipy/xsf#364), the frame switch requires y0 >= 1/512,
                        // or 1 - x0 <= y0 (rounding 1 - y0 then loses less than staying near
                        // x = 1, i.e. MACHEP / (1 - x0)), or x0 = 1 - MACHEP, the largest x tried.
                        if x0 > 0.75 && (y0 >= 1.0 / 512.0 || 1.0 - x0 <= y0.max(MACHEP)) {
                            rflg = !rflg;
                            (a, b, y0) = frame(rflg);
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
                // Unlike in xsf, the halving starts from `y` at `x`, instead of from a stale `y`
                // or the bisection fraction that the loop above stores in it.
                y = incbet(a, b, x);
                Label::Ihalve
            },
            Label::Done => break,
        };
    }

    // done:
    if rflg {
        (if x <= MACHEP { 1.0 - MACHEP } else { 1.0 - x }, x)
    } else {
        (x, 1.0 - x)
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
/// For example, with $10^3 < a < 10^5$ and $1/2 < b < 10$, the relative error is up to 5e-14 for
/// $10^{-300} \le y < 10^{-5}$, but for $a = 10^{10}$ it can be as large as 6e-5.
///
/// # See also
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function $\I_x(a, b)$
/// - [`betainccinv`]: Inverse of the complement of the regularized incomplete Beta function
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`scipy.special.betaincinv`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincinv.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17>.
/// [^CEPHES]: Cephes Math Library Release 2.4.
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
        incbi(a, b, y).0
    }
}

/// Inverse of the complement of the regularized incomplete Beta function, $y = 1 - \I_x(a, b)$
///
/// Computes $x \in \[0, 1\]$ such that [^DLMF]:
///
/// $$
/// y = 1 - \I_x(a, b) = \I_{1-x}(b, a) ,
/// $$
///
/// for $0 \le y \le 1$ and $a,b > 0$, with $\I_x(a, b)$ the regularized incomplete Beta function.
///
/// This function is the inverse survival function of the Beta distribution.
///
/// # Notes
/// This function uses the `incbi` Cephes routine [^CEPHES] for the root of $\I_{1-x}(b, a) = y$,
/// refined by safeguarded Newton steps, whereas [`scipy.special.betainccinv`][scipy] wraps the
/// Boost `ibetac_inv` routine [^BOOST]. For example, with $10^{-3} \le a, b \le 10^5$, the
/// relative error is up to 4e-13, for $a = 10^{-3}$, where the root is ill-conditioned.
/// Where $1 - \I_x(a, b)$ rounds to $y$ over a range of $x$, e.g. for a subnormal $y$, the result
/// can be anywhere in that range.
///
/// # See also
/// - [`betaincc`](crate::betaincc): Complement of the regularized incomplete Beta function
/// - [`betaincinv`]: Inverse of the regularized incomplete Beta function
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function $\I_x(a, b)$
/// - [`scipy.special.betainccinv`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainccinv.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17.E4>.
/// [^CEPHES]: Cephes Math Library Release 2.4.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(alias = "ibetac_inv")]
#[must_use]
pub fn betainccinv(a: f64, b: f64, y: f64) -> f64 {
    if a.is_nan() || b.is_nan() || y.is_nan() || a <= 0.0 || b <= 0.0 || !(0.0..=1.0).contains(&y) {
        return f64::NAN;
    }
    if y == 0.0 || y == 1.0 {
        return 1.0 - y;
    }
    if a.is_infinite() || b.is_infinite() {
        // As in SciPy, NaN unless a = 1 or b = 1, where the roots are `1 - y^(1/b)` and
        // `(1 - y)^(1/a)`
        return if a == 1.0 {
            0.0
        } else if b == 1.0 {
            1.0
        } else {
            f64::NAN
        };
    }
    let x = incbi(b, a, y).1;
    // `incbi` returns 1 if its root `t` rounds to 0, 0 for roots below eps in its reflected frame
    // or if `t` rounds to 1, and NaN if it fails, e.g. if `lgam(a + b)` overflows
    if x == 1.0 {
        x
    } else {
        betainccinv_refine(a, b, y, if x.is_nan() { 0.5 } else { x })
    }
}

/// Enough for the bisection alone to narrow `ln(x)` from `[ln(5e-324), 0]` to `eps`
const BETAINCCINV_MAXITER: usize = 64;

/// Refines the root `x` of `I_(1-x)(b, a) = y` by Newton steps in `ln(x)` on the log of the
/// smaller side, safeguarded by bisection in `ln(x)` (as `rtsafe` in Numerical Recipes)
///
/// The root from `incbi` can be off, as `1 - t` of its root `t` near 1 (by ~eps / x), or if the
/// rounding of `1 - y` in its reflected frame lands on a plateau of `I_x(a, b)`, and it is 0 for
/// roots below eps.
fn betainccinv_refine(a: f64, b: f64, y: f64, mut x: f64) -> f64 {
    let lbeta_ab = lbeta(a, b);
    // `g(x)` is decreasing, and >= 0 at `lo` and <= 0 at `hi`
    let (mut lo, mut hi) = (0.0, 1.0);
    let mut step = f64::INFINITY;
    for _ in 0..BETAINCCINV_MAXITER {
        // `g = ln(I_(1-x)(b, a) / y)`, or `ln((1 - y) / I_x(a, b))` where `1 - y` is exact, and
        // `p` the one of these that depends on `x`
        let (p, g) = if y < 0.5 {
            let p = betainc_xy(b, a, 1.0 - x, x);
            (p, (p / y).ln())
        } else {
            let p = betainc_xy(a, b, x, 1.0 - x);
            (p, ((1.0 - y) / p).ln())
        };
        if g > 0.0 {
            lo = x;
        } else if g < 0.0 {
            hi = x;
        } else if g == 0.0 {
            break;
        } else {
            return f64::NAN;
        }
        // the root is below the smallest subnormal
        if hi == f64::from_bits(1) {
            return 0.0;
        }
        // within rounding of the root
        if hi - lo <= f64::EPSILON * lo {
            break;
        }
        // the Newton step in `ln(x)`, with `x^a (1 - x)^(b-1) / B(a, b)` the `|dp/d(ln x)|`
        let newton = g * p / (a * x.ln() + (b - 1.0) * (-x).ln_1p() - lbeta_ab).exp();
        let mut xn = x * newton.exp();
        if newton.abs() <= f64::EPSILON {
            return xn;
        }
        // bisection if the step leaves the bracket, or does not halve once there is one
        if !(lo < xn && xn < hi) || (lo > 0.0 && hi < 1.0 && newton.abs() > 0.5 * step) {
            // `lo` is at least the smallest subnormal here
            xn = lo.max(f64::from_bits(1)).sqrt() * hi.sqrt();
            if (xn - x).abs() <= f64::EPSILON * x {
                return xn;
            }
        }
        step = (xn / x).ln().abs();
        x = xn;
    }
    x
}
