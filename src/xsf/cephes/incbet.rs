//! Translated into pure Rust from `xsf/cephes/incbet.h` (scipy/xsf@09c93c7).
//!
//! Cephes Math Library, Release 2.3:  March, 1995
//! Copyright 1984, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::utils::{MACHEP, MAXGAM, MAXLOG, MINLOG, SCIPY_EULER, SQRT_TAU};
use crate::xsf::cephes::{beta, lbeta, lgam_stirling_corr, psi};
use crate::xsf::faddeeva::erfcx;
use crate::xsf::log::log1pmx;

use core::f64::consts::{SQRT_2, TAU};

const INCBET_BIG: f64 = 4.503_599_627_370_496e15;
const INCBET_BIGINV: f64 = 2.220_446_049_250_313e-16;

/// If `a` and `b` are both below this, `I_x(a, b)` is `b / (a + b)` up to a relative
/// `745 max(a, b) < eps`
pub(super) const INCBET_TINY: f64 = 1e-3 * MACHEP;

/// If `m = a max(|ln(x)|, 1/b)` is at most this, `1 - incbet_apser` is `I_x(a, b)` up to its
/// neglected terms of ~`2.5 m^2 < eps / 40` (measured)
const INCBET_APSER_MAX: f64 = 1e-9;

/// `t / (a B(a, b))` if `a` or `b` is subnormal, where `B(a, b)` and `1 / a` can overflow
///
/// Unlike in xsf, which returns NaN or 0 then.
///
/// Here `1 / (a B(a, b)) = b / (a + b) Gamma(a + b + 1) / (Gamma(a + 1) Gamma(b + 1))`, and the
/// ratio of Gamma functions is 1 up to `min(a, b) (psi(a + b + 1) + gamma) < eps`. If `b / (a + b)`
/// is subnormal, `t` is divided first, so that the product keeps its digits.
fn incbet_subnormal_div_a_beta(a: f64, b: f64, t: f64) -> Option<f64> {
    if a.min(b) >= f64::MIN_POSITIVE {
        return None;
    }
    let r = b / (a + b);
    Some(if r >= f64::MIN_POSITIVE {
        t * r
    } else {
        t / (a + b) * b
    })
}

/// The sum of the power series `x^a / B(a, b) (1 / a + sum_(n >= 1) (1 - b)_n / n! x^n / (a + n))`
/// of `I_x(a, b)`, without its `1 / a` term
pub(crate) fn incbet_pseries_sum(a: f64, b: f64, x: f64) -> f64 {
    let mut u = (1.0 - b) * x;
    let mut v = u / (a + 1.0);
    let t1 = v;
    let mut t = u;
    let mut n = 2.0;
    let mut s = 0.0;
    let z = MACHEP * (1.0 / a);
    while v.abs() > z {
        u = (n - b) * x / n;
        t *= u;
        v = t / (a + n);
        s += v;
        n += 1.0;
    }
    s + t1
}

/// Power series for incomplete beta integral.
/// Use when b*x is small and x not too close to 1.
fn incbet_pseries(a: f64, b: f64, x: f64) -> f64 {
    // Unlike in xsf, `I_x(a, b) = 1 - O(a)` is computed as `1 - incbet_apser` for tiny `a`
    if a * x.ln().abs().max(1.0 / b) <= INCBET_APSER_MAX {
        return 1.0 - incbet_apser(a, b, x);
    }
    let ai = 1.0 / a;
    let mut s = incbet_pseries_sum(a, b, x);
    if let Some(p) = incbet_subnormal_div_a_beta(a, b, x.powf(a) * (a * s + 1.0)) {
        return p;
    }
    s += ai;

    let u = a * x.ln();
    if (a + b) < MAXGAM && u.abs() < MAXLOG {
        let t = 1.0 / beta(a, b);
        // `triple_product(s, t, pow(x, a))` in xsf, which only disables optimizations for gcc on
        // 32-bit x86, to prevent aggressive optimizations from degrading accuracy
        s = (s * t) * x.powf(a);
    } else {
        let t = -lbeta(a, b) + u + s.ln();
        // Unlike in xsf, also 0 if `t` is NaN, as `lbeta` can be `inf - inf` if `a` and `b` are
        // both huge, where `x^a` underflows
        s = if t >= MINLOG { t.exp() } else { 0.0 };
    }
    s
}

/// `1 - I_x(a, b)` for tiny `a`, `b x <= 1` and `x <= 0.95` (`apser` of TOMS 708 \[1\])
///
/// To first order in `a`, with `gamma` the Euler-Mascheroni constant,
///
/// ```text
/// 1 - I_x(a, b) = -a (ln(x) + psi(b) + gamma + sum_{j >= 1} (1 - b)_j x^j / (j j!)),
/// ```
///
/// where the neglected terms are smaller by a factor of ~`a max(|ln(x)|, 1/b)`.
///
/// \[1\] See `incbet_lnpow`.
fn incbet_apser(a: f64, b: f64, x: f64) -> f64 {
    let bx = b * x;
    let mut t = x - bx;
    // for b eps > 0.02, psi(b) is ln(b) up to eps / 0.04, and ln(b x) avoids the cancellation
    let ln_x_psi_b = if b * MACHEP <= 0.02 {
        x.ln() + psi(b)
    } else {
        bx.ln()
    };
    let c = ln_x_psi_b + SCIPY_EULER + t;
    let tol = 5.0 * MACHEP * c.abs();
    let mut j = 1.0;
    let mut s = 0.0;
    loop {
        j += 1.0;
        t *= x - bx / j;
        let aj = t / j;
        s += aj;
        if aj.abs() <= tol {
            break;
        }
    }
    -a * (c + s)
}

/// The iteration limit of the continued fractions, which bounds the run time
const INCBET_MAXITER: f64 = 1e6;

/// The iteration limit of the continued fractions at `(a, b, x)`
///
/// Unlike in xsf (scipy/xsf#358), where it is always 300. The continued fractions also stop at
/// NaN, which would otherwise persist until this limit.
fn incbet_maxiter(a: f64, b: f64, x: f64) -> f64 {
    (8.0 * ((a + b).sqrt() + 1.0 / (1.0 - x).sqrt())).clamp(300.0, INCBET_MAXITER)
}

/// `corr(a + b) - corr(a) - corr(b)`, the Stirling correction of `-ln(B(a, b))`, with `corr` the
/// correction to Stirling's formula of `ln(Gamma)` (`lgam_stirling_corr`)
pub(super) fn incbet_corr(a: f64, b: f64) -> f64 {
    lgam_stirling_corr(a + b) - lgam_stirling_corr(a) - lgam_stirling_corr(b)
}

/// `ln(x^a y^b / B(a, b))` for `a, b >= 13` and `x <= a / (a + b)`, with `y = 1 - x`
///
/// Unlike in xsf (scipy/xsf#358), where `incbet` computes this as
/// `a ln(x) + b ln(y) - lbeta(a, b)`.
///
/// Instead, Stirling's formula is used for the Gamma functions in `B(a, b)`, as in `brcomp` of
/// TOMS 708 \[1\]. With `c = a + b` and `d = c x - a = b - c y`, the large linear parts `d` and
/// `-d` of `a ln(c x / a) = a ln(1 + d / a)` and `b ln(c y / b) = b ln(1 - d / b)` then cancel
/// analytically:
///
/// ```text
/// ln(x^a y^b / B(a, b)) = a log1pmx(d / a) + b log1pmx(-d / b) + ln(a b / (2 pi c)) / 2
///                         + incbet_corr(a, b).
/// ```
///
/// \[1\] A. R. DiDonato and A. H. Morris, "Algorithm 708: Significant digit computation of the
/// incomplete beta function ratios", ACM Trans. Math. Softw. 18(3), 360-373 (1992).
pub(super) fn incbet_lnpow(a: f64, b: f64, x: f64, y: f64) -> f64 {
    let c = a + b;
    // without the rounding error of c
    let d = x * b - y * a;
    // 1 + d / a = c x / a can be close to 0, so its log is then computed directly
    let t = d / a;
    let phi_a = if t > -0.5 {
        log1pmx(t)
    } else {
        (c * x / a).ln() - t
    };
    a * phi_a + b * log1pmx(-d / b) + 0.5 * (a / c * (b / TAU)).ln() + incbet_corr(a, b)
}

/// `I_x(a, b)` for `a, b > 100` and `|lambda| <= 0.03 min(a, b)`, with `lambda = a (1 - x) - b x`
/// (`basym` of TOMS 708 \[1\])
///
/// This is an expansion in powers of `w0 = 1 / sqrt(min(a, b) (1 + min(a, b) / max(a, b)))`, with
/// leading term `erfc(z0) / 2`, where `z0^2 = -a log1pmx(-lambda / a) - b log1pmx(lambda / b)`.
///
/// \[1\] See `incbet_lnpow`.
fn incbet_basym(a: f64, b: f64, lambda: f64) -> f64 {
    // the maximum (even) number of terms
    const NUM: usize = 20;

    // `x` above the mean, e.g. within a few ulp of it after the swap in `incbet_xy`
    if lambda < 0.0 {
        return 1.0 - incbet_basym(b, a, -lambda);
    }
    let f = -a * log1pmx(-lambda / a) - b * log1pmx(lambda / b);
    // e^-f with the Stirling corrections of 1 / B(a, b)
    let t = (incbet_corr(a, b) - f).exp();
    if t == 0.0 {
        return 0.0;
    }
    let z0 = f.sqrt();
    let z2 = 2.0 * f;
    let (lo, hi) = (a.min(b), a.max(b));
    let h = lo / hi;
    let r0 = 1.0 / (1.0 + h);
    let r1 = (b - a) / hi;
    let w0 = 1.0 / (lo * (1.0 + h)).sqrt();

    let [mut a0, mut b0, mut c, mut d] = [[0.0; NUM + 1]; 4];
    a0[0] = 2.0 / 3.0 * r1;
    c[0] = -0.5 * a0[0];
    d[0] = -c[0];
    // TOMS 708's `j0` and `j1` times its final factor 2 / sqrt(pi)
    let mut j0 = 0.5 * erfcx(z0);
    let mut j1 = 1.0 / SQRT_TAU;
    let mut sum = j0 + d[0] * w0 * j1;
    let (mut s, mut hn, mut w, mut znm1, mut zn) = (1.0, 1.0, w0, SQRT_2 * z0, z2);
    for n in (2..=NUM).step_by(2) {
        hn *= h * h;
        a0[n - 1] = 2.0 * r0 * (1.0 + h * hn) / (n + 2) as f64;
        s += hn;
        a0[n] = 2.0 * r1 * s / (n + 3) as f64;
        for i in n..=n + 1 {
            let r = -0.5 * (i + 1) as f64;
            b0[0] = r * a0[0];
            for m in 2..=i {
                let bsum: f64 = (1..m)
                    .map(|j| (j as f64 * r - (m - j) as f64) * a0[j - 1] * b0[m - j - 1])
                    .sum();
                b0[m - 1] = r * a0[m - 1] + bsum / m as f64;
            }
            c[i - 1] = b0[i - 1] / (i + 1) as f64;
            let dsum: f64 = (1..i).map(|j| d[i - j - 1] * c[j - 1]).sum();
            d[i - 1] = -(dsum + c[i - 1]);
        }
        j0 = znm1 / SQRT_TAU + (n - 1) as f64 * j0;
        j1 = zn / SQRT_TAU + n as f64 * j1;
        znm1 *= z2;
        zn *= z2;
        w *= w0;
        let t0 = d[n - 1] * w * j0;
        w *= w0;
        let t1 = d[n] * w * j1;
        sum += t0 + t1;
        if t0.abs() + t1.abs() <= MACHEP * sum {
            break;
        }
    }
    t * sum
}

/// `x k1 k2 / (k3 k4)`, a term of the continued fractions
///
/// Unlike in xsf, the ratios are taken first if a product overflows (e.g. `k3 k4 ~ a^2` for
/// `a > 1.3e154`).
fn incbet_term(x: f64, k1: f64, k2: f64, k3: f64, k4: f64) -> f64 {
    let (num, den) = (x * k1 * k2, k3 * k4);
    if num.is_finite() && den.is_finite() {
        num / den
    } else {
        x * (k1 / k3) * (k2 / k4)
    }
}

/// Continued fraction expansion #1 for incomplete beta integral
fn incbcf(a: f64, b: f64, x: f64) -> f64 {
    incbet_cf(a, x, a + b, b - 1.0, 1.0, incbet_maxiter(a, b, x))
}

/// Continued fraction expansion #2 for incomplete beta integral
///
/// `xc = 1 - x` is passed as well, since `x` can be rounded (see `incbet_xy`).
fn incbd(a: f64, b: f64, x: f64, xc: f64) -> f64 {
    incbet_cf(a, x / xc, b - 1.0, a + b, -1.0, incbet_maxiter(a, b, x))
}

/// The continued fraction of `incbcf` and `incbd`, with iteration limit `nmax`
///
/// Unlike in xsf, `incbcf` and `incbd` share this loop (xsf duplicates it). They only differ in
/// their argument `z`, and in `k2` and `k6`, which are swapped and step by `step` in opposite
/// directions:
///
/// ```text
///           z        k2       k6       step
/// incbcf    x        a + b    b - 1     1
/// incbd     x / xc   b - 1    a + b    -1
/// ```
fn incbet_cf(a: f64, z: f64, mut k2: f64, mut k6: f64, step: f64, nmax: f64) -> f64 {
    let mut xk;
    let mut pk;
    let mut qk;
    let mut t;

    let mut k1 = a;
    let mut k3 = a;
    let mut k4 = a + 1.0;
    let mut k5 = 1.0;
    let mut k7 = a + 1.0;
    let mut k8 = a + 2.0;

    let mut pkm2 = 0.0;
    let mut qkm2 = 1.0;
    let mut pkm1 = 1.0;
    let mut qkm1 = 1.0;
    let mut ans = 1.0;
    let mut r = 1.0;
    let mut n = 0;
    let thresh = 3.0 * MACHEP;
    loop {
        xk = -incbet_term(z, k1, k2, k3, k4);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = incbet_term(z, k5, k6, k7, k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if qk != 0.0 {
            r = pk / qk;
        }
        if r != 0.0 {
            t = ((ans - r) / r).abs();
            ans = r;
        } else {
            t = 1.0;
        }
        // also stop at NaN (see `incbet_maxiter`)
        if t < thresh || r.is_nan() {
            break;
        }

        k1 += 1.0;
        k2 += step;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 -= step;
        k7 += 2.0;
        k8 += 2.0;

        if (qk.abs() + pk.abs()) > INCBET_BIG {
            pkm2 *= INCBET_BIGINV;
            pkm1 *= INCBET_BIGINV;
            qkm2 *= INCBET_BIGINV;
            qkm1 *= INCBET_BIGINV;
        }
        if (qk.abs() < INCBET_BIGINV) || (pk.abs() < INCBET_BIGINV) {
            pkm2 *= INCBET_BIG;
            pkm1 *= INCBET_BIG;
            qkm2 *= INCBET_BIG;
            qkm1 *= INCBET_BIG;
        }

        n += 1;
        if f64::from(n) >= nmax {
            break;
        }
    }

    ans
}

/// Incomplete beta integral
///
/// Returns incomplete beta integral of the arguments, evaluated
/// from zero to x.  The function is defined as
///
/// ```text
///                  x
///     -            -
///    | (a+b)      | |  a-1     b-1
///  -----------    |   t   (1-t)   dt.
///   -     -     | |
///  | (a) | (b)   -
///                 0
/// ```
///
/// The domain of definition is 0 <= x <= 1.  In this
/// implementation a and b are restricted to positive values.
/// The integral from x to 1 may be obtained by the symmetry
/// relation
///
/// ```text
///    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
/// ```
///
/// The integral is evaluated by a continued fraction expansion
/// or, when b*x is small, by a power series.
///
/// ACCURACY:
///
/// Tested at uniformly distributed random points (a,b,x) with a and b
/// in "domain" and x between 0 and 1.
///
/// ```text
///                                        Relative error
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0,5         10000       3.6e-15     4.9e-16
///    IEEE      0,85       250000       1.5e-13     1.7e-14
///    IEEE      0,1000        400       2.1e-13     3.8e-14
///    IEEE      0,10000       400       1.9e-13     3.8e-14
///    IEEE      0,100000      400       5.5e-13     7.9e-14
/// ```
///
/// Outputs smaller than the IEEE gradual underflow threshold
/// were excluded from these statistics.
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// incbet domain      x<0, x>1          NaN
/// incbet underflow                     0.0
/// ```
pub(super) fn incbet(aa: f64, bb: f64, xx: f64) -> f64 {
    incbet_xy(aa, bb, xx, 1.0 - xx)
}

/// `incbet(aa, bb, xx)`, given `yy = 1 - xx`
///
/// Unlike in xsf (scipy/xsf#358), `yy` is passed, which is exact also if `xx` is rounded (e.g. in
/// `fdtr`).
pub(super) fn incbet_xy(aa: f64, bb: f64, xx: f64, yy: f64) -> f64 {
    let mut flag = false;

    if aa <= 0.0 || bb <= 0.0 {
        // domerr
        return f64::NAN;
    }

    if (xx <= 0.0) || (yy <= 0.0) {
        if xx == 0.0 {
            return 0.0;
        }
        if yy == 0.0 {
            return 1.0;
        }
        // domerr
        return f64::NAN;
    }

    let t = 'done: {
        // Unlike in xsf, `I_x(a, b)` is `b / (a + b)` for tiny `a` and `b` (see `INCBET_TINY`)
        if aa.max(bb) < INCBET_TINY {
            break 'done bb / (aa + bb);
        }
        if (bb * xx) <= 1.0 && xx <= 0.95 {
            break 'done incbet_pseries(aa, bb, xx);
        }
        // Unlike in xsf, `I_x(a, b) = O(b)` is computed as `incbet_apser(b, a, 1 - x)` for tiny `b`,
        // i.e. `b <= eps min(1, a)` as in TOMS 708. Being off by a relative ~`b max(|ln(1 - x)|,
        // 1/a)`, it is also used up to `INCBET_APSER_MAX` where the continued fractions would be cut
        // off at `INCBET_MAXITER`, and are worse.
        if aa * yy <= 1.0
            && (bb <= MACHEP * aa.min(1.0)
                || bb * yy.ln().abs().max(1.0 / aa) <= INCBET_APSER_MAX
                    && incbet_maxiter(aa, bb, xx) >= INCBET_MAXITER)
        {
            break 'done incbet_apser(bb, aa, yy);
        }

        /* Reverse a and b if x is greater than the mean. */
        // Unlike in xsf (scipy/xsf#358), `yy` is also compared with its mean, as `xx` or the mean
        // can be rounded to 1.
        let past_mean = xx > (aa / (aa + bb)) || yy < bb / (aa + bb);
        // Unlike in xsf (scipy/xsf#358), they are also reversed if the power series of the
        // complement applies and its absolute error of ~(1 - ln(b)) eps, relative to its result of
        // at least ~b/5 (both measured), is below the ~0.1 eps / y of the continued fractions.
        let pseries_y = aa * yy <= 1.0 && 50.0 * (1.0 - bb.ln()) * yy < bb;
        flag = past_mean || pseries_y;
        let (a, b, xc, x) = if flag {
            (bb, aa, xx, yy)
        } else {
            (aa, bb, yy, xx)
        };

        if flag && (b * x) <= 1.0 && x <= 0.95 {
            break 'done incbet_pseries(a, b, x);
        }

        // Unlike in xsf (scipy/xsf#358), `incbet_basym` is used near the mean (as in TOMS 708),
        // where it is generally more accurate than the continued fractions.
        let lambda = a * xc - b * x;
        if a.min(b) > 100.0 && lambda <= 0.03 * a.min(b) {
            break 'done incbet_basym(a, b, lambda);
        }

        /* Choose expansion for better convergence. */
        let mut y = x * (a + b - 2.0) - (a - 1.0);
        // Unlike in xsf (scipy/xsf#358), `incbd` is also used if `x` is close to 1, where `incbcf`
        // amplifies the rounding error of `x` by ~`a / |a xc - b|`.
        let w = if y < 0.0 && x <= 0.95 {
            incbcf(a, b, x)
        } else {
            incbd(a, b, x, xc) / xc
        };

        /* Multiply w by the factor
         * a      b   _             _     _
         * x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

        y = a * x.ln();
        let mut t = b * xc.ln();
        if (a + b) < MAXGAM && y.abs() < MAXLOG && t.abs() < MAXLOG {
            t = xc.powf(b);
            t *= x.powf(a);
            if let Some(p) = incbet_subnormal_div_a_beta(a, b, t * w) {
                break 'done p;
            }
            t /= a;
            t *= w;
            t *= 1.0 / beta(a, b);
            break 'done t;
        }
        /* Resort to logarithms.  */
        // `incbet_lnpow` for `a, b >= 13`, where `lgam` uses Stirling's formula
        if a.min(b) >= 13.0 {
            y = incbet_lnpow(a, b, x, xc);
        } else if x < xc {
            // Unlike in xsf (scipy/xsf#358), the log of the larger (possibly rounded) of `x` and
            // `xc` is `ln_1p` of minus the smaller one.
            y += b * (-x).ln_1p() - lbeta(a, b);
        } else {
            y = a * (-xc).ln_1p() + t - lbeta(a, b);
        }
        y += (w / a).ln();
        if y < MINLOG { 0.0 } else { y.exp() }
    };

    // done:
    // Unlike in xsf, `1 - t` is not clamped to `1 - MACHEP` if `t <= MACHEP`, so that `incbi` can
    // tell `x` far past the root of `I_x(a, b) = 1 - MACHEP` from the root.
    if flag { 1.0 - t } else { t }
}

/// Regularized incomplete Beta function $\I_x(a,b)$
///
/// Computes the regularized incomplete Beta function, defined as [^DLMF]:
///
/// $$
/// \I_x(a,b)
/// = {\Gamma(a+b) \over \Gamma(a)\Gamma(b)} \int_0^x t^{a-1} (1-t)^{b-1} \dd t
/// = {\B_x(a,b) \over \B(a,b)} ,
/// $$
///
/// for $0 \le x \le 1$, with $\B$ the Beta function, and $\B_x$ the (non-regularized)
/// incomplete Beta function.
///
/// # Notes
/// This function uses the `incbet` Cephes routine [^CEPHES], making it less accurate than
/// [`scipy.special.betainc`][scipy], which wraps the Boost `ibeta` routine [^BOOST].
/// Especially for small `x` the accuracy may be poor.
///
/// # See also
/// - [`betaincinv`](crate::betaincinv): Inverse of the regularized incomplete Beta function
/// - [`betaincc`]: Complement of the regularized incomplete Beta function
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`scipy.special.betainc`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17>.
/// [^CEPHES]: Cephes Math Library Release 2.3.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(
    alias = "inc_beta",
    alias = "beta_inc",
    alias = "ibeta",
    alias = "incbet"
)]
#[must_use]
pub fn betainc(a: f64, b: f64, x: f64) -> f64 {
    betainc_xy(a, b, x, 1.0 - x)
}

/// `betainc(a, b, x)`, given `y = 1 - x`, which is exact also if `x` is rounded
pub(crate) fn betainc_xy(a: f64, b: f64, x: f64, y: f64) -> f64 {
    if (a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&x))
        || (a.is_nan() || b.is_nan())
        || (a == 0.0 && b == 0.0)
        || (a.is_infinite() && b.is_infinite())
    {
        f64::NAN
    } else if a == 0.0 || b.is_infinite() {
        if x > 0.0 { 1.0 } else { 0.0 }
    } else if a.is_infinite() || b == 0.0 {
        if y > 0.0 { 0.0 } else { 1.0 }
    } else {
        incbet_xy(a, b, x, y)
    }
}

/// Complement of the regularized incomplete Beta function, $1 - \I_x(a,b)$
///
/// Computes [^DLMF]:
///
/// $$
/// 1 - \I_x(a,b)
/// = \I_{1-x}(b,a)
/// = {\Gamma(a+b) \over \Gamma(a)\Gamma(b)} \int_x^1 t^{a-1} (1-t)^{b-1} \dd t ,
/// $$
///
/// for $0 \le x \le 1$, with $\I_x(a,b)$ the regularized incomplete Beta function.
///
/// This function is the survival function of the Beta distribution.
///
/// # Notes
/// This function uses the `incbet` Cephes routine [^CEPHES] for $\I_{1-x}(b,a)$, given both $x$
/// and $1 - x$, whereas [`scipy.special.betaincc`][scipy] wraps the Boost `ibetac` routine
/// [^BOOST].
///
/// # See also
/// - [`betainc`]: Regularized incomplete Beta function $\I_x(a,b)$
/// - [`betainccinv`](crate::betainccinv): Inverse of the complement of the regularized incomplete
///   Beta function
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`scipy.special.betaincc`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17.E4>.
/// [^CEPHES]: Cephes Math Library Release 2.3.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(alias = "ibetac")]
#[must_use]
pub fn betaincc(a: f64, b: f64, x: f64) -> f64 {
    // `betainc_xy` checks `1 - x`, which rounds to 1 for tiny negative `x`
    if x < 0.0 {
        f64::NAN
    } else {
        betainc_xy(b, a, 1.0 - x, x)
    }
}
