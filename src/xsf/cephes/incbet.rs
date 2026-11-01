//! Translated into pure Rust from `xsf/cephes/incbet.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library, Release 2.3:  March, 1995
//! Copyright 1984, 1995 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::cephes::{beta, lbeta, lgam_stirling_corr};
use crate::xsf::log::log1pmx;

use core::f64::consts::TAU;

// cephes/const.h
const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53
const MAXLOG: f64 = 709.782_712_893_384; // log(DBL_MAX)
const MINLOG: f64 = -745.133_219_101_941_2; // log(2**(-1075))
const MAXGAM: f64 = 171.624_376_956_302_7; // Largest x such that Gamma(x) is finite

const INCBET_BIG: f64 = 4.503_599_627_370_496e15;
const INCBET_BIGINV: f64 = 2.220_446_049_250_313e-16;

/// Power series for incomplete beta integral.
/// Use when b*x is small and x not too close to 1.
fn incbet_pseries(a: f64, b: f64, x: f64) -> f64 {
    let ai = 1.0 / a;
    let mut u = (1.0 - b) * x;
    let mut v = u / (a + 1.0);
    let t1 = v;
    let mut t = u;
    let mut n = 2.0;
    let mut s = 0.0;
    let z = MACHEP * ai;
    while v.abs() > z {
        u = (n - b) * x / n;
        t *= u;
        v = t / (a + n);
        s += v;
        n += 1.0;
    }
    s += t1;
    s += ai;

    u = a * x.ln();
    if (a + b) < MAXGAM && u.abs() < MAXLOG {
        t = 1.0 / beta(a, b);
        // `triple_product(s, t, pow(x, a))` in xsf, which only disables optimizations for gcc on
        // 32-bit x86, to prevent aggressive optimizations from degrading accuracy
        s = (s * t) * x.powf(a);
    } else {
        t = -lbeta(a, b) + u + s.ln();
        if t < MINLOG {
            s = 0.0;
        } else {
            s = t.exp();
        }
    }
    s
}

/// Unlike in xsf (scipy/xsf#358), the continued fractions are not limited to 300 iterations, which
/// is too few near the mean for large `a + b`, and for `x` near 1 with small `b`, where they need
/// O(sqrt(a + b)) and O(1/sqrt(1 - x)) iterations, respectively. For example,
/// `incbet(1e8, 1e8, 0.4999964644661029)` needs 2623 iterations and was 0.45307 instead of
/// 0.46017, and `incbet(1e4, 0.002, 0.9999996)` needs 6942 iterations and was 0.0097217 instead
/// of 0.0098510. The limit of 1e6 iterations bounds the run time.
fn incbet_maxiter(a: f64, b: f64, x: f64) -> f64 {
    (8.0 * ((a + b).sqrt() + 1.0 / (1.0 - x).sqrt())).clamp(300.0, 1e6)
}

/// `ln(x^a y^b / B(a, b))` for `a, b >= 13` and `x <= a / (a + b)`, with `y = 1 - x`
///
/// Unlike in xsf (scipy/xsf#358), where `incbet` computes this as
/// `a ln(x) + b ln(y) - lbeta(a, b)`, whose terms are ~`(a + b) ln(a + b)` and cancel, so that e.g.
/// `incbet(1e8, 2e8, 0.3333)` had a relative error of 4e-7 (now 1e-12). Instead, Stirling's formula
/// is used for the Gamma functions in `B(a, b)`, as in `brcomp` of TOMS 708 \[1\]. With
/// `c = a + b` and `d = c x - a = b - c y`, this gives `a ln(c x / a) = a ln(1 + d / a)` and
/// `b ln(c y / b) = b ln(1 - d / b)`, whose large linear parts `d` and `-d` cancel analytically,
/// so that
///
/// ```text
/// ln(x^a y^b / B(a, b)) = a log1pmx(d / a) + b log1pmx(-d / b) + ln(a b / (2 pi c)) / 2
///                         + corr(c) - corr(a) - corr(b),
/// ```
///
/// where `corr` is the correction to Stirling's formula (`lgam_stirling_corr`).
///
/// \[1\] A. R. DiDonato and A. H. Morris, "Algorithm 708: Significant digit computation of the
/// incomplete beta function ratios", ACM Trans. Math. Softw. 18(3), 360-373 (1992).
fn incbet_lnpow(a: f64, b: f64, x: f64, y: f64) -> f64 {
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
    a * phi_a + b * log1pmx(-d / b) + 0.5 * (a / c * (b / TAU)).ln() + lgam_stirling_corr(c)
        - lgam_stirling_corr(a)
        - lgam_stirling_corr(b)
}

/// Continued fraction expansion #1 for incomplete beta integral
#[allow(clippy::if_not_else)]
fn incbcf(a: f64, b: f64, x: f64) -> f64 {
    let mut xk;
    let mut pk;
    let mut qk;
    let mut t;

    let mut k1 = a;
    let mut k2 = a + b;
    let mut k3 = a;
    let mut k4 = a + 1.0;
    let mut k5 = 1.0;
    let mut k6 = b - 1.0;
    let mut k7 = k4;
    let mut k8 = a + 2.0;

    let mut pkm2 = 0.0;
    let mut qkm2 = 1.0;
    let mut pkm1 = 1.0;
    let mut qkm1 = 1.0;
    let mut ans = 1.0;
    let mut r = 1.0;
    let mut n = 0;
    let nmax = incbet_maxiter(a, b, x);
    let thresh = 3.0 * MACHEP;
    loop {
        xk = -(x * k1 * k2) / (k3 * k4);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (x * k5 * k6) / (k7 * k8);
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
        // Unlike in xsf (scipy/xsf#358), also stop once the result is NaN (which it then stays),
        // because the larger iteration limit would otherwise make this very slow for huge `a + b`.
        if t < thresh || r.is_nan() {
            break;
        }

        k1 += 1.0;
        k2 += 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 -= 1.0;
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

/// Continued fraction expansion #2 for incomplete beta integral
///
/// Unlike in xsf (scipy/xsf#358), `xc = 1 - x` is passed, since `x` can be rounded (see `incbet`).
#[allow(clippy::if_not_else)]
fn incbd(a: f64, b: f64, x: f64, xc: f64) -> f64 {
    let mut xk;
    let mut pk;
    let mut qk;
    let mut t;

    let mut k1 = a;
    let mut k2 = b - 1.0;
    let mut k3 = a;
    let mut k4 = a + 1.0;
    let mut k5 = 1.0;
    let mut k6 = a + b;
    let mut k7 = a + 1.0;
    let mut k8 = a + 2.0;

    let mut pkm2 = 0.0;
    let mut qkm2 = 1.0;
    let mut pkm1 = 1.0;
    let mut qkm1 = 1.0;
    let z = x / xc;
    let mut ans = 1.0;
    let mut r = 1.0;
    let mut n = 0;
    let nmax = incbet_maxiter(a, b, x);
    let thresh = 3.0 * MACHEP;
    loop {
        xk = -(z * k1 * k2) / (k3 * k4);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (z * k5 * k6) / (k7 * k8);
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
        // Unlike in xsf (scipy/xsf#358), also stop once the result is NaN (which it then stays),
        // because the larger iteration limit would otherwise make this very slow for huge `a + b`.
        if t < thresh || r.is_nan() {
            break;
        }

        k1 += 1.0;
        k2 -= 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 += 1.0;
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
///    IEEE      0,1000      30000       4.5e-12     6.1e-13
///    IEEE      0,10000    250000       8.4e-11     6.9e-12
///    IEEE      0,100000    10000       5.8e-10     4.4e-11
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
///
/// (Unlike in xsf, where the value returned on a domain error is documented as 0.0.)
pub(super) fn incbet(aa: f64, bb: f64, xx: f64) -> f64 {
    incbet_xy(aa, bb, xx, 1.0 - xx)
}

/// `incbet(aa, bb, xx)`, given `yy = 1 - xx`
///
/// Unlike in xsf (scipy/xsf#358), where only `xx` is passed, so that `1 - xx` is exact also if `xx`
/// is rounded (see `betainc_xy`).
#[allow(clippy::float_cmp)]
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
        if (bb * xx) <= 1.0 && xx <= 0.95 {
            break 'done incbet_pseries(aa, bb, xx);
        }

        let mut w = yy;

        /* Reverse a and b if x is greater than the mean. */
        let (a, b, xc, x) = if xx > (aa / (aa + bb)) {
            flag = true;
            (bb, aa, xx, w)
        } else {
            (aa, bb, w, xx)
        };

        if flag && (b * x) <= 1.0 && x <= 0.95 {
            break 'done incbet_pseries(a, b, x);
        }

        /* Choose expansion for better convergence. */
        let mut y = x * (a + b - 2.0) - (a - 1.0);
        if x > 0.95 {
            // Unlike in xsf (scipy/xsf#358), `incbd` is used with `xc` instead of `1 - x` if `x` is
            // close to 1 (as for the power series above), where `incbcf`, which only uses `x`,
            // amplifies its rounding error by about `a / |a xc - b|`.
            //
            // E.g. `incbet(15, 5e19, 6e-19)` was NaN instead of 0.99908, since `x = 1 - 6e-19`
            // rounds to 1.
            w = incbd(a, b, x, xc) / xc;
        } else if y < 0.0 {
            w = incbcf(a, b, x);
        } else {
            w = incbd(a, b, x, 1.0 - x) / xc;
        }

        /* Multiply w by the factor
         * a      b   _             _     _
         * x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

        y = a * x.ln();
        let mut t = b * xc.ln();
        if (a + b) < MAXGAM && y.abs() < MAXLOG && t.abs() < MAXLOG {
            t = xc.powf(b);
            t *= x.powf(a);
            t /= a;
            t *= w;
            t *= 1.0 / beta(a, b);
            break 'done t;
        }
        /* Resort to logarithms.  */
        // Unlike in xsf (scipy/xsf#358), `incbet_lnpow` is used for `a, b >= 13` (where `lgam` uses
        // Stirling's formula)
        if a.min(b) >= 13.0 {
            y = incbet_lnpow(a, b, x, xc);
        } else if x < xc {
            // Unlike in xsf (scipy/xsf#358), the log of the larger of `x` and `xc` is computed as
            // `ln_1p` of minus the smaller one, since the larger one can be rounded.
            //
            // E.g. `incbet(5, 5e19, 5e-20)` was 1 instead of 0.10882.
            y += b * (-x).ln_1p() - lbeta(a, b);
        } else {
            y = a * (-xc).ln_1p() + t - lbeta(a, b);
        }
        y += (w / a).ln();
        if y < MINLOG { 0.0 } else { y.exp() }
    };

    // done:
    if flag {
        if t <= MACHEP { 1.0 - MACHEP } else { 1.0 - t }
    } else {
        t
    }
}

/// `betainc(a, b, x)`, given `y = 1 - x` computed without cancellation by the caller
///
/// Unlike in xsf (scipy/xsf#358), where only `x` is passed: if `x` is close to 1 and above the
/// mean, where `incbet` swaps to `1 - x`, the rounding error of `x` is amplified (or `1 - x` even
/// rounds to 0), so `1 - betainc(b, a, y)` is used instead.
///
/// Below the mean, the continued fractions of `incbet` need O(1/sqrt(y)) iterations, over which
/// their relative rounding error grows to ~0.1 eps / y, e.g. `fdtr(1e7, 0.001, 0.5)` had a
/// relative error of 4e-8. So the complement is used there as well if `a y <= 1`, where
/// `betainc(b, a, y)` is a power series with an absolute error of ~(1 - ln(b)) eps, and the result
/// is at least ~b/5 (both measured), if its relative error ~5 (1 - ln(b)) eps / b is smaller.
///
/// Otherwise `y` is passed to `incbet` as well, e.g. `stdtr(1e20, -3)` was 0.5 instead of 0.0013499.
///
/// The mean is compared with `y`, since `x` and `a / (a + b)` can both round to 1, e.g.
/// `fdtrc(100, 1e20, 0.5)` was 1 instead of 0.99999.
pub(crate) fn betainc_xy(a: f64, b: f64, x: f64, y: f64) -> f64 {
    if x > 0.95 && a.is_finite() && b.is_finite() {
        let p = if y < b / (a + b) || (a * y <= 1.0 && 50.0 * (1.0 - b.ln()) * y < b) {
            1.0 - betainc(b, a, y)
        } else {
            incbet_xy(a, b, x, y)
        };
        // `incbet` can be NaN for huge `a` and `b` (or for `a` or `b` zero), and slightly
        // exceed 1
        if !p.is_nan() {
            return p.max(0.0);
        }
    }
    betainc(a, b, x)
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
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`scipy.special.betainc`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17>.
/// [^CEPHES]: Cephes Math Library Release 2.3. Translated into C++ by SciPy developers in 2024.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(
    alias = "inc_beta",
    alias = "beta_inc",
    alias = "ibeta",
    alias = "incbet"
)]
#[must_use]
pub fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if (a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&x))
        || (a.is_nan() || b.is_nan() || x.is_nan())
        || (a == 0.0 && b == 0.0)
        || (a.is_infinite() && b.is_infinite())
    {
        f64::NAN
    } else if a == 0.0 || b.is_infinite() {
        if x > 0.0 { 1.0 } else { 0.0 }
    } else if a.is_infinite() || b == 0.0 {
        if x < 1.0 { 0.0 } else { 1.0 }
    } else {
        incbet(a, b, x)
    }
}

#[cfg(test)]
mod tests {
    // based on scipy.special.tests.test_basic.TestBetaInc

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_betainc_a1_b1() {
        // betainc(1, 1, x) is x.
        let x = [0.0, 0.25, 1.0];
        assert_eq!(x.map(|xv| crate::betainc(1.0, 1.0, xv)), x);
    }

    #[test]
    fn test_betainc_nontrivial() {
        for &(a, b, x, p) in &[
            (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
            (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
            // scipy/scipy#3761:
            (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
            // scipy/scipy#4244 (NOTE: this results in a relative error of 5.134e-10)
            // (
            //     0.0002742794749792665,
            //     289206.03125,
            //     1.639984034231756e-56,
            //     0.9688708782196045,
            // ),
            // scipy/scipy#12796:
            (4.0, 99997.0, 0.000_194_784_157_889_212_1, 0.999_995),
        ] {
            let p1 = crate::betainc(a, b, x);
            // NOTE: The original `rtol = 1e-15` is too strict for Cephes
            crate::np_assert_allclose!(&[p1], &[p], rtol = 1e-14);
        }
    }

    #[test]
    fn test_betainc_slow_continued_fraction() {
        // reference values from mpmath

        // large a + b near the mean
        let p = crate::betainc(1e8, 1e8, 0.499_996_464_466_102_9);
        crate::np_assert_allclose!(&[p], &[0.460_172_162_871_263_6], rtol = 5e-12);

        // x near 1 with small b
        for &(a, b, x, p) in &[
            (1e4, 0.002, 0.999_999_6, 0.009_851_010_029_838_827),
            (116.3, 0.0008, 0.999_98, 0.004_385_278_106_520_217),
            (1000.0, 1e-6, 0.999_999, 6.332_019_221_868_309e-6),
        ] {
            crate::np_assert_allclose!(&[crate::betainc(a, b, x)], &[p], rtol = 1e-10);
        }

        // the result should stay within [0, 1] for huge a + b
        assert!((0.0..=1.0).contains(&crate::betainc(1e15, 1e15, 0.5)));
        assert!((0.0..=1.0).contains(&crate::betainc(2.0, 1e17, 2e-17)));
        assert!((0.0..=1.0).contains(&crate::fdtr(2e15, 2e15, 1.0)));
        assert!((0.0..=1.0).contains(&crate::fdtrc(2e15, 2e15, 1.0)));
    }

    #[test]
    fn test_betainc_large_ab() {
        // reference values from mpmath
        for &(a, b, x, p) in &[
            (1e8, 2e8, 0.3333, 0.110_334_398_548_217_36),
            (5e7, 8e7, 0.3846, 0.359_220_529_919_169_7),
            (3e8, 1e9, 0.2307, 1.562_168_507_020_367_6e-9),
            // tiny x with huge b, where 1 - x is rounded
            (13.0, 1e15, 1e-14, 0.208_443_523_605_136_1),
        ] {
            crate::np_assert_allclose!(&[crate::betainc(a, b, x)], &[p], rtol = 5e-12);
        }
    }

    #[test]
    fn test_betainc_huge_b() {
        // reference values from mpmath
        for &(a, b, x, p) in &[
            (2.0, 1e17, 2e-17, 0.593_994_150_290_162),
            (5.0, 5e19, 5e-20, 0.108_821_981_085_848_76),
            (15.0, 5e19, 6e-19, 0.999_079_317_603_851_4),
            (20.0, 1e12, 3e-11, 0.978_126_531_565_180_6),
        ] {
            crate::np_assert_allclose!(&[crate::betainc(a, b, x)], &[p], rtol = 1e-13);
        }
    }

    #[test]
    fn test_betainc_edge_cases() {
        for &((a, b, x), expected) in &[
            ((0.0, 0.0, 0.0), f64::NAN),
            ((0.0, 0.0, 0.5), f64::NAN),
            ((0.0, 0.0, 1.0), f64::NAN),
            ((f64::INFINITY, f64::INFINITY, 0.0), f64::NAN),
            ((f64::INFINITY, f64::INFINITY, 0.5), f64::NAN),
            ((f64::INFINITY, f64::INFINITY, 1.0), f64::NAN),
            ((0.0, 1.0, 0.0), 0.0),
            ((0.0, 1.0, 0.5), 1.0),
            ((0.0, 1.0, 1.0), 1.0),
            ((1.0, 0.0, 0.0), 0.0),
            ((1.0, 0.0, 0.5), 0.0),
            ((1.0, 0.0, 1.0), 1.0),
            ((0.0, f64::INFINITY, 0.0), 0.0),
            ((0.0, f64::INFINITY, 0.5), 1.0),
            ((0.0, f64::INFINITY, 1.0), 1.0),
            ((f64::INFINITY, 0.0, 0.0), 0.0),
            ((f64::INFINITY, 0.0, 0.5), 0.0),
            ((f64::INFINITY, 0.0, 1.0), 1.0),
            ((1.0, f64::INFINITY, 0.0), 0.0),
            ((1.0, f64::INFINITY, 0.5), 1.0),
            ((1.0, f64::INFINITY, 1.0), 1.0),
            ((f64::INFINITY, 1.0, 0.0), 0.0),
            ((f64::INFINITY, 1.0, 0.5), 0.0),
            ((f64::INFINITY, 1.0, 1.0), 1.0),
        ] {
            let observed = crate::betainc(a, b, x);
            crate::np_assert_equal!(&[observed], &[expected]);
        }
    }

    #[test]
    fn test_betainc_gh22682() {
        // scipy/scipy#22682: betainc returned incorrect results for tiny
        // single precision inputs. test that this is resolved
        let rtol = 1e-15;
        for &(a, b, x, reference) in &[
            (1e-20, 1e-21, 0.5, 0.090_909_090_909_090_9),
            (1e-15, 1e-16, 0.5, 0.090_909_090_909_090_91),
        ] {
            let res = crate::betainc(a, b, x);
            crate::np_assert_allclose!(&[res], &[reference], rtol = rtol);
        }
    }
}
