//! The complex `ndtr`, `log_ndtr`, `poisson_binom_{pmf,cdf}_all`, `cdf_cvm_inf`, and `cdf_cvm` are
//! translated into pure Rust from `xsf/stats.h` (xsf v0.2.2).
//!
//! Unlike in xsf (scipy/xsf#340), `cdf_cvm_inf(x)` returns `1` for `x > f64::MAX / 16` (incl.
//! `inf`), and `cdf_cvm(NaN, n)` returns `NaN`, instead of looping forever. Like in
//! `scipy.stats._hypotests`, the series also stop at a NaN term, so that e.g.
//! `cdf_cvm(1e4, 1_000_000)` returns `NaN` instead of looping forever.
//!
//! The `wilcoxon_{pmf,cdf}_all` and `wilcoxon_{pmf,cdf,sf}` functions are translated from
//! `xsf/stats.h` in xsf main. Unlike in xsf, a `k` that is out of range for `i64` (e.g. `inf`)
//! saturates, instead of being undefined behavior in the cast to `long long`.

use core::f64::consts::{FRAC_1_PI, FRAC_1_SQRT_2, LN_2, PI, TAU};

use num_complex::{Complex64, c64};

use crate::utils::{LN_PI, LN_TAU, SQRT_PI, SQRT2OPI};
use crate::xsf::bessel::{cyl_bessel_k_f64, cyl_bessel_ke_f64};
use crate::xsf::cephes::{betainc_xy, gamma, lgam, lgam1p};

pub trait StatsArg: crate::sealed::Sealed {
    fn ndtr(self) -> Self;
    fn log_ndtr(self) -> Self;
}

impl StatsArg for f64 {
    #[inline]
    fn ndtr(self) -> Self {
        crate::xsf::cephes::ndtr(self)
    }

    #[inline]
    fn log_ndtr(self) -> Self {
        log_ndtr_f64(self)
    }
}

impl StatsArg for Complex64 {
    #[inline]
    fn ndtr(self) -> Self {
        ndtr_c64(self)
    }

    #[inline]
    fn log_ndtr(self) -> Self {
        log_ndtr_c64(self)
    }
}

fn ndtr_c64(z: Complex64) -> Complex64 {
    crate::xsf::faddeeva::erfc_c64(-z * FRAC_1_SQRT_2, 0.0) * 0.5
}

/*
 * Log of the CDF of the normal distribution for double x.
 *
 * Let F(x) be the CDF of the standard normal distribution.
 * This implementation of log(F(x)) is based on the identities
 *
 *   F(x) = erfc(-x/√2)/2
 *        = 1 - erfc(x/√2)/2
 *
 * We use the first formula for x < -1, with erfc(z) replaced
 * by erfcx(z)*exp(-z**2) to ensure high precision for large
 * negative values when we take the logarithm:
 *
 *   log F(x) = log(erfc(-x/√2)/2)
 *            = log(erfcx(-x/√2)/2)*exp(-x**2/2))
 *            = log(erfcx(-x/√2)/2) - x**2/2
 *
 * For x >= -1, we use the second formula for F(x):
 *
 *   log F(x) = log(1 - erfc(x/√2)/2)
 *            = log1p(-erfc(x/√2)/2)
 */
fn log_ndtr_f64(x: f64) -> f64 {
    let t = x * FRAC_1_SQRT_2;
    if x < -1.0 {
        (crate::xsf::faddeeva::erfcx(-t) / 2.0).ln() - t * t
    } else {
        (-crate::xsf::cephes::erfc(t) / 2.0).ln_1p()
    }
}

/*
 * Log of the normal CDF for complex arguments.
 *
 * This is equivalent to log(ndtr(z)), but is more robust to overflow at $z\to\infty$.
 * This implementation uses $\erfc(z) = \exp(-z^2) w(iz)$ taking special care to select
 * the principal branch of the log function log( exp(-z^2) w(i z) )
 */
#[allow(non_snake_case)]
fn log_ndtr_c64(mut z: Complex64) -> Complex64 {
    if z.re > 6.0 {
        // Underflow. Close to the real axis, expand the log in log(1 - ndtr(-z)).
        let w = crate::xsf::faddeeva::erfc_c64(z * FRAC_1_SQRT_2, 0.0) * -0.5;
        if w.norm() < 1e-8 {
            return w;
        }
    }

    z *= -FRAC_1_SQRT_2;
    let x = z.re;
    let y = z.im;

    /* Compute the principal branch of $log(exp(-z^2))$, using the fact that
     * $log(e^t) = log|e^t| + i Arg(e^t)$, and that if $t = r + is$, then
     * $e^t = e^r (\cos(s) + i \sin(s))$.
     */
    let mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    let mIm_z2 = -2.0 * x * y; // Im(-z^2)

    let mut im = mIm_z2 % TAU;
    if im > PI {
        im -= TAU;
    }

    let val1 = c64(mRe_z2, im);

    let val2 = crate::xsf::log::clog(crate::xsf::faddeeva::w(c64(-y, x), 0.0));
    let result = val1 + val2 - LN_2;

    /* Again, select the principal branch: log(z) = log|z| + i arg(z), thus
     * the imaginary part of the result should belong to [-pi, pi].
     */
    im = result.im;
    if im >= PI {
        im -= TAU;
    }
    if im < -PI {
        im += TAU;
    }

    c64(result.re, im)
}

// Student's t

/// Student's t distribution cumulative distribution function
///
/// Rust implementation of [`scipy.special.stdtr`][stdtr].
///
/// # Notes
///
/// The CDF is given by
///
/// $$
/// \begin{align*}
/// F_\nu(x)
/// \&= \frac{1}{\nu \B \left( \frac{1}{2}, \frac{\nu}{2} \right)}
///     \int_{-\infty}^x g(u)^{\frac{\nu+1}{2}} \dd u \\\\
/// \&= 1 - \frac{1}{2} I_{g(x)}\left( \frac{1}{2}, \frac{\nu}{2} \right) ,
/// \end{align*}
/// $$
///
/// with $g(\square) = \frac{\nu}{\square^2 + \nu}$ a helper function, $\B(\cdot, \cdot)$ the Beta
/// function, and $I_\square(\cdot,\cdot)$ the regularized incomplete Beta function.
///
/// # See also
/// - [`stdtri`]: Inverse of the CDF
/// - [`beta`](crate::beta): Beta function
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function
///
/// [stdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtr.html
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
pub fn stdtr(nu: f64, x: f64) -> f64 {
    if nu <= 0.0 {
        f64::NAN
    } else if x == 0.0 {
        0.5
    } else if nu > 1e150 {
        // The terms of the continued fractions in `incbet` (~nu^2) overflow, and the 1/nu
        // correction phi(x) (x^3 + x) / (4 nu) is negligible, e.g. stdtr(1e200, -2.33) was 0.5
        // instead of 0.0099
        ndtr(x)
    } else if nu == 1.0 {
        0.5 + FRAC_1_PI * x.atan()
    } else if nu == 2.0 {
        f64::midpoint(1.0, x / (2.0 + x * x).sqrt())
    } else {
        let t2 = x * x;
        let w = nu / (nu + t2);
        let p = 0.5
            * if w < f64::MIN_POSITIVE && nu < 2.0 && x.is_finite() {
                // `w` underflows (and `t2` may overflow), so use the leading term of the power
                // series I_w(a, 1/2) ~ w^a / (a B(a, 1/2)) in log space, with w ~ nu / t^2 and
                // a B(a, 1/2) = Gamma(a + 1) sqrt(pi) / Gamma(a + 1/2), e.g. stdtr(0.0016, -3.8e274)
                // was 0 instead of 0.1807. For nu >= 2 the tail probability is subnormal at most.
                let a = 0.5 * nu;
                let lnw = nu.ln() - 2.0 * x.abs().ln();
                (a * lnw - lgam1p(a) + lgam(a + 0.5) - 0.5 * LN_PI).exp()
            } else {
                // `1 - w` is passed as well (see `betainc_xy`), e.g. stdtr(7.8e5, 5.9e-8) was 0.5
                // instead of 0.5 + 2.4e-8
                betainc_xy(0.5 * nu, 0.5, w, t2 / (nu + t2))
            };
        if x >= 0.0 { 1.0 - p } else { p }
    }
}

/// Inverse of [`stdtr`]
///
/// Rust implementation of [`scipy.special.stdtrit`][stdtrit] with comparable or better accuracy.
///
/// [stdtrit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtrit.html
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
pub fn stdtri(nu: f64, p: f64) -> f64 {
    if nu <= 0.0 || !(0.0..=1.0).contains(&p) {
        f64::NAN
    } else if p == 0.0 {
        f64::NEG_INFINITY
    } else if p == 1.0 {
        f64::INFINITY
    } else if p == 0.5 {
        0.0
    } else if nu.is_infinite() {
        ndtri(p)
    } else if nu == 1.0 {
        (PI * (p - 0.5)).tan()
    } else if nu == 2.0 {
        (p - 0.5) * (2.0 / (p * (1.0 - p))).sqrt()
    } else {
        let a = 0.5 * nu;
        let y = 2.0 * p.min(1.0 - p);
        if nu <= 1e3 && p > 0.25 && p < 0.75 {
            // Solve I_x(1/2, a) = |1 - 2p| (exact) for x = t^2 / (nu + t^2) = 1 - w, which is more
            // accurate than computing `1 - w` if x < 1/2, e.g. stdtri(10, 0.5 + 1e-10) was
            // 3.3e-8 instead of 2.6e-10
            let x = crate::betaincinv(0.5, a, (1.0 - 2.0 * p).abs());
            if x < 0.5 {
                let t = (nu * x / (1.0 - x)).sqrt();
                return if p < 0.5 { -t } else { t };
            }
        }
        let w = crate::betaincinv(a, 0.5, y);
        let mut t = (nu * (1.0 - w) / w).sqrt();
        if t.is_infinite() {
            // For small nu, `w` can underflow (or `nu / w` overflow), so invert the leading term
            // of the power series I_w(a, 1/2) ~ w^a / (a B(a, 1/2)) (see `stdtr`) instead:
            //
            // - t = sqrt(nu / w) = (y c)^(-1/2a), with c = a B(a, 1/2) / nu^a
            // - computed as a square, since y^(-1/4a) can't overflow if t doesn't
            //
            // e.g. stdtri(0.0017, 0.79) was inf instead of 8.6e219
            let e = -0.25 / a;
            let c = a * crate::beta(a, 0.5) / nu.powf(a);
            t = (y.powf(e) * c.powf(e)).powi(2);
        }
        if (nu > 1e3 || w > 0.5) && y >= 2.0 * f64::MIN_POSITIVE {
            // `1 - w` is inaccurate if w > 1/2, and `betaincinv` too for huge `a`, so start from
            // the Cornish-Fisher expansion if nu > 1e3, and (unless that's already exact) solve
            // stdtr(nu, -t) = y / 2 with Newton's method on ln(stdtr), e.g. stdtri(1e20, 0.01)
            // was -105.4 instead of -2.326
            let q = 0.5 * y;
            let mut err = 1.0;
            if nu > 1e3 {
                (t, err) = stdtri_cornish_fisher(nu, ndtri(q));
                t = -t;
            }
            let c = nu.sqrt() * crate::beta(a, 0.5);
            for _ in 0..20 {
                if err <= 4.0 * f64::EPSILON * t {
                    break;
                }
                let f = stdtr(nu, -t);
                let pdf = (-(a + 0.5) * (t * t / nu).ln_1p()).exp() / c;
                let dt = (f / q).ln() * f / pdf;
                t += dt;
                err = dt.abs();
            }
        }
        if p < 0.5 { -t } else { t }
    }
}

/// Cornish-Fisher expansion of the Student t quantile in terms of the normal quantile `z`, up to
/// order 1 / nu^4 (A&S 26.7.5), and the magnitude of that last term
fn stdtri_cornish_fisher(nu: f64, z: f64) -> (f64, f64) {
    let z2 = z * z;
    let g1 = (z2 + 1.0) / 4.0;
    let g2 = ((5.0 * z2 + 16.0) * z2 + 3.0) / 96.0;
    let g3 = (((3.0 * z2 + 19.0) * z2 + 17.0) * z2 - 15.0) / 384.0;
    let g4 = ((((79.0 * z2 + 776.0) * z2 + 1482.0) * z2 - 1920.0) * z2 - 945.0) / 92160.0;
    let t = z * (1.0 + (g1 + (g2 + (g3 + g4 / nu) / nu) / nu) / nu);
    (t, (z * g4 / nu.powi(4)).abs())
}

// Normal

/// CDF of the standard normal distribution, $\Phi(z)$
///
/// Corresponds to [`scipy.special.ndtr`][ndtr].
///
/// [ndtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtr.html
///
/// # Notes
///
/// The CDF is given by
///
/// $$
/// \begin{align*}
/// \Phi(z)
/// &= {1 \over \sqrt{2\pi}} \int_{-\infty}^z e^{-t^2 / 2} \dd t \\\\
/// &= {1 \over 2} + {1 \over 2} \erf \left( {z \over \sqrt{2}} \right) ,
/// \end{align*}
/// $$
///
/// with $\erf(z)$ the [error function](crate::erf).
///
/// # See also
/// - [`log_ndtr`]: $\ln \Phi(z)$
/// - [`ndtri`]: Normal quantile function $\Phi^{-1}(z)$, a.k.a. the probit function
/// - [`erf`](crate::erf): Error function $\erf(z)$
#[must_use]
#[inline]
pub fn ndtr<T: StatsArg>(z: T) -> T {
    z.ndtr()
}

/// Logarithm of [`ndtr`], $\ln \Phi(z)$
///
/// Corresponds to [`scipy.special.log_ndtr`][log_ndtr].
///
/// [log_ndtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_ndtr.html
///
/// # See also
/// - [`ndtr`]: Normal distribution function $\Phi(z)$
/// - [`erf`](crate::erf): Error function $\erf(z)$
/// - [`ndtri_exp`](crate::ndtri_exp): Inverse of `log_ndtr`
#[must_use]
#[inline]
pub fn log_ndtr<T: StatsArg>(z: T) -> T {
    z.log_ndtr()
}

/// Inverse of [`ndtr`], the probit function $\Phi^{-1}(p)$
///
/// Corresponds to [`scipy.special.ndtri`][ndtri].
///
/// [ndtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri.html
///
/// # Notes
///
/// The normal quantile function (probit function) is given by
///
/// $$
/// \Phi^{-1}(p) = \sqrt{2}\\,\erf^{-1}(2p-1) ,
/// $$
///
/// with $\erf^{-1}(z)$ the [inverse error function](crate::erfinv) and $p \in (0,1)$.
///
/// # See also
/// - [`ndtr`]: CDF of the standard normal distribution, $\Phi(z)$
/// - [`crate::erfinv`]: Inverse error function $\erf^{-1}(z)$
#[doc(alias = "probit")]
#[must_use]
#[inline]
pub fn ndtri(x: f64) -> f64 {
    crate::xsf::cephes::ndtri(x)
}

/// Mean of the normal distribution, given the CDF value `p`, the standard deviation `std`, and the
/// quantile `x`
///
/// Corresponds to [`scipy.special.nrdtrimn`][nrdtrimn], and `xsf::nrdtrimn` in xsf main.
///
/// [nrdtrimn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nrdtrimn.html
///
/// # See also
/// - [`ndtri`]: Inverse of [`ndtr`]
/// - [`nrdtrisd`]: Standard deviation of the normal distribution
#[must_use]
#[inline]
pub fn nrdtrimn(p: f64, std: f64, x: f64) -> f64 {
    if std.is_nan() || std <= 0.0 {
        return f64::NAN;
    }
    if p.is_nan() || p <= 0.0 || p >= 1.0 {
        return f64::NAN;
    }
    if x.is_nan() {
        return f64::NAN;
    }
    x - std * crate::xsf::cephes::ndtri(p)
}

/// Standard deviation of the normal distribution, given the mean `mean`, the CDF value `p`, and
/// the quantile `x`
///
/// Corresponds to [`scipy.special.nrdtrisd`][nrdtrisd], and `xsf::nrdtrisd` in xsf main.
///
/// [nrdtrisd]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nrdtrisd.html
///
/// # See also
/// - [`ndtri`]: Inverse of [`ndtr`]
/// - [`nrdtrimn`]: Mean of the normal distribution
#[must_use]
#[inline]
pub fn nrdtrisd(mean: f64, p: f64, x: f64) -> f64 {
    if mean.is_nan() || p.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if p <= 0.0 || p >= 1.0 {
        return f64::NAN;
    }
    (x - mean) / crate::xsf::cephes::ndtri(p)
}

fn von_mises_cdf_series(k: f64, x: f64, p: u32) -> f64 {
    let s = x.sin();
    let c = x.cos();
    let mut sn = (f64::from(p) * x).sin();
    let mut cn = (f64::from(p) * x).cos();
    let mut r = 0.0;
    let mut v = 0.0;
    for n in (1..p).rev() {
        let n = f64::from(n);
        let sn_new = sn * c - cn * s;
        let cn_new = cn * c + sn * s;
        sn = sn_new;
        cn = cn_new;
        r = k / (2.0 * n + k * r);
        v = r * (sn / n + v);
    }
    0.5 + x / (2.0 * PI) + v / PI
}

fn von_mises_cdf_normalapprox(k: f64, x: f64) -> f64 {
    let b = SQRT2OPI / crate::xsf::cephes::i0e(k);
    let z = b * (x / 2.0).sin();
    crate::xsf::cephes::ndtr(z)
}

/// CDF of the von Mises distribution with concentration `k`, extended periodically over `x`
///
/// For `k < 50`, this uses the backward-recursion series method of Hill \[1\], Algorithm 518,
/// with constants chosen for about 12 decimal digits. For `k >= 50`, this switches to a
/// large-concentration normal approximation.
///
/// Corresponds to `xsf::von_mises_cdf` in xsf main, which is used by
/// [`scipy.stats.vonmises`][vonmises].
///
/// [vonmises]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.vonmises.html
///
/// \[1\] G. W. Hill, "Algorithm 518: Incomplete Bessel Function I0.
///     The Von Mises Distribution \[S14\]", ACM Transactions on Mathematical
///     Software, 3(3), 279-284, 1977.
///     DOI: <https://doi.org/10.1145/355744.355753>
#[must_use]
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn von_mises_cdf(k: f64, x: f64) -> f64 {
    if k < 0.0 {
        return f64::NAN;
    }

    let ix = (x / (2.0 * PI)).round();
    let x = x - ix * 2.0 * PI;

    // These values should give 12 decimal digits
    let ck = 50.0;
    let (a1, a2, a3, a4) = (28.0, 0.5, 100.0, 5.0);
    let result = if k < ck {
        let p = (1.0 + a1 + a2 * k - a3 / (k + a4)) as u32;
        let result = von_mises_cdf_series(k, x, p);
        // `cxx::min(cxx::max(result, 0.0), 1.0)`, with `max(a, b) = a < b ? b : a` and
        // `min(a, b) = a < b ? a : b`, which maps NaN to 1 (e.g. for x = inf)
        let result = if result < 0.0 { 0.0 } else { result };
        if result < 1.0 { result } else { 1.0 }
    } else {
        von_mises_cdf_normalapprox(k, x)
    };
    result + ix
}

/// Owen's T function, $T(h, a)$
///
/// The function $T(h, a)$ gives the probability of the event $(X > h \wedge 0 < Y < a X)$ where
/// $X$ and $Y$ are independent standard normal random variables.
///
/// Corresponds to [`scipy.special.owens_t`][owens_t].
///
/// [owens_t]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.owens_t.html
///
/// # Definition
///
/// $$ T(h, a) = {1 \over 2\pi} \int_0^a {e^{-h^2 (1+x^2) / 2} \over 1+x^2} \dd x $$
///
/// # See also
/// - [`ndtr`]: CDF of the standard normal distribution, $\Phi(z)$
#[must_use]
#[inline]
pub fn owens_t(h: f64, a: f64) -> f64 {
    crate::xsf::cephes::owens_t(h, a)
}

// Kolmogorov

/// Survival function of the Kolmogorov distribution
///
/// Corresponds to [`scipy.special.kolmogorov`][kolmogorov].
///
/// [kolmogorov]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kolmogorov.html
///
/// # See also
/// - [`kolmogi`]: Inverse survival function
/// - [`kolmogc`]: Distribution function
/// - [`kolmogci`]: Quantile function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogorov(x: f64) -> f64 {
    crate::xsf::cephes::kolmogorov(x)
}

/// Inverse of [`kolmogorov`]
///
/// Corresponds to [`scipy.special.kolmogi`][kolmogi].
///
/// [kolmogi]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kolmogi.html
///
/// # See also
/// - [`kolmogorov`]: Survival function
/// - [`kolmogc`]: Distribution function
/// - [`kolmogci`]: Quantile function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogi(x: f64) -> f64 {
    crate::xsf::cephes::kolmogi(x)
}

/// CDF of the Kolmogorov distribution
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`kolmogorov`]: Survival function
/// - [`kolmogi`]: Inverse survival function
/// - [`kolmogci`]: Quantile function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogc(x: f64) -> f64 {
    crate::xsf::cephes::kolmogc(x)
}

/// Inverse of [`kolmogc`], the quantile function of the Kolmogorov distribution
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`kolmogorov`]: Survival function
/// - [`kolmogi`]: Inverse survival function
/// - [`kolmogc`]: Distribution function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogci(x: f64) -> f64 {
    crate::xsf::cephes::kolmogci(x)
}

/// Derivative of [`kolmogorov`]
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`kolmogorov`]: Survival function of the Kolmogorov distribution
/// - [`kolmogi`]: Inverse of the Kolmogorov distribution function
/// - [`kolmogc`]: CDF of the Kolmogorov distribution
/// - [`kolmogci`]: Inverse of the Kolmogorov CDF
#[must_use]
#[inline]
pub fn kolmogp(x: f64) -> f64 {
    crate::xsf::cephes::kolmogp(x)
}

// Kolmogorov-Smirnov

/// Survival function of the Kolmogorov-Smirnov distribution
///
/// Corresponds to [`scipy.special.smirnov`][smirnov].
///
/// [smirnov]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.smirnov.html
///
/// # See also
/// - [`smirnovi`]: Inverse survival function
/// - [`smirnovc`]: Distribution function
/// - [`smirnovp`]: Derivative of the survival function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnov(n: i32, x: f64) -> f64 {
    crate::xsf::cephes::smirnov(n, x)
}

/// Inverse of [`smirnov`]
///
/// Corresponds to [`scipy.special.smirnovi`][smirnovi].
///
/// [smirnovi]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.smirnovi.html
///
/// # See also
/// - [`smirnov`]: Survival function
/// - [`smirnovc`]: Distribution function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnovi(n: i32, q: f64) -> f64 {
    crate::xsf::cephes::smirnovi(n, q)
}

/// CDF of the Kolmogorov-Smirnov distribution
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`smirnov`]: Survival function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnovc(n: i32, x: f64) -> f64 {
    crate::xsf::cephes::smirnovc(n, x)
}

/// Inverse of [`smirnovc`]
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`smirnovc`]: Distribution function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnovci(n: i32, p: f64) -> f64 {
    crate::xsf::cephes::smirnovci(n, p)
}

/// Derivative of [`smirnov`]
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`smirnov`]: Survival function
#[must_use]
#[inline]
pub fn smirnovp(n: i32, x: f64) -> f64 {
    crate::xsf::cephes::smirnovp(n, x)
}

// Tukey-Lambda

/// CDF of the Tukey-Lambda distribution
///
/// Corresponds to [`scipy.special.tklmbda`][tklmbda].
///
/// [tklmbda]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.tklmbda.html
#[doc(alias = "tklmbda")]
#[must_use]
#[inline]
pub fn tukeylambdacdf(x: f64, lambda: f64) -> f64 {
    crate::xsf::cephes::tukeylambdacdf(x, lambda)
}

// Chi-squared

/// CDF of the Chi-squared distribution
///
/// Corresponds to [`scipy.special.chdtr`][chdtr].
///
/// [chdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtr.html
///
/// # See also
/// - [`chdtrc`]: Survival function
/// - [`chdtri`]: Quantile function
#[must_use]
#[inline]
pub fn chdtr(df: f64, x: f64) -> f64 {
    crate::xsf::cephes::chdtr(df, x)
}

/// Survival function of the Chi-squared distribution
///
/// Corresponds to [`scipy.special.chdtrc`][chdtrc].
///
/// [chdtrc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtrc.html
///
/// # See also
/// - [`chdtr`]: Distribution function
/// - [`chdtri`]: Quantile function
#[must_use]
#[inline]
pub fn chdtrc(df: f64, x: f64) -> f64 {
    crate::xsf::cephes::chdtrc(df, x)
}

/// Quantile function of the Chi-squared distribution
///
/// Corresponds to [`scipy.special.chdtri`][chdtri].
///
/// [chdtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtri.html
///
/// # See also
/// - [`chdtr`]: Distribution function
/// - [`chdtrc`]: Survival function
#[must_use]
#[inline]
pub fn chdtri(df: f64, y: f64) -> f64 {
    crate::xsf::cephes::chdtri(df, y)
}

// F

/// CDF of the F distribution
///
/// Corresponds to [`scipy.special.fdtr`][fdtr].
///
/// [fdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtr.html
///
/// # See also
/// - [`fdtrc`]: Survival function
/// - [`fdtri`]: Quantile function
#[must_use]
#[inline]
pub fn fdtr(a: f64, b: f64, x: f64) -> f64 {
    crate::xsf::cephes::fdtr(a, b, x)
}

/// Survival function of the F distribution
///
/// Corresponds to [`scipy.special.fdtrc`][fdtrc].
///
/// [fdtrc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtrc.html
///
/// # See also
/// - [`fdtr`]: Distribution function
/// - [`fdtri`]: Quantile function
#[must_use]
#[inline]
pub fn fdtrc(a: f64, b: f64, x: f64) -> f64 {
    crate::xsf::cephes::fdtrc(a, b, x)
}

/// Quantile function of the F distribution
///
/// Corresponds to [`scipy.special.fdtri`][fdtri].
///
/// [fdtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fdtri.html
///
/// # See also
/// - [`fdtr`]: Distribution function
/// - [`fdtrc`]: Survival function
#[must_use]
#[inline]
pub fn fdtri(a: f64, b: f64, y: f64) -> f64 {
    crate::xsf::cephes::fdtri(a, b, y)
}

// Gamma

/// CDF of the Gamma distribution
///
/// Corresponds to [`scipy.special.gdtr`][gdtr].
///
/// [gdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtr.html
///
/// # See also
/// - [`gdtrc`]: Survival function
/// - [`gdtrib`](crate::gdtrib): Inverse of the CDF w.r.t. `b`
#[must_use]
#[inline]
pub fn gdtr(a: f64, b: f64, x: f64) -> f64 {
    crate::xsf::cephes::gdtr(a, b, x)
}

/// Survival function of the Gamma distribution
///
/// Corresponds to [`scipy.special.gdtrc`][gdtrc].
///
/// [gdtrc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtrc.html
///
/// # See also
/// - [`gdtr`]: Cumulative distribution function
/// - [`gdtrib`](crate::gdtrib): Inverse of the CDF w.r.t. `b`
#[must_use]
#[inline]
pub fn gdtrc(a: f64, b: f64, x: f64) -> f64 {
    crate::xsf::cephes::gdtrc(a, b, x)
}

// Poisson

/// Poisson distribution function
///
/// Corresponds to [`scipy.special.pdtr`][pdtr].
///
/// [pdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtr.html
///
/// # See also
/// - [`pdtri`]: Quantile function
/// - [`pdtrc`]: Survival function
#[must_use]
#[inline]
pub fn pdtr(k: f64, m: f64) -> f64 {
    crate::xsf::cephes::pdtr(k, m)
}

/// Poisson quantile function
///
/// Corresponds to [`scipy.special.pdtri`][pdtri].
///
/// [pdtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtri.html
///
/// # See also
/// - [`pdtr`]: Cumulative distribution function
/// - [`pdtrc`]: Survival function
#[must_use]
#[inline]
pub fn pdtri(k: i32, y: f64) -> f64 {
    crate::xsf::cephes::pdtri(k, y)
}

/// Poisson survival function
///
/// Corresponds to [`scipy.special.pdtrc`][pdtrc].
///
/// [pdtrc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pdtrc.html
///
/// # See also
/// - [`pdtr`]: Cumulative distribution function
/// - [`pdtri`]: Quantile function
#[must_use]
#[inline]
pub fn pdtrc(k: f64, m: f64) -> f64 {
    crate::xsf::cephes::pdtrc(k, m)
}

// Binomial

/// Binomial distribution function
///
/// Corresponds to [`scipy.special.bdtr`][bdtr].
///
/// [bdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtr.html
///
/// # See also
/// - [`bdtrc`]: Survival function
/// - [`bdtri`]: Quantile function
#[must_use]
#[inline]
pub fn bdtr(k: f64, n: i32, p: f64) -> f64 {
    crate::xsf::cephes::bdtr(k, n, p)
}

/// Binomial survival function
///
/// Corresponds to [`scipy.special.bdtrc`][bdtrc].
///
/// [bdtrc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtrc.html
///
/// # See also
/// - [`bdtr`]: Distribution function
/// - [`bdtri`]: Quantile function
#[must_use]
#[inline]
pub fn bdtrc(k: f64, n: i32, p: f64) -> f64 {
    crate::xsf::cephes::bdtrc(k, n, p)
}

/// Binomial quantile function
///
/// Corresponds to [`scipy.special.bdtri`][bdtri].
///
/// [bdtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bdtri.html
///
/// # See also
/// - [`bdtr`]: Distribution function
/// - [`bdtrc`]: Survival function
#[must_use]
#[inline]
pub fn bdtri(k: f64, n: i32, y: f64) -> f64 {
    crate::xsf::cephes::bdtri(k, n, y)
}

// Poisson binomial

#[allow(clippy::needless_range_loop)]
fn poisson_binom_pmf_all_impl(p: &[f64], res: &mut [f64]) {
    let n = p.len();

    // Initialize output array with zeros.
    for i in 0..=n {
        res[i] = 0.0;
    }

    if n == 0 {
        res[0] = 1.0;
        return;
    }

    res[0] = 1.0 - p[0];
    res[1] = p[0];

    for i in 1..n {
        let p_i = p[i];
        let q_i = 1.0 - p_i;
        for j in (1..=i + 1).rev() {
            let tmp = res[j - 1] * p_i;
            res[j - 1] *= q_i;
            res[j] += tmp;
        }
    }
}

/// Full PMF of a Poisson-Binomial distribution
///
/// Used internally by [`scipy.stats.poisson_binom`][poisson_binom].
///
/// [poisson_binom]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.poisson_binom.html
///
/// `p` contains the success probabilities for the $n$ Bernoulli trials. The result at index `k`
/// contains the probability of observing $k$ successes for $k$ from $0$ to $n$.
///
/// # See also
/// - [`poisson_binom_cdf_all`]: Full CDF of a Poisson-Binomial distribution
#[must_use]
#[inline]
pub fn poisson_binom_pmf_all(p: &[f64]) -> Vec<f64> {
    let mut res = vec![0.0; p.len() + 1];
    poisson_binom_pmf_all_impl(p, &mut res);
    res
}

/// Full CDF of a Poisson-Binomial distribution
///
/// Used internally by [`scipy.stats.poisson_binom`][poisson_binom].
///
/// [poisson_binom]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.poisson_binom.html
///
/// `p` contains the success probabilities for the $n$ Bernoulli trials. The result at index `k`
/// contains the probability of observing at most $k$ successes for $k$ from $0$ to $n$.
///
/// # See also
/// - [`poisson_binom_pmf_all`]: Full PMF of a Poisson-Binomial distribution
#[must_use]
#[inline]
pub fn poisson_binom_cdf_all(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut res = vec![0.0; n + 1];

    poisson_binom_pmf_all_impl(p, &mut res);
    for i in 1..n {
        // `std::min(res(i) + res(i - 1), T(1.0))`, which propagates NaN
        let r = res[i] + res[i - 1];
        res[i] = if 1.0 < r { 1.0 } else { r };
    }
    res[n] = 1.0;
    res
}

fn take_from_pmf(pmf: &[f64], k: i64) -> f64 {
    let size = pmf.len();
    match usize::try_from(k) {
        Ok(k) if k < size => pmf[k],
        _ => 0.0,
    }
}

fn take_from_discrete_cdf(cdf: &[f64], k: i64) -> f64 {
    // Unlike in xsf, this returns 1 instead of reading out of bounds if `cdf` is empty.
    let size = cdf.len();
    if k < 0 {
        return 0.0;
    }
    match usize::try_from(k) {
        Ok(k) if k < size.saturating_sub(1) => cdf[k],
        _ => 1.0,
    }
}

// Wilcoxon signed-rank

/// Shared implementation used by PMF and CDF table generation.
#[allow(clippy::manual_midpoint, clippy::needless_range_loop)]
fn wilcoxon_pmf_all_impl(n: usize, res: &mut [f64]) {
    let max_k = n * (n + 1) / 2;

    for i in 0..=max_k {
        res[i] = 0.0;
    }
    res[0] = 1.0;

    for k in 1..=n {
        let curr_max = k * (k + 1) / 2;
        // Iterate backwards so res(i - k) still holds its previous-step value.
        for i in (0..=curr_max).rev() {
            res[i] = 0.5 * (res[i] + if i >= k { res[i - k] } else { 0.0 });
        }
    }
}

/// Full PMF of the Wilcoxon signed-rank statistic
///
/// Used internally by [`scipy.stats.wilcoxon`][wilcoxon].
///
/// [wilcoxon]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
///
/// The result has length $n (n + 1) / 2 + 1$, and contains the probability that the signed-rank
/// statistic of a sample of size $n$ equals $k$ at index `k`.
///
/// Corresponds to `xsf::wilcoxon_pmf_all` in xsf main.
///
/// # See also
/// - [`wilcoxon_cdf_all`]: Full CDF of the Wilcoxon signed-rank statistic
/// - [`wilcoxon_pmf`]: Look up the PMF in this table
#[must_use]
#[inline]
pub fn wilcoxon_pmf_all(n: usize) -> Vec<f64> {
    let mut res = vec![0.0; n * (n + 1) / 2 + 1];
    wilcoxon_pmf_all_impl(n, &mut res);
    res
}

/// Full CDF of the Wilcoxon signed-rank statistic
///
/// Used internally by [`scipy.stats.wilcoxon`][wilcoxon].
///
/// [wilcoxon]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
///
/// The result has length $n (n + 1) / 2 + 1$, and contains the probability that the signed-rank
/// statistic of a sample of size $n$ is at most $k$ at index `k`.
///
/// Corresponds to `xsf::wilcoxon_cdf_all` in xsf main.
///
/// # See also
/// - [`wilcoxon_pmf_all`]: Full PMF of the Wilcoxon signed-rank statistic
/// - [`wilcoxon_cdf`]: Look up the CDF in this table
/// - [`wilcoxon_sf`]: Look up the survival function in this table
#[must_use]
#[inline]
pub fn wilcoxon_cdf_all(n: usize) -> Vec<f64> {
    let out_size = n * (n + 1) / 2 + 1;
    let mut res = vec![0.0; out_size];

    wilcoxon_pmf_all_impl(n, &mut res);
    for i in 1..out_size - 1 {
        // `cxx::min(res(i) + res(i - 1), T(1.0))`, with `min(a, b) = a < b ? a : b`
        let r = res[i] + res[i - 1];
        res[i] = if r < 1.0 { r } else { 1.0 };
    }
    res[out_size - 1] = 1.0;
    res
}

/// PMF of the Wilcoxon signed-rank statistic, from a precomputed table
///
/// After truncating `k` towards zero, this returns `pmf[k]`, or `0` if `k` is out of bounds.
/// `pmf` is the table returned by [`wilcoxon_pmf_all`].
///
/// Corresponds to `xsf::wilcoxon_pmf` in xsf main.
///
/// # See also
/// - [`wilcoxon_pmf_all`]: Full PMF of the Wilcoxon signed-rank statistic
#[must_use]
#[inline]
#[allow(clippy::cast_possible_truncation)]
pub fn wilcoxon_pmf(k: f64, pmf: &[f64]) -> f64 {
    if k.is_nan() {
        return f64::NAN;
    }
    take_from_pmf(pmf, k as i64)
}

/// CDF of the Wilcoxon signed-rank statistic, from a precomputed table
///
/// After truncating `k` towards zero, this returns `cdf[k]`, `0` if `k < 0`, or `1` if `k` is at
/// least the last index. `cdf` is the table returned by [`wilcoxon_cdf_all`].
///
/// Corresponds to `xsf::wilcoxon_cdf` in xsf main.
///
/// # See also
/// - [`wilcoxon_cdf_all`]: Full CDF of the Wilcoxon signed-rank statistic
/// - [`wilcoxon_sf`]: Survival function
#[must_use]
#[inline]
#[allow(clippy::cast_possible_truncation)]
pub fn wilcoxon_cdf(k: f64, cdf: &[f64]) -> f64 {
    if k.is_nan() {
        return f64::NAN;
    }
    take_from_discrete_cdf(cdf, k as i64)
}

/// Survival function of the Wilcoxon signed-rank statistic, from a precomputed table
///
/// After truncating `k` towards zero, this returns `1 - cdf[k - 1]`, `1` if `k <= 0`, or `0` if
/// `k` exceeds the last index. `cdf` is the table returned by [`wilcoxon_cdf_all`].
///
/// Corresponds to `xsf::wilcoxon_sf` in xsf main.
///
/// # See also
/// - [`wilcoxon_cdf_all`]: Full CDF of the Wilcoxon signed-rank statistic
/// - [`wilcoxon_cdf`]: Cumulative distribution function
#[must_use]
#[inline]
#[allow(clippy::cast_possible_truncation)]
pub fn wilcoxon_sf(k: f64, cdf: &[f64]) -> f64 {
    if k.is_nan() {
        return f64::NAN;
    }
    // Unlike in xsf, `k - 1` saturates instead of overflowing for `k = -inf`.
    1.0 - take_from_discrete_cdf(cdf, (k as i64).saturating_sub(1))
}

// Negative Binomial

/// Negative binomial distribution function
///
/// Corresponds to [`scipy.special.nbdtr`][nbdtr].
///
/// [nbdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtr.html
///
/// # See also
/// - [`nbdtrc`]: Survival function
/// - [`nbdtri`]: Quantile function
#[must_use]
#[inline]
pub fn nbdtr(k: i32, n: i32, p: f64) -> f64 {
    crate::xsf::cephes::nbdtr(k, n, p)
}

/// Negative binomial survival function
///
/// Corresponds to [`scipy.special.nbdtrc`][nbdtrc].
///
/// [nbdtrc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtrc.html
///
/// # See also
/// - [`nbdtr`]: Distribution function
/// - [`nbdtri`]: Quantile function
#[must_use]
#[inline]
pub fn nbdtrc(k: i32, n: i32, p: f64) -> f64 {
    crate::xsf::cephes::nbdtrc(k, n, p)
}

/// Negative binomial quantile function
///
/// Corresponds to [`scipy.special.nbdtri`][nbdtri].
///
/// [nbdtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.nbdtri.html
///
/// # See also
/// - [`nbdtr`]: Distribution function
/// - [`nbdtrc`]: Survival function
#[must_use]
#[inline]
pub fn nbdtri(k: i32, n: i32, p: f64) -> f64 {
    crate::xsf::cephes::nbdtri(k, n, p)
}

/// CDF of the Cramér-von Mises test statistic (infinite sample limit)
///
/// Corresponds to `scipy.stats._hypotests._cdf_cvm_inf`, which is used in
/// [`scipy.stats.cramervonmises`][f1] and [`scipy.stats.cramervonmises_2samp`][f2]
///
/// [f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cramervonmises.html
/// [f2]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cramervonmises_2samp.html
///
/// Accurate for practical hypothesis testing but not expected to be accurate for large values of
/// `x`, e.g. `x > 4`, when the cdf is very close to `1`.
#[must_use]
#[inline]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn cdf_cvm_inf(x: f64) -> f64 {
    // CDF of the Cramer-von Mises test statistic (infinite sample limit).
    // Accurate for practical hypothesis testing but not expected to be accurate
    // for large values of x, e.g. x > 4, when the cdf is very close to 1.
    if x.is_nan() {
        return x;
    }
    if 16.0 * x == f64::INFINITY {
        // (unlike in scipy, scipy/xsf#340) instead of an infinite loop, as `q` below would be `0`
        return 1.0;
    }
    if x <= 0.0 {
        return 0.0;
    }

    let mut total = 0.0;
    let mut k = 0;

    loop {
        let k_f = f64::from(k);
        let u = (lgam(k_f + 0.5) - lgam(k_f + 1.0)).exp() / (PI.powf(1.5) * x.sqrt());
        let y = 4.0 * k_f + 1.0;
        let q = (y * y) / (16.0 * x);
        let b = cyl_bessel_k_f64(0.25, q);
        let z = u * y.sqrt() * (-q).exp() * b;

        total += z;
        // (unlike in xsf, scipy/xsf#340) like in scipy, so that a NaN term doesn't cause an
        // infinite loop
        if !(z.abs() >= 1e-7) {
            break;
        }
        k += 1;
    }

    total
}

// Helpers for psi1 computations
#[allow(non_snake_case)]
fn psi1_G(y: f64) -> f64 {
    // Helper implementing the G-type Bessel combination appearing in equation 1.10
    // in Csörgő, S. and Faraway, J. (1996).
    let z = y * y / 4.0;
    let b = cyl_bessel_k_f64(0.25, z) + cyl_bessel_k_f64(0.75, z);
    (-z).exp() * (y / 2.0).powf(1.5) * b / SQRT_PI
}

#[allow(non_snake_case)]
fn psi1_H(y: f64) -> f64 {
    // Helper implementing the H-type Bessel combination appearing in equation 1.10
    // in Csörgő, S. and Faraway, J. (1996).
    let z = y * y / 4.0;
    let b = 2.0 * cyl_bessel_k_f64(0.25, z) + 3.0 * cyl_bessel_k_f64(0.75, z)
        - cyl_bessel_k_f64(1.25, z);
    (-z).exp() * (y / 2.0).powf(2.5) * b / SQRT_PI
}

fn psi1_term(k: i32, x: f64) -> f64 {
    // Compute the k-th term of the psi1 series expansion in equation 1.10
    // of Csörgő, S. and Faraway, J. (1996).
    let k_f = f64::from(k);
    let m = 2.0 * k_f + 1.0;
    let sx = 2.0 * x.sqrt();
    let y1 = x.powf(0.75);
    let y2 = x.powf(1.25);
    let gamma_kp1_2 = gamma(k_f + 0.5);
    let gamma_kp3_2 = gamma(k_f + 1.5);
    let e1 = m * gamma_kp1_2 * psi1_G((4.0 * k_f + 3.0) / sx) / (9.0 * y1);
    let e2 = gamma_kp1_2 * psi1_H((4.0 * k_f + 1.0) / sx) / (72.0 * y2);
    let e3 = 2.0 * (m + 2.0) * gamma_kp3_2 * psi1_H((4.0 * k_f + 5.0) / sx) / (12.0 * y2);
    let e4 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k_f + 1.0) / sx) / (144.0 * y1);
    let e5 = 7.0 * m * gamma_kp1_2 * psi1_G((4.0 * k_f + 5.0) / sx) / (144.0 * y1);
    e1 + e2 + e3 + e4 + e5
}

#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn psi1_mod(x: f64) -> f64 {
    // psi1 is defined in equation 1.10 in Csörgő, S. and Faraway, J. (1996).
    // This implements a modified version by excluding the term V(x) / 12
    // (here: cdf_cvm_inf(x) / 12) to avoid evaluating cdf_cvm_inf(x) twice in cdf_cvm.
    //
    // Implementation based on MAPLE code of Julian Faraway and R code of the
    // function pCvM in the package goftest (v1.1.1), permission granted
    // by Adrian Baddeley. Main difference in the implementation: the code
    // here keeps adding terms of the series until the terms are small enough.
    let mut k = 0;
    let mut tot = 0.0;
    loop {
        let gamma_kp1 = gamma(f64::from(k) + 1.0);
        let z = -psi1_term(k, x) / (PI * gamma_kp1);
        tot += z;
        // (unlike in xsf, scipy/xsf#340) like in scipy, so that a NaN term doesn't cause an
        // infinite loop
        if !(z.abs() >= 1e-7) {
            break;
        }
        k += 1;
    }
    tot
}

/// CDF of the Cramér-von Mises statistic for a finite sample size `n`
///
/// Corresponds to `scipy.stats._hypotests._cdf_cvm`, which is used in
/// [`scipy.stats.cramervonmises`][f1] and [`scipy.stats.cramervonmises_2samp`][f2]
///
/// [f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cramervonmises.html
/// [f2]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cramervonmises_2samp.html
///
/// See equation 1.8 in Csörgő, S. and Faraway, J. (1996) for finite samples, 1.2 for the asymptotic
/// CDF.
///
/// For finite `n`, the approximation is less accurate near the support boundaries `[1/(12*n), n/3]`
/// and for larger `x` values where the CDF is close to `1`.
/// The returned value is clipped to `[0, 1]`.
#[must_use]
#[inline]
pub fn cdf_cvm(x: f64, n: i32) -> f64 {
    if x.is_nan() {
        // (unlike in scipy, scipy/xsf#340) instead of an infinite loop
        return x;
    }
    // Calculate the CDF of the Cramér-von Mises statistic for a finite sample
    // size n. If n=-1, use the asymptotic CDF (n=inf).
    if n == -1 {
        return cdf_cvm_inf(x);
    }
    let lower = 1.0 / (12.0 * f64::from(n));
    let upper = f64::from(n) / 3.0;
    if x <= lower {
        return 0.0;
    }
    if x >= upper {
        return 1.0;
    }
    let cdf = cdf_cvm_inf(x) * (1.0 + lower) + psi1_mod(x) / f64::from(n);
    // `clamp` propagates NaN like `std::min(std::max(cdf, 0.0), 1.0)`
    cdf.clamp(0.0, 1.0)
}

/// Logarithm of the generalized inverse Gaussian probability density function
///
/// Corresponds to `xsf::geninvgauss_logpdf` in xsf main, and to the `logpdf` of
/// [`scipy.stats.geninvgauss`][geninvgauss] with `loc=0` and `scale=1`.
///
/// [geninvgauss]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.geninvgauss.html
///
/// # Definition
///
/// $$ f(x; p, b) = {x^{p-1} \over 2 K_p(b)} \exp\left(-{b \over 2} \left(x + {1 \over x}\right)\right) $$
///
/// for $x > 0$, where $K_p$ is the modified Bessel function of the second kind.
///
/// # See also
/// - [`geninvgauss_pdf`]: the density itself
#[must_use]
pub fn geninvgauss_logpdf(x: f64, p: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return f64::NEG_INFINITY;
    }

    let z = cyl_bessel_ke_f64(p, b);
    if z.is_infinite() {
        return f64::NAN;
    }

    -LN_2 - z.ln() + b + (p - 1.0) * x.ln() - b * (x + 1.0 / x) / 2.0
}

/// Generalized inverse Gaussian probability density function
///
/// Corresponds to `xsf::geninvgauss_pdf` in xsf main, and to the `pdf` of
/// [`scipy.stats.geninvgauss`][geninvgauss] with `loc=0` and `scale=1`.
///
/// [geninvgauss]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.geninvgauss.html
///
/// # See also
/// - [`geninvgauss_logpdf`]: its logarithm, and the definition
#[must_use]
#[inline]
pub fn geninvgauss_pdf(x: f64, p: f64, b: f64) -> f64 {
    geninvgauss_logpdf(x, p, b).exp()
}

fn genhyperbolic_log_norming_constant(p: f64, a: f64, b: f64) -> f64 {
    let t1 = (a + b) * (a - b);
    let t2 = p * 0.5 * t1.ln();
    let t3 = 0.5 * LN_TAU;
    let t4 = (p - 0.5) * a.ln();
    let t5 = t1.sqrt();
    let t6 = cyl_bessel_ke_f64(p, t5).ln() - t5;

    t2 - t3 - t4 - t6
}

/// Logarithm of the generalized hyperbolic probability density function
///
/// Corresponds to `xsf::genhyperbolic_logpdf` in xsf main, and to the `logpdf` of
/// [`scipy.stats.genhyperbolic`][genhyperbolic] with `loc=0` and `scale=1`.
///
/// [genhyperbolic]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genhyperbolic.html
///
/// # Definition
///
/// $$
/// f(x; p, a, b) = {(a^2 - b^2)^{p/2} \over \sqrt{2\pi} a^{p-1/2} K_p\left(\sqrt{a^2 - b^2}\right)}
/// e^{b x} {K_{p-1/2}\left(a \sqrt{1 + x^2}\right) \over \left(\sqrt{1 + x^2}\right)^{1/2-p}}
/// $$
///
/// for $|b| < a$, where $K_p$ is the modified Bessel function of the second kind.
///
/// # See also
/// - [`genhyperbolic_pdf`]: the density itself
#[must_use]
pub fn genhyperbolic_logpdf(x: f64, p: f64, a: f64, b: f64) -> f64 {
    let t1 = genhyperbolic_log_norming_constant(p, a, b);
    let t2 = 1.0_f64.hypot(x);
    let t3 = (p - 0.5) * t2.ln();
    let t4 = cyl_bessel_ke_f64(p - 0.5, a * t2).ln() - a * t2;
    let t5 = b * x;

    t1 + t3 + t4 + t5
}

/// Generalized hyperbolic probability density function
///
/// Corresponds to `xsf::genhyperbolic_pdf` in xsf main, and to the `pdf` of
/// [`scipy.stats.genhyperbolic`][genhyperbolic] with `loc=0` and `scale=1`.
///
/// [genhyperbolic]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genhyperbolic.html
///
/// # See also
/// - [`genhyperbolic_logpdf`]: its logarithm, and the definition
#[must_use]
#[inline]
pub fn genhyperbolic_pdf(x: f64, p: f64, a: f64, b: f64) -> f64 {
    genhyperbolic_logpdf(x, p, a, b).exp()
}
