//! The complex `ndtr`, `log_ndtr`, `poisson_binom_{pmf,cdf}_all`, `cdf_cvm_inf`, and `cdf_cvm` are
//! translated into pure Rust from `xsf/stats.h` (xsf v0.2.2).
//!
//! Unlike in xsf, `cdf_cvm_inf(x)` returns `1` for `x > f64::MAX / 16` (incl. `inf`), and
//! `cdf_cvm(NaN, n)` returns `NaN`, instead of looping forever. Like in
//! `scipy.stats._hypotests`, the series also stop at a NaN term, so that e.g.
//! `cdf_cvm(1e4, 1_000_000)` returns `NaN` instead of looping forever.
//!
//! The `wilcoxon_{pmf,cdf}_all` and `wilcoxon_{pmf,cdf,sf}` functions are translated from
//! `xsf/stats.h` in xsf main. Unlike in xsf, a `k` that is out of range for `i64` (e.g. `inf`)
//! saturates, instead of being undefined behavior in the cast to `long long`.

use core::f64::consts::{FRAC_1_PI, FRAC_1_SQRT_2, LN_2, PI, TAU};

use num_complex::{Complex64, c64};

use crate::utils::LN_PI;
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
    } else if nu.is_infinite() {
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
        if p < 0.5 { -t } else { t }
    }
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
    const SQRT2OPI: f64 = 0.797_884_560_802_865_4; // sqrt(2/pi)
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
        // (unlike in scipy) instead of an infinite loop, as `q` below would be `0`
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
        // (unlike in xsf) like in scipy, so that a NaN term doesn't cause an infinite loop
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
    (-z).exp() * (y / 2.0).powf(1.5) * b / PI.sqrt()
}

#[allow(non_snake_case)]
fn psi1_H(y: f64) -> f64 {
    // Helper implementing the H-type Bessel combination appearing in equation 1.10
    // in Csörgő, S. and Faraway, J. (1996).
    let z = y * y / 4.0;
    let b = 2.0 * cyl_bessel_k_f64(0.25, z) + 3.0 * cyl_bessel_k_f64(0.75, z)
        - cyl_bessel_k_f64(1.25, z);
    (-z).exp() * (y / 2.0).powf(2.5) * b / PI.sqrt()
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
        // (unlike in xsf) like in scipy, so that a NaN term doesn't cause an infinite loop
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
        // (unlike in scipy) instead of an infinite loop
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
    let t3 = 0.5 * TAU.ln();
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

#[cfg(test)]
mod tests {
    #![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

    use num_complex::c64;

    #[test]
    fn test_stdtr() {
        const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
        const X: [f64; 9] = [-8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0];
        // values from scipy.special.stdtr
        const P_REF: [[f64; 9]; 6] = [
            [
                0.113_252_546_403_789,
                0.159_610_041_494_336,
                0.222_757_445_091_566,
                0.301_121_610_841_322,
                0.500_000_000_000_000,
                0.698_878_389_158_678,
                0.777_242_554_908_434,
                0.840_389_958_505_664,
                0.886_747_453_596_211,
            ],
            [
                0.039_583_424_160_566,
                0.077_979_130_377_369,
                0.147_583_617_650_433,
                0.250_000_000_000_000,
                0.500_000_000_000_000,
                0.750_000_000_000_000,
                0.852_416_382_349_567,
                0.922_020_869_622_631,
                0.960_416_575_839_434,
            ],
            [
                0.007_634_036_082_669,
                0.028_595_479_208_968,
                0.091_751_709_536_137,
                0.211_324_865_405_187,
                0.500_000_000_000_000,
                0.788_675_134_594_813,
                0.908_248_290_463_863,
                0.971_404_520_791_032,
                0.992_365_963_917_331,
            ],
            [
                0.002_038_288_793_893,
                0.014_004_228_005_073,
                0.069_662_984_279_422,
                0.195_501_109_477_885,
                0.500_000_000_000_000,
                0.804_498_890_522_115,
                0.930_337_015_720_578,
                0.985_995_771_994_927,
                0.997_961_711_206_107,
            ],
            [
                0.000_661_948_454_609,
                0.008_065_044_950_046,
                0.058_058_261_758_408,
                0.186_950_483_150_030,
                0.500_000_000_000_000,
                0.813_049_516_849_970,
                0.941_941_738_241_592,
                0.991_934_955_049_954,
                0.999_338_051_545_391,
            ],
            [
                0.000_005_887_471_395,
                0.001_259_166_312_368,
                0.036_694_017_385_370,
                0.170_446_566_151_030,
                0.500_000_000_000_000,
                0.829_553_433_848_970,
                0.963_305_982_614_630,
                0.998_740_833_687_632,
                0.999_994_112_528_605,
            ],
        ];
        for (i, &nu) in NU.iter().enumerate() {
            let p = X.map(|x| crate::stdtr(nu, x));
            crate::np_assert_allclose!(p, P_REF[i], atol = 1e-15);
        }
    }

    #[test]
    fn test_stdtr_extreme() {
        // reference values from mpmath
        for &(nu, x, p) in &[
            // x^2 << nu
            (7.8e5, 5.9e-8, 0.500_000_023_537_587),
            (7.8e5, -5.9e-8, 0.499_999_976_462_413),
            // large nu and 1 < x^2 <= 2
            (1e9, 1.1, 0.864_333_938_921_217_7),
            (1e8, -1.3, 0.096_800_486_083_800_24),
            // huge |x|
            (0.5, -1e200, 3.207_009_754_142_229e-101),
            (0.0016, -3.8e274, 0.180_685_979_862_586_23),
            (0.0016, 3.8e274, 0.819_314_020_137_413_8),
            (1e306, -1e200, 0.0),
            (1e306, 1e200, 1.0),
            (5e-324, f64::NEG_INFINITY, 0.0),
            (5e-324, f64::INFINITY, 1.0),
        ] {
            crate::np_assert_allclose!(&[crate::stdtr(nu, x)], &[p], rtol = 1e-13);
        }
    }

    #[test]
    fn test_stdtri() {
        const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
        const P: [f64; 7] = [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999];
        // values from Wolfram Alpha
        const X_REF: [[f64; 7]; 6] = [
            [
                -102_849.115_630_175_55,
                -1_028.491_010_471_621_9,
                -10.270_324_410_234_506,
                0.0,
                10.270_324_410_234_506,
                1_028.491_010_471_621_9,
                102_849.115_630_175_55,
            ],
            [
                -318.308_838_985_550_45,
                -31.820_515_953_773_96,
                -3.077_683_537_175_253_4,
                0.0,
                3.077_683_537_175_253_4,
                31.820_515_953_773_96,
                318.308_838_985_550_45,
            ],
            [
                -22.327_124_770_119_875,
                -6.964_556_734_283_274,
                -1.885_618_083_164_126_7,
                0.0,
                1.885_618_083_164_126_7,
                6.964_556_734_283_274,
                22.327_124_770_119_875,
            ],
            [
                -10.214_531_852_407_387,
                -4.540_702_858_568_133_6,
                -1.637_744_353_696_210_1,
                0.0,
                1.637_744_353_696_210_1,
                4.540_702_858_568_133_6,
                10.214_531_852_407_387,
            ],
            [
                -7.173_182_219_782_308_5,
                -3.746_947_387_979_197,
                -1.533_206_274_058_943_9,
                0.0,
                1.533_206_274_058_943_9,
                3.746_947_387_979_197,
                7.173_182_219_782_308_5,
            ],
            [
                -4.143_700_494_046_59,
                -2.763_769_458_112_696,
                -1.372_183_641_110_335_6,
                0.0,
                1.372_183_641_110_335_6,
                2.763_769_458_112_696,
                4.143_700_494_046_59,
            ],
        ];
        for (i, &nu) in NU.iter().enumerate() {
            let x = P.map(|p| crate::stdtri(nu, p));
            crate::np_assert_allclose!(x, X_REF[i], rtol = 1e-13);
        }
    }

    #[test]
    fn test_stdtri_huge_t() {
        // reference values from mpmath
        for &(nu, p, t) in &[
            (0.0017, 0.79, 8.561_742_739_129_646e219),
            (0.0017, 0.21, -8.561_742_739_128_98e219),
            (1.08, 5.3e-315, -3.489_628_880_599_662e290),
            (2.02, 1e-320, -1.860_555_975_007_462_2e158),
        ] {
            crate::np_assert_allclose!(&[crate::stdtri(nu, p)], &[t], rtol = 2e-13);
        }
    }

    #[test]
    fn test_ndtr_f64() {
        xsref::test("ndtr", "d-d", |x| crate::ndtr(x[0]));
    }

    #[test]
    fn test_ndtr_c64() {
        xsref::test("ndtr", "cd-cd", |x| crate::ndtr(c64(x[0], x[1])));
    }

    // Based on `scipy.special.tests.test_ndtr.TestLogNdtr.test_log_ndtr_moderate_le8`
    #[test]
    #[allow(clippy::approx_constant)]
    fn test_log_ndtr() {
        let x = [-0.75, -0.25, 0.0, 0.5, 1.5, 2.5, 3.0, 4.0, 5.0, 7.0, 8.0];
        let expected = [
            -1.484_448_229_919_656_2,
            -0.913_061_764_811_135_1,
            -0.693_147_180_559_945_3,
            -0.368_946_415_288_656_4,
            -0.069_143_455_612_233_98,
            -0.006_229_025_485_860_002,
            -0.001_350_809_964_748_193_8,
            -3.167_174_337_748_927e-05,
            -2.866_516_129_637_636e-07,
            -1.279_812_543_886_654e-12,
            -6.220_960_574_271_786e-16,
        ];
        let y = x.map(crate::log_ndtr);
        crate::np_assert_allclose!(&y, &expected, rtol = 1e-14);
    }

    #[test]
    fn test_ndtri() {
        xsref::test("ndtri", "d-d", |x| crate::ndtri(x[0]));
    }

    #[test]
    fn test_von_mises_cdf() {
        // mpmath: the integral of the density from -pi to x, plus the period count
        let ks = [1.0, 10.0, 0.1, 5.0, 49.0];
        let xs = [0.5, -1.0, 3.0, 7.0, 0.1];
        let expected = [
            0.664_076_474_604_993_1,
            0.001_410_801_706_404_327_5,
            0.979_653_423_965_240_7,
            1.935_915_236_709_023,
            0.757_381_643_056_756,
        ];
        let actual: Vec<f64> = ks
            .iter()
            .zip(xs)
            .map(|(&k, x)| crate::von_mises_cdf(k, x))
            .collect();
        xsref::np_assert_allclose!(&actual, &expected, atol = 5e-13);
        // normal approximation for k >= 50
        xsref::np_assert_allclose!(
            [crate::von_mises_cdf(100.0, 0.2)],
            [0.976_932_103_211_657],
            atol = 1e-5
        );

        assert!(crate::von_mises_cdf(-1.0, 0.5).is_nan());
        assert_eq!(crate::von_mises_cdf(1.0, f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn test_geninvgauss() {
        // (x, p, b), with the reference values from mpmath
        let args = [
            (0.01, 0.5, 1.0),
            (1.0, 0.5, 1.0),
            (5.0, 0.5, 1.0),
            (0.3, -2.5, 10.0),
            (7.0, 3.0, 0.25),
        ];
        let logpdf = [
            -47.621_353_440_210_626,
            -0.918_938_533_204_672_8,
            -3.323_657_489_421_723,
            -4.005_586_780_467_203,
            -3.924_725_429_203_324,
        ];
        let pdf = [
            2.081_176_820_202_832e-21,
            0.398_942_280_401_432_7,
            0.036_020_844_672_153_66,
            0.018_213_598_738_512_592,
            0.019_747_558_229_980_97,
        ];
        let actual = args.map(|(x, p, b)| crate::geninvgauss_logpdf(x, p, b));
        xsref::np_assert_allclose!(&actual, &logpdf, rtol = 1e-14);
        let actual = args.map(|(x, p, b)| crate::geninvgauss_pdf(x, p, b));
        xsref::np_assert_allclose!(&actual, &pdf, rtol = 1e-14);

        assert_eq!(crate::geninvgauss_logpdf(0.0, 0.5, 0.5), f64::NEG_INFINITY);
        assert_eq!(crate::geninvgauss_pdf(0.0, 0.5, 0.5), 0.0);
        assert_eq!(crate::geninvgauss_pdf(2e6, 50.0, 2.0), 0.0);
        assert!(crate::geninvgauss_logpdf(1.0, 0.5, 0.0).is_nan());
        assert!(crate::geninvgauss_logpdf(1.0, 0.5, -1.0).is_nan());
    }

    #[test]
    fn test_genhyperbolic() {
        // (x, p, a, b), with the reference values from mpmath
        let args = [
            (-3.0, 2.0, 3.0, 1.5),
            (0.0, 2.0, 3.0, 1.5),
            (1.2, 2.0, 3.0, 1.5),
            (4.0, -0.7, 0.5, -0.2),
            (-10.0, 1.0, 2.0, 0.0),
        ];
        let logpdf = [
            -11.470_885_563_926_888,
            -1.447_878_690_690_421_3,
            -0.982_289_498_087_753,
            -5.972_599_781_992_024,
            -18.825_827_120_241_21,
        ];
        let pdf = [
            1.042_936_092_158_763_5e-5,
            0.235_068_412_380_730_5,
            0.374_452_807_711_144_86,
            0.002_547_609_572_561_395,
            6.668_791_884_131_342e-9,
        ];
        let actual = args.map(|(x, p, a, b)| crate::genhyperbolic_logpdf(x, p, a, b));
        xsref::np_assert_allclose!(&actual, &logpdf, rtol = 1e-14);
        let actual = args.map(|(x, p, a, b)| crate::genhyperbolic_pdf(x, p, a, b));
        xsref::np_assert_allclose!(&actual, &pdf, rtol = 1e-14);

        // Cauchy limit: p = -1/2, a -> 0, b = 0
        let x = [-30.0, -1.0, 0.0, 0.5, 20.0];
        let actual = x.map(|x| crate::genhyperbolic_pdf(x, -0.5, f64::from(f32::EPSILON), 0.0));
        let expected = x.map(|x| 1.0 / (core::f64::consts::PI * (1.0 + x * x)));
        xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-6);
    }

    #[test]
    fn test_nrdtrimn_nrdtrisd() {
        // reference values from scipy.special
        xsref::np_assert_allclose!([crate::nrdtrimn(0.3, 2.0, 1.0)], [2.048_801_025_416_082]);
        xsref::np_assert_allclose!([crate::nrdtrisd(1.0, 0.7, 3.0)], [3.813_878_803_572_981]);
        assert!(crate::nrdtrimn(0.5, 0.0, 1.0).is_nan());
        assert!(crate::nrdtrisd(0.0, 1.0, 1.0).is_nan());
    }

    #[test]
    fn test_kolmogorov() {
        xsref::test("kolmogorov", "d-d", |x| crate::kolmogorov(x[0]));
    }

    #[test]
    fn test_kolmogc() {
        xsref::test("kolmogc", "d-d", |x| crate::kolmogc(x[0]));
    }

    #[test]
    fn test_kolmogi() {
        xsref::test("kolmogi", "d-d", |x| crate::kolmogi(x[0]));
    }

    #[test]
    fn test_kolmogci() {
        xsref::test("kolmogci", "d-d", |x| crate::kolmogci(x[0]));
    }

    #[test]
    fn test_kolmogp() {
        xsref::test("kolmogp", "d-d", |x| crate::kolmogp(x[0]));
    }

    #[test]
    fn test_smirnov() {
        xsref::test("smirnov", "p_d-d", |x| crate::smirnov(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovc() {
        xsref::test("smirnovc", "p_d-d", |x| crate::smirnovc(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovi() {
        xsref::test("smirnovi", "p_d-d", |x| crate::smirnovi(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovci() {
        xsref::test("smirnovci", "p_d-d", |x| {
            crate::smirnovci(x[0] as i32, x[1])
        });
    }

    #[test]
    fn test_smirnovp_large_n() {
        // in xsf this is -0.0, due to `int` overflow
        let (n, x, h) = (50_000, 2e-3, 1e-8);
        let expected = (crate::smirnov(n, x + h) - crate::smirnov(n, x - h)) / (2.0 * h);
        let actual = crate::smirnovp(n, x);
        assert!(
            ((actual - expected) / expected).abs() < 1e-6,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn test_smirnovp() {
        xsref::test("smirnovp", "p_d-d", |x| crate::smirnovp(x[0] as i32, x[1]));
    }

    #[test]
    fn test_owens_t() {
        xsref::test("owens_t", "d_d-d", |x| crate::owens_t(x[0], x[1]));
    }

    #[test]
    fn test_chdtr() {
        xsref::test("chdtr", "d_d-d", |x| crate::chdtr(x[0], x[1]));
    }

    #[test]
    fn test_chdtrc() {
        xsref::test("chdtrc", "d_d-d", |x| crate::chdtrc(x[0], x[1]));
    }

    #[test]
    fn test_chdtri() {
        xsref::test("chdtri", "d_d-d", |x| crate::chdtri(x[0], x[1]));
    }

    #[test]
    fn test_fdtr() {
        xsref::test("fdtr", "d_d_d-d", |x| crate::fdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtrc() {
        xsref::test("fdtrc", "d_d_d-d", |x| crate::fdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtr_fdtrc_w_near_1() {
        // reference values from mpmath
        let p = crate::fdtr(1.38, 0.112, 5.77e62);
        crate::np_assert_allclose!(&[p], &[0.999_743_687_449_994_5], rtol = 1e-13);
        let q = crate::fdtrc(0.00118, 1.46e11, 0.00286);
        crate::np_assert_allclose!(&[q], &[0.007_474_128_501_559_321], rtol = 1e-12);
        // `incbet` slightly exceeds 1 here (the true value is 4.6245e-10)
        let p = crate::fdtr(
            458_884.656_463_958_03,
            2.951_258_965_120_171e-11,
            1_072.009_775_175_712_8,
        );
        assert!((0.0..1e-9).contains(&p));
        let q = crate::fdtrc(
            2.951_258_965_120_171e-11,
            458_884.656_463_958_03,
            0.000_932_827_314_784_597_3,
        );
        assert!((0.0..1e-9).contains(&q));
    }

    #[test]
    fn test_fdtr_fdtrc_tiny_df() {
        // reference values from mpmath; one df is tiny and the other large, so that `w` is near 1
        // but below the mean
        for &(a, b, x, p) in &[
            (1e7, 0.001, 0.5, 0.003_160_968_634_795_953_3),
            (1e7, 0.002, 1.0, 0.006_312_353_390_409_258),
        ] {
            crate::np_assert_allclose!(&[crate::fdtr(a, b, x)], &[p], rtol = 1e-12);
        }
        // inaccurate either way, but the complement would be far worse
        let p = crate::fdtr(1e16, 1e-13, 6e-14);
        crate::np_assert_allclose!(&[p], &[1.462_771_840_746_950_6e-14], rtol = 0.02);
        for &(a, b, x, q) in &[
            (0.002, 1e7, 1.0, 0.006_312_353_390_409_258),
            (
                0.001_153_561_542_308_523_6,
                731_611.160_931_477_6,
                3.267_656_625_181_253,
                0.003_281_744_887_778_894_7,
            ),
        ] {
            crate::np_assert_allclose!(&[crate::fdtrc(a, b, x)], &[q], rtol = 1e-12);
        }
    }

    #[test]
    fn test_fdtr_bdtr_large_df() {
        // reference values from mpmath (see `incbet_lnpow`)
        let p = crate::fdtr(100.0, 1e20, 0.5);
        crate::np_assert_allclose!(&[p], &[6.953_305_247_616_099e-6], rtol = 1e-12);
        let p = crate::bdtr(3e8, 1_000_000_000, 0.3);
        crate::np_assert_allclose!(&[p], &[0.500_015_600_125_545], rtol = 1e-11);
    }

    #[test]
    fn test_bdtr_nbdtrc_large_n() {
        // reference values from mpmath
        let p = [
            crate::bdtr(0.0, 1_000_000_000, 1e-10),
            crate::bdtr(2.0, 2_000_000_000, 1e-9),
            crate::nbdtrc(2_000_000_000, 1, 1e-13),
            crate::nbdtrc(2_000_000_000, 3, 1e-9),
        ];
        let p_ref = [
            0.904_837_418_031_435_4,
            0.676_676_416_183_063_4,
            0.999_800_019_998_566_7,
            0.676_676_415_371_051_7,
        ];
        crate::np_assert_allclose!(&p, &p_ref, rtol = 1e-14);
    }

    #[test]
    fn test_fdtr_stdtr_huge_df() {
        // reference values from mpmath (see `betainc_xy`)
        let abx = [
            (100.0, 1e20, 0.5),
            (30.0, 1e20, 2.0),
            (1e20, 30.0, 2.0),
            (3.0, 1e14, 2.0),
        ];
        let p = [
            6.953_305_247_616_099e-6,
            0.999_079_317_603_851_4,
            0.989_739_572_087_657_4,
            0.888_389_774_905_280_1,
        ];
        let q = [
            0.999_993_046_694_752_4,
            9.206_823_961_486_663e-4,
            0.010_260_427_912_342_617,
            0.111_610_225_094_719_86,
        ];
        let rtol = 1e-13;
        crate::np_assert_allclose!(&abx.map(|(a, b, x)| crate::fdtr(a, b, x)), &p, rtol = rtol);
        crate::np_assert_allclose!(&abx.map(|(a, b, x)| crate::fdtrc(a, b, x)), &q, rtol = rtol);

        let p = [crate::stdtr(1e20, -3.0), crate::stdtr(1e20, 3.0)];
        let q = 0.001_349_898_031_630_094_6;
        crate::np_assert_allclose!(&p, &[q, 1.0 - q], rtol = 1e-14);
    }

    #[test]
    fn test_fdtri() {
        xsref::test("fdtri", "d_d_d-d", |x| crate::fdtri(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtr() {
        xsref::test("gdtr", "d_d_d-d", |x| crate::gdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtrc() {
        xsref::test("gdtrc", "d_d_d-d", |x| crate::gdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_pdtr() {
        xsref::test("pdtr", "d_d-d", |x| crate::pdtr(x[0], x[1]));
    }

    #[test]
    fn test_pdtrc() {
        xsref::test("pdtrc", "d_d-d", |x| crate::pdtrc(x[0], x[1]));
    }

    #[test]
    fn test_pdtri() {
        xsref::test("pdtri", "p_d-d", |x| crate::pdtri(x[0] as i32, x[1]));
    }

    #[test]
    fn test_bdtr() {
        xsref::test("bdtr", "d_p_d-d", |x| crate::bdtr(x[0], x[1] as i32, x[2]));
    }

    #[test]
    fn test_bdtri_small_y() {
        // mpmath: the p for which bdtr(96, 191, p) == 2.1012712563041036e-22
        let p = crate::bdtri(96.0, 191, 2.101_271_256_304_103_6e-22);
        assert!((p / 0.814_405_035_443_584 - 1.0).abs() < 1e-14, "p = {p}");
    }

    #[test]
    fn test_bdtrc() {
        xsref::test("bdtrc", "d_p_d-d", |x| {
            crate::bdtrc(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_bdtri() {
        xsref::test("bdtri", "d_p_d-d", |x| {
            crate::bdtri(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_poisson_binom_smoke() {
        // Subset of xsf's `test_poisson_binom.cpp`
        const K: [usize; 8] = [0, 7, 14, 21, 25, 26, 44, 51];
        const PMF_REF: [f64; 8] = [
            3.428_603_61e-26,
            2.575_893_01e-11,
            4.150_992_19e-5,
            4.122_831_29e-2,
            1.357_810_59e-1,
            1.357_810_59e-1,
            2.575_893_01e-11,
            3.428_603_61e-26,
        ];

        // np.linspace(1e-5, 1 - 1e-5, 51)
        let p: Vec<f64> = (0..51)
            .map(|i| 1e-5 + f64::from(i) * (1.0 - 2e-5) / 50.0)
            .collect();

        let pmf = crate::poisson_binom_pmf_all(&p);
        let cdf = crate::poisson_binom_cdf_all(&p);

        crate::np_assert_allclose!(K.map(|k| pmf[k]), PMF_REF, rtol = 1e-8);

        crate::np_assert_allclose!([cdf[25], cdf[51]], [0.5, 1.0], atol = 1e-15);
    }

    #[test]
    fn test_wilcoxon() {
        // Subset of xsf's `test_wilcoxon.cpp`
        assert_eq!(crate::wilcoxon_pmf_all(0), [1.0]);
        assert_eq!(crate::wilcoxon_pmf_all(2), [0.25; 4]);
        assert_eq!(
            crate::wilcoxon_pmf_all(5),
            [
                0.03125, 0.03125, 0.03125, 0.0625, 0.0625, 0.09375, 0.09375, 0.09375, 0.09375,
                0.09375, 0.09375, 0.0625, 0.0625, 0.03125, 0.03125, 0.03125,
            ]
        );
        assert_eq!(
            crate::wilcoxon_cdf_all(3),
            [0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 1.0]
        );

        let pmf = crate::wilcoxon_pmf_all(3);
        assert_eq!(crate::wilcoxon_pmf(3.0, &pmf), 0.25);
        assert_eq!(crate::wilcoxon_pmf(-1.0, &pmf), 0.0);
        assert_eq!(crate::wilcoxon_pmf(7.0, &pmf), 0.0);
        assert!(crate::wilcoxon_pmf(f64::NAN, &pmf).is_nan());

        let cdf = crate::wilcoxon_cdf_all(5);
        let k = [-1.0, -0.5, 5.0, 7.9, 9.0, 15.0, f64::INFINITY];
        assert_eq!(
            k.map(|k| crate::wilcoxon_cdf(k, &cdf)),
            [0.0, 0.03125, 0.3125, 0.5, 0.6875, 1.0, 1.0]
        );
        assert_eq!(
            k.map(|k| crate::wilcoxon_sf(k, &cdf)),
            [1.0, 1.0, 0.78125, 0.59375, 0.40625, 0.03125, 0.0]
        );
        assert_eq!(crate::wilcoxon_sf(f64::NEG_INFINITY, &cdf), 1.0);
        assert!(crate::wilcoxon_cdf(f64::NAN, &cdf).is_nan());
        assert!(crate::wilcoxon_sf(f64::NAN, &cdf).is_nan());
    }

    #[test]
    fn test_nbdtr() {
        xsref::test("nbdtr", "p_p_d-d", |x| {
            crate::nbdtr(x[0] as i32, x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtrc() {
        xsref::test("nbdtrc", "p_p_d-d", |x| {
            crate::nbdtrc(x[0] as i32, x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_cdf_cvm_inf_smoke() {
        assert!(crate::cdf_cvm_inf(f64::NAN).is_nan());
        assert!(crate::cdf_cvm_inf(-1.0).abs() <= f64::EPSILON);

        let result = crate::cdf_cvm_inf(0.5);
        assert!(result.is_finite());
        assert!(result > 0.0 && result < 1.0);
    }

    #[test]
    fn test_cdf_cvm_smoke() {
        const N: i32 = 10;
        let lower = 1.0 / (12.0 * f64::from(N));
        let upper = f64::from(N) / 3.0;

        assert!(crate::cdf_cvm(lower, N).abs() <= f64::EPSILON);
        assert!((crate::cdf_cvm(upper, N) - 1.0).abs() <= f64::EPSILON);

        let result = crate::cdf_cvm(0.5, N);
        assert!(result.is_finite());
        assert!((0.0..=1.0).contains(&result));
    }

    #[test]
    fn test_cdf_cvm_nonfinite() {
        assert_eq!(crate::cdf_cvm_inf(f64::INFINITY), 1.0);
        assert_eq!(crate::cdf_cvm_inf(f64::MAX), 1.0);
        assert!(crate::cdf_cvm_inf(f64::NAN).is_nan());
        for n in [-1, 5] {
            assert_eq!(crate::cdf_cvm(f64::INFINITY, n), 1.0);
            assert!(crate::cdf_cvm(f64::NAN, n).is_nan());
        }
        assert!(crate::cdf_cvm(1e4, 1_000_000).is_nan());
    }
}
