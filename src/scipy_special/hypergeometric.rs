use num_complex::Complex64;

use crate::utils::MAXLOG;
use crate::xsf::cephes::{ikv_asymptotic_uniform_sums, lgam_stirling_corr};

const LN_MAX: f64 = MAXLOG; // f64::MAX.ln()
const LN_MIN: f64 = -708.396_418_532_264_1; // f64::MIN_POSITIVE.ln()

/// Debye's expansion of $\Gamma(v) z^{(1-v)/2} I_{v-1}(2\sqrt{z})$ for real $z$ and large $|v|$
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_asy`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L59-L100>
///
/// Unlike in xsf (scipy/xsf#366), this also handles $z<0$ with $4|z| < (v-1)^2$, through Debye's
/// expansion of $J_{v-1}(2\sqrt{-z})$ (DLMF 10.19.3).
#[inline]
fn hyp0f1_asy(v: f64, z: f64) -> f64 {
    let arg = z.abs().sqrt();
    let v1 = (v - 1.0).abs();
    let x = 2.0 * arg / v1;
    let p1 = (1.0 + z.signum() * x * x).sqrt();
    let eta = p1 + x.ln() - crate::log1p(p1);

    // Unlike in xsf (scipy/xsf#366), the ~v1 ln(v1) terms of the exponent cancel analytically
    // (Stirling's formula), e.g. hyp0f1(1e12, 2e12) is 7.389 instead of 7.360:
    //
    //   zeta = ln(Gamma(v) arg^(1-v) e^(v1 eta) / sqrt(2 pi v1))                (v > 1)
    //        = -ln(Gamma(v) 2 sin(pi v1) arg^v1 e^(-v1 eta) / sqrt(2 pi v1))    (v < 1)
    //        = v1 (h - log1pmx(h)) + corr(v1),  h = (p1 - 1) / 2
    let h = z.signum() * x * x / (2.0 * (1.0 + p1));
    let zeta = v1 * (h - crate::log1pmx(h)) + lgam_stirling_corr(v1);
    let ln_sqrt_p1 = 0.5 * p1.ln();

    // Unlike in xsf (scipy/xsf#366), the large-v asymptotic correction has 11 terms (DLMF 10.41.9)
    // instead of 4 (DLMF 10.41.10)
    let (u_corr_i, u_corr_k) = ikv_asymptotic_uniform_sums(v1, 1.0 / p1);

    if v - 1.0 < 0.0 {
        // DLMF 10.27.2: I_{-v} = I_{v} + (2/pi) sin(pi*v) K_v
        // (or DLMF 10.2.3: J_{-v} = cos(pi*v) J_v - sin(pi*v) Y_v)
        let out = (-zeta - ln_sqrt_p1).exp() * u_corr_k;
        let s = 2.0 * crate::sinpi(v);
        if z < 0.0 {
            // the J_v term relative to the Y_v term, with cos(pi v1) = -cos(pi v)
            out * (1.0 - crate::cospi(v) / s * (2.0 * v1 * eta).exp() * u_corr_i / u_corr_k)
        } else {
            // with Gamma(v) = pi / (sin(pi v) Gamma(v1)), divided within the exponent
            let i = (2.0 * v1 * eta - zeta - ln_sqrt_p1 - s.abs().ln()).exp();
            out + i.copysign(s) * u_corr_i
        }
    } else {
        (zeta - ln_sqrt_p1).exp() * u_corr_i
    }
}

/// Real-valued kernel
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_real`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L31-L56>
#[inline]
fn hyp0f1_real(v: f64, z: f64) -> f64 {
    if v.is_nan() || z.is_nan() {
        f64::NAN
    } else if v <= 0.0 && v == v.floor() {
        // poles
        f64::NAN
    } else if z == 0.0 && v != 0.0 {
        // zeros
        1.0
    } else if z.abs() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)
        1.0 + z / v + z * z / (2.0 * v * (v + 1.0))
    } else if z.abs() <= v.abs().max(1.0) {
        // Unlike in xsf (scipy/xsf#366), the Maclaurin series is summed instead of the Bessel
        // functions, which lose digits for large |v| and near the poles, e.g.
        // hyp0f1(-7.000001, 0.1) is 0.985832546664 instead of 0.985832547540
        hyp0f1_series(v, z)
    } else if debye_is_accurate(v, z) {
        // Unlike in xsf (scipy/xsf#366), Debye's expansion is used instead of the Bessel
        // functions, whose prefactor loses ~eps v ln(v), e.g. hyp0f1(300, -1000) is
        // 0.035011885022349 instead of 0.035011885022344
        hyp0f1_asy(v, z)
    } else if z > 0.0 {
        let arg = z.sqrt();
        let arg_exp = crate::xlogy(1.0 - v, arg) + crate::gammaln(v);
        let bess_val = crate::xsf::cephes::iv(v - 1.0, 2.0 * arg);

        // Unlike in xsf (scipy/xsf#366), v = 1 doesn't fall back to `hyp0f1_asy`.
        if v != 1.0
            && (arg_exp > LN_MAX || bess_val == 0.0 || arg_exp < LN_MIN || bess_val.is_infinite())
        {
            // overflow or underflow
            hyp0f1_asy(v, z)
        } else {
            arg_exp.exp() * bess_val * crate::gammasgn(v)
        }
    } else {
        let arg = (-z).sqrt();
        let arg_pow = arg.powf(1.0 - v);
        let gamma_v = crate::gamma(v);
        let prefactor = arg_pow * gamma_v;
        let bess_val = crate::bessel_j(v - 1.0, 2.0 * arg);
        let out = prefactor * bess_val;
        let prefactor_ok = arg_pow.is_normal() && gamma_v.is_normal() && prefactor.is_normal();

        // Unlike in xsf (scipy/xsf#366), fall back if a factor over- or underflows
        if prefactor_ok && out.is_normal() {
            out
        } else if bess_val.is_normal() {
            (crate::xlogy(1.0 - v, arg) + crate::gammaln(v)).exp() * bess_val * crate::gammasgn(v)
        } else if 4.0 * -z < (v - 1.0).powi(2) && (v > 1.0 || !bess_val.is_finite()) {
            // Debye's expansion if J_{v-1}(2 sqrt(-z)) under- (or for v < 1, over-) flowed in its
            // non-oscillatory regime
            hyp0f1_asy(v, z)
        } else {
            out
        }
    }
}

/// Whether Debye's expansion in `hyp0f1_asy` is accurate (empirically to within ~50 eps)
///
/// As in cephes `iv`, this is for $|v-1| > 50$. For $z<0$ its terms grow like $p^{-3k}$, with
/// $p = \sqrt{1 - 4|z|/(v-1)^2}$, so then this is for $|v-1| p^3 > 50$.
fn debye_is_accurate(v: f64, z: f64) -> bool {
    let v1 = (v - 1.0).abs();
    let p2 = 1.0 + 4.0 * z.min(0.0) / (v1 * v1);
    v1 * p2 * p2.sqrt() > 50.0
}

/// Maclaurin series $\sum_n z^n / ((v)_n n!)$ (DLMF 16.2.1)
///
/// Only used for $|z| \le \max(|v|, 1)$, where the terms decrease rapidly, except near the poles
/// at $n \approx -v$ if $v < 0$.
fn hyp0f1_series(v: f64, z: f64) -> f64 {
    let mut term = 1.0;
    let mut sum = 1.0;
    for n in 0..1000 {
        let n = f64::from(n);
        term *= z / (v + n) / (n + 1.0);
        sum += term;
        // not before the term with the smallest |v + n|, e.g. 1 / (-7.000001)_8 ~ 1e6 / 7!
        if term.abs() <= f64::EPSILON * sum.abs() && (v + n > -0.5 || term == 0.0) {
            break;
        }
    }
    sum
}

/// Complex valued kernel
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_cmplx`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L106-L136>
#[inline]
fn hyp0f1_cmplx(v: f64, z: Complex64) -> Complex64 {
    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        Complex64::new(f64::NAN, f64::NAN)
    } else if v <= 0.0 && v == v.floor() {
        // poles
        f64::NAN.into()
    } else if z.re == 0.0 && z.im == 0.0 && v != 0.0 {
        // zeros
        1.0.into()
    } else if z.norm() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)

        // need to do computations in this order, for otherwise $v\approx -z \ll 1$
        // it can lose precision (as was reported for 32-bit linux, see scipy/scipy#6365)
        let t1 = 1.0 + z / v;
        let t2 = z * z / (2.0 * v * (v + 1.0));
        t1 + t2
    } else {
        let (arg, r) = if z.re > 0.0 {
            let arg = z.sqrt();
            let bess_val = crate::bessel_i(v - 1.0, 2.0 * arg);
            (arg, bess_val)
        } else {
            let arg = (-z).sqrt();
            let bess_val = crate::bessel_j(v - 1.0, 2.0 * arg);
            (arg, bess_val)
        };

        r * crate::gamma(v) * arg.powf(1.0 - v)
    }
}

pub trait HypergeometricArg: crate::sealed::Sealed {
    fn hyp0f0(self) -> Self;
    fn hyp1f0(self, a: f64) -> Self;
    fn hyp0f1(self, b: f64) -> Self;
}

impl HypergeometricArg for f64 {
    #[inline]
    fn hyp0f0(self) -> Self {
        self.exp()
    }

    #[inline]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self == 1.0 {
            f64::NAN
        } else {
            (1.0 - self).powf(-a)
        }
    }

    #[inline]
    fn hyp0f1(self, b: f64) -> Self {
        hyp0f1_real(b, self)
    }
}

impl HypergeometricArg for num_complex::Complex<f64> {
    #[inline]
    fn hyp0f0(self) -> Self {
        self.exp()
    }

    #[inline]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self.re == 1.0 && self.im == 0.0 {
            f64::NAN.into()
        } else {
            (1.0 - self).powf(-a)
        }
    }

    #[inline]
    fn hyp0f1(self, b: f64) -> Self {
        hyp0f1_cmplx(b, self)
    }
}

/// Hypergeometric function $_0F_0\left[\middle\| z\right]$ for real or complex $z$
///
/// This function is currently not implemented in SciPy, but is straightforward to evaluate.
/// Both [`f64`] and [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _0F_0\left[ \middle\| z\right] = \sum\_{n=0}^\infty {z^n \over n!},
/// $$
///
/// which is just the exponential function, $e^z$.
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function $_1F_0\left[a\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
#[inline]
pub fn hyp0f0<T: HypergeometricArg>(z: T) -> T {
    z.hyp0f0()
}

/// Hypergeometric function $_1F_0\left[a\middle\| z\right]$ for real or complex $z$
///
/// This function is currently not implemented in SciPy, but is straightforward to evaluate.
/// Both [`f64`] and [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _1F_0\left[a\middle\| z\right] = \sum\_{n=0}^\infty \rpow a n {z^n \over n!}
/// $$
///
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
/// It's also equal to $(1 - z)^{-a}$.
///
/// # See also
/// - [`hyp0f0`]: Hypergeometric function $_0F_0\left[\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
#[inline]
pub fn hyp1f0<T: HypergeometricArg>(a: f64, z: T) -> T {
    z.hyp1f0(a)
}

/// Confluent hypergeometric limit function $_0F_1\left[b\middle\| z\right]$ for real or complex $z$
///
/// This is a translation of the [`scipy.special.hyp0f1`][hyp0f1] Cython implementation into Rust,
/// as in [scipy/xsf][xsf], which uses AMOS (instead of Cephes) for $J_v$.
///
/// [hyp0f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp0f1.html
/// [xsf]: https://github.com/scipy/xsf
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _0F_1\left[b\middle\| z\right] = \sum\_{n=0}^\infty {1\over \rpow b n} {z^n \over n!}.
/// $$
///
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
///
/// It's also the limit as $a \to \infty$ of $\hyp 1 1 a b {z/a}$, and satisfies the
/// differential equation $z f\'\'(z) + b f\'(z) = f(z)$. See [^1] for more information.
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function, $_1F_0\left[a\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
/// [^1]: Weisstein, Eric W. "Confluent Hypergeometric Limit Function." From MathWorld -- A Wolfram
/// Resource. <https://mathworld.wolfram.com/ConfluentHypergeometricLimitFunction.html>
///
#[inline]
pub fn hyp0f1<T: HypergeometricArg>(b: f64, z: T) -> T {
    z.hyp0f1(b)
}
