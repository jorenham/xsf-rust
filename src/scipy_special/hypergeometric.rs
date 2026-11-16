//! `hyp0f1` is translated into pure Rust from `xsf/cpu/hyp0f1.h` (scipy/xsf@fa36f89).

use crate::utils::is_nonpos_int;
use crate::xsf::bessel::{cyl_bessel_i_c64, cyl_bessel_j_c64};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::consts::MAXLOG;
use crate::xsf::cephes::dd_real::{two_prod, two_sum};
use crate::xsf::cephes::{ikv_asymptotic_uniform_sums, lgam_stirling_corr};
use crate::xsf::complex::Num;
use crate::xsf::log::clog1p;

const LN_MIN: f64 = -708.396_418_532_264_1; // f64::MIN_POSITIVE.ln()

/// Debye's expansion of `Gamma(v) z^((1-v)/2) I_{v-1}(2 sqrt(z))` for real `z` and large `|v|`
///
/// Unlike in xsf (scipy/xsf#366), this also handles `z<0` with `4|z| < (v-1)^2`, through Debye's
/// expansion of `J_{v-1}(2 sqrt(-z))` (DLMF 10.19.3).
#[inline]
fn hyp0f1_asy(v: f64, z: f64) -> f64 {
    let arg = z.abs().sqrt();
    let v1 = (v - 1.0).abs();
    let x = 2.0 * arg / v1;
    let p1 = (1.0 + z.signum() * x * x).sqrt();
    let eta = p1 + x.ln() - crate::log1p(p1);

    // Unlike in xsf (scipy/xsf#366), the ~v1 ln(v1) terms of the exponent cancel analytically
    // (Stirling's formula):
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

/// `hyp0f1_asy` for complex `z`, if `I_{v-1}(2 sqrt(z))` over- or underflows
///
/// As in xsf, it is only used for `v > 1`, so without the `K` term of DLMF 10.27.2.
///
/// Unlike in xsf (scipy/xsf#557), this is `None` unless `|e^(v1 eta)|` is beyond `e^(+-MAXLOG/2)`,
/// so also near `(-inf, -(v-1)^2/4]`, where `|e^(v1 eta)| = 1` and the expansion is invalid.
fn hyp0f1_asy_c64(v: f64, z: C64) -> Option<C64> {
    let arg = z.sqrt();
    let v1 = v - 1.0;
    let x = 2.0 * arg / v1;
    let p1 = (1.0 + x * x).sqrt();

    // `v1 eta = v1 (1 + ln(arg / v1) + d)`, with `d = 2 h - ln(1 + h)` and `h = (p1 - 1) / 2`
    let h = x * x / (2.0 * (1.0 + p1));
    let v1_d = v1 * (2.0 * h - clog1p(h));

    // `|I_{v-1}(2 arg)|` is `|e^(v1 eta)|` up to a factor `(2 pi v1 |p1|)^(-1/2)`, for which
    // `MAXLOG / 2` leaves room
    if (v1_d.re + v1 * (1.0 + (arg.norm() / v1).ln())).abs() < MAXLOG / 2.0 {
        return None;
    }

    // Unlike in xsf (scipy/xsf#557), the ~v1 ln(v1) terms of the exponent cancel analytically, as
    // in `hyp0f1_asy`
    let zeta = v1_d + lgam_stirling_corr(v1) - 0.5 * p1.ln();

    // large-v asymptotic correction, DLMF 10.41.10
    let pp = 1.0 / p1;
    let p2 = pp * pp;
    let p4 = p2 * p2;
    let p6 = p4 * p2;
    let u1 = (3.0 - 5.0 * p2) * pp / 24.0;
    let u2 = (81.0 - 462.0 * p2 + 385.0 * p4) * p2 / 1152.0;
    let u3 = (30375.0 - 369_603.0 * p2 + 765_765.0 * p4 - 425_425.0 * p6) * pp * p2 / 414_720.0;
    let u_corr_i = 1.0 + u1 / v1 + u2 / (v1 * v1) + u3 / (v1 * v1 * v1);

    Some(zeta.exp() * u_corr_i)
}

/// `xsf::hyp0f1` for real `z`
#[inline]
fn hyp0f1_real(v: f64, z: f64) -> f64 {
    if v.is_nan() || z.is_nan() {
        f64::NAN
    } else if is_nonpos_int(v) {
        // poles
        f64::NAN
    } else if z == 0.0 && v != 0.0 {
        // zeros
        1.0
    } else if z.abs() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)
        1.0 + z / v + z * z / (2.0 * v * (v + 1.0))
    } else if z.abs() <= v.abs().max(1.0) {
        // Unlike in xsf (scipy/xsf#366), the Maclaurin series is summed here instead of using the
        // Bessel functions, which lose digits near the poles and for large |v|.
        hyp0f1_series(v, z)
    } else if debye_is_accurate(v, z) {
        // Unlike in xsf (scipy/xsf#366), Debye's expansion is used wherever it is accurate, as the
        // prefactor of the Bessel functions loses ~eps v ln(v) for large |v|.
        hyp0f1_asy(v, z)
    } else if z > 0.0 {
        let arg = z.sqrt();
        let arg_exp = crate::xlogy(1.0 - v, arg) + crate::gammaln(v);
        let bess_val = crate::xsf::cephes::iv(v - 1.0, 2.0 * arg);

        // Unlike in xsf (scipy/xsf#366), v = 1 doesn't fall back to `hyp0f1_asy`, which divides
        // by |v - 1|.
        if v != 1.0
            && (arg_exp > MAXLOG || bess_val == 0.0 || arg_exp < LN_MIN || bess_val.is_infinite())
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
        let mut bess_val = crate::bessel_j(v - 1.0, 2.0 * arg);
        // Unlike in xsf, the rounding error `d` of `arg` is corrected for, as it shifts the phase
        // of `J` by `2 d`, up to `arg eps`. By DLMF 10.6.6 and 10.17.3, `a^(1-v) J_{v-1}(2 a)` at
        // `a = arg + d` is about `arg^(1-v) (J_{v-1}(x) cos 2d - J_v(x) sin 2d)` for `x = 2 arg`,
        // exactly so to first order in `d`.
        let (sq, sq_err) = two_prod(arg, arg);
        let d = (-z - sq - sq_err) / (2.0 * arg);
        let (sin, cos) = (2.0 * d).sin_cos();
        let shift = -sin * crate::bessel_j(v, 2.0 * arg);
        if shift.is_finite() {
            bess_val = bess_val * cos + shift;
        }
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
/// As in cephes `iv`, this is for `|v-1| > 50`. For `z<0` its terms grow like `p^(-3k)`, with
/// `p = sqrt(1 - 4|z|/(v-1)^2)`, so then this is for `|v-1| p^3 > 50`.
fn debye_is_accurate(v: f64, z: f64) -> bool {
    let v1 = (v - 1.0).abs();
    let p2 = 1.0 + 4.0 * z.min(0.0) / (v1 * v1);
    v1 * p2 * p2.sqrt() > 50.0
}

/// Maclaurin series `sum_n z^n / ((v)_n n!)` (DLMF 16.2.1)
///
/// Only used for `|z| <= max(|v|, 1)`, where the terms decrease rapidly, except near the poles
/// at `n ~ -v` if `v < 0`.
fn hyp0f1_series<T: Num>(v: f64, z: T) -> T {
    let mut term = T::from(1.0);
    let mut sum = term;
    for n in 0..1000 {
        let n = f64::from(n);
        term = term * (z / (v + n) / (n + 1.0));
        sum = sum + term;
        // not before the term with the smallest |v + n|, e.g. 1 / (-7.000001)_8 ~ 1e6 / 7!
        // `sum` is scaled first, as its complex `abs` can overflow
        if term.abs() <= (sum * f64::EPSILON).abs() && (v + n > -0.5 || term.abs() == 0.0) {
            break;
        }
    }
    sum
}

/// `xsf::hyp0f1` for complex `z`
#[inline]
fn hyp0f1_cmplx(v: f64, z: C64) -> C64 {
    if v.is_nan() || z.is_nan() {
        C64::NAN
    } else if is_nonpos_int(v) {
        // poles
        f64::NAN.into()
    } else if z == C64::ZERO && v != 0.0 {
        // zeros
        C64::ONE
    } else if z.norm() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)

        // need to do computations in this order, for otherwise v ~ -z << 1
        // it can lose precision (as was reported for 32-bit linux, see scipy/scipy#6365)
        let t1 = 1.0 + z / v;
        let t2 = z * z / (2.0 * v * (v + 1.0));
        t1 + t2
    } else {
        // Unlike in xsf, the Maclaurin series is summed here like for real `z`, instead of using
        // the Bessel functions, which lose digits for large `|v|`, e.g. ~eps ln(Gamma(v)) through
        // the logarithmic prefactor; they remain for where the partial sums overflow (tiny `|v|`)
        if z.norm() <= v.abs().max(1.0) {
            let sum = hyp0f1_series(v, z);
            if sum.is_finite() {
                return sum;
            }
        }
        // `arg = sqrt(w)`, with `w = z` for `I` and `w = -z` for `J`
        let (w, s, bessel): (_, _, fn(f64, C64) -> C64) = if z.re > 0.0 {
            (z, 1.0, cyl_bessel_i_c64)
        } else {
            (-z, -1.0, cyl_bessel_j_c64)
        };
        let arg = w.sqrt();
        let mut r = bessel(v - 1.0, 2.0 * arg);

        // Unlike in xsf, the rounding error `d` of `arg` is corrected for (as for real `z`): it
        // shifts the argument of the Bessel function by `2 d`, up to `|arg| eps`. By DLMF 10.29.4
        // (10.6.6), `a^(1-v) I_{v-1}(2 a)` has the derivative `2 a^(1-v) I_v(2 a)` (with `J`:
        // `-2 a^(1-v) J_v(2 a)`), and `d = (w - arg^2) / (2 arg)`, with the residual exact up to
        // `O(eps^2 |w|)`.
        let (re2, re2_err) = two_prod(arg.re, arg.re);
        let (im2, im2_err) = two_prod(arg.im, arg.im);
        let (diff, diff_err) = two_sum(re2, -im2);
        let (cross, cross_err) = two_prod(2.0 * arg.re, arg.im);
        let res_re = (w.re - diff) - (diff_err + re2_err - im2_err);
        let res = c64(res_re, (w.im - cross) - cross_err);
        if res != C64::ZERO {
            let shift = s * res / arg * bessel(v, 2.0 * arg);
            if shift.is_finite() {
                r += shift;
            }
        }

        let (gamma_v, arg_pow) = (crate::gamma(v), arg.powf(1.0 - v));
        let out = r * gamma_v * arg_pow;
        let normal = gamma_v.is_normal() && arg_pow.norm().is_normal() && out.norm().is_normal();
        if !r.is_finite() || r == C64::ZERO {
            // the Bessel function over- or underflows, or rounds to 0 near a zero
            let asy = if v > 1.0 { hyp0f1_asy_c64(v, z) } else { None };
            asy.unwrap_or(out)
        } else if normal {
            out
        } else {
            // Gamma(v) or the power over- or underflows
            //
            // Unlike in xsf (scipy/xsf#557), this is also done if a factor or the product is
            // subnormal, and the phase is applied first, so that `e^Re(ln)` only overflows if the
            // result does.
            let m = r.re.abs().max(r.im.abs());
            let ln = arg.ln() * (1.0 - v) + crate::gammaln(v) + m.ln();
            (r / m * C64::cis(ln.im) * crate::gammasgn(v)).mul_exp(ln.re)
        }
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
        self.via(C64::exp)
    }

    #[inline]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self.re == 1.0 && self.im == 0.0 {
            f64::NAN.into()
        } else {
            self.via(|z| (1.0 - z).powf(-a))
        }
    }

    #[inline]
    fn hyp0f1(self, b: f64) -> Self {
        self.via(|z| hyp0f1_cmplx(b, z))
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
