//! `hyp0f1` is translated into pure Rust from `xsf/cpu/hyp0f1.h` (scipy/xsf@fa36f89).

use core::f64::consts::PI;

use crate::ext::ln_dd;
use crate::utils::is_nonpos_int;
use crate::xsf::bessel::{cyl_bessel_i_c64, cyl_bessel_j_c64};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::consts::MAXLOG;
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, two_prod, two_sum};
use crate::xsf::cephes::{ikv_asymptotic_uniform_sums, iv, lgam_stirling_corr};
use crate::xsf::complex::Num;
use crate::xsf::log::clog1p;

const LN_MIN: f64 = -708.396_418_532_264_1; // f64::MIN_POSITIVE.ln()

/// The `|v-1|` above which Debye's expansion is used, as in cephes `iv`
const DEBYE_MIN_V1: f64 = 50.0;

/// Debye's expansion of `Gamma(v) z^((1-v)/2) I_{v-1}(2 sqrt(z))` for real `z` and large `|v|`
///
/// Unlike in xsf (scipy/xsf#366), this also handles `z<0` with `4|z| < (v-1)^2`, through Debye's
/// expansion of `J_{v-1}(2 sqrt(-z))` (DLMF 10.19.3).
#[inline]
fn hyp0f1_asy(v: f64, z: f64) -> f64 {
    // Unlike in xsf, the exponents are in double-double, as they are about `ln(F)`, which amplifies
    // their rounding errors, e.g. of `x = 2 arg / v1`, with `arg = sqrt(|z|)`, or of `v1 = |v - 1|`
    let v1 = (Dd::new(v) - 1.0).abs();
    // `z / v1`, with a huge `z` halved first (exact, as is the doubling of the normal quotient), as
    // the division forms `v1 (z / v1)`, which may overflow near `f64::MAX`
    let z_v1 = if z.abs() > f64::MAX / 16.0 {
        Dd::new(0.5 * z) / v1 * 2.0
    } else {
        Dd::new(z) / v1
    };
    let x2 = z_v1 / v1 * 4.0; // sign(z) x^2
    let p1 = (x2 + 1.0).sqrt();

    // Unlike in xsf (scipy/xsf#366), the ~v1 ln(v1) terms of the exponent cancel analytically
    // (Stirling's formula):
    //
    //   zeta = ln(Gamma(v) arg^(1-v) e^(v1 eta) / sqrt(2 pi v1))                (v > 1)
    //        = -ln(Gamma(v) 2 sin(pi v1) arg^v1 e^(-v1 eta) / sqrt(2 pi v1))    (v < 1)
    //        = v1 d + corr(v1),  d = 2 h - ln(1 + h),  h = (p1 - 1) / 2
    let h = x2 / ((p1 + 1.0) * 2.0);
    let v1_d = (h * 2.0 - ln_dd(h + 1.0)) * v1;
    let zeta = v1_d + lgam_stirling_corr(v1.hi);
    // (small and not multiplied by `v1`, so `f64` suffices, up to the low word of `p1`)
    let ln_sqrt_p1 = Dd::new(0.5 * p1.hi.ln()) + 0.5 * p1.lo / p1.hi;

    // Unlike in xsf (scipy/xsf#366), the large-v asymptotic correction has 11 terms (DLMF 10.41.9)
    // instead of 4 (DLMF 10.41.10)
    let (u_corr_i, u_corr_k) = ikv_asymptotic_uniform_sums(v1.hi, 1.0 / p1.hi);

    if v - 1.0 < 0.0 {
        // DLMF 10.27.2: I_{-v} = I_{v} + (2/pi) sin(pi*v) K_v
        // (or DLMF 10.2.3: J_{-v} = cos(pi*v) J_v - sin(pi*v) Y_v)
        let out = exp_dd(-zeta - ln_sqrt_p1) * u_corr_k;
        let s = 2.0 * crate::sinpi(v);
        // `v1 eta = v1 (1 + ln(x / 2)) + v1 d`
        let v1_eta = (ln_dd(x2.abs() / 4.0) * 0.5 + 1.0) * v1 + v1_d;
        if z < 0.0 {
            // the J_v term relative to the Y_v term, with cos(pi v1) = -cos(pi v)
            out * (1.0 - crate::cospi(v) / s * exp_dd(v1_eta * 2.0) * u_corr_i / u_corr_k)
        } else {
            // with Gamma(v) = pi / (sin(pi v) Gamma(v1)), divided within the exponent where `e^e`
            // overflows
            let e = v1_eta * 2.0 - zeta - ln_sqrt_p1;
            let i = if e.hi < MAXLOG {
                exp_dd(e) / s.abs()
            } else {
                exp_dd(e - s.abs().ln())
            };
            out + i.copysign(s) * u_corr_i
        }
    } else {
        exp_dd(zeta - ln_sqrt_p1) * u_corr_i
    }
}

/// `e^x`, with the low word of `x` to first order (it is NaN where `x` overflowed)
fn exp_dd(x: Dd) -> f64 {
    let e = x.hi.exp();
    if e.is_normal() { e + e * x.lo } else { e }
}

/// `hyp0f1_asy` for complex `z`, if the Bessel function over- or underflows
///
/// Unlike in xsf, this is also used for `v < 1`, with the `K` term of DLMF 10.27.2 (as in
/// `hyp0f1_asy`), which overflows where `Gamma(v)` underflows.
///
/// Unlike in xsf (scipy/xsf#557), this is `None` unless `|e^(v1 eta)|` is beyond `e^(+-MAXLOG/2)`,
/// so also near `(-inf, -(v-1)^2/4]`, where `|e^(v1 eta)| = 1` and the expansion is invalid.
fn hyp0f1_asy_c64(v: f64, z: C64) -> Option<C64> {
    let arg = z.sqrt();
    let v1 = (v - 1.0).abs();
    let x = 2.0 * arg / v1;
    let p1 = (1.0 + x * x).sqrt();

    // `v1 eta = v1 (1 + ln(arg / v1) + d)`, with `d = 2 h - ln(1 + h)` and `h = (p1 - 1) / 2`
    let h = x * x / (2.0 * (1.0 + p1));
    let v1_d = v1 * (2.0 * h - clog1p(h));

    // `|I_{v1}(2 arg)|` and `|K_{v1}(2 arg)|` are `|e^(+-v1 eta)|` up to a factor of about
    // `v1^(-1/2)`, for which `MAXLOG / 2` leaves room
    if (v1_d.re + v1 * (1.0 + (arg.norm() / v1).ln())).abs() < MAXLOG / 2.0 {
        return None;
    }

    // Unlike in xsf (scipy/xsf#557), the ~v1 ln(v1) terms of the exponent cancel analytically, as
    // in `hyp0f1_asy`
    let corr = lgam_stirling_corr(v1);
    let zeta = v1_d + corr;
    let ln_sqrt_p1 = 0.5 * p1.ln();

    // Unlike in xsf, the large-v asymptotic correction has 11 terms instead of 4, as in `hyp0f1_asy`
    let (u_corr_i, u_corr_k) = ikv_asymptotic_uniform_sums(v1, 1.0 / p1);

    if v > 1.0 {
        return Some((zeta - ln_sqrt_p1).exp() * u_corr_i);
    }
    // DLMF 10.27.2 as in `hyp0f1_asy`, where `2 v1 eta - zeta` is the exponent of the `I` term
    let s = 2.0 * crate::sinpi(v);
    let k = (-zeta - ln_sqrt_p1).exp() * u_corr_k;
    let i_exp = 2.0 * v1 * (1.0 + (arg / v1).ln()) + v1_d - corr;
    let i = (i_exp - ln_sqrt_p1 - s.abs().ln()).exp() * u_corr_i;
    Some(k + i * s.signum())
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
    } else if z.abs() <= v.abs().max(1.0) {
        // Unlike in xsf (scipy/xsf#366), the Maclaurin series is summed here instead of using the
        // Bessel functions, which lose digits near the poles and for large |v|. Also for small
        // `z`, as a near-zero `v + n` can amplify later terms.
        hyp0f1_series(v, z)
    } else if debye_is_accurate(v, z) {
        // Unlike in xsf (scipy/xsf#366), Debye's expansion is used wherever it is accurate, as the
        // prefactor of the Bessel functions loses ~eps v ln(v) for large |v|.
        hyp0f1_asy(v, z)
    } else if z > 0.0 {
        let arg = z.sqrt();
        let arg_exp = crate::xlogy(1.0 - v, arg) + crate::gammaln(v);
        let x = 2.0 * arg;
        let bess_val = bessel_below(v, x, 1.0, || iv(v, x), |n| iv(n, x));

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
        let x = 2.0 * arg;
        // AMOS, as in xsf (SciPy uses Cephes `jv`)
        let j_v = crate::bessel_j(v, x);
        let mut bess_val = bessel_below(v, x, -1.0, || j_v, |n| crate::bessel_j(n, x));
        // Unlike in xsf, the rounding error `d` of `arg` is corrected for, as it shifts the phase
        // of `J` by `2 d`, up to `arg eps`. By DLMF 10.6.6 and 10.17.3, `a^(1-v) J_{v-1}(2 a)` at
        // `a = arg + d` is about `arg^(1-v) (J_{v-1}(x) cos 2d - J_v(x) sin 2d)` for `x = 2 arg`,
        // exactly so to first order in `d`.
        let (sq, sq_err) = two_prod(arg, arg);
        let d = (-z - sq - sq_err) / (2.0 * arg);
        let (sin, cos) = (2.0 * d).sin_cos();
        let shift = -sin * j_v;
        if shift.is_finite() {
            bess_val = bess_val * cos + shift;
        }
        let (out, prefactor_ok) = if (v - 1.0).abs() > DEBYE_MIN_V1 {
            (gamma_pow(v, arg, bess_val), true)
        } else {
            let (arg_pow, gamma_v) = (arg.powf(1.0 - v), crate::gamma(v));
            let prefactor = arg_pow * gamma_v;
            let normal = arg_pow.is_normal() && gamma_v.is_normal() && prefactor.is_normal();
            (prefactor * bess_val, normal)
        };

        // Unlike in xsf (scipy/xsf#366), fall back if a factor over- or underflows, with `J` in
        // the exponent, as the prefactor alone can overflow where the result doesn't
        if prefactor_ok && out.is_normal() {
            out
        } else if bess_val.is_normal() {
            let ln = crate::xlogy(1.0 - v, arg) + crate::gammaln(v) + bess_val.abs().ln();
            ln.exp().copysign(bess_val) * crate::gammasgn(v)
        } else if 4.0 * -z < (v - 1.0).powi(2) && (v > 1.0 || !bess_val.is_finite()) {
            // Debye's expansion if J_{v-1}(2 sqrt(-z)) under- (or for v < 1, over-) flowed in its
            // non-oscillatory regime
            hyp0f1_asy(v, z)
        } else {
            out
        }
    }
}

/// `Gamma(v) arg^(1-v) j`, for `|v-1| > DEBYE_MIN_V1`
///
/// Unlike in xsf, the ~v1 ln(v1) terms of its exponent cancel analytically (Stirling's formula, as
/// in `hyp0f1_asy`), as the power `arg^(1-v)` loses ~eps v1 ln(arg).
fn gamma_pow(v: f64, arg: f64, j: f64) -> f64 {
    // `ln(Gamma(v1 + 1) arg^-v1) = v1 (ln(v1 / arg) - 1) + ln(2 pi v1) / 2 + corr(v1)`
    let v1 = (Dd::new(v) - 1.0).abs();
    let e = (ln_dd(v1 / arg) - 1.0) * v1 + lgam_stirling_corr(v1.hi);
    let c = (2.0 * PI * v1.hi).sqrt();
    // DLMF 5.5.3: `Gamma(v) = pi v1 / (sin(pi v) Gamma(v1 + 1))` for `v < 1`
    let (e, f) = if v > 1.0 {
        (e, c)
    } else {
        (-e, PI * v1.hi / (crate::sinpi(v) * c))
    };
    let exp = exp_dd(e);
    let out = exp * f * j;
    if (exp.is_normal() && out.is_normal()) || !j.is_normal() {
        out
    } else {
        // `e^e` or a product over- or underflows, but maybe not the result
        let m = f * j;
        exp_dd(e + m.abs().ln()).copysign(m)
    }
}

/// Whether Debye's expansion in `hyp0f1_asy` is accurate, i.e. more so than the Bessel functions
///
/// This is for `|v-1| > DEBYE_MIN_V1`. For `z<0` its terms grow like `(|v-1| p^3)^-k`,
/// with `p = sqrt(1 - 4|z|/(v-1)^2)`, so then it is also for `|v-1| p^3 > 9 (1 + ln|v-1|)`, where
/// the error of its 11 terms (measured: up to ~1000 eps) is about that of the Bessel functions.
fn debye_is_accurate(v: f64, z: f64) -> bool {
    let v1 = (v - 1.0).abs();
    let p2 = 1.0 + 4.0 * z.min(0.0) / (v1 * v1);
    v1 > DEBYE_MIN_V1 && v1 * p2 * p2.sqrt() > 9.0 * (1.0 + v1.ln())
}

/// `B_{v-1}(x)` of the Bessel function `B = I` (with `s = 1`) or `J` (`s = -1`), given `B_v(x)`
///
/// Unlike in xsf, if `v - 1` rounds, this uses the recurrence (DLMF 10.29.1, 10.6.1) from the
/// orders `v` and `v + 1` where its estimated error is below a quarter of that of `b(v - 1)`,
/// which `sin(pi v)` amplifies near the poles of `Gamma(v)` (DLMF 10.27.2, 10.2.3).
fn bessel_below<T: Num>(v: f64, x: T, s: f64, b_v: impl FnOnce() -> T, b: impl Fn(f64) -> T) -> T {
    let mu = 1.0 - v;
    let d = two_sum(v, -1.0).1.abs();
    // The error of `b(v - 1)`: the rounding error `d` of the order, times `pi |cot(pi v)|`, times
    // the share of the `sin(pi v)` term, about `2 |sin(pi v)| e^(2 mu Re g)` (DLMF 10.41.3-4,
    // 10.19.3). The orders `v`, `v + 1` are exact for `v <= -1/2`, and below the turning point
    // `|x| = mu` (past it, the error of `B_v` and `B_{v+1}` grows with `|x|`).
    let amp = if v <= -0.5 && d > 0.0 && x.abs() < mu {
        let t = c64(x.real(), x.imag()) / mu;
        let w = (1.0 + s * t * t).sqrt();
        let g = ((1.0 + w) / t).norm().ln() - w.re;
        let sin = crate::sinpi(v).abs();
        let share = (2.0 * sin * (2.0 * mu * g).exp()).min(1.0);
        PI * d * crate::cospi(v).abs() / sin * share
    } else {
        0.0
    };
    let tol = v.abs() * f64::EPSILON;
    if amp > tol {
        let (t1, t2) = (T::from(2.0 * v) / x * b_v(), b(v + 1.0) * s);
        let r = t1 + t2;
        let kappa = (t1.abs() + t2.abs()) / r.abs();
        // The (also rough) error of `r`: its cancellation `kappa` times the error of `B_v` and
        // `B_{v+1}` (empirically up to `|v| eps`), plus the `d |ln(x/2)|` that the rounded
        // exponent `1 - v` of `(x/2)^(1-v)` costs. False if `r` overflows or underflows to zero.
        if 4.0 * (kappa * tol + d * (x.abs() / 2.0).ln().abs()) < amp {
            return r;
        }
    }
    b(v - 1.0)
}

/// Maclaurin series `sum_n z^n / ((v)_n n!)` (DLMF 16.2.1)
///
/// Only used for `|z| <= max(|v|, 1)`, where the terms decrease rapidly, except near the poles
/// at `n ~ -v` if `v < 0`.
fn hyp0f1_series<T: Num>(v: f64, z: T) -> T {
    // a tiny `v + n` can overflow the terms of a finite result: then the sum is retried times
    // `2^-SHIFT`, but only then, as that underflows tiny results, e.g. for `v = -z << 1`
    const SHIFT: i32 = 600;
    let sum = hyp0f1_series_scaled(v, z, 0);
    if sum.is_finite() {
        sum
    } else {
        hyp0f1_series_scaled(v, z, -SHIFT).ldexp(SHIFT)
    }
}

/// `2^shift` times the Maclaurin series
fn hyp0f1_series_scaled<T: Num>(v: f64, z: T, shift: i32) -> T {
    let one = T::from(1.0).ldexp(shift);
    let mut term = if shift == 0 {
        one * (z / v)
    } else {
        one / v * z
    };
    // the first two terms as `(v + z) / v`, which is exact where they cancel, e.g. `v ~ -z`, but
    // not if `v + z` overflows, nor on the scaled retry, where the scaled `z.im` can underflow
    let mut sum = match T::from(v) + z {
        vz if shift == 0 && vz.is_finite() => vz / v,
        _ => one + term,
    };
    // the larger part is cheaper than the complex `abs`, and as good for stopping
    let size = |x: T| x.real().abs().max(x.imag().abs());
    for n in 1..1000 {
        let n = f64::from(n);
        // `1 / (v + n)` first only where `z / (v + n)` would overflow, as `term` is scaled down
        term = if shift == 0 {
            term * (z / (v + n) / (n + 1.0))
        } else {
            term / (v + n) * z / (n + 1.0)
        };
        sum = sum + term;
        // not before the term with the smallest |v + n|, e.g. 1 / (-7.000001)_8 ~ 1e6 / 7!
        if size(term) <= size(sum) * (f64::EPSILON / 2.0) && (v + n > -0.5 || size(term) == 0.0) {
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
    } else {
        // Unlike in xsf, the Maclaurin series is summed here like for real `z`, instead of using
        // the Bessel functions, which lose digits for large `|v|`, e.g. ~eps ln(Gamma(v)) through
        // the logarithmic prefactor. An overflow of the series is kept, e.g. for `v = 5e-324`, as
        // within its rounding of `f64::MAX` the Bessel result is too rough to tell if it overflows.
        if z.norm() <= v.abs().max(1.0) {
            let sum = hyp0f1_series(v, z);
            if !sum.is_nan() {
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
        let x = 2.0 * arg;
        let r_v = bessel(v, x);
        let mut r = bessel_below(v, x, s, || r_v, |n| bessel(n, x));

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
            let shift = s * res / arg * r_v;
            if shift.is_finite() {
                r += shift;
            }
        }

        let (gamma_v, arg_pow) = (crate::gamma(v), arg.powf(1.0 - v));
        let out = r * gamma_v * arg_pow;
        let normal = gamma_v.is_normal() && arg_pow.norm().is_normal() && out.norm().is_normal();
        if !r.is_finite() || r == C64::ZERO {
            // the Bessel function over- or underflows, or rounds to 0 near a zero
            let asy = if v == 1.0 { None } else { hyp0f1_asy_c64(v, z) };
            asy.unwrap_or_else(|| {
                // AMOS fails (NaN) for huge `v` or `|z|`. For `v - 1 > DEBYE_MIN_V1` on the axis,
                // the real path (with `|J| <= 1`, DLMF 10.14.1) tells if the result underflows;
                // elsewhere its phase can be lost, so this keeps the NaN.
                if r.is_nan() && z.im == 0.0 && z.re < 0.0 && v - 1.0 > DEBYE_MIN_V1 {
                    let f = hyp0f1_real(v, z.re);
                    if f == 0.0 {
                        return f.into();
                    }
                }
                out
            })
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
/// This function is not available in SciPy. Both [`f64`] and
/// [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Definition
///
/// $$
/// _0F_0\left[ \middle\| z\right] = \sum\_{n=0}^\infty {z^n \over n!} = e^z
/// $$
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function $_1F_0\left[a\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
#[inline]
pub fn hyp0f0<T: HypergeometricArg>(z: T) -> T {
    z.hyp0f0()
}

/// Hypergeometric function $_1F_0\left[a\middle\| z\right]$ for real or complex $z$
///
/// This function is not available in SciPy. Both [`f64`] and
/// [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Definition
///
/// $$
/// _1F_0\left[a\middle\| z\right] = \sum\_{n=0}^\infty \rpow a n {z^n \over n!} = (1 - z)^{-a}
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)), where the
/// series converges for $\abs z < 1$, or terminates. It is NaN at $z = 1$ if $a > 0$, and for
/// `f64` $z > 1$ if $a$ isn't an integer.
///
/// # See also
/// - [`hyp0f0`]: Hypergeometric function $_0F_0\left[\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
#[inline]
pub fn hyp1f0<T: HypergeometricArg>(a: f64, z: T) -> T {
    z.hyp1f0(a)
}

/// Confluent hypergeometric limit function $_0F_1\left[b\middle\| z\right]$ for real or complex $z$
///
/// Corresponds to [`scipy.special.hyp0f1`][hyp0f1].
///
/// [hyp0f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp0f1.html
///
/// # Definition
///
/// $$
/// _0F_1\left[b\middle\| z\right] = \sum\_{n=0}^\infty {1\over \rpow b n} {z^n \over n!}
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1]. It is NaN
/// if $b$ is a non-positive integer, where it has a pole, or $-\infty$, and if an argument is NaN.
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function, $_1F_0\left[a\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
#[inline]
pub fn hyp0f1<T: HypergeometricArg>(b: f64, z: T) -> T {
    z.hyp0f1(b)
}
