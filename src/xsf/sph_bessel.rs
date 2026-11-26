//! Translated into pure Rust from `xsf/cpu/sph_bessel.h` (scipy/xsf@316ab68).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the order `n` is always 64-bit (C++ `long` is 32-bit on Windows).

use core::f64::consts::{FRAC_PI_2, LOG2_E, PI};

use num_complex::Complex64;

use crate::utils::{neg_if_odd, neg_one_pow};
use crate::xsf::amos::MACH;
use crate::xsf::bessel::{
    cyl_bessel_i_c64, cyl_bessel_j_c64, cyl_bessel_j_f64, cyl_bessel_k_c64, cyl_bessel_k_f64,
    cyl_bessel_y_c64,
};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::{LOG2, frexp_f64};
use crate::xsf::complex::Num;

pub trait SphBesselArg: crate::sealed::Sealed {
    fn sph_bessel_j(self, n: i64) -> Self;
    fn sph_bessel_y(self, n: i64) -> Self;
    fn sph_bessel_i(self, n: i64) -> Self;
    fn sph_bessel_k(self, n: i64) -> Self;
    fn sph_bessel_j_prime(self, n: i64) -> Self;
    fn sph_bessel_y_prime(self, n: i64) -> Self;
    fn sph_bessel_i_prime(self, n: i64) -> Self;
    fn sph_bessel_k_prime(self, n: i64) -> Self;
}

impl SphBesselArg for f64 {
    #[inline]
    fn sph_bessel_j(self, n: i64) -> f64 {
        sph_bessel_j_f64(n, self)
    }
    #[inline]
    fn sph_bessel_y(self, n: i64) -> f64 {
        sph_bessel_y_f64(n, self)
    }
    #[inline]
    fn sph_bessel_i(self, n: i64) -> f64 {
        sph_bessel_i_f64(n, self)
    }
    #[inline]
    fn sph_bessel_k(self, n: i64) -> f64 {
        sph_bessel_k_f64(n, self)
    }
    #[inline]
    fn sph_bessel_j_prime(self, n: i64) -> f64 {
        sph_bessel_j_jac_f64(n, self)
    }
    #[inline]
    fn sph_bessel_y_prime(self, n: i64) -> f64 {
        sph_bessel_y_jac_f64(n, self)
    }
    #[inline]
    fn sph_bessel_i_prime(self, n: i64) -> f64 {
        sph_bessel_i_jac_f64(n, self)
    }
    #[inline]
    fn sph_bessel_k_prime(self, n: i64) -> f64 {
        sph_bessel_k_jac_f64(n, self)
    }
}

impl SphBesselArg for Complex64 {
    #[inline]
    fn sph_bessel_j(self, n: i64) -> Self {
        self.via(|z| sph_bessel_j_c64(n, z))
    }
    #[inline]
    fn sph_bessel_y(self, n: i64) -> Self {
        self.via(|z| sph_bessel_y_c64(n, z))
    }
    #[inline]
    fn sph_bessel_i(self, n: i64) -> Self {
        self.via(|z| sph_bessel_i_c64(n, z))
    }
    #[inline]
    fn sph_bessel_k(self, n: i64) -> Self {
        self.via(|z| sph_bessel_k_c64(n, z))
    }
    #[inline]
    fn sph_bessel_j_prime(self, n: i64) -> Self {
        self.via(|z| sph_bessel_j_jac_c64(n, z))
    }
    #[inline]
    fn sph_bessel_y_prime(self, n: i64) -> Self {
        self.via(|z| sph_bessel_y_jac_c64(n, z))
    }
    #[inline]
    fn sph_bessel_i_prime(self, n: i64) -> Self {
        self.via(|z| sph_bessel_i_jac_c64(n, z))
    }
    #[inline]
    fn sph_bessel_k_prime(self, n: i64) -> Self {
        self.via(|z| sph_bessel_k_jac_c64(n, z))
    }
}

// Translated from `xsf/sph_bessel.h`

/// `n + 1 / static_cast<T>(2)`
fn half(n: i64) -> f64 {
    n as f64 + 0.5
}

/// `static_cast<T>(n + 1UL)`
fn np1(n: i64) -> f64 {
    (n as u64).wrapping_add(1) as f64
}

/// `std::pow(-1, n)`
fn pow_m1(n: i64) -> f64 {
    (-1.0_f64).powf(n as f64)
}

/// Leading term `z^n / (2n + 1)!!` of `j_n(z)` and `i_n(z)` as `z -> 0`,
/// <https://dlmf.nist.gov/10.52.E1>
///
/// The extra stop if the product overflows (returning `inf`) is for `sph_bessel_j_series`.
fn sph_bessel_tiny<T: Num>(n: i64, z: T) -> T {
    let mut res = T::from(1.0);
    let mut k: i64 = 1;
    while k <= n && res.abs() != 0.0 && res.is_finite() {
        res = res * (z / (2.0 * k as f64 + 1.0));
        k += 1;
    }
    res
}

/// `j_n'(z)` and `i_n'(z)` for tiny `z`, or `None`
fn sph_bessel_jac_tiny<T: Num>(n: i64, z: T) -> Option<T> {
    // Leading term n z^(n - 1) / (2n + 1)!! of the derivative of https://dlmf.nist.gov/10.52.E1,
    // as DLMF 10.51.2 fails when j_n(z) underflows; the first correction term is below eps / 2 here.
    // |z| < 1 keeps the product from overflowing for huge n, as in sph_bessel_j.
    let (nf, abs_z) = (n as f64, z.abs());
    (abs_z < 1.0 && abs_z * abs_z < (nf + 0.5) * f64::EPSILON)
        .then(|| sph_bessel_tiny(n - 1, z) * (nf / (2.0 * nf + 1.0)))
}

/// The power series of `j_n(x)` (DLMF 10.53.1) if `x^2 < 2n + 3`, so that each term is at most
/// half the previous one, and if the product for the leading term doesn't overflow (else `j_n(x)`
/// underflows).
fn sph_bessel_j_series(n: i64, x: f64) -> Option<f64> {
    let nf = n as f64;
    if x * x >= 2.0 * nf + 3.0 {
        return None;
    }
    let (mut term, mut sum, mut k) = (1.0_f64, 1.0, 1.0);
    while term.abs() > f64::EPSILON * sum {
        term *= -x * x / (2.0 * k * (2.0 * (nf + k) + 1.0));
        sum += term;
        k += 1.0;
    }
    let lead = sph_bessel_tiny(n, x);
    lead.is_finite().then_some(lead * sum)
}

/// `f_n'(x) = f_{n-1}(x) - (n + 1) f_n(x) / x` (DLMF 10.51.2, 10.51.5) for `f = j` (`c = -1`)
/// or `f = i` (`c = 1`).
fn sph_bessel_ji_jac(f: fn(i64, f64) -> f64, c: f64, n: i64, x: f64) -> f64 {
    let nf = n as f64;

    // Unlike in xsf, the sign is from f_n'(-x) = (-1)^(n+1) f_n'(x), so that it's also right if
    // the product underflows.
    if let Some(res) = sph_bessel_jac_tiny(n, x.abs()) {
        let sign = if x < 0.0 && n % 2 == 0 { -1.0 } else { 1.0 };
        return sign * res;
    }

    // Unlike in xsf, use (2n + 1) f_n' = n f_{n-1} + c (n + 1) f_{n+1} (DLMF 10.51.1, 10.51.4)
    // if f_n(x) underflows, since (n + 1) f_n(x) / x then loses precision.
    let (f_m, f_n) = (f(n.wrapping_sub(1), x), f(n, x));
    if f_n.abs() < f64::MIN_POSITIVE {
        let f_p = f(n.saturating_add(1), x);
        return (nf * f_m + c * np1(n) * f_p) / (2.0 * nf + 1.0);
    }

    f_m - np1(n) * f_n / x
}

/// `std::sqrt(M_PI_2 / z)` for complex `z`
fn sqrt_pi_2_div(z: C64) -> C64 {
    (FRAC_PI_2 / z).sqrt()
}

/// `out`, or unlike in xsf, `fallback()` (see [`sph_bessel_overflow`]) where `out` is not finite
fn finite_or(out: C64, fallback: impl FnOnce() -> Option<C64>) -> C64 {
    if out.is_finite() {
        return out;
    }
    fallback().unwrap_or(out)
}

/// `(f_n(z), f_n'(z))` for complex `z` where they overflow, by the upward recurrence
/// `f_(k+1) = (2k + 1) f_k / z + s f_(k-1)` from `g_k = z exp(x) f_k` at `k = -1, 0`, and
/// `f_n' = -s f_(n-1) - (n + 1) f_n / z`, for `y_n` (`s = -1`, `x = -|Im z|`, DLMF 10.51.1 and
/// 10.51.2) and `k_n` (`s = 1`, `x = Re z`, DLMF 10.51.4 and 10.51.5).
///
/// `None` for `n > N_MAX`, as AMOS for `z = 0` and past its range, and where the rounding errors
/// grow by more than `GROWTH_MAX`, about `exp(min(-2x, n (n + 1) / |z|))` (a measured heuristic),
/// as the other solution grows faster for `x < 0`. A part is NaN unless its error interval
/// overflows, or is finite with at least half of its digits correct, as an overflowing part can
/// mask it, e.g. `Re k_2(1e-120 i) = 1e-241` with `Im k_2 = inf`.
fn sph_bessel_overflow(n: i64, z: C64, s: f64, g_m1: C64, g_0: C64, x: f64) -> Option<(C64, C64)> {
    // Past `2^+-E_MAX`, any nonzero part of the result over- or underflows
    const E_MAX: i64 = 1 << 13;
    // The limits on the steps and on the error growth
    const N_MAX: i64 = 1 << 20;
    const GROWTH_MAX: f64 = 1024.0;
    // A bound on the roundings per step relative to `eps |a_k|`, also of `1 / z`
    const ROUNDINGS: f64 = 8.0;

    let r = z.norm();
    let growth = (-2.0 * x).min(n as f64 * np1(n) / r);
    if !((0..=N_MAX).contains(&n) && r > 0.0 && r <= 0.5 / MACH.0 && growth <= GROWTH_MAX.ln()) {
        return None;
    }

    // The binary exponent of the larger part
    let exp2 = |g: C64| frexp_f64(g.linf_norm()).1;
    // `1 / z = 2^-m / w`, as `1 / z` overflows for tiny `z`
    let m = exp2(z).min(0);
    let rw = C64::ONE / z.ldexp(-m);
    // The partwise bound `|a b|` for `a >= 0`
    let abs_mul = |a: C64, b: C64| {
        let (re, im) = (b.re.abs(), b.im.abs());
        c64(a.re * re + a.im * im, a.re * im + a.im * re)
    };
    // An underflow, in units of `2^de`
    let tiny = |de: i64| {
        let t = f64::MIN_POSITIVE.ldexp(de.clamp(-E_MAX, 0) as i32);
        c64(t, t)
    };
    // `(f_0, f_1)` in units of `2^(e - m)`, rescaled so that all parts are below 1, and the new `e`
    let rescale = |f_0: C64, f_1: C64, e: i64| {
        let de = exp2(f_0).max(exp2(f_1));
        (f_0.ldexp(-de), f_1.ldexp(-de), e + i64::from(de - m))
    };

    // `g_k` in units of `2^e`, so that only the result overflows (all `n` steps, as each one turns
    // the phase), and its partwise bound `a_k >= |g_k|` plus the underflow, in units of `2^e_a` as
    // it can grow much faster
    let [a_m1, a_0] = [g_m1, g_0].map(|g| c64(g.re.abs(), g.im.abs()) + tiny(0));
    let init = ((g_m1, g_0, 0), (a_m1, a_0, 0));
    let ((g_m1, g_0, e), (a_m1, a_0, e_a)) =
        (0..n).fold(init, |((g_m1, g_0, e), (a_m1, a_0, e_a)), k| {
            let (c, t) = ((2 * k + 1) as f64, tiny(e - e_a));
            (
                rescale(g_0.ldexp(m), c * g_0 * rw + s * g_m1.ldexp(m), e),
                rescale(a_0.ldexp(m) + t, c * abs_mul(a_0, rw) + a_m1.ldexp(m) + t, e_a),
            )
        });

    // `exp(-x) = 2^q exp(t)` for `|t| <= ln(2) / 2` (reduced in double-double), so that only the
    // final `ldexp` can overflow
    let q = (-x * LOG2_E).round();
    let exp_t = (-f64::from(LOG2 * q + x)).exp();
    // The error is the smaller of the normwise `tol |f|` and the partwise `tol_a |a / z|`
    let [tol, tol_a] = [GROWTH_MAX, ROUNDINGS].map(|c| c * np1(n) * f64::EPSILON);
    // `p 2^e exp(t)` for a part, with its own scale so that a subnormal one rounds only once
    let part = |p: f64, err: f64, e: i64| {
        // `y_n(iy)` is real or imaginary
        if s < 0.0 && z.re == 0.0 && p == 0.0 {
            return p;
        }
        let (p, ep) = frexp_f64(p);
        let (err, e) = (err.ldexp(-ep), (e + i64::from(ep)).clamp(-E_MAX, E_MAX) as i32);
        let [lo, hi] = [-err, err].map(|d| ((p.abs() + d) * exp_t).ldexp(e).is_finite());
        if p.abs() >= err && (!lo || hi && err <= f64::EPSILON.sqrt() * p.abs()) {
            (p * exp_t).ldexp(e)
        } else {
            f64::NAN
        }
    };
    // `f = g 2^e exp(-x) / z`, with the error from `a 2^e_a`
    let d_a = e_a - e;
    let scale = |g: C64, a: C64, e: i64| {
        let de = exp2(g);
        let a =
            (abs_mul(a, rw) + tiny(-d_a)).ldexp((d_a - i64::from(de)).clamp(-E_MAX, E_MAX) as i32);
        let g = g.ldexp(-de) * rw;
        let min = tol * g.linf_norm();
        let e = e + i64::from(de - m) + q as i64;
        c64(part(g.re, min.min(tol_a * a.re), e), part(g.im, min.min(tol_a * a.im), e))
    };
    let dg = -s * g_m1.ldexp(m) - np1(n) * g_0 * rw;
    let da = a_m1.ldexp(m) + np1(n) * abs_mul(a_0, rw);
    Some((scale(g_0, a_0, e), scale(dg, da, e - i64::from(m))))
}

/// `f_n(x)` for `n >= 0` and `x > 0`, from `f_0 = s0` and `f_1 = s1` by the upward recurrence
/// `f_(k+1) = (2k + 1) f_k / x - f_(k-1)` of `j_n` and `y_n` (DLMF 10.51.1), which stops at the
/// first infinite `f_k`.
///
/// Unlike in xsf, `j_n` and `y_n` share this implementation (xsf duplicates it).
fn sph_bessel_recur(n: i64, x: f64, mut s0: f64, mut s1: f64) -> f64 {
    if n == 0 {
        return s0;
    }

    // If s1 overflows (tiny x, only for y_n), so does y_n(x); if s0 does too, the recurrence would
    // give inf - inf.
    if n == 1 || s1.is_infinite() {
        return s1;
    }

    let mut sn = 0.0;
    let mut i: i64 = 0;
    while i < n - 1 {
        sn = (2 * i + 3) as f64 * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if sn.is_infinite() {
            // Overflow occurred already: terminate recurrence.
            return sn;
        }
        i += 1;
    }

    sn
}

fn sph_bessel_j_f64(n: i64, x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }

    if n < 0 {
        return f64::NAN;
    }

    // For real negative x, the principal-branch product
    // sqrt(M_PI_2/x) * J_{n+1/2}(x) introduces a spurious sign relative
    // to the (entire) function j_n. Use the reflection identity
    // j_n(-x) = (-1)^n * j_n(x) explicitly (DLMF 10.47.14).
    if x < 0.0 {
        let r = sph_bessel_j_f64(n, -x);
        return neg_if_odd(n, r);
    }

    if x.is_infinite() {
        return 0.0;
    }

    if x == 0.0 {
        if n == 0 {
            return 1.0;
        }

        return 0.0;
    }

    // Avoid underflow of J_{n + 1/2}(x) as in sph_bessel_i; the first correction term,
    // -x^2 / (4n + 6) (https://dlmf.nist.gov/10.53.E1), is below eps / 4 here.
    // x < 1 keeps the product from overflowing for huge n.
    if x < 1.0 && x * x < (n as f64 + 1.5) * f64::EPSILON {
        return sph_bessel_tiny(n, x);
    }

    if n > 0 && n as f64 >= x {
        let jv = cyl_bessel_j_f64(half(n), x);
        // Unlike in xsf, use the power series if J_{n + 1/2}(x) underflows.
        if jv.abs() < f64::MIN_POSITIVE {
            if let Some(res) = sph_bessel_j_series(n, x) {
                return res;
            }
        }
        return (FRAC_PI_2 / x).sqrt() * jv;
    }

    let s0 = x.sin() / x;
    sph_bessel_recur(n, x, s0, (s0 - x.cos()) / x)
}

fn sph_bessel_j_c64(n: i64, z: C64) -> C64 {
    if z.is_nan() {
        return z;
    }

    if n < 0 {
        return C64::NAN_RE;
    }

    if z.re.is_infinite() {
        // https://dlmf.nist.gov/10.52.E3
        if z.is_real() {
            return C64::ZERO;
        }

        return c64(1.0, 1.0) * f64::INFINITY;
    }

    if z == C64::ZERO {
        if n == 0 {
            return C64::ONE;
        }

        return C64::ZERO;
    }

    // For purely real negative z, the principal-branch product
    // sqrt(M_PI_2/z) * J_{n+1/2}(z) introduces a spurious sign relative
    // to the (entire) function j_n. Use the reflection identity
    // j_n(-z) = (-1)^n * j_n(z) explicitly (DLMF 10.47.14).
    if z.is_real() && z.re < 0.0 {
        let r = sph_bessel_j_c64(n, c64(-z.re, 0.0));
        return neg_if_odd(n, r);
    }

    let out = sqrt_pi_2_div(z) * cyl_bessel_j_c64(half(n), z);
    if z.is_real() {
        return c64(out.re, 0.0); // Small imaginary part is spurious
    }

    out
}

fn sph_bessel_j_jac_f64(n: i64, z: f64) -> f64 {
    if n == 0 {
        return -sph_bessel_j_f64(1, z);
    }

    sph_bessel_ji_jac(sph_bessel_j_f64, -1.0, n, z)
}

fn sph_bessel_j_jac_c64(n: i64, z: C64) -> C64 {
    if n == 0 {
        return -sph_bessel_j_c64(1, z);
    }

    if let Some(res) = sph_bessel_jac_tiny(n, z) {
        return res;
    }

    // DLMF 10.51.2
    sph_bessel_j_c64(n.wrapping_sub(1), z) - c64(np1(n), 0.0) * sph_bessel_j_c64(n, z) / z
}

fn sph_bessel_y_f64(n: i64, x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }

    if n < 0 {
        return f64::NAN;
    }

    if x < 0.0 {
        let r = sph_bessel_y_f64(n, -x);
        return -neg_if_odd(n, r);
    }

    if x.is_infinite() {
        return 0.0;
    }

    if x == 0.0 {
        return f64::NEG_INFINITY;
    }

    let s0 = -x.cos() / x;
    sph_bessel_recur(n, x, s0, (s0 - x.sin()) / x)
}

fn sph_bessel_y_c64(n: i64, z: C64) -> C64 {
    if z.is_nan() {
        return z;
    }

    if n < 0 {
        return C64::NAN_RE;
    }

    if z == C64::ZERO {
        // https://dlmf.nist.gov/10.52.E2
        return C64::NAN_RE;
    }

    if z.re.is_infinite() {
        // https://dlmf.nist.gov/10.52.E3
        if z.is_real() {
            return C64::ZERO;
        }

        return c64(1.0, 1.0) * f64::INFINITY;
    }

    if z.is_real() {
        return c64(sph_bessel_y_f64(n, z.re), 0.0);
    }

    let out = sqrt_pi_2_div(z) * cyl_bessel_y_c64(half(n), z);
    finite_or(out, || sph_bessel_y_overflow(n, z).map(|f| f.0))
}

/// [`sph_bessel_overflow`] of `y_n`, from `y_(-1) = j_0` and `y_0`
fn sph_bessel_y_overflow(n: i64, z: C64) -> Option<(C64, C64)> {
    // `sin z` and `cos z` times `exp(-|y|)`, as they overflow for `|y| > 710`
    let (y, (sin_x, cos_x)) = (z.im.abs(), z.re.sin_cos());
    // `cosh(y) exp(-|y|)` and `sinh(y) exp(-|y|)`
    let (ch, sh) = ((1.0 + (-2.0 * y).exp()) / 2.0, (-(-2.0 * y).exp_m1() / 2.0).copysign(z.im));
    let (sin_z, cos_z) = (c64(sin_x * ch, cos_x * sh), c64(cos_x * ch, -sin_x * sh));
    sph_bessel_overflow(n, z, -1.0, sin_z, -cos_z, -y)
}

fn sph_bessel_y_jac_f64(n: i64, x: f64) -> f64 {
    if n == 0 {
        return -sph_bessel_y_f64(1, x);
    }

    // DLMF 10.51.2; y_{n-1}(x) is negligible if the second term overflows, so avoid inf - inf.
    let res = -np1(n) * sph_bessel_y_f64(n, x) / x;
    if res.is_infinite() {
        return res;
    }

    sph_bessel_y_f64(n.wrapping_sub(1), x) + res
}

fn sph_bessel_y_jac_c64(n: i64, x: C64) -> C64 {
    // Real z: avoid nan + nan j where y_n(z) overflows; complex 0 still gives nan.
    if x.is_real() && x.re != 0.0 {
        return c64(sph_bessel_y_jac_f64(n, x.re), 0.0);
    }

    if n == 0 {
        return -sph_bessel_y_c64(1, x);
    }

    let out =
        sph_bessel_y_c64(n.wrapping_sub(1), x) - c64(np1(n), 0.0) * sph_bessel_y_c64(n, x) / x;
    finite_or(out, || sph_bessel_y_overflow(n, x).map(|f| f.1))
}

fn sph_bessel_i_f64(n: i64, x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }

    if n < 0 {
        return f64::NAN;
    }

    if x < 0.0 {
        let r = sph_bessel_i_f64(n, -x);
        return neg_if_odd(n, r);
    }

    if x == 0.0 {
        // https://dlmf.nist.gov/10.52.E1
        if n == 0 {
            return 1.0;
        }
        return 0.0;
    }

    if x.is_infinite() {
        // https://dlmf.nist.gov/10.49.E8
        if x == f64::NEG_INFINITY {
            return pow_m1(n) * f64::INFINITY;
        }

        return f64::INFINITY;
    }

    // Computing i_n through I_{n + 1/2} can underflow before multiplication
    // by the spherical normalization.  In the tiny-x region, use the leading
    // term of the power series i_n(x) = x^n / (2n + 1)!! * (1 + O(x^2)).
    // The cutoff bounds the omitted first correction by roughly eps / 4.
    if x * x < (n as f64 + 1.5) * f64::EPSILON {
        return sph_bessel_tiny(n, x);
    }

    (FRAC_PI_2 / x).sqrt() * crate::bessel_i(half(n), x)
}

fn sph_bessel_i_c64(n: i64, z: C64) -> C64 {
    if z.is_nan() {
        return z;
    }

    if n < 0 {
        return C64::NAN_RE;
    }

    if z.norm() == 0.0 {
        // https://dlmf.nist.gov/10.52.E1
        if n == 0 {
            return C64::ONE;
        }

        return C64::ZERO;
    }

    if z.re.is_infinite() || z.im.is_infinite() {
        // https://dlmf.nist.gov/10.52.E5
        if z.is_real() {
            if z.re == f64::NEG_INFINITY {
                return c64(pow_m1(n) * f64::INFINITY, 0.0);
            }

            return C64::INFINITY;
        }

        return C64::NAN_RE;
    }

    // i_n(-z) = (-1)^n * i_n(z); see sph_bessel_j note above.
    if z.is_real() && z.re < 0.0 {
        let r = sph_bessel_i_c64(n, c64(-z.re, 0.0));
        return neg_if_odd(n, r);
    }

    sqrt_pi_2_div(z) * cyl_bessel_i_c64(half(n), z)
}

fn sph_bessel_i_jac_f64(n: i64, z: f64) -> f64 {
    if n == 0 {
        return sph_bessel_i_f64(1, z);
    }

    sph_bessel_ji_jac(sph_bessel_i_f64, 1.0, n, z)
}

fn sph_bessel_i_jac_c64(n: i64, z: C64) -> C64 {
    if n == 0 {
        return sph_bessel_i_c64(1, z);
    }

    if let Some(res) = sph_bessel_jac_tiny(n, z) {
        return res;
    }

    sph_bessel_i_c64(n.wrapping_sub(1), z) - c64(np1(n), 0.0) * sph_bessel_i_c64(n, z) / z
}

fn sph_bessel_k_f64(n: i64, z: f64) -> f64 {
    if z.is_nan() {
        return z;
    }

    if n < 0 {
        return f64::NAN;
    }

    if z == 0.0 {
        return f64::INFINITY;
    }

    if z < 0.0 {
        // https://dlmf.nist.gov/10.47#E17
        // and
        // https://dlmf.nist.gov/10.47.E11
        let sign = neg_one_pow(n);
        return -(PI * sph_bessel_i_f64(n, -z) + sign * sph_bessel_k_f64(n, -z));
    }

    if z.is_infinite() {
        // https://dlmf.nist.gov/10.52.E6
        if z == f64::INFINITY {
            return 0.0;
        }

        return f64::NEG_INFINITY;
    }

    (FRAC_PI_2 / z).sqrt() * cyl_bessel_k_f64(half(n), z)
}

fn sph_bessel_k_c64(n: i64, z: C64) -> C64 {
    if z.is_nan() {
        return z;
    }

    if n < 0 {
        return C64::NAN_RE;
    }

    if z.norm() == 0.0 {
        return C64::NAN_RE;
    }

    if z.re.is_infinite() || z.im.is_infinite() {
        // https://dlmf.nist.gov/10.52.E6
        if z.is_real() {
            if z.re == f64::INFINITY {
                return C64::ZERO;
            }

            return C64::NEG_INFINITY;
        }

        return C64::NAN_RE;
    }

    if z.is_real() {
        return c64(sph_bessel_k_f64(n, z.re), 0.0);
    }

    let out = sqrt_pi_2_div(z) * cyl_bessel_k_c64(half(n), z);
    finite_or(out, || sph_bessel_k_overflow(n, z).map(|f| f.0))
}

/// [`sph_bessel_overflow`] of `k_n`, from `k_(-1) = k_0`
fn sph_bessel_k_overflow(n: i64, z: C64) -> Option<(C64, C64)> {
    let g = FRAC_PI_2 * C64::cis(-z.im);
    sph_bessel_overflow(n, z, 1.0, g, g, z.re)
}

fn sph_bessel_k_jac_f64(n: i64, x: f64) -> f64 {
    if n == 0 {
        return -sph_bessel_k_f64(1, x);
    }

    -sph_bessel_k_f64(n.wrapping_sub(1), x) - np1(n) * sph_bessel_k_f64(n, x) / x
}

fn sph_bessel_k_jac_c64(n: i64, x: C64) -> C64 {
    // Real z: avoid nan + nan j where k_n(z) overflows; complex 0 still gives nan.
    if x.is_real() && x.re != 0.0 {
        return c64(sph_bessel_k_jac_f64(n, x.re), 0.0);
    }

    if n == 0 {
        return -sph_bessel_k_c64(1, x);
    }

    let out =
        -sph_bessel_k_c64(n.wrapping_sub(1), x) - c64(np1(n), 0.0) * sph_bessel_k_c64(n, x) / x;
    finite_or(out, || sph_bessel_k_overflow(n, x).map(|f| f.1))
}

/// Spherical Bessel function of the first kind, $j_n(z)$
///
/// Corresponds to [`scipy.special.spherical_jn(n, z)`][jn].
///
/// [jn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_jn.html
///
/// # See also
/// - [`sph_bessel_j_prime`]: derivative $j_n\'(z)$
/// - [`bessel_j`](crate::bessel_j): Bessel function $J_v(z)$
#[doc(alias = "spherical_jn")]
pub fn sph_bessel_j<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_j(n)
}

/// Derivative of [`sph_bessel_j`], $j_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_jn(n, z, derivative=True)`][jn].
///
/// [jn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_jn.html
///
/// # See also
/// - [`sph_bessel_j`]: spherical Bessel function $j_n(z)$
/// - [`bessel_j_prime`](crate::bessel_j_prime): derivative of Bessel function $J_v\'(z)$
#[doc(alias = "spherical_jnp")]
pub fn sph_bessel_j_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_j_prime(n)
}

/// Spherical Bessel function of the second kind, $y_n(z)$
///
/// Corresponds to [`scipy.special.spherical_yn(n, z)`][yn].
///
/// [yn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_yn.html
///
/// # See also
/// - [`sph_bessel_y_prime`]: derivative $y_n\'(z)$
/// - [`bessel_y`](crate::bessel_y): Bessel function $Y_v(z)$
#[doc(alias = "spherical_yn")]
pub fn sph_bessel_y<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_y(n)
}

/// Derivative of [`sph_bessel_y`], $y_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_yn(n, z, derivative=True)`][yn].
///
/// [yn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_yn.html
///
/// # See also
/// - [`sph_bessel_y`]: spherical Bessel function $y_n(z)$
/// - [`bessel_y_prime`](crate::bessel_y_prime): derivative of Bessel function $Y_v\'(z)$
#[doc(alias = "spherical_ynp")]
pub fn sph_bessel_y_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_y_prime(n)
}

/// Modified spherical Bessel function of the first kind, $i_n(z)$
///
/// Corresponds to [`scipy.special.spherical_in(n, z)`][in_].
///
/// [in_]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_in.html
///
/// # See also
/// - [`sph_bessel_i_prime`]: derivative $i_n\'(z)$
/// - [`bessel_i`](crate::bessel_i): modified Bessel function $I_v(z)$
#[doc(alias = "spherical_in")]
pub fn sph_bessel_i<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_i(n)
}

/// Derivative of [`sph_bessel_i`], $i_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_in(n, z, derivative=True)`][in_].
///
/// [in_]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_in.html
///
/// # See also
/// - [`sph_bessel_i`]: spherical modified Bessel function $i_n(z)$
/// - [`bessel_i_prime`](crate::bessel_i_prime): derivative of modified Bessel function $I_v\'(z)$
#[doc(alias = "spherical_inp")]
pub fn sph_bessel_i_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_i_prime(n)
}

/// Modified spherical Bessel function of the second kind, $k_n(z)$
///
/// Corresponds to [`scipy.special.spherical_kn(n, z)`][kn].
///
/// [kn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_kn.html
///
/// # See also
/// - [`sph_bessel_k_prime`]: derivative $k_n\'(z)$
/// - [`bessel_k`](crate::bessel_k): modified Bessel function $K_v(z)$
#[doc(alias = "spherical_kn")]
pub fn sph_bessel_k<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_k(n)
}

/// Derivative of [`sph_bessel_k`], $k_n\'(z)$
///
/// Corresponds to [`scipy.special.spherical_kn(n, z, derivative=True)`][kn].
///
/// [kn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_kn.html
///
/// # See also
/// - [`sph_bessel_k`]: spherical modified Bessel function $k_n(z)$
/// - [`bessel_k_prime`](crate::bessel_k_prime): derivative of modified Bessel function $K_v\'(z)$
#[doc(alias = "spherical_knp")]
pub fn sph_bessel_k_prime<T: SphBesselArg>(n: i64, z: T) -> T {
    z.sph_bessel_k_prime(n)
}
