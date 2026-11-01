//! Translated into pure Rust from `xsf/sph_bessel.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the order `n` is always 64-bit (C++ `long` is 32-bit on Windows), and the loop
//! counters of the recurrences are also 64-bit (instead of `int`, which could overflow).

use core::f64::consts::{FRAC_PI_2, PI};

use num_complex::{Complex64, c64};

use crate::xsf::bessel::{
    cyl_bessel_i_c64, cyl_bessel_j_c64, cyl_bessel_j_f64, cyl_bessel_k_c64, cyl_bessel_k_f64,
    cyl_bessel_y_c64,
};
use crate::xsf::complex::{cdiv, cmul, csqrt};

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
        sph_bessel_j_c64(n, self)
    }
    #[inline]
    fn sph_bessel_y(self, n: i64) -> Self {
        sph_bessel_y_c64(n, self)
    }
    #[inline]
    fn sph_bessel_i(self, n: i64) -> Self {
        sph_bessel_i_c64(n, self)
    }
    #[inline]
    fn sph_bessel_k(self, n: i64) -> Self {
        sph_bessel_k_c64(n, self)
    }
    #[inline]
    fn sph_bessel_j_prime(self, n: i64) -> Self {
        sph_bessel_j_jac_c64(n, self)
    }
    #[inline]
    fn sph_bessel_y_prime(self, n: i64) -> Self {
        sph_bessel_y_jac_c64(n, self)
    }
    #[inline]
    fn sph_bessel_i_prime(self, n: i64) -> Self {
        sph_bessel_i_jac_c64(n, self)
    }
    #[inline]
    fn sph_bessel_k_prime(self, n: i64) -> Self {
        sph_bessel_k_jac_c64(n, self)
    }
}

// Translated from `xsf/sph_bessel.h`

/// `n + 1 / static_cast<T>(2)`
#[allow(clippy::cast_precision_loss)]
fn half(n: i64) -> f64 {
    n as f64 + 0.5
}

/// `static_cast<T>(n + 1)`
#[allow(clippy::cast_precision_loss)]
fn np1(n: i64) -> f64 {
    n.wrapping_add(1) as f64
}

/// `std::pow(-1, n)`
#[allow(clippy::cast_precision_loss)]
fn pow_m1(n: i64) -> f64 {
    (-1.0_f64).powf(n as f64)
}

/// The leading term `x^n / (2n + 1)!!` of the power series of both `j_n(x)` and `i_n(x)`, if the
/// first omitted correction `x^2 / (4n + 6)` is below roughly `eps / 4`.
#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
fn sph_bessel_tiny_x(n: i64, x: f64) -> Option<f64> {
    if x * x >= (n as f64 + 1.5) * f64::EPSILON {
        return None;
    }
    Some(sph_bessel_lead(n, x))
}

/// `x^n / (2n + 1)!!`, or `inf` if the product overflows
#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
fn sph_bessel_lead(n: i64, x: f64) -> f64 {
    let mut res = 1.0_f64;
    let mut k: i64 = 1;
    while k <= n && res != 0.0 && res.is_finite() {
        res *= x / (2.0 * k as f64 + 1.0);
        k += 1;
    }
    res
}

/// The power series of `j_n(x)` (DLMF 10.53.1, `s = -1`) or `i_n(x)` (DLMF 10.53.3, `s = 1`), if
/// `x^2 < 2n + 3`, so that each term is at most half the previous one.
///
/// Then the leading term can only overflow (for `x > 1e3`) if the result underflows.
#[allow(clippy::cast_precision_loss)]
fn sph_bessel_series(n: i64, x: f64, s: f64) -> Option<f64> {
    let nf = n as f64;
    if x * x >= 2.0 * nf + 3.0 {
        return None;
    }
    let lead = sph_bessel_lead(n, x);
    if lead.is_infinite() {
        return None;
    }
    let mut term = 1.0_f64;
    let mut sum = 1.0;
    let mut k = 1.0;
    while term.abs() > f64::EPSILON * sum {
        term *= s * x * x / (2.0 * k * (2.0 * (nf + k) + 1.0));
        sum += term;
        k += 1.0;
    }
    Some(lead * sum)
}

/// The leading term `n x^(n-1) / (2n + 1)!!` of the power series of both `j_n'(x)` and
/// `i_n'(x)` for `n >= 1`, if `x` is tiny enough for `sph_bessel_tiny_x(n - 1, x)`.
#[allow(clippy::cast_precision_loss)]
fn sph_bessel_prime_tiny_x(n: i64, x: f64) -> Option<f64> {
    if n < 1 || x.is_nan() {
        return None;
    }
    let nf = n as f64;
    // f_n'(-x) = (-1)^(n+1) f_n'(x), which also gives the right sign if the product underflows
    let sign = if x < 0.0 && n % 2 == 0 { -1.0 } else { 1.0 };
    sph_bessel_tiny_x(n - 1, x.abs()).map(|res| sign * nf * res / (2.0 * nf + 1.0))
}

/// `std::sqrt(M_PI_2 / z)` for complex `z`
fn sqrt_pi_2_div(z: Complex64) -> Complex64 {
    csqrt(cdiv(c64(FRAC_PI_2, 0.0), z))
}

#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
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
        return if n % 2 == 0 { r } else { -r };
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

    // Unlike in xsf (scipy/xsf#359), the leading power series term is also used for j_n (see
    // sph_bessel_i_f64), since J_{n + 1/2}(x) underflows for tiny x, e.g. j_1(1.166e-222) was 0
    // instead of 3.89e-223.
    if let Some(res) = sph_bessel_tiny_x(n, x) {
        return res;
    }

    if n > 0 && n as f64 >= x {
        let jv = cyl_bessel_j_f64(half(n), x);
        // Unlike in xsf, use the power series if J_{n + 1/2}(x) underflowed (AMOS returns 0 below
        // about 4.5e-305), e.g. j_50(3.6e-5) was 0 instead of 2.37e-303.
        if jv.abs() < f64::MIN_POSITIVE {
            if let Some(res) = sph_bessel_series(n, x, -1.0) {
                return res;
            }
        }
        return (FRAC_PI_2 / x).sqrt() * jv;
    }

    let mut s0 = x.sin() / x;
    if n == 0 {
        return s0;
    }

    let mut s1 = (s0 - x.cos()) / x;
    if n == 1 {
        return s1;
    }

    let mut sn = 0.0;
    let mut i: i64 = 0;
    while i < n - 1 {
        sn = (2 * i + 3) as f64 * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if sn.is_infinite() {
            // Overflow occurred already : terminate recurrence.
            return sn;
        }
        i += 1;
    }

    sn
}

#[allow(clippy::float_cmp)]
fn sph_bessel_j_c64(n: i64, z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return z;
    }

    if n < 0 {
        return c64(f64::NAN, 0.0);
    }

    if z.re.is_infinite() {
        // https://dlmf.nist.gov/10.52.E3
        if z.im == 0.0 {
            return c64(0.0, 0.0);
        }

        return c64(1.0, 1.0) * f64::INFINITY;
    }

    if z.re == 0.0 && z.im == 0.0 {
        if n == 0 {
            return c64(1.0, 0.0);
        }

        return c64(0.0, 0.0);
    }

    // For purely real negative z, the principal-branch product
    // sqrt(M_PI_2/z) * J_{n+1/2}(z) introduces a spurious sign relative
    // to the (entire) function j_n. Use the reflection identity
    // j_n(-z) = (-1)^n * j_n(z) explicitly (DLMF 10.47.14).
    if z.im == 0.0 && z.re < 0.0 {
        let r = sph_bessel_j_c64(n, c64(-z.re, 0.0));
        return if n % 2 == 0 { r } else { -r };
    }

    let out = cmul(sqrt_pi_2_div(z), cyl_bessel_j_c64(half(n), z));
    if z.im == 0.0 {
        return c64(out.re, 0.0); // Small imaginary part is spurious
    }

    out
}

#[allow(clippy::cast_precision_loss)]
fn sph_bessel_j_jac_f64(n: i64, z: f64) -> f64 {
    if n == 0 {
        return -sph_bessel_j_f64(1, z);
    }

    if z == 0.0 {
        // DLMF 10.51.2 doesn't work, so use 10.51.1 to get the exact value
        if n == 1 {
            return 1.0 / 3.0;
        }

        return 0.0;
    }

    // Unlike in xsf (scipy/xsf#373), use the leading power series term n x^(n-1) / (2n + 1)!! for
    // tiny x, e.g. j_2'(7.5e-179) was 2.5e-179 instead of 1e-179.
    if let Some(res) = sph_bessel_prime_tiny_x(n, z) {
        return res;
    }

    // Unlike in xsf, use DLMF 10.51.1 if j_n(x) underflows, since (n + 1) j_n(x) / x then loses
    // precision, e.g. j_74'(0.0037) was 0 instead of 3.66e-307.
    let jn = sph_bessel_j_f64(n, z);
    if jn.abs() < f64::MIN_POSITIVE {
        let jm = n as f64 * sph_bessel_j_f64(n - 1, z);
        let jp = np1(n) * sph_bessel_j_f64(n.saturating_add(1), z);
        return (jm - jp) / (2.0 * n as f64 + 1.0);
    }

    // DLMF 10.51.2
    sph_bessel_j_f64(n - 1, z) - np1(n) * jn / z
}

fn sph_bessel_j_jac_c64(n: i64, z: Complex64) -> Complex64 {
    if n == 0 {
        return -sph_bessel_j_c64(1, z);
    }

    if z == c64(0.0, 0.0) {
        // DLMF 10.51.2 doesn't work, so use 10.51.1 to get the exact value
        if n == 1 {
            return cdiv(c64(1.0, 0.0), c64(3.0, 0.0));
        }

        return c64(0.0, 0.0);
    }

    // DLMF 10.51.2
    sph_bessel_j_c64(n - 1, z) - cdiv(cmul(c64(np1(n), 0.0), sph_bessel_j_c64(n, z)), z)
}

#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
fn sph_bessel_y_f64(n: i64, x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }

    if n < 0 {
        return f64::NAN;
    }

    if x < 0.0 {
        let r = sph_bessel_y_f64(n, -x);
        return if n % 2 == 1 { r } else { -r };
    }

    if x.is_infinite() {
        return 0.0;
    }

    if x == 0.0 {
        return f64::NEG_INFINITY;
    }

    let mut s0 = -x.cos() / x;
    if n == 0 {
        return s0;
    }

    let mut s1 = (s0 - x.sin()) / x;
    if n == 1 {
        return s1;
    }

    let mut sn = 0.0;
    let mut idx: i64 = 0;
    while idx < n - 1 {
        sn = (2 * idx + 3) as f64 * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if sn.is_infinite() {
            // Overflow occurred already: terminate recurrence.
            return sn;
        }
        idx += 1;
    }

    sn
}

#[allow(clippy::float_cmp)]
fn sph_bessel_y_c64(n: i64, z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return z;
    }

    if n < 0 {
        return c64(f64::NAN, 0.0);
    }

    if z.re == 0.0 && z.im == 0.0 {
        // https://dlmf.nist.gov/10.52.E2
        return c64(f64::NAN, 0.0);
    }

    if z.re.is_infinite() {
        // https://dlmf.nist.gov/10.52.E3
        if z.im == 0.0 {
            return c64(0.0, 0.0);
        }

        return c64(1.0, 1.0) * f64::INFINITY;
    }

    // y_n(-z) = (-1)^(n+1) * y_n(z); see sph_bessel_j note above.
    if z.im == 0.0 && z.re < 0.0 {
        let r = sph_bessel_y_c64(n, c64(-z.re, 0.0));
        return if n % 2 == 0 { -r } else { r };
    }

    cmul(sqrt_pi_2_div(z), cyl_bessel_y_c64(half(n), z))
}

fn sph_bessel_y_jac_f64(n: i64, x: f64) -> f64 {
    if n == 0 {
        return -sph_bessel_y_f64(1, x);
    }

    // Unlike in xsf (scipy/xsf#373), if y_n(x) is infinite, return -(n + 1) y_n(x) / x instead of
    // inf - inf = NaN, e.g. y_8'(7e-46) was NaN instead of inf.
    //
    // Here x = -0 is treated as +0, as in y_n(x).
    let yn = sph_bessel_y_f64(n, x);
    if yn.is_infinite() {
        return if x < 0.0 { yn } else { -yn };
    }

    sph_bessel_y_f64(n - 1, x) - np1(n) * yn / x
}

fn sph_bessel_y_jac_c64(n: i64, x: Complex64) -> Complex64 {
    if n == 0 {
        return -sph_bessel_y_c64(1, x);
    }

    sph_bessel_y_c64(n - 1, x) - cdiv(cmul(c64(np1(n), 0.0), sph_bessel_y_c64(n, x)), x)
}

#[allow(clippy::float_cmp)]
fn sph_bessel_i_f64(n: i64, x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }

    if n < 0 {
        return f64::NAN;
    }

    if x < 0.0 {
        let r = sph_bessel_i_f64(n, -x);
        return if n % 2 == 0 { r } else { -r };
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
    if let Some(res) = sph_bessel_tiny_x(n, x) {
        return res;
    }

    let iv = crate::bessel_i(half(n), x);
    // Unlike in xsf, use the power series if I_{n + 1/2}(x) underflowed (cephes iv returns 0 if
    // K_{n + 3/2}(x) overflows), e.g. i_39(6.36e-7) was 0 instead of 2.71e-301.
    if iv.abs() < f64::MIN_POSITIVE {
        if let Some(res) = sph_bessel_series(n, x, 1.0) {
            return res;
        }
    }
    (FRAC_PI_2 / x).sqrt() * iv
}

#[allow(clippy::float_cmp)]
fn sph_bessel_i_c64(n: i64, z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return z;
    }

    if n < 0 {
        return c64(f64::NAN, 0.0);
    }

    if z.norm() == 0.0 {
        // https://dlmf.nist.gov/10.52.E1
        if n == 0 {
            return c64(1.0, 0.0);
        }

        return c64(0.0, 0.0);
    }

    if z.re.is_infinite() || z.im.is_infinite() {
        // https://dlmf.nist.gov/10.52.E5
        if z.im == 0.0 {
            if z.re == f64::NEG_INFINITY {
                return c64(pow_m1(n) * f64::INFINITY, 0.0);
            }

            return c64(f64::INFINITY, 0.0);
        }

        return c64(f64::NAN, 0.0);
    }

    // i_n(-z) = (-1)^n * i_n(z); see sph_bessel_j note above.
    if z.im == 0.0 && z.re < 0.0 {
        let r = sph_bessel_i_c64(n, c64(-z.re, 0.0));
        return if n % 2 == 0 { r } else { -r };
    }

    cmul(sqrt_pi_2_div(z), cyl_bessel_i_c64(half(n), z))
}

#[allow(clippy::cast_precision_loss)]
fn sph_bessel_i_jac_f64(n: i64, z: f64) -> f64 {
    if n == 0 {
        return sph_bessel_i_f64(1, z);
    }

    if z == 0.0 {
        if n == 1 {
            return 1. / 3.;
        }
        return 0.0;
    }

    // Unlike in xsf (scipy/xsf#373), use the leading power series term for tiny x (see
    // sph_bessel_j_jac_f64), e.g. i_2'(7.5e-179) was 2.5e-179 instead of 1e-179.
    if let Some(res) = sph_bessel_prime_tiny_x(n, z) {
        return res;
    }

    // Unlike in xsf, use DLMF 10.51.4 if i_n(x) underflows (see sph_bessel_j_jac_f64).
    let i_n = sph_bessel_i_f64(n, z);
    if i_n.abs() < f64::MIN_POSITIVE {
        let im = n as f64 * sph_bessel_i_f64(n - 1, z);
        let ip = np1(n) * sph_bessel_i_f64(n.saturating_add(1), z);
        return (im + ip) / (2.0 * n as f64 + 1.0);
    }

    sph_bessel_i_f64(n - 1, z) - np1(n) * i_n / z
}

fn sph_bessel_i_jac_c64(n: i64, z: Complex64) -> Complex64 {
    if n == 0 {
        return sph_bessel_i_c64(1, z);
    }

    if z == c64(0.0, 0.0) {
        if n == 1 {
            return c64(1. / 3., 0.0);
        }
        return c64(0.0, 0.0);
    }

    sph_bessel_i_c64(n - 1, z) - cdiv(cmul(c64(np1(n), 0.0), sph_bessel_i_c64(n, z)), z)
}

#[allow(clippy::float_cmp)]
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
        let sign = if n.abs() % 2 == 0 { 1.0 } else { -1.0 };
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

#[allow(clippy::float_cmp)]
fn sph_bessel_k_c64(n: i64, z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return z;
    }

    if n < 0 {
        return c64(f64::NAN, 0.0);
    }

    if z.norm() == 0.0 {
        return c64(f64::NAN, 0.0);
    }

    if z.re.is_infinite() || z.im.is_infinite() {
        // https://dlmf.nist.gov/10.52.E6
        if z.im == 0.0 {
            if z.re == f64::INFINITY {
                return c64(0.0, 0.0);
            }

            return c64(f64::NEG_INFINITY, 0.0);
        }

        return c64(f64::NAN, 0.0);
    }

    if z.im == 0.0 {
        return c64(sph_bessel_k_f64(n, z.re), 0.0);
    }

    cmul(sqrt_pi_2_div(z), cyl_bessel_k_c64(half(n), z))
}

fn sph_bessel_k_jac_f64(n: i64, x: f64) -> f64 {
    if n == 0 {
        return -sph_bessel_k_f64(1, x);
    }

    -sph_bessel_k_f64(n - 1, x) - np1(n) * sph_bessel_k_f64(n, x) / x
}

fn sph_bessel_k_jac_c64(n: i64, x: Complex64) -> Complex64 {
    if n == 0 {
        return -sph_bessel_k_c64(1, x);
    }

    -sph_bessel_k_c64(n - 1, x) - cdiv(cmul(c64(np1(n), 0.0), sph_bessel_k_c64(n, x)), x)
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
