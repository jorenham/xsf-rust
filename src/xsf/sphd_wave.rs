//! The spheroidal wave functions are translated into pure Rust from `xsf/sphd_wave.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

use crate::xsf::specfun::{self, Status, try_zeroed};

/// Characteristic value of prolate spheroidal function
///
/// Corresponds to [`scipy.special.pro_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_cv.html
///
/// # See also
/// - [`oblate_segv`]: Characteristic value of oblate spheroidal function
/// - [`prolate_aswfa`]: Prolate spheroidal angular function of the 1st kind for precomputed
///   characteristic value
/// - [`prolate_radial1`]: Prolate spheroidal radial function of the 1st kind for precomputed
///   characteristic value
/// - [`prolate_radial2`]: Prolate spheroidal radial function of the 2nd kind for precomputed
///   characteristic value
#[doc(alias = "pro_cv")]
#[must_use]
#[inline]
pub fn prolate_segv(m: u64, n: u64, c: f64) -> f64 {
    let (m, n) = (m as f64, n as f64);
    if (n < m) || ((n - m) > 198.0) {
        return f64::NAN;
    }
    segv(m, n, c, 1).unwrap_or(f64::NAN)
}

/// Characteristic value of oblate spheroidal function
///
/// Corresponds to [`scipy.special.obl_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_cv.html
///
/// # See also
/// - [`prolate_segv`]: Characteristic value of prolate spheroidal function
/// - [`oblate_aswfa`]: Oblate spheroidal angular function of the 1st kind for precomputed
///   characteristic value
/// - [`oblate_radial1`]: Oblate spheroidal radial function of the 1st kind for precomputed
///   characteristic value
/// - [`oblate_radial2`]: Oblate spheroidal radial function of the 2nd kind for precomputed
///   characteristic value
#[doc(alias = "obl_cv")]
#[must_use]
#[inline]
pub fn oblate_segv(m: u64, n: u64, c: f64) -> f64 {
    let (m, n) = (m as f64, n as f64);
    if (n < m) || ((n - m) > 198.0) {
        return f64::NAN;
    }
    segv(m, n, c, -1).unwrap_or(f64::NAN)
}

/// Prolate spheroidal angular function of the 1st kind and its derivative
///
/// Computes the prolate spheroidal angular function of the 1st kind and its derivative (w.r.t.
/// $x$).
///
/// Corresponds to [`scipy.special.pro_ang1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_ang1.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument ($|x| < 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`prolate_aswfa`]: Prolate spheroidal angular function for precomputed characteristic value
/// - [`oblate_aswfa_nocv`]: Oblate spheroidal angular function of the 1st kind
#[doc(alias = "pro_ang1")]
#[must_use]
#[inline]
pub fn prolate_aswfa_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    aswfa_nocv(m as f64, n as f64, c, x, 1)
}

/// Oblate spheroidal angular function of the 1st kind and its derivative
///
/// Computes the oblate spheroidal angular function of the 1st kind and its derivative (w.r.t. $x$).
///
/// Corresponds to [`scipy.special.obl_ang1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_ang1.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument ($|x| < 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`oblate_aswfa`]: Oblate spheroidal angular function for precomputed characteristic value
/// - [`oblate_radial1_nocv`]: Oblate spheroidal radial function of the 1st kind
/// - [`oblate_radial2_nocv`]: Oblate spheroidal radial function of the 2nd kind
#[doc(alias = "obl_ang1")]
#[must_use]
#[inline]
pub fn oblate_aswfa_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    aswfa_nocv(m as f64, n as f64, c, x, -1)
}

/// Prolate spheroidal radial function of the 1st kind and its derivative
///
/// Computes the prolate spheroidal radial function of the 1st kind and its derivative (w.r.t. $x$).
///
/// Corresponds to [`scipy.special.pro_rad1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad1.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`prolate_aswfa_nocv`]: Prolate spheroidal angular function of the 1st kind
/// - [`prolate_radial1`]: Prolate spheroidal radial function for precomputed characteristic value
/// - [`prolate_radial2_nocv`]: Prolate spheroidal radial function of the 2nd kind
#[doc(alias = "pro_rad1")]
#[must_use]
#[inline]
pub fn prolate_radial1_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    radial_nocv(m as f64, n as f64, c, x, 1, 1)
}

/// Oblate spheroidal radial function of the 1st kind and its derivative
///
/// Computes the oblate spheroidal radial function of the 1st kind and its derivative (w.r.t. $x$).
///
/// Corresponds to [`scipy.special.obl_rad1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad1.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`oblate_aswfa_nocv`]: Oblate spheroidal angular function of the 1st kind
/// - [`oblate_radial1`]: Oblate spheroidal radial function of the 1st kind for precomputed
///   characteristic value
/// - [`oblate_radial2_nocv`]: Oblate spheroidal radial function of the 2nd kind
#[doc(alias = "obl_rad1")]
#[must_use]
#[inline]
pub fn oblate_radial1_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    radial_nocv(m as f64, n as f64, c, x, -1, 1)
}

/// Prolate spheroidal radial function of the 2nd kind and its derivative
///
/// Computes the prolate spheroidal radial function of the 2nd kind and its derivative (w.r.t. $x$).
///
/// Corresponds to [`scipy.special.pro_rad2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad2.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`prolate_aswfa_nocv`]: Prolate spheroidal angular function of the 1st kind
/// - [`prolate_radial1_nocv`]: Prolate spheroidal radial function of the 1st kind
/// - [`prolate_radial2`]: Prolate spheroidal radial function of the 2nd kind for precomputed
///   characteristic value
#[doc(alias = "pro_rad2")]
#[must_use]
#[inline]
pub fn prolate_radial2_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    radial_nocv(m as f64, n as f64, c, x, 1, 2)
}

/// Oblate spheroidal radial function of the 2nd kind and its derivative
///
/// Computes the oblate spheroidal radial function of the 2nd kind and its derivative (w.r.t. $x$).
///
/// Corresponds to [`scipy.special.obl_rad2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad2.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`oblate_aswfa_nocv`]: Oblate spheroidal angular function of the 1st kind
/// - [`oblate_radial1_nocv`]: Oblate spheroidal radial function of the 1st kind
/// - [`oblate_radial2`]: Oblate spheroidal radial function of the 2nd kind for precomputed
///   characteristic value
#[doc(alias = "obl_rad2")]
#[must_use]
#[inline]
pub fn oblate_radial2_nocv(m: u64, n: u64, c: f64, x: f64) -> (f64, f64) {
    radial_nocv(m as f64, n as f64, c, x, -1, 2)
}

/// Prolate spheroidal angular function for precomputed characteristic value
///
/// Computes the prolate spheroidal angular function of the 1st kind and its derivative (w.r.t.
/// $x$). Requires pre-computed characteristic value.
///
/// Corresponds to [`scipy.special.pro_ang1_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_ang1_cv.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument ($|x| < 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`prolate_aswfa_nocv`]: Prolate spheroidal angular function of the 1st kind
#[doc(alias = "pro_ang1_cv")]
#[must_use]
#[inline]
pub fn prolate_aswfa(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    aswfa_cv(m as f64, n as f64, c, cv, x, 1)
}

/// Oblate spheroidal angular function for precomputed characteristic value
///
/// Computes the oblate spheroidal angular function of the 1st kind and its derivative (w.r.t. $x$).
/// Requires pre-computed characteristic value.
///
/// Corresponds to [`scipy.special.obl_ang1_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_ang1_cv.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument ($|x| < 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`oblate_aswfa_nocv`]: Oblate spheroidal angular function of the 1st kind
#[doc(alias = "obl_ang1_cv")]
#[must_use]
#[inline]
pub fn oblate_aswfa(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    aswfa_cv(m as f64, n as f64, c, cv, x, -1)
}

/// Prolate spheroidal radial function of the 1st kind for precomputed characteristic value
///
/// Computes the prolate spheroidal radial function of the 1st kind and its derivative (w.r.t. $x$).
/// Requires pre-computed characteristic value.
///
/// Corresponds to [`scipy.special.pro_rad1_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad1_cv.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`prolate_radial1_nocv`]: Prolate spheroidal radial function of the 1st kind
/// - [`prolate_radial2`]: Prolate spheroidal radial function of the 2nd kind for precomputed
///   characteristic value
#[doc(alias = "pro_rad1_cv")]
#[must_use]
#[inline]
pub fn prolate_radial1(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    radial_cv(m as f64, n as f64, c, cv, x, 1, 1)
}

/// Oblate spheroidal radial function of the 1st kind for precomputed characteristic value
///
/// Computes the oblate spheroidal radial function of the 1st kind and its derivative (w.r.t. $x$).
/// Requires pre-computed characteristic value.
///
/// Corresponds to [`scipy.special.obl_rad1_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad1_cv.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`oblate_radial1_nocv`]: Oblate spheroidal radial function of the 1st kind
/// - [`oblate_radial2`]: Oblate spheroidal radial function of the 2nd kind for precomputed
///   characteristic value
#[doc(alias = "obl_rad1_cv")]
#[must_use]
#[inline]
pub fn oblate_radial1(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    radial_cv(m as f64, n as f64, c, cv, x, -1, 1)
}

/// Prolate spheroidal radial function of the 2nd kind for precomputed characteristic value
///
/// Computes the prolate spheroidal radial function of the 2nd kind and its derivative (w.r.t. $x$).
/// Requires pre-computed characteristic value.
///
/// Corresponds to [`scipy.special.pro_rad2_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pro_rad2_cv.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`prolate_radial1`]: Prolate spheroidal radial function of the 1st kind for precomputed
///   characteristic value
/// - [`prolate_radial2_nocv`]: Prolate spheroidal radial function of the 2nd kind
#[doc(alias = "pro_rad2_cv")]
#[must_use]
#[inline]
pub fn prolate_radial2(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    radial_cv(m as f64, n as f64, c, cv, x, 1, 2)
}

/// Oblate spheroidal radial function of the 2nd kind for precomputed characteristic value
///
/// Computes the oblate spheroidal radial function of the 2nd kind and its derivative (w.r.t. $x$).
/// Requires pre-computed characteristic value.
///
/// Corresponds to [`scipy.special.obl_rad2_cv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.obl_rad2_cv.html
///
/// # Arguments
/// - `m`: Mode parameter ($m \ge 0$)
/// - `n`: Mode parameter ($n \ge m$)
/// - `c`: Spheroidal parameter
/// - `cv`: Characteristic value
/// - `x`: Real argument ($x > 1$)
///
/// # Returns
/// - $s$: Value of the function
/// - $s\'$: Value of the derivative w.r.t. $x$
///
/// # See also
/// - [`oblate_radial1`]: Oblate spheroidal radial function of the 1st kind for precomputed
///   characteristic value
/// - [`oblate_radial2_nocv`]: Oblate spheroidal radial function of the 2nd kind
#[doc(alias = "obl_rad2_cv")]
#[must_use]
#[inline]
pub fn oblate_radial2(m: u64, n: u64, c: f64, cv: f64, x: f64) -> (f64, f64) {
    radial_cv(m as f64, n as f64, c, cv, x, -1, 2)
}

/// Unlike in xsf, the functions return NaN instead of indexing out of bounds of the fixed-size
/// buffers in `specfun` (or of overflowing `int`), which happens if `0.5 (n - m) + c` is not
/// within `(lo, hi)`, or for the radial functions, if `m + n + 2 c` is roughly larger than 200 or
/// `|c x| >= 2^31 / 1.1`.
fn out_of_bounds(m: f64, n: f64, c: f64, lo: f64, hi: f64) -> bool {
    // `sckb` uses `max(c, 1e-10)`
    let h = 0.5 * (n - m);
    !(lo < h + c && h + c.max(0.0) < hi) || m >= 1e8
}

/// `specfun::segv`, or `None` if an allocation fails
fn segv(m: f64, n: f64, c: f64, kd: i32) -> Option<f64> {
    // `nm = 10 + (int)(0.5 (n - m) + c)` indexes buffers of size 300
    if out_of_bounds(m, n, c, -1e9, 291.0) {
        return Some(f64::NAN);
    }
    let int_m = m as i32;
    let int_n = n as i32;
    let mut eg = try_zeroed((n - m + 2.0) as usize)?;
    let mut cv = 0.0;
    let status = specfun::segv(int_m, int_n, c, kd, &mut cv, &mut eg);
    (status != Status::NoMemory).then_some(cv)
}

/// `specfun::aswfa`, i.e. `prolate_aswfa` (`kd = 1`) or `oblate_aswfa` (`kd = -1`)
fn aswfa(m: f64, n: f64, c: f64, cv: f64, x: f64, kd: i32) -> (f64, f64) {
    // `nm = 25 + (int)(0.5 (n - m) + c)` indexes buffers of size 200 in `sdmn` and `sckb`
    if out_of_bounds(m, n, c, -1e9, 175.0) {
        return (f64::NAN, f64::NAN);
    }
    let (mut s1f, mut s1d) = (f64::NAN, f64::NAN);
    if specfun::aswfa(x, m as i32, n as i32, c, kd, cv, &mut s1f, &mut s1d) == Status::NoMemory {
        return (f64::NAN, f64::NAN);
    }
    (s1f, s1d)
}

fn aswfa_cv(m: f64, n: f64, c: f64, cv: f64, x: f64, kd: i32) -> (f64, f64) {
    if (x >= 1.0) || (x <= -1.0) || (n < m) {
        return (f64::NAN, f64::NAN);
    }
    aswfa(m, n, c, cv, x, kd)
}

fn aswfa_nocv(m: f64, n: f64, c: f64, x: f64, kd: i32) -> (f64, f64) {
    if (x >= 1.0) || (x <= -1.0) || (n < m) || ((n - m) > 198.0) {
        return (f64::NAN, f64::NAN);
    }
    let Some(cv) = segv(m, n, c, kd) else {
        return (f64::NAN, f64::NAN);
    };
    aswfa(m, n, c, cv, x, kd)
}

/// `specfun::rswfp` (`kd = 1`) or `specfun::rswfo` (`kd = -1`), returning `(r1f, r1d)` if `kf = 1`,
/// or `(r2f, r2d)` if `kf = 2`
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn radial(m: f64, n: f64, c: f64, cv: f64, x: f64, kd: i32, kf: i32) -> (f64, f64) {
    // as in `aswfa`, and `kmn` requires `nn = 25 + (int)(0.5 (n - m) + c) + m >= m + 2`
    if out_of_bounds(m, n, c, -24.0, 175.0) {
        return (f64::NAN, f64::NAN);
    }
    // `sphj` (`rmn1`) and `sphy` (`rmn2l`) write up to index `2 nm + m` in buffers of size 252,
    // with `nm = 25 + (n - m) / 2 + (int)c`, and `msta1(c x, ..)` (`rmn1`) converts `1.1 |c x|` to
    // `int`
    let nm = 25.0 + ((n - m) / 2.0).floor() + c.trunc();
    if 2.0 * nm + m > 251.0 || !((c * x).abs() < 1e9) {
        return (f64::NAN, f64::NAN);
    }
    // `kmn` (`rmn2so`, `rmn2sp`) writes up to index `25 + (int)(0.5 (n - m) + c) + m` in buffers
    // of size 200 or 201, and `lpmns` and `lqmns` (`rmn2sp`) up to `(int)(2 (nm + c) + m) + 1`
    // in buffers of size 252
    if kf == 2
        && (25.0 + (0.5 * (n - m) + c).trunc() + m > 199.0
            || (2.0 * (nm - c.trunc() + c) + m).trunc() > 250.0)
    {
        return (f64::NAN, f64::NAN);
    }
    let rswf = if kd == 1 {
        specfun::rswfp
    } else {
        specfun::rswfo
    };
    let (mut r1f, mut r1d, mut r2f, mut r2d) = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    if rswf(
        m as i32, n as i32, c, x, cv, kf, &mut r1f, &mut r1d, &mut r2f, &mut r2d,
    ) == Status::NoMemory
    {
        return (f64::NAN, f64::NAN);
    }
    if kf == 1 { (r1f, r1d) } else { (r2f, r2d) }
}

/// The domain of `x`: `x > 1` for the prolate (`kd = 1`), and `x >= 0` for the oblate functions
fn radial_x_invalid(x: f64, kd: i32) -> bool {
    if kd == 1 { x <= 1.0 } else { x < 0.0 }
}

fn radial_cv(m: f64, n: f64, c: f64, cv: f64, x: f64, kd: i32, kf: i32) -> (f64, f64) {
    if radial_x_invalid(x, kd) || (n < m) {
        return (f64::NAN, f64::NAN);
    }
    radial(m, n, c, cv, x, kd, kf)
}

fn radial_nocv(m: f64, n: f64, c: f64, x: f64, kd: i32, kf: i32) -> (f64, f64) {
    if radial_x_invalid(x, kd) || (n < m) || ((n - m) > 198.0) {
        return (f64::NAN, f64::NAN);
    }
    let Some(cv) = segv(m, n, c, kd) else {
        return (f64::NAN, f64::NAN);
    };
    radial(m, n, c, cv, x, kd, kf)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

    #[test]
    fn test_prolate_segv() {
        xsref::test("prolate_segv", "d_d_d-d", |x| {
            crate::prolate_segv(x[0] as u64, x[1] as u64, x[2])
        });
    }

    // https://github.com/scipy/xsref/issues/9
    // #[test]
    // fn test_oblate_segv() {
    //     testing::test("oblate_segv", "d_d_d-d", |x| {
    //         crate::oblate_segv(x[0] as u64, x[1] as u64, x[2])
    //     });
    // }

    #[test]
    fn test_prolate_aswfa_nocv() {
        xsref::test("prolate_aswfa_nocv", "d_d_d_d-d_d", |x| {
            crate::prolate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_aswfa_nocv() {
        xsref::test("oblate_aswfa_nocv", "d_d_d_d-d_d", |x| {
            crate::oblate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_radial1_nocv() {
        xsref::test("prolate_radial1_nocv", "d_d_d_d-d_d", |x| {
            crate::prolate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_oblate_radial1_nocv() {
        xsref::test("oblate_radial1_nocv", "d_d_d_d-d_d", |x| {
            crate::oblate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_radial2_nocv() {
        xsref::test("prolate_radial2_nocv", "d_d_d_d-d_d", |x| {
            crate::prolate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
        });
    }

    #[test]
    fn test_prolate_aswfa() {
        xsref::test("prolate_aswfa", "d_d_d_d_d-d_d", |x| {
            crate::prolate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_aswfa() {
        xsref::test("oblate_aswfa", "d_d_d_d_d-d_d", |x| {
            crate::oblate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial1() {
        xsref::test("prolate_radial1", "d_d_d_d_d-d_d", |x| {
            crate::prolate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_oblate_radial1() {
        xsref::test("oblate_radial1", "d_d_d_d_d-d_d", |x| {
            crate::oblate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial2() {
        xsref::test("prolate_radial2", "d_d_d_d_d-d_d", |x| {
            crate::prolate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
        });
    }

    #[test]
    fn test_prolate_radial1_nocv_out_of_bounds() {
        // the 252-element buffers of `rmn1` would overflow
        let (r1f, r1d) = crate::prolate_radial1_nocv(15, 174, 30.2, 1.077);
        assert!(r1f.is_nan() && r1d.is_nan());
    }

    // The xsref `oblate_radial2` cases are computed with scipy, which has the `qstar` and `gmn`
    // bugs that are fixed in xsf-rust.
    #[test]
    fn test_oblate_radial2_wronskian() {
        // `R1 R2' - R1' R2 = 1 / (c (x^2 + 1))`
        for (m, n, c, x) in [
            (0, 0, 1.0, 0.0),
            (1, 1, 1.0, 0.0),
            (0, 1, 1.0, 0.5),
            (1, 1, 1.0, 0.5),
            (1, 2, 1.0, 0.5),
            (2, 3, 2.0, 0.2),
            (3, 6, 0.7, 0.05),
            (5, 8, 3.0, 0.4),
        ] {
            let (r1f, r1d) = crate::oblate_radial1_nocv(m, n, c, x);
            let (r2f, r2d) = crate::oblate_radial2_nocv(m, n, c, x);
            let w = (r1f * r2d - r1d * r2f) * c * (x * x + 1.0);
            assert!((w - 1.0).abs() < 1e-13, "{m} {n} {c} {x}: {w}");

            let cv = crate::oblate_segv(m, n, c);
            assert_eq!(crate::oblate_radial2(m, n, c, cv, x), (r2f, r2d));
        }
    }

    #[test]
    fn test_prolate_segv_c0_large_n() {
        // `n (n + 1)`, without `int` overflow
        assert_eq!(crate::prolate_segv(46341, 46341, 0.0), 46341.0 * 46342.0);
    }
}
