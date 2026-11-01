use core::ops::{Add, Mul};

use core::f64::consts::{FRAC_2_PI, PI};

use num_complex::{Complex, Complex64, c64};
use num_traits::ToPrimitive;

use crate::utils::PISQ_6;
use crate::xsf::amos::{self, ierr_to_sferr, set_error_and_nan};
use crate::xsf::specfun;
use crate::xsf::trig::{cephes_cospi, cephes_sinpi};

// the arithmetic bounds are needed for `bessel_diff_formula` in `scipy_special::bessel_prime`
pub trait BesselArg:
    crate::sealed::Sealed + Copy + From<f64> + Add<Output = Self> + Mul<Output = Self>
{
    fn bessel_j(self, v: f64) -> Self;
    fn bessel_je(self, v: f64) -> Self;
    fn bessel_y(self, v: f64) -> Self;
    fn bessel_ye(self, v: f64) -> Self;
    fn bessel_i(self, v: f64) -> Self;
    fn bessel_ie(self, v: f64) -> Self;
    fn bessel_k(self, v: f64) -> Self;
    fn bessel_ke(self, v: f64) -> Self;
    fn hankel_1(self, v: f64) -> Complex<f64>;
    fn hankel_1e(self, v: f64) -> Complex<f64>;
    fn hankel_2(self, v: f64) -> Complex<f64>;
    fn hankel_2e(self, v: f64) -> Complex<f64>;
}

impl BesselArg for f64 {
    #[inline]
    fn bessel_j(self, v: f64) -> Self {
        cyl_bessel_j_f64(v, self)
    }
    #[inline]
    fn bessel_je(self, v: f64) -> Self {
        cyl_bessel_je_f64(v, self)
    }
    #[inline]
    fn bessel_y(self, v: f64) -> Self {
        cyl_bessel_y_f64(v, self)
    }
    #[inline]
    fn bessel_ye(self, v: f64) -> Self {
        cyl_bessel_ye_f64(v, self)
    }
    #[inline]
    fn bessel_i(self, v: f64) -> Self {
        crate::xsf::cephes::iv(v, self)
    }
    #[inline]
    fn bessel_ie(self, v: f64) -> Self {
        cyl_bessel_ie_f64(v, self)
    }
    #[inline]
    fn bessel_k(self, v: f64) -> Self {
        cyl_bessel_k_f64(v, self)
    }
    #[inline]
    fn bessel_ke(self, v: f64) -> Self {
        cyl_bessel_ke_f64(v, self)
    }
    #[inline]
    fn hankel_1(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self.into(), 1, 1)
    }
    #[inline]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self.into(), 2, 1)
    }
    #[inline]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self.into(), 1, 2)
    }
    #[inline]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self.into(), 2, 2)
    }
}

impl BesselArg for Complex<f64> {
    #[inline]
    fn bessel_j(self, v: f64) -> Self {
        cyl_bessel_j_c64(v, self)
    }
    #[inline]
    fn bessel_je(self, v: f64) -> Self {
        cyl_bessel_je_c64(v, self)
    }
    #[inline]
    fn bessel_y(self, v: f64) -> Self {
        cyl_bessel_y_c64(v, self)
    }
    #[inline]
    fn bessel_ye(self, v: f64) -> Self {
        cyl_bessel_ye_c64(v, self)
    }
    #[inline]
    fn bessel_i(self, v: f64) -> Self {
        cyl_bessel_i_c64(v, self)
    }
    #[inline]
    fn bessel_ie(self, v: f64) -> Self {
        cyl_bessel_ie_c64(v, self)
    }
    #[inline]
    fn bessel_k(self, v: f64) -> Self {
        cyl_bessel_k_c64(v, self)
    }
    #[inline]
    fn bessel_ke(self, v: f64) -> Self {
        cyl_bessel_ke_c64(v, self)
    }
    #[inline]
    fn hankel_1(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self, 1, 1)
    }
    #[inline]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self, 2, 1)
    }
    #[inline]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self, 1, 2)
    }
    #[inline]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        cyl_hankel_c64(v, self, 2, 2)
    }
}

// Translated from `xsf/bessel.h` (xsf v0.2.2)

fn rotate(z: Complex64, v: f64) -> Complex64 {
    let c = cephes_cospi(v);
    let s = cephes_sinpi(v);

    c64(c * z.re - s * z.im, s * z.re + c * z.im)
}

fn rotate_jy(j: Complex64, y: Complex64, v: f64) -> Complex64 {
    let c = cephes_cospi(v);
    let s = cephes_sinpi(v);

    c * j - s * y
}

#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
fn reflect_jy(jy: &mut Complex64, v: f64) -> bool {
    /* NB: Y_v may be huge near negative integers -- so handle exact
     *     integers carefully
     */
    if v != v.floor() {
        return false;
    }

    let i = (v - 16384.0 * (v / 16384.0).floor()) as i32;
    if i & 1 != 0 {
        *jy = -*jy;
    }
    true
}

#[allow(clippy::float_cmp)]
fn reflect_i(_ik: &mut Complex64, v: f64) -> bool {
    if v != v.floor() {
        return false;
    }

    true /* I is symmetric for integer v */
}

fn rotate_i(i: Complex64, k: Complex64, v: f64) -> Complex64 {
    let s = (v * PI).sin() * FRAC_2_PI;

    i + s * k
}

fn cyl_bessel_ie_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 2;
    let mut sign = 1;
    let mut nz;
    let mut ierr = 0;
    let mut cy = c64(f64::NAN, f64::NAN);
    let mut cy_k = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    nz = amos::besi(z, v, kode, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);

    if sign == -1 && !reflect_i(&mut cy, v) {
        nz = amos::besk(z, v, kode, n, core::slice::from_mut(&mut cy_k), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_k);
        /* adjust scaling to match zbesi */
        cy_k = rotate(cy_k, -z.im / PI);
        if z.re > 0.0 {
            cy_k.re *= (-2.0 * z.re).exp();
            cy_k.im *= (-2.0 * z.re).exp();
        }
        /* v -> -v */
        cy = rotate_i(cy, cy_k, v);
    }

    cy
}

#[allow(clippy::float_cmp)]
fn cyl_bessel_ie_f64(v: f64, x: f64) -> f64 {
    if v != v.floor() && x < 0.0 {
        return f64::NAN;
    }

    cyl_bessel_ie_c64(v, c64(x, 0.0)).re
}

fn cyl_bessel_ke_c64(mut v: f64, z: Complex64) -> Complex64 {
    let mut cy = c64(f64::NAN, f64::NAN);
    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }

    if v < 0.0 {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }

    let n = 1;
    let kode = 2;
    let mut ierr = 0;
    let nz = amos::besk(z, v, kode, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if ierr == 2 && z.re >= 0.0 && z.im == 0.0 {
        /* overflow */
        cy = c64(f64::INFINITY, 0.0);
    }

    cy
}

#[allow(clippy::float_cmp)]
pub(crate) fn cyl_bessel_ke_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return f64::INFINITY;
    }

    let ke = cyl_bessel_ke_c64(v, c64(x, 0.0)).re;
    if ke == f64::INFINITY {
        return cyl_bessel_k_overflow(cyl_bessel_ke_c64, v, x);
    }
    ke
}

/// Unlike in xsf (scipy/xsf#372), retry an AMOS overflow (K_v(x) >~ 1e304) for |v| > 1 with the
/// recurrence K_v(x) = K_{v-2}(x) + 2 (v - 1) / x K_{v-1}(x) (DLMF 10.29.1), which also holds for
/// e^x K_v(x), e.g. `bessel_k(1.5, 1e-205)` was inf instead of 3.96e307.
fn cyl_bessel_k_overflow(k: fn(f64, Complex64) -> Complex64, v: f64, x: f64) -> f64 {
    let v = v.abs();
    if v <= 1.0 {
        return f64::INFINITY;
    }
    let z = c64(x, 0.0);
    k(v - 2.0, z).re + 2.0 * (v - 1.0) / x * k(v - 1.0, z).re
}

#[allow(clippy::float_cmp)]
pub(crate) fn cyl_bessel_i_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 1;
    let mut sign = 1.0;
    let mut nz;
    let mut ierr = 0;
    let mut cy = c64(f64::NAN, f64::NAN);
    let mut cy_k = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }
    if v < 0.0 {
        v = -v;
        sign = -1.0;
    }
    nz = amos::besi(z, v, kode, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if ierr == 2 {
        /* overflow */
        if z.im == 0.0 && (z.re >= 0.0 || v == v.floor()) {
            if z.re < 0.0 && v / 2.0 != (v / 2.0).floor() {
                cy.re = f64::NEG_INFINITY;
            } else {
                cy.re = f64::INFINITY;
            }
            cy.im = 0.0;
        } else {
            cy = cyl_bessel_ie_c64(v * sign, z);
            cy.re *= f64::INFINITY;
            cy.im *= f64::INFINITY;
        }
    }

    if sign == -1.0 && !reflect_i(&mut cy, v) {
        nz = amos::besk(z, v, kode, n, core::slice::from_mut(&mut cy_k), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_k);
        cy = rotate_i(cy, cy_k, v);
    }

    cy
}

pub(crate) fn cyl_bessel_k_c64(mut v: f64, z: Complex64) -> Complex64 {
    let mut cy = c64(f64::NAN, f64::NAN);
    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }

    if v < 0.0 {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }

    let n = 1;
    let kode = 1;
    let mut ierr = 0;
    let nz = amos::besk(z, v, kode, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if ierr == 2 && z.re >= 0.0 && z.im == 0.0 {
        /* overflow */
        cy = c64(f64::INFINITY, 0.0);
    }

    cy
}

#[allow(clippy::float_cmp)]
pub(crate) fn cyl_bessel_k_f64(v: f64, z: f64) -> f64 {
    if z < 0.0 {
        return f64::NAN;
    }

    if z == 0.0 {
        return f64::INFINITY;
    }

    if z > 710.0 * (1.0 + v.abs()) {
        /* Underflow. See uniform expansion https://dlmf.nist.gov/10.41
         * This condition is not a strict bound (it can underflow earlier),
         * rather, we are here working around a restriction in AMOS.
         */
        return 0.0;
    }

    let k = cyl_bessel_k_c64(v, c64(z, 0.0)).re;
    if k == f64::INFINITY {
        return cyl_bessel_k_overflow(cyl_bessel_k_c64, v, z);
    }
    k
}

fn cyl_bessel_je_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 2;
    let mut nz;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy_j = c64(f64::NAN, f64::NAN);
    let mut cy_y = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy_j;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    nz = amos::besj(z, v, kode, n, core::slice::from_mut(&mut cy_j), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_j);
    if sign == -1 && !reflect_jy(&mut cy_j, v) {
        nz = amos::besy(z, v, kode, n, core::slice::from_mut(&mut cy_y), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_y);
        cy_j = rotate_jy(cy_j, cy_y, v);
    }
    cy_j
}

#[allow(clippy::float_cmp)]
fn cyl_bessel_je_f64(v: f64, x: f64) -> f64 {
    if v != v.floor() && x < 0.0 {
        return f64::NAN;
    }

    cyl_bessel_je_c64(v, c64(x, 0.0)).re
}

fn cyl_bessel_ye_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 2;
    let mut nz;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy_y = c64(f64::NAN, f64::NAN);
    let mut cy_j = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy_y;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    nz = amos::besy(z, v, kode, n, core::slice::from_mut(&mut cy_y), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_y);
    if ierr == 2 && z.re >= 0.0 && z.im == 0.0 {
        /* overflow */
        cy_y.re = f64::INFINITY;
        cy_y.im = 0.0;
    }

    if sign == -1 && !reflect_jy(&mut cy_y, v) {
        nz = amos::besj(z, v, kode, n, core::slice::from_mut(&mut cy_j), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_j);
        cy_y = rotate_jy(cy_y, cy_j, -v);
    }
    cy_y
}

fn cyl_bessel_ye_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    cyl_bessel_ye_c64(v, c64(x, 0.0)).re
}

/// Unlike in xsf, the four Hankel functions share this implementation (xsf duplicates it), with
/// `kode = 2` for the exponentially scaled ones and `m` the kind.
#[allow(clippy::float_cmp)]
fn cyl_hankel_c64(mut v: f64, z: Complex64, kode: i32, m: i32) -> Complex64 {
    let n = 1;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }
    if kode == 1 && m == 2 && v == 0.0 && z.re == 0.0 && z.im == 0.0 {
        cy.re = f64::NAN;
        cy.im = f64::INFINITY;
        return cy;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let nz = amos::besh(z, v, kode, m, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if sign == -1 {
        cy = rotate(cy, if m == 1 { v } else { -v });
    }
    cy
}

pub(crate) fn cyl_bessel_j_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 1;
    let mut nz;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy_j = c64(f64::NAN, f64::NAN);
    let mut cy_y = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy_j;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    nz = amos::besj(z, v, kode, n, core::slice::from_mut(&mut cy_j), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_j);
    if ierr == 2 {
        /* overflow */
        cy_j = cyl_bessel_je_c64(v, z);
        cy_j.re *= f64::INFINITY;
        cy_j.im *= f64::INFINITY;
    }

    if sign == -1 && !reflect_jy(&mut cy_j, v) {
        nz = amos::besy(z, v, kode, n, core::slice::from_mut(&mut cy_y), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_y);
        cy_j = rotate_jy(cy_j, cy_y, v);
    }
    cy_j
}

#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
pub(crate) fn cyl_bessel_j_f64(v: f64, x: f64) -> f64 {
    if f64::from(v as i32) != v && x < 0.0 {
        return f64::NAN;
    }

    let res = cyl_bessel_j_c64(v, c64(x, 0.0));
    if res.re.is_nan() {
        /* AMOS returned NaN, possibly due to overflow */
        return crate::xsf::cephes::jv(v, x);
    }

    res.re
}

#[allow(clippy::float_cmp)]
pub(crate) fn cyl_bessel_y_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 1;
    let mut nz;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy_y = c64(f64::NAN, f64::NAN);
    let mut cy_j = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy_y;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }

    if z.re == 0.0 && z.im == 0.0 {
        /* overflow */
        cy_y.re = f64::NEG_INFINITY;
        cy_y.im = 0.0;
    } else {
        nz = amos::besy(z, v, kode, n, core::slice::from_mut(&mut cy_y), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_y);
        if ierr == 2 && z.re >= 0.0 && z.im == 0.0 {
            /* overflow */
            cy_y.re = f64::NEG_INFINITY;
            cy_y.im = 0.0;
        }
    }

    if sign == -1 && !reflect_jy(&mut cy_y, v) {
        nz = amos::besj(z, v, kode, n, core::slice::from_mut(&mut cy_j), &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy_j);
        cy_y = rotate_jy(cy_y, cy_j, -v);
    }
    cy_y
}

pub(crate) fn cyl_bessel_y_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    let res = cyl_bessel_y_c64(v, c64(x, 0.0));
    if res.re.is_nan() {
        return crate::xsf::cephes::yv(v, x);
    }

    res.re
}

// Bessel J

/// Bessel function of the first kind, $J_v(z)$
///
/// Corresponds to [`scipy.special.jv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html
///
/// # See also
/// - [`bessel_j0`]: faster version of $J_0$ for real argument
/// - [`bessel_j1`]: faster version of $J_1$ for real argument
/// - [`bessel_je`]: $J_v$ without the leading exponential factor
/// - [`bessel_j_prime`](crate::bessel_j_prime): derivative $J_v\'(z)$
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "jv", alias = "cyl_bessel_j")]
#[must_use]
#[inline]
pub fn bessel_j<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_j(v)
}

/// Bessel function of the first kind of order 0, $J_0(x)$
///
/// Corresponds to [`scipy.special.j0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j0.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$ of real order and complex argument
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "j0", alias = "cyl_bessel_j0")]
#[must_use]
#[inline]
pub fn bessel_j0(x: f64) -> f64 {
    crate::xsf::cephes::j0(x)
}

/// Bessel function of the first kind of order 1, $J_1(x)$
///
/// Corresponds to [`scipy.special.j1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.j1.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$ of real order and complex argument
/// - [`sph_bessel_j`](crate::sph_bessel_j): spherical Bessel function $j_n(x)$
#[doc(alias = "j1", alias = "cyl_bessel_j1")]
#[must_use]
#[inline]
pub fn bessel_j1(x: f64) -> f64 {
    crate::xsf::cephes::j1(x)
}

/// Exponentially scaled Bessel function of the first kind
///
/// Corresponds to [`scipy.special.jve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jve.html
///
/// # See also
/// - [`bessel_j`]: Bessel function $J_v(z)$
#[doc(alias = "jve", alias = "cyl_bessel_je")]
#[must_use]
#[inline]
pub fn bessel_je<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_je(v)
}

// Bessel Y (Neumann)

/// Bessel function of the second kind of order 0, $Y_0(x)$
///
/// Corresponds to [`scipy.special.y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y0.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$ of real order and complex argument
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "y0", alias = "cyl_bessel_y0", alias = "cyl_neumann_0")]
#[must_use]
#[inline]
pub fn bessel_y0(x: f64) -> f64 {
    crate::xsf::cephes::y0(x)
}

/// Bessel function of the second kind of order 1, $Y_1(x)$
///
/// Corresponds to [`scipy.special.y1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.y1.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$ of real order and complex argument
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "y1", alias = "cyl_bessel_y1", alias = "cyl_neumann_1")]
#[must_use]
#[inline]
pub fn bessel_y1(x: f64) -> f64 {
    crate::xsf::cephes::y1(x)
}

/// Bessel function of the second kind, $Y_v(z)$
///
/// Corresponds to [`scipy.special.yv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yv.html
///
/// # See also
/// - [`bessel_y0`]: faster version of $Y_0$ for real argument
/// - [`bessel_y1`]: faster version of $Y_1$ for real argument
/// - [`bessel_ye`]: $Y_v$ without the leading exponential factor
/// - [`bessel_y_prime`](crate::bessel_y_prime): derivative $Y_v\'(z)$
/// - [`sph_bessel_y`](crate::sph_bessel_y): spherical Bessel function $y_n(x)$
#[doc(alias = "yv", alias = "cyl_neumann", alias = "cyl_bessel_y")]
#[must_use]
#[inline]
pub fn bessel_y<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_y(v)
}

/// Exponentially scaled Bessel function of the second kind
///
/// Corresponds to [`scipy.special.yve`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yve.html
///
/// # See also
/// - [`bessel_y`]: Bessel function $Y_v(z)$
#[doc(alias = "yve", alias = "cyl_neumann_e", alias = "cyl_bessel_ye")]
#[must_use]
#[inline]
pub fn bessel_ye<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ye(v)
}

// Bessel I

/// Modified Bessel function of the first kind of order 0, $I_0(x)$
///
/// Corresponds to [`scipy.special.i0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$ of real order and complex argument
/// - [`bessel_i0e`]: exponentially scaled $I_0$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "i0", alias = "cyl_bessel_i0")]
#[must_use]
#[inline]
pub fn bessel_i0(x: f64) -> f64 {
    crate::xsf::cephes::i0(x)
}

/// Exponentially scaled modified Bessel function of the first kind of order 0, $e^{-|x|} I_0(x)$
///
/// Corresponds to [`scipy.special.i0e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i0e.html
///
/// # See also
/// - [`bessel_i0`]: modified Bessel function $I_0(x)$
/// - [`bessel_ie`]: exponentially scaled $I_v(z)$
#[doc(alias = "i0e", alias = "cyl_bessel_i0e")]
#[must_use]
#[inline]
pub fn bessel_i0e(x: f64) -> f64 {
    crate::xsf::cephes::i0e(x)
}

/// Modified Bessel function of the first kind of order 1, $I_1(x)$
///
/// Corresponds to [`scipy.special.i1`][i1].
///
/// [i1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$ of real order and complex argument
/// - [`bessel_i1e`]: exponentially scaled $I_1$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "i1", alias = "cyl_bessel_i1")]
#[must_use]
#[inline]
pub fn bessel_i1(x: f64) -> f64 {
    crate::xsf::cephes::i1(x)
}

/// Exponentially scaled modified Bessel function of the first kind of order 1, $e^{-|x|} I_1(x)$
///
/// Corresponds to [`scipy.special.i1e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1e.html
///
/// # See also
/// - [`bessel_i1`]: modified Bessel function $I_1(x)$
/// - [`bessel_ie`]: exponentially scaled $I_v(z)$
#[doc(alias = "i1e", alias = "cyl_bessel_i1e")]
#[must_use]
#[inline]
pub fn bessel_i1e(x: f64) -> f64 {
    crate::xsf::cephes::i1e(x)
}

/// Modified Bessel function of the first kind, $I_v(z)$
///
/// Corresponds to [`scipy.special.iv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iv.html
///
/// # See also
/// - [`bessel_i0`]: faster version of $I_0$ for real argument
/// - [`bessel_i1`]: faster version of $I_1$ for real argument
/// - [`bessel_ie`]: $I_v$ without the leading exponential factor
/// - [`bessel_i_prime`](crate::bessel_i_prime): derivative $I_v\'(z)$
/// - [`sph_bessel_i`](crate::sph_bessel_i): spherical modified Bessel function $i_n(x)$
#[doc(alias = "iv", alias = "cyl_bessel_i")]
#[must_use]
#[inline]
pub fn bessel_i<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_i(v)
}

/// Exponentially scaled modified Bessel function of the first kind
///
/// Corresponds to [`scipy.special.ive`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ive.html
///
/// # See also
/// - [`bessel_i`]: modified Bessel function $I_v(z)$
#[doc(alias = "ive", alias = "cyl_bessel_ie")]
#[must_use]
#[inline]
pub fn bessel_ie<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ie(v)
}

// Bessel K

/// Modified Bessel function of the second kind of order 0, $K_0(x)$
///
/// Corresponds to [`scipy.special.k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$ of real order and complex argument
/// - [`bessel_k0e`]: exponentially scaled $K_0$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "k0", alias = "cyl_bessel_k0")]
#[must_use]
#[inline]
pub fn bessel_k0(x: f64) -> f64 {
    crate::xsf::cephes::k0(x)
}

/// Exponentially scaled modified Bessel function of the second kind of order 0, $e^x K_0(x)$
///
/// Corresponds to [`scipy.special.k0e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k0e.html
///
/// # See also
/// - [`bessel_k0`]: modified Bessel function $K_0(x)$
/// - [`bessel_ke`]: exponentially scaled $K_v(z)$
#[doc(alias = "k0e", alias = "cyl_bessel_k0e")]
#[must_use]
#[inline]
pub fn bessel_k0e(x: f64) -> f64 {
    crate::xsf::cephes::k0e(x)
}

/// Modified Bessel function of the second kind of order 1, $K_1(x)$
///
/// Corresponds to [`scipy.special.k1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$ of real order and complex argument
/// - [`bessel_k1e`]: exponentially scaled $K_1$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "k1", alias = "cyl_bessel_k1")]
#[must_use]
#[inline]
pub fn bessel_k1(x: f64) -> f64 {
    crate::xsf::cephes::k1(x)
}

/// Exponentially scaled modified Bessel function of the second kind of order 1, $e^x K_1(x)$
///
/// Corresponds to [`scipy.special.k1e`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1e.html
///
/// # See also
/// - [`bessel_k1`]: modified Bessel function $K_1(x)$
/// - [`bessel_ke`]: exponentially scaled $K_v(z)$
#[doc(alias = "k1e", alias = "cyl_bessel_k1e")]
#[must_use]
#[inline]
pub fn bessel_k1e(x: f64) -> f64 {
    crate::xsf::cephes::k1e(x)
}

/// Modified Bessel function of the second kind, $K_v(z)$
///
/// Corresponds to [`scipy.special.kv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kv.html
///
/// # See also
/// - [`bessel_k0`]: faster version of $K_0$ for real argument
/// - [`bessel_k1`]: faster version of $K_1$ for real argument
/// - [`bessel_ke`]: $K_v$ without the leading exponential factor
/// - [`bessel_k_prime`](crate::bessel_k_prime): derivative $K_v\'(z)$
/// - [`sph_bessel_k`](crate::sph_bessel_k): spherical modified Bessel function $k_n(x)$
#[doc(alias = "kv", alias = "cyl_bessel_k")]
#[must_use]
#[inline]
pub fn bessel_k<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_k(v)
}

/// Exponentially scaled modified Bessel function of the second kind
///
/// Corresponds to [`scipy.special.kve`][kve].
///
/// [kve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kve.html
///
/// # See also
/// - [`bessel_k`]: modified Bessel function $K_v(z)$
#[doc(alias = "kve", alias = "cyl_bessel_ke")]
#[must_use]
#[inline]
pub fn bessel_ke<T: BesselArg>(v: f64, z: T) -> T {
    z.bessel_ke(v)
}

// Zeros

/// Compute $N$ zeros of Bessel functions $J_v(x)$, $J_v\'(x)$, $Y_v(x)$, and $Y_v\'(x)$
///
/// Corresponds to [`scipy.special.jnyn_zeros`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jnyn_zeros.html
///
/// # Panics
/// - Panics if `N` is greater than 1200
#[doc(alias = "jnyn_zeros")]
#[must_use]
#[inline]
pub fn bessel_zeros<const N: usize>(v: u32) -> [[f64; N]; 4] {
    assert!(N <= 1200, "N must be at most 1200");

    let (mut j, mut jp) = ([f64::NAN; N], [f64::NAN; N]);
    let (mut y, mut yp) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 || v > 100_000_000 {
        // like in scipy (where `nt` must be positive), since `jyzo` writes at least 1 zero, and
        // (unlike in scipy) avoids the `int` overflows in `msta1` and `msta2` for huge `v`, where
        // the zeros are NaN anyway (as for `v >~ 5000`)
        return [j, jp, y, yp];
    }
    specfun::jyzo(
        v.try_into().unwrap(),
        N.try_into().unwrap(),
        &mut j,
        &mut jp,
        &mut y,
        &mut yp,
    );
    [j, jp, y, yp]
}

// Hankel 1

/// Hankel function of the first kind, $H_v^{(1)}(z)$
///
/// Corresponds to [`scipy.special.hankel1`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1.html
///
/// # See also
/// - [`hankel_1e`]: $H_v^{(1)}$ without the leading exponential factor
/// - [`hankel_1_prime`](crate::hankel_1_prime): $n$th derivative, ${d^n\over dz^n}H_v^{(1)}(z)$
/// - [`hankel_2`]: Hankel function of the second kind $H_v^{(2)}(z)$
#[doc(alias = "h1v", alias = "hankel1", alias = "cyl_hankel_1")]
#[must_use]
#[inline]
pub fn hankel_1<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_1(v)
}

/// Exponentially scaled Hankel function of the first kind
///
/// Corresponds to [`scipy.special.hankel1e`][h1ve].
///
/// [h1ve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel1e.html
///
/// # See also
/// - [`hankel_1`]: Hankel function $H_v^{(1)}(z)$
#[doc(alias = "h1ve", alias = "hankel1e", alias = "cyl_hankel_1e")]
#[must_use]
#[inline]
pub fn hankel_1e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_1e(v)
}

// Hankel 2

/// Hankel function of the second kind, $H_v^{(2)}(z)$
///
/// Corresponds to [`scipy.special.hankel2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2.html
///
/// # See also
/// - [`hankel_2e`]: $H_v^{(2)}$ without the leading exponential factor
/// - [`hankel_2_prime`](crate::hankel_2_prime): $n$th derivative, ${d^n\over dz^n}H_v^{(2)}(z)$
/// - [`hankel_1`]: Hankel function of the first kind $H_v^{(1)}(z)$
#[doc(alias = "h2v", alias = "hankel2", alias = "cyl_hankel_2")]
#[must_use]
#[inline]
pub fn hankel_2<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_2(v)
}

/// Exponentially scaled Hankel function of the second kind
///
/// Corresponds to [`scipy.special.hankel2e`][h2ve].
///
/// [h2ve]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hankel2e.html
///
/// # See also
/// - [`hankel_2`]: Hankel function $H_v^{(2)}(z)$
#[doc(alias = "h2ve", alias = "hankel2e", alias = "cyl_hankel_2e")]
#[must_use]
#[inline]
pub fn hankel_2e<T: BesselArg>(v: f64, z: T) -> Complex<f64> {
    z.hankel_2e(v)
}

// Integrals

/// Weighted integral of the Bessel function of the first kind,
/// $\int_0^1 x^\lambda \mathop{J}_v(2ax) \dd x$
///
/// Corresponds to [`scipy.special.besselpoly`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.besselpoly.html
#[must_use]
#[inline]
pub fn besselpoly(a: f64, lambda: f64, v: f64) -> f64 {
    crate::xsf::cephes::besselpoly(a, lambda, v)
}

/* Integrals of bessel functions */

#[allow(
    clippy::assign_op_pattern,
    clippy::float_cmp,
    clippy::unreadable_literal
)]
fn itika(x: f64) -> (f64, f64) {
    // =======================================================
    // Purpose: Integrate modified Bessel functions I0(t) and
    //          K0(t) with respect to t from 0 to x
    // Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
    // Output:  TI --- Integration of I0(t) from 0 to x
    //          TK --- Integration of K0(t) from 0 to x
    // =======================================================

    const A: [f64; 10] = [
        0.625,
        1.0078125,
        2.5927734375,
        9.1868591308594,
        4.1567974090576e+1,
        2.2919635891914e+2,
        1.491504060477e+3,
        1.1192354495579e+4,
        9.515939374212e+4,
        9.0412425769041e+5,
    ];
    let el = 0.5772156649015329;

    let mut ti;
    let mut tk = 0.0;
    let mut r;
    let x2;

    if x == 0.0 {
        return (0.0, 0.0);
    } else if x < 20.0 {
        x2 = x * x;
        ti = 1.0;
        r = 1.0;
        for k in 1..=50 {
            r = 0.25 * r * (f64::from(2 * k) - 1.0) / (f64::from(2 * k) + 1.0) / f64::from(k * k)
                * x2;
            ti += r;
            if (r / ti).abs() < 1.0e-12 {
                break;
            }
        }
        ti *= x;
    } else {
        x2 = 0.0;
        ti = 1.0;
        r = 1.0;
        for k in 1..=10 {
            r /= x;
            ti += A[k - 1] * r;
        }
        let rc1 = 1.0 / (2.0 * PI * x).sqrt();
        ti = rc1 * x.exp() * ti;
    }

    if x < 12.0 {
        let e0 = el + (x / 2.0).ln();
        let mut b1 = 1.0 - e0;
        let mut b2 = 0.0;
        let mut rs = 0.0;
        r = 1.0;
        let mut tw = 0.0;
        for k in 1..=50 {
            r = 0.25 * r * (f64::from(2 * k) - 1.0) / (f64::from(2 * k) + 1.0) / f64::from(k * k)
                * x2;
            b1 += r * (1.0 / f64::from(2 * k + 1) - e0);
            rs += 1.0 / f64::from(k);
            b2 += r * rs;
            tk = b1 + b2;
            if ((tk - tw) / tk).abs() < 1.0e-12 {
                break;
            }
            tw = tk;
        }
        tk *= x;
    } else {
        tk = 1.0;
        r = 1.0;
        for k in 1..=10 {
            r = -r / x;
            tk += A[k - 1] * r;
        }
        let rc2 = (PI / (2.0 * x)).sqrt();
        tk = PI / 2.0 - rc2 * tk * (-x).exp();
    }
    (ti, tk)
}

#[allow(clippy::float_cmp, clippy::manual_midpoint, clippy::unreadable_literal)]
fn itjya(x: f64) -> (f64, f64) {
    let mut a = [0.0; 18];
    let el = 0.5772156649015329;
    let eps = 1.0e-12;

    let tj;
    let ty;
    if x == 0.0 {
        tj = 0.0;
        ty = 0.0;
    } else if x <= 20.0 {
        let x2 = x * x;
        let mut tj_ = x;
        let mut r = x;
        for ki in 1..61 {
            let k = f64::from(ki);
            r = -0.25 * r * (2.0 * k - 1.0) / (2.0 * k + 1.0) / f64::from(ki * ki) * x2;
            tj_ += r;
            if r.abs() < tj_.abs() * eps {
                break;
            }
        }
        tj = tj_;
        let ty1 = (el + (x / 2.0).ln()) * tj;
        let mut rs = 0.0;
        let mut ty2 = 1.0;
        r = 1.0;
        for ki in 1..61 {
            let k = f64::from(ki);
            r = -0.25 * r * (2.0 * k - 1.0) / (2.0 * k + 1.0) / f64::from(ki * ki) * x2;
            rs += 1.0 / k;
            let r2 = r * (rs + 1.0 / (2.0 * k + 1.0));
            ty2 += r2;
            if r2.abs() < ty2.abs() * eps {
                break;
            }
        }
        ty = (ty1 - x * ty2) * 2.0 / PI;
    } else {
        let mut a0 = 1.0;
        let mut a1 = 5.0 / 8.0;
        a[0] = a1;

        for k in 1_u32..=16 {
            let kf = f64::from(k);
            let af = (1.5 * (kf + 0.5) * (kf + 5.0 / 6.0) * a1
                - 0.5 * (kf + 0.5) * (kf + 0.5) * (kf - 0.5) * a0)
                / (kf + 1.0);
            a[k as usize] = af;
            a0 = a1;
            a1 = af;
        }
        let mut bf = 1.0;
        let mut r = 1.0;
        for k in 1..=8 {
            r = -r / (x * x);
            bf += a[2 * k - 1] * r;
        }
        let mut bg = a[0] / x;
        r = 1.0 / x;
        for k in 1..=8 {
            r = -r / (x * x);
            bg += a[2 * k] * r;
        }
        let xp = x + 0.25 * PI;
        let rc = (2.0 / (PI * x)).sqrt();
        tj = 1.0 - rc * (bf * xp.cos() + bg * xp.sin());
        ty = rc * (bg * xp.cos() - bf * xp.sin());
    }
    (tj, ty)
}

#[allow(
    clippy::float_cmp,
    clippy::needless_late_init,
    clippy::unreadable_literal
)]
fn ittika(x: f64) -> (f64, f64) {
    // =========================================================
    // Purpose: Integrate [I0(t)-1]/t with respect to t from 0
    //          to x, and K0(t)/t with respect to t from x to ∞
    // Input :  x   --- Variable in the limits  ( x ≥ 0 )
    // Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
    //          TTK --- Integration of K0(t)/t from x to ∞
    // =========================================================

    const C: [f64; 8] = [
        1.625,
        4.1328125,
        1.45380859375e1,
        6.553353881835e1,
        3.6066157150269e2,
        2.3448727161884e3,
        1.7588273098916e4,
        1.4950639538279e5,
    ];

    let el = 0.5772156649015329;
    let mut tti;
    let ttk;
    let mut r;

    if x == 0.0 {
        return (0.0, 1.0e+300);
    }
    if x < 40.0 {
        tti = 1.0;
        r = 1.0;
        for k in 2..=50 {
            r = 0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            tti += r;
            if (r / tti).abs() < 1.0e-12 {
                break;
            }
        }
        tti *= 0.125 * x * x;
    } else {
        tti = 1.0;
        r = 1.0;
        for k in 1..=8 {
            r /= x;
            tti += C[k - 1] * r;
        }
        let rc = x * (2.0 * PI * x).sqrt();
        tti = tti * x.exp() / rc;
    }
    if x <= 12.0 {
        let e0 = (0.5 * (x / 2.0).ln() + el) * (x / 2.0).ln() + PI * PI / 24.0 + 0.5 * el * el;
        let mut b1 = 1.5 - (el + (x / 2.0).ln());
        let mut rs = 1.0;
        r = 1.0;
        for k in 2..=50 {
            r = 0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            rs += 1.0 / f64::from(k);
            let r2 = r * (rs + 1.0 / (2.0 * f64::from(k)) - (el + (x / 2.0).ln()));
            b1 += r2;
            if (r2 / b1).abs() < 1.0e-12 {
                break;
            }
        }
        ttk = e0 - 0.125 * x * x * b1;
    } else {
        let mut ttk_ = 1.0;
        r = 1.0;
        for k in 1..=8 {
            r = -r / x;
            ttk_ += C[k - 1] * r;
        }
        let rc = x * (FRAC_2_PI * x).sqrt();
        ttk = ttk_ * (-x).exp() / rc;
    }
    (tti, ttk)
}

#[allow(
    clippy::assign_op_pattern,
    clippy::float_cmp,
    clippy::unreadable_literal
)]
fn ittjya(x: f64) -> (f64, f64) {
    // =========================================================
    // Purpose: Integrate [1-J0(t)]/t with respect to t from 0
    //          to x, and Y0(t)/t with respect to t from x to ∞
    // Input :  x   --- Variable in the limits  ( x ≥ 0 )
    // Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
    //          TTY --- Integration of Y0(t)/t from x to ∞
    // =========================================================

    let el = 0.5772156649015329;

    let mut ttj;
    let tty;
    let mut r;

    if x == 0.0 {
        ttj = 0.0;
        tty = -1.0e+300;
    } else if x <= 20.0 {
        ttj = 1.0;
        r = 1.0;
        for k in 2..=100 {
            r = -0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            ttj += r;
            if r.abs() < ttj.abs() * 1.0e-12 {
                break;
            }
        }
        ttj *= 0.125 * x * x;
        let e0 = 0.5 * (PISQ_6 - el * el) - (0.5 * (x / 2.0).ln() + el) * (x / 2.0).ln();
        let mut b1 = el + (x / 2.0).ln() - 1.5;
        let mut rs = 1.0;
        r = -1.0;
        for k in 2..=100 {
            r = -0.25 * r * (f64::from(k) - 1.0) / f64::from(k * k * k) * x * x;
            rs += 1.0 / f64::from(k);
            let r2 = r * (rs + 1.0 / (2.0 * f64::from(k)) - (el + (x / 2.0).ln()));
            b1 += r2;
            if r2.abs() < b1.abs() * 1.0e-12 {
                break;
            }
        }
        tty = FRAC_2_PI * (e0 + 0.125 * x * x * b1);
    } else {
        let a0 = (2.0 / (PI * x)).sqrt();
        let mut bj0 = 0.0;
        let mut by0 = 0.0;
        let mut bj1 = 0.0;
        let mut by1 = 0.0;
        for l in 0..=1 {
            let l = f64::from(l);
            let vt = 4.0 * l * l;
            let mut px = 1.0;
            r = 1.0;
            // `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM
            // only does when optimizing
            for k in 1..=14 {
                let k = f64::from(k);
                let (p0, p1) = (4.0 * k - 3.0, 4.0 * k - 1.0);
                r = -0.0078125 * r * (vt - p0 * p0) / (x * k) * (vt - p1 * p1)
                    / ((2.0 * k - 1.0) * x);
                px += r;
                if r.abs() < px.abs() * 1.0e-12 {
                    break;
                }
            }
            let mut qx = 1.0;
            r = 1.0;
            for k in 1..=14 {
                let k = f64::from(k);
                let (q0, q1) = (4.0 * k - 1.0, 4.0 * k + 1.0);
                r = -0.0078125 * r * (vt - q0 * q0) / (x * k) * (vt - q1 * q1)
                    / ((2.0 * k + 1.0) * x);
                qx += r;
                if r.abs() < qx.abs() * 1.0e-12 {
                    break;
                }
            }
            qx = 0.125 * (vt - 1.0) / x * qx;
            let xk = x - (0.25 + 0.5 * l) * PI;
            bj1 = a0 * (px * xk.cos() - qx * xk.sin());
            by1 = a0 * (px * xk.sin() + qx * xk.cos());
            if l == 0.0 {
                bj0 = bj1;
                by0 = by1;
            }
        }
        let t = 2.0 / x;
        let mut g0 = 1.0;
        let mut r0 = 1.0;
        for k in 1..=10 {
            let k = f64::from(k);
            r0 = -r0 * k * k * t * t;
            g0 += r0;
        }
        let mut g1 = 1.0;
        let mut r1 = 1.0;
        for k in 1..=10 {
            let k = f64::from(k);
            r1 = -r1 * k * (k + 1.0) * t * t;
            g1 += r1;
        }
        ttj = 2.0 * g1 * bj0 / (x * x) - g0 * bj1 / x + el + (x / 2.0).ln();
        tty = 2.0 * g1 * by0 / (x * x) - g0 * by1 / x;
    }
    (ttj, tty)
}

/// Integrals of Bessel functions of order 0
///
/// Corresponds to [`scipy.special.itj0y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itj0y0.html
///
/// # Returns
/// - $\int_0^x J_0(t) \dd t$
/// - $\int_0^x Y_0(t) \dd t$
///
/// # See also
/// - [`it2j0y0`]
/// - [`bessel_j0`]
/// - [`bessel_y0`]
#[doc(alias = "itj0y0")]
#[must_use]
#[inline]
pub fn it1j0y0(mut x: f64) -> (f64, f64) {
    /* int(j0(t),t=0..x) */
    /* int(y0(t),t=0..x) */
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (mut j0int, mut y0int) = itjya(x);
    if flag {
        j0int = -j0int;
        y0int = f64::NAN; /* domain error */
    }
    (j0int, y0int)
}

/// Integrals related to Bessel functions of the first kind of order 0
///
/// Corresponds to [`scipy.special.it2j0y0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2j0y0.html
///
/// # Returns
/// - $\int_0^x {1 - J_0(t) \over t} \dd t$
/// - $\int_x^\infty {Y_0(t) \over t} \dd t$
///
/// # See also
/// - [`it1j0y0`]
/// - [`bessel_j0`]
/// - [`bessel_y0`]
#[must_use]
#[inline]
pub fn it2j0y0(mut x: f64) -> (f64, f64) {
    /* int((1-j0(t))/t,t=0..x) */
    /* int(y0(t)/t,t=x..inf) */
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (j0int, mut y0int) = ittjya(x);
    if flag {
        y0int = f64::NAN; /* domain error */
    }
    (j0int, y0int)
}

/// Integrals of modified Bessel functions of order 0
///
/// Corresponds to [`scipy.special.iti0k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iti0k0.html
///
/// # Returns
/// - $\int_0^x I_0(t) \dd t$
/// - $\int_0^x K_0(t) \dd t$
///
/// # See also
/// - [`it2i0k0`]
/// - [`bessel_i0`]
/// - [`bessel_k0`]
#[doc(alias = "iti0k0")]
#[must_use]
#[inline]
pub fn it1i0k0(mut x: f64) -> (f64, f64) {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (mut i0int, mut k0int) = itika(x);
    if flag {
        i0int = -i0int;
        k0int = f64::NAN; /* domain error */
    }
    (i0int, k0int)
}

/// Integrals related to modified Bessel functions of order 0.
///
/// Corresponds to [`scipy.special.it2i0k0`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.it2i0k0.html
///
/// # Returns
/// - $\int_0^x {I_0(t) - 1 \over t} \dd t$
/// - $\int_x^\infty {K_0(t) \over t} \dd t$
///
/// # See also
/// - [`it1i0k0`]
/// - [`bessel_i0`]
/// - [`bessel_k0`]
#[must_use]
#[inline]
pub fn it2i0k0(mut x: f64) -> (f64, f64) {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let (i0int, mut k0int) = ittika(x);
    if flag {
        k0int = f64::NAN; /* domain error */
    }
    (i0int, k0int)
}

// Riccati-Bessel

#[allow(
    clippy::assign_op_pattern,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::similar_names
)]
fn rctj(x: f64, nm: &mut i32, rj: &mut [f64], dj: &mut [f64]) {
    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the first
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of jn(x)  ( n = 0,1,2,... )
    // Output:  RJ(n) --- x·jn(x)
    //          DJ(n) --- [x·jn(x)]'
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    // ========================================================

    // NOTE: rj[0] and rj[1] are always written (for |x| >= 1e-100), so buffers of length 1
    // (n == 0) are out of bounds; riccati_j uses length 2 for N == 1. Otherwise, only
    // rj[0..=*nm] / dj[0..=*nm] are written, with *nm <= n.
    //
    // NOTE: for x NaN, +-inf, or |x| >= 2^31 / 1.1, msta1/msta2 hit their `(int)` UB (see there),
    // so `*nm` may differ from C++, and the `+ 1` / `+ 5` in msta1 overflow (panic in debug builds)
    // for x +-inf and |x| >= ~2^31.
    let n = rj.len() as i32 - 1;

    *nm = n;
    if x.abs() < 1.0e-100 {
        for k in 0..=n {
            rj[k as usize] = 0.0;
            dj[k as usize] = 0.0;
        }
        dj[0] = 1.0;
        return;
    }
    rj[0] = x.sin();
    rj[1] = rj[0] / x - x.cos();
    let rj0 = rj[0];
    let rj1 = rj[1];
    let mut f = 0.0;

    if n >= 2 {
        let mut m = specfun::msta1(x, 200);
        if m < n {
            *nm = m;
        } else {
            m = specfun::msta2(x, n, 15);
        }

        let mut f0 = 0.0;
        let mut f1 = 1.0e-100;

        for k in (0..=m).rev() {
            f = (2.0 * f64::from(k) + 3.0) * f1 / x - f0;
            if k <= *nm {
                rj[k as usize] = f;
            }
            f0 = f1;
            f1 = f;
        }
        let cs = if rj0.abs() > rj1.abs() {
            rj0 / f
        } else {
            rj1 / f0
        };
        for k in 0..=*nm {
            rj[k as usize] = cs * rj[k as usize];
        }
    }
    dj[0] = x.cos();
    for k in 1..=*nm {
        dj[k as usize] = f64::from(-k) * rj[k as usize] / x + rj[(k - 1) as usize];
    }
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]
fn rcty(x: f64, nm: &mut i32, ry: &mut [f64], dy: &mut [f64]) {
    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the second
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of yn(x)
    // Output:  RY(n) --- x·yn(x)
    //          DY(n) --- [x·yn(x)]'
    //          NM --- Highest order computed
    // ========================================================

    // NOTE: ry[0], ry[1] and dy[0], dy[1] are always written (for x >= 1e-60), so buffers of
    // length 1 (n == 0) are out of bounds; riccati_y uses length 2 for N == 1.
    let n = ry.len() as i32 - 1;

    let mut k;
    *nm = n;
    if x < 1.0e-60 {
        for k in 0..=n {
            ry[k as usize] = -1.0e+300;
            dy[k as usize] = 1.0e+300;
        }
        ry[0] = -1.0;
        dy[0] = 0.0;
        return;
    }

    ry[0] = -x.cos();
    ry[1] = ry[0] / x - x.sin();
    let mut rf0 = ry[0];
    let mut rf1 = ry[1];

    k = 2;
    while k <= n {
        let rf2 = (2.0 * f64::from(k) - 1.0) * rf1 / x - rf0;
        if rf2.abs() > 1.0e+300 {
            break;
        }
        ry[k as usize] = rf2;
        rf0 = rf1;
        rf1 = rf2;
        k += 1;
    }

    *nm = k - 1;
    dy[0] = x.sin();
    for k in 1..=*nm {
        dy[k as usize] = f64::from(-k) * ry[k as usize] / x + ry[(k - 1) as usize];
    }
}

/// Compute Riccati-Bessel function of the first kind and derivatives for the first $N$ orders
///
/// The Riccati-Bessel function of the first kind is defined as $x \mathop{j}_n(x)$, where $j_n$ is
/// the spherical Bessel function of the first kind of order $n$.
///
/// This function computes the value and first derivative of the function for all orders up to $N$.
///
/// Corresponds to [`scipy.special.riccati_jn`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_jn.html
///
/// # Returns
/// - $\[j_0(x), \ldots, j_{N-1}(x)\]$
/// - $\[j_0\'(x), \ldots, j_{N-1}\'(x)\]$
///
/// # See also
/// - [`riccati_y`]
/// - [`bessel_j`]
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
#[doc(alias = "rctj", alias = "riccati_jn")]
#[must_use]
#[inline]
pub fn riccati_j<const N: usize>(x: f64) -> ([f64; N], [f64; N]) {
    let (mut rj, mut dj) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 {
        return (rj, dj);
    }
    if N == 1 {
        // like in scipy, since `rctj` writes at least 2 terms
        let (rj2, dj2) = riccati_j::<2>(x);
        rj.copy_from_slice(&rj2[..N]);
        dj.copy_from_slice(&dj2[..N]);
        return (rj, dj);
    }

    if N > 2 && (x.is_nan() || x.abs() > 1e9) {
        // (unlike in scipy) avoids the `int` casts of `|x|` in `msta1` and `msta2`, which are
        // undefined for NaN and `|x| >= 2**31 / 1.1`
        return (rj, dj);
    }

    let mut nm = 0;
    rctj(x, &mut nm, &mut rj, &mut dj);
    (rj, dj)
}

/// Compute Riccati-Bessel function of the second kind and derivatives for the first $N$ orders
///
/// The Riccati-Bessel function of the second kind is defined here as $x \mathop{y}_n(x)$, where
/// $y_n$ is the spherical Bessel function of the second kind of order $n$.
/// *Note that this is in contrast to a common convention that includes a minus sign in the
/// definition.*
///
/// This function computes the value and first derivative of the function for all orders up to $N$.
///
/// Corresponds to [`scipy.special.riccati_yn`][scipy].
///
/// # Returns
/// - $\[y_0(x), \ldots, y_{N-1}(x)\]$
/// - $\[y_0\'(x), \ldots, y_{N-1}\'(x)\]$
///
/// # See also
/// - [`riccati_j`]
/// - [`bessel_y`]
/// - [`sph_bessel_y`](crate::sph_bessel_y)
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.riccati_yn.html
#[doc(alias = "rcty", alias = "riccati_yn")]
#[must_use]
#[inline]
pub fn riccati_y<const N: usize>(x: f64) -> ([f64; N], [f64; N]) {
    let (mut ry, mut dy) = ([f64::NAN; N], [f64::NAN; N]);
    if N == 0 {
        return (ry, dy);
    }
    if N == 1 {
        // like in scipy, since `rcty` writes at least 2 terms
        let (ry2, dy2) = riccati_y::<2>(x);
        ry.copy_from_slice(&ry2[..N]);
        dy.copy_from_slice(&dy2[..N]);
        return (ry, dy);
    }

    let mut nm = 0;
    rcty(x, &mut nm, &mut ry, &mut dy);
    (ry, dy)
}

/// Evaluate the Jahnke-Emden Lambda function $\Lambda_v(x)$ and its derivatives
///
/// The Jahnke-Emden Lambda function is defined as
///
/// $$\Lambda_v(x) = \Gamma(v+1) \left({2 \over x}\right)^v J_v(x)$$
///
/// where $\Gamma$ is the [Gamma function](crate::gamma) and $J_v$ is the
/// [Bessel function of the first kind](crate::bessel_j).
///
/// Corresponds to [`scipy.special.lmbda`][scipy] in SciPy, and uses translations of
/// `xsf::specfun::lamn` and `xsf::specfun::lamv`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lmbda.html
///
/// # Arguments
/// - `v`: Order of the Lambda function as a non-negative real number
/// - `x`: Value at which to evaluate the function and derivatives
///
/// # Returns
/// - `vl`: Values of $\Lambda_{v_i}(x)$ for
///   $v_i = v-\lfloor v\rfloor, v-\lfloor v\rfloor+1, \ldots, v$
/// - `dl`: Derivatives $\Lambda_{v_i}\'(x)$ for
///   $v_i = v-\lfloor v\rfloor, v-\lfloor v\rfloor+1, \ldots, v$
///
/// Both vectors have length $\lfloor v\rfloor+1$, i.e. `(v as usize) + 1`.
/// They are all NaN if `x` is NaN or if $|x| > 10^9$.
///
/// # See also
/// - [`bessel_j`]
/// - [`gamma`](crate::gamma)
/// - [`sph_bessel_j`](crate::sph_bessel_j)
///
/// # Panics
/// - Panics if `v` is negative
/// - Panics if `v` is too large to convert to `c_int`
#[doc(alias = "lmbda")]
#[must_use]
pub fn jahnke_emden_lambda<V: Into<f64>>(v: V, x: f64) -> (Vec<f64>, Vec<f64>) {
    // based on https://github.com/scipy/scipy/blob/51dfbcc/scipy/special/_basic.py#L2008-L2055
    let v: f64 = v.into();
    assert!(v >= 0.0, "v must be non-negative");

    let n = v.to_usize().unwrap();
    // checked before allocating, since `lamv` casts `v` to `int` as well
    assert!(core::ffi::c_int::try_from(n).is_ok(), "v is too large");
    if x.is_nan() || x.abs() > 1e9 {
        // (unlike in scipy) avoids out-of-bounds indexing in `lamn` and `lamv`, where `int` casts
        // of `|x|` (e.g. in `msta1` and `msta2`) are undefined for NaN and `|x| >= 2**31 / 1.1`
        return (vec![f64::NAN; n + 1], vec![f64::NAN; n + 1]);
    }
    let v0 = v - n.to_f64().unwrap();
    let n1 = if n < 1 { 1 } else { n };
    let v1 = v0 + n1.to_f64().unwrap();

    let size = v1.to_usize().unwrap() + 1;
    let (mut vl, mut dl) = (vec![f64::NAN; size], vec![f64::NAN; size]);
    #[allow(clippy::float_cmp)]
    if v == n.to_f64().unwrap() {
        let mut nm = 0;
        let n = v1.to_usize().unwrap().try_into().unwrap();
        specfun::lamn(n, x, &mut nm, &mut vl, &mut dl);
    } else {
        let mut vm = 0.0;
        specfun::lamv(v1, x, &mut vm, &mut vl, &mut dl);
    }

    vl.resize(n + 1, f64::NAN);
    dl.resize(n + 1, f64::NAN);
    (vl, dl)
}
