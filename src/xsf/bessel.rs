use core::ops::{Add, Mul};

use core::f64::consts::PI;

use num_complex::{Complex, Complex64, c64};
use num_traits::ToPrimitive;

use crate::xsf::amos::{self, ierr_to_sferr, set_error_and_nan};
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
        cyl_hankel_1_c64(v, self.into())
    }
    #[inline]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        cyl_hankel_1e_c64(v, self.into())
    }
    #[inline]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        cyl_hankel_2_c64(v, self.into())
    }
    #[inline]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        cyl_hankel_2e_c64(v, self.into())
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
        cyl_hankel_1_c64(v, self)
    }
    #[inline]
    fn hankel_1e(self, v: f64) -> Complex<f64> {
        cyl_hankel_1e_c64(v, self)
    }
    #[inline]
    fn hankel_2(self, v: f64) -> Complex<f64> {
        cyl_hankel_2_c64(v, self)
    }
    #[inline]
    fn hankel_2e(self, v: f64) -> Complex<f64> {
        cyl_hankel_2e_c64(v, self)
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
    let s = (v * PI).sin() * (2.0 / PI);

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
fn cyl_bessel_ke_f64(v: f64, x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    if x == 0.0 {
        return f64::INFINITY;
    }

    cyl_bessel_ke_c64(v, c64(x, 0.0)).re
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

    cyl_bessel_k_c64(v, c64(z, 0.0)).re
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

fn cyl_hankel_1e_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 2;
    let m = 1;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let nz = amos::besh(z, v, kode, m, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if sign == -1 {
        cy = rotate(cy, v);
    }
    cy
}

fn cyl_hankel_2e_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 2;
    let m = 2;
    let mut ierr = 0;
    let mut sign = 1;

    let mut cy = c64(f64::NAN, f64::NAN);
    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }

    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let nz = amos::besh(z, v, kode, m, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if sign == -1 {
        cy = rotate(cy, -v);
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
        return unsafe { crate::ffi::xsf::cephes_jv(v, x) };
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

fn cyl_hankel_1_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 1;
    let m = 1;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }
    if v < 0.0 {
        v = -v;
        sign = -1;
    }
    let nz = amos::besh(z, v, kode, m, n, core::slice::from_mut(&mut cy), &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cy);
    if sign == -1 {
        cy = rotate(cy, v);
    }
    cy
}

#[allow(clippy::float_cmp)]
fn cyl_hankel_2_c64(mut v: f64, z: Complex64) -> Complex64 {
    let n = 1;
    let kode = 1;
    let m = 2;
    let mut ierr = 0;
    let mut sign = 1;
    let mut cy = c64(f64::NAN, f64::NAN);

    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return cy;
    }
    if v == 0.0 && z.re == 0.0 && z.im == 0.0 {
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
        cy = rotate(cy, -v);
    }
    cy
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
    if N == 0 {
        // like in scipy (where `nt` must be positive), since `jyzo` writes at least 1 zero
        return [j, jp, y, yp];
    }
    unsafe {
        crate::ffi::xsf::jyzo(
            v.try_into().unwrap(),
            N.try_into().unwrap(),
            j.as_mut_ptr(),
            jp.as_mut_ptr(),
            y.as_mut_ptr(),
            yp.as_mut_ptr(),
        );
    }
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
pub fn it1j0y0(x: f64) -> (f64, f64) {
    crate::utils::out2(|j0int, y0int| unsafe { crate::ffi::xsf::it1j0y0(x, j0int, y0int) })
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
pub fn it2j0y0(x: f64) -> (f64, f64) {
    crate::utils::out2(|j0int, y0int| unsafe { crate::ffi::xsf::it2j0y0(x, j0int, y0int) })
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
pub fn it1i0k0(x: f64) -> (f64, f64) {
    crate::utils::out2(|i0int, k0int| unsafe { crate::ffi::xsf::it1i0k0(x, i0int, k0int) })
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
pub fn it2i0k0(x: f64) -> (f64, f64) {
    crate::utils::out2(|i0int, k0int| unsafe { crate::ffi::xsf::it2i0k0(x, i0int, k0int) })
}

// Riccati-Bessel

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

    let _ = unsafe { crate::ffi::xsf::rctj(N, x, rj.as_mut_ptr(), dj.as_mut_ptr()) };
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

    let _ = unsafe { crate::ffi::xsf::rcty(N, x, ry.as_mut_ptr(), dy.as_mut_ptr()) };
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
/// Corresponds to [`scipy.special.lmbda`][scipy] in SciPy, and calls the FFI functions
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
    unsafe {
        #[allow(clippy::float_cmp)]
        if v == n.to_f64().unwrap() {
            crate::ffi::xsf::lamn(
                v1.to_usize().unwrap().try_into().unwrap(),
                x,
                vl.as_mut_ptr(),
                dl.as_mut_ptr(),
            );
        } else {
            crate::ffi::xsf::lamv(v1, x, vl.as_mut_ptr(), dl.as_mut_ptr());
        }
    };

    vl.resize(n + 1, f64::NAN);
    dl.resize(n + 1, f64::NAN);
    (vl, dl)
}

// Tests

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;
    use num_complex::c64;

    // Bessel J

    #[test]
    fn test_bessel_j0() {
        xsref::test("cyl_bessel_j0", "d-d", |x| crate::bessel_j0(x[0]));
    }

    #[test]
    fn test_bessel_j1() {
        xsref::test("cyl_bessel_j1", "d-d", |x| crate::bessel_j1(x[0]));
    }

    #[test]
    fn test_bessel_j_f64() {
        xsref::test("cyl_bessel_j", "d_d-d", |x| crate::bessel_j(x[0], x[1]));
    }

    #[test]
    fn test_bessel_j_c64() {
        xsref::test("cyl_bessel_j", "d_cd-cd", |x| {
            crate::bessel_j(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_je_f64() {
        xsref::test("cyl_bessel_je", "d_d-d", |x| crate::bessel_je(x[0], x[1]));
    }

    #[test]
    fn test_bessel_je_c64() {
        xsref::test("cyl_bessel_je", "d_cd-cd", |x| {
            crate::bessel_je(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel Y

    #[test]
    fn test_bessel_y0() {
        xsref::test("cyl_bessel_y0", "d-d", |x| crate::bessel_y0(x[0]));
    }

    #[test]
    fn test_bessel_y1() {
        xsref::test("cyl_bessel_y1", "d-d", |x| crate::bessel_y1(x[0]));
    }

    #[test]
    fn test_bessel_y_f64() {
        xsref::test("cyl_bessel_y", "d_d-d", |x| crate::bessel_y(x[0], x[1]));
    }

    #[test]
    fn test_bessel_y_c64() {
        xsref::test("cyl_bessel_y", "d_cd-cd", |x| {
            crate::bessel_y(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ye_f64() {
        xsref::test("cyl_bessel_ye", "d_d-d", |x| crate::bessel_ye(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ye_c64() {
        xsref::test("cyl_bessel_ye", "d_cd-cd", |x| {
            crate::bessel_ye(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel I

    #[test]
    fn test_bessel_i0() {
        xsref::test("cyl_bessel_i0", "d-d", |x| crate::bessel_i0(x[0]));
    }

    #[test]
    fn test_bessel_i0e() {
        xsref::test("cyl_bessel_i0e", "d-d", |x| crate::bessel_i0e(x[0]));
    }

    #[test]
    fn test_bessel_i1() {
        xsref::test("cyl_bessel_i1", "d-d", |x| crate::bessel_i1(x[0]));
    }

    #[test]
    fn test_bessel_i1e() {
        xsref::test("cyl_bessel_i1e", "d-d", |x| crate::bessel_i1e(x[0]));
    }

    #[test]
    fn test_bessel_i_f64() {
        xsref::test("cyl_bessel_i", "d_d-d", |x| crate::bessel_i(x[0], x[1]));
    }

    #[test]
    fn test_bessel_i_c64() {
        xsref::test("cyl_bessel_i", "d_cd-cd", |x| {
            crate::bessel_i(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ie_f64() {
        xsref::test("cyl_bessel_ie", "d_d-d", |x| crate::bessel_ie(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ie_c64() {
        xsref::test("cyl_bessel_ie", "d_cd-cd", |x| {
            crate::bessel_ie(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel K

    #[test]
    fn test_bessel_k0() {
        xsref::test("cyl_bessel_k0", "d-d", |x| crate::bessel_k0(x[0]));
    }

    #[test]
    fn test_bessel_k0e() {
        xsref::test("cyl_bessel_k0e", "d-d", |x| crate::bessel_k0e(x[0]));
    }

    #[test]
    fn test_bessel_k1() {
        xsref::test("cyl_bessel_k1", "d-d", |x| crate::bessel_k1(x[0]));
    }

    #[test]
    fn test_bessel_k1e() {
        xsref::test("cyl_bessel_k1e", "d-d", |x| crate::bessel_k1e(x[0]));
    }

    #[test]
    fn test_bessel_k_f64() {
        xsref::test("cyl_bessel_k", "d_d-d", |x| crate::bessel_k(x[0], x[1]));
    }

    #[test]
    fn test_bessel_k_c64() {
        xsref::test("cyl_bessel_k", "d_cd-cd", |x| {
            crate::bessel_k(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_bessel_ke_f64() {
        xsref::test("cyl_bessel_ke", "d_d-d", |x| crate::bessel_ke(x[0], x[1]));
    }

    #[test]
    fn test_bessel_ke_c64() {
        xsref::test("cyl_bessel_ke", "d_cd-cd", |x| {
            crate::bessel_ke(x[0], c64(x[1], x[2]))
        });
    }

    // Hankel

    #[test]
    fn test_hankel_1_c64() {
        xsref::test("cyl_hankel_1", "d_cd-cd", |x| {
            crate::hankel_1(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_1e_c64() {
        xsref::test("cyl_hankel_1e", "d_cd-cd", |x| {
            crate::hankel_1e(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_2_c64() {
        xsref::test("cyl_hankel_2", "d_cd-cd", |x| {
            crate::hankel_2(x[0], c64(x[1], x[2]))
        });
    }

    #[test]
    fn test_hankel_2e_c64() {
        xsref::test("cyl_hankel_2e", "d_cd-cd", |x| {
            crate::hankel_2e(x[0], c64(x[1], x[2]))
        });
    }

    // Bessel zeros

    /// Based on `scipy.special.tests.test_basic.TestBessel.test_jnyn_zeros`
    #[test]
    fn test_bessel_zeros() {
        let jnz = crate::bessel_zeros::<5>(1);
        crate::np_assert_allclose!(
            jnz[0].as_ref(),
            &[3.83171, 7.01559, 10.17347, 13.32369, 16.47063],
            atol = 1.5e-5
        );
        crate::np_assert_allclose!(
            jnz[1].as_ref(),
            &[1.84118, 5.33144, 8.53632, 11.70600, 14.86359],
            atol = 1.5e-5
        );
        crate::np_assert_allclose!(
            jnz[2].as_ref(),
            &[2.19714, 5.42968, 8.59601, 11.74915, 14.89744],
            atol = 1.5e-5
        );
        crate::np_assert_allclose!(
            jnz[3].as_ref(),
            &[3.68302, 6.94150, 10.12340, 13.28576, 16.44006],
            atol = 1.5e-5
        );
    }

    // Bessel integrals

    #[test]
    fn test_it1j0y0() {
        xsref::test("it1j0y0", "d-d_d", |x| crate::it1j0y0(x[0]));
    }

    #[test]
    fn test_it2j0y0() {
        // This xsref table contains incorrect values for the 4th integral,
        // see https://github.com/scipy/xsref/issues/10.
        // So we instead use values that were verified using a Wolfram notebook.

        // xsref::test("it2j0y0", "d-d_d", |x| crate::it2j0y0(x[0]));

        let xs = [-1.0, 10.0, -10.0, 1.0, 0.2];
        let expect_0 = [
            0.121_165_246_995_068_71,
            2.177_866_420_093_336,
            2.177_866_420_093_336,
            0.121_165_246_995_068_71,
            0.004_993_754_627_460_185,
        ];
        let expect_1 = [
            f64::NAN,
            -0.022_987_933_564_659_65,
            f64::NAN,
            0.395_272_901_699_293_2,
            -0.434_230_670_112_316_34,
        ];
        np_assert_allclose!(xs.map(|x| crate::it2j0y0(x).0), expect_0);
        np_assert_allclose!(xs.map(|x| crate::it2j0y0(x).1), expect_1);
    }

    #[test]
    fn test_it1i0k0() {
        xsref::test("it1i0k0", "d-d_d", |x| crate::it1i0k0(x[0]));
    }

    #[test]
    fn test_it2i0k0() {
        xsref::test("it2i0k0", "d-d_d", |x| crate::it2i0k0(x[0]));
    }

    #[test]
    fn test_besselpoly() {
        xsref::test("besselpoly", "d_d_d-d", |x| {
            crate::besselpoly(x[0], x[1], x[2])
        });
    }

    // Riccati-Bessel

    /// Based on `scipy.special.tests.test_basic.TestRiccati.test_riccati_jn`
    #[test]
    fn test_riccati_j() {
        // N, x = 2, 0.2
        const N: usize = 2;
        let x: f64 = 0.2;
        // S = np.empty((N, N))
        let mut s = [[0.0; N], [0.0; N]];
        // for n in range(N):
        #[allow(clippy::needless_range_loop)]
        for n in 0..N {
            // j = special.spherical_jn(n, x)
            let j = crate::sph_bessel_j(n.try_into().unwrap(), x);
            // jp = special.spherical_jn(n, x, derivative=True)
            let jp = crate::sph_bessel_j_prime(n.try_into().unwrap(), x);
            // S[0,n] = x*j
            s[0][n] = x * j;
            // S[1,n] = x*jp + j
            s[1][n] = x * jp + j;
        }

        // assert_allclose(S, special.riccati_jn(n, x), atol=1.5e-8, rtol=0)
        let (jn, jnp) = crate::riccati_j::<N>(x);
        crate::np_assert_allclose!(&s[0], &jn, atol = 1.5e-8);
        crate::np_assert_allclose!(&s[1], &jnp, atol = 1.5e-8);
    }

    /// Based on `scipy.special.tests.test_basic.TestRiccati.test_riccati_yn`
    #[test]
    fn test_riccati_y() {
        // N, x = 2, 0.2
        const N: usize = 2;
        let x: f64 = 0.2;
        // S = np.empty((N, N))
        let mut c = [[0.0; N], [0.0; N]];
        // for n in range(N):
        #[allow(clippy::needless_range_loop)]
        for n in 0..N {
            // j = special.spherical_jn(n, x)
            let j = crate::sph_bessel_y(n.try_into().unwrap(), x);
            // jp = special.spherical_jn(n, x, derivative=True)
            let jp = crate::sph_bessel_y_prime(n.try_into().unwrap(), x);
            // S[0,n] = x*j
            c[0][n] = x * j;
            // S[1,n] = x*jp + j
            c[1][n] = x * jp + j;
        }

        // assert_allclose(S, special.riccati_jn(n, x), atol=1.5e-8, rtol=0)
        let (yn, ynp) = crate::riccati_y::<N>(x);
        crate::np_assert_allclose!(&c[0], &yn, atol = 1.5e-8);
        crate::np_assert_allclose!(&c[1], &ynp, atol = 1.5e-8);
    }

    #[test]
    fn test_riccati_n1() {
        // computed with N = 2, since `rctj` and `rcty` write at least 2 terms
        let (j, jp) = crate::riccati_j::<2>(0.2);
        assert_eq!(crate::riccati_j::<1>(0.2), ([j[0]], [jp[0]]));
        let (y, yp) = crate::riccati_y::<2>(0.2);
        assert_eq!(crate::riccati_y::<1>(0.2), ([y[0]], [yp[0]]));
    }

    // Based on `scipy.special.tests.test_basic.TestLambda.test_lmbda`
    #[test]
    fn test_jahnke_emden_lambda() {
        // lam = special.lmbda(1,.1)
        let lam = crate::jahnke_emden_lambda(1, 0.1);
        // lamr = (
        //     array([special.jn(0,.1), 2*special.jn(1,.1)/.1]),
        //     array([special.jvp(0,.1), -2*special.jv(1,.1)/.01 + 2*special.jvp(1,.1)/.1])
        // )
        let lamr = (
            [crate::bessel_j(0.0, 0.1), 2e1 * crate::bessel_j(1.0, 0.1)],
            [
                crate::bessel_j_prime(0.0, 0.1, 1),
                -2e2 * crate::bessel_j(1.0, 0.1) + 2e1 * crate::bessel_j_prime(1.0, 0.1, 1),
            ],
        );
        // assert_allclose(lam, lamr, atol=1.5e-8, rtol=0)
        crate::np_assert_allclose!(&lam.0, &lamr.0, atol = 1.5e-8);
        crate::np_assert_allclose!(&lam.1, &lamr.1, atol = 1.5e-8);
    }

    #[test]
    fn test_bessel_zeros_n0() {
        assert_eq!(crate::bessel_zeros::<0>(0), [[]; 4]);
    }

    #[test]
    #[should_panic(expected = "v is too large")]
    fn test_jahnke_emden_lambda_huge_v() {
        let _ = crate::jahnke_emden_lambda(2_147_483_648.5, 1.0);
    }

    #[test]
    fn test_jahnke_emden_lambda_nan() {
        // instead of out-of-bounds indexing in `lamn` and `lamv`
        for v in [2.0, 2.25] {
            for x in [f64::NAN, f64::INFINITY, -1e10] {
                let (vl, dl) = crate::jahnke_emden_lambda(v, x);
                assert_eq!((vl.len(), dl.len()), (3, 3));
                assert!(vl.iter().chain(&dl).all(|y| y.is_nan()), "v = {v}, x = {x}");
            }
        }
    }
}
