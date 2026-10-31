//! Translated into pure Rust from `xsf/ellint_carlson/ellint_argcheck.hh` (xsf main).

use core::num::FpCategory;

use num_complex::Complex64;

use super::typing::Num;

#[inline]
pub(super) fn invalid_rerr(r: f64, upper: f64) -> bool {
    r <= 0.0 || r > upper
}

#[inline]
pub(super) fn too_small(x: f64) -> bool {
    x == 0.0 || x.classify() == FpCategory::Subnormal
}

#[inline]
pub(super) fn too_small_cplx(z: Complex64) -> bool {
    too_small(z.re) && too_small(z.im)
}

/* Argument (phase) angle is NOT +/- pi, and is not indeterminate value,
 * checked without computing the actual argument. */
#[inline]
pub(super) fn ph_good(x: f64) -> bool {
    /* also handles x being NaN, which compares to false, hence is the
     * correct result. */
    x >= 0.0
}

#[inline]
pub(super) fn ph_good_cplx(z: Complex64) -> bool {
    match z.im.classify() {
        FpCategory::Normal | FpCategory::Subnormal => {
            let r = z.re;
            r.is_finite() || (r.is_infinite() && (r > 0.0))
        },
        FpCategory::Zero => ph_good(z.re),
        FpCategory::Nan => false,
        FpCategory::Infinite => z.re.is_finite(),
    }
}

/* "let two of the variables x, y, z be nonzero and conjugate complex with
 * phase less in magnitude than pi and the third variable be real and
 * nonnegative." */
#[inline]
pub(super) fn r1conj2<T: Num>(x: T, y: T, z: T) -> bool {
    x.isconj(y)
        && !(x.too_small() || y.too_small())
        && too_small(z.imag())
        && (z.real() >= 0.0)
        && x.ph_good()
        && y.ph_good()
}
