//! Translated into pure Rust from `xsf/ellint_carlson/ellint_typing.hh` (xsf main).
//!
//! The C++ templates are generic over the real or complex floating-point type `T`, with
//! `RT = typing::decplx_t<T>` the underlying real type. Here `T` is generic over [`Num`], which
//! is implemented for `f64` and `Complex64`, and `RT` is always `f64`. The C++ overloads that are
//! selected by template matching on real or complex `T` (`typing::real_only` and
//! `typing::cplx_only`) are the methods of [`Num`].

use core::ops::{Add, Div, Mul, Neg, Sub};

use num_complex::Complex64;

use super::{argcheck, arithmetic, rg, rj};
use crate::xsf::complex::{cdiv, cmul, csqrt};

/// The arithmetic of `T` in C++, i.e. `double` or `std::complex<double>`
///
/// Note that `T + RT` and `T * RT` are component-wise in C++, which differs from `T + T(RT)` if
/// the imaginary part is `-0.0`.
pub(crate) trait Num:
    Copy
    + PartialEq
    + From<f64>
    + Add<Output = Self>
    + Add<f64, Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    /// `typing::nan<T>()`
    const NAN: Self;
    /// `typing::huge<T>()`
    const HUGE: Self;

    /// `T * T`
    fn prod(self, rhs: Self) -> Self;
    /// `T / T`
    fn quot(self, rhs: Self) -> Self;
    /// `std::sqrt(x)`
    fn sqrt(self) -> Self;
    /// `std::abs(x)`
    fn abs(self) -> f64;
    /// `std::real(x)`
    fn real(self) -> f64;
    /// `std::imag(x)`
    fn imag(self) -> f64;

    /// `argcheck::too_small(x)`
    fn too_small(self) -> bool;
    /// `argcheck::ph_good(x)`
    fn ph_good(self) -> bool;
    /// `argcheck::isinf(x)`
    fn isinf(self) -> bool;
    /// `argcheck::isconj(x, y)`
    fn isconj(self, y: Self) -> bool;

    /// `arithmetic::ndot2(x, y, n)`
    fn ndot2(x: &[Self], y: &[Self], n: usize) -> Self;
    /// `arithmetic::dcomp_horner(x, poly, degree)`
    fn dcomp_horner(self, poly: &[f64], degree: usize) -> Self;
    /// `arithmetic::aux::rg_dot2_acc(fac, term, acc, cor)`
    fn rg_dot2_acc(fac: f64, term: Self, acc: &mut Self, cor: &mut Self);
    /// `rjimpl::safe_atan_sqrt_div(x)`
    fn safe_atan_sqrt_div(self) -> Self;
}

impl Num for f64 {
    const NAN: Self = f64::NAN;
    const HUGE: Self = f64::INFINITY;

    #[inline]
    fn prod(self, rhs: Self) -> Self {
        self * rhs
    }
    #[inline]
    fn quot(self, rhs: Self) -> Self {
        self / rhs
    }
    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    #[inline]
    fn abs(self) -> f64 {
        f64::abs(self)
    }
    #[inline]
    fn real(self) -> f64 {
        self
    }
    #[inline]
    fn imag(self) -> f64 {
        0.0
    }

    #[inline]
    fn too_small(self) -> bool {
        argcheck::too_small(self)
    }
    #[inline]
    fn ph_good(self) -> bool {
        argcheck::ph_good(self)
    }
    #[inline]
    fn isinf(self) -> bool {
        self.is_infinite()
    }
    #[inline]
    fn isconj(self, y: Self) -> bool {
        // `std::conj(y)` is `std::complex<double>(y, -0.0)` for real `y`, so `x - std::conj(y)`
        // has a zero imaginary part
        argcheck::too_small(self - y)
    }

    #[inline]
    fn ndot2(x: &[Self], y: &[Self], n: usize) -> Self {
        arithmetic::ndot2_real(x, y, n)
    }
    #[inline]
    fn dcomp_horner(self, poly: &[f64], degree: usize) -> Self {
        arithmetic::dcomp_horner_real(self, poly, degree)
    }
    #[inline]
    fn rg_dot2_acc(fac: f64, term: Self, acc: &mut Self, cor: &mut Self) {
        arithmetic::fdot2_acc(fac, term, acc, cor);
    }
    #[inline]
    fn safe_atan_sqrt_div(self) -> Self {
        rj::safe_atan_sqrt_div(self)
    }
}

impl Num for Complex64 {
    /* Complexified specialisation, with the complex NaN and infinity as
     * return values. */
    const NAN: Self = Complex64::new(f64::NAN, f64::NAN);
    const HUGE: Self = Complex64::new(f64::INFINITY, 0.0);

    #[inline]
    fn prod(self, rhs: Self) -> Self {
        cmul(self, rhs)
    }
    #[inline]
    fn quot(self, rhs: Self) -> Self {
        cdiv(self, rhs)
    }
    #[inline]
    fn sqrt(self) -> Self {
        csqrt(self)
    }
    #[inline]
    fn abs(self) -> f64 {
        self.norm()
    }
    #[inline]
    fn real(self) -> f64 {
        self.re
    }
    #[inline]
    fn imag(self) -> f64 {
        self.im
    }

    #[inline]
    fn too_small(self) -> bool {
        argcheck::too_small_cplx(self)
    }
    #[inline]
    fn ph_good(self) -> bool {
        argcheck::ph_good_cplx(self)
    }
    #[inline]
    fn isinf(self) -> bool {
        self.re.is_infinite() || self.im.is_infinite()
    }
    #[inline]
    fn isconj(self, y: Self) -> bool {
        argcheck::too_small_cplx(self - y.conj())
    }

    #[inline]
    fn ndot2(x: &[Self], y: &[Self], n: usize) -> Self {
        arithmetic::ndot2_cplx(x, y, n)
    }
    #[inline]
    fn dcomp_horner(self, poly: &[f64], degree: usize) -> Self {
        arithmetic::dcomp_horner_cplx(self, poly, degree)
    }
    #[inline]
    fn rg_dot2_acc(fac: f64, term: Self, acc: &mut Self, cor: &mut Self) {
        rg::rg_dot2_acc_cplx(fac, term, acc, cor);
    }
    #[inline]
    fn safe_atan_sqrt_div(self) -> Self {
        rj::safe_atan_sqrt_div_cplx(self)
    }
}
