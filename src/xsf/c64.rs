//! A complex number with the arithmetic of `std::complex<double>` in C++.

use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_complex::Complex;

use crate::utils::cmul_exp;
use crate::xsf::complex::{catan, cdiv, cmul, csqrt, rdiv, rsub};
use crate::xsf::evalpoly::cevalpoly;
use crate::xsf::exp::cexp;
use crate::xsf::log::clog;
use crate::xsf::trig::{ccos, ccosh, csin, csinh};

/// A complex number with the arithmetic of `std::complex<double>` in C++ (GCC), so that ported
/// code can use the same operators:
///
/// - `z * w` is [`cmul`] and `z / w` is [`cdiv`], also for a real numerator (`x / w`);
/// - a real operand scales both parts in `z * x`, `x * z` and `z / x`, and only changes the real
///   part in `z + x`, `x + z` and `z - x`;
/// - `x - z` has the imaginary part `-z.im`.
///
/// The public API uses [`num_complex::Complex`], whose `*` and `/` are the naive formulas;
/// convert with [`ViaC64`] at the boundary.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct C64 {
    pub(crate) re: f64,
    pub(crate) im: f64,
}

/// `re + im i`
#[inline]
pub(crate) const fn c64(re: f64, im: f64) -> C64 {
    C64 { re, im }
}

impl From<Complex<f64>> for C64 {
    #[inline]
    fn from(z: Complex<f64>) -> Self {
        c64(z.re, z.im)
    }
}

impl From<f64> for C64 {
    #[inline]
    fn from(re: f64) -> Self {
        c64(re, 0.0)
    }
}

impl From<C64> for Complex<f64> {
    #[inline]
    fn from(z: C64) -> Self {
        Self::new(z.re, z.im)
    }
}

impl C64 {
    /// The imaginary unit
    pub(crate) const I: Self = c64(0.0, 1.0);
    /// NaN + NaN i
    pub(crate) const NAN: Self = c64(f64::NAN, f64::NAN);
    /// 1 + 0 i
    pub(crate) const ONE: Self = c64(1.0, 0.0);
    /// 0 + 0 i
    pub(crate) const ZERO: Self = c64(0.0, 0.0);

    /// `std::abs(z)`
    #[inline]
    pub(crate) fn norm(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// `std::atan(z)`
    #[inline]
    pub(crate) fn atan(self) -> Self {
        self.via(catan)
    }

    /// `std::conj(z)`
    #[inline]
    pub(crate) fn conj(self) -> Self {
        c64(self.re, -self.im)
    }

    /// `std::cos(z)`
    #[inline]
    pub(crate) fn cos(self) -> Self {
        self.via(ccos)
    }

    /// `std::cosh(z)`
    #[inline]
    pub(crate) fn cosh(self) -> Self {
        self.via(ccosh)
    }

    /// `xsf::cospi(z)`
    #[inline]
    pub(crate) fn cospi(self) -> Self {
        self.via(crate::cospi)
    }

    /// `xsf::cevalpoly(coeffs, degree, z)`
    #[inline]
    pub(crate) fn evalpoly(self, coeffs: &[f64]) -> Self {
        self.via(|z| cevalpoly(coeffs, z))
    }

    /// `std::exp(z)`
    #[inline]
    pub(crate) fn exp(self) -> Self {
        self.via(cexp)
    }

    /// `std::isfinite(z.real()) && std::isfinite(z.imag())`
    #[inline]
    pub(crate) fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }

    /// `std::isnan(z.real()) || std::isnan(z.imag())`
    #[inline]
    pub(crate) fn is_nan(self) -> bool {
        self.re.is_nan() || self.im.is_nan()
    }

    /// `std::log(z)`
    #[inline]
    pub(crate) fn ln(self) -> Self {
        self.via(clog)
    }

    /// `z * std::exp(x)` without overflowing prematurely, as in glibc's complex `exp`
    #[inline]
    pub(crate) fn mul_exp(self, x: f64) -> Self {
        self.via(|z| cmul_exp(z, x))
    }

    /// `std::sin(z)`
    #[inline]
    pub(crate) fn sin(self) -> Self {
        self.via(csin)
    }

    /// `std::sinh(z)`
    #[inline]
    pub(crate) fn sinh(self) -> Self {
        self.via(csinh)
    }

    /// `xsf::sinpi(z)`
    #[inline]
    pub(crate) fn sinpi(self) -> Self {
        self.via(crate::sinpi)
    }

    /// `std::sqrt(z)`
    #[inline]
    pub(crate) fn sqrt(self) -> Self {
        self.via(csqrt)
    }

    /// `f(z)` for a function `f` of [`num_complex::Complex`]
    #[inline]
    pub(crate) fn via(self, f: impl FnOnce(Complex<f64>) -> Complex<f64>) -> Self {
        f(self.into()).into()
    }

    /// `f(z, w)` for a function `f` of [`num_complex::Complex`]
    #[inline]
    fn via2(self, w: Self, f: fn(Complex<f64>, Complex<f64>) -> Complex<f64>) -> Self {
        f(self.into(), w.into()).into()
    }
}

/// The conversion at the public API
pub(crate) trait ViaC64 {
    /// `f(z)` for a function `f` of [`C64`]
    fn via(self, f: impl FnOnce(C64) -> C64) -> Self;
}

impl ViaC64 for Complex<f64> {
    #[inline]
    fn via(self, f: impl FnOnce(C64) -> C64) -> Self {
        f(self.into()).into()
    }
}

impl Neg for C64 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        c64(-self.re, -self.im)
    }
}

/// `C64 op C64`, `C64 op f64`, `f64 op C64`, and the two assignment forms of `op`
macro_rules! impl_op {
    (
        $Op:ident::$op:ident, $OpAssign:ident::$op_assign:ident,
        |$z:ident, $w:ident| $zw:expr, |$x:ident| $zx:expr, $xw:expr
    ) => {
        impl $Op for C64 {
            type Output = Self;

            #[inline]
            fn $op(self, $w: Self) -> Self {
                let $z = self;
                $zw
            }
        }

        impl $Op<f64> for C64 {
            type Output = Self;

            #[inline]
            fn $op(self, $x: f64) -> Self {
                let $z = self;
                $zx
            }
        }

        impl $Op<C64> for f64 {
            type Output = C64;

            #[inline]
            fn $op(self, $w: C64) -> C64 {
                let $x = self;
                $xw
            }
        }

        impl $OpAssign for C64 {
            #[inline]
            fn $op_assign(&mut self, w: Self) {
                *self = $Op::$op(*self, w);
            }
        }

        impl $OpAssign<f64> for C64 {
            #[inline]
            fn $op_assign(&mut self, x: f64) {
                *self = $Op::$op(*self, x);
            }
        }
    };
}

impl_op!(
    Add::add,
    AddAssign::add_assign,
    |z, w| c64(z.re + w.re, z.im + w.im),
    |x| c64(z.re + x, z.im),
    c64(x + w.re, w.im)
);
impl_op!(
    Sub::sub,
    SubAssign::sub_assign,
    |z, w| c64(z.re - w.re, z.im - w.im),
    |x| c64(z.re - x, z.im),
    w.via(|w| rsub(x, w))
);
impl_op!(
    Mul::mul,
    MulAssign::mul_assign,
    |z, w| z.via2(w, cmul),
    |x| c64(z.re * x, z.im * x),
    c64(x * w.re, x * w.im)
);
impl_op!(
    Div::div,
    DivAssign::div_assign,
    |z, w| z.via2(w, cdiv),
    |x| c64(z.re / x, z.im / x),
    w.via(|w| rdiv(x, w))
);

#[cfg(test)]
mod tests {
    use super::{C64, c64};
    use num_complex::Complex;

    fn bits(z: C64) -> [u64; 2] {
        [z.re.to_bits(), z.im.to_bits()]
    }

    #[test]
    fn test_real_operand() {
        // only the real part changes, and the zero keeps its sign
        assert_eq!(bits(c64(1.0, -0.0) + 2.0), bits(c64(3.0, -0.0)));
        assert_eq!(bits(2.0 + c64(1.0, -0.0)), bits(c64(3.0, -0.0)));
        assert_eq!(bits(c64(1.0, -0.0) - 2.0), bits(c64(-1.0, -0.0)));
        // x - z is x + (-z)
        assert_eq!(bits(2.0 - c64(1.0, 0.0)), bits(c64(1.0, -0.0)));

        // both parts are scaled, without the cross terms of a complex product or quotient
        let z = c64(1.0, f64::INFINITY);
        assert_eq!(bits(z * 0.5), bits(c64(0.5, f64::INFINITY)));
        assert_eq!(bits(0.5 * z), bits(c64(0.5, f64::INFINITY)));
        assert_eq!(bits(z / 2.0), bits(c64(0.5, f64::INFINITY)));
    }

    #[test]
    fn test_complex_operand() {
        assert_eq!(c64(1.0, 2.0) + c64(3.0, -5.0), c64(4.0, -3.0));
        assert_eq!(c64(1.0, 2.0) - c64(3.0, -5.0), c64(-2.0, 7.0));
        assert_eq!(bits(-c64(1.0, -0.0)), bits(c64(-1.0, 0.0)));
        assert_eq!(c64(1.0, 2.0) * c64(3.0, -5.0), c64(13.0, 1.0));
        assert_eq!(c64(13.0, 1.0) / c64(3.0, -5.0), c64(1.0, 2.0));
        assert_eq!(1.0 / c64(0.0, 2.0), c64(0.0, -0.5));

        // C99 Annex G: an infinity times a nonzero finite number is an infinity
        let inf = f64::INFINITY;
        assert_eq!(c64(inf, inf) * c64(1.0, 0.0), c64(inf, inf));
        // no overflow of |w|^2
        assert_eq!(c64(1e200, 1e200) / c64(1e200, 1e200), c64(1.0, 0.0));
    }

    #[test]
    fn test_assign() {
        let mut z = c64(1.0, 2.0);
        z += c64(1.0, 1.0);
        z -= 1.0;
        z *= c64(0.0, 1.0);
        z /= 2.0;
        assert_eq!(z, c64(-1.5, 0.5));
    }

    #[test]
    fn test_conversion() {
        // the bits are kept, also of a NaN and of the zeros
        let nan = f64::from_bits(0x7ff8_0000_0000_1234);
        let z = Complex::new(-0.0, nan);
        let w = Complex::from(C64::from(z));
        assert_eq!(
            [w.re.to_bits(), w.im.to_bits()],
            [z.re.to_bits(), nan.to_bits()]
        );
    }
}
