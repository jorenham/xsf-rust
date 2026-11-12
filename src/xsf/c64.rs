//! A complex number with the arithmetic of `std::complex<double>` in C++.

use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_complex::Complex;

use crate::utils::{cmul_exp, pow2};
use crate::xsf::complex::{catan, cpow, csqrt};
use crate::xsf::evalpoly::cevalpoly_c64;
use crate::xsf::exp::cexp;
use crate::xsf::log::clog;
use crate::xsf::trig::{ccos, ccosh, ccospi, csin, csinh, csinpi};

/// A complex number with the arithmetic of `std::complex<double>` in C++ (GCC), so that ported
/// code can use the same operators:
///
/// - `z * w` is [`mul`] and `z / w` is [`div`], also for a real numerator (`x / w`);
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
        catan(self)
    }

    /// `std::conj(z)`
    #[inline]
    pub(crate) fn conj(self) -> Self {
        c64(self.re, -self.im)
    }

    /// `std::cos(z)`
    #[inline]
    pub(crate) fn cos(self) -> Self {
        ccos(self)
    }

    /// `std::cosh(z)`
    #[inline]
    pub(crate) fn cosh(self) -> Self {
        ccosh(self)
    }

    /// `xsf::cospi(z)`
    #[inline]
    pub(crate) fn cospi(self) -> Self {
        ccospi(self)
    }

    /// `xsf::cevalpoly(coeffs, degree, z)`
    #[inline]
    pub(crate) fn evalpoly(self, coeffs: &[f64]) -> Self {
        cevalpoly_c64(coeffs, self)
    }

    /// `std::exp(z)`
    #[inline]
    pub(crate) fn exp(self) -> Self {
        cexp(self)
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

    /// `std::abs(z.real()) + std::abs(z.imag())`
    #[inline]
    pub(crate) fn l1_norm(self) -> f64 {
        self.re.abs() + self.im.abs()
    }

    /// `std::log(z)`
    #[inline]
    pub(crate) fn ln(self) -> Self {
        clog(self)
    }

    /// `z * std::exp(x)` without overflowing prematurely, as in glibc's complex `exp`
    #[inline]
    pub(crate) fn mul_exp(self, x: f64) -> Self {
        cmul_exp(self, x)
    }

    /// `z.real() * z.real() + z.imag() * z.imag()`
    #[inline]
    pub(crate) fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// `std::pow(z, x)`
    #[inline]
    pub(crate) fn powf(self, x: f64) -> Self {
        cpow(self, x)
    }

    /// `std::sin(z)`
    #[inline]
    pub(crate) fn sin(self) -> Self {
        csin(self)
    }

    /// `std::sinh(z)`
    #[inline]
    pub(crate) fn sinh(self) -> Self {
        csinh(self)
    }

    /// `xsf::sinpi(z)`
    #[inline]
    pub(crate) fn sinpi(self) -> Self {
        csinpi(self)
    }

    /// `std::sqrt(z)`
    #[inline]
    pub(crate) fn sqrt(self) -> Self {
        csqrt(self)
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
    // as in C++, the imaginary part is `-w.im` instead of `0.0 - w.im`
    c64(x - w.re, -w.im)
);
impl_op!(
    Mul::mul,
    MulAssign::mul_assign,
    |z, w| mul(z, w),
    |x| c64(z.re * x, z.im * x),
    c64(x * w.re, x * w.im)
);
impl_op!(
    Div::div,
    DivAssign::div_assign,
    |z, w| div(z, w),
    |x| c64(z.re / x, z.im / x),
    // as in C++, i.e. `(x + 0i) / w`
    div(c64(x, 0.0), w)
);

/// Complex division `a / b` using Smith's algorithm, which avoids the overflow and underflow of
/// `|b|^2` in the naive formula. Tiny operands are first scaled by the same power of 2, so that
/// the intermediate results don't lose precision by becoming subnormal, and finite operands are
/// halved if a sum of the algorithm overflows.
///
/// This replaces the complex division in C++, i.e. GCC's `__divdc3`. Like `__divdc3`, it recovers
/// the infinities and zeros that would otherwise be computed as NaN + NaN i (C99 Annex G).
///
/// R. L. Smith, "Algorithm 116: Complex division", Communications of the ACM 5(8), 1962.
fn div(mut a: C64, mut b: C64) -> C64 {
    const TWO_M400: f64 = pow2(-400);
    const TWO_600: f64 = pow2(600);

    if b.re.abs().max(b.im.abs()) < TWO_M400 {
        a *= TWO_600;
        b *= TWO_600;
    }
    // the denominator and the two numerators
    let smith = |a: C64, b: C64| {
        if b.re.abs() >= b.im.abs() {
            let r = b.im / b.re;
            [b.re + b.im * r, a.re + a.im * r, a.im - a.re * r]
        } else {
            let r = b.re / b.im;
            [b.re * r + b.im, a.re * r + a.im, a.im * r - a.re]
        }
    };
    let mut sums = smith(a, b);
    if sums.iter().any(|s| s.is_infinite()) && a.is_finite() && b.is_finite() {
        // halved, the sums are at most twice their largest term
        sums = smith(0.5 * a, 0.5 * b);
    }
    let [den, x, y] = sums;
    let w = c64(x / den, y / den);

    // Recover infinities and zeros that computed as NaN+iNaN; the only
    // cases are nonzero/zero, infinite/finite, and finite/infinite.
    if w.re.is_nan() && w.im.is_nan() {
        let (C64 { re: a, im: b }, C64 { re: c, im: d }) = (a, b);
        if c == 0.0 && d == 0.0 && (!a.is_nan() || !b.is_nan()) {
            let inf = f64::INFINITY.copysign(c);
            return c64(inf * a, inf * b);
        } else if (a.is_infinite() || b.is_infinite()) && c.is_finite() && d.is_finite() {
            let a = box_inf(a);
            let b = box_inf(b);
            return c64(
                f64::INFINITY * (a * c + b * d),
                f64::INFINITY * (b * c - a * d),
            );
        } else if (c.is_infinite() || d.is_infinite()) && a.is_finite() && b.is_finite() {
            let c = box_inf(c);
            let d = box_inf(d);
            return c64(0.0 * (a * c + b * d), 0.0 * (b * c - a * d));
        }
    }
    w
}

/// Complex multiplication `a * b`, like GCC's `__muldc3` (C99 Annex G)
///
/// This replaces the complex multiplication in C++. For finite results it is the same as the naive
/// formula (i.e. `num_complex`'s `*`), but it recovers the infinities that would otherwise be
/// computed as NaN + NaN i, e.g. `(inf + inf i) * (1 + 0i)`.
#[inline]
fn mul(a: C64, b: C64) -> C64 {
    let (
        C64 {
            re: mut a,
            im: mut b,
        },
        C64 {
            re: mut c,
            im: mut d,
        },
    ) = (a, b);
    let (ac, bd, ad, bc) = (a * c, b * d, a * d, b * c);
    let (mut x, mut y) = (ac - bd, ad + bc);

    if x.is_nan() && y.is_nan() {
        // Recover infinities that computed as NaN+iNaN.
        let mut recalc = false;
        if a.is_infinite() || b.is_infinite() {
            // z is infinite.  "Box" the infinity and change NaNs in
            // the other factor to 0.
            a = box_inf(a);
            b = box_inf(b);
            c = nan_to_zero(c);
            d = nan_to_zero(d);
            recalc = true;
        }
        if c.is_infinite() || d.is_infinite() {
            // w is infinite.  "Box" the infinity and change NaNs in
            // the other factor to 0.
            c = box_inf(c);
            d = box_inf(d);
            a = nan_to_zero(a);
            b = nan_to_zero(b);
            recalc = true;
        }
        if !recalc && (ac.is_infinite() || bd.is_infinite() || ad.is_infinite() || bc.is_infinite())
        {
            // Recover infinities from overflow by changing NaNs to 0.
            a = nan_to_zero(a);
            b = nan_to_zero(b);
            c = nan_to_zero(c);
            d = nan_to_zero(d);
            recalc = true;
        }
        if recalc {
            x = f64::INFINITY * (a * c - b * d);
            y = f64::INFINITY * (a * d + b * c);
        }
    }
    c64(x, y)
}

/// `isnan(x) ? copysign(0, x) : x`, used by [`mul`]
fn nan_to_zero(x: f64) -> f64 {
    if x.is_nan() { 0.0_f64.copysign(x) } else { x }
}

/// `copysign(isinf(x) ? 1 : 0, x)`, used to "box" infinities in C99 Annex G
fn box_inf(x: f64) -> f64 {
    (if x.is_infinite() { 1.0_f64 } else { 0.0 }).copysign(x)
}

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
    fn test_div_overflow() {
        // exact, where a sum of Smith's algorithm overflows
        let max = f64::MAX;
        let z = c64(max, max);
        assert_eq!(z / z, c64(1.0, 0.0));
        assert_eq!(z / c64(2.0, 2.0), c64(0.5 * max, 0.0));
        assert_eq!(c64(1.0, 0.0) / z, c64(0.5 / max, -0.5 / max));

        // not halved otherwise, which would lose the sign of the zero
        let w = c64(-1.0, 0.0) / c64(max, -(2.0_f64.powi(-50)));
        assert_eq!(w.im.to_bits(), (-0.0_f64).to_bits());
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
