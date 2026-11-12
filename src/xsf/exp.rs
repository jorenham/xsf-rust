//! `expm1`, `exp2`, and `exp10`: these use the implementations from the Rust standard library,
//! which (at least on Linux) are more accurate than the Cephes ones in xsf. Only the complex
//! `expm1` is translated into pure Rust from `xsf/exp.h` (scipy/xsf@b1eb270).

use num_complex::Complex;

use crate::utils::{cmul_exp, sin_cos_glibc};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::trig::cosm1;

/// Exponential of complex `z`
///
/// This replaces `std::exp(std::complex<double>)` from C++. For finite `z`, it is computed like
/// glibc's `cexp`, which, unlike `num_complex::Complex::exp`, doesn't overflow prematurely when
/// only `exp(z.re)` overflows (instead of NaN, the imaginary part of `exp(x + 0i)` is `0`). The
/// sign of the zero imaginary part of `exp(-inf + i(+-inf | NaN))` also follows glibc.
pub(crate) fn cexp(z: C64) -> C64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        if z.re == f64::NEG_INFINITY && !z.im.is_finite() {
            // like glibc (`num_complex` returns `0 + 0i` here)
            return c64(0.0, 0.0_f64.copysign(z.im));
        }
        // See the note in cexpm1(z) below.
        return Complex::from(z).exp().into();
    }

    let (sinix, cosix) = sin_cos_glibc(z.im);
    cmul_exp(c64(cosix, sinix), z.re)
}

// cexpm1(z) = cexp(z) - 1
//
// The imaginary part of this is easily computed via exp(z.real)*sin(z.imag)
// The real part is difficult to compute when there is cancellation e.g. when
// z.real = -log(cos(z.imag)).  There isn't a way around this problem  that
// doesn't involve computing exp(z.real) and/or cos(z.imag) to higher
// precision.
fn cexpm1(z: C64) -> C64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        // NOTE(xsf-rust): `num_complex::Complex::exp` follows libc++'s special cases, which differ
        // from glibc's `cexp` (used by libstdc++) only in the sign of the zero imaginary part of
        // exp(-inf + i(+-inf | NaN)).
        return (Complex::from(z).exp() - 1.0).into();
    }

    let mut ezr = 0.0;
    let x = if z.re <= -40.0 {
        -1.0
    } else {
        ezr = z.re.exp_m1();
        ezr * z.im.cos() + cosm1(z.im)
    };

    // don't compute exp(zr) too, unless necessary
    let y = if z.re > -1.0 {
        (ezr + 1.0) * z.im.sin()
    } else {
        z.re.exp() * z.im.sin()
    };

    c64(x, y)
}

pub trait ExpArg: crate::sealed::Sealed {
    fn expm1(self) -> Self;
}

impl ExpArg for f64 {
    #[inline]
    fn expm1(self) -> Self {
        self.exp_m1()
    }
}

impl ExpArg for Complex<f64> {
    #[inline]
    fn expm1(self) -> Self {
        self.via(cexpm1)
    }
}

/// $e^x - 1$ for real or complex input
///
/// Corresponds to [`scipy.special.expm1`][expm1]
///
/// [expm1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expm1.html
///
/// # See also
/// - [`exp2`] for $2^x$
/// - [`exp10`] for $10^x$
#[doc(alias = "exp_m1")]
#[must_use]
#[inline]
pub fn expm1<T: ExpArg>(z: T) -> T {
    z.expm1()
}

/// $2^x$
///
/// Corresponds to [`scipy.special.exp2`][exp2]
///
/// [exp2]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp2.html
///
/// # See also
/// - [`expm1`] for $e^x - 1$
/// - [`exp10`] for $10^x$
#[must_use]
#[inline]
pub fn exp2(x: f64) -> f64 {
    x.exp2()
}

/// $10^x$
///
/// Corresponds to [`scipy.special.exp10`][exp10]
///
/// [exp10]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp10.html
///
/// # See also
/// - [`expm1`] for $e^x - 1$
/// - [`exp2`] for $2^x$
#[must_use]
#[inline]
pub fn exp10(x: f64) -> f64 {
    10.0_f64.powf(x)
}

#[cfg(test)]
mod tests {
    use crate::xsf::c64::c64;

    fn assert_close(actual: f64, expected: f64, rtol: f64) {
        assert!(
            (actual - expected).abs() <= rtol * expected.abs(),
            "actual: {actual:e}, expected: {expected:e}"
        );
    }

    #[test]
    fn test_cexp_overflow() {
        use crate::xsf::exp::cexp;

        // exp(z.re) overflows, but not all of exp(z) does (reference values from mpmath)
        let w = cexp(c64(709.9, 1.5));
        assert_close(w.re, 1.429_883_248_952_389_1e307, 1e-15);
        assert_eq!(w.im, f64::INFINITY);
        let w = cexp(c64(1450.0, 5e-324));
        assert_eq!(w.re, f64::INFINITY);
        assert_close(w.im, 2.635_016_970_633_634_3e306, 1e-15);

        // real z: inf + 0i instead of inf + NaN i (once or twice scaled, and z.re > 3 * 709)
        for re in [800.0, 1500.0, 3000.0] {
            let w = cexp(c64(re, -0.0));
            assert_eq!(w.re, f64::INFINITY);
            assert_eq!(w.im.to_bits(), (-0.0_f64).to_bits());
        }
    }
}
