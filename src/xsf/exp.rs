//! `expm1`, `exp2`, and `exp10`: these use the implementations from the Rust standard library,
//! which (at least on Linux) are more accurate than the Cephes ones in xsf. Only the complex
//! `expm1` is translated into pure Rust from `xsf/exp.h` (xsf v0.2.2).

use num_complex::Complex;

use crate::xsf::trig::cosm1;

/// Exponential of complex `z`
///
/// This replaces `std::exp(std::complex<double>)` from C++. For finite `z`, it is computed like
/// glibc's `cexp`, which, unlike `num_complex::Complex::exp`, doesn't overflow prematurely when
/// only `exp(z.re)` overflows (instead of NaN, the imaginary part of `exp(x + 0i)` is `0`).
pub(crate) fn cexp(z: Complex<f64>) -> Complex<f64> {
    // (int) ((DBL_MAX_EXP - 1) * ln(2))
    const T: f64 = 709.0;

    if !z.re.is_finite() || !z.im.is_finite() {
        // See the note in cexpm1(z) below.
        return z.exp();
    }

    let (mut sinix, mut cosix) = if z.im.abs() > f64::MIN_POSITIVE {
        z.im.sin_cos()
    } else {
        (z.im, 1.0)
    };
    let mut x = z.re;
    if x > T {
        let exp_t = T.exp();
        x -= T;
        sinix *= exp_t;
        cosix *= exp_t;
        if x > T {
            x -= T;
            sinix *= exp_t;
            cosix *= exp_t;
        }
    }
    if x > T {
        // overflow (the original z.re > 3 T)
        Complex::new(f64::MAX * cosix, f64::MAX * sinix)
    } else {
        let exp_val = x.exp();
        Complex::new(exp_val * cosix, exp_val * sinix)
    }
}

// cexpm1(z) = cexp(z) - 1
//
// The imaginary part of this is easily computed via exp(z.real)*sin(z.imag)
// The real part is difficult to compute when there is cancellation e.g. when
// z.real = -log(cos(z.imag)).  There isn't a way around this problem  that
// doesn't involve computing exp(z.real) and/or cos(z.imag) to higher
// precision.
#[allow(clippy::needless_late_init)]
fn cexpm1(z: Complex<f64>) -> Complex<f64> {
    if !z.re.is_finite() || !z.im.is_finite() {
        // NOTE(xsf-rust): `num_complex::Complex::exp` follows libc++'s special cases, which differ
        // from glibc's `cexp` (used by libstdc++) only in the sign of the zero imaginary part of
        // exp(-inf + i(+-inf | NaN)).
        return z.exp() - 1.0;
    }

    let x;
    let mut ezr = 0.0;
    if z.re <= -40.0 {
        x = -1.0;
    } else {
        ezr = z.re.exp_m1();
        x = ezr * z.im.cos() + cosm1(z.im);
    }

    // don't compute exp(zr) too, unless necessary
    let y = if z.re > -1.0 {
        (ezr + 1.0) * z.im.sin()
    } else {
        z.re.exp() * z.im.sin()
    };

    Complex::new(x, y)
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
        cexpm1(self)
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
    use num_complex::c64;

    fn assert_close(actual: f64, expected: f64, rtol: f64) {
        assert!(
            (actual - expected).abs() <= rtol * expected.abs(),
            "actual: {actual:e}, expected: {expected:e}"
        );
    }

    #[test]
    fn test_expm1_f64() {
        xsref::test("expm1", "d-d", |x| crate::expm1(x[0]));
    }

    #[test]
    fn test_expm1_c64() {
        xsref::test("expm1", "cd-cd", |x| crate::expm1(c64(x[0], x[1])));
    }

    #[test]
    fn test_exp2_f64() {
        xsref::test("exp2", "d-d", |x| crate::exp2(x[0]));
    }

    #[test]
    fn test_exp10_f64() {
        xsref::test("exp10", "d-d", |x| crate::exp10(x[0]));
    }

    #[test]
    fn test_expm1_c64_edges() {
        // small |z|, where exp(z) - 1 would lose precision
        let w = crate::expm1(c64(1e-10, 1e-10));
        assert_close(w.re, 1e-10, 1e-15);
        assert_close(w.im, 1e-10 + 1e-20, 1e-15);

        // each branch: re <= -40, -40 < re <= -1, and re > -1
        for re in [-41.0_f64, -40.0, -1.5, -1.0, -0.25, 2.0] {
            let w = crate::expm1(c64(re, 2.0));
            assert_close(w.re, re.exp() * 2.0_f64.cos() - 1.0, 1e-14);
            assert_close(w.im, re.exp() * 2.0_f64.sin(), 1e-14);
        }
        assert_eq!(crate::expm1(c64(-40.0, 1.0)).re, -1.0);
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
