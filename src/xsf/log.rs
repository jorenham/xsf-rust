//! Translated into pure Rust from `xsf/log.h` and `xsf/cephes/unity.h` (xsf v0.2.2). The real
//! `log1p` uses `f64::ln_1p` from the Rust standard library instead, which (at least on Linux) is
//! more accurate than the Cephes one in xsf.

use core::f64::consts::LN_2;

use num_complex::Complex;

use crate::xsf::cephes::dd_real::DoubleDouble;

// cephes/const.h
const MAXITER: u64 = 500;
const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53

/// Natural logarithm of complex `z`
///
/// This replaces `std::log(std::complex<double>)` from C++. Unlike `num_complex::Complex::ln`, it
/// doesn't lose precision near the unit circle (the real part is computed in double-double
/// precision there) or for extreme magnitudes (the input is scaled by a power of two), like
/// glibc's `clog`. For non-finite or zero `z`, the return values are as in C99 Annex G.
#[allow(clippy::many_single_char_names)]
pub(crate) fn clog(z: Complex<f64>) -> Complex<f64> {
    let (x, y) = (z.re, z.im);
    let m = x.abs().max(y.abs());

    let re = if m < f64::MIN_POSITIVE {
        // scale up subnormal inputs by 2^54
        let s = f64::from_bits(0x4350_0000_0000_0000);
        (x * s).hypot(y * s).ln() - 54.0 * LN_2
    } else if m > f64::MAX / 2.0 {
        // scale down to avoid overflow in hypot
        (x * 0.5).hypot(y * 0.5).ln() + LN_2
    } else {
        let r = x.hypot(y);
        if (0.5..=2.0).contains(&r) {
            // ln(r) = log1p(x^2 + y^2 - 1) / 2, with x^2 + y^2 - 1 in double-double precision
            let (dx, dy) = (DoubleDouble::new(x), DoubleDouble::new(y));
            0.5 * f64::from(dx * dx + dy * dy + -1.0).ln_1p()
        } else {
            r.ln()
        }
    };
    Complex::new(re, y.atan2(x))
}

fn clog1p_ddouble(zr: f64, zi: f64) -> Complex<f64> {
    let r = DoubleDouble::new(zr);
    let i = DoubleDouble::new(zi);
    let two = DoubleDouble::new(2.0);

    let rsqr = r * r;
    let isqr = i * i;
    let rtwo = two * r;
    let mut absm1 = rsqr + isqr;
    absm1 = absm1 + rtwo;

    let x = 0.5 * f64::from(absm1).ln_1p();
    let y = zi.atan2(zr + 1.0);
    Complex::new(x, y)
}

// log(z + 1) = log(x + 1 + 1j*y)
//             = log(sqrt((x+1)**2 + y**2)) + 1j*atan2(y, x+1)
//
// Using atan2(y, x+1) for the imaginary part is always okay.  The real part
// needs to be calculated more carefully.  For |z| large, the naive formula
// log(z + 1) can be used.  When |z| is small, rewrite as
//
// log(sqrt((x+1)**2 + y**2)) = 0.5*log(x**2 + 2*x +1 + y**2)
//       = 0.5 * log1p(x**2 + y**2 + 2*x)
//       = 0.5 * log1p(hypot(x,y) * (hypot(x, y) + 2*x/hypot(x,y)))
//
// This expression suffers from cancellation when x < 0 and
// y = +/-sqrt(2*fabs(x)). To get around this cancellation problem, we use
// double-double precision when necessary.
fn clog1p(mut z: Complex<f64>) -> Complex<f64> {
    if !z.re.is_finite() || !z.im.is_finite() {
        z += 1.0;
        return clog(z);
    }

    let zr = z.re;
    let zi = z.im;

    if zi == 0.0 && zr >= -1.0 {
        return Complex::new(zr.ln_1p(), 0.0);
    }

    let az = z.norm();
    if az < 0.707 {
        let azi = zi.abs();
        if zr < 0.0 && (-zr - azi * azi / 2.0).abs() / (-zr) < 0.5 {
            return clog1p_ddouble(zr, zi);
        }
        let x = 0.5 * (az * (az + 2.0 * zr / az)).ln_1p();
        let y = zi.atan2(zr + 1.0);
        return Complex::new(x, y);
    }

    z += 1.0;
    clog(z)
}

pub trait LogArg: crate::sealed::Sealed {
    fn xsf_log1p(self) -> Self;
    fn xsf_xlogy(self, x: Self) -> Self;
    fn xsf_xlog1py(self, x: Self) -> Self;
}

impl LogArg for f64 {
    #[inline]
    fn xsf_log1p(self) -> Self {
        self.ln_1p()
    }

    #[inline]
    fn xsf_xlogy(self, x: Self) -> Self {
        let y = self;
        if x == 0.0 && !y.is_nan() {
            return 0.0;
        }

        x * y.ln()
    }

    #[inline]
    fn xsf_xlog1py(self, x: Self) -> Self {
        let y = self;
        if x == 0.0 && !y.is_nan() {
            return 0.0;
        }

        x * y.ln_1p()
    }
}

// NOTE(xsf-rust): the complex multiplication of num-complex doesn't do the C99 Annex G inf/NaN
// recovery of C/C++ (e.g. `__muldc3`), so e.g. `xlogy(1+1j, inf+NaNj)` is `NaN+NaNj` here,
// instead of `inf+infj`.
impl LogArg for Complex<f64> {
    #[inline]
    fn xsf_log1p(self) -> Self {
        clog1p(self)
    }

    #[inline]
    fn xsf_xlogy(self, x: Self) -> Self {
        let y = self;
        if x == Complex::ZERO && !y.is_nan() {
            return Complex::ZERO;
        }

        x * clog(y)
    }

    #[inline]
    fn xsf_xlog1py(self, x: Self) -> Self {
        let y = self;
        if x == Complex::ZERO && !y.is_nan() {
            return Complex::ZERO;
        }

        x * clog1p(y)
    }
}

/// $ \log(1+z)$ for real or complex input
///
/// Corresponds to [`scipy.special.log1p`][log1p].
///
/// [log1p]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log1p.html
///
/// # See also
/// - [`expm1`](crate::expm1)
/// - [`cosm1`](crate::cosm1)
/// - [`log1pmx`]
#[doc(alias = "ln_1p", alias = "log_1p")]
#[must_use]
#[inline]
pub fn log1p<T: LogArg>(z: T) -> T {
    z.xsf_log1p()
}

/// Compute $\log(1+x)-x$ for real input
///
/// Has no analogue in `scipy.special`.
///
/// # See also
/// - [`log1p`]
#[must_use]
#[inline]
#[allow(clippy::cast_precision_loss)]
pub fn log1pmx(x: f64) -> f64 {
    if x.abs() < 0.5 {
        let mut xfac = x;
        let mut res = 0.0;

        for n in 2..MAXITER {
            xfac *= -x;
            let term = xfac / n as f64;
            res += term;
            if term.abs() < MACHEP * res.abs() {
                break;
            }
        }
        res
    } else {
        x.ln_1p() - x
    }
}

/// Compute $x \log(y)$ for real or complex input
///
/// Corresponds to [`scipy.special.xlogy`][xlogy].
///
/// [xlogy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.xlogy.html
///
/// # See also
/// - [`xlog1py`]
#[doc(alias = "x_ln_y", alias = "x_log_y")]
#[must_use]
#[inline]
pub fn xlogy<T: LogArg>(x: T, y: T) -> T {
    y.xsf_xlogy(x)
}

/// Compute $x \log(1+y)$ for real or complex input
///
/// Corresponds to [`scipy.special.xlog1py`][xlog1py].
///
/// [xlog1py]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.xlog1py.html
///
/// # See also
/// - [`xlogy`]
#[doc(alias = "x_ln_1py", alias = "x_log_1py")]
#[must_use]
#[inline]
pub fn xlog1py<T: LogArg>(x: T, y: T) -> T {
    y.xsf_xlog1py(x)
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_log1p_f64() {
        xsref::test("log1p", "d-d", |x| crate::log1p(x[0]));
    }

    #[test]
    fn test_log1p_c64() {
        xsref::test("log1p", "cd-cd", |x| crate::log1p(c64(x[0], x[1])));
    }

    #[test]
    fn test_log1pmx_f64() {
        xsref::test("log1pmx", "d-d", |x| crate::log1pmx(x[0]));
    }

    #[test]
    fn test_xlogy_f64() {
        xsref::test("xlogy", "d_d-d", |x| crate::xlogy(x[0], x[1]));
    }

    #[test]
    fn test_xlogy_c64() {
        xsref::test("xlogy", "cd_cd-cd", |x| {
            crate::xlogy(c64(x[0], x[1]), c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_xlog1py_f64() {
        xsref::test("xlog1py", "d_d-d", |x| crate::xlog1py(x[0], x[1]));
    }

    #[test]
    fn test_xlog1py_c64() {
        xsref::test("xlog1py", "cd_cd-cd", |x| {
            crate::xlog1py(c64(x[0], x[1]), c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_clog() {
        use super::clog;

        // near the unit circle: |z|^2 = 1 + 2^-40 exactly, so ln|z| = 2^-41 - 2^-82 + ...
        let w = clog(c64(1.0, 2f64.powi(-20)));
        crate::np_assert_allclose!([w.re], [2f64.powi(-41) - 2f64.powi(-82)], rtol = 1e-15);
        // extreme magnitudes
        let w = clog(c64(f64::MAX, f64::MAX));
        crate::np_assert_allclose!([w.re], [710.129_286_483_663_9], rtol = 1e-15);
        let w = clog(c64(f64::from_bits(1), f64::from_bits(1)));
        crate::np_assert_allclose!([w.re], [-744.093_498_331_101_4], rtol = 1e-15);
    }

    #[test]
    fn test_log1p_c64_double_double() {
        // |1 + z|^2 = 1 + 2^-50 exactly
        let w = crate::log1p(c64(-2f64.powi(-25), 2f64.powi(-12)));
        crate::np_assert_allclose!([w.re], [0.5 * 2f64.powi(-50).ln_1p()], rtol = 1e-15);
    }

    #[test]
    fn test_xlogy_xlog1py_zero() {
        assert_eq!(crate::xlogy(0.0, f64::INFINITY), 0.0);
        assert_eq!(crate::xlog1py(0.0, -1.0), 0.0);
        assert_eq!(crate::xlogy(c64(0.0, 0.0), c64(0.0, 0.0)), c64(0.0, 0.0));
        assert!(crate::xlogy(0.0, f64::NAN).is_nan());
    }
}
