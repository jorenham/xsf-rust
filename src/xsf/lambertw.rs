//! Translated into pure Rust from `xsf/lambertw.h` (xsf v0.2.2).
//!
//! Translated from Cython into C++ by SciPy developers in 2023.
//! Original header with Copyright information appears below.
//!
//! Implementation of the Lambert W function \[1\]. Based on MPMath
//! Implementation \[2\], and documentation \[3\].
//!
//! Copyright: Yosef Meller, 2009
//! Author email: mellerf@netvision.net.il
//!
//! Distributed under the same license as SciPy
//!
//! References:
//! \[1\] On the Lambert W function, Adv. Comp. Math. 5 (1996) 329-359,
//!     available online: <https://web.archive.org/web/20230123211413/https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf>
//! \[2\] mpmath source code,
//!     <https://github.com/mpmath/mpmath/blob/c5939823669e1bcce151d89261b802fe0d8978b4/mpmath/functions/functions.py#L435-L461>
//! \[3\] <https://web.archive.org/web/20230504171447/https://mpmath.org/doc/current/functions/powers.html#lambert-w-function>
//!
//! TODO: use a series expansion when extremely close to the branch point
//! at `-1/e` and make sure that the proper branch is chosen there.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted. Complex division
//! uses Smith's algorithm (see [`cdiv`]), which, like GCC's `__divdc3`, avoids the spurious
//! underflow of the naive formula for tiny `z` with `k != 0`.

use crate::xsf::evalpoly::cevalpoly;
use crate::xsf::log::clog;

use core::f64::consts::{E, PI};
use num_complex::Complex;

const EXPN1: f64 = 0.367_879_441_171_442_33; // exp(-1)
const OMEGA: f64 = 0.567_143_290_409_783_8; // W(1, 0)

/// Complex division `a / b` using Smith's algorithm, which avoids the overflow and underflow of
/// `|b|^2` in the naive formula. Tiny operands are first scaled by the same power of 2, so that
/// the intermediate results don't lose precision by becoming subnormal.
///
/// R. L. Smith, "Algorithm 116: Complex division", Communications of the ACM 5(8), 1962.
#[inline]
fn cdiv(mut a: Complex<f64>, mut b: Complex<f64>) -> Complex<f64> {
    const TWO_M400: f64 = f64::from_bits((1023 - 400) << 52); // 2^-400
    const TWO_600: f64 = f64::from_bits((1023 + 600) << 52); // 2^600

    if b.re.abs().max(b.im.abs()) < TWO_M400 {
        a *= TWO_600;
        b *= TWO_600;
    }
    if b.re.abs() >= b.im.abs() {
        let r = b.im / b.re;
        let den = b.re + b.im * r;
        Complex::new((a.re + a.im * r) / den, (a.im - a.re * r) / den)
    } else {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        Complex::new((a.re * r + a.im) / den, (a.im * r - a.re) / den)
    }
}

fn lambertw_branchpt(z: Complex<f64>) -> Complex<f64> {
    // Series for W(z, 0) around the branch point; see 4.22 in [1].
    let coeffs = [-1.0 / 3.0, 1.0, -1.0];
    let p = (2.0 * (E * z + 1.0)).sqrt();

    cevalpoly(&coeffs, p)
}

fn lambertw_pade0(z: Complex<f64>) -> Complex<f64> {
    // (3, 2) Pade approximation for W(z, 0) around 0.
    let num = [12.851_063_829_787_234, 12.340_425_531_914_894, 1.0];
    let denom = [32.531_914_893_617_02, 14.340_425_531_914_894, 1.0];

    /* This only gets evaluated close to 0, so we don't need a more
     * careful algorithm that avoids overflow in the numerator for
     * large z. */
    cdiv(z * cevalpoly(&num, z), cevalpoly(&denom, z))
}

#[allow(clippy::cast_precision_loss)]
fn lambertw_asy(z: Complex<f64>, k: isize) -> Complex<f64> {
    /* Compute the W function using the first two terms of the
     * asymptotic series. See 4.20 in [1].
     */
    let w = clog(z) + 2.0 * PI * k as f64 * Complex::I;
    w - clog(w)
}

/// Lambert W function.
///
/// The Lambert W function `W(z)` is defined as the inverse function of `w * exp(w)`. In other
/// words, the value of `W(z)` is such that `z = W(z) * exp(W(z))` for any complex number `z`.
///
/// The Lambert W function is a multivalued function with infinitely many branches. Each branch
/// gives a separate solution of the equation `z = w exp(w)`. Here, the branches are indexed by the
/// integer `k`.
///
/// Corresponds to [`scipy.special.lambertw`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lambertw.html
#[doc(alias = "lambert_w")]
#[must_use]
#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
pub fn lambertw(z: Complex<f64>, k: isize, tol: f64) -> Complex<f64> {
    if z.is_nan() {
        return z;
    }
    if z.re == f64::INFINITY {
        return z + 2.0 * PI * k as f64 * Complex::I;
    }
    if z.re == f64::NEG_INFINITY {
        return -z + (2.0 * PI * k as f64 + PI) * Complex::I;
    }
    if z == Complex::ZERO {
        if k == 0 {
            return z;
        }
        return Complex::new(f64::NEG_INFINITY, 0.0);
    }
    if z == Complex::ONE && k == 0 {
        // Split out this case because the asymptotic series blows up
        return Complex::new(OMEGA, 0.0);
    }

    let absz = z.norm();
    // Get an initial guess for Halley's method
    let mut w = if k == 0 {
        if (z + EXPN1).norm() < 0.3 {
            lambertw_branchpt(z)
        } else if -1.0 < z.re && z.re < 1.5 && z.im.abs() < 1.0 && -2.5 * z.im.abs() - 0.2 < z.re {
            /* Empirically determined decision boundary where the Pade
             * approximation is more accurate. */
            lambertw_pade0(z)
        } else {
            lambertw_asy(z, k)
        }
    } else if k == -1 {
        if absz <= EXPN1 && z.im == 0.0 && z.re < 0.0 {
            Complex::new((-z.re).ln(), 0.0)
        } else {
            lambertw_asy(z, k)
        }
    } else {
        lambertw_asy(z, k)
    };

    // Halley's method; see 5.9 in [1]
    if w.re >= 0.0 {
        // Rearrange the formula to avoid overflow in exp
        for _ in 0..100 {
            let ew = (-w).exp();
            let wewz = w - z * ew;
            let wn = w - cdiv(wewz, w + 1.0 - cdiv((w + 2.0) * wewz, 2.0 * w + 2.0));
            if (wn - w).norm() <= tol * wn.norm() {
                return wn;
            }
            w = wn;
        }
    } else {
        for _ in 0..100 {
            let ew = w.exp();
            let wew = w * ew;
            let wewz = wew - z;
            let wn = w - cdiv(wewz, wew + ew - cdiv((w + 2.0) * wewz, 2.0 * w + 2.0));
            if (wn - w).norm() <= tol * wn.norm() {
                return wn;
            }
            w = wn;
        }
    }

    Complex::new(f64::NAN, f64::NAN)
}

#[cfg(test)]
mod tests {
    use num_complex::c64;
    use num_traits::ToPrimitive;

    #[test]
    fn test_lambertw_c64() {
        xsref::test("lambertw", "cd_p_d-cd", |x| {
            crate::lambertw(c64(x[0], x[1]), x[2].to_isize().unwrap(), x[3])
        });
    }
}
