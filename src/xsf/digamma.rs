//! Translated into pure Rust from `xsf/digamma.h` (xsf v0.2.2).
//!
//! Unlike in xsf (scipy/xsf#347), the number of forward recurrence steps for complex `z` with
//! `z.re < 0` is clamped to `n >= 0`, e.g. `digamma(-1e-300 + 16i)` was `digamma(z + 1)`.

use core::f64::consts::PI;

use num_complex::Complex;

use crate::xsf::complex::{cdiv, rsub};
use crate::xsf::log::clog;

// All of the following were computed with mpmath
/// Location of the positive root
const DIGAMMA_POSROOT: f64 = 1.461_632_144_968_362_2;
/// Value of the positive root
const DIGAMMA_POSROOTVAL: f64 = -9.241_265_521_729_427e-17;
/// Location of the negative root
const DIGAMMA_NEGROOT: f64 = -0.504_083_008_264_455_4;
/// Value of the negative root
const DIGAMMA_NEGROOTVAL: f64 = 7.289_763_902_976_895e-17;

fn digamma_zeta_series(z: f64, root: f64, rootval: f64) -> f64 {
    let mut res = rootval;
    let mut coeff = -1.0;

    let z = z - root;
    for n in 1..100 {
        coeff *= -z;
        let term = coeff * crate::xsf::cephes::zeta(f64::from(n + 1), root);
        res += term;
        if term.abs() < f64::EPSILON * res.abs() {
            break;
        }
    }
    res
}

/// `digamma_zeta_series` for complex `z`
fn digamma_zeta_series_c(z: Complex<f64>, root: f64, rootval: f64) -> Complex<f64> {
    let mut res = Complex::from(rootval);
    let mut coeff = Complex::new(-1.0, 0.0);

    let z = z - root;
    for n in 1..100 {
        coeff *= -z;
        let term = coeff * crate::xsf::cephes::zeta(f64::from(n + 1), root);
        res += term;
        if term.norm() < f64::EPSILON * res.norm() {
            break;
        }
    }
    res
}

fn digamma_forward_recurrence(z: Complex<f64>, psiz: Complex<f64>, n: i32) -> Complex<f64> {
    /* Compute digamma(z + n) using digamma(z) using the recurrence
     * relation
     *
     * digamma(z + 1) = digamma(z) + 1/z.
     *
     * See https://dlmf.nist.gov/5.5#E2 */
    let mut res = psiz;

    for k in 0..n {
        res += cdiv(Complex::ONE, z + f64::from(k));
    }
    res
}

fn digamma_backward_recurrence(z: Complex<f64>, psiz: Complex<f64>, n: i32) -> Complex<f64> {
    /* Compute digamma(z - n) using digamma(z) and a recurrence relation. */
    let mut res = psiz;

    for k in 1..=n {
        res -= cdiv(Complex::ONE, z - f64::from(k));
    }
    res
}

fn digamma_asymptotic_series(z: Complex<f64>) -> Complex<f64> {
    /* Evaluate digamma using an asymptotic series. See
     *
     * https://dlmf.nist.gov/5.11#E2 */
    const BERNOULLI2K: [f64; 16] = [
        1.666_666_666_666_666_6e-1,
        -3.333_333_333_333_333e-2,
        2.380_952_380_952_380_8e-2,
        -3.333_333_333_333_333e-2,
        7.575_757_575_757_576e-2,
        -2.531_135_531_135_531e-1,
        1.166_666_666_666_666_7,
        -7.092_156_862_745_098,
        5.497_117_794_486_215_6e1,
        -5.291_242_424_242_424e2,
        6.192_123_188_405_797e3,
        -8.658_025_311_355_312e4,
        1.425_517_166_666_666_7e6,
        -2.729_823_106_781_609_4e7,
        6.015_808_739_006_424e8,
        -1.511_631_576_709_215_7e10,
    ];
    let rzz = cdiv(cdiv(Complex::ONE, z), z);
    let mut zfac = Complex::ONE;

    if !(z.re.is_finite() && z.im.is_finite()) {
        /* Check for infinity (or nan) and return early.
         * Result of division by complex infinity is implementation dependent.
         * and has been observed to vary between C++ stdlib and CUDA stdlib.
         */
        return clog(z);
    }

    let mut res = clog(z) - cdiv(0.5.into(), z);

    for (k, &b2k) in (1_u32..17).zip(&BERNOULLI2K) {
        zfac *= rzz;
        let term = -b2k * zfac / (2.0 * f64::from(k));
        res += term;
        if term.norm() < f64::EPSILON * res.norm() {
            break;
        }
    }
    res
}

/// `xsf::digamma` for real `z`
fn digamma_f64(z: f64) -> f64 {
    /* Wrap Cephes' psi to take advantage of the series expansion around
     * the smallest negative zero.
     */
    if (z - DIGAMMA_NEGROOT).abs() < 0.3 {
        return digamma_zeta_series(z, DIGAMMA_NEGROOT, DIGAMMA_NEGROOTVAL);
    }
    crate::xsf::cephes::psi(z)
}

/// `xsf::digamma` for complex `z`
#[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
fn digamma_c64(mut z: Complex<f64>) -> Complex<f64> {
    /*
     * Compute the digamma function for complex arguments. The strategy
     * is:
     *
     * - Around the two zeros closest to the origin (posroot and negroot)
     * use a Taylor series with precomputed zero order coefficient.
     * - If close to the origin, use a recurrence relation to step away
     * from the origin.
     * - If close to the negative real axis, use the reflection formula
     * to move to the right halfplane.
     * - If |z| is large (> 16), use the asymptotic series.
     * - If |z| is small, use a recurrence relation to make |z| large
     * enough to use the asymptotic series.
     */
    let mut absz = z.norm();
    let mut res = Complex::ZERO;
    /* Use the asymptotic series for z away from the negative real axis
     * with abs(z) > smallabsz. */
    let smallabsz = 16.0;
    /* Use the reflection principle for z with z.real < 0 that are within
     * smallimag of the negative real axis.
     * int smallimag = 6  # unused below except in a comment */

    if z.re <= 0.0 && z == Complex::from(z.re.ceil()) {
        // Poles
        return Complex::new(f64::NAN, f64::NAN);
    }
    if (z - DIGAMMA_NEGROOT).norm() < 0.3 {
        // First negative root.
        return digamma_zeta_series_c(z, DIGAMMA_NEGROOT, DIGAMMA_NEGROOTVAL);
    }

    if z.re < 0.0 && z.im.abs() < smallabsz {
        /* Reflection formula for digamma. See
         *
         *https://dlmf.nist.gov/5.5#E4
         */
        res = cdiv(-PI * crate::cospi(z), crate::sinpi(z));
        z = rsub(1.0, z);
        absz = z.norm();
    }

    if absz < 0.5 {
        /* Use one step of the recurrence relation to step away from
         * the pole. */
        res = cdiv((-1.0).into(), z);
        z += 1.0;
        absz = z.norm();
    }

    if (z - DIGAMMA_POSROOT).norm() < 0.5 {
        res += digamma_zeta_series_c(z, DIGAMMA_POSROOT, DIGAMMA_POSROOTVAL);
    } else if absz > smallabsz {
        res += digamma_asymptotic_series(z);
    } else if z.re >= 0.0 {
        let n = (smallabsz - absz).trunc() + 1.0;
        let init = digamma_asymptotic_series(z + n);
        res += digamma_backward_recurrence(z + n, init, n as i32);
    } else {
        // z.real() < 0, absz < smallabsz, and z.imag() > smallimag
        // (unlike in xsf, scipy/xsf#347) clamped to n >= 0
        let n = ((smallabsz - absz).trunc() - 1.0).max(0.0);
        let init = digamma_asymptotic_series(z - n);
        res += digamma_forward_recurrence(z - n, init, n as i32);
    }
    res
}

pub trait DigammaArg: crate::sealed::Sealed {
    fn digamma(self) -> Self;
}

impl DigammaArg for f64 {
    #[inline]
    fn digamma(self) -> Self {
        digamma_f64(self)
    }
}

impl DigammaArg for num_complex::Complex<f64> {
    #[inline]
    fn digamma(self) -> Self {
        digamma_c64(self)
    }
}

/// Digamma function $\psi(z)$ for real or complex input
///
/// This function is defined as [^DLMF]:
///
/// $$
/// \psi(z) = {\Gamma\'(z) \over \Gamma(z)} ,
/// $$
///
/// where $\Gamma$ is the [Gamma function](crate::gamma).
///
/// Corresponds to [`scipy.special.digamma`][digamma] in scipy
///
/// [digamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.digamma.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/5.2.E2>.
///
/// # See also
/// - [`polygamma`](crate::polygamma): Polygamma function (n-th derivative)
#[doc(alias = "psi")]
#[must_use]
#[inline]
pub fn digamma<T: DigammaArg>(x: T) -> T {
    x.digamma()
}
