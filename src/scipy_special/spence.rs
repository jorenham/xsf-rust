//! NOTE: The Cephes library only provides a real implementation of the Spence function.
//! For complex argument, `scipy.special.spence` uses the Cython implementation at
//! `scipy/special/_spence.pxd`.
//!
//! Cython converts the real operand of a complex subtraction to complex, so `1 - z` is
//! `C64::ONE - z` here, with the imaginary part `0.0 - z.im`.

use num_complex::Complex64;

use crate::utils::FRAC_PI_SQ_6;
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::zlog1_c64;

// Relative tolerance for the series
const TOL: f64 = f64::EPSILON;

/// Compute Spence's function for complex arguments. The strategy is:
/// - If z is close to 0, use a series centered at 0.
/// - If z is far away from 1, use the reflection formula
///
/// spence(z) = -spence(z/(z - 1)) - pi**2/6 - ln(z - 1)**2/2
///
/// to move close to 1. See \[1\].
/// - If z is close to 1, use a series centered at 1.
#[inline]
fn cspence(z: C64) -> C64 {
    // if zabs(z) < 0.5:
    //     # This step isn't necessary, but this series converges faster.
    //     return cspence_series0(z)
    // elif zabs(1 - z) > 1:
    //     return -cspence_series1(z/(z - 1)) - FRAC_PI_SQ_6 - 0.5*zlog1(z - 1)**2
    // else:
    //     return cspence_series1(z)

    if z.is_nan() {
        // like in xsf, instead of iterating 500 times
        C64::NAN
    } else if z.norm() < 0.5 {
        // This step isn't necessary, but this series converges faster.
        cspence_series0(z)
    } else if (C64::ONE - z).norm() > 1.0 {
        let zm1 = z - 1.0;
        let log_zm1 = zlog1_c64(zm1);
        -cspence_series1(z / zm1) - FRAC_PI_SQ_6 - 0.5 * (log_zm1 * log_zm1)
    } else {
        cspence_series1(z)
    }
}

/// A series centered at z = 0; see <https://functions.wolfram.com/10.07.06.0005.02>
#[inline]
fn cspence_series0(z: C64) -> C64 {
    if z == C64::ZERO {
        return c64(FRAC_PI_SQ_6, 0.0);
    }
    // cdef:
    //     int n
    //     double complex zfac = 1
    //     double complex sum1 = 0
    //     double complex sum2 = 0
    //     double complex term1, term2
    let mut zfac = C64::ONE;
    let mut sum1 = C64::ZERO;
    let mut sum2 = C64::ZERO;

    // for n in range(1, 500):
    for n in 1..500 {
        // zfac *= z
        zfac *= z;

        // term1 = zfac/n**2
        let term1 = zfac / f64::from(n * n);
        // sum1 += term1
        sum1 += term1;

        // term2 = zfac/n
        let term2 = zfac / f64::from(n);
        // sum2 += term2
        sum2 += term2;

        // if zabs(term1) <= TOL*zabs(sum1) and zabs(term2) <= TOL*zabs(sum2):
        //     break
        if term1.norm() <= TOL * sum1.norm() && term2.norm() <= TOL * sum2.norm() {
            break;
        }
    }

    // return FRAC_PI_SQ_6 - sum1 + zlog1(z)*sum2
    c64(FRAC_PI_SQ_6, 0.0) - sum1 + zlog1_c64(z) * sum2
}

/// A series centered at z = 1 which enjoys faster convergence than the Taylor series.
/// See Ginsberg, Zaborowski, "The Dilogarithm Function of a Real Argument"].
/// The number of terms used comes from bounding the absolute tolerance at the edge of the radius of
/// convergence where the sum is O(1).
#[inline]
fn cspence_series1(z: C64) -> C64 {
    // z = 1 - z
    let z = C64::ONE - z;

    // if z == 1:
    //     return 0
    if z == C64::ZERO {
        return z;
    }

    // cdef:
    //     int n
    //     double complex zfac = 1
    //     double complex res = 0
    //     double complex term, zz
    let mut zfac = C64::ONE;
    let mut res = C64::ZERO;

    // for n in range(1, 500):
    //     if zabs(term) <= TOL*zabs(res):
    //         break
    for n in 1..500 {
        // zfac *= z
        zfac *= z;
        // # Do the divisions one at a time to guard against overflow
        // term = ((zfac/n**2)/(n + 1)**2)/(n + 2)**2
        let n = f64::from(n);
        let term = zfac / n.powi(2) / (n + 1.0).powi(2) / (n + 2.0).powi(2);
        // res += term
        res += term;
        if term.norm() <= TOL * res.norm() {
            break;
        }
    }

    // zz = z**2
    let zz = z * z;

    // res *= 4*zz
    res *= 4.0 * zz;
    // res += 4*z + 5.75*zz + 3*(1 - zz)*zlog1(1 - z)
    res += 4.0 * z + 5.75 * zz + 3.0 * (C64::ONE - zz) * zlog1_c64(C64::ONE - z);
    // res /= 1 + 4*z + zz
    res /= 1.0 + 4.0 * z + zz;
    // return res
    res
}

pub trait SpenceArg: crate::sealed::Sealed {
    fn spence(self) -> Self;
}

impl SpenceArg for f64 {
    #[inline]
    fn spence(self) -> Self {
        crate::xsf::cephes::spence(self)
    }
}

impl SpenceArg for Complex64 {
    #[inline]
    fn spence(self) -> Self {
        self.via(cspence)
    }
}

/// Spence's function for real or complex argument, also known as the dilogarithm
///
/// $$ \int_1^z {\log(t) \over 1 - t} \dd t $$
///
/// Corresponds to [`scipy.special.spence`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spence.html
#[doc(alias = "dilogarithm", alias = "li2")]
#[must_use]
#[inline]
pub fn spence<T: SpenceArg>(z: T) -> T {
    z.spence()
}
