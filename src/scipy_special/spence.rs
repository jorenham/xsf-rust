//! NOTE: The Cephes library only provides a real implementation of the Spence function.
//! For complex argument, `scipy.special.spence` uses the Cython implementation at
//! `scipy/special/_spence.pxd`.

use num_complex::{Complex64, c64};
use num_traits::Zero;

// Relative tolerance for the series
const TOL: f64 = f64::EPSILON;

// π^2 / 6
#[allow(clippy::excessive_precision)]
const PISQ_6: f64 = 1.644_934_066_848_226_436_5;

/// Compute Spence's function for complex arguments. The strategy is:
/// - If z is close to 0, use a series centered at 0.
/// - If z is far away from 1, use the reflection formula
///
/// spence(z) = -spence(z/(z - 1)) - pi**2/6 - ln(z - 1)**2/2
///
/// to move close to 1. See [1].
/// - If z is close to 1, use a series centered at 1.
#[inline]
fn cspence(z: Complex64) -> Complex64 {
    // if zabs(z) < 0.5:
    //     # This step isn't necessary, but this series converges faster.
    //     return cspence_series0(z)
    // elif zabs(1 - z) > 1:
    //     return -cspence_series1(z/(z - 1)) - PISQ_6 - 0.5*zlog1(z - 1)**2
    // else:
    //     return cspence_series1(z)

    if z.norm() < 0.5 {
        // This step isn't necessary, but this series converges faster.
        cspence_series0(z)
    } else if (1.0 - z).norm() > 1.0 {
        let zm1 = z - 1.0;
        -cspence_series1(z / zm1) - PISQ_6 - 0.5 * zm1.ln().powi(2)
    } else {
        cspence_series1(z)
    }
}

/// A series centered at z = 0; see <https://functions.wolfram.com/10.07.06.0005.02>
#[inline]
fn cspence_series0(z: Complex64) -> Complex64 {
    if z.is_zero() {
        return c64(PISQ_6, 0.0);
    }
    // cdef:
    //     int n
    //     double complex zfac = 1
    //     double complex sum1 = 0
    //     double complex sum2 = 0
    //     double complex term1, term2
    let mut zfac = c64(1.0, 0.0);
    let mut sum1 = c64(0.0, 0.0);
    let mut sum2 = c64(0.0, 0.0);

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

    // return PISQ_6 - sum1 + zlog1(z)*sum2
    PISQ_6 - sum1 + z.ln() * sum2
}

/// A series centered at z = 1 which enjoys faster convergence than the Taylor series.
/// See Ginsberg, Zaborowski, "The Dilogarithm Function of a Real Argument"].
/// The number of terms used comes from bounding the absolute tolerance at the edge of the radius of
/// convergence where the sum is O(1).
#[inline]
fn cspence_series1(z: Complex64) -> Complex64 {
    // z = 1 - z
    let z = 1.0 - z;

    // if z == 1:
    //     return 0
    if z.is_zero() {
        return z;
    }

    // cdef:
    //     int n
    //     double complex zfac = 1
    //     double complex res = 0
    //     double complex term, zz
    let mut zfac = c64(1.0, 0.0);
    let mut res = c64(0.0, 0.0);

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
    res += 4.0 * z + 5.75 * zz + 3.0 * (1.0 - zz) * (1.0 - z).ln();
    // res /= 1 + 4*z + zz
    res /= 1.0 + 4.0 * z + zz;
    // return res
    res
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait SpenceArg: sealed::Sealed {
    fn spence(self) -> Self;
}

impl SpenceArg for f64 {
    #[inline]
    fn spence(self) -> Self {
        unsafe { crate::ffi::xsf::spence(self) }
    }
}

impl SpenceArg for Complex64 {
    #[inline]
    fn spence(self) -> Self {
        cspence(self)
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

#[cfg(test)]
mod tests {
    use num_complex::c64;
    use xsref::np_assert_allclose;

    #[test]
    fn test_spence_xsref_f64() {
        xsref::test("spence", "d-d", |x| crate::spence(x[0]));
    }

    #[test]
    fn test_spence_xsref_c64() {
        xsref::test("spence", "d-d", |x| crate::spence(c64(x[0], 0.0)).re);
    }

    #[test]
    fn test_spence_manual_c64() {
        // values calculated using Wolfram|Alpha with PolyLog[2, 1-z]
        let zs = [
            c64(-1.0, -1.0),
            c64(-1.0, 1.0),
            c64(0.0, -1.0),
            c64(0.0, 1.0),
            c64(1.0, -1.0),
            c64(1.0, 1.0),
        ];
        let expect = [
            c64(1.186_688_537_000_058, 2.407_740_769_345_772),
            c64(1.186_688_537_000_058, -2.407_740_769_345_772),
            c64(0.616_850_275_068_085, 1.460_362_116_753_12),
            c64(0.616_850_275_068_085, -1.460_362_116_753_12),
            c64(-0.205_616_758_356_028_3, 0.915_965_594_177_219),
            c64(-0.205_616_758_356_028_3, -0.915_965_594_177_219),
        ];
        np_assert_allclose!(zs.map(crate::spence), expect, rtol = 1e-13);
    }
}
