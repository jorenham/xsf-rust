//! The real `digamma` is translated into pure Rust from `xsf/digamma.h` (xsf v0.2.2).

// All of the following were computed with mpmath
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
        unsafe { crate::ffi::xsf::digamma(self) }
    }
}

/// Digamma function for real or complex input
///
/// Corresponds to [`scipy.special.digamma`][digamma] in scipy
///
/// See also [`polygamma`](crate::polygamma).
///
/// [digamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.digamma.html
///
/// # See also
/// - [`polygamma`](crate::polygamma): Polygamma function (n-th derivative)
#[doc(alias = "psi")]
#[must_use]
#[inline]
pub fn digamma<T: DigammaArg>(x: T) -> T {
    x.digamma()
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_digamma_f64() {
        xsref::test("digamma", "d-d", |x| crate::digamma(x[0]));
    }

    #[test]
    fn test_digamma_near_negative_pole() {
        // just above a negative integer (reference values from mpmath)
        for (x, expected) in [
            (-0.999_999, -999_999.577_184_264_3),
            (-2.999_999_9, -9_999_998.760_247_82),
        ] {
            let rel = ((crate::digamma(x) - expected) / expected).abs();
            assert!(rel < 1e-15, "x = {x}: rel = {rel:e}");
        }
    }

    #[test]
    fn test_digamma_c64() {
        use num_complex::c64;

        xsref::test("digamma", "cd-cd", |x| crate::digamma(c64(x[0], x[1])));
    }
}
