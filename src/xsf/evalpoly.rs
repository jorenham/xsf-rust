//! Translated into pure Rust from `xsf/evalpoly.h` (xsf v0.2.2): translated from Cython into C++
//! by SciPy developers in 2024, original author: Josh Wilson, 2016.
//!
//! References
//! ----------
//! [1] Knuth, "The Art of Computer Programming, Volume II"

use num_complex::Complex;

/// `xsf::cevalpoly`, with `degree = coeffs.len() - 1`
fn xsf_cevalpoly(coeffs: &[f64], z: Complex<f64>) -> Complex<f64> {
    /* Evaluate a polynomial with real coefficients at a complex point.
     *
     * Uses equation (3) in section 4.6.4 of [1]. Note that it is more
     * efficient than Horner's method.
     */
    let mut a = coeffs[0];
    let mut b = coeffs[1];
    let r = 2.0 * z.re;
    let s = z.norm_sqr();

    for &c in &coeffs[2..] {
        let tmp = b;
        b = (-s).mul_add(a, c);
        a = r.mul_add(a, tmp);
    }

    z * a + b
}

/// Evaluate polynomials
///
/// All of the coefficients are stored in reverse order, i.e. if the polynomial is:
///
/// $$
/// u_n x^n + u_{n-1} x^{n-1} + \ldots + u_0
/// $$
///
/// then `coeffs[0]` = $u_n$, `coeffs[1]` = $u_{n-1}$, …, `coeffs[n]` = $u_0$
///
/// # Arguments
///
/// - `coeffs`: Polynomial coefficients in reverse order
/// - `z`: Complex value at which to evaluate the polynomial
///
/// # Returns
/// - `p(z)`: Value of the polynomial evaluated at `z`
#[doc(alias = "evalpoly", alias = "polynomial")]
#[must_use]
#[inline]
pub fn cevalpoly(coeffs: &[f64], z: Complex<f64>) -> Complex<f64> {
    match coeffs {
        [] => xsf_cevalpoly(&[0.0, 0.0], z),
        &[c] => xsf_cevalpoly(&[0.0, c], z),
        _ => xsf_cevalpoly(coeffs, z),
    }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_cevalpoly_0() {
        // p(z) = 0
        let y = crate::cevalpoly(&[], c64(2.0, 3.0));
        assert_eq!(y, c64(0.0, 0.0));
    }

    #[test]
    fn test_cevalpoly_1() {
        // p(z) = 5
        let y = crate::cevalpoly(&[5.0], c64(2.0, 3.0));
        // p(2+3i) = 5
        assert_eq!(y, c64(5.0, 0.0));
    }

    #[test]
    fn test_cevalpoly_2() {
        // p(z) = 2z + 3
        let y = crate::cevalpoly(&[2.0, 3.0], c64(1.0, 1.0));
        // p(1+i) = 5 + 2i
        assert_eq!(y, c64(5.0, 2.0));
    }

    #[test]
    fn test_cevalpoly_4() {
        // p(z) = z^3 + 2z^2 + 3z + 4
        let y = crate::cevalpoly(&[1.0, 2.0, 3.0, 4.0], c64(1.0, 1.0));
        // p(1+i) = 5 + 9i
        assert_eq!(y, c64(5.0, 9.0));
    }
}
