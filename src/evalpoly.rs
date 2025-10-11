use crate::ffi;

/// Evaluate polynomials
///
/// All of the coefficients are stored in reverse order, i.e. if the polynomial is:
///
/// uₙxⁿ + uₙ₋₁xⁿ⁻¹ + … + u₀
///
/// then `coeffs[0]` = uₙ, `coeffs[1]` = uₙ₋₁, …, `coeffs[n]` = u₀
///
/// # Arguments
///
/// - `coeffs`: Polynomial coefficients in reverse order
/// - `z`: Complex value at which to evaluate the polynomial
///
/// # Returns
/// - `p(z)`: Value of the polynomial evaluated at `z`
pub fn cevalpoly(coeffs: &[f64], z: num_complex::Complex<f64>) -> num_complex::Complex<f64> {
    let degree = (coeffs.len() as i32) - 1;
    if degree == -1 {
        unsafe { ffi::cevalpoly([0.0, 0.0].as_ptr(), 1, z.into()) }
    } else if degree == 0 {
        unsafe { ffi::cevalpoly([0.0, coeffs[0]].as_ptr(), 1, z.into()) }
    } else {
        unsafe { ffi::cevalpoly(coeffs.as_ptr(), degree, z.into()) }
    }
    .into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::c64;

    #[test]
    fn test_cevalpoly_0() {
        // p(z) = 0
        let y = cevalpoly(&[], c64(2.0, 3.0));
        assert_eq!(y.re, 0.0);
        assert_eq!(y.im, 0.0);
    }

    #[test]
    fn test_cevalpoly_1() {
        // p(z) = 5
        let y = cevalpoly(&[5.0], c64(2.0, 3.0));
        // p(2+3i) = 5
        assert_eq!(y.re, 5.0);
        assert_eq!(y.im, 0.0);
    }

    #[test]
    fn test_cevalpoly_2() {
        // p(z) = 2z + 3
        let y = cevalpoly(&[2.0, 3.0], c64(1.0, 1.0));
        // p(1+i) = 5 + 2i
        assert_eq!(y.re, 5.0);
        assert_eq!(y.im, 2.0);
    }
}
