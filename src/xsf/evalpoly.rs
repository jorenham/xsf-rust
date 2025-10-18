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
#[doc(alias = "evalpoly")]
#[doc(alias = "polynomial")]
pub fn cevalpoly(coeffs: &[f64], z: num_complex::Complex<f64>) -> num_complex::Complex<f64> {
    let degree = (coeffs.len() as i32) - 1;
    if degree == -1 {
        unsafe { crate::ffi::xsf::cevalpoly([0.0, 0.0].as_ptr(), 1, z) }
    } else if degree == 0 {
        unsafe { crate::ffi::xsf::cevalpoly([0.0, coeffs[0]].as_ptr(), 1, z) }
    } else {
        unsafe { crate::ffi::xsf::cevalpoly(coeffs.as_ptr(), degree, z) }
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
}
