use num_complex::Complex;

/// Lambert W function.
///
/// The Lambert W function `W(z)` is defined as the inverse function of `w * exp(w)`. In other
/// words, the value of `W(z)` is such that `z = W(z) * exp(W(z))` for any complex number `z`.
///
/// The Lambert W function is a multivalued function with infinitely many branches. Each branch
/// gives a separate solution of the equation `z = w exp(w)`. Here, the branches are indexed by the
/// integer `k`.
#[doc(alias = "lambert_w")]
pub fn lambertw(z: Complex<f64>, k: isize, tol: f64) -> Complex<f64> {
    unsafe { crate::ffi::xsf::lambertw(z, k as core::ffi::c_long, tol) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_lambertw_c64() {
        xsref::test("lambertw", "cd_p_d-cd", |x| {
            crate::lambertw(c64(x[0], x[1]), x[2] as isize, x[3])
        });
    }
}
