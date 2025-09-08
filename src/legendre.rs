use crate::bindings;
use std::os::raw::c_int;

use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait LegendreArg: sealed::Sealed {
    fn legendre_p(self, n: c_int) -> Self;
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self;
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> Self;
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> Self;
}

impl LegendreArg for f64 {
    #[inline(always)]
    fn legendre_p(self, n: c_int) -> f64 {
        unsafe { bindings::legendre_p(n, self) }
    }
    #[inline(always)]
    fn sph_legendre_p(self, n: c_int, m: c_int) -> f64 {
        unsafe { bindings::sph_legendre_p(n, m, self) }
    }
    #[inline(always)]
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> f64 {
        unsafe { bindings::assoc_legendre_p_0(n, m, self, branch_cut) }
    }
    #[inline(always)]
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> f64 {
        unsafe { bindings::assoc_legendre_p_1(n, m, self, branch_cut) }
    }
}

impl LegendreArg for Complex<f64> {
    #[inline(always)]
    fn legendre_p(self, n: c_int) -> Complex<f64> {
        unsafe { bindings::legendre_p_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Complex<f64> {
        unsafe { bindings::sph_legendre_p_1(n, m, self.into()) }.into()
    }
    #[inline(always)]
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> Complex<f64> {
        unsafe { bindings::assoc_legendre_p_0_1(n, m, self.into(), branch_cut) }.into()
    }
    #[inline(always)]
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> Complex<f64> {
        unsafe { bindings::assoc_legendre_p_1_1(n, m, self.into(), branch_cut) }.into()
    }
}

/// Legendre polynomial of degree n
pub fn legendre_p<T: LegendreArg>(n: i32, z: T) -> T {
    z.legendre_p(n as c_int)
}

/// Spherical Legendre polynomial of degree n and order m
pub fn sph_legendre_p<T: LegendreArg>(n: i32, m: i32, z: T) -> T {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// Associated Legendre polynomial of the first kind
pub fn assoc_legendre_p<T: LegendreArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2)
}

/// Normalized associated Legendre polynomial of the first kind
pub fn assoc_legendre_p_norm<T: LegendreArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p_norm(n as c_int, m as c_int, 2)
}

#[cfg(test)]
mod tests {
    // there are no xsref tables, so we use smoke tests instead

    use super::*;
    use num_complex::c64;
    use std::f64::consts;

    const I: Complex<f64> = Complex::new(0.0, 1.0);

    // legendre_p

    #[test]
    fn test_legendre_p_f64() {
        assert_eq!(legendre_p(0, 0.0), 1.0);
        assert_eq!(legendre_p(1, 0.0), 0.0);
        assert_eq!(legendre_p(2, 0.0), -0.5);
    }

    #[test]
    fn test_legendre_p_c64() {
        assert_eq!(legendre_p(0, I), c64(1.0, 0.0));
        assert_eq!(legendre_p(1, I), I);
        assert_eq!(legendre_p(2, I), c64(-2.0, 0.0));
    }

    // sph_legendre_p

    #[test]
    fn test_sph_legendre_p_f64() {
        assert_eq!(sph_legendre_p(0, 0, 0.0), 0.28209479177387814);
        assert_eq!(sph_legendre_p(1, 0, 0.0), 0.48860251190291987);
        assert_eq!(sph_legendre_p(1, 1, 0.0), 0.0);
    }

    #[test]
    fn test_sph_legendre_p_c64() {
        assert!((sph_legendre_p(0, 0, I).norm() - 0.2820947917738782).abs() < f64::EPSILON);
        assert!((sph_legendre_p(1, 0, I).norm() - 0.7539530742394804).abs() < f64::EPSILON);
        assert!((sph_legendre_p(1, 1, I).norm() - 0.4060251368556634).abs() < f64::EPSILON);
    }

    // assoc_legendre_p

    #[test]
    fn test_assoc_legendre_p_f64() {
        assert_eq!(assoc_legendre_p(0, 0, 0.0), 1.0);
        assert_eq!(assoc_legendre_p(1, 0, 0.0), 0.0);
        assert_eq!(assoc_legendre_p(0, 1, 0.0), 0.0);
        assert_eq!(assoc_legendre_p(1, 1, 0.0), -1.0);
    }

    #[test]
    fn test_assoc_legendre_p_c64() {
        assert_eq!(assoc_legendre_p(0, 0, I), c64(1.0, 0.0));
        assert_eq!(assoc_legendre_p(1, 0, I), I);
        assert_eq!(assoc_legendre_p(0, 1, I), c64(0.0, 0.0));
        assert_eq!(assoc_legendre_p(1, 1, I), c64(-consts::SQRT_2, 0.0));
    }

    // assoc_legendre_p_norm

    #[test]
    fn test_assoc_legendre_p_norm_f64() {
        assert!((assoc_legendre_p_norm(0, 0, 0.0) - consts::FRAC_1_SQRT_2).abs() < f64::EPSILON);
        assert_eq!(assoc_legendre_p_norm(1, 0, 0.0), 0.0);
        assert_eq!(assoc_legendre_p_norm(0, 1, 0.0), 0.0);
        assert_eq!(
            assoc_legendre_p_norm(1, 1, 0.0),
            -0.8660254037844386 // -sqrt(3 / 4)
        );
    }

    #[test]
    fn test_assoc_legendre_p_norm_c64() {
        assert!((assoc_legendre_p_norm(0, 0, I).re - consts::FRAC_1_SQRT_2).abs() < f64::EPSILON);
        assert_eq!(
            assoc_legendre_p_norm(1, 0, I),
            c64(0.0, 1.224744871391589) // sqrt(3 / 2) i
        );
        assert_eq!(assoc_legendre_p_norm(0, 1, I), c64(0.0, 0.0));
        assert_eq!(
            assoc_legendre_p_norm(1, 1, I),
            c64(-1.2247448713915892, 0.0) // -sqrt(3 / 2)
        );
    }
}
