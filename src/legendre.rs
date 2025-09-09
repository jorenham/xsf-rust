use crate::bindings;
use core::ffi::c_int;

use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait LegendrePArg: sealed::Sealed + Sized {
    fn legendre_p(self, n: c_int) -> Self;
    fn legendre_p_all<const N: usize>(self) -> [Self; N];
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self;
    fn assoc_legendre_p(self, n: c_int, m: c_int, branch_cut: c_int) -> Self;
    fn assoc_legendre_p_norm(self, n: c_int, m: c_int, branch_cut: c_int) -> Self;
}

impl LegendrePArg for f64 {
    #[inline(always)]
    fn legendre_p(self, n: c_int) -> f64 {
        unsafe { bindings::legendre_p(n, self) }
    }
    #[inline(always)]
    fn legendre_p_all<const N: usize>(self) -> [Self; N] {
        let mut pn = [0.0; N];
        unsafe {
            bindings::legendre_p_all(N - 1, self, pn.as_mut_ptr());
        }
        pn
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

impl LegendrePArg for Complex<f64> {
    #[inline(always)]
    fn legendre_p(self, n: c_int) -> Complex<f64> {
        unsafe { bindings::legendre_p_1(n, self.into()) }.into()
    }
    #[inline(always)]
    fn legendre_p_all<const N: usize>(self) -> [Self; N] {
        let mut pn = bindings::complex_zeros::<N>();
        unsafe {
            bindings::legendre_p_all_1(N - 1, self.into(), pn.as_mut_ptr());
        }
        pn.map(|c| c.into())
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

pub trait LegendreQArg: sealed::Sealed + Sized {
    fn legendre_q_all<const N: usize>(self) -> ([Self; N], [Self; N]);
}

impl LegendreQArg for f64 {
    #[inline(always)]
    fn legendre_q_all<const N: usize>(self) -> ([Self; N], [Self; N]) {
        let mut qn = [0.0; N];
        let mut qd = [0.0; N];

        unsafe {
            bindings::lqn(N - 1, self, qn.as_mut_ptr(), qd.as_mut_ptr());
        }

        (qn, qd)
    }
}

impl LegendreQArg for Complex<f64> {
    #[inline(always)]
    fn legendre_q_all<const N: usize>(self) -> ([Self; N], [Self; N]) {
        let mut cqn = bindings::complex_zeros::<N>();
        let mut cqd = bindings::complex_zeros::<N>();

        unsafe {
            bindings::lqn_1(N - 1, self.into(), cqn.as_mut_ptr(), cqd.as_mut_ptr());
        }

        (cqn.map(|c| c.into()), cqd.map(|c| c.into()))
    }
}

/// Legendre polynomial of degree n
pub fn legendre_p<T: LegendrePArg>(n: i32, z: T) -> T {
    z.legendre_p(n as c_int)
}

/// Sequence of Legendre functions of the 1st kind
///
/// # Example
///
/// ```
/// use xsf::legendre_p_all;
///
/// // Compute P_n(z) for degrees 0, 1, 2, 3 (N=4)
/// let pn: [f64; 4] = legendre_p_all::<_, 4>(0.5);
/// ```
pub fn legendre_p_all<T: LegendrePArg, const N: usize>(z: T) -> [T; N] {
    z.legendre_p_all()
}

/// Spherical Legendre polynomial of degree n and order m
pub fn sph_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// Associated Legendre polynomial of the first kind
pub fn assoc_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2)
}

/// Normalized associated Legendre polynomial of the first kind
pub fn assoc_legendre_p_norm<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p_norm(n as c_int, m as c_int, 2)
}

/// Sequence of Legendre functions of the 2nd kind
///
/// Compute sequence of Legendre functions of the second kind, *Qn(z)* and derivatives for all
/// degrees from *0* to *n* (inclusive).
///
/// The array size `N` must be specified as a const generic parameter and determines the maximum
/// degree computed (*n* = `N` - 1).
///
/// # Example
///
/// ```
/// use xsf::legendre_q_all;
///
/// // Compute Q_n(z) and Q'_n(z) for degrees 0, 1, 2, 3, 4 (N=5)
/// let (qn, qn_deriv) = legendre_q_all::<_, 5>(0.5);
/// ```
///
/// For complex input
///
/// ```
/// use num_complex::c64;
/// use xsf::legendre_q_all;
///
/// let (qn, qn_deriv) = legendre_q_all::<_, 3>(c64(0.5, 0.3));
/// ```
#[doc(alias = "lqn")]
pub fn legendre_q_all<T: LegendreQArg, const N: usize>(z: T) -> ([T; N], [T; N]) {
    z.legendre_q_all()
}

#[cfg(test)]
mod tests {
    // there are no xsref tables, so we use smoke tests instead

    use super::*;
    use num_complex::c64;
    use std::f64::consts;

    const I: Complex<f64> = Complex::new(0.0, 1.0);
    const LN_3: f64 = 1.098_612_288_668_109_8_f64;

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

    // legendre_p_all

    #[test]
    fn test_legendre_p_all_f64() {
        assert_eq!(legendre_p_all::<_, 3>(0.0), [1.0, 0.0, -0.5]);
    }

    #[test]
    fn test_legendre_p_all_c64() {
        assert_eq!(
            legendre_p_all::<_, 3>(I),
            [c64(1.0, 0.0), I, c64(-2.0, 0.0)]
        );
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

    #[test]
    fn test_legendre_q_all_f64() {
        let (qn, qd) = legendre_q_all::<_, 5>(0.5);

        assert_eq!(qn.len(), 5);
        assert_eq!(qd.len(), 5);

        assert_eq!(qn[0], LN_3 * 0.5);
        assert_eq!(qn[1], LN_3 * 0.25 - 1.0);

        assert_eq!(qd[0], 4.0 / 3.0);
        assert_eq!(qd[1], 2.0 / 3.0 + 0.5 * LN_3);
    }

    #[test]
    fn test_legendre_q_all_c64() {
        let (qn, qd) = legendre_q_all::<_, 5>(c64(0.0, 1.0));

        assert_eq!(qn.len(), 5);
        assert_eq!(qd.len(), 5);

        assert_eq!(qn[0], c64(0.0, consts::FRAC_PI_4));
        assert_eq!(qn[1], c64(-1.0 - consts::FRAC_PI_4, 0.0));
    }
}
