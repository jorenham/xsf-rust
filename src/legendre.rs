use crate::bindings;
use alloc::vec::Vec;
use core::ffi::c_int;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait LegendrePArg: sealed::Sealed + Sized {
    fn legendre_p(self, n: c_int) -> Self;
    fn legendre_p_all(self, n: usize) -> Vec<Self>;
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self;
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>>;
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self;
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>>;
}

impl LegendrePArg for f64 {
    #[inline(always)]
    fn legendre_p(self, n: c_int) -> f64 {
        unsafe { bindings::legendre_p(n, self) }
    }

    #[inline(always)]
    fn legendre_p_all(self, n: usize) -> Vec<Self> {
        let mut pn = alloc::vec![0.0; n + 1];
        unsafe {
            bindings::legendre_p_all(n, self, pn.as_mut_ptr());
        }
        pn
    }

    #[inline(always)]
    fn sph_legendre_p(self, n: c_int, m: c_int) -> f64 {
        unsafe { bindings::sph_legendre_p(n, m, self) }
    }

    #[inline(always)]
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = alloc::vec![0.0; ni * nj];
        unsafe {
            bindings::sph_legendre_p_all(n, m, self, pnm.as_mut_ptr());
        }
        (0..ni)
            .map(|i| pnm[i * nj..(i + 1) * nj].to_vec())
            .collect()
    }

    #[inline(always)]
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> f64 {
        unsafe {
            if norm {
                bindings::assoc_legendre_p_1(n, m, self, bc)
            } else {
                bindings::assoc_legendre_p_0(n, m, self, bc)
            }
        }
    }

    #[inline(always)]
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = alloc::vec![0.0; ni * nj];
        unsafe {
            if norm {
                bindings::assoc_legendre_p_all_1(n, m, self, bc, pnm.as_mut_ptr());
            } else {
                bindings::assoc_legendre_p_all_0(n, m, self, bc, pnm.as_mut_ptr());
            }
        }
        (0..ni)
            .map(|i| pnm[i * nj..(i + 1) * nj].to_vec())
            .collect()
    }
}

impl LegendrePArg for Complex<f64> {
    #[inline(always)]
    fn legendre_p(self, n: c_int) -> Complex<f64> {
        unsafe { bindings::legendre_p_1(n, self.into()) }.into()
    }

    #[inline(always)]
    fn legendre_p_all(self, n: usize) -> Vec<Self> {
        let mut pn = bindings::complex_zeros(n + 1);
        unsafe {
            bindings::legendre_p_all_1(n, self.into(), pn.as_mut_ptr());
        }
        bindings::cvec_into(pn)
    }

    #[inline(always)]
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Complex<f64> {
        unsafe { bindings::sph_legendre_p_1(n, m, self.into()) }.into()
    }

    #[inline(always)]
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = bindings::complex_zeros(ni * nj);
        unsafe {
            bindings::sph_legendre_p_all_1(n, m, self.into(), pnm.as_mut_ptr());
        }
        let pnm = bindings::cvec_into(pnm);
        (0..ni)
            .map(|i| pnm[i * nj..(i + 1) * nj].to_vec())
            .collect()
    }

    #[inline(always)]
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Complex<f64> {
        unsafe {
            if norm {
                bindings::assoc_legendre_p_1_1(n, m, self.into(), bc)
            } else {
                bindings::assoc_legendre_p_0_1(n, m, self.into(), bc)
            }
        }
        .into()
    }

    #[inline(always)]
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = bindings::complex_zeros(ni * nj);
        unsafe {
            if norm {
                bindings::assoc_legendre_p_all_1_1(n, m, self.into(), bc, pnm.as_mut_ptr());
            } else {
                bindings::assoc_legendre_p_all_0_1(n, m, self.into(), bc, pnm.as_mut_ptr());
            }
        }
        let pnm = bindings::cvec_into(pnm);
        (0..ni)
            .map(|i| pnm[i * nj..(i + 1) * nj].to_vec())
            .collect()
    }
}

pub trait LegendreQArg: sealed::Sealed + Sized {
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>);
}

impl LegendreQArg for f64 {
    #[inline(always)]
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
        let mut qn = alloc::vec![0.0; n + 1];
        let mut qd = alloc::vec![0.0; n + 1];

        unsafe {
            bindings::lqn(n, self, qn.as_mut_ptr(), qd.as_mut_ptr());
        }

        (qn, qd)
    }
}

impl LegendreQArg for Complex<f64> {
    #[inline(always)]
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
        let mut cqn: Vec<_> = bindings::complex_zeros(n + 1);
        let mut cqd: Vec<_> = bindings::complex_zeros(n + 1);

        unsafe {
            bindings::lqn_1(n, self.into(), cqn.as_mut_ptr(), cqd.as_mut_ptr());
        }

        (bindings::cvec_into(cqn), bindings::cvec_into(cqd))
    }
}

/// Legendre polynomial of degree n
pub fn legendre_p<T: LegendrePArg>(n: i32, z: T) -> T {
    z.legendre_p(n as c_int)
}

/// All Legendre polynomials of the 1st kind up to degree `n`
///
/// Output length is `n + 1`. The entry at `j` corresponds to degree `j` for all `0 <= j <= n`.
pub fn legendre_p_all<T: LegendrePArg>(n: usize, z: T) -> Vec<T> {
    z.legendre_p_all(n)
}

/// Spherical Legendre polynomial of degree `n` and order `m`
pub fn sph_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// All spherical Legendre polynomials of the 1st kind up to degree `n` and order `m`
///
/// Output shape is `(n + 1, 2 * m + 1)`. The entry at `(j, i)` corresponds to degree `j` and
/// order `i` for all  `0 <= j <= n` and `-m <= i <= m`.
pub fn sph_legendre_p_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.sph_legendre_p_all(n, m)
}

/// Associated Legendre polynomial of the 1st kind
pub fn assoc_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2, false)
}

/// All associated Legendre polynomials of the 1st kind up to degree `n` and order `m`
pub fn assoc_legendre_p_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.assoc_legendre_p_all(n, m, 2, false)
}

/// Normalized associated Legendre polynomial of the 1st kind
pub fn assoc_legendre_p_norm<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2, true)
}

/// All normalized associated Legendre polynomials of the 1st kind up to degree `n` and order `m`
pub fn assoc_legendre_p_norm_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.assoc_legendre_p_all(n, m, 2, true)
}

/// All Legendre polynomials of the 2nd kind up to degree `n`
///
/// Output length is `n + 1`. The entry at `j` corresponds to degree `j` for all `0 <= j <= n`.
#[doc(alias = "lqn")]
pub fn legendre_q_all<T: LegendreQArg>(n: usize, z: T) -> (Vec<T>, Vec<T>) {
    z.legendre_q_all(n)
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
        assert_eq!(legendre_p_all(2, 0.0), vec![1.0, 0.0, -0.5]);
    }

    #[test]
    fn test_legendre_p_all_c64() {
        assert_eq!(legendre_p_all(2, I), vec![c64(1.0, 0.0), I, c64(-2.0, 0.0)]);
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

    // sph_legendre_p_all

    #[test]
    fn test_sph_legendre_p_all_f64() {
        let pnm = sph_legendre_p_all(1, 1, 0.0);

        assert_eq!(pnm.len(), 2);
        assert_eq!(pnm[0].len(), 3);
        assert_eq!(pnm[1].len(), 3);

        assert_eq!(pnm[0], vec![0.28209479177387814, 0.0, 0.0]);
        assert_eq!(pnm[1], vec![0.48860251190291987, 0.0, 0.0]);
    }

    #[test]
    fn test_sph_legendre_p_all_c64() {
        let pnm = sph_legendre_p_all(1, 1, I);
        assert_eq!(pnm.len(), 2);
        assert_eq!(pnm[0].len(), 3);
        assert_eq!(pnm[1].len(), 3);

        assert!((pnm[0][0].norm() - 0.2820947917738782).abs() < f64::EPSILON);
        assert!((pnm[1][0].norm() - 0.7539530742394804).abs() < f64::EPSILON);
        assert!((pnm[1][1].norm() - 0.4060251368556634).abs() < f64::EPSILON);
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

    // assoc_legendre_p_all

    #[test]
    fn test_assoc_legendre_p_all_f64() {
        let pnm = assoc_legendre_p_all(1, 1, 0.0);

        assert_eq!(pnm.len(), 2);
        assert_eq!(pnm[0].len(), 3);
        assert_eq!(pnm[1].len(), 3);

        assert_eq!(pnm[0], vec![1.0, 0.0, 0.0]);
        assert_eq!(pnm[1], vec![0.0, -1.0, 0.5]);
    }

    #[test]
    fn test_assoc_legendre_p_all_c64() {
        let pnm = assoc_legendre_p_all(1, 1, I);

        assert_eq!(pnm.len(), 2);
        assert_eq!(pnm[0].len(), 3);
        assert_eq!(pnm[1].len(), 3);

        assert_eq!(pnm[0], vec![c64(1.0, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)]);
        assert_eq!(
            pnm[1],
            vec![
                c64(0.0, 1.0),
                c64(-consts::SQRT_2, 0.0),
                c64(consts::FRAC_1_SQRT_2, 0.0)
            ]
        );
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
        let (qn, qd) = legendre_q_all(4, 0.5);

        assert_eq!(qn.len(), 5);
        assert_eq!(qd.len(), 5);

        assert_eq!(qn[0], LN_3 * 0.5);
        assert_eq!(qn[1], LN_3 * 0.25 - 1.0);

        assert_eq!(qd[0], 4.0 / 3.0);
        assert_eq!(qd[1], 2.0 / 3.0 + 0.5 * LN_3);
    }

    #[test]
    fn test_legendre_q_all_c64() {
        let (qn, qd) = legendre_q_all(4, c64(0.0, 1.0));

        assert_eq!(qn.len(), 5);
        assert_eq!(qd.len(), 5);

        assert_eq!(qn[0], c64(0.0, consts::FRAC_PI_4));
        assert_eq!(qn[1], c64(-1.0 - consts::FRAC_PI_4, 0.0));
    }
}
