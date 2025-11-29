use crate::utils;
use core::ffi::c_int;

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
    #[inline]
    fn legendre_p(self, n: c_int) -> Self {
        unsafe { crate::ffi::xsf::legendre_p(n, self) }
    }

    #[inline]
    fn legendre_p_all(self, n: usize) -> Vec<Self> {
        let mut pn = vec![f64::NAN; n + 1];
        unsafe {
            crate::ffi::xsf::legendre_p_all(n, self, pn.as_mut_ptr());
        }
        pn
    }

    #[inline]
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self {
        unsafe { crate::ffi::xsf::sph_legendre_p(n, m, self) }
    }

    #[inline]
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = vec![f64::NAN; ni * nj];
        unsafe {
            crate::ffi::xsf::sph_legendre_p_all(n, m, self, pnm.as_mut_ptr());
        }
        utils::vec_to_vecvec(&pnm, ni, nj, false)
    }

    #[inline]
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self {
        unsafe {
            if norm {
                crate::ffi::xsf::assoc_legendre_p_1(n, m, self, bc)
            } else {
                crate::ffi::xsf::assoc_legendre_p_0(n, m, self, bc)
            }
        }
    }

    #[inline]
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = vec![f64::NAN; ni * nj];
        unsafe {
            if norm {
                crate::ffi::xsf::assoc_legendre_p_all_1(n, m, self, bc, pnm.as_mut_ptr());
            } else {
                crate::ffi::xsf::assoc_legendre_p_all_0(n, m, self, bc, pnm.as_mut_ptr());
            }
        }
        utils::vec_to_vecvec(&pnm, ni, nj, false)
    }
}

impl LegendrePArg for num_complex::Complex<f64> {
    #[inline]
    fn legendre_p(self, n: c_int) -> Self {
        unsafe { crate::ffi::xsf::legendre_p_1(n, self) }
    }

    #[inline]
    fn legendre_p_all(self, n: usize) -> Vec<Self> {
        let mut pn = vec![num_complex::Complex::new(f64::NAN, 0.0); n + 1];
        unsafe {
            crate::ffi::xsf::legendre_p_all_1(n, self, pn.as_mut_ptr());
        }
        pn
    }

    #[inline]
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self {
        unsafe { crate::ffi::xsf::sph_legendre_p_1(n, m, self) }
    }

    #[inline]
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = vec![num_complex::Complex::new(f64::NAN, 0.0); ni * nj];
        unsafe {
            crate::ffi::xsf::sph_legendre_p_all_1(n, m, self, pnm.as_mut_ptr());
        }
        utils::vec_into_vecvec(pnm, ni, nj, false)
    }

    #[inline]
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self {
        unsafe {
            if norm {
                crate::ffi::xsf::assoc_legendre_p_1_1(n, m, self, bc)
            } else {
                crate::ffi::xsf::assoc_legendre_p_0_1(n, m, self, bc)
            }
        }
    }

    #[inline]
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>> {
        let (ni, nj) = (n + 1, 2 * m + 1);
        let mut pnm = vec![num_complex::Complex::new(f64::NAN, 0.0); ni * nj];
        unsafe {
            if norm {
                crate::ffi::xsf::assoc_legendre_p_all_1_1(n, m, self, bc, pnm.as_mut_ptr());
            } else {
                crate::ffi::xsf::assoc_legendre_p_all_0_1(n, m, self, bc, pnm.as_mut_ptr());
            }
        }
        utils::vec_into_vecvec(pnm, ni, nj, false)
    }
}

pub trait LegendreQArg: sealed::Sealed + Sized {
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>);
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>);
}

impl LegendreQArg for f64 {
    #[inline]
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
        let mut qn = vec![f64::NAN; n + 1];
        let mut qd = vec![f64::NAN; n + 1];

        unsafe {
            crate::ffi::xsf::lqn(n, self, qn.as_mut_ptr(), qd.as_mut_ptr());
        }
        (qn, qd)
    }

    #[inline]
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>) {
        let (ni, nj) = (m + 1, n + 1);
        let mut qm = vec![f64::NAN; ni * nj];
        let mut qd = vec![f64::NAN; ni * nj];

        unsafe {
            crate::ffi::xsf::lqmn(m, n, self, qm.as_mut_ptr(), qd.as_mut_ptr());
        }

        (
            utils::vec_to_vecvec(&qm, ni, nj, true),
            utils::vec_to_vecvec(&qd, ni, nj, true),
        )
    }
}

impl LegendreQArg for num_complex::Complex<f64> {
    #[inline]
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
        let mut cqn = vec![num_complex::Complex::new(f64::NAN, 0.0); n + 1];
        let mut cqd = vec![num_complex::Complex::new(f64::NAN, 0.0); n + 1];

        unsafe {
            crate::ffi::xsf::lqn_1(n, self, cqn.as_mut_ptr(), cqd.as_mut_ptr());
        }
        (cqn, cqd)
    }

    #[inline]
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>) {
        let (ni, nj) = (m + 1, n + 1);
        let mut cqm = vec![num_complex::Complex::new(f64::NAN, 0.0); ni * nj];
        let mut cqd = vec![num_complex::Complex::new(f64::NAN, 0.0); ni * nj];

        unsafe {
            crate::ffi::xsf::lqmn_1(m, n, self, cqm.as_mut_ptr(), cqd.as_mut_ptr());
        }

        (
            utils::vec_into_vecvec(cqm, ni, nj, true),
            utils::vec_into_vecvec(cqd, ni, nj, true),
        )
    }
}

/// Legendre polynomial of degree `n`
#[doc(alias = "eval_legendre")]
#[inline]
pub fn legendre_p<T: LegendrePArg>(n: i32, z: T) -> T {
    z.legendre_p(n as c_int)
}

/// All Legendre polynomials of the 1st kind
///
/// Output length is `n + 1`. The entry at `j` corresponds to degree `j` in `0..=n`.
#[doc(alias = "lpn", alias = "clpn")]
#[inline]
pub fn legendre_p_all<T: LegendrePArg>(n: usize, z: T) -> Vec<T> {
    z.legendre_p_all(n)
}

/// Spherical Legendre polynomial of degree `n` and order `m`
#[inline]
pub fn sph_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// All spherical Legendre polynomials of the 1st kind
///
/// Output shape is `(n+1, 2m+1)`. The entry at `(j, i)` corresponds to degree `j` in `0..=n` and
/// order `i` in `-m..=m`.
#[inline]
pub fn sph_legendre_p_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.sph_legendre_p_all(n, m)
}

/// Associated Legendre polynomial of the 1st kind
#[doc(alias = "lpmv")]
#[inline]
pub fn assoc_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2, false)
}

/// All associated Legendre polynomials of the 1st kind
#[doc(alias = "lpmn")]
#[inline]
pub fn assoc_legendre_p_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.assoc_legendre_p_all(n, m, 2, false)
}

/// Normalized associated Legendre polynomial of the 1st kind
#[inline]
pub fn assoc_legendre_p_norm<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2, true)
}

/// All normalized associated Legendre polynomials of the 1st kind
#[inline]
pub fn assoc_legendre_p_norm_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.assoc_legendre_p_all(n, m, 2, true)
}

/// All Legendre polynomials of the 2nd kind and their derivatives
///
/// Output lengths are `n + 1`. The entry at `j` corresponds to degree `j` in `0..=n`.
#[doc(alias = "lqn")]
#[inline]
pub fn legendre_q_all<T: LegendreQArg>(n: usize, z: T) -> (Vec<T>, Vec<T>) {
    z.legendre_q_all(n)
}

/// All associated Legendre polynomials of the 2nd kind and their derivatives
///
/// Computes the associated Legendre function of the second kind of order `m` and degree `n`, `Qmn(z)`,
/// and its derivative, `Qmn'(z)`. Returns two arrays of size `(n+1, m+1)` containing `Qmn(z)` and
/// `Qmn'(z)` for all degrees from `0..=n` and orders from `0..=m`.
#[doc(alias = "lqmn")]
#[inline]
pub fn assoc_legendre_q_all<T: LegendreQArg>(
    n: usize,
    m: usize,
    z: T,
) -> (Vec<Vec<T>>, Vec<Vec<T>>) {
    z.assoc_legendre_q_all(n, m)
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    // there are no xsref tables, so we use smoke tests instead

    use core::f64::consts;
    use num_complex::c64;
    const I: num_complex::Complex<f64> = num_complex::Complex::new(0.0, 1.0);
    const LN_3: f64 = 1.098_612_288_668_109_8_f64;

    // legendre_p

    #[test]
    fn test_legendre_p_f64() {
        assert_eq!(
            [
                crate::legendre_p(0, 0.0),
                crate::legendre_p(1, 0.0),
                crate::legendre_p(2, 0.0),
            ],
            [1.0, 0.0, -0.5],
        );
    }

    #[test]
    fn test_legendre_p_c64() {
        assert_eq!(
            [
                crate::legendre_p(0, I),
                crate::legendre_p(1, I),
                crate::legendre_p(2, I),
            ],
            [c64(1.0, 0.0), I, c64(-2.0, 0.0)],
        );
    }

    // legendre_p_all

    #[test]
    fn test_legendre_p_all_f64() {
        assert_eq!(crate::legendre_p_all(2, 0.0), vec![1.0, 0.0, -0.5]);
    }

    #[test]
    fn test_legendre_p_all_c64() {
        assert_eq!(
            crate::legendre_p_all(2, I),
            vec![c64(1.0, 0.0), I, c64(-2.0, 0.0)],
        );
    }

    // sph_legendre_p

    #[test]
    fn test_sph_legendre_p_f64() {
        assert_eq!(
            [
                crate::sph_legendre_p(0, 0, 0.0),
                crate::sph_legendre_p(1, 0, 0.0),
                crate::sph_legendre_p(1, 1, 0.0),
            ],
            [0.282_094_791_773_878_14, 0.488_602_511_902_919_87, 0.0],
        );
    }

    #[test]
    fn test_sph_legendre_p_c64() {
        crate::np_assert_allclose!(
            &[
                crate::sph_legendre_p(0, 0, I),
                crate::sph_legendre_p(1, 0, I),
                crate::sph_legendre_p(1, 1, I),
            ],
            &[
                c64(0.282_094_791_773_878_2, 0.0),
                c64(0.753_953_074_239_480_4, 0.0),
                c64(-0.406_025_136_855_663_4, 0.0),
            ],
            atol = f64::EPSILON
        );
    }

    // sph_legendre_p_all

    #[test]
    fn test_sph_legendre_p_all_f64() {
        assert_eq!(
            crate::sph_legendre_p_all(1, 1, 0.0),
            vec![
                vec![0.282_094_791_773_878_14, 0.0, 0.0],
                vec![0.488_602_511_902_919_87, 0.0, 0.0],
            ]
        );
    }

    #[test]
    fn test_sph_legendre_p_all_c64() {
        let pnm = crate::sph_legendre_p_all(1, 1, I);
        assert_eq!(pnm.len(), 2);

        crate::np_assert_allclose!(
            &pnm[0],
            &[
                c64(0.282_094_791_773_878_2, 0.0),
                c64(0.0, 0.0),
                c64(0.0, 0.0)
            ],
            atol = f64::EPSILON
        );
        crate::np_assert_allclose!(
            &pnm[1],
            &[
                c64(0.753_953_074_239_480_4, 0.0),
                c64(-0.406_025_136_855_663_4, 0.0),
                c64(0.406_025_136_855_663_4, 0.0),
            ],
            atol = f64::EPSILON
        );
    }

    // assoc_legendre_p

    #[test]
    fn test_assoc_legendre_p_f64() {
        assert_eq!(
            [
                crate::assoc_legendre_p(0, 0, 0.0),
                crate::assoc_legendre_p(1, 0, 0.0),
                crate::assoc_legendre_p(0, 1, 0.0),
                crate::assoc_legendre_p(1, 1, 0.0),
            ],
            [1.0, 0.0, 0.0, -1.0],
        );
    }

    #[test]
    fn test_assoc_legendre_p_c64() {
        assert_eq!(
            [
                crate::assoc_legendre_p(0, 0, I),
                crate::assoc_legendre_p(1, 0, I),
                crate::assoc_legendre_p(0, 1, I),
                crate::assoc_legendre_p(1, 1, I),
            ],
            [c64(1.0, 0.0), I, c64(0.0, 0.0), c64(-consts::SQRT_2, 0.0)],
        );
    }

    // assoc_legendre_p_all

    #[test]
    fn test_assoc_legendre_p_all_f64() {
        assert_eq!(
            crate::assoc_legendre_p_all(1, 1, 0.0),
            vec![vec![1.0, 0.0, 0.0], vec![0.0, -1.0, 0.5]],
        );
    }

    #[test]
    fn test_assoc_legendre_p_all_c64() {
        assert_eq!(
            crate::assoc_legendre_p_all(1, 1, I),
            vec![
                vec![c64(1.0, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)],
                vec![
                    c64(0.0, 1.0),
                    c64(-consts::SQRT_2, 0.0),
                    c64(consts::FRAC_1_SQRT_2, 0.0),
                ],
            ]
        );
    }

    // assoc_legendre_p_norm

    #[test]
    fn test_assoc_legendre_p_norm_f64() {
        crate::np_assert_allclose!(
            &[
                crate::assoc_legendre_p_norm(0, 0, 0.0),
                crate::assoc_legendre_p_norm(1, 0, 0.0),
                crate::assoc_legendre_p_norm(0, 1, 0.0),
                crate::assoc_legendre_p_norm(1, 1, 0.0),
            ],
            &[
                consts::FRAC_1_SQRT_2,
                0.0,
                0.0,
                -0.866_025_403_784_438_6_f64, // -sqrt(3) / 2
            ],
            atol = f64::EPSILON
        );
    }

    #[test]
    fn test_assoc_legendre_p_norm_c64() {
        crate::np_assert_allclose!(
            &[
                crate::assoc_legendre_p_norm(0, 0, I),
                crate::assoc_legendre_p_norm(1, 0, I),
                crate::assoc_legendre_p_norm(0, 1, I),
                crate::assoc_legendre_p_norm(1, 1, I),
            ],
            &[
                c64(consts::FRAC_1_SQRT_2, 0.0),
                c64(0.0, 1.224_744_871_391_589),
                c64(0.0, 0.0),
                c64(-1.224_744_871_391_589_2, 0.0),
            ],
            atol = f64::EPSILON
        );
    }

    // legendre_q_all

    #[test]
    fn test_legendre_q_all_f64() {
        let (qn, qd) = crate::legendre_q_all(4, 0.5);

        assert_eq!(qn.len(), 5);
        assert_eq!(qd.len(), 5);

        assert_eq!(qn[..2], [LN_3 * 0.5, LN_3 * 0.25 - 1.0]);
        assert_eq!(qd[..2], [4.0 / 3.0, 2.0 / 3.0 + 0.5 * LN_3]);
    }

    #[test]
    fn test_legendre_q_all_c64() {
        let (qn, qd) = crate::legendre_q_all(4, c64(0.0, 1.0));

        assert_eq!(qn.len(), 5);
        assert_eq!(qd.len(), 5);

        assert_eq!(
            qn[..2],
            [
                c64(0.0, consts::FRAC_PI_4),
                c64(-1.0 - consts::FRAC_PI_4, 0.0),
            ]
        );
    }

    // assoc_legendre_q_all

    #[test]
    fn assoc_test_legendre_q_all_f64() {
        assert_eq!(
            crate::assoc_legendre_q_all(1, 2, 0.0),
            (
                vec![
                    vec![0.0, -1.0, 0.0], // order 0,1,2 for degree 0
                    vec![-1.0, 0.0, 2.0], // order 0,1,2 for degree 1
                ],
                vec![
                    vec![1.0, 0.0, 2.0],  // order 0,1,2 for degree 0
                    vec![0.0, -2.0, 0.0], // order 0,1,2 for degree 1
                ],
            ),
        );
    }

    #[test]
    fn assoc_test_legendre_q_all_c64() {
        assert_eq!(
            crate::assoc_legendre_q_all(1, 2, c64(0.0, 0.0)),
            (
                vec![
                    vec![c64(0.0, 0.0), c64(-1.0, 0.0), c64(0.0, 0.0)],
                    vec![c64(-1.0, 0.0), c64(0.0, 0.0), c64(2.0, 0.0)],
                ],
                vec![
                    vec![c64(1.0, 0.0), c64(0.0, 0.0), c64(2.0, 0.0)],
                    vec![c64(0.0, 0.0), c64(-2.0, 0.0), c64(0.0, 0.0)],
                ],
            ),
        );
    }
}
