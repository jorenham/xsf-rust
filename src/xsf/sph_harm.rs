use core::ffi::c_int;
use num_complex::Complex64;

/// Spherical harmonics
///
/// Corresponds to [`scipy.special.sph_harm_y`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html
///
/// # Definition
///
/// $$
/// Y_n^m(\theta, \phi) =
/// \sqrt{ {2n + 1 \over 4\pi} {(n - m)! \over (n + m)!} } \\
/// P_n^m(\cos \theta) \\
/// e^{i m \phi}
/// $$
///
/// where $P_n^m$ are the (unnormalized) associated Legendre polynomials.
///
/// With SciPy's convention, the first several spherical harmonics are
///
/// $$
/// \begin{align*}
/// Y_0^0(\theta, \phi)    &=  \sqrt{1 \over 4\pi} \\\\
/// Y_1^{-1}(\theta, \phi) &=  \sqrt{3 \over 8\pi} \\ \sin \theta \\ e^{-i \phi} \\\\
/// Y_1^0(\theta, \phi)    &=  \sqrt{3 \over 4\pi} \\ \cos \theta \\\\
/// Y_1^1(\theta, \phi)    &= -\sqrt{3 \over 8\pi} \\ \sin \theta \\ e^{ i \phi}
/// \end{align*}
/// $$
///
/// # Panics
/// - if `n` is greater than [`c_int::MAX`](core::ffi::c_int::MAX)
/// - if `m` is greater than [`c_int::MAX`](core::ffi::c_int::MAX)
///
/// # See also
/// - [`sph_harm_y_all`]
/// - [`assoc_legendre_p`](crate::assoc_legendre_p)
#[must_use]
#[inline]
pub fn sph_harm_y(n: usize, m: isize, theta: f64, phi: f64) -> Complex64 {
    assert!(n <= c_int::MAX as usize);
    assert!(m.abs() <= c_int::MAX as isize);
    unsafe { crate::ffi::xsf::sph_harm_y(n.try_into().unwrap(), m.try_into().unwrap(), theta, phi) }
}

/// All spherical harmonics up to the specified degree $n$ and order $m$
///
/// Corresponds to [`scipy.special.sph_harm_y_all`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y_all.html
///
/// # Notes
///
/// Output shape is `(n + 1, 2 * m + 1)`. The entry at `(j, i)` corresponds to degree $j$ and
/// order $i$ for all $0 \leq j \leq n$ and $-m \leq i \leq m$.
///
/// Note: Unlike Python, Rust uses only positive array indices. The mapping from array index
/// to spherical harmonic order is:
/// - Index $0$ to $m$: orders $0, 1, 2, \ldots, m$
/// - Index $m+1$ to $2*m$: orders $-m, -(m-1), \ldots, -1$
///
/// # See also
/// - [`sph_harm_y`]
#[must_use]
#[inline]
pub fn sph_harm_y_all(n: usize, m: usize, theta: f64, phi: f64) -> Vec<Vec<Complex64>> {
    let (nr, nc) = (n + 1, 2 * m + 1);
    let mut out = vec![f64::NAN.into(); nr * nc];
    unsafe { crate::ffi::xsf::sph_harm_y_all(n, m, theta, phi, out.as_mut_ptr()) };
    crate::utils::vec_into_vecvec(out, nr, nc, false)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;
    use num_complex::c64;
    use num_traits::ToPrimitive;

    /// Ported from `scipy.special.tests.test_basic.TestSphericalHarmonics.test_all`
    #[test]
    fn test_sph_harm_y_vs_all() {
        let cases = [
            [7, 1],
            [7, 4],
            [7, 5],
            [10, 1],
            [10, 4],
            [10, 5],
            [10, 9],
            [50, 1],
            [50, 4],
            [50, 5],
            [50, 9],
            [50, 14],
        ];
        let args = [
            (0.0, 0.0),
            (PI / 4.0, PI / 6.0),
            (PI / 2.0, PI / 3.0),
            (3.0 * PI / 4.0, PI),
            (PI, 2.0 * PI / 3.0),
        ];

        for [n_max, m_max] in cases {
            for (theta, phi) in &args {
                test_sph_harm_y_vs_all_single(n_max, m_max, *theta, *phi);
            }
        }
    }

    fn test_sph_harm_y_vs_all_single(n_max: usize, m_max: usize, theta: f64, phi: f64) {
        const RTOL: f64 = 1e-5;
        const ATOL: f64 = 1e-8;

        fn col_to_order(col: usize, m_max: usize) -> isize {
            if col <= m_max {
                col.to_isize().unwrap()
            } else {
                -(2 * m_max + 1 - col).to_isize().unwrap()
            }
        }

        let y_all = crate::sph_harm_y_all(n_max, m_max, theta, phi);

        // Check shape
        assert_eq!(y_all.len(), n_max + 1);
        assert_eq!(y_all[0].len(), 2 * m_max + 1);

        // Check each entry matches the corresponding individual call
        for (n, y_row) in y_all.iter().enumerate().take(n_max + 1) {
            let expected_row: Vec<_> = (0..=2 * m_max)
                .map(|col_idx| crate::sph_harm_y(n, col_to_order(col_idx, m_max), theta, phi))
                .collect();
            crate::np_assert_allclose!(y_row, &expected_row, rtol = RTOL, atol = ATOL);
        }
    }

    // Basic smoke tests for functionality
    #[test]
    fn test_sph_harm_y_basic() {
        // Test Y_0^0 at theta=0, phi=0 should be sqrt(1/(4*pi))
        let y00 = crate::sph_harm_y(0, 0, 0.0, 0.0);
        let expected = c64((4.0 * PI).sqrt().recip(), 0.0);
        assert!((y00 - expected).norm() < 1e-14);
    }

    #[test]
    fn test_sph_harm_y_all_basic() {
        let y_all = crate::sph_harm_y_all(1, 1, PI / 2.0, 0.0);

        // Check shape: (n+1, 2*m+1) = (2, 3)
        assert_eq!(y_all.len(), 2);
        assert_eq!(y_all[0].len(), 3);
        assert_eq!(y_all[1].len(), 3);

        // Y_0^0 should be non-zero and is at index 0
        assert!(y_all[0][0].norm() > 0.0);
    }

    #[test]
    fn test_scipy_compatibility() {
        const ATOL: f64 = 1.5e-8;

        // Verify that our implementation matches scipy.special.sph_harm_y_all
        let y_all = crate::sph_harm_y_all(3, 1, PI / 4.0, PI / 6.0);

        // Expected scipy output with reasonable tolerance for floating-point differences
        let expected = [
            vec![c64(0.282_094_79, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)],
            vec![
                c64(0.345_494_15, 0.0),
                c64(-0.211_571_09, -0.122_150_63),
                c64(0.211_571_09, -0.122_150_63),
            ],
            vec![
                c64(0.157_695_78, 0.0),
                c64(-0.334_523_27, -0.193_137_1),
                c64(0.334_523_27, -0.193_137_1),
            ],
            vec![
                c64(-0.131_937_76, 0.0),
                c64(-0.296_859_95, -0.171_392_17),
                c64(0.296_859_95, -0.171_392_17),
            ],
        ];

        assert_eq!(y_all.len(), expected.len());

        for (row_actual, row_expected) in y_all.iter().zip(expected.iter()) {
            crate::np_assert_allclose!(&row_actual, &row_expected, atol = ATOL);
        }
    }
}
