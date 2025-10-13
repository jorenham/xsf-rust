use alloc::{vec, vec::Vec};
use core::ffi::c_int;
use num_complex::Complex;

/// Spherical harmonics
pub fn sph_harm_y(n: usize, m: isize, theta: f64, phi: f64) -> Complex<f64> {
    assert!(n <= c_int::MAX as usize);
    assert!(m.abs() <= c_int::MAX as isize);
    unsafe { crate::ffi::xsf::sph_harm_y(n as c_int, m as c_int, theta, phi) }.into()
}

/// All spherical harmonics up to the specified degree `n` and order `m`
///
/// Output shape is `(n + 1, 2 * m + 1)`. The entry at `(j, i)` corresponds to degree `j` and
/// order `i` for all `0 <= j <= n` and `-m <= i <= m`.
///
/// Note: Unlike Python, Rust uses only positive array indices. The mapping from array index
/// to spherical harmonic order is:
/// - Index `0` to `m`: orders `0, 1, 2, ..., m`
/// - Index `m+1` to `2*m`: orders `-m, -(m-1), ..., -1`
pub fn sph_harm_y_all(n: usize, m: usize, theta: f64, phi: f64) -> Vec<Vec<Complex<f64>>> {
    let (nr, nc) = (n + 1, 2 * m + 1);
    let mut out = vec![f64::NAN.into(); nr * nc];
    unsafe { crate::ffi::xsf::sph_harm_y_all(n, m, theta, phi, out.as_mut_ptr()) };
    crate::utils::vec_into_vecvec(out, nr, nc, false)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;
    use num_complex::c64;

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

        let y_all = crate::sph_harm_y_all(n_max, m_max, theta, phi);

        // Check shape
        assert_eq!(y_all.len(), n_max + 1);
        assert_eq!(y_all[0].len(), 2 * m_max + 1);

        fn col_to_order(col: usize, m_max: usize) -> isize {
            if col <= m_max {
                col as isize
            } else {
                -((2 * m_max + 1 - col) as isize)
            }
        }

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
            vec![c64(0.28209479, 0.0), c64(0.0, 0.0), c64(0.0, 0.0)],
            vec![
                c64(0.34549415, 0.0),
                c64(-0.21157109, -0.12215063),
                c64(0.21157109, -0.12215063),
            ],
            vec![
                c64(0.15769578, 0.0),
                c64(-0.33452327, -0.1931371),
                c64(0.33452327, -0.1931371),
            ],
            vec![
                c64(-0.13193776, 0.0),
                c64(-0.29685995, -0.17139217),
                c64(0.29685995, -0.17139217),
            ],
        ];

        assert_eq!(y_all.len(), expected.len());

        for (row_actual, row_expected) in y_all.iter().zip(expected.iter()) {
            crate::np_assert_allclose!(&row_actual, &row_expected, atol = ATOL);
        }
    }
}
