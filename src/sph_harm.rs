use crate::{bindings, utils};
use alloc::vec::Vec;
use core::ffi::c_int;
use num_complex::Complex;

/// Spherical harmonics
pub fn sph_harm_y(n: usize, m: isize, theta: f64, phi: f64) -> Complex<f64> {
    assert!(n <= c_int::MAX as usize);
    assert!(m.abs() <= c_int::MAX as isize);
    unsafe { bindings::sph_harm_y(n as c_int, m as c_int, theta, phi) }.into()
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
    let mut res = bindings::complex_zeros(nr * nc);
    unsafe { bindings::sph_harm_y_all(n, m, theta, phi, res.as_mut_ptr()) };
    utils::vec_to_vecvec(bindings::cvec_into(res), nr, nc, false)
}

#[cfg(test)]
mod tests {
    use super::*;
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

        let y_all = sph_harm_y_all(n_max, m_max, theta, phi);

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
            for (col_idx, &y_from_all) in y_row.iter().enumerate() {
                let m = col_to_order(col_idx, m_max);
                let expected = sph_harm_y(n, m, theta, phi);
                let error = (expected - y_from_all).norm();
                let tolerance = ATOL + RTOL * expected.norm().max(y_from_all.norm());
                assert!(error <= tolerance);
            }
        }
    }

    // Basic smoke tests for functionality
    #[test]
    fn test_sph_harm_y_basic() {
        // Test Y_0^0 at theta=0, phi=0 should be sqrt(1/(4*pi))
        let y00 = sph_harm_y(0, 0, 0.0, 0.0);
        let expected = c64((4.0 * PI).sqrt().recip(), 0.0);
        assert!((y00 - expected).norm() < 1e-14);
    }

    #[test]
    fn test_sph_harm_y_all_basic() {
        let y_all = sph_harm_y_all(1, 1, PI / 2.0, 0.0);

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
        let y_all = sph_harm_y_all(3, 1, PI / 4.0, PI / 6.0);

        // Expected scipy output with reasonable tolerance for floating-point differences
        let expected = [
            vec![(0.28209479, 0.0), (0.0, 0.0), (0.0, 0.0)],
            vec![
                (0.34549415, 0.0),
                (-0.21157109, -0.12215063),
                (0.21157109, -0.12215063),
            ],
            vec![
                (0.15769578, 0.0),
                (-0.33452327, -0.1931371),
                (0.33452327, -0.1931371),
            ],
            vec![
                (-0.13193776, 0.0),
                (-0.29685995, -0.17139217),
                (0.29685995, -0.17139217),
            ],
        ];

        assert_eq!(y_all.len(), expected.len());

        for (rust_row, expected_row) in y_all.iter().zip(expected.iter()) {
            assert_eq!(rust_row.len(), expected_row.len());

            for (rust_val, &(exp_re, exp_im)) in rust_row.iter().zip(expected_row.iter()) {
                let expected_val = Complex::new(exp_re, exp_im);
                let diff = (rust_val - expected_val).norm();
                assert!(diff < ATOL);
            }
        }
    }
}
