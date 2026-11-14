//! Eigenvalues and eigenvectors of a real symmetric tridiagonal matrix `T`, given its diagonal `d`
//! and off-diagonal `e` (`e.len() + 1 == d.len()`)
//!
//! Like in LAPACK, `T` is first split into unreduced blocks, which are then scaled separately.

use core::iter::once;
use core::ops::Range;

use crate::utils::EPSILON_SQ;
use crate::xsf::cephes::dd_real::{frexp_f64, ldexp_f64};

/// Maximum number of QL iterations per eigenvalue, like in EISPACK's `tql1`
const MAX_ITER: usize = 30;

/// The largest magnitude of the entries
fn max_abs(d: &[f64], e: &[f64]) -> f64 {
    d.iter().chain(e).fold(0.0, |m, x| x.abs().max(m))
}

/// The diagonal and off-diagonal of the block `b` of `T`
fn block<'a>(d: &'a [f64], e: &'a [f64], b: &Range<usize>) -> (&'a [f64], &'a [f64]) {
    (&d[b.clone()], &e[b.start..b.end - 1])
}

/// The ranges of the unreduced blocks of `T`
///
/// Like in LAPACK's `dsteqr`, it is split where `|e_i| <= eps sqrt(|d_i| |d_{i + 1}|)`, which is
/// relative, so that the small eigenvalues of graded matrices stay accurate.
fn blocks(d: &[f64], e: &[f64]) -> Vec<Range<usize>> {
    let mut parts = Vec::new();
    let mut start = 0;
    for (i, ei) in e.iter().enumerate() {
        if ei.abs() <= f64::EPSILON * d[i].abs().sqrt() * d[i + 1].abs().sqrt() {
            parts.push(start..i + 1);
            start = i + 1;
        }
    }
    parts.push(start..d.len());
    parts
}

/// `d` and `e` scaled by `2^-k`, and `k`, so that their largest magnitude is in
/// `[2^-405, 2^509)`, as in LAPACK's `dsteqr`
///
/// This avoids over- and underflow, e.g. of `e^2`, without needlessly pushing small entries
/// towards underflow.
fn scaled(d: &[f64], e: &[f64]) -> (Vec<f64>, Vec<f64>, i32) {
    let (_, m) = frexp_f64(max_abs(d, e));
    let k = m - m.clamp(-404, 509);
    let scale = |x: &[f64]| x.iter().map(|&xi| ldexp_f64(xi, -k)).collect();
    (scale(d), scale(e), k)
}

/// The eigenvalues of an unreduced block scaled by `2^-k`, unsorted, and `k`, by the implicit QL
/// method with Wilkinson shifts (EISPACK's `tql1` and `tql2`, LAPACK's `dsteqr`), or `None` if
/// it doesn't converge
///
/// The rotations are applied to the rows of the row-major `n x n` matrix `z` (if not empty), so
/// that for `z = I`, its column `j` becomes the eigenvector of eigenvalue `j`.
fn ql(d: &[f64], e: &[f64], z: &mut [f64]) -> Option<(Vec<f64>, i32)> {
    let n = d.len();
    let (mut d, e, k) = scaled(d, e);
    // `e[i]` couples `d[i]` and `d[i + 1]`, and `e[n - 1] = 0`
    let mut e: Vec<f64> = e.into_iter().chain(once(0.0)).collect();
    for l in 0..n {
        for iter in 0.. {
            // split off the block `l..=m` at the first negligible `e[m]`, relative like in
            // `blocks`, or if `e[m]^2` underflows (as in `dsteqr`)
            let m = (l..n - 1)
                .find(|&m| {
                    e[m] * e[m] <= EPSILON_SQ * d[m].abs() * d[m + 1].abs() + f64::MIN_POSITIVE
                })
                .unwrap_or(n - 1);
            if m == l {
                break;
            }
            if iter == MAX_ITER {
                return None;
            }

            // the Wilkinson shift is the eigenvalue of the top 2x2 block closest to `d[l]`
            let g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut g = d[m] - d[l] + e[l] / (g + g.hypot(1.0).copysign(g));

            // chase the bulge up from `m` to `l` with Givens rotations
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut deflated = false;
            for i in (l..m).rev() {
                let f = s * e[i];
                let b = c * e[i];
                let r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    // `e[i + 1]` underflowed, so restart on the split blocks
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    deflated = true;
                    break;
                }
                s = f / r;
                c = g / r;
                let g1 = d[i + 1] - p;
                let r = (d[i] - g1) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g1 + p;
                g = c * r - b;
                for row in z.chunks_exact_mut(n) {
                    let zi1 = row[i + 1];
                    row[i + 1] = s * row[i] + c * zi1;
                    row[i] = c * row[i] - s * zi1;
                }
            }
            if !deflated {
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        }
    }
    Some((d, k))
}

/// Whether all entries are finite, after checking that the sizes match
fn check(d: &[f64], e: &[f64]) -> bool {
    assert!(e.len() + 1 == d.len(), "invalid matrix size");
    d.iter().chain(e).all(|x| x.is_finite())
}

/// All eigenvalues in ascending order
///
/// They are NaN if an entry isn't finite, or if the QL iteration doesn't converge.
///
/// # Panics
/// - if `e.len() + 1 != d.len()`
pub(crate) fn eigvals(d: &[f64], e: &[f64]) -> Vec<f64> {
    let n = d.len();
    if !check(d, e) {
        return vec![f64::NAN; n];
    }
    let mut w = Vec::with_capacity(n);
    for b in blocks(d, e) {
        let (db, eb) = block(d, e, &b);
        let Some((wb, k)) = ql(db, eb, &mut []) else {
            return vec![f64::NAN; n];
        };
        w.extend(wb.into_iter().map(|x| ldexp_f64(x, k)));
    }
    w.sort_by(f64::total_cmp);
    w
}

/// A unit eigenvector, for the eigenvalue closest to `lambda`
///
/// It is NaN if an entry or `lambda` isn't finite, or if the QL iteration doesn't converge.
///
/// # Panics
/// - if `e.len() + 1 != d.len()`
pub(crate) fn eigvec(d: &[f64], e: &[f64], lambda: f64) -> Vec<f64> {
    let n = d.len();
    let mut x = vec![f64::NAN; n];
    if !check(d, e) {
        return x;
    }
    let mut closest = f64::INFINITY;
    let mut best = None;
    for b in blocks(d, e) {
        let (db, eb) = block(d, e, &b);
        let m = b.len();
        let mut z = vec![0.0; m * m];
        z.iter_mut().step_by(m + 1).for_each(|zii| *zii = 1.0);
        let Some((w, k)) = ql(db, eb, &mut z) else {
            return vec![f64::NAN; n];
        };
        // a quarter of the distance, which can't overflow, since `|w| <= 3 max(|d|, |e|)`
        for (j, wj) in w.into_iter().enumerate() {
            let distance = (ldexp_f64(wj, k - 2) - 0.25 * lambda).abs();
            if distance < closest {
                closest = distance;
                best = Some((b.clone(), wj));
                x.fill(0.0);
                for (xk, row) in x[b.clone()].iter_mut().zip(z.chunks_exact(m)) {
                    *xk = row[j];
                }
            }
        }
    }
    if let Some((b, w)) = best {
        let (db, eb) = block(d, e, &b);
        let (db, eb, _) = scaled(db, eb);
        refine(&db, &eb, w, &mut x[b]);
    }
    x
}

/// The unit eigenvector `x` of eigenvalue `w` of an unreduced block after two steps of inverse
/// iteration (like LAPACK's `dstein`), since the small components of the QL eigenvector are only
/// accurate relative to its norm
///
/// `(T - w I) y = x` is solved with partial pivoting (as in `dgtsv`), where a zero pivot becomes
/// `eps max(|d|, |e|)`, and a step that overflows is discarded.
fn refine(d: &[f64], e: &[f64], w: f64, x: &mut [f64]) {
    let n = d.len();
    let tiny = f64::EPSILON * max_abs(d, e);
    for _ in 0..2 {
        // `U` with diagonal `u`, and superdiagonals `u1` and `u2`
        let mut u: Vec<f64> = d.iter().map(|di| di - w).collect();
        let mut u1 = e.to_vec();
        let mut u2 = vec![0.0; n.saturating_sub(2)];
        let mut y = x.to_vec();
        for (i, &l) in e.iter().enumerate() {
            if u[i].abs() >= l.abs() {
                let f = l / u[i];
                u[i + 1] -= f * u1[i];
                y[i + 1] -= f * y[i];
            } else {
                // swap the rows `i` and `i + 1`
                let f = u[i] / l;
                (u[i], u[i + 1], u1[i]) = (l, u1[i] - f * u[i + 1], u[i + 1]);
                if i + 2 < n {
                    (u2[i], u1[i + 1]) = (u1[i + 1], -f * u1[i + 1]);
                }
                y.swap(i, i + 1);
                y[i + 1] -= f * y[i];
            }
        }
        for i in (0..n).rev() {
            let above = if i + 1 < n { u1[i] * y[i + 1] } else { 0.0 };
            let above2 = if i + 2 < n { u2[i] * y[i + 2] } else { 0.0 };
            let pivot = if u[i] == 0.0 { tiny } else { u[i] };
            y[i] = (y[i] - above - above2) / pivot;
        }
        let max = y.iter().fold(0.0, |m, yi| yi.abs().max(m));
        if !max.is_finite() {
            return;
        }
        let norm = y
            .iter()
            .map(|yi| (yi / max) * (yi / max))
            .sum::<f64>()
            .sqrt()
            * max;
        for (xi, yi) in x.iter_mut().zip(y) {
            *xi = yi / norm;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const N: usize = 9;

    /// The second-difference matrix `tridiag(-1, 2, -1)` of size `N`, with eigenvalues
    /// `2 - 2 cos(k pi / (N + 1))` and eigenvectors `sin(j k pi / (N + 1))`
    fn laplacian() -> ([f64; N], [f64; N - 1]) {
        ([2.0; N], [-1.0; N - 1])
    }

    #[test]
    fn test_eigvals_laplacian() {
        let (d, e) = laplacian();
        let actual = eigvals(&d, &e);
        for (k, lambda) in actual.into_iter().enumerate() {
            let expected = 2.0 - 2.0 * ((k + 1) as f64 * PI / (N + 1) as f64).cos();
            assert!(
                (lambda - expected).abs() <= 1e-15,
                "{k}: {lambda} != {expected}"
            );
        }
    }

    #[test]
    fn test_eigvec_laplacian() {
        let (d, e) = laplacian();
        for k in 0..N {
            let x = eigvec(&d, &e, eigvals(&d, &e)[k]);
            let theta = (k + 1) as f64 * PI / (N + 1) as f64;
            let y: Vec<f64> = (1..=N).map(|j| (j as f64 * theta).sin()).collect();
            let norm = y.iter().map(|yj| yj * yj).sum::<f64>().sqrt();
            let sign = x[0].signum() * y[0].signum();
            for (xj, yj) in x.iter().zip(&y) {
                assert!((xj - sign * yj / norm).abs() <= 4e-15, "{k}: {x:?}");
            }
        }
    }

    #[test]
    fn test_eigvals_legendre() {
        // Golub-Welsch: the roots `0, +-sqrt(3/5)` of the Legendre polynomial `P_3`
        let e = [1.0 / 3.0_f64.sqrt(), 2.0 / 15.0_f64.sqrt()];
        let r = (3.0_f64 / 5.0).sqrt();
        let actual = eigvals(&[0.0; 3], &e);
        for (x, y) in actual.iter().zip([-r, 0.0, r]) {
            assert!((x - y).abs() <= 2e-16, "{actual:?}");
        }
    }

    #[test]
    fn test_eigvals_diagonal() {
        // unsorted, and already split into 1x1 blocks
        assert_eq!(eigvals(&[3.0, 0.0, -1.0], &[0.0, 0.0]), [-1.0, 0.0, 3.0]);
        assert_eq!(eigvals(&[5.0], &[]), [5.0]);
    }

    /// `max |(T - lambda I) x|`, relative to `max(|d|, |e|)`, for a unit `x`
    fn residual(d: &[f64], e: &[f64], lambda: f64, x: &[f64]) -> f64 {
        let norm = x.iter().map(|xi| xi * xi).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() <= 1e-15, "{x:?}");
        let anorm = max_abs(d, e);
        (0..d.len())
            .map(|i| {
                let below = if i > 0 { e[i - 1] * x[i - 1] } else { 0.0 };
                let above = if i + 1 < d.len() {
                    e[i] * x[i + 1]
                } else {
                    0.0
                };
                (below + (d[i] - lambda) * x[i] + above).abs()
            })
            .fold(0.0, f64::max)
            / anorm
    }

    #[test]
    fn test_eigvec_residual() {
        let cases: [(Vec<f64>, Vec<f64>, f64); 8] = [
            // `d - lambda = 0`
            (vec![1.0; 3], vec![2.0; 2], 1.0),
            // `e^2` underflows without scaling
            (vec![0.0; 3], vec![1e-200; 2], 2.0_f64.sqrt() * 1e-200),
            // the eigenvalues `1 +- 1e-20` are equal in `f64`
            (vec![1.0; 2], vec![1e-20], 1.0),
            // `lambda = 0` is a multiple eigenvalue, of the blocks split off at the tiny `e_i`
            (vec![0.0; 5], vec![1e-320, 1e-320, 1e-150, 1.0], 0.0),
            (
                vec![0.0; 7],
                vec![1e-150, 1e-320, 1.0, 1e-300, 1.0, 1e-150],
                0.0,
            ),
            (
                vec![0.0; 7],
                vec![1e-150, 1.0, 1e-300, 1.0, 1e-320, 1e-150],
                0.0,
            ),
            // `lambda = 1` is a double eigenvalue, up to `3e-16`
            (vec![0.0; 5], vec![3e-16, 1.0, 3e-16, 1.0], 1.0),
            (vec![0.0; 5], vec![3e-16, 1.0, 3e-16, 1.0], -1.0),
        ];
        for (d, e, lambda) in cases {
            let x = eigvec(&d, &e, lambda);
            assert!(residual(&d, &e, lambda, &x) <= 1e-15, "{d:?} {e:?}: {x:?}");
        }
    }

    #[test]
    fn test_eigvec_max() {
        // the eigenvalue `(1 - 1e-17) f64::MAX`, and the residual of `T / 2^1024`
        let d = [-4.128_240_709_468_263_5e307, 4.128_240_709_468_263_5e307];
        let e = [1.749_650_506_123_434_2e308];
        let x = eigvec(&d, &e, f64::MAX);
        let scale = |v: &[f64]| v.iter().map(|&vi| ldexp_f64(vi, -1024)).collect::<Vec<_>>();
        let lambda = ldexp_f64(f64::MAX, -1024);
        assert!(
            residual(&scale(&d), &scale(&e), lambda, &x) <= 1e-15,
            "{x:?}"
        );
    }

    #[test]
    fn test_eigvec_scales() {
        // blocks of very different scales, whose distances to `lambda` overflow when unscaled
        assert_eq!(eigvec(&[1e-300, -1e308], &[0.0], 1e308), [1.0, 0.0]);
        assert_eq!(eigvec(&[1e-300, -1e308], &[0.0], -1e308), [0.0, 1.0]);
    }

    #[test]
    fn test_eigvec_nan() {
        assert!(
            eigvec(&[1.0, 1.0], &[1.0], f64::NAN)
                .into_iter()
                .all(f64::is_nan)
        );
        assert!(
            eigvec(&[1.0, f64::NAN], &[1.0], 1.0)
                .into_iter()
                .all(f64::is_nan)
        );
    }

    /// `eigvals(d, e)` to within 2 ulp
    fn assert_eigvals(d: &[f64], e: &[f64], expected: &[f64]) {
        let actual = eigvals(d, e);
        for (a, b) in actual.iter().zip(expected) {
            assert!((a - b).abs() <= 2.0 * f64::EPSILON * b.abs(), "{actual:?}");
        }
    }

    #[test]
    fn test_eigvals_split() {
        // each block is scaled separately
        assert_eigvals(&[1e308, 1e-100], &[0.0], &[1e-100, 1e308]);
        assert_eigvals(
            &[1e308, 0.0, 0.0],
            &[0.0, 1e-100],
            &[-1e-100, 1e-100, 1e308],
        );
    }

    #[test]
    fn test_graded() {
        // `e_0` isn't negligible relative to `d_0`, with `det(T) = 1e-40`
        let (d, e) = ([2e-40, 1.0], [1e-20]);
        assert_eigvals(&d, &e, &[1e-40, 1.0]);
        let x = eigvec(&d, &e, 1e-40);
        assert!((x[1] / x[0] + 1e-20).abs() <= 1e-35, "{x:?}");
    }

    #[test]
    fn test_eigvals_unscaled() {
        // scaling to `[0.5, 1)` would underflow `e^2`, and deflate the `0` block
        assert_eigvals(
            &[0.0, 0.0, 1e150],
            &[1e-150, 1e-150],
            &[-1e-150, 1e-150, 1e150],
        );
    }

    #[test]
    fn test_eigvals_extreme() {
        // over- and underflow without scaling
        assert_eigvals(&[1e308; 2], &[1e307], &[9e307, 1.1e308]);
        assert_eigvals(&[0.0; 2], &[1.5e308], &[-1.5e308, 1.5e308]);
        assert_eigvals(&[0.0; 2], &[5e-324], &[-5e-324, 5e-324]);
    }

    #[test]
    fn test_eigvals_nan() {
        for (d, e) in [([1.0, f64::NAN], 1.0), ([1.0, 1.0], f64::INFINITY)] {
            assert!(eigvals(&d, &[e]).into_iter().all(f64::is_nan));
        }
    }
}
