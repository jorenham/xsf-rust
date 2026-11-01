//! Translated into pure Rust from `xsf/sph_harm.h` (xsf v0.2.2).

use core::ffi::c_int;
use num_complex::{Complex64, c64};

use crate::xsf::exp::cexp;
use crate::xsf::legendre::{Mdspan, int, sph_legendre_p_for_each_n, sph_legendre_p_for_each_n_m};

fn sph_harm_y_next(m: i32, phi: f64, p: f64, res: &mut Complex64) {
    *res = cexp(c64(0.0, 1.0) * f64::from(m) * phi) * p;
}

fn sph_harm_y_for_each_n(
    n: i32,
    m: i32,
    theta: f64,
    phi: f64,
    res: &mut Complex64,
    mut f: impl FnMut(i32, i32, &Complex64),
) {
    let mut p = [0.0; 2];
    sph_legendre_p_for_each_n(n, m, theta, &mut p, |n, p| {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

fn sph_harm_y_for_each_n_m(
    n: i32,
    m: i32,
    theta: f64,
    phi: f64,
    res: &mut Complex64,
    mut f: impl FnMut(i32, i32, &Complex64),
) {
    let mut p = [0.0; 2];
    sph_legendre_p_for_each_n_m(n, m, theta, &mut p, |n, m, p| {
        sph_harm_y_next(m, phi, p[1], res);

        f(n, m, res);
    });
}

fn sph_harm_y_t(n: i32, m: i32, theta: f64, phi: f64) -> Complex64 {
    let mut res_n = c64(0.0, 0.0);
    sph_harm_y_for_each_n(n, m, theta, phi, &mut res_n, |_, _, _| {});

    res_n
}

fn sph_harm_y_all_t(theta: f64, phi: f64, mut res: Mdspan<'_, Complex64, 2>) {
    let n_max = int(res.extent(0)) - 1;
    let m_max = (int(res.extent(1)) - 1) / 2;

    let mut res_n_m = c64(0.0, 0.0);
    sph_harm_y_for_each_n_m(n_max, m_max, theta, phi, &mut res_n_m, |n, m, res_n_m| {
        if m >= 0 {
            res[(n, m)] = *res_n_m;
        } else {
            res[(n, m + 2 * m_max + 1)] = *res_n_m;
        }
    });
}

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
    sph_harm_y_t(n.try_into().unwrap(), m.try_into().unwrap(), theta, phi)
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
    sph_harm_y_all_t(theta, phi, Mdspan::new(&mut out, [nr, nc]));
    crate::utils::vec_into_vecvec(out, nr, nc, false)
}
