//! Translated into pure Rust from `xsf/cpu/legendre.h` (scipy/xsf@a0c50d0) and `xsf/cpu/recur.h`
//! (scipy/xsf@b1eb270), for `T` either `double` or `std::complex<double>`. The dual numbers of
//! `xsf/cpu/dual.h`, which xsf uses for derivatives, are not translated.
//!
//! Unlike in xsf:
//! - The integer products in the (normalized) recurrence coefficients are computed in `f64`
//!   instead of `int` (scipy/scipy#26574).
//! - Other `int` arithmetic that can overflow wraps (like C++ does in practice), instead of
//!   panicking in debug builds.
//! - The real `lqmn` uses `log(|x| - 1)` (like the complex `lqmn`) instead of `log(x - 1)`
//!   (scipy/xsf#368).
//! - `lqmn` keeps the scale of each value of its rescaled backward recurrences apart, instead of
//!   rescaling the stored values, which can become subnormal (scipy/xsf#571).
//! - The complex `lqn` and `lqmn` take the imaginary parts for a tiny `Im z` from a lifted one,
//!   since those of their seeds can be subnormal.
//! - The complex `sin` and `cos` (for complex `theta`) follow glibc only for finite arguments.

use core::ops::{Index, IndexMut};

use core::ffi::c_int;
use num_complex::Complex64;

use crate::utils;
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble, frexp_f64, ldexp_f64};
use crate::xsf::complex::Num;
use crate::xsf::specfun::{SPECFUN_BIG, SPECFUN_BIG_LOG2, SPECFUN_INF};

/// A mutable view of a row-major array, like `std::mdspan` with `R` extents and `int` indices
pub(crate) struct Mdspan<'a, T, const R: usize> {
    data: &'a mut [T],
    extents: [usize; R],
}

impl<'a, T, const R: usize> Mdspan<'a, T, R> {
    pub(crate) fn new(data: &'a mut [T], extents: [usize; R]) -> Self {
        assert_eq!(data.len(), extents.iter().product::<usize>());
        Self { data, extents }
    }

    pub(crate) fn extent(&self, r: usize) -> usize {
        self.extents[r]
    }

    fn size(&self) -> usize {
        self.data.len()
    }
}

impl<T> Index<i32> for Mdspan<'_, T, 1> {
    type Output = T;
    fn index(&self, i: i32) -> &T {
        &self.data[i as usize]
    }
}

impl<T> IndexMut<i32> for Mdspan<'_, T, 1> {
    fn index_mut(&mut self, i: i32) -> &mut T {
        &mut self.data[i as usize]
    }
}

impl<T> Index<(i32, i32)> for Mdspan<'_, T, 2> {
    type Output = T;
    fn index(&self, (i, j): (i32, i32)) -> &T {
        &self.data[i as usize * self.extents[1] + j as usize]
    }
}

impl<T> IndexMut<(i32, i32)> for Mdspan<'_, T, 2> {
    fn index_mut(&mut self, (i, j): (i32, i32)) -> &mut T {
        &mut self.data[i as usize * self.extents[1] + j as usize]
    }
}

#[allow(clippy::cast_possible_wrap)]
pub(crate) const fn int(x: usize) -> i32 {
    x as i32
}

/// The `int` product `static_cast<T>(a * b * ...)`, computed in `i64` (instead of `int`) to avoid
/// overflow, so that, unlike a product of `f64`, a zero product is `+0.0`
fn int_product(factors: &[i32]) -> f64 {
    factors
        .iter()
        .fold(1_i64, |p, &x| p.wrapping_mul(i64::from(x))) as f64
}

// recur.h

fn dot<T: Num, const N: usize>(x: &[T; N], y: &[T; N]) -> T {
    let mut res = T::from(0.0);
    for i in 0..N {
        res = res + x[i] * y[i];
    }

    res
}

fn forward_recur_shift_left<T: Copy, const K: usize>(res: &mut [T; K]) {
    for k in 1..K {
        res[k - 1] = res[k];
    }
}

fn forward_recur_rotate_left<T: Copy, const K: usize>(res: &mut [T; K]) {
    let tmp = res[0];
    forward_recur_shift_left(res);
    res[K - 1] = tmp;
}

/// Compute a forward recurrence that depends on K previous values.
///
/// - `first`: begin iterator
/// - `last`: end iterator
/// - `r`: recurrence
/// - `res`: values, initialised to the leading K values
/// - `f`: a function to be called as `f(it, res)`
fn forward_recur<T: Num, const K: usize>(
    first: i32,
    last: i32,
    r: impl Fn(i32, &mut [T; K]),
    res: &mut [T; K],
    mut f: impl FnMut(i32, &[T; K]),
) {
    let mut it = first;
    while it.wrapping_sub(first) != int(K) && it != last {
        forward_recur_rotate_left(res);

        f(it, res);
        it = it.wrapping_add(1);
    }

    if last.wrapping_sub(first) > int(K) {
        while it != last {
            let mut coef = [T::from(0.0); K];
            r(it, &mut coef);

            let tmp = dot(&coef, res);
            forward_recur_shift_left(res);
            res[K - 1] = tmp;

            f(it, res);
            it = it.wrapping_add(1);
        }
    }
}

fn backward_recur<T: Num, const K: usize>(
    first: i32,
    last: i32,
    r: impl Fn(i32, &mut [T; K]),
    res: &mut [T; K],
    mut f: impl FnMut(i32, &[T; K]),
) {
    let mut it = first;
    while it.wrapping_sub(first).wrapping_abs() != int(K) && it != last {
        forward_recur_rotate_left(res);

        f(it, res);
        it = it.wrapping_sub(1);
    }

    if last.wrapping_sub(first).wrapping_abs() > int(K) {
        while it != last {
            let mut coef = [T::from(0.0); K];
            r(it, &mut coef);

            let tmp = dot(&coef, res);
            forward_recur_shift_left(res);
            res[K - 1] = tmp;

            f(it, res);
            it = it.wrapping_sub(1);
        }
    }
}

// legendre.h
//
// The `_t` suffix of the templates with the same name as a public function avoids a name clash.

fn legendre_p_initializer_n<T: Num>(z: T, res: &mut [T; 2]) {
    res[0] = T::from(1.0);
    res[1] = z;
}

fn legendre_p_recurrence_n<T: Num>(z: T, n: i32, res: &mut [T; 2]) {
    let n = f64::from(n);
    let fac0 = -T::from(n - 1.0) / T::from(n);
    let fac1 = T::from(2.0 * n - 1.0) / T::from(n);

    res[0] = fac0;
    res[1] = fac1 * z;
}

/// Compute the Legendre polynomial of degree n.
///
/// - `n`: degree of the polynomial
/// - `z`: argument of the polynomial, either real or complex
/// - `f`: a function to be called as `callback(j, res)` for `0 <= j <= n`
/// - `res`: value and derivatives of the polynomial
fn legendre_p_for_each_n<T: Num>(n: i32, z: T, res: &mut [T; 2], f: impl FnMut(i32, &[T; 2])) {
    legendre_p_initializer_n(z, res);

    let re_n = |n, res: &mut _| legendre_p_recurrence_n(z, n, res);
    forward_recur(0, n.wrapping_add(1), re_n, res, f);
}

/// Compute the Legendre polynomial of degree n.
///
/// - `n`: degree of the polynomial
/// - `z`: argument of the polynomial, either real or complex
///
/// Returns the value of the polynomial
fn legendre_p_t<T: Num>(n: i32, z: T) -> T {
    let mut res_n = [T::from(0.0); 2];
    legendre_p_for_each_n(n, z, &mut res_n, |_, _| {});

    res_n[1]
}

/// Compute all Legendre polynomials of degree j, where 0 <= j <= n.
///
/// - `z`: argument of the polynomials, either real or complex
/// - `res`: a view into a multidimensional array with element type T and size n + 1 to store the
///   value of each polynomial
fn legendre_p_all_t<T: Num>(z: T, mut res: Mdspan<'_, T, 1>) {
    let n = int(res.size()) - 1;

    let mut res_n = [T::from(0.0); 2];
    legendre_p_for_each_n(n, z, &mut res_n, |n, res_n| res[n] = res_n[1]);
}

// `NORM` replaces the `assoc_legendre_unnorm_policy` and `assoc_legendre_norm_policy` types

fn assoc_legendre_p_initializer_m_abs_m<const NORM: bool, T: Num>(
    m_signbit: bool,
    z: T,
    branch_cut: i32,
    res: &mut [T; 2],
) {
    let mut w;
    if branch_cut == 3 {
        w = (z - T::from(1.0)).sqrt() * (z + T::from(1.0)).sqrt(); // form for this branch cut
    } else {
        w = -(T::from(1.0) - z * z).sqrt(); // form for this branch cut
        if m_signbit {
            w = -w;
        }
    }

    if NORM {
        res[0] = T::from(1.0) / T::from(2.0).sqrt();
        res[1] = T::from(3.0).sqrt() * w / T::from(2.0);
    } else {
        res[0] = T::from(1.0);
        res[1] = w;

        if m_signbit {
            res[1] = res[1] / 2.0;
        }
    }
}

// other square roots can be avoided if each iteration increments by 2
fn assoc_legendre_p_recurrence_m_abs_m<const NORM: bool, T: Num>(
    z: T,
    branch_cut: i32,
    m: i32,
    res: &mut [T; 2],
) {
    let branch_cut_sign = if branch_cut == 3 {
        T::from(-1.0)
    } else {
        T::from(1.0)
    };
    let m_abs = f64::from(m.wrapping_abs());

    let fac = if NORM {
        branch_cut_sign
            * (T::from((2.0 * m_abs + 1.0) * (2.0 * m_abs - 1.0))
                / T::from(4.0 * m_abs * (m_abs - 1.0)))
            .sqrt()
    } else if m < 0 {
        branch_cut_sign / T::from((2.0 * m_abs) * (2.0 * m_abs - 2.0))
    } else {
        branch_cut_sign * T::from((2.0 * m_abs - 1.0) * (2.0 * m_abs - 3.0))
    };

    res[0] = fac * (T::from(1.0) - z * z);
    res[1] = T::from(0.0);
}

fn assoc_legendre_p_for_each_m_abs_m<const NORM: bool, T: Num>(
    m: i32,
    z: T,
    branch_cut: i32,
    res: &mut [T; 2],
    f: impl FnMut(i32, &[T; 2]),
) {
    let m_signbit = m < 0;

    assoc_legendre_p_initializer_m_abs_m::<NORM, T>(m_signbit, z, branch_cut, res);

    let re_m_abs_m =
        |m, res: &mut _| assoc_legendre_p_recurrence_m_abs_m::<NORM, T>(z, branch_cut, m, res);
    if m >= 0 {
        forward_recur(0, m.wrapping_add(1), re_m_abs_m, res, f);
    } else {
        backward_recur(0, m.wrapping_sub(1), re_m_abs_m, res, f);
    }
}

// Compute the associated Legendre polynomial of degree n and order n.
//
// We need to be careful with complex arithmetic, in particular the square roots
// should not be modified. This is because the sign bit of a real or imaginary part,
// even if it is equal to zero, can affect the branch cut.

fn assoc_legendre_p_initializer_n<const NORM: bool, T: Num>(
    m: i32,
    z: T,
    res_m_abs_m: T,
    res: &mut [T; 2],
) {
    let m_abs = f64::from(m.wrapping_abs());
    let fac = if NORM {
        T::from(2.0 * m_abs + 3.0).sqrt()
    } else {
        T::from(2.0 * (m_abs + 1.0) - 1.0) / T::from(m_abs + 1.0 - f64::from(m))
    };

    res[0] = res_m_abs_m;
    res[1] = fac * z * res_m_abs_m;
}

fn assoc_legendre_p_recurrence_n<const NORM: bool, T: Num>(m: i32, z: T, n: i32, res: &mut [T; 2]) {
    let (n, m) = (f64::from(n), f64::from(m));
    let (fac0, fac1);
    if NORM {
        let den = T::from((2.0 * n - 3.0) * (n * n - m * m));
        fac0 = -(T::from((2.0 * n + 1.0) * ((n - 1.0) * (n - 1.0) - m * m)) / den).sqrt();
        fac1 = (T::from((2.0 * n + 1.0) * (4.0 * (n - 1.0) * (n - 1.0) - 1.0)) / den).sqrt();
    } else {
        fac0 = -T::from(n + m - 1.0) / T::from(n - m);
        fac1 = T::from(2.0 * n - 1.0) / T::from(n - m);
    }

    res[0] = fac0;
    res[1] = fac1 * z;
}

fn assoc_legendre_p_pm1<const NORM: bool, T: Num>(n: i32, m: i32, z: T, res: &mut T) {
    if m == 0 {
        if z.real() >= 0.0 {
            *res = T::from(1.0);
        } else {
            *res = T::from(utils::neg_one_pow(n));
        }
        // Apply normalization for m=0
        if NORM {
            *res = *res * (T::from(2.0 * f64::from(n) + 1.0) / T::from(2.0)).sqrt();
        }
    } else {
        *res = T::from(0.0);
    }
}

/// Compute the associated Legendre polynomial of degree n and order m.
///
/// This is the overload that starts from the given `res_m_abs_m`.
///
/// - `n`: degree of the polynomial
/// - `m`: order of the polynomial
/// - `type`: specifies the branch cut of the polynomial, either 1, 2, or 3
/// - `z`: argument of the polynomial, either real or complex
/// - `callback`: a function to be called as `callback(j, m, type, z, p, p_prev, args...)` for
///   `0 <= j <= n`
/// - `args`: arguments to forward to the callback
///
/// Returns the value of the polynomial
fn assoc_legendre_p_for_each_n_from<const NORM: bool, T: Num>(
    n: i32,
    m: i32,
    z: T,
    _branch_cut: i32,
    res_m_abs_m: T,
    res: &mut [T; 2],
    mut f: impl FnMut(i32, &[T; 2]),
) {
    res[0] = T::from(0.0);
    res[1] = T::from(0.0);

    let m_abs = m.wrapping_abs();
    if m_abs > n {
        for j in 0..=n {
            f(j, res);
        }
    } else {
        for j in 0..m_abs {
            f(j, res);
        }

        if z.real().abs() == 1.0 && z.imag() == 0.0 {
            for j in m_abs..=n {
                forward_recur_shift_left(res);
                assoc_legendre_p_pm1::<NORM, T>(j, m, z, &mut res[1]);

                f(j, res);
            }
        } else {
            assoc_legendre_p_initializer_n::<NORM, T>(m, z, res_m_abs_m, res);

            let re_n = |n, res: &mut _| assoc_legendre_p_recurrence_n::<NORM, T>(m, z, n, res);
            forward_recur(m_abs, n.wrapping_add(1), re_n, res, f);
        }
    }
}

fn assoc_legendre_p_for_each_n<const NORM: bool, T: Num>(
    n: i32,
    m: i32,
    z: T,
    branch_cut: i32,
    res: &mut [T; 2],
    f: impl FnMut(i32, &[T; 2]),
) {
    assoc_legendre_p_for_each_m_abs_m::<NORM, T>(m, z, branch_cut, res, |_, _| {});

    let res_m_abs_m = res[1];
    assoc_legendre_p_for_each_n_from::<NORM, T>(n, m, z, branch_cut, res_m_abs_m, res, f);
}

fn assoc_legendre_p_for_each_n_m<const NORM: bool, T: Num>(
    n: i32,
    m: i32,
    z: T,
    branch_cut: i32,
    res: &mut [T; 2],
    mut f: impl FnMut(i32, i32, &[T; 2]),
) {
    let mut res_m_abs_m = [T::from(0.0); 2];
    assoc_legendre_p_for_each_m_abs_m::<NORM, T>(
        m,
        z,
        branch_cut,
        &mut res_m_abs_m,
        |m, res_m_abs_m| {
            res[0] = res_m_abs_m[1];

            assoc_legendre_p_for_each_n_from::<NORM, T>(
                n,
                m,
                z,
                branch_cut,
                res_m_abs_m[1],
                res,
                |n, res_n| {
                    f(n, m, res_n);
                },
            );
        },
    );
    assoc_legendre_p_for_each_m_abs_m::<NORM, T>(
        m.wrapping_neg(),
        z,
        branch_cut,
        &mut res_m_abs_m,
        |m, res_m_abs_m| {
            res[0] = res_m_abs_m[1];

            assoc_legendre_p_for_each_n_from::<NORM, T>(
                n,
                m,
                z,
                branch_cut,
                res_m_abs_m[1],
                res,
                |n, res_n| {
                    f(n, m, res_n);
                },
            );
        },
    );
}

/// Compute the associated Legendre polynomial of degree n and order m.
///
/// - `n`: degree of the polynomial
/// - `m`: order of the polynomial
/// - `type`: specifies the branch cut of the polynomial, either 1, 2, or 3
/// - `z`: argument of the polynomial, either real or complex
///
/// Returns the value of the polynomial
fn assoc_legendre_p_t<const NORM: bool, T: Num>(n: i32, m: i32, z: T, branch_cut: i32) -> T {
    let mut res_n = [T::from(0.0); 2];
    assoc_legendre_p_for_each_n::<NORM, T>(n, m, z, branch_cut, &mut res_n, |_, _| {});

    res_n[1]
}

/// Compute all associated Legendre polynomials of degree j and order i, where 0 <= j <= n and
/// -m <= i <= m.
///
/// - `type`: specifies the branch cut of the polynomial, either 1, 2, or 3
/// - `z`: argument of the polynomial, either real or complex
/// - `res`: a view into the output with element type T and extents (2 * m + 1, n + 1)
///
/// Returns the value of the polynomial
fn assoc_legendre_p_all_t<const NORM: bool, T: Num>(
    z: T,
    branch_cut: i32,
    mut res: Mdspan<'_, T, 2>,
) {
    let n = int(res.extent(0)) - 1;
    let m = (int(res.extent(1)) - 1) / 2;

    let mut p = [T::from(0.0); 2];
    let extent1 = int(res.extent(1));
    assoc_legendre_p_for_each_n_m::<NORM, T>(n, m, z, branch_cut, &mut p, |n, m, res_n_m| {
        if m >= 0 {
            res[(n, m)] = res_n_m[1];
        } else {
            res[(n, m + extent1)] = res_n_m[1];
        }
    });
}

fn sph_legendre_p_initializer_m_abs_m<T: Num>(m_signbit: bool, theta_sin: T, res: &mut [T; 2]) {
    let fac0 = T::from(1.0) / (T::from(2.0) * T::from(utils::SQRT_PI));
    let mut fac1 = -T::from(3.0).sqrt() / (T::from(2.0) * T::from(utils::SQRT_TAU_DOWN));
    if m_signbit {
        fac1 = -fac1;
    }

    res[0] = fac0;
    res[1] = fac1 * theta_sin.abs();
}

fn sph_legendre_p_recurrence_m_abs_m<T: Num>(theta_sin: T, m: i32, res: &mut [T; 2]) {
    let m_abs = f64::from(m.wrapping_abs());

    let fac = (T::from((2.0 * m_abs + 1.0) * (2.0 * m_abs - 1.0))
        / T::from(4.0 * m_abs * (m_abs - 1.0)))
    .sqrt();

    res[0] = fac * theta_sin * theta_sin;
    res[1] = T::from(0.0);
}

fn sph_legendre_p_for_each_m_abs_m<T: Num>(
    m: i32,
    theta: T,
    res: &mut [T; 2],
    f: impl FnMut(i32, &[T; 2]),
) {
    let m_signbit = m < 0;

    let theta_sin = theta.sin();

    sph_legendre_p_initializer_m_abs_m(m_signbit, theta_sin, res);

    let re_m_abs_m = |m, res: &mut _| sph_legendre_p_recurrence_m_abs_m(theta_sin, m, res);
    if m >= 0 {
        forward_recur(0, m.wrapping_add(1), re_m_abs_m, res, f);
    } else {
        backward_recur(0, m.wrapping_sub(1), re_m_abs_m, res, f);
    }
}

fn sph_legendre_p_initializer_n<T: Num>(m: i32, theta_cos: T, res_m_abs_m: T, res: &mut [T; 2]) {
    let fac = T::from(2.0 * f64::from(m.wrapping_abs()) + 3.0).sqrt();

    res[0] = res_m_abs_m;
    res[1] = fac * theta_cos * res_m_abs_m;
}

fn sph_legendre_p_recurrence_n<T: Num>(m: i32, theta_cos: T, n: i32, res: &mut [T; 2]) {
    let (n, m) = (f64::from(n), f64::from(m));
    let den = T::from((2.0 * n - 3.0) * (n * n - m * m));
    let fac0 = -(T::from((2.0 * n + 1.0) * ((n - 1.0) * (n - 1.0) - m * m)) / den).sqrt();
    let fac1 = (T::from((2.0 * n + 1.0) * (4.0 * (n - 1.0) * (n - 1.0) - 1.0)) / den).sqrt();

    res[0] = fac0;
    res[1] = fac1 * theta_cos;
}

/// Compute the spherical Legendre polynomial of degree n and order m.
///
/// This is the overload that starts from the given `res_m_abs_m`.
///
/// - `n`: degree of the polynomial
/// - `m`: order of the polynomial
/// - `theta`: z = cos(theta) argument of the polynomial, either real or complex
/// - `callback`: a function to be called as `callback(j, m, type, z, p, p_prev, args...)` for
///   `0 <= j <= n`
/// - `args`: arguments to forward to the callback
///
/// Returns the value of the polynomial
fn sph_legendre_p_for_each_n_from<T: Num>(
    n: i32,
    m: i32,
    theta: T,
    res_m_abs_m: T,
    res: &mut [T; 2],
    mut f: impl FnMut(i32, &[T; 2]),
) {
    res[0] = T::from(0.0);
    res[1] = T::from(0.0);

    let m_abs = m.wrapping_abs();
    if m_abs > n {
        for j in 0..=n {
            f(j, res);
        }
    } else {
        for j in 0..m_abs {
            f(j, res);
        }

        let theta_cos = theta.cos();

        sph_legendre_p_initializer_n(m, theta_cos, res_m_abs_m, res);

        let re_n = |n, res: &mut _| sph_legendre_p_recurrence_n(m, theta_cos, n, res);
        forward_recur(m_abs, n.wrapping_add(1), re_n, res, f);
    }
}

fn sph_legendre_p_for_each_n<T: Num>(
    n: i32,
    m: i32,
    theta: T,
    res: &mut [T; 2],
    f: impl FnMut(i32, &[T; 2]),
) {
    sph_legendre_p_for_each_m_abs_m(m, theta, res, |_, _| {});

    let res_m_abs_m = res[1];
    sph_legendre_p_for_each_n_from(n, m, theta, res_m_abs_m, res, f);
}

pub(crate) fn sph_legendre_p_for_each_n_m<T: Num>(
    n: i32,
    m: i32,
    theta: T,
    res: &mut [T; 2],
    mut f: impl FnMut(i32, i32, &[T; 2]),
) {
    let mut res_m_abs_m = [T::from(0.0); 2];
    sph_legendre_p_for_each_m_abs_m(m, theta, &mut res_m_abs_m, |m, res_m_abs_m| {
        res[0] = res_m_abs_m[1];

        sph_legendre_p_for_each_n_from(n, m, theta, res_m_abs_m[1], res, |n, res_n| {
            f(n, m, res_n);
        });
    });
    sph_legendre_p_for_each_m_abs_m(m.wrapping_neg(), theta, &mut res_m_abs_m, |m, res_m_abs_m| {
        res[0] = res_m_abs_m[1];

        sph_legendre_p_for_each_n_from(n, m, theta, res_m_abs_m[1], res, |n, res_n| {
            f(n, m, res_n);
        });
    });
}

pub(crate) fn sph_legendre_p_t<T: Num>(n: i32, m: i32, theta: T) -> T {
    let mut res_n = [T::from(0.0); 2];
    sph_legendre_p_for_each_n(n, m, theta, &mut res_n, |_, _| {});

    res_n[1]
}

fn sph_legendre_p_all_t<T: Num>(theta: T, mut res: Mdspan<'_, T, 2>) {
    let n_max = int(res.extent(0)) - 1;
    let m_max = (int(res.extent(1)) - 1) / 2;

    let mut res_n_m = [T::from(0.0); 2];
    sph_legendre_p_for_each_n_m(n_max, m_max, theta, &mut res_n_m, |n, m, res_n_m| {
        if m >= 0 {
            res[(n, m)] = res_n_m[1];
        } else {
            res[(n, m + 2 * m_max + 1)] = res_n_m[1];
        }
    });
}

/// `Q_0(x) = ln |(1 + x) / (1 - x)| / 2` without cancellation, for real `x != +-1`
fn lq0(x: f64) -> f64 {
    let ax = x.abs();
    let t = if ax < 1.0 {
        2.0 * ax / (1.0 - ax)
    } else {
        2.0 / (ax - 1.0)
    };
    (t.ln_1p() / 2.0).copysign(x)
}

/// `Q_0(z) = ln(ls (1 + z) / (1 - z)) / 2` for complex `z`, with `ls = -1` for `|z| > 1`, else 1
fn q0_c64(z: C64) -> C64 {
    if z.is_nan() {
        return C64::NAN;
    }
    if z.re.is_infinite() || z.im.is_infinite() {
        return C64::ZERO;
    }
    // Unlike in xsf, subtract the logs when the quotient can overflow or underflow.
    // Here |z| rounds to 1 and the argument difference is on the principal branch.
    if z.re.abs() == 1.0 && z.im.is_subnormal() {
        return ((1.0 + z).ln() - (1.0 - z).ln()) * 0.5;
    }
    // the quotient loses digits for `|z| < 1/2` and `|z| > 2`, where `Q_0(z) = atanh(w)`, with
    // `w = z` or `1 / z`, so there `Re atanh(w) = ln_1p(t) / 4` with `t = 4 Re w / |1 - w|^2`, and
    // `Im atanh(w) = atan2(2 Im w, 1 - |w|^2) / 2`
    let a = z.norm();
    let w = if a > 2.0 { 1.0 / z } else { z };
    let t = 4.0 * w.re / (1.0 - w).norm_sqr();
    if (0.5..=2.0).contains(&a) {
        let ls = if a > 1.0 { -1.0 } else { 1.0 };
        let mut q = ((1.0 + z) * ls / (1.0 - z)).ln() / 2.0;
        // in between, only a tiny `Re w` loses digits
        if (-0.5..=1.0).contains(&t) {
            q.re = t.ln_1p() / 4.0;
        }
        return q;
    }
    c64(t.ln_1p() / 4.0, (2.0 * w.im).atan2(1.0 - w.norm_sqr()) / 2.0)
}

/// Replaces the imaginary parts of `Q` and `Q'` (by `lq`, `lqn_c64` or `lqmn_c64`, up to degree
/// plus order `n`) at `z` with a tiny `y = Im z` by those at `z` with `y` lifted by a power of 2
/// to about `2^-36 d`, scaled back, with `d` the distance of `Re z` to `+-1` over `n + 1`, or
/// `sqrt|Re z|` if smaller and `Re z != 0`
///
/// The imaginary parts of the seeds are about `y Q'(Re z)`, which can be subnormal while those
/// of the results are not. Up to the lifted `y`, the imaginary parts are linear in `y` to about
/// `2^-70 y` times the size of `Q'` near `Re z`, and the real parts do not depend on `y` to double
/// precision, so they are kept. So are the imaginary parts that overflow in the lifted run,
/// unless they are NaN.
fn lq_lift_im<const R: usize>(
    lq: fn(C64, Mdspan<'_, C64, R>, Mdspan<'_, C64, R>),
    z: C64,
    n: i32,
    q: &mut Mdspan<'_, C64, R>,
    dq: &mut Mdspan<'_, C64, R>,
) {
    let x = z.re.abs();
    let mut d = (x - 1.0).abs() / f64::from(n + 1);
    if x != 0.0 {
        // the lifted `y^2` stays far below `|x|`, since the parts linear in `x` cancel to it
        d = d.min(x.sqrt());
    }
    let y = z.im.abs();
    // also not for `d = 0` (where the expansion does not hold), and NaN
    if !(0.0 < y && y < d * const { utils::pow2(-37) } && d.is_finite()) {
        return;
    }
    // at least 1, and the lifted `y` is at least `2^-37 d`, so that it is not lifted again
    let s = frexp_f64(d).1 - 36 - frexp_f64(y).1;
    let len = q.data.len();
    let mut lifted = vec![C64::ZERO; 2 * len];
    let (lq_, ldq) = lifted.split_at_mut(len);
    lq(c64(z.re, ldexp_f64(z.im, s)), Mdspan::new(lq_, q.extents), Mdspan::new(ldq, dq.extents));
    for (c, l) in q.data.iter_mut().chain(dq.data.iter_mut()).zip(lifted) {
        if l.im.is_finite() || c.im.is_nan() {
            c.im = ldexp_f64(l.im, -s);
        }
    }
}

// ====================================================
// Purpose: Compute Legendre functions Qn(x) & Qn'(x)
// Input :  x  --- Argument of Qn(x)
//          n  --- Degree of Qn(x)  ( n = 0,1,2,…)
// Output:  QN(n) --- Qn(x)
//          QD(n) --- Qn'(x)
// ====================================================

fn lqn(x: f64, mut qn: Mdspan<'_, f64, 1>, mut qd: Mdspan<'_, f64, 1>) {
    const EPS: f64 = 1.0e-16;
    let n = int(qn.size()) - 1;

    let (mut q0, mut q1, mut qf, mut qr, mut qf0, mut qf1, mut qf2);

    if x.abs() == 1.0 {
        for k in 0..=n {
            qn[k] = SPECFUN_INF;
            qd[k] = SPECFUN_INF;
        }
        return;
    }

    if x.abs() <= 1.021 {
        q0 = lq0(x);
        q1 = x * q0 - 1.0;
        qn[0] = q0;
        qn[1] = q1;
        qd[0] = 1.0 / (1.0 - x * x);
        qd[1] = qn[0] + x * qd[0];

        for k in 2..=n {
            let k_ = f64::from(k);
            qf = ((2.0 * k_ - 1.0) * x * q1 - (k_ - 1.0) * q0) / k_;
            qn[k] = qf;
            qd[k] = (qn[k - 1] - x * qf) * k_ / (1.0 - x * x);
            q0 = q1;
            q1 = qf;
        }
    } else {
        let mut qc1 = if n == 0 { 0.0 } else { 1.0 / x };
        let mut qc2 = 1.0 / x;

        for j in 1..=n {
            let j_ = f64::from(j);
            qc2 *= j_ / ((2.0 * j_ + 1.0) * x);
            if j == n - 1 {
                qc1 = qc2;
            }
        }

        for (degree, prefactor) in [(n - 1, qc1), (n, qc2)] {
            // The series can overflow when its prefactor has already underflowed.
            if prefactor == 0.0 {
                qn[degree] = prefactor;
                continue;
            }
            let nl = f64::from(degree) + 1.0;
            qf = 1.0;
            qr = 1.0;

            // Just above |x| = 1.021 the terms decay like 1/x^2, taking up to about
            // 3800 terms for n near 1000.
            for k in 1..=5000 {
                let k = f64::from(k);
                qr = qr * (0.5 * nl + k - 1.0) * (0.5 * (nl - 1.0) + k)
                    / ((nl + k - 0.5) * k * x * x);
                qf += qr;
                if (qr / qf).abs() < EPS {
                    break;
                }
            }

            qn[degree] = qf * prefactor;
        }

        qf2 = qn[n];
        qf1 = qn[n - 1];

        for k in (2..=n).rev() {
            let k_ = f64::from(k);
            qf0 = ((2.0 * k_ - 1.0) * x * qf1 - k_ * qf2) / (k_ - 1.0);
            qn[k - 2] = qf0;
            qf2 = qf1;
            qf1 = qf0;
        }

        qd[0] = 1.0 / (1.0 - x * x);

        for k in 1..=n {
            qd[k] = f64::from(k) * (qn[k - 1] - x * qn[k]) / (1.0 - x * x);
        }
    }
}

// ==================================================
// Purpose: Compute the Legendre functions Qn(z) and
//          their derivatives Qn'(z) for a complex
//          argument
// Input :  x --- Real part of z
//          y --- Imaginary part of z
//          n --- Degree of Qn(z), n = 0,1,2,...
// Output:  CQN(n) --- Qn(z)
//          CQD(n) --- Qn'(z)
// ==================================================

fn lqn_c64(z: C64, mut cqn: Mdspan<'_, C64, 1>, mut cqd: Mdspan<'_, C64, 1>) {
    let n = int(cqn.size()) - 1;

    let (cq0, cq1);
    let mut cqf0 = C64::ZERO;
    let (mut cqf1, mut cqf2);

    // Unlike in xsf, require zero imaginary part for the pole at 1.
    if z.re == 1.0 && z.is_real() {
        for k in 0..=n {
            cqn[k] = c64(SPECFUN_INF, 0.0);
            cqd[k] = c64(SPECFUN_INF, 0.0);
        }
        return;
    }

    cq0 = q0_c64(z);
    cq1 = z * cq0 - 1.0;

    cqn[0] = cq0;
    cqn[1] = cq1;

    if z.norm() < 1.0001 {
        cqf0 = cq0;
        cqf1 = cq1;
        for k in 2..=n {
            let k_ = f64::from(k);
            cqf2 = (z * (2.0 * k_ - 1.0) * cqf1 - cqf0 * (k_ - 1.0)) / k_;
            cqn[k] = cqf2;
            cqf0 = cqf1;
            cqf1 = cqf2;
        }
    } else {
        let km = if z.norm() > 1.1 {
            40_i32.wrapping_add(n)
        } else {
            // Unlike in xsf, this uses `|z| - 1` instead of `|z - 1|`, which is not positive for
            // `z` far from 1 and left `cqf0` at zero (as in `lqmn_large`)
            (f64::from(40_i32.wrapping_add(n)) * (-1.0 - 1.8 * (z.norm() - 1.0).ln()).floor())
                as i32
        };

        cqf2 = C64::ZERO;
        cqf1 = C64::ONE;
        for k in (0..=km).rev() {
            let k_ = f64::from(k);
            cqf0 = (z * (2.0 * k_ + 3.0) * cqf1 - cqf2 * (k_ + 2.0)) / (k_ + 1.0);
            if k <= n {
                cqn[k] = cqf0;
            }
            // Unlike in xsf, this rescales by a power of 2 (as xsf's `lqmn` does), since `cqf0`
            // grows like `(2 z)^k` and overflowed for large `|z|`
            if cqf0.l1_norm() > SPECFUN_BIG {
                (cqf0, cqf1) = (cqf0 / SPECFUN_BIG, cqf1 / SPECFUN_BIG);
                for j in k.min(n)..=n {
                    cqn[j] /= SPECFUN_BIG;
                }
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }
        for k in 1..=n {
            cqn[k] *= cq0 / cqf0;
        }
        cqn[0] = cq0;
    }
    cqd[0] = (cqn[1] - z * cqn[0]) / (z * z - 1.0);

    for k in 1..=n {
        let k_ = f64::from(k);
        cqd[k] = (z * k_ * cqn[k] - cqn[k - 1] * k_) / (z * z - 1.0);
    }

    // Unlike in xsf, the imaginary parts for a tiny `Im z` are those at a lifted one, since the
    // imaginary parts of the seeds of the recurrence can be subnormal
    lq_lift_im(lqn_c64, z, n, &mut cqn, &mut cqd);

    // For a tiny real coordinate, use the first-order terms at z = i*y. The dominant parts
    // are even in Re z, so their values from the recurrence already include them to rounding
    // accuracy. The small part of Q is x times Q', and Legendre's equation gives x times Q''.
    let (x, y) = (z.re, z.im);
    let d = (y.abs() / f64::from(n + 1)).min(1.0);
    let denom = 1.0 + y * y;
    if x != 0.0 && x.abs() < d * const { utils::pow2(-37) } && y.abs() > 1.1 && denom.is_finite() {
        for k in 0..=n {
            let degree = f64::from(k);
            let (q, dq) = (&mut cqn[k], &mut cqd[k]);
            let (small_q, small_dq, large_q, large_dq, cross) = if k % 2 == 0 {
                (&mut q.re, &mut dq.im, q.im, dq.re, y)
            } else {
                (&mut q.im, &mut dq.re, q.re, dq.im, -y)
            };
            let q_linear = x * large_dq;
            let numerator = 2.0 * cross * large_dq - degree * (degree + 1.0) * large_q;
            let dq_linear = x * (numerator / denom);
            if q_linear.is_finite() && (q_linear != 0.0 || *small_q != 0.0) {
                *small_q = q_linear;
            }
            if dq_linear.is_finite() && (dq_linear != 0.0 || *small_dq != 0.0) {
                *small_dq = dq_linear;
            }
        }
    }
}

/// The `Q_n^m(z)` of `lqmn` and `lqmn_c64` for `|z| >= 1.0001`, from `Q_0^0(z) = q0` and
/// `zq = sqrt(z^2 - 1)`, with the backward recurrences in the degree of the orders 0 and 1 from a
/// start degree `km` and the recurrence in the order; xsf writes these out separately for real and
/// complex `z`, and the backward recurrences also for each of the orders 0 and 1
fn lqmn_large<T: Num>(z: T, zq: T, q0: T, mut qm: Mdspan<'_, T, 2>) {
    let m = int(qm.extent(0)) - 1;
    let n = int(qm.extent(1)) - 1;
    // Unlike in xsf (scipy/xsf#368), `log(|x| - 1)` for real `x` (as for complex `z`) instead of
    // `log(x - 1)`
    let km = if z.abs() > 1.1 {
        40_i32.wrapping_add(m).wrapping_add(n)
    } else {
        40_i32
            .wrapping_add(m)
            .wrapping_add(n)
            .wrapping_mul((-1.0 - 1.8 * (z.abs() - 1.0).ln()) as i32)
    };
    let l1 = |q: T| q.real().abs() + q.imag().abs();
    let ldexp = |q: T, e: i32| if e == 0 { q } else { q.ldexp(e) };

    let mut cs = [T::from(0.0); 2];
    let mut ex = vec![[0; 2]; n as usize + 1];
    for (i, qi0) in [(0, q0), (1, T::from(-1.0) / zq)] {
        let i_ = f64::from(i);
        let (mut qf0, mut qf1, mut qf2) = (T::from(0.0), T::from(1.0), T::from(0.0));
        let mut sc = 0;
        for k in (0..=km).rev() {
            let k_ = f64::from(k);
            qf0 = (z * (2.0 * k_ + 3.0) * qf1 - qf2 * (k_ + 2.0 - i_)) / (k_ + 1.0 + i_);
            if k <= n {
                qm[(i, k)] = qf0;
                ex[k as usize][i as usize] = sc;
            }
            // rescaled by powers of 2, since the recurrence grows like `(2 z)^k`, which overflows
            // for large `|z|` and `m + n`
            //
            // Unlike in xsf (scipy/xsf#571), `ex` keeps the scale of each stored value, instead of
            // rescaling the stored values, which can become subnormal
            if l1(qf0) > SPECFUN_BIG {
                (qf0, qf1) = (qf0 / SPECFUN_BIG, qf1 / SPECFUN_BIG);
                sc += 1;
            }
            qf2 = qf1;
            qf1 = qf0;
        }

        for e in &mut ex {
            e[i as usize] = SPECFUN_BIG_LOG2 * (e[i as usize] - sc);
        }
        cs[i as usize] = qi0 / qf0;
    }

    // `s c 2^k`, with `s` scaled to [0.5, 1) first, so that neither factor is subnormal
    let seed = |s: T, c: T, k: i32| {
        if k == 0 {
            return s * c;
        }
        let a = -frexp_f64(l1(s)).1;
        ldexp(s, a) * ldexp(c, k - a)
    };
    for j in 0..=n {
        // Unlike in xsf, the normalization of the orders 0 and 1 is deferred to here, where it is
        // lifted by `2^e` if they would be subnormal, since the recurrence in the order can raise
        // them back into range
        let o = ex[j as usize];
        let seeds = |e: i32| {
            [0, 1].map(|i| match (i, j) {
                // as is for complex `z`, since the normalization can cost the digits of a part
                (0, 0) if T::COMPLEX => ldexp(q0, e),
                _ => seed(qm[(i, j)], cs[i as usize], e + o[i as usize]),
            })
        };
        let [mut q0, mut q1] = seeds(0);
        let mut e = 0;
        if l1(q0).max(l1(q1)) < f64::MIN_POSITIVE {
            // to about `MIN_POSITIVE / EPSILON` for the larger one, where all its digits are normal
            let k_min = frexp_f64(f64::MIN_POSITIVE / f64::EPSILON).1;
            let k = |i: i32| frexp_f64(l1(qm[(i, j)])).1 + frexp_f64(l1(cs[i as usize])).1;
            e = (0..2)
                .filter(|&i| l1(qm[(i, j)]) != 0.0)
                .map(|i| k(i) + o[i as usize])
                .max()
                .map_or(0, |k| (k_min - k).max(0));
            [q0, q1] = seeds(e);
        }
        qm[(0, j)] = ldexp(q0, -e);
        qm[(1, j)] = ldexp(q1, -e);
        let lifted = e != 0;
        for i in 0..=(m - 2) {
            let (i_, j_) = (f64::from(i), f64::from(j));
            let mut qf = z * -(2.0 * (i_ + 1.0)) / zq * q1 + q0 * ((j_ - i_) * (j_ + i_ + 1.0));
            // the lifted values can overflow before they are scaled back, and their products
            // before they cancel, also once the lift is undone
            if lifted && l1(qf) > SPECFUN_BIG {
                (qf, q1) = (qf / SPECFUN_BIG, q1 / SPECFUN_BIG);
                e -= SPECFUN_BIG_LOG2;
            }
            qm[(i + 2, j)] = ldexp(qf, -e);
            q0 = q1;
            q1 = qf;
        }
    }
}

// Use the leading term in 1/x. Keep its exponent separate because raising
// the order can bring an underflowed value back into range.
fn lqmn_extreme_real(x: f64, mut qm: Mdspan<'_, f64, 2>, mut qd: Mdspan<'_, f64, 2>) {
    #[expect(
        clippy::cast_possible_wrap,
        reason = "the binary64 mantissa width fits i32"
    )]
    const MIN_SUB_EXP: i32 = f64::MIN_EXP - f64::MANTISSA_DIGITS as i32;

    let (xm, xe) = frexp_f64(x);
    let xe = i64::from(xe);
    let inv = DoubleDouble::ONE / xm;
    let (mut seed, e) = inv.frexp();
    let mut seed_e = i64::from(e) - xe;
    let order_sign = -x.signum();
    // A value mantissa is normalized; a derivative multiplies two such mantissas.
    // Their magnitudes lie between 1/4 and 1.
    let scaled = |a: DoubleDouble, e: i64| {
        if e > i64::from(f64::MIN_EXP) {
            return ldexp_f64(a.hi, e.min(i64::from(i32::MAX)) as i32);
        }
        if e < i64::from(MIN_SUB_EXP - 1) {
            return 0.0_f64.copysign(a.hi);
        }
        // Round in units of the least subnormal. At an exact midpoint
        // the first omitted term increases the magnitude.
        let units = ldexp_f64(a.hi, (e - i64::from(MIN_SUB_EXP)) as i32);
        let rounded = if units.fract().abs() == 0.5
            && a.lo != 0.0
            && a.lo.is_sign_negative() != a.hi.is_sign_negative()
        {
            units.trunc()
        } else {
            units.round()
        };
        ldexp_f64(rounded, MIN_SUB_EXP)
    };
    for n in 0..int(qm.extent(1)) {
        let next = f64::from(n) + 1.0;
        let (slope, e) = (inv * -next).frexp();
        let slope_e = i64::from(e) - xe;
        let (mut value, mut exp) = (seed, seed_e);
        for m in 0..int(qm.extent(0)) {
            qm[(m, n)] = scaled(value, exp);
            qd[(m, n)] = scaled(value * slope, exp + slope_e);
            let (v, e) = (value * (order_sign * (next + f64::from(m)))).frexp();
            (value, exp) = (v, exp + i64::from(e));
        }
        let (v, e) = (seed * -slope / (2.0 * next + 1.0)).frexp();
        (seed, seed_e) = (v, seed_e + slope_e + i64::from(e));
    }
}

// ==========================================================
// Purpose: Compute the associated Legendre functions of the
//          second kind, Qmn(x) and Qmn'(x)
// Input :  x  --- Argument of Qmn(x)
//          m  --- Order of Qmn(x)  ( m = 0,1,2,… )
//          n  --- Degree of Qmn(x) ( n = 0,1,2,… )
//          mm --- Physical dimension of QM and QD
// Output:  QM(m,n) --- Qmn(x)
//          QD(m,n) --- Qmn'(x)
// ==========================================================

fn lqmn(x: f64, mut qm: Mdspan<'_, f64, 2>, mut qd: Mdspan<'_, f64, 2>) {
    // Leave room for 2*k+3 and the scale in the backward recurrence.
    #[expect(clippy::cast_possible_wrap, reason = "the i32 bit width fits i32")]
    const LARGE_X: f64 = utils::pow2(f64::MAX_EXP - (i32::BITS as i32 + 1) - SPECFUN_BIG_LOG2);

    let m = int(qm.extent(0)) - 1;
    let n = int(qm.extent(1)) - 1;

    if x.abs() == 1.0 {
        for i in 0..=m {
            for j in 0..=n {
                qm[(i, j)] = SPECFUN_INF;
                qd[(i, j)] = SPECFUN_INF;
            }
        }
        return;
    }
    // Unlike in xsf, use the leading terms before the backward recurrence can overflow.
    if x.abs() >= LARGE_X && x.is_finite() {
        lqmn_extreme_real(x, qm, qd);
        return;
    }

    let mut ls = 1.0;
    if x.abs() > 1.0 {
        ls = -1.0;
    }
    let xs = ls * (1.0 - x * x);
    let xq = xs.sqrt();
    let q0 = lq0(x);
    if x.abs() < 1.0001 {
        qm[(0, 0)] = q0;
        qm[(0, 1)] = x * q0 - 1.0;
        qm[(1, 0)] = -1.0 / xq;
        qm[(1, 1)] = -ls * xq * (q0 + x / (1. - x * x));
        for i in 0..=1 {
            for j in 2..=n {
                let (i_, j_) = (f64::from(i), f64::from(j));
                qm[(i, j)] = ((2.0 * j_ - 1.) * x * qm[(i, j - 1)]
                    - (j_ + i_ - 1.0) * qm[(i, j - 2)])
                    / (j_ - i_);
            }
        }
        /* 15 */
        for i in 2..=m {
            for j in 0..=n {
                let (i_, j_) = (f64::from(i), f64::from(j));
                qm[(i, j)] = -2.0 * (i_ - 1.0) * x / xq * qm[(i - 1, j)]
                    - ls * (j_ + i_ - 1.0) * (j_ - i_ + 2.0) * qm[(i - 2, j)];
            }
        }
    } else {
        lqmn_large(x, xq, q0, Mdspan::new(qm.data, qm.extents));
    }

    let exterior = ls < 0.0 && x.is_finite();
    let derivative = |value: f64, i, j| {
        // Unlike in xsf, restore the derivative's exterior parity after underflow.
        if value == 0.0 && exterior {
            utils::neg_if_odd(if x < 0.0 { j } else { i }, -0.0)
        } else {
            value
        }
    };
    qd[(0, 0)] = ls / xs;
    for j in 1..=n {
        let value = ls * f64::from(j) * (qm[(0, j - 1)] - x * qm[(0, j)]) / xs;
        qd[(0, j)] = derivative(value, 0, j);
    }

    for i in 1..=m {
        for j in 0..=n {
            let (i_, j_) = (f64::from(i), f64::from(j));
            let value =
                ls * i_ * x / xs * qm[(i, j)] + (i_ + j_) * (j_ - i_ + 1.) / xq * qm[(i - 1, j)];
            qd[(i, j)] = derivative(value, i, j);
        }
    }
}

// =======================================================
// Purpose: Compute the associated Legendre functions of
//          the second kind, Qmn(z) and Qmn'(z), for a
//          complex argument
// Input :  x  --- Real part of z
//          y  --- Imaginary part of z
//          m  --- Order of Qmn(z)  ( m = 0,1,2,… )
//          n  --- Degree of Qmn(z) ( n = 0,1,2,… )
//          mm --- Physical dimension of CQM and CQD
// Output:  CQM(m,n) --- Qmn(z)
//          CQD(m,n) --- Qmn'(z)
// =======================================================

fn lqmn_c64(z: C64, mut cqm: Mdspan<'_, C64, 2>, mut cqd: Mdspan<'_, C64, 2>) {
    let m = int(cqm.extent(0)) - 1;
    let n = int(cqm.extent(1)) - 1;

    if z.re.abs() == 1.0 && z.is_real() {
        for i in 0..=m {
            for j in 0..=n {
                cqm[(i, j)] = c64(SPECFUN_INF, 0.0);
                cqd[(i, j)] = c64(SPECFUN_INF, 0.0);
            }
        }

        return;
    }

    let xc = z.norm();
    let mut ls = 0;
    // Unlike in xsf, use the same branch as q0_c64 when the norm equals one.
    if z.is_real() || xc <= 1.0 {
        ls = 1;
    }
    if xc > 1.0 {
        ls = -1;
    }
    let zs = (1.0 - z * z) * f64::from(ls);
    let zq = zs.sqrt();

    let cq0 = q0_c64(z);
    if xc < 1.0001 {
        cqm[(0, 0)] = cq0;
        cqm[(1, 0)] = -1.0 / zq;
        cqm[(0, 1)] = z * cq0 - 1.0;
        // Unlike in xsf, use the exterior square-root factor for Q_1^1.
        let root = if ls < 0 { zq } else { -zq };
        cqm[(1, 1)] = root * (cq0 + z / (1.0 - z * z));

        for i in 0..=1 {
            for j in 2..=n {
                let (i_, j_) = (f64::from(i), f64::from(j));
                cqm[(i, j)] = (z * (2.0 * j_ - 1.0) * cqm[(i, j - 1)]
                    - cqm[(i, j - 2)] * (j_ + i_ - 1.0))
                    / (j_ - i_);
            }
        }

        for i in 2..=m {
            for j in 0..=n {
                let i_ = f64::from(i);
                cqm[(i, j)] = z * (-2.0 * (i_ - 1.0)) / zq * cqm[(i - 1, j)]
                    - cqm[(i - 2, j)] * int_product(&[ls, j + i - 1, j - i + 2]);
            }
        }
    } else {
        lqmn_large(z, zq, cq0, Mdspan::new(cqm.data, cqm.extents));
    }

    cqd[(0, 0)] = f64::from(ls) / zs;
    for j in 1..=n {
        cqd[(0, j)] = (cqm[(0, j - 1)] - z * cqm[(0, j)]) * (f64::from(ls) * f64::from(j)) / zs;
    }

    for i in 1..=m {
        for j in 0..=n {
            cqd[(i, j)] = z * f64::from(ls * i) / zs * cqm[(i, j)]
                + int_product(&[i + j, j - i + 1]) / zq * cqm[(i - 1, j)];
        }
    }

    // Unlike in xsf, the imaginary parts for a tiny `Im z` are those at a lifted one, since those
    // of the seeds of the order recurrence can be subnormal, while those of the orders it raises
    // back into range are not
    lq_lift_im(lqmn_c64, z, m + n, &mut cqm, &mut cqd);
}

// Rust API

pub trait LegendrePArg: crate::sealed::Sealed + Sized {
    fn legendre_p(self, n: c_int) -> Self;
    fn legendre_p_all(self, n: usize) -> Vec<Self>;
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self;
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>>;
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self;
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>>;
}

fn legendre_p_all_vec<T: Num + Into<P>, P>(z: T, n: usize) -> Vec<P> {
    let mut pn = vec![T::from(f64::NAN); n + 1];
    legendre_p_all_t(z, Mdspan::new(&mut pn, [n + 1]));
    pn.into_iter().map(Into::into).collect()
}

fn sph_legendre_p_all_vec<T: Num + Into<P>, P>(z: T, n: usize, m: usize) -> Vec<Vec<P>> {
    let extents = [n + 1, 2 * m + 1];
    let mut pnm = vec![T::from(f64::NAN); extents[0] * extents[1]];
    sph_legendre_p_all_t(z, Mdspan::new(&mut pnm, extents));
    utils::vec_to_vecvec(&pnm, extents[1])
}

fn assoc_legendre_p_dyn<T: Num>(z: T, n: c_int, m: c_int, bc: c_int, norm: bool) -> T {
    if norm {
        assoc_legendre_p_t::<true, T>(n, m, z, bc)
    } else {
        assoc_legendre_p_t::<false, T>(n, m, z, bc)
    }
}

fn assoc_legendre_p_all_vec<T: Num + Into<P>, P>(
    z: T,
    n: usize,
    m: usize,
    bc: c_int,
    norm: bool,
) -> Vec<Vec<P>> {
    let extents = [n + 1, 2 * m + 1];
    let mut pnm = vec![T::from(f64::NAN); extents[0] * extents[1]];
    let res = Mdspan::new(&mut pnm, extents);
    if norm {
        assoc_legendre_p_all_t::<true, T>(z, bc, res);
    } else {
        assoc_legendre_p_all_t::<false, T>(z, bc, res);
    }
    utils::vec_to_vecvec(&pnm, extents[1])
}

/// `$T` is the public type, and `$K` the type of the computation
macro_rules! impl_legendre_p_arg {
    ($T:ty, $K:ty) => {
        impl LegendrePArg for $T {
            #[inline]
            fn legendre_p(self, n: c_int) -> Self {
                legendre_p_t(n, <$K>::from(self)).into()
            }

            #[inline]
            fn legendre_p_all(self, n: usize) -> Vec<Self> {
                legendre_p_all_vec(<$K>::from(self), n)
            }

            #[inline]
            fn sph_legendre_p(self, n: c_int, m: c_int) -> Self {
                sph_legendre_p_t(n, m, <$K>::from(self)).into()
            }

            #[inline]
            fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>> {
                sph_legendre_p_all_vec(<$K>::from(self), n, m)
            }

            #[inline]
            fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self {
                assoc_legendre_p_dyn(<$K>::from(self), n, m, bc, norm).into()
            }

            #[inline]
            fn assoc_legendre_p_all(
                self,
                n: usize,
                m: usize,
                bc: c_int,
                norm: bool,
            ) -> Vec<Vec<Self>> {
                assoc_legendre_p_all_vec(<$K>::from(self), n, m, bc, norm)
            }
        }
    };
}

impl_legendre_p_arg!(f64, f64);
impl_legendre_p_arg!(Complex64, C64);

pub trait LegendreQArg: crate::sealed::Sealed + Sized {
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>);
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>);
}

/// `$T` is the public type, and `$K` the type of the computation
macro_rules! impl_legendre_q_arg {
    ($T:ty, $K:ty, $lqn:ident, $lqmn:ident) => {
        impl LegendreQArg for $T {
            #[inline]
            fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
                // like in scipy, since `lqn` writes at least 2 terms
                let n1 = n.max(1);

                let mut qn = vec![<$K>::from(f64::NAN); n1 + 1];
                let mut qd = vec![<$K>::from(f64::NAN); n1 + 1];

                $lqn(
                    <$K>::from(self),
                    Mdspan::new(&mut qn, [n1 + 1]),
                    Mdspan::new(&mut qd, [n1 + 1]),
                );
                let head = |q: Vec<$K>| q.into_iter().take(n + 1).map(Into::into).collect();
                (head(qn), head(qd))
            }

            #[inline]
            fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>) {
                // like in scipy, since `lqmn` writes at least 2 x 2 terms
                let (mm, nn) = (m.max(1), n.max(1));

                let extents = [mm + 1, nn + 1];
                let mut qm = vec![<$K>::from(f64::NAN); extents[0] * extents[1]];
                let mut qd = vec![<$K>::from(f64::NAN); extents[0] * extents[1]];

                $lqmn(
                    <$K>::from(self),
                    Mdspan::new(&mut qm, extents),
                    Mdspan::new(&mut qd, extents),
                );

                // [order][degree] -> [degree][order], without the padding of `mm` and `nn`
                let by_degree = |q: &[$K]| -> Vec<Vec<$T>> {
                    (0..=n)
                        .map(|j| (0..=m).map(|i| q[i * (nn + 1) + j].into()).collect())
                        .collect()
                };
                (by_degree(&qm), by_degree(&qd))
            }
        }
    };
}

impl_legendre_q_arg!(f64, f64, lqn, lqmn);
impl_legendre_q_arg!(Complex64, C64, lqn_c64, lqmn_c64);

/// Legendre polynomial of degree `n`
///
/// Corresponds to [`scipy.special.legendre_p`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre_p.html
///
/// # See also
/// - [`legendre_p_all`]: All Legendre polynomials up to degree `n`
/// - [`eval_legendre`](crate::eval_legendre): $P_n(z)$, also for non-integer $n$
#[doc(alias = "eval_legendre")]
#[inline]
pub fn legendre_p<T: LegendrePArg>(n: i32, z: T) -> T {
    z.legendre_p(n as c_int)
}

/// All Legendre polynomials of the 1st kind
///
/// Output length is `n + 1`. The entry at `j` corresponds to degree `j` in `0..=n`.
///
/// Corresponds to [`scipy.special.legendre_p_all`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre_p_all.html
///
/// # See also
/// - [`legendre_p`]
#[doc(alias = "lpn", alias = "clpn")]
#[inline]
pub fn legendre_p_all<T: LegendrePArg>(n: usize, z: T) -> Vec<T> {
    z.legendre_p_all(n)
}

/// Spherical Legendre polynomial of degree `n` and order `m`
///
/// For complex `z` and odd `m` with $|m| \le n$, the factor $\sin z$ is taken as its modulus
/// $|\sin z|$, as in xsf, so the result is not analytic in `z`.
///
/// Corresponds to [`scipy.special.sph_legendre_p`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_legendre_p.html
///
/// # See also
/// - [`sph_legendre_p_all`]: All spherical Legendre polynomials up to degree `n` and order `m`
/// - [`sph_harm_y`](crate::sph_harm_y): Spherical harmonics
#[inline]
pub fn sph_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// All spherical Legendre polynomials of the 1st kind
///
/// Output shape is `(n+1, 2m+1)`. The entry at `(j, i)` corresponds to degree `j` in `0..=n` and
/// order `i` in `0..=m`, or `i - (2m + 1)` in `-m..=-1` for `i > m`.
///
/// For complex `z`, the entries of odd order `k` with $|k| \le j$ take a factor $\sin z$ as its
/// modulus $|\sin z|$, as in xsf, so they are not analytic in `z`.
///
/// Corresponds to [`scipy.special.sph_legendre_p_all`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_legendre_p_all.html
///
/// # See also
/// - [`sph_legendre_p`]
#[inline]
pub fn sph_legendre_p_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.sph_legendre_p_all(n, m)
}

/// Associated Legendre polynomial of the 1st kind
///
/// Corresponds to [`scipy.special.assoc_legendre_p`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_legendre_p.html
///
/// # See also
/// - [`assoc_legendre_p_norm`]
/// - [`assoc_legendre_p_all`]
#[doc(alias = "lpmv")]
#[inline]
pub fn assoc_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2, false)
}

/// All associated Legendre polynomials of the 1st kind
///
/// Corresponds to [`scipy.special.assoc_legendre_p_all`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_legendre_p_all.html
///
/// # See also
/// - [`assoc_legendre_p`]
/// - [`assoc_legendre_p_norm_all`]
#[doc(alias = "lpmn")]
#[inline]
pub fn assoc_legendre_p_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.assoc_legendre_p_all(n, m, 2, false)
}

/// Normalized associated Legendre polynomial of the 1st kind
///
/// Corresponds to [`scipy.special.assoc_legendre_p`][scipy] with `norm=True`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_legendre_p.html
///
/// # See also
/// - [`assoc_legendre_p`]
/// - [`assoc_legendre_p_norm_all`]
#[inline]
pub fn assoc_legendre_p_norm<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.assoc_legendre_p(n as c_int, m as c_int, 2, true)
}

/// All normalized associated Legendre polynomials of the 1st kind
///
/// Corresponds to [`scipy.special.assoc_legendre_p_all`][scipy] with `norm=True`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.assoc_legendre_p_all.html
///
/// # See also
/// - [`assoc_legendre_p`]
/// - [`assoc_legendre_p_norm`]
#[inline]
pub fn assoc_legendre_p_norm_all<T: LegendrePArg>(n: usize, m: usize, z: T) -> Vec<Vec<T>> {
    z.assoc_legendre_p_all(n, m, 2, true)
}

/// All Legendre functions of the 2nd kind and their derivatives
///
/// Output lengths are `n + 1`. The entry at `j` corresponds to degree `j` in `0..=n`.
///
/// Corresponds to [`scipy.special.lqn`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lqn.html
///
/// # See also
/// - [`assoc_legendre_q_all`]
/// - [`legendre_p_all`]
#[doc(alias = "lqn")]
#[inline]
pub fn legendre_q_all<T: LegendreQArg>(n: usize, z: T) -> (Vec<T>, Vec<T>) {
    z.legendre_q_all(n)
}

/// All associated Legendre functions of the 2nd kind and their derivatives
///
/// Computes the associated Legendre function of the second kind of order `m` and degree `n`, `Qmn(z)`,
/// and its derivative, `Qmn'(z)`. Returns two arrays of size `(n+1, m+1)` containing `Qmn(z)` and
/// `Qmn'(z)` for all degrees from `0..=n` and orders from `0..=m`.
///
/// Corresponds to [`scipy.special.lqmn`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lqmn.html
///
/// # See also
/// - [`assoc_legendre_p_all`]
/// - [`legendre_p_all`]
/// - [`legendre_q_all`]
#[doc(alias = "lqmn")]
#[inline]
pub fn assoc_legendre_q_all<T: LegendreQArg>(
    n: usize,
    m: usize,
    z: T,
) -> (Vec<Vec<T>>, Vec<Vec<T>>) {
    z.assoc_legendre_q_all(n, m)
}
