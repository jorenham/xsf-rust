//! Translated into pure Rust from `xsf/legendre.h` (except for `lqn` and `lqmn`) and `xsf/recur.h`
//! (xsf v0.2.2), for `T` either `double` or `std::complex<double>`. The dual numbers of
//! `xsf/dual.h`, which xsf uses for derivatives, are not translated.
//!
//! Unlike in xsf:
//! - The integer products in the (normalized) recurrence coefficients, e.g. `(2n+1)(4(n-1)^2-1)`,
//!   are computed in `f64`, which is exact whenever the `int` product in xsf doesn't overflow. In
//!   xsf, they overflow for `n >= 646`, so that e.g. `sph_legendre_p(700, 0, 0.5)` is NaN.
//! - Other `int` arithmetic that can overflow wraps (like C++ does in practice), instead of
//!   panicking in debug builds.
//! - The complex `sin` and `cos` (for complex `theta`) follow glibc only for finite arguments.

use core::f64::consts::PI;
use core::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

use core::ffi::c_int;
use num_complex::Complex64;

use crate::utils;
use crate::xsf::complex::{cdiv, cmul, csqrt};

/// The arithmetic of `T` in C++, i.e. `double` or `std::complex<double>`
pub(crate) trait Num:
    Copy
    + From<f64>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    /// `T * T`
    fn prod(self, rhs: Self) -> Self;
    /// `T / T`
    fn quot(self, rhs: Self) -> Self;
    fn sqrt(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn abs(self) -> f64;
    fn real(self) -> f64;
    fn imag(self) -> f64;
}

impl Num for f64 {
    #[inline]
    fn prod(self, rhs: Self) -> Self {
        self * rhs
    }
    #[inline]
    fn quot(self, rhs: Self) -> Self {
        self / rhs
    }
    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    #[inline]
    fn sin(self) -> Self {
        f64::sin(self)
    }
    #[inline]
    fn cos(self) -> Self {
        f64::cos(self)
    }
    #[inline]
    fn abs(self) -> f64 {
        f64::abs(self)
    }
    #[inline]
    fn real(self) -> f64 {
        self
    }
    #[inline]
    fn imag(self) -> f64 {
        0.0
    }
}

impl Num for Complex64 {
    #[inline]
    fn prod(self, rhs: Self) -> Self {
        cmul(self, rhs)
    }
    #[inline]
    fn quot(self, rhs: Self) -> Self {
        cdiv(self, rhs)
    }
    #[inline]
    fn sqrt(self) -> Self {
        csqrt(self)
    }
    #[inline]
    fn sin(self) -> Self {
        crate::xsf::trig::csin(self)
    }
    #[inline]
    fn cos(self) -> Self {
        crate::xsf::trig::ccos(self)
    }
    #[inline]
    fn abs(self) -> f64 {
        self.norm()
    }
    #[inline]
    fn real(self) -> f64 {
        self.re
    }
    #[inline]
    fn imag(self) -> f64 {
        self.im
    }
}

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
    #[allow(clippy::cast_sign_loss)]
    fn index(&self, i: i32) -> &T {
        &self.data[i as usize]
    }
}

impl<T> IndexMut<i32> for Mdspan<'_, T, 1> {
    #[allow(clippy::cast_sign_loss)]
    fn index_mut(&mut self, i: i32) -> &mut T {
        &mut self.data[i as usize]
    }
}

impl<T> Index<(i32, i32)> for Mdspan<'_, T, 2> {
    type Output = T;
    #[allow(clippy::cast_sign_loss)]
    fn index(&self, (i, j): (i32, i32)) -> &T {
        &self.data[i as usize * self.extents[1] + j as usize]
    }
}

impl<T> IndexMut<(i32, i32)> for Mdspan<'_, T, 2> {
    #[allow(clippy::cast_sign_loss)]
    fn index_mut(&mut self, (i, j): (i32, i32)) -> &mut T {
        &mut self.data[i as usize * self.extents[1] + j as usize]
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub(crate) const fn int(x: usize) -> i32 {
    x as i32
}

// recur.h

fn dot<T: Num, const N: usize>(x: &[T; N], y: &[T; N]) -> T {
    let mut res = T::from(0.0);
    for i in 0..N {
        res = res + x[i].prod(y[i]);
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
    let fac0 = (-T::from(n - 1.0)).quot(T::from(n));
    let fac1 = T::from(2.0 * n - 1.0).quot(T::from(n));

    res[0] = fac0;
    res[1] = fac1.prod(z);
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
        w = (z - T::from(1.0)).sqrt().prod((z + T::from(1.0)).sqrt()); // form for this branch cut
    } else {
        w = -(T::from(1.0) - z.prod(z)).sqrt(); // form for this branch cut
        if m_signbit {
            w = -w;
        }
    }

    if NORM {
        res[0] = T::from(1.0).quot(T::from(2.0).sqrt());
        res[1] = T::from(3.0).sqrt().prod(w).quot(T::from(2.0));
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
        branch_cut_sign.prod(
            T::from((2.0 * m_abs + 1.0) * (2.0 * m_abs - 1.0))
                .quot(T::from(4.0 * m_abs * (m_abs - 1.0)))
                .sqrt(),
        )
    } else if m < 0 {
        branch_cut_sign.quot(T::from((2.0 * m_abs) * (2.0 * m_abs - 2.0)))
    } else {
        branch_cut_sign.prod(T::from((2.0 * m_abs - 1.0) * (2.0 * m_abs - 3.0)))
    };

    res[0] = fac.prod(T::from(1.0) - z.prod(z));
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
        T::from(2.0 * (m_abs + 1.0) - 1.0).quot(T::from(m_abs + 1.0 - f64::from(m)))
    };

    res[0] = res_m_abs_m;
    res[1] = fac.prod(z).prod(res_m_abs_m);
}

fn assoc_legendre_p_recurrence_n<const NORM: bool, T: Num>(m: i32, z: T, n: i32, res: &mut [T; 2]) {
    let (n, m) = (f64::from(n), f64::from(m));
    let (fac0, fac1);
    if NORM {
        let den = T::from((2.0 * n - 3.0) * (n * n - m * m));
        fac0 = -T::from((2.0 * n + 1.0) * ((n - 1.0) * (n - 1.0) - m * m))
            .quot(den)
            .sqrt();
        fac1 = T::from((2.0 * n + 1.0) * (4.0 * (n - 1.0) * (n - 1.0) - 1.0))
            .quot(den)
            .sqrt();
    } else {
        fac0 = (-T::from(n + m - 1.0)).quot(T::from(n - m));
        fac1 = T::from(2.0 * n - 1.0).quot(T::from(n - m));
    }

    res[0] = fac0;
    res[1] = fac1.prod(z);
}

fn assoc_legendre_p_pm1<const NORM: bool, T: Num>(n: i32, m: i32, z: T, res: &mut T) {
    if m == 0 {
        if z.real() >= 0.0 {
            *res = T::from(1.0);
        } else {
            *res = T::from((-1.0_f64).powi(n));
        }
        // Apply normalization for m=0
        if NORM {
            *res = res.prod(T::from(2.0 * f64::from(n) + 1.0).quot(T::from(2.0)).sqrt());
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
#[allow(clippy::float_cmp)]
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
    let fac0 = T::from(1.0).quot(T::from(2.0).prod(T::from(PI).sqrt()));
    let mut fac1 =
        (-T::from(3.0).sqrt()).quot(T::from(2.0).prod(T::from(2.0).prod(T::from(PI)).sqrt()));
    if m_signbit {
        fac1 = -fac1;
    }

    res[0] = fac0;
    res[1] = fac1 * theta_sin.abs();
}

fn sph_legendre_p_recurrence_m_abs_m<T: Num>(theta_sin: T, m: i32, res: &mut [T; 2]) {
    let m_abs = f64::from(m.wrapping_abs());

    let fac = T::from((2.0 * m_abs + 1.0) * (2.0 * m_abs - 1.0))
        .quot(T::from(4.0 * m_abs * (m_abs - 1.0)))
        .sqrt();

    res[0] = fac.prod(theta_sin).prod(theta_sin);
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
    res[1] = fac.prod(theta_cos).prod(res_m_abs_m);
}

fn sph_legendre_p_recurrence_n<T: Num>(m: i32, theta_cos: T, n: i32, res: &mut [T; 2]) {
    let (n, m) = (f64::from(n), f64::from(m));
    let den = T::from((2.0 * n - 3.0) * (n * n - m * m));
    let fac0 = -T::from((2.0 * n + 1.0) * ((n - 1.0) * (n - 1.0) - m * m))
        .quot(den)
        .sqrt();
    let fac1 = T::from((2.0 * n + 1.0) * (4.0 * (n - 1.0) * (n - 1.0) - 1.0))
        .quot(den)
        .sqrt();

    res[0] = fac0;
    res[1] = fac1.prod(theta_cos);
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

pub(crate) fn sph_legendre_p_for_each_n<T: Num>(
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
    sph_legendre_p_for_each_m_abs_m(
        m.wrapping_neg(),
        theta,
        &mut res_m_abs_m,
        |m, res_m_abs_m| {
            res[0] = res_m_abs_m[1];

            sph_legendre_p_for_each_n_from(n, m, theta, res_m_abs_m[1], res, |n, res_n| {
                f(n, m, res_n);
            });
        },
    );
}

fn sph_legendre_p_t<T: Num>(n: i32, m: i32, theta: T) -> T {
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

// Rust API

pub trait LegendrePArg: crate::sealed::Sealed + Sized {
    fn legendre_p(self, n: c_int) -> Self;
    fn legendre_p_all(self, n: usize) -> Vec<Self>;
    fn sph_legendre_p(self, n: c_int, m: c_int) -> Self;
    fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>>;
    fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self;
    fn assoc_legendre_p_all(self, n: usize, m: usize, bc: c_int, norm: bool) -> Vec<Vec<Self>>;
}

fn legendre_p_all_vec<T: Num>(z: T, n: usize) -> Vec<T> {
    let mut pn = vec![T::from(f64::NAN); n + 1];
    legendre_p_all_t(z, Mdspan::new(&mut pn, [n + 1]));
    pn
}

fn sph_legendre_p_all_vec<T: Num>(z: T, n: usize, m: usize) -> Vec<Vec<T>> {
    let extents = [n + 1, 2 * m + 1];
    let mut pnm = vec![T::from(f64::NAN); extents[0] * extents[1]];
    sph_legendre_p_all_t(z, Mdspan::new(&mut pnm, extents));
    utils::vec_into_vecvec(pnm, extents[0], extents[1], false)
}

fn assoc_legendre_p_dyn<T: Num>(z: T, n: c_int, m: c_int, bc: c_int, norm: bool) -> T {
    if norm {
        assoc_legendre_p_t::<true, T>(n, m, z, bc)
    } else {
        assoc_legendre_p_t::<false, T>(n, m, z, bc)
    }
}

fn assoc_legendre_p_all_vec<T: Num>(
    z: T,
    n: usize,
    m: usize,
    bc: c_int,
    norm: bool,
) -> Vec<Vec<T>> {
    let extents = [n + 1, 2 * m + 1];
    let mut pnm = vec![T::from(f64::NAN); extents[0] * extents[1]];
    let res = Mdspan::new(&mut pnm, extents);
    if norm {
        assoc_legendre_p_all_t::<true, T>(z, bc, res);
    } else {
        assoc_legendre_p_all_t::<false, T>(z, bc, res);
    }
    utils::vec_into_vecvec(pnm, extents[0], extents[1], false)
}

macro_rules! impl_legendre_p_arg {
    ($T:ty) => {
        impl LegendrePArg for $T {
            #[inline]
            fn legendre_p(self, n: c_int) -> Self {
                legendre_p_t(n, self)
            }

            #[inline]
            fn legendre_p_all(self, n: usize) -> Vec<Self> {
                legendre_p_all_vec(self, n)
            }

            #[inline]
            fn sph_legendre_p(self, n: c_int, m: c_int) -> Self {
                sph_legendre_p_t(n, m, self)
            }

            #[inline]
            fn sph_legendre_p_all(self, n: usize, m: usize) -> Vec<Vec<Self>> {
                sph_legendre_p_all_vec(self, n, m)
            }

            #[inline]
            fn assoc_legendre_p(self, n: c_int, m: c_int, bc: c_int, norm: bool) -> Self {
                assoc_legendre_p_dyn(self, n, m, bc, norm)
            }

            #[inline]
            fn assoc_legendre_p_all(
                self,
                n: usize,
                m: usize,
                bc: c_int,
                norm: bool,
            ) -> Vec<Vec<Self>> {
                assoc_legendre_p_all_vec(self, n, m, bc, norm)
            }
        }
    };
}

impl_legendre_p_arg!(f64);
impl_legendre_p_arg!(Complex64);

pub trait LegendreQArg: crate::sealed::Sealed + Sized {
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>);
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>);
}

impl LegendreQArg for f64 {
    #[inline]
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
        // like in scipy, since `lqn` writes at least 2 terms
        let n1 = n.max(1);

        let mut qn = vec![f64::NAN; n1 + 1];
        let mut qd = vec![f64::NAN; n1 + 1];

        unsafe {
            crate::ffi::xsf::lqn(n1, self, qn.as_mut_ptr(), qd.as_mut_ptr());
        }
        qn.truncate(n + 1);
        qd.truncate(n + 1);
        (qn, qd)
    }

    #[inline]
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>) {
        // like in scipy, since `lqmn` writes at least 2 x 2 terms
        let (mm, nn) = (m.max(1), n.max(1));

        let (ni, nj) = (mm + 1, nn + 1);
        let mut qm = vec![f64::NAN; ni * nj];
        let mut qd = vec![f64::NAN; ni * nj];

        unsafe {
            crate::ffi::xsf::lqmn(mm, nn, self, qm.as_mut_ptr(), qd.as_mut_ptr());
        }

        (
            utils::vecvec_truncate(utils::vec_to_vecvec(&qm, ni, nj, true), n + 1, m + 1),
            utils::vecvec_truncate(utils::vec_to_vecvec(&qd, ni, nj, true), n + 1, m + 1),
        )
    }
}

impl LegendreQArg for num_complex::Complex<f64> {
    #[inline]
    fn legendre_q_all(self, n: usize) -> (Vec<Self>, Vec<Self>) {
        // like in scipy, since `lqn` writes at least 2 terms
        let n1 = n.max(1);

        let mut cqn = vec![num_complex::Complex::new(f64::NAN, 0.0); n1 + 1];
        let mut cqd = vec![num_complex::Complex::new(f64::NAN, 0.0); n1 + 1];

        unsafe {
            crate::ffi::xsf::lqn_1(n1, self, cqn.as_mut_ptr(), cqd.as_mut_ptr());
        }
        cqn.truncate(n + 1);
        cqd.truncate(n + 1);
        (cqn, cqd)
    }

    #[inline]
    fn assoc_legendre_q_all(self, n: usize, m: usize) -> (Vec<Vec<Self>>, Vec<Vec<Self>>) {
        // like in scipy, since `lqmn` writes at least 2 x 2 terms
        let (mm, nn) = (m.max(1), n.max(1));

        let (ni, nj) = (mm + 1, nn + 1);
        let mut cqm = vec![num_complex::Complex::new(f64::NAN, 0.0); ni * nj];
        let mut cqd = vec![num_complex::Complex::new(f64::NAN, 0.0); ni * nj];

        unsafe {
            crate::ffi::xsf::lqmn_1(mm, nn, self, cqm.as_mut_ptr(), cqd.as_mut_ptr());
        }

        (
            utils::vecvec_truncate(utils::vec_into_vecvec(cqm, ni, nj, true), n + 1, m + 1),
            utils::vecvec_truncate(utils::vec_into_vecvec(cqd, ni, nj, true), n + 1, m + 1),
        )
    }
}

/// Legendre polynomial of degree `n`
///
/// Corresponds to [`scipy.special.eval_legendre`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_legendre.html
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
/// Corresponds to [`scipy.special.sph_legendre_p`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_legendre_p.html
#[inline]
pub fn sph_legendre_p<T: LegendrePArg>(n: i32, m: i32, z: T) -> T {
    z.sph_legendre_p(n as c_int, m as c_int)
}

/// All spherical Legendre polynomials of the 1st kind
///
/// Output shape is `(n+1, 2m+1)`. The entry at `(j, i)` corresponds to degree `j` in `0..=n` and
/// order `i` in `-m..=m`.
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
/// Corresponds to [`scipy.special.assoc_legendre_p`][scipy].
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
/// Corresponds to [`scipy.special.assoc_legendre_p_all`][scipy].
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

/// All Legendre polynomials of the 2nd kind and their derivatives
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

/// All associated Legendre polynomials of the 2nd kind and their derivatives
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

    #[test]
    fn test_legendre_q_all_n0() {
        // computed with n = 1, since `lqn` writes at least 2 terms
        let (qn, qd) = crate::legendre_q_all(1, 1.5);
        assert_eq!(
            crate::legendre_q_all(0, 1.5),
            (qn[..1].to_vec(), qd[..1].to_vec())
        );

        let z = c64(1.5, 0.0);
        let (qn, qd) = crate::legendre_q_all(1, z);
        assert_eq!(
            crate::legendre_q_all(0, z),
            (qn[..1].to_vec(), qd[..1].to_vec())
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

    #[test]
    fn test_assoc_legendre_q_all_n0_m0() {
        // computed with n = 1 and m = 1, since `lqmn` writes at least 2 x 2 terms
        macro_rules! check {
            ($z:expr) => {
                let (qm, qd) = crate::assoc_legendre_q_all(1, 1, $z);
                for (n, m) in [(0, 0), (0, 1), (1, 0)] {
                    let sub = |q: &[Vec<_>]| q[..=n].iter().map(|r| r[..=m].to_vec()).collect();
                    assert_eq!(crate::assoc_legendre_q_all(n, m, $z), (sub(&qm), sub(&qd)));
                }
            };
        }
        check!(1.5);
        check!(0.5);
        check!(c64(1.5, 0.0));
        check!(c64(0.3, 0.4));
    }

    // deviations from xsf

    #[test]
    fn test_large_n_normalized() {
        // NaN in xsf because of `int` overflow, for `n >= 646`
        let p = crate::sph_legendre_p(700, 0, 0.5);
        let expected = -0.337_154_106_114_745_7; // mpmath
        assert!(((p - expected) / expected).abs() < 1e-13, "{p}");

        let p = crate::assoc_legendre_p_norm(700, 0, 0.3);
        let expected = 0.802_058_211_710_411_2; // mpmath
        assert!(((p - expected) / expected).abs() < 1e-13, "{p}");
    }
}
