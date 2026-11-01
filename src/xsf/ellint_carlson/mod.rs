//! Carlson symmetric elliptic integrals, translated into pure Rust from
//! `xsf/cpu/ellint_carlson.h` and `xsf/ellint_carlson/*.hh` (xsf main).
//!
//! The C++ functions return an `ExitStatus` and write the result to an output argument; here they
//! return both as a tuple. `set_error` is a no-op in our build of xsf, so the status is discarded.

mod argcheck;
mod arith_aux;
mod arithmetic;
mod common;
mod rc;
mod rd;
mod rf;
mod rg;
mod rj;
mod typing;

use num_complex::Complex64;

const ELLINT_CARLSON_RERR: f64 = 5e-16;

pub trait CarlsonArg: crate::sealed::Sealed + Sized {
    fn elliprc(self, y: Self) -> Self;
    fn elliprd(self, y: Self, z: Self) -> Self;
    fn elliprf(self, y: Self, z: Self) -> Self;
    fn elliprg(self, y: Self, z: Self) -> Self;
    fn elliprj(self, y: Self, z: Self, p: Self) -> Self;
}

macro_rules! impl_carlson_arg {
    ($T:ty) => {
        impl CarlsonArg for $T {
            #[inline]
            fn elliprc(self, y: Self) -> Self {
                rc::rc(self, y, ELLINT_CARLSON_RERR).1
            }
            #[inline]
            fn elliprd(self, y: Self, z: Self) -> Self {
                rd::rd(self, y, z, ELLINT_CARLSON_RERR).1
            }
            #[inline]
            fn elliprf(self, y: Self, z: Self) -> Self {
                rf::rf(self, y, z, ELLINT_CARLSON_RERR).1
            }
            #[inline]
            fn elliprg(self, y: Self, z: Self) -> Self {
                rg::rg(self, y, z, ELLINT_CARLSON_RERR).1
            }
            #[inline]
            fn elliprj(self, y: Self, z: Self, p: Self) -> Self {
                rj::rj(self, y, z, p, ELLINT_CARLSON_RERR, false).1
            }
        }
    };
}

impl_carlson_arg!(f64);
impl_carlson_arg!(Complex64);

/// Degenerate symmetric elliptic integral $R_C(x, y)$ for real or complex input
///
/// $$
/// R_C(x, y) = {1 \over 2} \int_0^\infty (t + x)^{-1/2} (t + y)^{-1} \dd t
/// = R_F(x, y, y)
/// $$
///
/// For $y < 0$, the Cauchy principal value is returned.
///
/// Corresponds to [`scipy.special.elliprc`][elliprc] in SciPy
///
/// [elliprc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprc.html
///
/// # See also
/// - [`elliprf`]: Completely-symmetric elliptic integral of the first kind
#[must_use]
#[inline]
pub fn elliprc<T: CarlsonArg>(x: T, y: T) -> T {
    x.elliprc(y)
}

/// Symmetric elliptic integral of the second kind $R_D(x, y, z)$ for real or complex input
///
/// $$
/// R_D(x, y, z) = {3 \over 2} \int_0^\infty
/// \left[(t + x) (t + y)\right]^{-1/2} (t + z)^{-3/2} \dd t
/// = R_J(x, y, z, z)
/// $$
///
/// Corresponds to [`scipy.special.elliprd`][elliprd] in SciPy
///
/// [elliprd]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprd.html
///
/// # See also
/// - [`elliprj`]: Symmetric elliptic integral of the third kind
#[must_use]
#[inline]
pub fn elliprd<T: CarlsonArg>(x: T, y: T, z: T) -> T {
    x.elliprd(y, z)
}

/// Completely-symmetric elliptic integral of the first kind $R_F(x, y, z)$ for real or complex
/// input
///
/// $$
/// R_F(x, y, z) = {1 \over 2} \int_0^\infty \left[(t + x) (t + y) (t + z)\right]^{-1/2} \dd t
/// $$
///
/// Corresponds to [`scipy.special.elliprf`][elliprf] in SciPy
///
/// [elliprf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprf.html
///
/// # See also
/// - [`elliprg`]: Completely-symmetric elliptic integral of the second kind
/// - [`ellipk`](crate::ellipk): Complete elliptic integral of the first kind, $K(m) = R_F(0, 1 - m,
///   1)$
#[must_use]
#[inline]
pub fn elliprf<T: CarlsonArg>(x: T, y: T, z: T) -> T {
    x.elliprf(y, z)
}

/// Completely-symmetric elliptic integral of the second kind $R_G(x, y, z)$ for real or complex
/// input
///
/// $$
/// R_G(x, y, z) = {1 \over 4} \int_0^\infty \left[(t + x) (t + y) (t + z)\right]^{-1/2}
/// \left({x \over t + x} + {y \over t + y} + {z \over t + z}\right) t \dd t
/// $$
///
/// Corresponds to [`scipy.special.elliprg`][elliprg] in SciPy
///
/// [elliprg]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprg.html
///
/// # See also
/// - [`elliprf`]: Completely-symmetric elliptic integral of the first kind
/// - [`ellipe`](crate::ellipe): Complete elliptic integral of the second kind,
///   $E(m) = 2 R_G(0, 1 - m, 1)$
#[must_use]
#[inline]
pub fn elliprg<T: CarlsonArg>(x: T, y: T, z: T) -> T {
    x.elliprg(y, z)
}

/// Symmetric elliptic integral of the third kind $R_J(x, y, z, p)$ for real or complex input
///
/// $$
/// R_J(x, y, z, p) = {3 \over 2} \int_0^\infty
/// \left[(t + x) (t + y) (t + z)\right]^{-1/2} (t + p)^{-1} \dd t
/// $$
///
/// For real $x, y, z \ge 0$ and $p < 0$, the Cauchy principal value is returned.
///
/// Corresponds to [`scipy.special.elliprj`][elliprj] in SciPy
///
/// [elliprj]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprj.html
///
/// # See also
/// - [`elliprd`]: Symmetric elliptic integral of the second kind
#[must_use]
#[inline]
pub fn elliprj<T: CarlsonArg>(x: T, y: T, z: T, p: T) -> T {
    x.elliprj(y, z, p)
}
