//! The subset of `xsf/tools.h` (xsf v0.2.2) that is needed by `iv_ratio`.
//!
//! Building blocks for implementing special functions
//!
//! The C++ generators (callable objects that produce a value each time they are called) are
//! `FnMut` closures here, and the C++ template parameter `T` of the working type is generic over
//! [`Real`], which is implemented for `f64` and `DoubleDouble`.

use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::xsf::cephes::dd_real::DoubleDouble;

/// The working type of the series and continued fraction evaluators
pub(crate) trait Real:
    Copy
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    /// `T(x)`
    fn from_f64(x: f64) -> Self;
    /// `static_cast<double>(x)`
    fn to_f64(self) -> f64;
    /// `abs(x)`
    fn abs(self) -> Self;
    /// `fma(x, y, z)`
    fn fma(self, y: Self, z: Self) -> Self;
}

impl Real for f64 {
    #[inline]
    fn from_f64(x: f64) -> Self {
        x
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self
    }

    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }

    #[inline]
    fn fma(self, y: Self, z: Self) -> Self {
        self.mul_add(y, z)
    }
}

impl Real for DoubleDouble {
    #[inline]
    fn from_f64(x: f64) -> Self {
        Self::new(x)
    }

    #[inline]
    fn to_f64(self) -> f64 {
        f64::from(self)
    }

    #[inline]
    fn abs(self) -> Self {
        DoubleDouble::abs(self)
    }

    #[inline]
    fn fma(self, y: Self, z: Self) -> Self {
        DoubleDouble::fma(self, y, z)
    }
}

/* Performs one step of Kahan summation. */
#[inline]
fn kahan_step<T: Real>(sum: &mut T, comp: &mut T, x: T) {
    let y = x - *comp;
    let t = *sum + y;
    *comp = (t - *sum) - y;
    *sum = t;
}

/// Evaluates an infinite series using Kahan summation.
///
/// Denote the series by
///
/// ```text
///   S = a[0] + a[1] + a[2] + ...
/// ```
///
/// And for n = 0, 1, 2, ..., denote its n-th partial sum by
///
/// ```text
///   S[n] = a[0] + a[1] + ... + a[n]
/// ```
///
/// This function computes S\[0\], S\[1\], ... until a\[n\] is sufficiently
/// small or if the maximum number of terms have been evaluated.
///
/// Parameters
/// ----------
/// - `g`:
///   Generator that yields the sequence of values a\[1\], a\[2\], a\[3\], ...
/// - `tol`:
///   Relative tolerance for convergence.  Specifically, stop iteration
///   as soon as `abs(a[n]) <= tol * abs(S[n])` for some n >= 1.
/// - `max_terms`:
///   Maximum number of terms after a\[0\] to evaluate.  It should be set
///   large enough such that the convergence criterion is guaranteed
///   to have been satisfied within that many terms if there is no
///   rounding error.
/// - `init_val`:
///   a\[0\].  The type of this parameter (T) is used for intermediary computations as well as
///   the result.
///
/// Return Value
/// ------------
/// If the convergence criterion is satisfied by some `n <= max_terms`,
/// returns `(S[n], n)`.  Otherwise, returns `(S[max_terms], 0)`.
pub(crate) fn series_eval_kahan<T: Real>(
    mut g: impl FnMut() -> T,
    tol: T,
    max_terms: u64,
    init_val: T,
) -> (T, u64) {
    let mut sum = init_val;
    let mut comp = T::from_f64(0.0);
    for i in 0..max_terms {
        let term = g();
        kahan_step(&mut sum, &mut comp, term);
        if term.abs() <= tol * sum.abs() {
            return (sum, i + 1);
        }
    }
    (sum, 0)
}

/// Generator that yields the difference of successive convergents of a
/// continued fraction.
///
/// Let f\[n\] denote the n-th convergent of a continued fraction:
///
/// ```text
///                 a[1]   a[2]       a[n]
///   f[n] = b[0] + ------ ------ ... ----
///                 b[1] + b[2] +     b[n]
/// ```
///
/// with f\[0\] = b\[0\].  This generator yields the sequence of values
/// f\[1\]-f\[0\], f\[2\]-f\[1\], f\[3\]-f\[2\], ...
///
/// Constructor Arguments
/// ---------------------
/// - `cf`:
///   Generator that yields the terms of the continued fraction as (numerator, denominator)
///   pairs, starting from (a\[1\], b\[1\]).
///
///   The constructed object always eagerly retrieves the next term
///   of the continued fraction.  Specifically, (a\[1\], b\[1\]) is
///   retrieved upon construction, and (a\[n\], b\[n\]) is retrieved after
///   (n-1) calls of `call()`.
///
/// Type Arguments
/// --------------
/// - `T`:
///   Type in which computations are performed and results are turned.
///
/// Remarks
/// -------
/// The series is computed using the recurrence relation described in \[1\].
/// Let v\[n\], n >= 1 denote the terms of the series.  Then
///
/// ```text
///   v[1] = a[1] / b[1]
///   v[n] = v[n-1] * r[n-1], n >= 2
/// ```
///
/// where
///
/// ```text
///                              -(a[n] + a[n] * r[n-1])
///   r[1] = 0,  r[n] = ------------------------------------------, n >= 2
///                      (a[n] + a[n] * r[n-1]) + (b[n] * b[n-1])
/// ```
///
/// No error checking is performed.  The caller must ensure that all terms
/// are finite and that intermediary computations do not trigger floating
/// point exceptions such as overflow.
///
/// The numerical stability of this method depends on the characteristics
/// of the continued fraction being evaluated.
///
/// Reference
/// ---------
/// \[1\] Gautschi, W. (1967). “Computational Aspects of Three-Term
///     Recurrence Relations.” SIAM Review, 9(1):24-82.
pub(crate) struct ContinuedFractionSeriesGenerator<G, T> {
    cf: G, // continued fraction generator
    v: T,  // v[n] == f[n] - f[n-1], n >= 1
    r: T,  // r[1] = 0, r[n] = v[n]/v[n-1], n >= 2
    b: T,  // last denominator, i.e. b[n-1]
}

impl<G: FnMut() -> (T, T), T: Real> ContinuedFractionSeriesGenerator<G, T> {
    fn new(mut cf: G) -> Self {
        // init()
        let (num, denom) = cf();
        let a = num;
        let b = denom;
        Self {
            cf,
            v: a / b,
            r: T::from_f64(0.0),
            b,
        }
    }

    /// `operator()`
    pub(crate) fn call(&mut self) -> T {
        let v = self.v;
        self.advance();
        v
    }

    fn advance(&mut self) {
        let (num, denom) = (self.cf)();
        let a = num;
        let b = denom;
        let p = a + a * self.r;
        let q = p + b * self.b;
        self.r = -p / q;
        self.v = self.v * self.r;
        self.b = b;
    }
}

/// Converts a continued fraction into a series whose terms are the
/// difference of its successive convergents.
///
/// See [`ContinuedFractionSeriesGenerator`] for details.
pub(crate) fn continued_fraction_series<G: FnMut() -> (T, T), T: Real>(
    cf: G,
) -> ContinuedFractionSeriesGenerator<G, T> {
    ContinuedFractionSeriesGenerator::new(cf)
}
