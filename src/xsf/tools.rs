//! The subset of `xsf/tools.h` (xsf v0.2.2) that is needed by `iv_ratio`, `gdtrib`, and the
//! complex `hyp2f1`.
//!
//! Building blocks for implementing special functions
//!
//! The C++ generators (callable objects that produce a value each time they are called) are
//! `FnMut` closures here, and the C++ template parameter `T` of the working type is generic over
//! [`Real`], which is implemented for `f64` and `DoubleDouble`.
//!
//! It also has `int_x86`, which some ports use for undefined `double` -> `int` conversions.

use core::ops::{Add, Div, Mul, Neg, Sub};

use num_complex::Complex;

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

// Series evaluators.
//
// Since these are only used by the complex `hyp2f1`, they are only implemented for complex `T`.
// The `func_name` parameter is omitted, since it is only passed to `set_error`, which is a no-op
// in our build of xsf.

/// Sum an infinite series to a given precision.
///
/// - `g`: a generator of terms for the series.
/// - `init_val`: A starting value that terms are added to. This argument determines the type of
///   the result.
/// - `tol`: relative tolerance for stopping criterion.
/// - `max_terms`: The maximum number of terms to add before giving up and declaring
///   non-convergence.
pub(crate) fn series_eval(
    mut g: impl FnMut() -> Complex<f64>,
    init_val: Complex<f64>,
    tol: f64,
    max_terms: u64,
) -> Complex<f64> {
    let mut result = init_val;
    for _ in 0..max_terms {
        let term = g();
        result += term;
        if term.norm() < result.norm() * tol {
            return result;
        }
    }
    // Exceeded max terms without converging. Return NaN.
    Complex::new(f64::NAN, f64::NAN)
}

/// Sum a fixed number of terms from a series.
///
/// - `g`: a generator of terms for the series.
/// - `init_val`: A starting value that terms are added to. This argument determines the type of
///   the result.
/// - `num_terms`: The number of terms from the series to sum.
pub(crate) fn series_eval_fixed_length(
    mut g: impl FnMut() -> Complex<f64>,
    init_val: Complex<f64>,
    num_terms: u64,
) -> Complex<f64> {
    let mut result = init_val;
    for _ in 0..num_terms {
        result += g();
    }
    result
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

/// Find initial bracket for a bracketing scalar root finder. A valid bracket is a pair of points a < b for
/// which the signs of f(a) and f(b) differ. If f(x0) = 0, where x0 is the initial guess, this bracket finder
/// will return the bracket (x0, x0). It is expected that the rootfinder will check if the bracket
/// endpoints are roots.
///
/// This is a private function intended specifically for the situation where
/// the goal is to invert a CDF function F for a parametrized family of distributions with respect to one
/// parameter, when the other parameters are known, and where F is monotonic with respect to the unknown parameter.
#[allow(clippy::too_many_arguments, clippy::float_cmp)]
pub(crate) fn bracket_root_for_cdf_inversion(
    func: impl Fn(f64) -> f64,
    x0: f64,
    xmin: f64,
    xmax: f64,
    step0_left: f64,
    step0_right: f64,
    factor_left: f64,
    factor_right: f64,
    increasing: bool,
    maxiter: u64,
) -> (f64, f64, f64, f64, i32) {
    let y0 = func(x0);

    if y0 == 0.0 {
        // Initial guess is correct.
        return (x0, x0, y0, y0, 0);
    }

    let y0_sgn = y0.is_sign_negative();

    let search_left;
    /* The frontier is the new leading endpoint of the expanding bracket. The
     * interior endpoint trails behind the frontier. In each step, the old frontier
     * endpoint becomes the new interior endpoint. */
    let (mut interior, mut frontier, mut y_interior, mut y_frontier);
    let (mut y_interior_sgn, mut y_frontier_sgn);
    let (boundary, factor);
    if (increasing && y0 < 0.0) || (!increasing && y0 > 0.0) {
        /* If func is increasing  and func(x_right) < 0 or if func is decreasing and
         *  f(y_right) > 0, we should expand the bracket to the right. */
        (interior, y_interior) = (x0, y0);
        frontier = x0 + step0_right;
        y_interior_sgn = y0_sgn;
        search_left = false;
        boundary = xmax;
        factor = factor_right;
    } else {
        /* Otherwise we move and expand the bracket to the left. */
        (interior, y_interior) = (x0, y0);
        frontier = x0 + step0_left;
        y_interior_sgn = y0_sgn;
        search_left = true;
        boundary = xmin;
        factor = factor_left;
    }

    let mut reached_boundary = false;
    for _ in 0..maxiter {
        y_frontier = func(frontier);
        y_frontier_sgn = y_frontier.is_sign_negative();
        if y_frontier_sgn != y_interior_sgn || (y_frontier == 0.0) {
            /* Stopping condition, func evaluated at endpoints of bracket has opposing signs,
             * meeting requirement for bracketing root finder. (Or endpoint has reached a
             * zero.) */
            if search_left {
                /* Ensure we return an interval (a, b) with a < b. */
                core::mem::swap(&mut interior, &mut frontier);
                core::mem::swap(&mut y_interior, &mut y_frontier);
            }
            return (interior, frontier, y_interior, y_frontier, 0);
        }
        if reached_boundary {
            /* We've reached a boundary point without finding a root . */
            return (
                f64::NAN,
                f64::NAN,
                f64::NAN,
                f64::NAN,
                if search_left { 1 } else { 2 },
            );
        }
        let step = (frontier - interior) * factor;
        interior = frontier;
        y_interior = y_frontier;
        y_interior_sgn = y_frontier_sgn;
        frontier += step;
        if (search_left && frontier <= boundary) || (!search_left && frontier >= boundary) {
            /* If the frontier has reached the boundary, set a flag so the algorithm will know
             * not to search beyond this point. */
            frontier = boundary;
            reached_boundary = true;
        }
    }
    /* Failed to converge within maxiter iterations. If maxiter is sufficiently high and
     * factor_left and factor_right are set appropriately, this should only happen due to
     * a bug in this function. Limiting the number of iterations is a defensive programming measure. */
    (f64::NAN, f64::NAN, f64::NAN, f64::NAN, 3)
}

/// Find root of a scalar function using Chandrupatla's algorithm
#[allow(
    clippy::too_many_arguments,
    clippy::float_cmp,
    clippy::many_single_char_names
)]
pub(crate) fn find_root_chandrupatla(
    func: impl Fn(f64) -> f64,
    mut x1: f64,
    mut x2: f64,
    mut f1: f64,
    mut f2: f64,
    rtol: f64,
    atol: f64,
    maxiter: u64,
) -> (f64, i32) {
    if f1 == 0.0 {
        return (x1, 0);
    }
    if f2 == 0.0 {
        return (x2, 0);
    }
    let mut t = 0.5;
    let (mut x3, mut f3);
    for _ in 0..maxiter {
        let x = x1 + t * (x2 - x1);
        let f = func(x);
        if f.is_sign_negative() == f1.is_sign_negative() {
            x3 = x1;
            x1 = x;
            f3 = f1;
            f1 = f;
        } else {
            x3 = x2;
            x2 = x1;
            x1 = x;
            f3 = f2;
            f2 = f1;
            f1 = f;
        }
        let (xm, fm);
        if f2.abs() < f1.abs() {
            xm = x2;
            fm = f2;
        } else {
            xm = x1;
            fm = f1;
        }
        let tol = 2.0 * rtol * xm.abs() + 0.5 * atol;
        let tl = tol / (x2 - x1).abs();
        if tl > 0.5 || fm == 0.0 {
            return (xm, 0);
        }
        let xi = (x1 - x2) / (x3 - x2);
        let phi = (f1 - f2) / (f3 - f2);
        let fl = 1.0 - (1.0 - xi).sqrt();
        let fh = xi.sqrt();

        if (fl < phi) && (phi < fh) {
            t = (f1 / (f2 - f1)) * (f3 / (f2 - f3))
                + (f1 / (f3 - f1)) * (f2 / (f3 - f2)) * ((x3 - x1) / (x2 - x1));
        } else {
            t = 0.5;
        }
        t = t.max(tl).min(1.0 - tl);
    }
    (f64::NAN, 1)
}

/// `(int)x` as compiled by g++ for x86-64 (`cvttsd2si`), i.e. INT_MIN if `x` is NaN or out of
/// range, where it is undefined behavior in C++
#[allow(clippy::cast_possible_truncation)]
pub(crate) fn int_x86(x: f64) -> i32 {
    if x > -2_147_483_649.0 && x < 2_147_483_648.0 {
        x as i32
    } else {
        i32::MIN
    }
}
