//! Translated into pure Rust from `xsf/iv_ratio.h` (xsf v0.2.2, and `iv_ratioinv` from xsf main).
//!
//! Numerically stable computation of iv(v, x) / iv(v-1, x) and its inverse
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use core::ops::Mul;

use crate::xsf::cephes::dd_real::{DoubleDouble, frexp_f64, ldexp_f64};
use crate::xsf::complex::atanh;
use crate::xsf::tools::{
    Real, bracket_root_for_cdf_inversion, continued_fraction_series, find_root_chandrupatla,
    series_eval_kahan,
};

/// Generates the "tail" of Perron's continued fraction for iv(v,x)/iv(v-1,x).
///
/// The Perron continued fraction is studied in \[1\].  It is given by
///
/// ```text
///         iv(v, x)      x    -(2v+1)x   -(2v+3)x   -(2v+5)x
///   R := --------- = ------ ---------- ---------- ---------- ...
///        iv(v-1,x)   x+2v + 2(v+x)+1 + 2(v+x)+2 + 2(v+x)+3 +
/// ```
///
/// Given a suitable constant c, the continued fraction may be rearranged
/// into the following form to avoid premature floating point overflow:
///
/// ```text
///        xc                -(2vc+c)(xc) -(2vc+3c)(xc) -(2vc+5c)(xc)
///   R = -----,  fc = 2vc + ------------ ------------- ------------- ...
///       xc+fc              2(vc+xc)+c + 2(vc+xc)+2c + 2(vc+xc)+3c +
/// ```
///
/// This class generates the fractions of fc after 2vc.
///
/// \[1\] Gautschi, W. and Slavik, J. (1978). "On the computation of modified
///     Bessel function ratios." Mathematics of Computation, 32(143):865-875.
struct IvRatioCFTailGenerator<T> {
    a0: T,
    as_: T, // a[k] == a0 + as*k, k >= 1
    b0: T,
    bs: T,  // b[k] == b0 + bs*k, k >= 1
    k: u64, // current index
}

impl<T: Real> IvRatioCFTailGenerator<T>
where
    f64: Mul<T, Output = T>,
{
    fn new(vc: T, xc: T, c: T) -> Self {
        Self {
            a0: -(2.0 * vc - c) * xc,
            as_: -2.0 * c * xc,
            b0: 2.0 * (vc + xc),
            bs: c,
            k: 0,
        }
    }

    /// `operator()`
    #[allow(clippy::cast_precision_loss)]
    fn call(&mut self) -> (T, T) {
        self.k += 1;
        (
            T::from_f64(self.k as f64).fma(self.as_, self.a0),
            T::from_f64(self.k as f64).fma(self.bs, self.b0),
        )
    }
}

/// Computes f(v, x) using Perron's continued fraction.
///
/// T specifies the working type.  This allows the function to perform
/// calculations in a higher precision, such as double-double, even if
/// the return type is hardcoded to be double.
fn iv_ratio_cf<T: Real>(v: f64, x: f64, complement: bool) -> (f64, u64)
where
    f64: Mul<T, Output = T>,
{
    let (_, e) = frexp_f64(v.max(x));
    let c = T::from_f64(ldexp_f64(1.0, 2 - e)); // rescaling multiplier
    let vc = v * c;
    let xc = x * c;

    let mut cf = IvRatioCFTailGenerator::new(vc, xc, c);
    let mut g = continued_fraction_series(|| cf.call());
    let (fc, terms) = series_eval_kahan(|| g.call(), T::from_f64(f64::EPSILON), 1000, 2.0 * vc);

    let ret = (if complement { fc } else { xc }) / (xc + fc);
    (ret.to_f64(), terms)
}

/// Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind
#[must_use]
#[allow(clippy::float_cmp)]
pub fn iv_ratio(v: f64, x: f64) -> f64 {
    if v.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if v < 0.5 || x < 0.0 {
        return f64::NAN;
    }
    if v.is_infinite() && x.is_infinite() {
        // There is not a unique limit as both v and x tends to infinity.
        return f64::NAN;
    }
    if x == 0.0 {
        return x; // keep sign of x because iv_ratio is an odd function
    }
    if v.is_infinite() {
        return 0.0;
    }
    if x.is_infinite() {
        return 1.0;
    }
    if v == 0.5 {
        // Closed-form solution for v = 0.5: iv_ratio(0.5, x) = tanh(x)
        return x.tanh();
    }

    let (ret, terms) = iv_ratio_cf::<f64>(v, x, false);
    if terms == 0 {
        // failed to converge; should not happen
        return f64::NAN;
    }
    ret
}

/// Compute `1 - iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind
#[must_use]
#[allow(clippy::float_cmp)]
pub fn iv_ratio_c(v: f64, x: f64) -> f64 {
    if v.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if v < 0.5 || x < 0.0 {
        return f64::NAN;
    }
    if v.is_infinite() && x.is_infinite() {
        // There is not a unique limit as both v and x tends to infinity.
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if v.is_infinite() {
        return 1.0;
    }
    if x.is_infinite() {
        return 0.0;
    }

    if v >= 1.0 {
        // Numerical experiments show that evaluating the Perron c.f.
        // in double precision is sufficiently accurate if v >= 1.
        let (ret, terms) = iv_ratio_cf::<f64>(v, x, true);
        if terms == 0 {
            // failed to converge; should not happen
            return f64::NAN;
        }
        ret
    } else if v > 0.5 {
        // double-double arithmetic is needed for 0.5 < v < 1 to
        // achieve relative error on the scale of machine precision.
        let (ret, terms) = iv_ratio_cf::<DoubleDouble>(v, x, true);
        if terms == 0 {
            // failed to converge; should not happen
            return f64::NAN;
        }
        ret
    } else {
        // The previous branch (v > 0.5) also works for v == 0.5, but
        // the closed-form formula "1 - tanh(x)" is more efficient.
        let t = (-2.0 * x).exp();
        (2.0 * t) / (1.0 + t)
    }
}

/// Inverse of [`iv_ratio`] with respect to `x`
///
/// Compute x such that `iv(v,x)/iv(v-1,x) = r` for v >= 0.5 and 0 <= r <= 1.
///
/// Corresponds to `xsf::iv_ratioinv` in xsf main.
#[must_use]
#[allow(clippy::float_cmp)]
pub fn iv_ratioinv(v: f64, r: f64) -> f64 {
    if v.is_nan() || r.is_nan() {
        return f64::NAN;
    }
    if !v.is_finite() || v < 0.5 || r < 0.0 || r > 1.0 {
        // iv_ratioinv is only defined for v >= 0.5
        return f64::NAN;
    }
    if r == 0.0 {
        return 0.0;
    }
    if r == 1.0 {
        return f64::INFINITY;
    }
    if v == 0.5 {
        return atanh(r);
    }
    // Algorithm description: use Chandrupatla's method to find the root of
    // f(x) = iv_ratio(v, x) - r (or f(x) = 1 - iv_ratio_c(v, x) - r if r > 0.5).
    // Bounds from Amos (1974), "Computation of Modified Bessel Functions and Their Ratios".
    // That paper defines the ratio as iv(v+1, x)/iv(v, x), but we define it as iv(v, x)/iv(v-1, x).
    // Substituting v -> v-1 in eq. (9) and eq. (16) yields the following bounds for x:
    //   x/(v-1/2 + sqrt(x^2 + (v+1/2)^2)) <= r <= x/(v-1 + sqrt(x^2 + (v+1)^2))
    //
    // Inverting (inequality direction flips):
    //   lower x (from upper r-bound): r*[(v-1) + sqrt((v-1)^2 + 4v*(1-r^2))] / (1-r^2)
    //   upper x (from lower r-bound): r*[(v-0.5) + sqrt((v-0.5)^2 + 2v*(1-r^2))] / (1-r^2)

    let r_c = 1.0 - r;
    let one_minus_r_sq = r_c * (1.0 + r); // 1 - r^2

    let vm1 = v - 1.0;
    let lower_bound = r * (vm1 + (vm1 * vm1 + 4.0 * v * one_minus_r_sq).sqrt()) / one_minus_r_sq;

    let vm05 = v - 0.5;
    let mut upper_bound =
        r * (vm05 + (vm05 * vm05 + 2.0 * v * one_minus_r_sq).sqrt()) / one_minus_r_sq;

    // For small r both bounds converge to 2v*r (== true x). Ensure bracket has width.
    if upper_bound <= lower_bound {
        upper_bound = 2.0 * lower_bound;
    }

    let func = |x: f64| {
        if r <= 0.5 {
            return iv_ratio(v, x) - r;
        }
        r_c - iv_ratio_c(v, x)
    };

    let mut xl = lower_bound;
    let mut xr = upper_bound;
    let mut f_xl = func(xl);
    let mut f_xr = func(xr);

    if f_xl * f_xr > 0.0 {
        // The Amos bounds are theoretically correct but may yield an invalid bracket
        // due to precision issues in iv_ratio or iv_ratio_c. This happens especially for very small or large r.
        // Fallback: use bracket_root_for_cdf_inversion to find a valid bracket.
        // Bracketing parameters taken from gdtrib, only difference: our function is increasing
        let (b_xl, b_xr, b_f_xl, b_f_xr, bracket_status) = bracket_root_for_cdf_inversion(
            func,
            1.0,
            f64::MIN_POSITIVE,
            f64::MAX,
            -0.875,
            7.0,
            0.125,
            8.0,
            true,
            342,
        );
        if bracket_status == 1 {
            return 0.0;
        }
        if bracket_status == 2 {
            return f64::INFINITY;
        }
        if bracket_status >= 3 {
            return f64::NAN;
        }
        xl = b_xl;
        xr = b_xr;
        f_xl = b_f_xl;
        f_xr = b_f_xr;
    }

    let (result, root_status) =
        find_root_chandrupatla(func, xl, xr, f_xl, f_xr, f64::EPSILON, 1e-308, 1000);
    if root_status != 0 {
        // Root finding failed. This should never happen.
        return f64::NAN;
    }
    result
}
