//! Translated into pure Rust from `xsf/iv_ratio.h` (xsf v0.2.2).
//!
//! Numerically stable computation of iv(v, x) / iv(v-1, x)
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the first line of the header comment is corrected; it said `iv(v+1, x) / iv(v, x)`.

use core::ops::Mul;

use crate::xsf::cephes::dd_real::{DoubleDouble, frexp_f64, ldexp_f64};
use crate::xsf::tools::{Real, continued_fraction_series, series_eval_kahan};

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

#[cfg(test)]
mod tests {
    #[test]
    fn test_iv_ratio() {
        xsref::test("iv_ratio", "d_d-d", |x| crate::iv_ratio(x[0], x[1]));
        assert_eq!(crate::iv_ratio(0.5, 3.0), 3.0_f64.tanh());
    }

    #[test]
    fn test_iv_ratio_c() {
        xsref::test("iv_ratio_c", "d_d-d", |x| crate::iv_ratio_c(x[0], x[1]));
    }
}
