//! Exact references: terminating series, Chu-Vandermonde, and conjugate symmetry

// exact values, e.g. integer parameters
#![allow(clippy::float_cmp)]

use std::ops::RangeInclusive;

use num_complex::{Complex64, c64};
use proptest::prelude::*;

use super::is_pole;
use super::strategies::moderate;
use super::symmetry::close;

/// The largest degree `n` of the terminating series
const MAX_N: u32 = 20;

/// The unit roundoff
const U: f64 = f64::EPSILON / 2.0;

/// Multiples `m / 8` of 1/8, e.g. -3, 2.5 or 0.125, so that `t + k` is exact for `|m| <= 256`
fn eighths(m: RangeInclusive<i32>) -> impl Strategy<Value = f64> {
    m.prop_map(|m| f64::from(m) / 8.0)
}

/// `eighths` for all `|m| <= 256`
fn param() -> impl Strategy<Value = f64> {
    eighths(-256..=256)
}

/// An upper bound on `m u / (1 - m u)`, the relative error of `m` roundings (Higham, Accuracy
/// and Stability of Numerical Algorithms, Lemma 3.1)
fn rounding_bound(m: u32) -> f64 {
    let mu = f64::from(m) * U;
    (mu / (1.0 - mu)).next_up()
}

/// A double-double `hi + lo`, with the algorithms of Joldes, Muller and Popescu, "Tight and
/// rigorous error bounds for basic building blocks of double-word arithmetic" (2017)
#[derive(Clone, Copy)]
struct Dd {
    hi: f64,
    lo: f64,
}

impl Dd {
    /// A bound on the relative error of each operation, above the `3 u^2 + 13 u^3` of `add`
    const ERROR: f64 = 4.0 * U * U;

    /// `a + b` exactly (2Sum)
    fn two_sum(a: f64, b: f64) -> Self {
        let hi = a + b;
        let b_hi = hi - a;
        Self {
            hi,
            lo: (a - (hi - b_hi)) + (b - b_hi),
        }
    }

    /// `a + b` exactly if `|a| >= |b|` (Fast2Sum)
    fn fast_two_sum(a: f64, b: f64) -> Self {
        let hi = a + b;
        Self {
            hi,
            lo: b - (hi - a),
        }
    }

    /// AccurateDWPlusDW
    fn add(self, y: Self) -> Self {
        let s = Self::two_sum(self.hi, y.hi);
        let t = Self::two_sum(self.lo, y.lo);
        let v = Self::fast_two_sum(s.hi, s.lo + t.hi);
        Self::fast_two_sum(v.hi, t.lo + v.lo)
    }

    /// DWTimesFP3, with the exact error of the product by `mul_add`
    fn mul(self, y: f64) -> Self {
        let hi = self.hi * y;
        let lo = self.hi.mul_add(y, -hi);
        Self::fast_two_sum(hi, self.lo.mul_add(y, lo))
    }

    /// DWDivFP3
    fn div(self, y: f64) -> Self {
        let hi = self.hi / y;
        let p_hi = hi * y;
        let p_lo = hi.mul_add(y, -p_hi);
        let delta = ((self.hi - p_hi) - p_lo) + self.lo;
        Self::fast_two_sum(hi, delta / y)
    }

    /// `self * 2^e`, exact if it does not underflow
    fn scale(self, e: i32) -> Self {
        let s = 2_f64.powi(e);
        Self {
            hi: self.hi * s,
            lo: self.lo * s,
        }
    }

    fn abs(self) -> Self {
        let sign = 1_f64.copysign(self.hi);
        Self {
            hi: sign * self.hi,
            lo: sign * self.lo,
        }
    }
}

/// The terminating `pFq(a; b; x)`, and a tolerance for its evaluation
///
/// It terminates at the least `-a_i = n`, so that a `b_j = -m` with `m >= n` is no pole
/// (DLMF 15.2.5). A pole is a reject.
///
/// The tolerance allows the error bound `rounding_bound(n (p + q + 3))` times the sum `A` of the
/// absolute terms of the plain sum of
/// `t_{k+1} = t_k (a_1 + k) ... (a_p + k) x / ((b_1 + k) ... (b_q + k) (k + 1))`, which rounds
/// `p + q + 2` times per term and once per addition (`a_i + k` and `b_j + k` are exact here).
/// It adds one rounding for the returned reference, and the error of the reference: at most
/// `Dd::ERROR A` for each of its `n (p + q + 5)` operations.
///
/// The terms of the reference are double-doubles `s 2^e` with `|s|` in [1/2, 2], so that they do
/// not underflow. It leaves out a term below `u^2 <= u^2 A`, and adding a rescaled term loses
/// less than that.
fn terminating(a: &[f64], b: &[f64], x: f64) -> Result<(f64, f64), TestCaseError> {
    let Some(degree) = a
        .iter()
        .filter(|ai| is_pole(**ai))
        .map(|ai| -ai)
        .reduce(f64::min)
    else {
        return Err(TestCaseError::reject("not a polynomial"));
    };
    if b.iter().any(|bj| is_pole(*bj) && -bj < degree) {
        return Err(TestCaseError::reject("a pole"));
    }
    let n = degree as u32;
    let one = Dd { hi: 1.0, lo: 0.0 };
    let (mut s, mut e, mut sum, mut abs_sum) = (one, 0, one, one);
    for k in 0..n {
        let k = f64::from(k);
        s = a.iter().fold(s, |s, ai| s.mul(ai + k)).mul(x);
        s = b.iter().fold(s, |s, bj| s.div(bj + k)).div(k + 1.0);
        if s.hi == 0.0 {
            // x = 0
            break;
        }
        let shift = s.hi.abs().log2().floor() as i32;
        s = s.scale(-shift);
        e += shift;
        if 2_f64.powi(e + 1) > U * U {
            let t = s.scale(e);
            sum = sum.add(t);
            abs_sum = abs_sum.add(t.abs());
        }
    }
    let (p, q) = (a.len() as u32, b.len() as u32);
    let bound = rounding_bound(n * (p + q + 3) + 1) + Dd::ERROR * f64::from(n * (p + q + 5));
    // `next_up` rounds up the bound, and `A` past the error of `abs_sum`
    Ok((sum.hi, (bound.next_up() * abs_sum.hi.next_up()).next_up()))
}

/// `actual` within `tol` of `expected`, unless `tol` allows no correct digit
fn within<T: Into<Complex64>>(
    (actual, (expected, tol)): (T, (f64, f64)),
) -> Result<(), TestCaseError> {
    prop_assume!(tol < expected.abs(), "no correct digit");
    let actual = actual.into();
    let err = (actual - expected).norm();
    prop_assert!(err <= tol, "{actual:e} != {expected:e} +- {tol:e}");
    Ok(())
}

property!(
    within,
    hyp2f1_terminating(n in 0..=MAX_N, b in param(), c in param(), x in moderate()) => {
        let a = -f64::from(n);
        // known: n >= 10 with b < 0 and x > 1, where x is negated (#2258), and (1 - x)^n for a
        // negative non-integer b = c, n >= 9 and x < -1 (#2262)
        let x = if n >= 10 && b < 0.0 && x > 1.0 { -x } else { x };
        prop_assume!(!(b == c && c < 0.0 && c.fract() != 0.0 && n >= 9 && x < -1.0), "known");
        (xsf::hyp2f1(a, b, c, x), terminating(&[a, b], &[c], x)?)
    }
);
property!(
    within,
    hyp2f1_terminating_complex(n in 0..=MAX_N, b in param(), c in param(), x in moderate()) => {
        let a = -f64::from(n);
        (xsf::hyp2f1(a, b, c, c64(x, 0.0)), terminating(&[a, b], &[c], x)?)
    }
);

/// `a` and `b` without the parameters that they share, which cancel
fn cancel(a: &[f64], b: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    b.retain(|bj| {
        a.iter()
            .position(|ai| ai == bj)
            .map(|i| a.swap_remove(i))
            .is_none()
    });
    (a, b)
}

property!(
    within,
    hyp3f2_terminating(
        n in 0..=MAX_N, a2 in param(), a3 in param(), b1 in param(), b2 in param(), x in moderate()
    ) => {
        let a1 = -f64::from(n);
        let (a, b) = cancel(&[a1, a2, a3], &[b1, b2]);
        (xsf::hyp3f2(a1, a2, a3, b1, b2, x), terminating(&a, &b, x)?)
    }
);

/// Chu-Vandermonde, `2F1(-n, b; c; 1) = (c - b)_n / (c)_n` (DLMF 15.4.24), as a product of `n`
/// quotients of exact factors, and a tolerance for its evaluation
///
/// The tolerance allows the error bound `rounding_bound(2 n)` of this product relative to the
/// exact value, and as much for its evaluation, where `|exact| <= |value| / (1 - gamma_2n)`
/// and `gamma_k / (1 - gamma_k) = gamma_2k / 2`, with `gamma_k = rounding_bound(k)`.
fn chu_vandermonde(n: u32, b: f64, c: f64) -> (f64, f64) {
    let value = (0..n)
        .map(f64::from)
        .map(|k| (c - b + k) / (c + k))
        .product::<f64>();
    (value, (rounding_bound(4 * n) * value.abs()).next_up())
}

property!(
    within,
    hyp2f1_chu_vandermonde(n in 0..=MAX_N, b in param(), c in param())
        if !(is_pole(c) && -c < f64::from(n)) // a pole of (c)_n
        => (xsf::hyp2f1(-f64::from(n), b, c, 1.0), chu_vandermonde(n, b, c))
);

/// `f(conj(z)) = conj(f(z))` exactly, where NaN equals NaN and -0 equals 0
fn conjugate_symmetric([w_conj, w]: [Complex64; 2]) -> Result<(), TestCaseError> {
    let w = w.conj();
    prop_assert!(close(w_conj.re, w.re, 0.0) && close(w_conj.im, w.im, 0.0), "{w_conj:e} != {w:e}");
    Ok(())
}

/// `property!(conjugate_symmetric, ...)` for `xsf::$name(params.., z)` off the real axis, where a
/// branch cut would depend on the sign of zero
macro_rules! conjugate {
    ($name:ident($($param:ident),*)) => {
        property!(
            conjugate_symmetric,
            $name($($param,)* re, im) if im.is_finite() && im != 0.0
                => [xsf::$name($($param,)* c64(re, -im)), xsf::$name($($param,)* c64(re, im))]
        );
    };
}

// not hyp1f1, which rounds the phase of a power of z after adding that of a gamma function, so
// that its asymmetry grows with |z|, and takes time linear in a huge a or b (#2259)
conjugate!(hyp0f1(b));
conjugate!(hyp2f1(a, b, c));
conjugate!(bessel_j(v));
conjugate!(bessel_y(v));
conjugate!(bessel_i(v));
conjugate!(bessel_k(v));
conjugate!(gamma());
conjugate!(loggamma());
conjugate!(rgamma());
conjugate!(digamma());
conjugate!(erf());
conjugate!(erfc());
conjugate!(expi());
conjugate!(exp1());
conjugate!(spence());
// not riemann_zeta, which takes time linear in Im z (#2260)
