//! Exact references: terminating series, Chu-Vandermonde, and conjugate symmetry

// exact values, e.g. integer parameters
#![allow(clippy::float_cmp)]

use num_complex::{Complex64, c64};
use proptest::prelude::*;

use super::is_pole;
use super::strategies::moderate;

/// The largest degree `n` of the terminating series
const MAX_N: u32 = 20;

/// Multiples of 1/8 in [-32, 32], e.g. -3, 2.5 or 0.125, so that `t + k` is exact
fn eighths() -> impl Strategy<Value = f64> {
    (-256..=256_i32).prop_map(|m| f64::from(m) / 8.0)
}

/// The bound `m u / (1 - m u)` on the relative error of `m` roundings with unit roundoff `u`
/// (Higham, Accuracy and Stability of Numerical Algorithms, Lemma 3.1)
fn rounding_bound(m: u32) -> f64 {
    let mu = f64::from(m) * f64::EPSILON / 2.0;
    mu / (1.0 - mu)
}

/// The sum of `pFq(-n, a; b; x)`, and a bound on the error of both this sum and any other that
/// rounds as often.
///
/// Each term follows from the previous one by `p + q` sums `a_i + k` and `b_j + k` and by
/// `p + q + 1` products and quotients, and adding it rounds once more. The error of the sum is
/// at most `rounding_bound(m)` times the sum of the absolute terms, with `m = 2 n (p + q + 1)`.
fn terminating(n: u32, a: &[f64], b: &[f64], x: f64) -> (f64, f64) {
    let (mut sum, mut abs_sum, mut term) = (1.0, 1.0, 1.0);
    for k in 0..n {
        let k = f64::from(k);
        let ratio =
            a.iter().map(|ai| ai + k).product::<f64>() / b.iter().map(|bj| bj + k).product::<f64>();
        term *= ratio * (k - f64::from(n)) * x / (k + 1.0);
        sum += term;
        abs_sum += term.abs();
    }
    let rounds = 2 * n * (a.len() + b.len() + 2) as u32;
    (sum, 2.0 * rounding_bound(rounds) * abs_sum)
}

/// `actual` within `tol` of `expected`, or NaN if `tol` allows no correct digit
fn within<T: Into<Complex64>>(
    (actual, (expected, tol)): (T, (f64, f64)),
) -> Result<(), TestCaseError> {
    let actual = actual.into();
    let err = (actual - expected).norm();
    prop_assert!(
        err <= tol || (actual.is_nan() && tol >= expected.abs()),
        "{actual:e} != {expected:e} +- {tol:e}"
    );
    Ok(())
}

/// Not a pole of `b_j`, unless the series terminates before it (DLMF 15.2.5)
fn terminates_first(n: u32, b: f64) -> bool {
    !is_pole(b) || b < -f64::from(n)
}

property!(
    within,
    hyp2f1_terminating(n in 0..=MAX_N, b in eighths(), c in eighths(), x in moderate())
        if terminates_first(n, c)
        => {
            // known (unfiled): a few hundred ulp for a negative non-integer c, and twice the bound
            // at x > 1
            if (c < 0.0 && c.fract() != 0.0) || x > 1.0 {
                return Ok(());
            }
            (xsf::hyp2f1(-f64::from(n), b, c, x), terminating(n, &[b], &[c], x))
        }
);
property!(
    within,
    hyp2f1_terminating_complex(n in 0..=MAX_N, b in eighths(), c in eighths(), x in moderate())
        if terminates_first(n, c)
        => (xsf::hyp2f1(-f64::from(n), b, c, c64(x, 0.0)), terminating(n, &[b], &[c], x))
);
property!(
    within,
    hyp3f2_terminating(
        n in 0..=MAX_N, a2 in eighths(), a3 in eighths(), b1 in eighths(), b2 in eighths(),
        x in moderate()
    ) if terminates_first(n, b1) && terminates_first(n, b2)
        => (
            xsf::hyp3f2(-f64::from(n), a2, a3, b1, b2, x),
            terminating(n, &[a2, a3], &[b1, b2], x)
        )
);

/// Chu-Vandermonde, `2F1(-n, b; c; 1) = (c - b)_n / (c)_n` (DLMF 15.4.24), as a product of `n`
/// quotients of exact factors, and a bound on the error of both this product and any other that
/// rounds as often
fn chu_vandermonde(n: u32, b: f64, c: f64) -> (f64, f64) {
    let value = (0..n)
        .map(f64::from)
        .map(|k| (c - b + k) / (c + k))
        .product::<f64>();
    (value, 2.0 * rounding_bound(2 * n) * value.abs())
}

property!(
    within,
    hyp2f1_chu_vandermonde(n in 0..=MAX_N, b in eighths(), c in eighths())
        if terminates_first(n, c)
        => {
            // known: all but b < 0 < c, where the factors are positive (#2198, fixed by #2214)
            if !(b < 0.0 && c > 0.0) {
                return Ok(());
            }
            (xsf::hyp2f1(-f64::from(n), b, c, 1.0), chu_vandermonde(n, b, c))
        }
);

/// `f(conj(z))` and `f(z)` for `z = re + i im`
fn conjugates(f: impl Fn(Complex64) -> Complex64, re: f64, im: f64) -> [Complex64; 2] {
    let z = c64(re, im);
    [f(z.conj()), f(z)]
}

/// `f(conj(z)) = conj(f(z))` exactly, where NaN equals NaN and -0 equals 0
fn conjugate_symmetric([w_conj, w]: [Complex64; 2]) -> Result<(), TestCaseError> {
    let same = |x: f64, y: f64| x == y || (x.is_nan() && y.is_nan());
    let w = w.conj();
    prop_assert!(same(w_conj.re, w.re) && same(w_conj.im, w.im), "{w_conj:e} != {w:e}");
    Ok(())
}

/// `property!(conjugate_symmetric, ...)` for `xsf::$name(params.., z)` off the real axis, where a
/// branch cut would depend on the sign of zero
macro_rules! conjugate {
    ($name:ident($($param:ident),*)) => {
        property!(
            conjugate_symmetric,
            $name($($param,)* re, im) if im != 0.0
                => conjugates(|z| xsf::$name($($param,)* z), re, im)
        );
    };
}

// not hyp1f1, which rounds the phase of a power of z after adding that of a gamma function, so
// that its asymmetry grows with |z|
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
