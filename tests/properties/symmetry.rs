//! Symmetry in the parameters, e.g. B(a, b) = B(b, a)

// exact parameters, e.g. c = a
#![allow(clippy::float_cmp)]

use proptest::prelude::*;

use super::strategies::{arg, balanced, moderate};
use super::{NEAR_POLE, config, int_distance, is_pole, pole_distance};

/// Within `tol` relative to a finite `expected`, or both the same inf, or both NaN
pub(super) fn close(actual: f64, expected: f64, tol: f64) -> bool {
    actual.to_bits() == expected.to_bits()
        || (actual.is_nan() && expected.is_nan())
        || (expected.is_finite() && (actual - expected).abs() <= tol * expected.abs())
}

/// A pole of Gamma, or within `NEAR_POLE` of a negative one, e.g. -3 + 1e-13
fn at_pole(x: f64) -> bool {
    is_pole(x) || (x.round() < 0.0 && int_distance(x) < NEAR_POLE)
}

/// `hyp2f1(a, b, c, x) = hyp2f1(b, a, c, x)`, up to half the digits for ill-conditioned parameters
fn hyp2f1_symmetric(a: f64, b: f64, c: f64, x: f64) -> Result<(), TestCaseError> {
    // known: asymmetries near the poles of c
    prop_assume!(!at_pole(c));
    let (ab, ba) = (xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x));
    prop_assert!(close(ab, ba, f64::EPSILON.sqrt()), "{ab:e} != {ba:e}");
    Ok(())
}

/// The largest distance to the integers of the parameters of the known asymmetries #2326 and
/// #2333, e.g. 1.7e-7
const NEAR_DEGENERATE: f64 = 1e-6;

/// The known asymmetries of `hyp2f1` for a nearly integer `c - a - b`
fn known_hyp2f1_balanced(a: f64, b: f64, c: f64, x: f64) -> bool {
    let tiny = |t: f64| t != 0.0 && t.abs() <= NEAR_DEGENERATE;
    let near = |t: f64| int_distance(t) <= NEAR_DEGENERATE;
    // wrong for x > 0.9 if p is next to a pole, and c < 0, the other and c minus it next to
    // integers, and c - a - b to 0 or -1 (#2326)
    let n = (c - a - b).round();
    let integers = |p: f64, q: f64| {
        pole_distance(p) <= NEAR_DEGENERATE && near(q) && near(c - q) && (c - q).round() < 0.0
    };
    // wrong for x > 0.9 or x < -100 (seen from -288) if c - a - b is next to -1, for a tiny a or b
    // (#2333), or one next to a negative pole (#2343)
    let shifted = |p: f64| p.fract() != 0.0 && pole_distance(p) <= NEAR_DEGENERATE;
    // NaN in one order at x = -f64::MAX, for a tiny a or b (#2381)
    (x == -f64::MAX && (tiny(a) || tiny(b)))
        // NaN or infinite in one order for x > 1, where c - a - b is an integer in only one order
        // of a and b (#2341)
        || (x > 1.0 && ((c - a - b).fract() == 0.0) != ((c - b - a).fract() == 0.0))
        || ((-1.0..=0.0).contains(&n) && c < 0.0 && x > 0.9 && (integers(a, b) || integers(b, a)))
        || (n == -1.0 && near(c - a - b) && (x > 0.9 || x < -100.0) && (shifted(a) || shifted(b)))
        // 1 in one order for c = a or c = b in (-1e-13, 0) and x > -1 (#2342)
        || ((a == c || b == c) && (-1e-13..0.0).contains(&c) && x > -1.0)
}

proptest! {
    #![proptest_config(config())]

    #[test]
    fn hyp2f1(a in moderate(), b in moderate(), c in moderate(), x in arg()) {
        hyp2f1_symmetric(a, b, c, x)?;
    }

    // a nearly integer c - a - b
    #[test]
    fn hyp2f1_balanced([a, b, c] in balanced(), x in arg()) {
        if !known_hyp2f1_balanced(a, b, c, x) {
            hyp2f1_symmetric(a, b, c, x)?;
        }
    }

    #[test]
    fn beta(a in arg(), b in arg()) {
        // e.g. beta(-0, -29) = -inf but beta(-29, -0) = inf
        prop_assume!(!(is_pole(a) && is_pole(b)));
        let (ab, ba) = (xsf::beta(a, b), xsf::beta(b, a));
        prop_assert!(close(ab, ba, 0.0), "{ab:e} != {ba:e}");
    }
}
