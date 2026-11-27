//! Symmetry in the parameters, e.g. B(a, b) = B(b, a)

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
    // usually NaN for x > 1, but not where a c - a - b within about 1e-13 of an integer counts as
    // one (see its docs), which can depend on the order of a and b
    prop_assume!(!(x > 1.0 && int_distance(c - a - b) <= NEAR_POLE));
    let (ab, ba) = (xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x));
    prop_assert!(close(ab, ba, f64::EPSILON.sqrt()), "{ab:e} != {ba:e}");
    Ok(())
}

/// The largest distance to the integers of the parameters of the known asymmetries #2326 and
/// #2333, e.g. 1.7e-7
const NEAR_DEGENERATE: f64 = 1e-6;

/// The known asymmetries of `hyp2f1` for a nearly integer `c - a - b`
fn known_hyp2f1_balanced(a: f64, b: f64, c: f64, x: f64) -> bool {
    // NaN for c < 0 at a tiny x > 0, where a term ratio of the power series underflows (#2332)
    let ratio = (a * b * x / c).abs();
    let underflow = c < 0.0 && x > 0.0 && (x < f64::MIN_POSITIVE || ratio < f64::MIN_POSITIVE);
    // wrong for x > 0.9 or x < -1 if one of a and b is next to a pole, and c minus the other next
    // to a negative integer, the same or one less (#2326, #2333)
    let near = |t: f64| int_distance(t) <= NEAR_DEGENERATE;
    let degenerate = |p: f64, q: f64| {
        pole_distance(p) <= NEAR_DEGENERATE && near(c - q) && (c - q).round() < 0.0
    };
    let n = (c - a - b).round();
    let transformed = x > 0.9 || x < -1.0;
    // 0 or NaN at x = -inf, and inf at -f64::MAX (#2325)
    x <= -f64::MAX
        || underflow
        || (transformed && (-1.0..=0.0).contains(&n) && (degenerate(a, b) || degenerate(b, a)))
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
