//! Symmetry in the parameters, e.g. B(a, b) = B(b, a)

use proptest::prelude::*;

use super::strategies::{MAX_OFFSET, arg, balanced, moderate};
use super::{config, int_distance, is_pole};

/// Within `tol` relative to a finite `expected`, or both the same inf, or both NaN
pub(super) fn close(actual: f64, expected: f64, tol: f64) -> bool {
    actual.to_bits() == expected.to_bits()
        || (actual.is_nan() && expected.is_nan())
        || (expected.is_finite() && (actual - expected).abs() <= tol * expected.abs())
}

/// A pole of Gamma, or within 1e-12 of a negative one, e.g. -3 + 1e-13
fn near_pole(x: f64) -> bool {
    is_pole(x) || (x.round() < 0.0 && int_distance(x) < 1e-12)
}

/// `hyp2f1(a, b, c, x) = hyp2f1(b, a, c, x)`, up to half the digits for ill-conditioned parameters
fn hyp2f1_symmetric(a: f64, b: f64, c: f64, x: f64) -> Result<(), TestCaseError> {
    // known: asymmetries near the poles of c
    prop_assume!(!near_pole(c));
    // known: a value or NaN depending on the order for x > 1 near an integer c - a - b (#1687)
    prop_assume!(!(x > 1.0 && int_distance(c - a - b) <= f64::EPSILON.sqrt()));
    let (ab, ba) = (xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x));
    prop_assert!(close(ab, ba, f64::EPSILON.sqrt()), "{ab:e} != {ba:e}");
    Ok(())
}

proptest! {
    #![proptest_config(config())]

    #[test]
    fn hyp2f1(a in moderate(), b in moderate(), c in moderate(), x in arg()) {
        hyp2f1_symmetric(a, b, c, x)?;
    }

    // a nearly integer c - a - b; known: 0 or NaN at x = -inf and inf at -f64::MAX (#2325), and a wrong value for c
    // next to a negative integer (#2326)
    #[test]
    fn hyp2f1_balanced([a, b, c] in balanced(), x in arg()) {
        if x > -f64::MAX && !(c < 0.0 && int_distance(c) <= MAX_OFFSET) {
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
