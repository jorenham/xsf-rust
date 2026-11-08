//! Symmetry in the parameters, e.g. B(a, b) = B(b, a)

use proptest::prelude::*;

use super::config;
use super::strategies::{arg, moderate};

/// Within `tol` relative to a finite `expected`, or both the same inf, or both NaN
fn close(actual: f64, expected: f64, tol: f64) -> bool {
    actual.to_bits() == expected.to_bits()
        || (actual.is_nan() && expected.is_nan())
        || (expected.is_finite() && (actual - expected).abs() <= tol * expected.abs())
}

/// A pole of Gamma, e.g. -3 or -0
fn is_pole(x: f64) -> bool {
    x <= 0.0 && x.fract() == 0.0
}

/// A pole of Gamma, or within 1e-12 of a negative one, e.g. -3 + 1e-13
fn near_pole(x: f64) -> bool {
    is_pole(x) || (x.round() < 0.0 && (x - x.round()).abs() < 1e-12)
}

proptest! {
    #![proptest_config(config())]

    // known: asymmetries for x < -1 or x > 0.9, e.g. hyp2f1(156.000001, 0.5, 43, -96) is NaN
    #[test]
    fn hyp2f1(a in moderate(), b in moderate(), c in moderate(), x in -1.0..0.9_f64) {
        // known: also near the poles of c
        prop_assume!(!near_pole(c));
        let (ab, ba) = (xsf::hyp2f1(a, b, c, x), xsf::hyp2f1(b, a, c, x));
        // half the digits, for ill-conditioned parameters
        prop_assert!(close(ab, ba, f64::EPSILON.sqrt()), "{ab:e} != {ba:e}");
    }

    #[test]
    fn beta(a in arg(), b in arg()) {
        // e.g. beta(-0, -29) = -inf but beta(-29, -0) = inf
        prop_assume!(!(is_pole(a) && is_pole(b)));
        let (ab, ba) = (xsf::beta(a, b), xsf::beta(b, a));
        prop_assert!(close(ab, ba, 0.0), "{ab:e} != {ba:e}");
    }
}
