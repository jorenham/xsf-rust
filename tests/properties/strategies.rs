use std::ops::Range;

use proptest::num::f64::ANY;
use proptest::prelude::*;

/// Log-uniform magnitudes 2^e with a random sign, e.g. -3.7e-5 or 1.2e9
fn log_uniform(exponents: Range<f64>) -> impl Strategy<Value = f64> {
    (any::<bool>(), exponents).prop_map(|(neg, e)| if neg { -e.exp2() } else { e.exp2() })
}

/// Integers in [-200, 200), some of them perturbed, e.g. -3 or 2 + 1e-9
fn near_int() -> impl Strategy<Value = f64> {
    let delta = prop_oneof![Just(0.0), -1e-6..1e-6, -1e-14..1e-14];
    (-200_i32..200, delta).prop_map(|(n, d)| f64::from(n) + d)
}

/// Boundaries and branch points, e.g. -0.0, 1 - 2^-30 or f64::MAX
fn special() -> impl Strategy<Value = f64> {
    let points = vec![
        0.0,
        -0.0,
        0.5,
        -0.5,
        1.0,
        -1.0,
        2.0,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ];
    prop_oneof![
        prop::sample::select(points),
        (1..53_i32).prop_map(|k| 1.0 - (-f64::from(k)).exp2()),
        (1..53_i32).prop_map(|k| -1.0 - (-f64::from(k)).exp2()),
    ]
}

/// Any `f64`, but mostly the interesting ones
pub(super) fn arg() -> impl Strategy<Value = f64> {
    prop_oneof![ANY, log_uniform(-64.0..64.0), near_int(), special()]
}

/// Magnitudes up to 2^8, or near an integer
pub(super) fn moderate() -> impl Strategy<Value = f64> {
    prop_oneof![log_uniform(-64.0..8.0), near_int()]
}

/// Not NaN
pub(super) fn real() -> impl Strategy<Value = f64> {
    arg().prop_filter("not NaN", |x| !x.is_nan())
}

/// In [0, inf]
pub(super) fn nonnegative() -> impl Strategy<Value = f64> {
    real().prop_map(f64::abs)
}

/// In (0, inf)
pub(super) fn positive() -> impl Strategy<Value = f64> {
    nonnegative().prop_filter("finite and positive", |x| x.is_finite() && *x > 0.0)
}

/// In [0, 1], e.g. 0.3, 2e-15 or 1
pub(super) fn unit() -> impl Strategy<Value = f64> {
    let ends = prop::sample::select(vec![0.0, 1.0]);
    prop_oneof![0.0..=1.0, log_uniform(-64.0..0.0).prop_map(f64::abs), ends]
}
