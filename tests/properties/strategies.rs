use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, PI};
use std::ops::{Range, RangeInclusive};

use proptest::num::f64::ANY;
use proptest::prelude::*;

/// Log-uniform magnitudes 2^e with a random sign, e.g. -3.7e-5 or 1.2e9
fn log_uniform(exponents: Range<f64>) -> impl Strategy<Value = f64> {
    (any::<bool>(), exponents).prop_map(|(neg, e)| if neg { -e.exp2() } else { e.exp2() })
}

/// Integers in `range`, some of them perturbed, e.g. -3 or 2 + 1e-9
fn near_int(range: Range<i32>) -> impl Strategy<Value = f64> {
    let delta = prop_oneof![Just(0.0), -1e-6..1e-6, -1e-14..1e-14];
    (range, delta).prop_map(|(n, d)| f64::from(n) + d)
}

/// The largest `|offset()|`
pub(super) const MAX_OFFSET: f64 = 1.0 / 256.0;

/// Small offsets from `-MAX_OFFSET` to `MAX_OFFSET`, down to `2^-60`, e.g. 3e-9 or -1e-15
fn offset() -> impl Strategy<Value = f64> {
    log_uniform(-60.0..MAX_OFFSET.log2())
}

/// Within `MAX_OFFSET` of `-n` for `n` in `poles`, e.g. -3 + 1e-9 or 4e-17
fn near_pole(poles: RangeInclusive<i32>) -> impl Strategy<Value = f64> {
    (poles, offset()).prop_map(|(n, d)| -f64::from(n) + d)
}

/// Subnormal or near-overflow magnitudes with a random sign, e.g. -5e-324 or 2^900
fn extreme() -> impl Strategy<Value = f64> {
    let subnormal = (any::<bool>(), 1..1_u64 << (f64::MANTISSA_DIGITS - 1))
        .prop_map(|(neg, m)| f64::from_bits(m | u64::from(neg) << 63));
    let max_exp = f64::from(f64::MAX_EXP);
    prop_oneof![subnormal, log_uniform(max_exp / 2.0..max_exp)]
}

/// The radii of the switches between the transformations of `hyp2f1`, 1/2, 1 and 2, perturbed
/// by up to `2^-4` down to an ulp, e.g. 1 - 2^-30 or -0.5000001
fn near_threshold() -> impl Strategy<Value = f64> {
    let radii = prop::sample::select(vec![0.5, 1.0, 2.0, -0.5, -1.0, -2.0]);
    (radii, log_uniform(-f64::from(f64::MANTISSA_DIGITS)..-4.0)).prop_map(|(r, d)| r * (1.0 + d))
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
    prop_oneof![
        4 => ANY,
        4 => log_uniform(-64.0..64.0),
        4 => near_int(-200..200),
        4 => special(),
        1 => near_pole(0..=200),
        1 => extreme(),
        1 => near_threshold(),
    ]
}

/// The imaginary part of `z`, with extra `+0` and `-0`, the sides of the branch cuts
pub(super) fn imag() -> impl Strategy<Value = f64> {
    prop_oneof![8 => arg(), 1 => Just(0.0), 1 => Just(-0.0)]
}

/// `(re, im)` next to the circles `|z| = 1/2, 1, 2` and the rays of angle 0, pi/3, pi/2 and pi,
/// e.g. (0.5, 0.8660254) for `exp(i pi/3)`
pub(super) fn near_threshold_z() -> impl Strategy<Value = (f64, f64)> {
    let angle = prop_oneof![
        prop::sample::select(vec![0.0, FRAC_PI_3, FRAC_PI_2, PI]),
        log_uniform(-64.0..1.5),
    ];
    (near_threshold(), angle).prop_map(|(r, t)| (r * t.cos(), r * t.sin()))
}

/// `[a, b, c]` with `c - a - b` an integer in [-4, 4], or within `MAX_OFFSET` of it, e.g.
/// [0.3, 1.2, 2.5 + 1e-9]
pub(super) fn balanced() -> impl Strategy<Value = [f64; 3]> {
    let delta = prop_oneof![Just(0.0), offset()];
    (moderate(), moderate(), -4..=4_i32, delta)
        .prop_map(|(a, b, n, d)| [a, b, a + b + f64::from(n) + d])
}

/// Magnitudes up to 2^8, or near an integer
pub(super) fn moderate() -> impl Strategy<Value = f64> {
    prop_oneof![
        4 => log_uniform(-64.0..8.0),
        4 => near_int(-200..200),
        2 => near_pole(0..=200),
    ]
}

/// Magnitudes up to 2^4, or near an integer
pub(super) fn small() -> impl Strategy<Value = f64> {
    prop_oneof![4 => log_uniform(-64.0..4.0), 4 => near_int(-16..16), 2 => near_pole(0..=16)]
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
