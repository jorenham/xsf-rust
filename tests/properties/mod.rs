//! Property tests with random inputs from `proptest`.
//!
//! The seed is fixed unless `PROPTEST_RNG_SEED` is set, and the number of cases is
//! `PROPTEST_CASES` (default 256). Hangs are only caught with `PROPTEST_TIMEOUT=<ms>`.

#[macro_use]
mod macros;
mod strategies;

mod continuity;
mod exact;
mod no_panic;
mod ranges;
mod symmetry;

use proptest::test_runner::{Config, RngSeed};

/// The default config, but with a fixed seed, and as many rejects as cases
fn config() -> Config {
    let default = Config::default();
    let rng_seed = match default.rng_seed {
        RngSeed::Random => RngSeed::Fixed(0),
        seed @ RngSeed::Fixed(_) => seed,
    };
    Config {
        rng_seed,
        max_global_rejects: default.cases,
        ..default
    }
}

/// A pole of Gamma, e.g. -3 or -0
fn is_pole(x: f64) -> bool {
    x <= 0.0 && x.fract() == 0.0
}

/// A polynomial parameter of degree 0 or 1, i.e. 0 or -1
fn is_linear(x: f64) -> bool {
    (-1.0..=0.0).contains(&x) && x.fract() == 0.0
}

/// The distance to the nearest integer, e.g. 0.25 for 2.75
fn int_distance(x: f64) -> f64 {
    (x - x.round()).abs()
}

/// Within `sqrt(eps)` of an integer, e.g. 3 + 1e-12, where terms with opposite poles at the
/// integers cancel more than half the digits
fn almost_int(x: f64) -> bool {
    int_distance(x) <= f64::EPSILON.sqrt()
}

/// The distance to the nearest pole of Gamma, e.g. 1e-13 for -3 + 1e-13, or 0.75 for 0.75
fn pole_distance(x: f64) -> f64 {
    (x - x.round().min(0.0)).abs()
}

/// The largest `pole_distance` of the inputs next to a pole, e.g. -3 + 1e-13 or -1e-15, and the
/// largest `int_distance` of those next to an integer
const NEAR_POLE: f64 = 1e-12;
