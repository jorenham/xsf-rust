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

/// The distance to the nearest integer, e.g. 0.25 for 2.75
fn int_distance(x: f64) -> f64 {
    (x - x.round()).abs()
}
