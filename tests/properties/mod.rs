//! Property tests with random inputs from `proptest`.
//!
//! The seed is fixed unless `PROPTEST_RNG_SEED` is set, and the number of cases is
//! `PROPTEST_CASES` (default 256). Hangs are only caught with `PROPTEST_TIMEOUT=<ms>`.

#[macro_use]
mod macros;
mod strategies;

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
