//! Distribution functions in [0, 1], and not NaN, for valid parameters

use proptest::prelude::*;

use super::strategies::{nonnegative, positive, real, unit};

/// Known: `betainc` and `fdtr` can be NaN if both a and b are above this
const HUGE: f64 = 1.0 / f64::EPSILON;

/// Known: `gammainc` and `betainc` can be up to ~300 ulp (7e-14) above 1 for a < 1
const KNOWN_ERROR: f64 = 1e-12;

fn in_unit_interval(f: f64) -> Result<(), TestCaseError> {
    prop_assert!((0.0..=1.0).contains(&f), "{f:e}");
    Ok(())
}

/// In [0, 1 + `KNOWN_ERROR`]
fn nearly_in_unit_interval(f: f64) -> Result<(), TestCaseError> {
    prop_assert!((0.0..=1.0 + KNOWN_ERROR).contains(&f), "{f:e}");
    Ok(())
}

property!(nearly_in_unit_interval, gammainc(a in positive(), x in nonnegative()));
property!(in_unit_interval, gammaincc(a in positive(), x in nonnegative()));
property!(
    nearly_in_unit_interval,
    betainc(a in positive(), b in positive(), x in unit()) if a.min(b) <= HUGE
);

property!(
    nearly_in_unit_interval,
    fdtr(a in positive(), b in positive(), x in nonnegative()) if a.min(b) <= HUGE
);
property!(in_unit_interval, stdtr(nu in positive(), x in real()));
property!(nearly_in_unit_interval, chdtr(nu in positive(), x in nonnegative()));
property!(nearly_in_unit_interval, gdtr(a in positive(), b in positive(), x in nonnegative()));
// pdtr(inf, inf) is inf / inf
property!(
    in_unit_interval,
    pdtr(k in nonnegative(), m in nonnegative()) if k.is_finite() || m.is_finite()
);
property!(
    in_unit_interval,
    bdtr(k in 0..1000_i32, n in 0..1000_i32, p in unit()) => xsf::bdtr(f64::from(k.min(n)), n, p)
);
property!(in_unit_interval, nbdtr(k in 0..1000_i32, n in 1..1000_i32, p in unit()));

property!(in_unit_interval, kolmogorov(x in nonnegative()));
property!(in_unit_interval, smirnov(n in 1..i32::MAX, x in unit()));
property!(in_unit_interval, ndtr(x in real()));
