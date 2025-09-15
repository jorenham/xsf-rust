use crate::bindings;

#[inline(always)]
fn incbet(a: f64, b: f64, x: f64) -> f64 {
    unsafe { bindings::incbet(a, b, x) }
}

/// Regularized incomplete beta function
///
/// Note: This is a custom wrapper around `cephes::incbet`, analogous to how `scipy.special.betainc`
/// wraps `boost::math::ibeta` in `scipy/special/boost_special_functions.h`. The Cephes
/// implementation tends to be less accurate than Boost's, especially for small `x`.
pub fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if a.is_nan() || b.is_nan() || x.is_nan() {
        return f64::NAN;
    }

    if a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&x) {
        return f64::NAN;
    }

    // In limiting cases, SciPy treats `betainc` as a two parameter family
    // of functions of a single variable `x`, rather than as a function of
    // three variables `a`, `b`, `x`. The limit `(a, b) -> (a0, b0)` of
    // `betainc(a, b, x)` is treated as the pointwise limit in `x`.

    if (a == 0.0 && b == 0.0) || (a.is_infinite() && b.is_infinite()) {
        // In the limit (a, b) -> (0+, 0+), the Beta distribution converges
        // to a Bernoulli(p) distribution, where p depends on the path in
        // which (a, b) approaches (0+, 0+).
        // e.g. if a = t*b then the limiting distribution will be
        // Bernoulli(t / (t + 1)). The a = 0, b = 0 case is thus indeterminate.
        // A similar statement can be made for the limit (a, b) -> (inf, inf).
        return f64::NAN;
    }

    if a == 0.0 || b.is_infinite() {
        // Distribution in the limit a -> 0+, b > 0 is a point distribution
        // at x = 0. The same is true in the limit b -> inf for fixed a.
        return if x > 0.0 { 1.0 } else { 0.0 };
    }

    if a.is_infinite() || b == 0.0 {
        // Distribution in the limit b -> 0+, a > 0 is a point distribution
        // at x = 1. The same is true in the limit a -> inf for fixed b.
        return if x < 1.0 { 0.0 } else { 1.0 };
    }

    // Cephes does not require exception handling like Boost does.
    incbet(a, b, x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{np_assert_allclose, np_assert_equal};

    // based on scipy.special.tests.test_basic.TestBetaInc

    #[test]
    fn test_betainc_a1_b1() {
        // betainc(1, 1, x) is x.
        let x = [0.0, 0.25, 1.0];
        assert_eq!(x.map(|xv| betainc(1.0, 1.0, xv)), x);
    }

    #[test]
    fn test_betainc_and_inverses() {
        for &(a, b, x, p) in &[
            (2.0, 4.0, 0.3138101704556974, 0.5),
            (0.0342, 171.0, 1e-10, 0.5526991690180709),
            // scipy/scipy#3761:
            (0.0342, 171.0, 8.42313169354797e-21, 0.25),
            // scipy/scipy#4244 (NOTE: this results in a relative error of 5.134e-10)
            // (
            //     0.0002742794749792665,
            //     289206.03125,
            //     1.639984034231756e-56,
            //     0.9688708782196045,
            // ),
            // scipy/scipy#12796:
            (4.0, 99997.0, 0.0001947841578892121, 0.999995),
        ] {
            let p1 = betainc(a, b, x);
            // NOTE: The original `rtol = 1e-15` is too strict for Cephes
            np_assert_allclose(&[p1], &[p], 1e-14, 0.0);
        }
    }

    #[test]
    fn test_betainc_edge_cases() {
        for &((a, b, x), expected) in &[
            ((0.0, 0.0, 0.0), f64::NAN),
            ((0.0, 0.0, 0.5), f64::NAN),
            ((0.0, 0.0, 1.0), f64::NAN),
            ((f64::INFINITY, f64::INFINITY, 0.0), f64::NAN),
            ((f64::INFINITY, f64::INFINITY, 0.5), f64::NAN),
            ((f64::INFINITY, f64::INFINITY, 1.0), f64::NAN),
            ((0.0, 1.0, 0.0), 0.0),
            ((0.0, 1.0, 0.5), 1.0),
            ((0.0, 1.0, 1.0), 1.0),
            ((1.0, 0.0, 0.0), 0.0),
            ((1.0, 0.0, 0.5), 0.0),
            ((1.0, 0.0, 1.0), 1.0),
            ((0.0, f64::INFINITY, 0.0), 0.0),
            ((0.0, f64::INFINITY, 0.5), 1.0),
            ((0.0, f64::INFINITY, 1.0), 1.0),
            ((f64::INFINITY, 0.0, 0.0), 0.0),
            ((f64::INFINITY, 0.0, 0.5), 0.0),
            ((f64::INFINITY, 0.0, 1.0), 1.0),
            ((1.0, f64::INFINITY, 0.0), 0.0),
            ((1.0, f64::INFINITY, 0.5), 1.0),
            ((1.0, f64::INFINITY, 1.0), 1.0),
            ((f64::INFINITY, 1.0, 0.0), 0.0),
            ((f64::INFINITY, 1.0, 0.5), 0.0),
            ((f64::INFINITY, 1.0, 1.0), 1.0),
        ] {
            let observed = betainc(a, b, x);
            np_assert_equal(&[observed], &[expected]);
        }
    }

    #[test]
    fn test_betainc_gh22682() {
        // scipy/scipy#22682: betainc returned incorrect results for tiny
        // single precision inputs. test that this is resolved
        let rtol = 1e-15;
        for &(a, b, x, reference) in &[
            (1e-20, 1e-21, 0.5, 0.0909090909090909),
            (1e-15, 1e-16, 0.5, 0.09090909090909091),
        ] {
            let res = betainc(a, b, x);
            np_assert_allclose(&[res], &[reference], rtol, 0.0);
        }
    }
}
