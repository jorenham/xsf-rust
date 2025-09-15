use crate::bindings;

#[inline(always)]
fn incbi(a: f64, b: f64, y: f64) -> f64 {
    unsafe { bindings::incbi(a, b, y) }
}

/// Inverse of the regularized incomplete beta function
///
///
/// Note: The Cephes backend is less accurate than the `scipy.special` Boost backend, especially
/// for small `y` and large `a - b` differences.
///
/// See also: [`betainc`]
pub fn betaincinv(a: f64, b: f64, y: f64) -> f64 {
    if a.is_nan() || b.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&y) {
        return f64::NAN;
    }

    // Cephes does not require exception handling like Boost does.
    incbi(a, b, y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::np_assert_allclose;

    // based on scipy.special.tests.test_basic.TestBetaInc

    #[test]
    fn test_betaincinv_a1_b1() {
        // betaincinv(1, 1, x) is x.
        let x = [0.0, 0.25, 1.0];
        assert_eq!(x.map(|xv| betaincinv(1.0, 1.0, xv)), x);
    }

    #[test]
    fn test_betaincinv_nontrivial() {
        for &(a, b, x, p) in &[
            (2.0, 4.0, 0.3138101704556974, 0.5),
            (0.0342, 171.0, 1e-10, 0.5526991690180709),
            // scipy/scipy#3761:
            (0.0342, 171.0, 8.42313169354797e-21, 0.25),
            // scipy/scipy#4244 (NOTE: this results in a relative error of 1.872e-6)
            // (
            //     0.0002742794749792665,
            //     289206.03125,
            //     1.639984034231756e-56,
            //     0.9688708782196045,
            // ),
            // scipy/scipy#12796 (NOTE: this results in a relative error of 1.787e2)
            // (4.0, 99997.0, 0.0001947841578892121, 0.999995),
        ] {
            let x1 = betaincinv(a, b, p);
            np_assert_allclose(&[x1], &[x], 5e-13, 0.0);
        }
    }

    #[test]
    fn test_betaincinv_tiny_y() {
        // Test with extremely small y values.  This test includes
        // a regression test for an issue in the boost code;
        // see https://github.com/boostorg/math/issues/961
        //
        // The reference values were computed with mpmath. For example,
        //
        //   from mpmath import mp
        //   mp.dps = 1000
        //   a = 14.208308325339239
        //   p = 7.703145458496392e-307
        //   x = mp.findroot(lambda t: mp.betainc(a, a, 0, t,
        //                                        regularized=True) - p,
        //                   x0=8.566e-23)
        //   print(float(x))

        for &(a, b, y, ref_) in &[
            (
                14.208308325339239,
                14.208308325339239,
                7.703145458496392e-307,
                8.566004561846704e-23,
            ),
            (14.0, 14.5, 1e-280, 2.9343915006642424e-21),
            (3.5, 15.0, 4e-95, 1.3290751429289227e-28),
            (10.0, 1.25, 2e-234, 3.982659092143654e-24),
            // NOTE: This results in a relative error of 3.806e-11
            // (4.0, 99997.0, 5e-88, 3.309800566862242e-27),
        ] {
            let x = betaincinv(a, b, y);
            np_assert_allclose(&[x], &[ref_], 1e-14, 0.0);
        }
    }

    #[test]
    fn test_betaincinv_gh21426() {
        // Test for scipy/scipy#21426: betaincinv must not return NaN
        let a = 5.0;
        let x = 0.5;
        let result = betaincinv(a, a, x);
        np_assert_allclose(&[result], &[0.5], 10.0 * f64::EPSILON, 0.0);
    }
}
