/// Inverse of the regularized incomplete Beta function, $y = \I_x(a, b)$
///
/// Computes $x \in \[0, 1\]$ such that [^DLMF]:
///
/// $$
/// y = \I_x(a, b) = {1 \over \B(a, b)} \int_0^x t^{a-1} (1-t)^{b-1} \dd t ,
/// $$
///
/// for $0 \le y \le 1$ and $a,b \ge 0$, with $\B(a, b)$ the Beta function.
///
/// This function is the quantile function (inverse CDF) of the Beta distribution.
///
/// # Notes
/// This functions wraps the `incbi` Cephes routine [^CEPHES], making it less accurate than
/// [`scipy.special.betaincinv`][scipy], which wraps the Boost `ibeta_inv` routine [^BOOST].
/// Especially for small `y` and large `a - b` the accuracy may be poor.
///
/// # See also
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function $\I_x(a, b)$
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(a, b)$
/// - [`scipy.special.betaincinv`][scipy]: Corresponding function in SciPy
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincinv.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17>.
/// [^CEPHES]: Cephes Math Library Release 2.4. Translated into C++ by SciPy developers in 2024.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(
    alias = "inc_beta_inv",
    alias = "beta_inc_inv",
    alias = "ibeta_inv",
    alias = "incbi"
)]
#[must_use]
#[inline]
pub fn betaincinv(a: f64, b: f64, y: f64) -> f64 {
    if a.is_nan() || b.is_nan() || y.is_nan() || a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&y) {
        f64::NAN
    } else {
        unsafe { crate::ffi::xsf::incbi(a, b, y) }
    }
}

#[cfg(test)]
mod tests {
    //! based on scipy.special.tests.test_basic.TestBetaInc

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_betaincinv_a1_b1() {
        // betaincinv(1, 1, x) is x.
        let x = [0.0, 0.25, 1.0];
        assert_eq!(x.map(|xv| crate::betaincinv(1.0, 1.0, xv)), x);
    }

    #[test]
    fn test_betaincinv_nontrivial() {
        for &(a, b, x, p) in &[
            (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
            (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
            // scipy/scipy#3761:
            (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
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
            let x1 = crate::betaincinv(a, b, p);
            crate::np_assert_allclose!(&[x1], &[x], rtol = 5e-13);
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
                14.208_308_325_339_239,
                14.208_308_325_339_239,
                7.703_145_458_496_392e-307,
                8.566_004_561_846_704e-023,
            ),
            (14.0, 14.5, 1e-280, 2.934_391_500_664_242_4e-21),
            (3.5, 15.0, 4e-95, 1.329_075_142_928_922_7e-28),
            (10.0, 1.25, 2e-234, 3.982_659_092_143_654e-24),
            // NOTE: This results in a relative error of 3.806e-11
            // (4.0, 99997.0, 5e-88, 3.309800566862242e-27),
        ] {
            let x = crate::betaincinv(a, b, y);
            crate::np_assert_allclose!(&[x], &[ref_], rtol = 1e-14);
        }
    }

    #[test]
    fn test_betaincinv_gh21426() {
        // Test for scipy/scipy#21426: betaincinv must not return NaN
        let a = 5.0;
        let x = 0.5;
        let result = crate::betaincinv(a, a, x);
        crate::np_assert_allclose!(&[result], &[0.5], rtol = 10.0 * f64::EPSILON);
    }
}
