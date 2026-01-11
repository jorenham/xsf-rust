//! Translated from `scipy/special/_spfun_stats.py`

use core::f64::consts::PI;

/// Log of multivariate gamma, $\ln \Gamma_d(a)$, sometimes called the generalized gamma
///
/// Pure rust translation of [`scipy.special.multigammaln`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.multigammaln.html
///
/// # Definition
///
/// The multibariate gamma function of dimension $d$ for real $a$ is
///
/// $$
/// \Gamma_d(a) = \pi^{d (d-1) \over 4} \prod_{j=0}^{d-1} \Gamma\left(a - {j \over 2}\right)
/// $$
///
/// # Panics
/// - if $a \le (d - 1) / 2$
///
/// # See also
/// - [`gamma`](crate::gamma): Gamma function $\Gamma$
#[must_use]
#[inline]
pub fn multigammaln(a: f64, d: u32) -> f64 {
    let d_f64 = f64::from(d);
    let a_min = 0.5 * (d_f64 - 1.0);
    assert!(
        a > a_min,
        "condition a ({a}) > 0.5 * (d-1) ({a_min}) not met"
    );

    // res = (d * (d-1) * 0.25) * np.log(np.pi)
    // res += np.sum(loggam([(a - (j - 1.)/2) for j in range(1, d+1)]), axis=0)
    // return res
    d_f64 * (d_f64 - 1.0) * 0.25 * PI.ln()
        + (0..d)
            .map(|j| unsafe { crate::ffi::xsf::gammaln(a - f64::from(j) / 2.0) })
            .sum::<f64>()
}

#[cfg(test)]
mod tests {
    //! Translated from `scipy.special.tests.test_spfun_stats.TestMultiGammaLn` at
    //! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_spfun_stats.py>

    use core::f64::consts::PI;

    use crate::np_assert_allclose;

    /// TestMultiGammaLn.test1
    #[test]
    #[allow(clippy::float_cmp)]
    fn test_multigammaln_1() {
        // A test of the identity:
        //     Gamma_1(a) = Gamma(a)

        // np.random.seed(1234)
        // a = np.abs(np.random.randn())
        let a = 0.471_435_163_732_493_06;
        // assert_array_equal(multigammaln(a, 1), gammaln(a))
        assert_eq!(crate::multigammaln(a, 1), crate::gammaln(a));
    }

    /// TestMultiGammaLn.test2
    #[test]
    fn test_multigammaln_2() {
        // A test of the identity
        //     Gamma_2(a) = sqrt(pi) * Gamma(a) * Gamma(a - 0.5)

        // a = np.array([2.5, 10.0])
        let a = [2.5, 10.0];
        // result = multigammaln(a, 2)
        let result = a.map(|a| crate::multigammaln(a, 2));
        // expected = np.log(np.sqrt(np.pi)) + gammaln(a) + gammaln(a - 0.5)
        let expected = a.map(|a| PI.sqrt().ln() + crate::gammaln(a) + crate::gammaln(a - 0.5));
        // assert_allclose(result, expected, atol=1.5e-7, rtol=0)
        np_assert_allclose!(result, expected, atol = 1.5e-7, rtol = 0.0);
    }
}
