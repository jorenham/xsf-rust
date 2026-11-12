//! Translated from `scipy.special.tests.test_spfun_stats.TestMultiGammaLn` at
//! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_spfun_stats.py>

use core::f64::consts::PI;
use xsref::np_assert_allclose;

/// TestMultiGammaLn.test1
#[test]
fn test_multigammaln_1() {
    // A test of the identity:
    //     Gamma_1(a) = Gamma(a)

    // np.random.seed(1234)
    // a = np.abs(np.random.randn())
    let a = 0.471_435_163_732_493_06;
    // assert_array_equal(multigammaln(a, 1), gammaln(a))
    assert_eq!(xsf::multigammaln(a, 1), xsf::gammaln(a));
}

/// TestMultiGammaLn.test2
#[test]
fn test_multigammaln_2() {
    // A test of the identity
    //     Gamma_2(a) = sqrt(pi) * Gamma(a) * Gamma(a - 0.5)

    // a = np.array([2.5, 10.0])
    let a = [2.5, 10.0];
    // result = multigammaln(a, 2)
    let result = a.map(|a| xsf::multigammaln(a, 2));
    // expected = np.log(np.sqrt(np.pi)) + gammaln(a) + gammaln(a - 0.5)
    let expected = a.map(|a| PI.sqrt().ln() + xsf::gammaln(a) + xsf::gammaln(a - 0.5));
    // assert_allclose(result, expected, atol=1.5e-7, rtol=0)
    np_assert_allclose!(result, expected, atol = 1.5e-7, rtol = 0.0);
}
