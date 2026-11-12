//! Translated from `scipy/special/_spfun_stats.py`

use crate::utils::LN_PI;

/// Log of multivariate gamma, $\ln \Gamma_d(a)$, sometimes called the generalized gamma
///
/// Pure rust translation of [`scipy.special.multigammaln`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.multigammaln.html
///
/// # Definition
///
/// The multivariate gamma function of dimension $d$ for real $a$ is
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
    d_f64 * (d_f64 - 1.0) * 0.25 * LN_PI
        + (0..d)
            .map(|j| crate::gammaln(a - f64::from(j) / 2.0))
            .sum::<f64>()
}
