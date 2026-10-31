//! Translated from <https://github.com/scipy/scipy/blob/38ecfef/scipy/special/_basic.py#L1539-L1579>

use crate::xsf::cephes::zeta_euler_maclaurin;

/// Polygamma function ψ<sup>(n)</sup>(x)
///
/// The n-th derivative of the [`digamma`](crate::digamma) function.
///
/// Corresponds to [`scipy.special.polygamma`][polygamma] in scipy
///
/// [polygamma]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.polygamma.html
///
/// # See also
/// - [`digamma`](crate::digamma): Digamma function (0th derivative)
#[must_use]
#[inline]
pub fn polygamma(n: u32, x: f64) -> f64 {
    if n == 0 {
        crate::digamma(x)
    } else {
        let n1p = f64::from(n) + 1.0;
        let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
        let (g, z) = (crate::gamma(n1p), crate::zeta(n1p, x));
        // Unlike in SciPy, if `n!` overflows or `zeta(n + 1, x)` underflows, this uses
        // `n! zeta(n + 1, x) = (n - 1)! x^-n sum_k n / (x + k) (x / (x + k))^n` in log space,
        // e.g. `polygamma(30, 1e11)` was -0 instead of -8.84e-300.
        if x > 0.0 && x.is_finite() && (g.is_infinite() || z < f64::MIN_POSITIVE) {
            let nf = f64::from(n);
            let s = zeta_euler_maclaurin(n1p, x, |a| nf / a * (x / a).powf(nf));
            return sign * (crate::gammaln(nf) - nf * x.ln() + s.ln()).exp();
        }
        sign * g * z
    }
}

#[cfg(test)]
mod tests {
    //! Translated from `scipy.special.tests.TestPolygamma` at
    //! <https://github.com/scipy/scipy/blob/38ecfef/scipy/special/tests/test_basic.py#L4104-L4127>

    use crate::np_assert_allclose;

    #[test]
    fn test_polygamma() {
        // poly2 = special.polygamma(2, 1)
        let poly2 = crate::polygamma(2, 1.0);
        // poly3 = special.polygamma(3, 1)
        let poly3 = crate::polygamma(3, 1.0);
        // assert_allclose(poly2, -2.4041138063, atol=1.5e-10, rtol=0)
        np_assert_allclose!([poly2], [-2.404_113_806_3], atol = 1.5e-10, rtol = 0.0);
        // assert_allclose(poly3, 6.4939394023, atol=1.5e-10, rtol=0)
        np_assert_allclose!([poly3], [6.493_939_402_3], atol = 1.5e-10, rtol = 0.0);

        // # Test polygamma(0, x) == psi(x)
        // x = [2, 3, 1.1e14]
        let x = [2.0, 3.0, 1.1e14];
        // assert_allclose(special.polygamma(0, x), special.psi(x),
        //                 atol=1.5e-7, rtol=0)
        np_assert_allclose!(
            x.map(|xi| crate::polygamma(0, xi)),
            x.map(crate::digamma),
            atol = 1.5e-7,
            rtol = 0.0
        );

        // # Test broadcasting
        // n = [0, 1, 2]
        let n: [u32; 3] = [0, 1, 2];
        // x = [0.5, 1.5, 2.5]
        let x = [0.5, 1.5, 2.5];
        // expected = [-1.9635100260214238, 0.93480220054467933,
        //             -0.23620405164172739]
        let expected = [
            -1.963_510_026_021_423_8,
            0.934_802_200_544_679_3,
            -0.236_204_051_641_727_4,
        ];
        // assert_allclose(special.polygamma(n, x), expected, atol=1.5e-7, rtol=0)
        np_assert_allclose!(
            n.iter()
                .zip(x.iter())
                .map(|(&ni, &xi)| crate::polygamma(ni, xi))
                .collect::<Vec<f64>>(),
            expected,
            atol = 1.5e-7,
            rtol = 0.0
        );
    }

    #[test]
    fn test_polygamma_large_n() {
        // n! overflows or zeta(n + 1, x) underflows (reference values from mpmath)
        for (n, x, expected) in [
            (30, 1e11, -8.841_761_995_065_966e-300),
            (142, 227.031_298_928_706_3, -6.966_070_108_353_208e-92),
            (171, 1.02, 4.116_693_769_939_060_4e307),
            (175, 1000.0, 7.004_133_681_486_838e-210),
        ] {
            let rel = (crate::polygamma(n, x) - expected).abs() / expected.abs();
            assert!(rel < 5e-13, "polygamma({n}, {x}): rel = {rel:e}");
        }
        // subnormal results
        assert_eq!(crate::polygamma(170, 4_960.382_505_394_897), -5e-324);
        assert_eq!(crate::polygamma(1000, 773.186_605_102_581_9), -5e-324);
    }
}
