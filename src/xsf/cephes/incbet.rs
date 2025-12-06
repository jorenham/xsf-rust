/// Regularized incomplete Beta function $\I_x(a,b)$
///
/// Computes the regularized incomplete Beta function, defined as [^DLMF]:
///
/// $$
/// \I_x(a,b)
/// = {\Gamma(a+b) \over \Gamma(a)\Gamma(b)} \int_0^x t^{a-1} (1-t)^{b-1} \dd t
/// = {\B_x(a,b) \over \B(a,b)} ,
/// $$
///
/// for $0 \le x \le 1$, with $\B$ the Beta function, and $\B_x$ the (non-regularized)
/// incomplete Beta function.
///
/// # Notes
/// This function wraps the `incbet` Cephes routine [^CEPHES], making it less accurate than
/// [`scipy.special.betainc`][scipy], which wraps the Boost `ibeta` routine [^BOOST].
/// Especially for small `x` the accuracy may be poor.
///
/// # See also
/// - [`betaincinv`](crate::betaincinv): Inverse of the regularized incomplete Beta function
/// - [`beta`](crate::beta): Beta function $\B(a, b)$
/// - [`gamma`](crate::gamma): Gamma function $\Gamma(x)$
/// - [`scipy.special.betainc`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/8.17>.
/// [^CEPHES]: Cephes Math Library Release 2.4. Translated into C++ by SciPy developers in 2024.
/// [^BOOST]: The Boost Developers. “Boost C++ Libraries”, <https://www.boost.org>.
#[doc(
    alias = "inc_beta",
    alias = "beta_inc",
    alias = "ibeta",
    alias = "incbet"
)]
#[must_use]
pub fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if (a < 0.0 || b < 0.0 || !(0.0..=1.0).contains(&x))
        || (a.is_nan() || b.is_nan() || x.is_nan())
        || (a == 0.0 && b == 0.0)
        || (a.is_infinite() && b.is_infinite())
    {
        f64::NAN
    } else if a == 0.0 || b.is_infinite() {
        if x > 0.0 { 1.0 } else { 0.0 }
    } else if a.is_infinite() || b == 0.0 {
        if x < 1.0 { 0.0 } else { 1.0 }
    } else {
        unsafe { crate::ffi::xsf::incbet(a, b, x) }
    }
}

#[cfg(test)]
mod tests {
    // based on scipy.special.tests.test_basic.TestBetaInc

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_betainc_a1_b1() {
        // betainc(1, 1, x) is x.
        let x = [0.0, 0.25, 1.0];
        assert_eq!(x.map(|xv| crate::betainc(1.0, 1.0, xv)), x);
    }

    #[test]
    fn test_betainc_nontrivial() {
        for &(a, b, x, p) in &[
            (2.0, 4.0, 0.313_810_170_455_697_4, 0.5),
            (0.0342, 171.0, 1e-10, 0.552_699_169_018_070_9),
            // scipy/scipy#3761:
            (0.0342, 171.0, 8.423_131_693_547_97e-21, 0.25),
            // scipy/scipy#4244 (NOTE: this results in a relative error of 5.134e-10)
            // (
            //     0.0002742794749792665,
            //     289206.03125,
            //     1.639984034231756e-56,
            //     0.9688708782196045,
            // ),
            // scipy/scipy#12796:
            (4.0, 99997.0, 0.000_194_784_157_889_212_1, 0.999_995),
        ] {
            let p1 = crate::betainc(a, b, x);
            // NOTE: The original `rtol = 1e-15` is too strict for Cephes
            crate::np_assert_allclose!(&[p1], &[p], rtol = 1e-14);
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
            let observed = crate::betainc(a, b, x);
            crate::np_assert_equal!(&[observed], &[expected]);
        }
    }

    #[test]
    fn test_betainc_gh22682() {
        // scipy/scipy#22682: betainc returned incorrect results for tiny
        // single precision inputs. test that this is resolved
        let rtol = 1e-15;
        for &(a, b, x, reference) in &[
            (1e-20, 1e-21, 0.5, 0.090_909_090_909_090_9),
            (1e-15, 1e-16, 0.5, 0.090_909_090_909_090_91),
        ] {
            let res = crate::betainc(a, b, x);
            crate::np_assert_allclose!(&[res], &[reference], rtol = rtol);
        }
    }
}
