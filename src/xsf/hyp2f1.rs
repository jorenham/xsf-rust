mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait Hyp2F1Arg: sealed::Sealed {
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self;
}

impl Hyp2F1Arg for f64 {
    #[inline(always)]
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self {
        unsafe { crate::ffi::xsf::hyp2f1(self, a, b, c) }
    }
}

impl Hyp2F1Arg for num_complex::Complex<f64> {
    #[inline(always)]
    fn hyp2f1(self, a: f64, b: f64, c: f64) -> Self {
        unsafe { crate::ffi::xsf::hyp2f1_1(a, b, c, self) }
    }
}

/// Gauss' hypergeometric function $_2F_1\[a_1,a_2;\\,b\\,\rvert\\,z\]$ for real or complex $z$
///
/// Corresponds to [`scipy.special.hyp2f1`][hyp2f1] in SciPy, and accepts both `f64` and
/// `num_complex::Complex<f64>` inputs for `z`.
///
/// [hyp2f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp2f1.html
///
/// # Notes
///
/// This function is defined for $\|z\| < 1$ as
///
/// $$
/// \hyp{2}{1}{a_1,\\, a_2}{b}{\Big\|\\, z} = \sum\_{n=0}^\infty
///     \frac{\rpow{a_1}{n} \ \rpow{a_2}{n} }{\rpow{b}{n} }
///     \frac{z^n}{n!}
/// ,
/// $$
///
/// and defined on the rest of the complex $z$-plane by analytic continuation [^1].
/// Here $\rpow{\square}{n}$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
/// When $n$ is a non-negative integer the result is a polynomial of degree $n$.
///
/// The implementation for complex values of $z$ is described in [^2], except for $z$ in the region
/// defined by
///
/// $$
/// 0.9 \le \| z \| < 1.1 , \\
/// \|1 - z \| \ge 0.9 , \\
/// \mathrm{Re}(z) \ge 0
/// $$
///
/// in which the implementation follows [^4].
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/15.2>
/// [^2]: S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996
/// [^3]: Cephes Mathematical Functions Library, <http://www.netlib.org/cephes/>
/// [^4]: J.L. Lopez and N.M. Temme, "New series expansions of the Gauss hypergeometric function",
/// Adv Comput Math 39, 349-365 (2013). <https://doi.org/10.1007/s10444-012-9283-y>
///
/// # See also
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function, $_0F_1\[b\\,\rvert\\,z\]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function,
///   $\hyp{1}{1}{a}{b}{\big\|\\,z}$
///
pub fn hyp2f1<T: Hyp2F1Arg>(a: f64, b: f64, c: f64, z: T) -> T {
    z.hyp2f1(a, b, c)
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_hyp2f1_f64() {
        crate::xsref::test("hyp2f1", "d_d_d_d-d", |x| {
            crate::hyp2f1(x[0], x[1], x[2], x[3])
        });
    }

    #[test]
    fn test_hyp2f1_c64() {
        crate::xsref::test("hyp2f1", "d_d_d_cd-cd", |x| {
            crate::hyp2f1(x[0], x[1], x[2], c64(x[3], x[4]))
        });
    }
}
