mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait ExpIntArg: sealed::Sealed {
    fn expi(self) -> Self;
    fn exp1(self) -> Self;
}

impl ExpIntArg for f64 {
    #[inline]
    fn expi(self) -> Self {
        unsafe { crate::ffi::xsf::expi(self) }
    }

    #[inline]
    fn exp1(self) -> Self {
        unsafe { crate::ffi::xsf::exp1(self) }
    }
}

impl ExpIntArg for num_complex::Complex<f64> {
    #[inline]
    fn expi(self) -> Self {
        unsafe { crate::ffi::xsf::expi_1(self) }
    }

    #[inline]
    fn exp1(self) -> Self {
        unsafe { crate::ffi::xsf::exp1_1(self) }
    }
}

/// Exponential integral $E_1$ for real or complex input
///
/// Corresponds to [`scipy.special.exp1`][exp1]
///
/// [exp1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp1.html
///
/// # Definition
/// For complex $z \neq 0$ the exponential integral can be defined as
///
/// $$ E_1(z) = \int_z^\infty \frac{e^{-t}}{t} \dd t $$
///
/// where the path of the integral does not cross the negative real axis or pass through the origin.
///
/// # See also
/// - [`scaled_exp1`]: scaled version of $E_1$
/// - [`expn`](crate::expn): generalization of $E_1$
/// - [`expi`]: exponential integral $Ei$
#[must_use]
#[inline]
pub fn exp1<T: ExpIntArg>(z: T) -> T {
    z.exp1()
}

/// Exponential integral $Ei$ for real or complex input
///
/// # Definition
/// For real $x$, the exponential integral is defined as
///
/// $$ Ei(x) = \int_{-\infty}^x \frac{e^t}{t} \dd t $$
///
/// For x > 0 the integral is understood as a Cauchy principal value.
///
/// It is extended to the complex plane by analytic continuation of the function on the interval
/// $(0, \infty)$. The complex variant has a branch cut on the negative real axis.
///
/// # See also
/// - [`exp1`]: exponential integral $E_1$
/// - [`expn`](crate::expn): generalization of $E_1$
#[must_use]
#[inline]
pub fn expi<T: ExpIntArg>(z: T) -> T {
    z.expi()
}

/// Scaled version of the exponential integral $E_1$ for real input
///
/// There exists no corresponding function in `scipy.special`.
///
/// # Definition
/// This function computes $F(x)$, where $F$ is the factor remaining in $E_1(x)$
/// when $e^{-x}/x$ is factored out. That is:
///
/// $$ E_1(x) = x^{-1} e^{-x} F(x) $$
///
/// or
///
/// $$ F(x) = x e^{x} E_1(x) $$
///
/// The function is defined for real $x \geq 0$. For $x < 0$, NaN is returned.
///
/// $F$ has the properties:
///
/// - $F(0) = 0$
/// - $F(x)$ is increasing on $[0, \infty)$.
/// - $F(x) = 1$ in the limit as $x \to \infty$.
///
/// # See also
/// - [`exp1`]: exponential integral $E_1$
/// - [`expi`]: exponential integral $Ei$
#[must_use]
#[inline]
pub fn scaled_exp1(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::scaled_exp1(x) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_exp1_f64() {
        xsref::test("exp1", "d-d", |x| crate::exp1(x[0]));
    }

    #[test]
    fn test_exp1_c64() {
        xsref::test("exp1", "cd-cd", |x| crate::exp1(c64(x[0], x[1])));
    }

    #[test]
    fn test_expi_f64() {
        xsref::test("expi", "d-d", |x| crate::expi(x[0]));
    }

    #[test]
    fn test_expi_c64() {
        xsref::test("expi", "cd-cd", |x| crate::expi(c64(x[0], x[1])));
    }

    #[test]
    fn test_scaled_exp1_f64() {
        xsref::test("scaled_exp1", "d-d", |x| crate::scaled_exp1(x[0]));
    }
}
