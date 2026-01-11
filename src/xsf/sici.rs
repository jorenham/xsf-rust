use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex64 {}
}

pub trait SiciArg: sealed::Sealed + Sized {
    fn sici(self) -> (Self, Self);
    fn shichi(self) -> (Self, Self);
}

impl SiciArg for f64 {
    #[inline]
    fn sici(self) -> (Self, Self) {
        let (mut si, mut ci) = (f64::NAN, f64::NAN);
        unsafe {
            crate::ffi::xsf::sici(self, &raw mut si, &raw mut ci);
        }
        (si, ci)
    }

    #[inline]
    fn shichi(self) -> (Self, Self) {
        let (mut shi, mut chi) = (f64::NAN, f64::NAN);
        unsafe {
            crate::ffi::xsf::shichi(self, &raw mut shi, &raw mut chi);
        }
        (shi, chi)
    }
}

impl SiciArg for num_complex::Complex64 {
    #[inline]
    fn sici(self) -> (Self, Self) {
        let (mut si, mut ci) = (Complex::new(f64::NAN, 0.0), Complex::new(f64::NAN, 0.0));
        unsafe {
            crate::ffi::xsf::sici_1(self, &raw mut si, &raw mut ci);
        }
        (si, ci)
    }

    #[inline]
    fn shichi(self) -> (Self, Self) {
        let (mut shi, mut chi) = (Complex::new(f64::NAN, 0.0), Complex::new(f64::NAN, 0.0));
        unsafe {
            crate::ffi::xsf::shichi_1(self, &raw mut shi, &raw mut chi);
        }
        (shi, chi)
    }
}

/// Sine and cosine integrals.
///
/// Corresponds to [`scipy.special.sici`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sici.html
///
/// # Definitions
///
/// $$
/// \begin{align*}
/// \Si(z) &= \int_0^z {\sin t \over t} \dd t \\\\
/// \Ci(z) &= \gamma + \ln(z) + \int_0^z {\cos t - 1 \over t} \dd t
/// \end{align*}
/// $$
///
/// where $\gamma$ is Euler's constant and $\ln$ is the principal branch of the logarithm.
///
/// # Arguments
/// - $z$ - real- or complex-valued argument
///
/// # Returns
/// - $\Si(z)$ - Sine integral
/// - $\Ci(z)$ - Cosine integral
///
/// # See also
/// - [`shichi`]: $\Shi(z)$ and $\Chi(z)$
/// - [`exp1`](crate::exp1): Exponential integral $E_1$
/// - [`expi`](crate::expi): Exponential integral $Ei$
#[inline]
pub fn sici<T: SiciArg>(z: T) -> (T, T) {
    z.sici()
}

/// Hyperbolic sine and cosine integrals.
///
/// Corresponds to [`scipy.special.shichi`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.shichi.html
///
/// # Definitions
///
/// $$
/// \begin{align*}
/// \Shi(z) &= \int_0^z {\sinh t \over t} \dd t \\\\
/// \Chi(z) &= \gamma + \ln(z) + \int_0^z {\cosh t - 1 \over t} \dd t
/// \end{align*}
/// $$
///
/// where $\gamma$ is Euler's constant and $\ln$ is the principal branch of the logarithm.
///
/// # Arguments
/// - $z$ - real- or complex-valued argument
///
/// # Returns
/// - $\Shi(z)$ - Hyperbolic sine integral
/// - $\Chi(z)$ - Hyperbolic cosine integral
///
/// # See also
/// - [`sici`]: $\Si(z)$ and $\Ci(z)$
/// - [`exp1`](crate::exp1): Exponential integral $E_1$
/// - [`expi`](crate::expi): Exponential integral $Ei$
#[inline]
pub fn shichi<T: SiciArg>(z: T) -> (T, T) {
    z.shichi()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_sici_f64() {
        xsref::test("sici", "d-d_d", |x| crate::sici(x[0]));
    }

    #[test]
    fn test_sici_c64() {
        xsref::test("sici", "cd-cd_cd", |x| crate::sici(c64(x[0], x[1])));
    }

    #[test]
    fn test_shichi_f64() {
        xsref::test("shichi", "d-d_d", |x| crate::shichi(x[0]));
    }

    #[test]
    fn test_shichi_c64() {
        xsref::test("shichi", "cd-cd_cd", |x| crate::shichi(c64(x[0], x[1])));
    }
}
