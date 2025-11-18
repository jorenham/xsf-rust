use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait Hyp1F1Arg: sealed::Sealed {
    fn hyp1f1(&self, a: f64, b: f64) -> Self;
}

impl Hyp1F1Arg for f64 {
    #[inline(always)]
    fn hyp1f1(&self, a: f64, b: f64) -> Self {
        unsafe { crate::ffi::xsf::hyp1f1(a, b, *self) }
    }
}

impl Hyp1F1Arg for Complex<f64> {
    #[inline(always)]
    fn hyp1f1(&self, a: f64, b: f64) -> Self {
        unsafe { crate::ffi::xsf::hyp1f1_1(a, b, *self) }
    }
}

/// Bernoulli numbers B<sub>0</sub>, ..., B<sub>N-1</sub>
///
/// Corresponds to [`scipy.special.bernoulli`][scipy-bern] in SciPy, and calls the FFI function
/// `xsf::specfun::bernob`.
///
/// # Examples
/// ```
/// use xsf::bernoulli;
/// assert_eq!(bernoulli::<0>(), []);
/// assert_eq!(bernoulli::<1>(), [1.0]);
/// assert_eq!(bernoulli::<2>(), [1.0, -0.5]);
/// assert_eq!(bernoulli::<4>(), [1.0, -0.5, 1.0 / 6.0, 0.0]);
/// assert_eq!(bernoulli::<1000>()[999], 0.0);
/// ```
///
/// # See also
/// - [`euler`]: Euler numbers E<sub>0</sub>, ..., E<sub>N-1</sub>
///
/// [scipy-bern]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bernoulli.html
pub fn bernoulli<const N: usize>() -> [f64; N] {
    let mut bn = [0.0; N];
    if N < 3 {
        // `bernob` requires `N >= 3` (i.e. `n >= 2`), so handle small `N` manually
        if N >= 1 {
            bn[0] = 1.0;
            if N >= 2 {
                bn[1] = -0.5;
            }
        }
    } else {
        unsafe { crate::ffi::xsf::bernob((N - 1) as i32, bn.as_mut_ptr()) };
    };
    bn
}

/// Euler numbers E<sub>0</sub>, ..., E<sub>N-1</sub>
///
/// Corresponds to [`scipy.special.euler`][scipy-euler] in SciPy, and calls the FFI function
/// `xsf::specfun::eulerb`.
///
/// # Examples
/// ```
/// use xsf::euler;
/// assert_eq!(euler::<0>(), []);
/// assert_eq!(euler::<1>(), [1.0]);
/// assert_eq!(euler::<2>(), [1.0, 0.0]);
/// assert_eq!(euler::<4>(), [1.0, 0.0, -1.0, 0.0]);
/// assert_eq!(euler::<1000>()[999], 0.0);
/// ```
///
/// # See also
/// - [`bernoulli`]: Bernoulli numbers B<sub>0</sub>, ..., B<sub>N-1</sub>
///
/// [scipy-euler]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.euler.html
pub fn euler<const N: usize>() -> [f64; N] {
    let mut en = [0.0; N];
    if N < 3 {
        // `eulerb` requires `N >= 3` (i.e. `n >= 2`), so handle small `N` manually
        if N >= 1 {
            en[0] = 1.0;
        }
    } else {
        unsafe { crate::ffi::xsf::eulerb((N - 1) as i32, en.as_mut_ptr()) };
    };
    en
}

/// Kummer's Confluent hypergeometric function $_1F_1$
///
/// Corresponds to [`scipy.special.hyp1f1`][hyp1f1] in SciPy, and accepts both `f64` and
/// `num_complex::Complex<f64>` inputs for `z`.
///
/// [hyp1f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp1f1.html
///
/// # Notes
///
/// The confluent hypergeometric function is defined by the series
///
/// $$
/// M(a,b,z) = \hyp 1 1 a b z = \sum\_{n=0}^\infty  {\rpow a n \over \rpow b n} {z^n \over n!}
/// $$
///
/// See [^DLMF] for more details.
/// Here $\rpow{\square}{n}$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.2#E2>
///
/// # See also
/// - [`hypu`](crate::hypu): Tricomi's confluent hypergeometric function $U(a,b,x)$
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function,
///   $_0F_1\left[b\middle\| z\right]$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
pub fn hyp1f1<T: Hyp1F1Arg>(a: f64, b: f64, z: T) -> T {
    z.hyp1f1(a, b)
}

#[inline(always)]
fn xsf_hypu(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::hypu(a, b, x) }
}

/// Tricomi's confluent hypergeometric function $U(a,b,x)$
///
/// Corresponds to [`scipy.special.hyperu`][hyperu] in SciPy.
///
/// [hyperu]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyperu.html
///
/// # Notes
///
/// It is defined as the solution to the equation
///
/// $$
/// x \frac{\dd^2 w}{\dd x^2} + (b - x) \frac{\dd w}{\dd x} - aw = 0
/// $$
///
/// which satisfies the property $U(a,b,x) \sim x^{-a}$ as $x \to \infty$.
/// See [^DLMF] for more details.
///
/// # See also
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function $M(a,b,z)$
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.2#E6>
#[doc(alias = "hyperu")]
pub fn hypu(a: f64, b: f64, x: f64) -> f64 {
    if a.is_nan() || b.is_nan() || x.is_nan() {
        f64::NAN
    } else if x < 0.0 {
        // Domain error
        f64::NAN
    } else if x == 0.0 {
        if b > 1.0 {
            // Singular. DLMF 13.2.16-18
            f64::INFINITY
        } else {
            // DLMF 13.2.14-15 and 13.2.19-21
            unsafe { crate::ffi::xsf::poch(1.0 - b + a, -a) }
        }
    } else if x < 1.0 && b == 1.0 && a > -0.25 && a < 0.3 {
        // DLMF 13.3.7. Fixes scipy/scipy#15650
        let a1 = a + 1.0;
        (x + a + a1) * xsf_hypu(a1, 1.0, x) - a1 * a1 * xsf_hypu(a + 2.0, 1.0, x)
    } else {
        xsf_hypu(a, b, x)
    }
}

/// Associated Legendre function for `|x| ≤ 1`
pub fn pmv(m: i64, v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::pmv(m as f64, v, x) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_hypu() {
        // the table is called "hyperu" instead of "hypu"
        crate::xsref::test("hyperu", "d_d_d-d", |x| crate::hypu(x[0], x[1], x[2]));
    }

    #[test]
    fn test_hyp1f1() {
        crate::xsref::test("hyp1f1", "d_d_cd-cd", |x| {
            crate::hyp1f1(x[0], x[1], c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_pmv() {
        crate::xsref::test("pmv", "d_d_d-d", |x| crate::pmv(x[0] as i64, x[1], x[2]));
    }
}
