use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait Hyp1F1Arg: sealed::Sealed {
    fn hyp1f1(self, a: f64, b: f64) -> Self;
}

impl Hyp1F1Arg for f64 {
    #[inline(always)]
    fn hyp1f1(self, a: f64, b: f64) -> Self {
        unsafe { crate::ffi::xsf::hyp1f1(a, b, self) }
    }
}

impl Hyp1F1Arg for Complex<f64> {
    #[inline(always)]
    fn hyp1f1(self, a: f64, b: f64) -> Self {
        unsafe { crate::ffi::xsf::hyp1f1_1(a, b, self.into()) }.into()
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
/// [scipy-bern]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bernoulli.html
pub fn bernoulli<const N: usize>() -> [f64; N] {
    let mut out = [0.0; N];
    if N < 3 {
        if N >= 1 {
            out[0] = 1.0;
            if N >= 2 {
                out[1] = -0.5;
            }
        }
    } else {
        unsafe { crate::ffi::xsf::bernob((N - 1) as i32, out.as_mut_ptr()) };
    };
    out
}

/// Confluent hypergeometric function `1F1(a; b; z)` for real or complex `z`
pub fn hyp1f1<T: Hyp1F1Arg>(a: f64, b: f64, z: T) -> T {
    z.hyp1f1(a, b)
}

/// Confluent hypergeometric function `U(a,b,x)` for `x > 0`
#[doc(alias = "hyperu")]
pub fn hypu(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::hypu(a, b, x) }
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
