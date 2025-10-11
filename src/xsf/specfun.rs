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

/// Confluent hypergeometric function `1F1(a; b; z)` for real or complex `z`
pub fn hyp1f1<T: Hyp1F1Arg>(a: f64, b: f64, z: T) -> T {
    z.hyp1f1(a, b)
}

#[doc(alias = "huperu")]
/// Confluent hypergeometric function `U(a,b,x)` for `x > 0`
pub fn hypu(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::hypu(a, b, x) }
}

/// Associated Legendre function for `|x| ≤ 1`
pub fn pmv(m: i64, v: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::pmv(m as f64, v, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;
    use num_complex::{Complex, c64};

    #[test]
    fn test_hypu() {
        // the table is called "hyperu" instead of "hypu"
        xsref::test::<f64, _>("hyperu", "d_d_d-d", |x: &[f64]| hypu(x[0], x[1], x[2]));
    }

    #[test]
    fn test_hyp1f1() {
        xsref::test::<Complex<f64>, _>("hyp1f1", "d_d_cd-cd", |x: &[f64]| {
            hyp1f1(x[0], x[1], c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_pmv() {
        xsref::test::<f64, _>("pmv", "d_d_d-d", |x: &[f64]| pmv(x[0] as i64, x[1], x[2]));
    }
}
