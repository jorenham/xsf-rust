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
        unsafe { crate::ffi::xsf::hyp2f1_1(a, b, c, self.into()) }.into()
    }
}

/// Gauss hypergeometric function `2F1(a, b; c; z)` for real or complex `z`
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
