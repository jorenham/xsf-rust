use core::ffi::c_int;
use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait StatsArg: sealed::Sealed {
    fn ndtr(self) -> Self;
    fn log_ndtr(self) -> Self;
}

impl StatsArg for f64 {
    #[inline(always)]
    fn ndtr(self) -> f64 {
        unsafe { crate::ffi::xsf::ndtr(self) }
    }

    #[inline(always)]
    fn log_ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::log_ndtr(self) }
    }
}

impl StatsArg for Complex<f64> {
    #[inline(always)]
    fn ndtr(self) -> Complex<f64> {
        unsafe { crate::ffi::xsf::ndtr_1(self.into()) }.into()
    }

    #[inline(always)]
    fn log_ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::log_ndtr_1(self.into()) }.into()
    }
}

// Normal

/// Normal distribution function `F(z)` for real or complex `z`
pub fn ndtr<T: StatsArg>(z: T) -> T {
    z.ndtr()
}

/// Log of [`ndtr`] for real or complex argument
pub fn log_ndtr<T: StatsArg>(z: T) -> T {
    z.log_ndtr()
}

/// Inverse of [`ndtr`]
pub fn ndtri(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::ndtri(x) }
}

/// Owen's T function
pub fn owens_t(h: f64, a: f64) -> f64 {
    unsafe { crate::ffi::xsf::owens_t(h, a) }
}

// Kolmogorov

/// Kolmogorov survival function
pub fn kolmogorov(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogorov(x) }
}

/// Inverse of [`kolmogorov`]
pub fn kolmogi(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogi(x) }
}

/// Kolmogorov distribution function
pub fn kolmogc(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogc(x) }
}

/// Inverse of [`kolmogc`]
pub fn kolmogci(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogci(x) }
}

/// Derivative of [`kolmogorov`]
pub fn kolmogp(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogp(x) }
}

// Kolmogorov-Smirnov

/// Kolmogorov-Smirnov survival function
pub fn smirnov(n: c_int, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnov(n, x) }
}

/// Kolmogorov-Smirnov distribution function
pub fn smirnovc(n: c_int, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovc(n, x) }
}

/// Inverse of [`smirnov`]
pub fn smirnovi(n: c_int, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovi(n, x) }
}

/// Inverse of [`smirnovc`]
pub fn smirnovci(n: c_int, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovci(n, x) }
}

/// Derivative of [`smirnov`]
pub fn smirnovp(n: c_int, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovp(n, x) }
}

// Tukey-Lambda

/// Tukey-Lambda distribution function
pub fn tukeylambdacdf(x: f64, lmbda: f64) -> f64 {
    unsafe { crate::ffi::xsf::tukeylambdacdf(x, lmbda) }
}

// Chi-squared

/// Chi-squared distribution function
pub fn chdtr(df: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::chdtr(df, x) }
}

/// Chi-squared survival function
pub fn chdtrc(df: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::chdtrc(df, x) }
}

/// Chi-squared quantile function
pub fn chdtri(df: f64, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::chdtri(df, y) }
}

// F

/// F distribution function
pub fn fdtr(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::fdtr(a, b, x) }
}

/// F survival function
pub fn fdtrc(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::fdtrc(a, b, x) }
}

/// F quantile function
pub fn fdtri(a: f64, b: f64, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::fdtri(a, b, y) }
}

// Gamma

/// Gamma distribution function
pub fn gdtr(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gdtr(a, b, x) }
}

/// Gamma survival function
pub fn gdtrc(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gdtrc(a, b, x) }
}

// Poisson

/// Poisson distribution function
pub fn pdtr(k: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtr(k, m) }
}

/// Poisson quantile function
pub fn pdtri(k: c_int, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtri(k, y) }
}

/// Poisson survival function
pub fn pdtrc(k: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtrc(k, m) }
}

// Binomial

/// Binomial distribution function
pub fn bdtr(k: f64, n: c_int, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtr(k, n, p) }
}

/// Binomial survival function
pub fn bdtrc(k: f64, n: c_int, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtrc(k, n, p) }
}

/// Binomial quantile function
pub fn bdtri(k: f64, n: c_int, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtri(k, n, y) }
}

// Negative Binomial

/// Negative binomial distribution function
pub fn nbdtr(k: c_int, n: c_int, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtr(k, n, p) }
}

/// Negative binomial survival function
pub fn nbdtrc(k: c_int, n: c_int, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtrc(k, n, p) }
}

/// Negative binomial quantile function
pub fn nbdtri(k: c_int, n: c_int, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtri(k, n, p) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;
    use num_complex::{Complex, c64};

    #[test]
    fn test_ndtr_f64() {
        testing::test::<f64, _>("ndtr", "d-d", |x: &[f64]| ndtr(x[0]));
    }

    // TODO: no log_ndtr xsref tables => needs manual smoketests

    #[test]
    fn test_ndtr_c64() {
        testing::test::<Complex<f64>, _>("ndtr", "cd-cd", |x: &[f64]| ndtr(c64(x[0], x[1])));
    }

    #[test]
    fn test_ndtri() {
        testing::test::<f64, _>("ndtri", "d-d", |x: &[f64]| ndtri(x[0]));
    }

    #[test]
    fn test_kolmogorov() {
        testing::test::<f64, _>("kolmogorov", "d-d", |x: &[f64]| kolmogorov(x[0]));
    }

    #[test]
    fn test_kolmogc() {
        testing::test::<f64, _>("kolmogc", "d-d", |x: &[f64]| kolmogc(x[0]));
    }

    #[test]
    fn test_kolmogi() {
        testing::test::<f64, _>("kolmogi", "d-d", |x: &[f64]| kolmogi(x[0]));
    }

    #[test]
    fn test_kolmogci() {
        testing::test::<f64, _>("kolmogci", "d-d", |x: &[f64]| kolmogci(x[0]));
    }

    #[test]
    fn test_kolmogp() {
        testing::test::<f64, _>("kolmogp", "d-d", |x: &[f64]| kolmogp(x[0]));
    }

    #[test]
    fn test_smirnov() {
        testing::test::<f64, _>("smirnov", "p_d-d", |x: &[f64]| smirnov(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovc() {
        testing::test::<f64, _>("smirnovc", "p_d-d", |x: &[f64]| smirnovc(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovi() {
        testing::test::<f64, _>("smirnovi", "p_d-d", |x: &[f64]| smirnovi(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovci() {
        testing::test::<f64, _>("smirnovci", "p_d-d", |x: &[f64]| {
            smirnovci(x[0] as i32, x[1])
        });
    }

    #[test]
    fn test_smirnovp() {
        testing::test::<f64, _>("smirnovp", "p_d-d", |x: &[f64]| smirnovp(x[0] as i32, x[1]));
    }

    #[test]
    fn test_owens_t() {
        testing::test::<f64, _>("owens_t", "d_d-d", |x: &[f64]| owens_t(x[0], x[1]));
    }

    #[test]
    fn test_chdtr() {
        testing::test::<f64, _>("chdtr", "d_d-d", |x: &[f64]| chdtr(x[0], x[1]));
    }

    #[test]
    fn test_chdtrc() {
        testing::test::<f64, _>("chdtrc", "d_d-d", |x: &[f64]| chdtrc(x[0], x[1]));
    }

    #[test]
    fn test_chdtri() {
        testing::test::<f64, _>("chdtri", "d_d-d", |x: &[f64]| chdtri(x[0], x[1]));
    }

    #[test]
    fn test_fdtr() {
        testing::test::<f64, _>("fdtr", "d_d_d-d", |x: &[f64]| fdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtrc() {
        testing::test::<f64, _>("fdtrc", "d_d_d-d", |x: &[f64]| fdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtri() {
        testing::test::<f64, _>("fdtri", "d_d_d-d", |x: &[f64]| fdtri(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtr() {
        testing::test::<f64, _>("gdtr", "d_d_d-d", |x: &[f64]| gdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtrc() {
        testing::test::<f64, _>("gdtrc", "d_d_d-d", |x: &[f64]| gdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_pdtr() {
        testing::test::<f64, _>("pdtr", "d_d-d", |x: &[f64]| pdtr(x[0], x[1]));
    }

    #[test]
    fn test_pdtrc() {
        testing::test::<f64, _>("pdtrc", "d_d-d", |x: &[f64]| pdtrc(x[0], x[1]));
    }

    #[test]
    fn test_pdtri() {
        testing::test::<f64, _>("pdtri", "p_d-d", |x: &[f64]| pdtri(x[0] as i32, x[1]));
    }

    #[test]
    fn test_bdtr() {
        testing::test::<f64, _>("bdtr", "d_p_d-d", |x: &[f64]| bdtr(x[0], x[1] as i32, x[2]));
    }

    #[test]
    fn test_bdtrc() {
        testing::test::<f64, _>("bdtrc", "d_p_d-d", |x: &[f64]| {
            bdtrc(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_bdtri() {
        testing::test::<f64, _>("bdtri", "d_p_d-d", |x: &[f64]| {
            bdtri(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtr() {
        testing::test::<f64, _>("nbdtr", "p_p_d-d", |x: &[f64]| {
            nbdtr(x[0] as i32, x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtrc() {
        testing::test::<f64, _>("nbdtrc", "p_p_d-d", |x: &[f64]| {
            nbdtrc(x[0] as i32, x[1] as i32, x[2])
        });
    }
}
