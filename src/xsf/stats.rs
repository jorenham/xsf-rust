use core::ffi::c_int;

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
    fn ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::ndtr(self) }
    }

    #[inline(always)]
    fn log_ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::log_ndtr(self) }
    }
}

impl StatsArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn ndtr(self) -> Self {
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
pub fn smirnov(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnov(n as c_int, x) }
}

/// Kolmogorov-Smirnov distribution function
pub fn smirnovc(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovc(n as c_int, x) }
}

/// Inverse of [`smirnov`]
pub fn smirnovi(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovi(n as c_int, x) }
}

/// Inverse of [`smirnovc`]
pub fn smirnovci(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovci(n as c_int, x) }
}

/// Derivative of [`smirnov`]
pub fn smirnovp(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovp(n as c_int, x) }
}

// Tukey-Lambda

/// Tukey-Lambda distribution function
#[doc(alias = "tklmbda")]
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
pub fn pdtri(k: i32, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtri(k as c_int, y) }
}

/// Poisson survival function
pub fn pdtrc(k: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtrc(k, m) }
}

// Binomial

/// Binomial distribution function
pub fn bdtr(k: f64, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtr(k, n as c_int, p) }
}

/// Binomial survival function
pub fn bdtrc(k: f64, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtrc(k, n as c_int, p) }
}

/// Binomial quantile function
pub fn bdtri(k: f64, n: i32, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtri(k, n as c_int, y) }
}

// Negative Binomial

/// Negative binomial distribution function
pub fn nbdtr(k: i32, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtr(k as c_int, n as c_int, p) }
}

/// Negative binomial survival function
pub fn nbdtrc(k: i32, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtrc(k as c_int, n as c_int, p) }
}

/// Negative binomial quantile function
pub fn nbdtri(k: i32, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtri(k as c_int, n as c_int, p) }
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[test]
    fn test_ndtr_f64() {
        crate::xsref::test("ndtr", "d-d", |x| crate::ndtr(x[0]));
    }

    // TODO(@jorenham): manual `log_ndtr` tests (jorenham/xsf-rust#97)

    #[test]
    fn test_ndtr_c64() {
        crate::xsref::test("ndtr", "cd-cd", |x| crate::ndtr(c64(x[0], x[1])));
    }

    #[test]
    fn test_ndtri() {
        crate::xsref::test("ndtri", "d-d", |x| crate::ndtri(x[0]));
    }

    #[test]
    fn test_kolmogorov() {
        crate::xsref::test("kolmogorov", "d-d", |x| crate::kolmogorov(x[0]));
    }

    #[test]
    fn test_kolmogc() {
        crate::xsref::test("kolmogc", "d-d", |x| crate::kolmogc(x[0]));
    }

    #[test]
    fn test_kolmogi() {
        crate::xsref::test("kolmogi", "d-d", |x| crate::kolmogi(x[0]));
    }

    #[test]
    fn test_kolmogci() {
        crate::xsref::test("kolmogci", "d-d", |x| crate::kolmogci(x[0]));
    }

    #[test]
    fn test_kolmogp() {
        crate::xsref::test("kolmogp", "d-d", |x| crate::kolmogp(x[0]));
    }

    #[test]
    fn test_smirnov() {
        crate::xsref::test("smirnov", "p_d-d", |x| crate::smirnov(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovc() {
        crate::xsref::test("smirnovc", "p_d-d", |x| crate::smirnovc(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovi() {
        crate::xsref::test("smirnovi", "p_d-d", |x| crate::smirnovi(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovci() {
        crate::xsref::test("smirnovci", "p_d-d", |x| {
            crate::smirnovci(x[0] as i32, x[1])
        });
    }

    #[test]
    fn test_smirnovp() {
        crate::xsref::test("smirnovp", "p_d-d", |x| crate::smirnovp(x[0] as i32, x[1]));
    }

    #[test]
    fn test_owens_t() {
        crate::xsref::test("owens_t", "d_d-d", |x| crate::owens_t(x[0], x[1]));
    }

    #[test]
    fn test_chdtr() {
        crate::xsref::test("chdtr", "d_d-d", |x| crate::chdtr(x[0], x[1]));
    }

    #[test]
    fn test_chdtrc() {
        crate::xsref::test("chdtrc", "d_d-d", |x| crate::chdtrc(x[0], x[1]));
    }

    #[test]
    fn test_chdtri() {
        crate::xsref::test("chdtri", "d_d-d", |x| crate::chdtri(x[0], x[1]));
    }

    #[test]
    fn test_fdtr() {
        crate::xsref::test("fdtr", "d_d_d-d", |x| crate::fdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtrc() {
        crate::xsref::test("fdtrc", "d_d_d-d", |x| crate::fdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtri() {
        crate::xsref::test("fdtri", "d_d_d-d", |x| crate::fdtri(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtr() {
        crate::xsref::test("gdtr", "d_d_d-d", |x| crate::gdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtrc() {
        crate::xsref::test("gdtrc", "d_d_d-d", |x| crate::gdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_pdtr() {
        crate::xsref::test("pdtr", "d_d-d", |x| crate::pdtr(x[0], x[1]));
    }

    #[test]
    fn test_pdtrc() {
        crate::xsref::test("pdtrc", "d_d-d", |x| crate::pdtrc(x[0], x[1]));
    }

    #[test]
    fn test_pdtri() {
        crate::xsref::test("pdtri", "p_d-d", |x| crate::pdtri(x[0] as i32, x[1]));
    }

    #[test]
    fn test_bdtr() {
        crate::xsref::test("bdtr", "d_p_d-d", |x| crate::bdtr(x[0], x[1] as i32, x[2]));
    }

    #[test]
    fn test_bdtrc() {
        crate::xsref::test("bdtrc", "d_p_d-d", |x| {
            crate::bdtrc(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_bdtri() {
        crate::xsref::test("bdtri", "d_p_d-d", |x| {
            crate::bdtri(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtr() {
        crate::xsref::test("nbdtr", "p_p_d-d", |x| {
            crate::nbdtr(x[0] as i32, x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtrc() {
        crate::xsref::test("nbdtrc", "p_p_d-d", |x| {
            crate::nbdtrc(x[0] as i32, x[1] as i32, x[2])
        });
    }
}
