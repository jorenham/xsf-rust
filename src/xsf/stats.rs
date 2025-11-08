use core::f64;
use core::f64::consts::{FRAC_1_PI, PI};
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
        unsafe { crate::ffi::xsf::ndtr_1(self) }
    }

    #[inline(always)]
    fn log_ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::log_ndtr_1(self) }
    }
}

// Student's t

/// Student's t distribution cumulative distribution function
///
/// Rust implementation of [`scipy.special.stdtr`][stdtr].
///
/// # Notes
///
/// The CDF is given by
///
/// $$
/// \begin{align*}
/// F_\nu(x)
/// \&= \frac{1}{\nu \B \left( \frac{1}{2}, \frac{\nu}{2} \right)}
///     \int_{-\infty}^x g(u)^{\frac{\nu+1}{2}} \dd u \\\\
/// \&= 1 - \frac{1}{2} I_{g(x)}\left( \frac{1}{2}, \frac{\nu}{2} \right) ,
/// \end{align*}
/// $$
///
/// with $g(\square) = \frac{\nu}{\square^2 + \nu}$ a helper function, $\B(\cdot, \cdot)$ the Beta
/// function, and $I_\square(\cdot,\cdot)$ the regularized incomplete Beta function.
///
/// # See also
/// - [`stdtri`](stdtri): Inverse of the CDF
/// - [`beta`](crate::beta): Beta function
/// - [`betainc`](crate::betainc): Regularized incomplete Beta function
///
/// [stdtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtr.html
pub fn stdtr(nu: f64, x: f64) -> f64 {
    if nu <= 0.0 {
        f64::NAN
    } else if x == 0.0 {
        0.5
    } else if nu.is_infinite() {
        ndtr(x)
    } else if nu == 1.0 {
        0.5 + FRAC_1_PI * x.atan()
    } else if nu == 2.0 {
        0.5 * (1.0 + x / (2.0 + x * x).sqrt())
    } else {
        let p = 0.5 * crate::betainc(0.5 * nu, 0.5, nu / (nu + x * x));
        if x >= 0.0 { 1.0 - p } else { p }
    }
}

/// Inverse of [`stdtr`]
///
/// Rust implementation of [`scipy.special.stdtrit`][stdtrit] with comparable or better accuracy.
///
/// [stdtrit]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stdtrit.html
pub fn stdtri(nu: f64, p: f64) -> f64 {
    if nu <= 0.0 || !(0.0..=1.0).contains(&p) {
        f64::NAN
    } else if p == 0.0 {
        f64::NEG_INFINITY
    } else if p == 1.0 {
        f64::INFINITY
    } else if p == 0.5 {
        0.0
    } else if nu.is_infinite() {
        ndtri(p)
    } else if nu == 1.0 {
        (PI * (p - 0.5)).tan()
    } else if nu == 2.0 {
        (p - 0.5) * (2.0 / (p * (1.0 - p))).sqrt()
    } else {
        let x = crate::betaincinv(0.5 * nu, 0.5, 2.0 * p.min(1.0 - p));
        let x = (nu * (1.0 - x) / x).sqrt();
        if p < 0.5 { -x } else { x }
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
    use core::f64;
    use num_complex::c64;

    #[test]
    fn test_stdtr() {
        const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
        const X: [f64; 9] = [-8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0];
        // values from scipy.special.stdtr
        const P_REF: [[f64; 9]; 6] = [
            [
                0.113252546403789,
                0.159610041494336,
                0.222757445091566,
                0.301121610841322,
                0.5,
                0.698878389158678,
                0.777242554908434,
                0.840389958505664,
                0.886747453596211,
            ],
            [
                0.039583424160566,
                0.077979130377369,
                0.147583617650433,
                0.25,
                0.5,
                0.75,
                0.852416382349567,
                0.922020869622631,
                0.960416575839434,
            ],
            [
                0.007634036082669,
                0.028595479208968,
                0.091751709536137,
                0.211324865405187,
                0.5,
                0.788675134594813,
                0.908248290463863,
                0.971404520791032,
                0.992365963917331,
            ],
            [
                0.002038288793893,
                0.014004228005073,
                0.069662984279422,
                0.195501109477885,
                0.5,
                0.804498890522115,
                0.930337015720578,
                0.985995771994927,
                0.997961711206107,
            ],
            [
                0.000661948454609,
                0.008065044950046,
                0.058058261758408,
                0.18695048315003,
                0.5,
                0.81304951684997,
                0.941941738241592,
                0.991934955049954,
                0.999338051545391,
            ],
            [
                0.000005887471395,
                0.001259166312368,
                0.03669401738537,
                0.17044656615103,
                0.5,
                0.82955343384897,
                0.96330598261463,
                0.998740833687632,
                0.999994112528605,
            ],
        ];
        for (i, &nu) in NU.iter().enumerate() {
            let p = X.map(|x| crate::stdtr(nu, x));
            crate::np_assert_allclose!(p, P_REF[i], atol = 1e-15);
        }
    }

    #[test]
    fn test_stdtri() {
        const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
        const P: [f64; 7] = [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999];
        // values from Wolfram Alpha
        const X_REF: [[f64; 7]; 6] = [
            [
                -102849.11563017555,
                -1028.4910104716219,
                -10.270324410234506,
                0.0,
                10.270324410234506,
                1028.4910104716219,
                102849.11563017555,
            ],
            [
                -318.30883898555045,
                -31.82051595377396,
                -3.0776835371752534,
                0.0,
                3.0776835371752534,
                31.82051595377396,
                318.30883898555045,
            ],
            [
                -22.327124770119875,
                -6.964556734283274,
                -1.8856180831641267,
                0.0,
                1.8856180831641267,
                6.964556734283274,
                22.327124770119875,
            ],
            [
                -10.214531852407387,
                -4.5407028585681336,
                -1.6377443536962101,
                0.0,
                1.6377443536962101,
                4.5407028585681336,
                10.214531852407387,
            ],
            [
                -7.1731822197823085,
                -3.746947387979197,
                -1.5332062740589439,
                0.0,
                1.5332062740589439,
                3.746947387979197,
                7.1731822197823085,
            ],
            [
                -4.14370049404659,
                -2.763769458112696,
                -1.3721836411103356,
                0.0,
                1.3721836411103356,
                2.763769458112696,
                4.14370049404659,
            ],
        ];
        for (i, &nu) in NU.iter().enumerate() {
            let x = P.map(|p| crate::stdtri(nu, p));
            crate::np_assert_allclose!(x, X_REF[i], rtol = 1e-13);
        }
    }

    #[test]
    fn test_ndtr_f64() {
        crate::xsref::test("ndtr", "d-d", |x| crate::ndtr(x[0]));
    }

    #[test]
    fn test_ndtr_c64() {
        crate::xsref::test("ndtr", "cd-cd", |x| crate::ndtr(c64(x[0], x[1])));
    }

    // Based on `scipy.special.tests.test_ndtr.TestLogNdtr.test_log_ndtr_moderate_le8`
    #[test]
    fn test_log_ndtr() {
        let x = [-0.75, -0.25, 0.0, 0.5, 1.5, 2.5, 3.0, 4.0, 5.0, 7.0, 8.0];
        let expected = [
            -1.4844482299196562,
            -0.9130617648111351,
            #[allow(clippy::approx_constant)]
            -0.6931471805599453,
            -0.3689464152886564,
            -0.06914345561223398,
            -0.006229025485860002,
            -0.0013508099647481938,
            -3.167174337748927e-05,
            -2.866516129637636e-07,
            -1.279812543886654e-12,
            -6.220960574271786e-16,
        ];
        let y = x.map(crate::log_ndtr);
        crate::np_assert_allclose!(&y, &expected, rtol = 1e-14);
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
