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
    #[inline]
    fn ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::ndtr(self) }
    }

    #[inline]
    fn log_ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::log_ndtr(self) }
    }
}

impl StatsArg for num_complex::Complex<f64> {
    #[inline]
    fn ndtr(self) -> Self {
        unsafe { crate::ffi::xsf::ndtr_1(self) }
    }

    #[inline]
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
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
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
#[must_use]
#[inline]
#[allow(clippy::float_cmp)]
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

/// CDF of the standard normal distribution, $\Phi(z)$
///
/// Corresponds to [`scipy.special.ndtr`][ndtr].
///
/// [ndtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtr.html
///
/// # Notes
///
/// The CDF is given by
///
/// $$
/// \begin{align*}
/// \Phi(z)
/// &= {1 \over \sqrt{2\pi}} \int_{-\infty}^z e^{-t^2 / 2} \dd t \\\\
/// &= {1 \over 2} + {1 \over 2} \erf \left( {z \over \sqrt{2}} \right) ,
/// \end{align*}
/// $$
///
/// with $\erf(z)$ the [error function](crate::erf).
///
/// # See also
/// - [`log_ndtr`]: $\ln \Phi(z)$
/// - [`ndtri`]: Normal quantile function $\Phi^{-1}(z)$, a.k.a. the probit function
/// - [`erf`](crate::erf): Error function $\erf(z)$
#[must_use]
#[inline]
pub fn ndtr<T: StatsArg>(z: T) -> T {
    z.ndtr()
}

/// Logarithm of [`ndtr`], $\ln \Phi(z)$
///
/// Corresponds to [`scipy.special.log_ndtr`][log_ndtr].
///
/// [log_ndtr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.log_ndtr.html
///
/// # See also
/// - [`ndtr`]: Normal distribution function $\Phi(z)$
/// - [`crate::erf`]: Error function $\erf(z)$
#[must_use]
#[inline]
pub fn log_ndtr<T: StatsArg>(z: T) -> T {
    z.log_ndtr()
}

/// Inverse of [`ndtr`], the probit function $\Phi^{-1}(p)$
///
/// Corresponds to [`scipy.special.ndtri`][ndtri].
///
/// [ndtri]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri.html
///
/// # Notes
///
/// The normal quantile function (probit function) is given by
///
/// $$
/// \Phi^{-1}(p) = \sqrt{2}\\,\erf^{-1}(2p-1) ,
/// $$
///
/// with $\erf^{-1}(z)$ the [inverse error function](crate::erfinv) and $p \in (0,1)$.
///
/// # See also
/// - [`ndtr`]: CDF of the standard normal distribution, $\Phi(z)$
/// - [`crate::erfinv`]: Inverse error function $\erf^{-1}(z)$
#[must_use]
#[inline]
#[doc(alias = "probit")]
pub fn ndtri(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::ndtri(x) }
}

/// Owen's T function, $T(h, a)$
///
/// The function $T(h, a)$ gives the probability of the event $(X > h \wedge 0 < Y < a X)$ where
/// $X$ and $Y$ are independent standard normal random variables.
///
/// Corresponds to [`scipy.special.owens_t`][owens_t].
///
/// [owens_t]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.owens_t.html
///
/// # Definition
///
/// $$ T(h, a) = {1 \over 2\pi} \int_0^a {e^{-h^2 (1+x^2) / 2} \over 1+x^2} \dd x $$
///
/// # See also
/// - [`ndtr`]: CDF of the standard normal distribution, $\Phi(z)$
#[must_use]
#[inline]
pub fn owens_t(h: f64, a: f64) -> f64 {
    unsafe { crate::ffi::xsf::owens_t(h, a) }
}

// Kolmogorov

/// Survival function of the Kolmogorov distribution
///
/// Corresponds to [`scipy.special.kolmogorov`][kolmogorov].
///
/// [kolmogorov]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kolmogorov.html
///
/// # See also
/// - [`kolmogi`]: Inverse survival function
/// - [`kolmogc`]: Cumulative distribution function
/// - [`kolmogci`]: Quantile function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogorov(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogorov(x) }
}

/// Inverse of [`kolmogorov`]
///
/// Corresponds to [`scipy.special.kolmogi`][kolmogi].
///
/// [kolmogi]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kolmogi.html
///
/// # See also
/// - [`kolmogorov`]: Survival function
/// - [`kolmogc`]: Cumulative distribution function
/// - [`kolmogci`]: Quantile function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogi(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogi(x) }
}

/// CDF of the Kolmogorov distribution
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`kolmogorov`]: Survival function
/// - [`kolmogi`]: Inverse survival function
/// - [`kolmogci`]: Quantile function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogc(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogc(x) }
}

/// Inverse of [`kolmogc`], the quantile function of the Kolmogorov distribution
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`kolmogorov`]: Survival function
/// - [`kolmogi`]: Inverse survival function
/// - [`kolmogc`]: Cumulative distribution function
/// - [`kolmogp`]: Derivative of the survival function
#[must_use]
#[inline]
pub fn kolmogci(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogci(x) }
}

/// Derivative of [`kolmogorov`]
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`kolmogorov`]: Survival function of the Kolmogorov distribution
/// - [`kolmogi`]: Inverse of the Kolmogorov distribution function
/// - [`kolmogc`]: CDF of the Kolmogorov distribution
/// - [`kolmogci`]: Inverse of the Kolmogorov CDF
#[must_use]
#[inline]
pub fn kolmogp(x: f64) -> f64 {
    unsafe { crate::ffi::xsf::kolmogp(x) }
}

// Kolmogorov-Smirnov

/// Survival function of the Kolmogorov-Smirnov distribution
///
/// Corresponds to [`scipy.special.smirnov`][smirnov].
///
/// [smirnov]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.smirnov.html
///
/// # See also
/// - [`smirnovi`]: Inverse survival function
/// - [`smirnovc`]: Cumulative distribution function
/// - [`smirnovp`]: Derivative of the survival function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnov(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnov(n as c_int, x) }
}

/// Inverse of [`smirnov`]
///
/// Corresponds to [`scipy.special.smirnovi`][smirnovi].
///
/// [smirnovi]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.smirnovi.html
///
/// # See also
/// - [`smirnov`]: Survival function
/// - [`smirnovc`]: Cumulative distribution function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnovi(n: i32, q: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovi(n as c_int, q) }
}

/// CDF of the Kolmogorov-Smirnov distribution
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`smirnov`]: Survival function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnovc(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovc(n as c_int, x) }
}

/// Inverse of [`smirnovc`]
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`smirnovc`]: Cumulative distribution function
/// - [`smirnovci`]: Quantile function
#[must_use]
#[inline]
pub fn smirnovci(n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovci(n as c_int, p) }
}

/// Derivative of [`smirnov`]
///
/// Does not have a direct counterpart in SciPy.
///
/// # See also
/// - [`smirnov`]: Survival function
#[must_use]
#[inline]
pub fn smirnovp(n: i32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::smirnovp(n as c_int, x) }
}

// Tukey-Lambda

/// CDF of the Tukey-Lambda distribution
///
/// Corresponds to [`scipy.special.tklmbda`][tklmbda].
///
/// [tklmbda]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.tklmbda.html
#[doc(alias = "tklmbda")]
#[must_use]
#[inline]
pub fn tukeylambdacdf(x: f64, lambda: f64) -> f64 {
    unsafe { crate::ffi::xsf::tukeylambdacdf(x, lambda) }
}

// Chi-squared

/// Chi-squared distribution function
#[must_use]
#[inline]
pub fn chdtr(df: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::chdtr(df, x) }
}

/// Chi-squared survival function
#[must_use]
#[inline]
pub fn chdtrc(df: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::chdtrc(df, x) }
}

/// Chi-squared quantile function
#[must_use]
#[inline]
pub fn chdtri(df: f64, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::chdtri(df, y) }
}

// F

/// F distribution function
#[must_use]
#[inline]
pub fn fdtr(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::fdtr(a, b, x) }
}

/// F survival function
#[must_use]
#[inline]
pub fn fdtrc(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::fdtrc(a, b, x) }
}

/// F quantile function
#[must_use]
#[inline]
pub fn fdtri(a: f64, b: f64, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::fdtri(a, b, y) }
}

// Gamma

/// Gamma distribution function
#[must_use]
#[inline]
pub fn gdtr(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gdtr(a, b, x) }
}

/// Gamma survival function
#[must_use]
#[inline]
pub fn gdtrc(a: f64, b: f64, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::gdtrc(a, b, x) }
}

// Poisson

/// Poisson distribution function
#[must_use]
#[inline]
pub fn pdtr(k: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtr(k, m) }
}

/// Poisson quantile function
#[must_use]
#[inline]
pub fn pdtri(k: i32, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtri(k as c_int, y) }
}

/// Poisson survival function
#[must_use]
#[inline]
pub fn pdtrc(k: f64, m: f64) -> f64 {
    unsafe { crate::ffi::xsf::pdtrc(k, m) }
}

// Binomial

/// Binomial distribution function
#[must_use]
#[inline]
pub fn bdtr(k: f64, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtr(k, n as c_int, p) }
}

/// Binomial survival function
#[must_use]
#[inline]
pub fn bdtrc(k: f64, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtrc(k, n as c_int, p) }
}

/// Binomial quantile function
#[must_use]
#[inline]
pub fn bdtri(k: f64, n: i32, y: f64) -> f64 {
    unsafe { crate::ffi::xsf::bdtri(k, n as c_int, y) }
}

// Negative Binomial

/// Negative binomial distribution function
#[must_use]
#[inline]
pub fn nbdtr(k: i32, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtr(k as c_int, n as c_int, p) }
}

/// Negative binomial survival function
#[must_use]
#[inline]
pub fn nbdtrc(k: i32, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtrc(k as c_int, n as c_int, p) }
}

/// Negative binomial quantile function
#[must_use]
#[inline]
pub fn nbdtri(k: i32, n: i32, p: f64) -> f64 {
    unsafe { crate::ffi::xsf::nbdtri(k as c_int, n as c_int, p) }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

    use num_complex::c64;

    #[test]
    fn test_stdtr() {
        const NU: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 4.0, 10.0];
        const X: [f64; 9] = [-8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0];
        // values from scipy.special.stdtr
        const P_REF: [[f64; 9]; 6] = [
            [
                0.113_252_546_403_789,
                0.159_610_041_494_336,
                0.222_757_445_091_566,
                0.301_121_610_841_322,
                0.500_000_000_000_000,
                0.698_878_389_158_678,
                0.777_242_554_908_434,
                0.840_389_958_505_664,
                0.886_747_453_596_211,
            ],
            [
                0.039_583_424_160_566,
                0.077_979_130_377_369,
                0.147_583_617_650_433,
                0.250_000_000_000_000,
                0.500_000_000_000_000,
                0.750_000_000_000_000,
                0.852_416_382_349_567,
                0.922_020_869_622_631,
                0.960_416_575_839_434,
            ],
            [
                0.007_634_036_082_669,
                0.028_595_479_208_968,
                0.091_751_709_536_137,
                0.211_324_865_405_187,
                0.500_000_000_000_000,
                0.788_675_134_594_813,
                0.908_248_290_463_863,
                0.971_404_520_791_032,
                0.992_365_963_917_331,
            ],
            [
                0.002_038_288_793_893,
                0.014_004_228_005_073,
                0.069_662_984_279_422,
                0.195_501_109_477_885,
                0.500_000_000_000_000,
                0.804_498_890_522_115,
                0.930_337_015_720_578,
                0.985_995_771_994_927,
                0.997_961_711_206_107,
            ],
            [
                0.000_661_948_454_609,
                0.008_065_044_950_046,
                0.058_058_261_758_408,
                0.186_950_483_150_030,
                0.500_000_000_000_000,
                0.813_049_516_849_970,
                0.941_941_738_241_592,
                0.991_934_955_049_954,
                0.999_338_051_545_391,
            ],
            [
                0.000_005_887_471_395,
                0.001_259_166_312_368,
                0.036_694_017_385_370,
                0.170_446_566_151_030,
                0.500_000_000_000_000,
                0.829_553_433_848_970,
                0.963_305_982_614_630,
                0.998_740_833_687_632,
                0.999_994_112_528_605,
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
                -102_849.115_630_175_55,
                -1_028.491_010_471_621_9,
                -10.270_324_410_234_506,
                0.0,
                10.270_324_410_234_506,
                1_028.491_010_471_621_9,
                102_849.115_630_175_55,
            ],
            [
                -318.308_838_985_550_45,
                -31.820_515_953_773_96,
                -3.077_683_537_175_253_4,
                0.0,
                3.077_683_537_175_253_4,
                31.820_515_953_773_96,
                318.308_838_985_550_45,
            ],
            [
                -22.327_124_770_119_875,
                -6.964_556_734_283_274,
                -1.885_618_083_164_126_7,
                0.0,
                1.885_618_083_164_126_7,
                6.964_556_734_283_274,
                22.327_124_770_119_875,
            ],
            [
                -10.214_531_852_407_387,
                -4.540_702_858_568_133_6,
                -1.637_744_353_696_210_1,
                0.0,
                1.637_744_353_696_210_1,
                4.540_702_858_568_133_6,
                10.214_531_852_407_387,
            ],
            [
                -7.173_182_219_782_308_5,
                -3.746_947_387_979_197,
                -1.533_206_274_058_943_9,
                0.0,
                1.533_206_274_058_943_9,
                3.746_947_387_979_197,
                7.173_182_219_782_308_5,
            ],
            [
                -4.143_700_494_046_59,
                -2.763_769_458_112_696,
                -1.372_183_641_110_335_6,
                0.0,
                1.372_183_641_110_335_6,
                2.763_769_458_112_696,
                4.143_700_494_046_59,
            ],
        ];
        for (i, &nu) in NU.iter().enumerate() {
            let x = P.map(|p| crate::stdtri(nu, p));
            crate::np_assert_allclose!(x, X_REF[i], rtol = 1e-13);
        }
    }

    #[test]
    fn test_ndtr_f64() {
        xsref::test("ndtr", "d-d", |x| crate::ndtr(x[0]));
    }

    #[test]
    fn test_ndtr_c64() {
        xsref::test("ndtr", "cd-cd", |x| crate::ndtr(c64(x[0], x[1])));
    }

    // Based on `scipy.special.tests.test_ndtr.TestLogNdtr.test_log_ndtr_moderate_le8`
    #[test]
    #[allow(clippy::approx_constant)]
    fn test_log_ndtr() {
        let x = [-0.75, -0.25, 0.0, 0.5, 1.5, 2.5, 3.0, 4.0, 5.0, 7.0, 8.0];
        let expected = [
            -1.484_448_229_919_656_2,
            -0.913_061_764_811_135_1,
            -0.693_147_180_559_945_3,
            -0.368_946_415_288_656_4,
            -0.069_143_455_612_233_98,
            -0.006_229_025_485_860_002,
            -0.001_350_809_964_748_193_8,
            -3.167_174_337_748_927e-05,
            -2.866_516_129_637_636e-07,
            -1.279_812_543_886_654e-12,
            -6.220_960_574_271_786e-16,
        ];
        let y = x.map(crate::log_ndtr);
        crate::np_assert_allclose!(&y, &expected, rtol = 1e-14);
    }

    #[test]
    fn test_ndtri() {
        xsref::test("ndtri", "d-d", |x| crate::ndtri(x[0]));
    }

    #[test]
    fn test_kolmogorov() {
        xsref::test("kolmogorov", "d-d", |x| crate::kolmogorov(x[0]));
    }

    #[test]
    fn test_kolmogc() {
        xsref::test("kolmogc", "d-d", |x| crate::kolmogc(x[0]));
    }

    #[test]
    fn test_kolmogi() {
        xsref::test("kolmogi", "d-d", |x| crate::kolmogi(x[0]));
    }

    #[test]
    fn test_kolmogci() {
        xsref::test("kolmogci", "d-d", |x| crate::kolmogci(x[0]));
    }

    #[test]
    fn test_kolmogp() {
        xsref::test("kolmogp", "d-d", |x| crate::kolmogp(x[0]));
    }

    #[test]
    fn test_smirnov() {
        xsref::test("smirnov", "p_d-d", |x| crate::smirnov(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovc() {
        xsref::test("smirnovc", "p_d-d", |x| crate::smirnovc(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovi() {
        xsref::test("smirnovi", "p_d-d", |x| crate::smirnovi(x[0] as i32, x[1]));
    }

    #[test]
    fn test_smirnovci() {
        xsref::test("smirnovci", "p_d-d", |x| {
            crate::smirnovci(x[0] as i32, x[1])
        });
    }

    #[test]
    fn test_smirnovp() {
        xsref::test("smirnovp", "p_d-d", |x| crate::smirnovp(x[0] as i32, x[1]));
    }

    #[test]
    fn test_owens_t() {
        xsref::test("owens_t", "d_d-d", |x| crate::owens_t(x[0], x[1]));
    }

    #[test]
    fn test_chdtr() {
        xsref::test("chdtr", "d_d-d", |x| crate::chdtr(x[0], x[1]));
    }

    #[test]
    fn test_chdtrc() {
        xsref::test("chdtrc", "d_d-d", |x| crate::chdtrc(x[0], x[1]));
    }

    #[test]
    fn test_chdtri() {
        xsref::test("chdtri", "d_d-d", |x| crate::chdtri(x[0], x[1]));
    }

    #[test]
    fn test_fdtr() {
        xsref::test("fdtr", "d_d_d-d", |x| crate::fdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtrc() {
        xsref::test("fdtrc", "d_d_d-d", |x| crate::fdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_fdtri() {
        xsref::test("fdtri", "d_d_d-d", |x| crate::fdtri(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtr() {
        xsref::test("gdtr", "d_d_d-d", |x| crate::gdtr(x[0], x[1], x[2]));
    }

    #[test]
    fn test_gdtrc() {
        xsref::test("gdtrc", "d_d_d-d", |x| crate::gdtrc(x[0], x[1], x[2]));
    }

    #[test]
    fn test_pdtr() {
        xsref::test("pdtr", "d_d-d", |x| crate::pdtr(x[0], x[1]));
    }

    #[test]
    fn test_pdtrc() {
        xsref::test("pdtrc", "d_d-d", |x| crate::pdtrc(x[0], x[1]));
    }

    #[test]
    fn test_pdtri() {
        xsref::test("pdtri", "p_d-d", |x| crate::pdtri(x[0] as i32, x[1]));
    }

    #[test]
    fn test_bdtr() {
        xsref::test("bdtr", "d_p_d-d", |x| crate::bdtr(x[0], x[1] as i32, x[2]));
    }

    #[test]
    fn test_bdtrc() {
        xsref::test("bdtrc", "d_p_d-d", |x| {
            crate::bdtrc(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_bdtri() {
        xsref::test("bdtri", "d_p_d-d", |x| {
            crate::bdtri(x[0], x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtr() {
        xsref::test("nbdtr", "p_p_d-d", |x| {
            crate::nbdtr(x[0] as i32, x[1] as i32, x[2])
        });
    }

    #[test]
    fn test_nbdtrc() {
        xsref::test("nbdtrc", "p_p_d-d", |x| {
            crate::nbdtrc(x[0] as i32, x[1] as i32, x[2])
        });
    }
}
