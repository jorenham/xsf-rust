use core::f64::consts::TAU;
use num_complex::Complex64;

use crate::utils::MAXLOG;

const LN_MAX: f64 = MAXLOG; // f64::MAX.ln()
const LN_MIN: f64 = -708.396_418_532_264_1; // f64::MIN_POSITIVE.ln()

/// Asymptotic expansion for $I_{v-1}(2\sqrt{z}) \Gamma(v)$ for real $z>0$ and $v \to +\infty$.
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_asy`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L59-L100>
///
/// Unlike in xsf (scipy/xsf#366), this also handles $z<0$ with $4|z| < (v-1)^2$ through Debye's expansion of
/// $J_{v-1}(2\sqrt{-z})$ (DLMF 10.19.3), i.e. with $\sqrt{1-x^2}$ instead of $\sqrt{1+x^2}$.
#[inline]
fn hyp0f1_asy(v: f64, z: f64) -> f64 {
    let arg = z.abs().sqrt();
    let v1 = (v - 1.0).abs();
    let x = 2.0 * arg / v1;
    let p1 = (1.0 + z.signum() * x * x).sqrt();
    let eta = p1 + x.ln() - crate::log1p(p1);

    let mut arg_exp_i = -0.5 * p1.ln();
    arg_exp_i -= 0.5 * (TAU * v1).ln();
    arg_exp_i += crate::gammaln(v);
    let gs = crate::gammasgn(v);

    let mut arg_exp_k = arg_exp_i;
    arg_exp_i += v1 * eta;
    arg_exp_k -= v1 * eta;

    // large-v asymptotic correction, DLMF 10.41.10
    let pp = 1.0 / p1;
    let p2 = pp * pp;
    let p4 = p2 * p2;
    let p6 = p4 * p2;
    let u1 = (3.0 - 5.0 * p2) * pp / 24.0;
    let u2 = (81.0 - 462.0 * p2 + 385.0 * p4) * p2 / 1152.0;
    let u3 = (30375.0 - 369_603.0 * p2 + 765_765.0 * p4 - 425_425.0 * p6) * pp * p2 / 414_720.0;
    let u_corr_i = 1.0 + u1 / v1 + u2 / (v1 * v1) + u3 / (v1 * v1 * v1);

    let out = (arg_exp_i - crate::xlogy(v1, arg)).exp() * gs * u_corr_i;

    if v - 1.0 < 0.0 {
        // DLMF 10.27.2: I_{-v} = I_{v} + (2/pi) sin(pi*v) K_v
        // (or DLMF 10.2.3: J_{-v} = cos(pi*v) J_v - sin(pi*v) Y_v, where J_v is negligible if
        // J_{-v} overflows)
        let out = if z < 0.0 { 0.0 } else { out };
        let u_corr_k = 1.0 - u1 / v1 + u2 / (v1 * v1) - u3 / (v1 * v1 * v1);
        out + (arg_exp_k + crate::xlogy(v1, arg)).exp() * gs * 2.0 * crate::sinpi(v1) * u_corr_k
    } else {
        out
    }
}

/// Real-valued kernel
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_real`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L31-L56>
#[inline]
#[allow(clippy::float_cmp)]
fn hyp0f1_real(v: f64, z: f64) -> f64 {
    if v.is_nan() || z.is_nan() {
        f64::NAN
    } else if v <= 0.0 && v == v.floor() {
        // poles
        f64::NAN
    } else if z == 0.0 && v != 0.0 {
        // zeros
        1.0
    } else if z.abs() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)
        1.0 + z / v + z * z / (2.0 * v * (v + 1.0))
    } else if z > 0.0 {
        let arg = z.sqrt();
        let arg_exp = crate::xlogy(1.0 - v, arg) + crate::gammaln(v);
        let bess_val = crate::xsf::cephes::iv(v - 1.0, 2.0 * arg);

        if arg_exp > LN_MAX || bess_val == 0.0 || arg_exp < LN_MIN || bess_val.is_infinite() {
            // overflow or underflow
            hyp0f1_large_v(v, z)
        } else {
            arg_exp.exp() * bess_val * crate::gammasgn(v)
        }
    } else {
        let arg = (-z).sqrt();
        let arg_pow = arg.powf(1.0 - v);
        let gamma_v = crate::gamma(v);
        let prefactor = arg_pow * gamma_v;
        let bess_val = crate::bessel_j(v - 1.0, 2.0 * arg);
        let out = prefactor * bess_val;
        let prefactor_ok = arg_pow.is_normal() && gamma_v.is_normal() && prefactor.is_normal();

        // Unlike in xsf (scipy/xsf#366), fall back if a factor over- or underflows (for |v| >~ 100), e.g.
        // hyp0f1(200, -1) would be NaN instead of 0.995 (Gamma(200) = inf and J_199(2) = 0)
        if prefactor_ok && out.is_normal() {
            out
        } else if -z > v.abs() && bess_val.is_normal() {
            (crate::xlogy(1.0 - v, arg) + crate::gammaln(v)).exp() * bess_val * crate::gammasgn(v)
        } else if -z <= v.abs()
            || (4.0 * -z < (v - 1.0).powi(2) && (v > 1.0 || !bess_val.is_finite()))
        {
            // the series, or Debye's expansion if J_{v-1}(2 sqrt(-z)) under- (or for v < 1,
            // over-) flowed in its non-oscillatory regime
            hyp0f1_large_v(v, z)
        } else {
            out
        }
    }
}

/// Fallback for large $|v|$, if the evaluation in terms of Bessel functions over- or underflows
///
/// Unlike in xsf (scipy/xsf#366), the Maclaurin series is summed if $|z| \le |v|$, since `hyp0f1_asy` loses roughly
/// $\epsilon |v| \ln|v|$ in its exponent, e.g. `hyp0f1(1.3e14, 2.3e8)` would be $e$ instead of
/// 1.0000018.
fn hyp0f1_large_v(v: f64, z: f64) -> f64 {
    if z.abs() <= v.abs() {
        hyp0f1_series(v, z)
    } else {
        hyp0f1_asy(v, z)
    }
}

/// Maclaurin series $\sum_n z^n / ((v)_n n!)$ (DLMF 16.2.1)
///
/// Only used for $|z| \le |v|$ with large $|v|$, where the terms decrease rapidly (up to the
/// negligible ones with $n > -v$ if $v < 0$), and cancel by at most roughly $e^2$.
fn hyp0f1_series(v: f64, z: f64) -> f64 {
    let mut term = 1.0;
    let mut sum = 1.0;
    for n in 0..1000 {
        let n = f64::from(n);
        term *= z / (v + n) / (n + 1.0);
        sum += term;
        if term.abs() <= f64::EPSILON * sum.abs() {
            break;
        }
    }
    sum
}

/// Complex valued kernel
///
/// Corresponds to `scipy.special._hyp0f1._hyp0f1_cmplx`:
/// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/_hyp0f1.pxd#L106-L136>
#[inline]
#[allow(clippy::float_cmp)]
fn hyp0f1_cmplx(v: f64, z: Complex64) -> Complex64 {
    if v.is_nan() || z.re.is_nan() || z.im.is_nan() {
        Complex64::new(f64::NAN, f64::NAN)
    } else if v <= 0.0 && v == v.floor() {
        // poles
        f64::NAN.into()
    } else if z.re == 0.0 && z.im == 0.0 && v != 0.0 {
        // zeros
        1.0.into()
    } else if z.norm() < 1e-6 * (1.0 + v.abs()) {
        // both v and z small: truncate the Taylor series at O(z**2)

        // need to do computations in this order, for otherwise $v\approx -z \ll 1$
        // it can lose precision (as was reported for 32-bit linux, see scipy/scipy#6365)
        let t1 = 1.0 + z / v;
        let t2 = z * z / (2.0 * v * (v + 1.0));
        t1 + t2
    } else {
        let (arg, r) = if z.re > 0.0 {
            let arg = z.sqrt();
            let bess_val = crate::bessel_i(v - 1.0, 2.0 * arg);
            (arg, bess_val)
        } else {
            let arg = (-z).sqrt();
            let bess_val = crate::bessel_j(v - 1.0, 2.0 * arg);
            (arg, bess_val)
        };

        r * crate::gamma(v) * arg.powf(1.0 - v)
    }
}

pub trait HypergeometricArg: crate::sealed::Sealed {
    fn hyp0f0(self) -> Self;
    fn hyp1f0(self, a: f64) -> Self;
    fn hyp0f1(self, b: f64) -> Self;
}

impl HypergeometricArg for f64 {
    #[inline]
    fn hyp0f0(self) -> Self {
        self.exp()
    }

    #[inline]
    #[allow(clippy::float_cmp)]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self == 1.0 {
            f64::NAN
        } else {
            (1.0 - self).powf(-a)
        }
    }

    #[inline]
    fn hyp0f1(self, b: f64) -> Self {
        hyp0f1_real(b, self)
    }
}

impl HypergeometricArg for num_complex::Complex<f64> {
    #[inline]
    fn hyp0f0(self) -> Self {
        self.exp()
    }

    #[inline]
    #[allow(clippy::float_cmp)]
    fn hyp1f0(self, a: f64) -> Self {
        if a > 0.0 && self.re == 1.0 && self.im == 0.0 {
            f64::NAN.into()
        } else {
            (1.0 - self).powf(-a)
        }
    }

    #[inline]
    fn hyp0f1(self, b: f64) -> Self {
        hyp0f1_cmplx(b, self)
    }
}

/// Hypergeometric function $_0F_0\left[\middle\| z\right]$ for real or complex $z$
///
/// This function is currently not implemented in SciPy, but is straightforward to evaluate.
/// Both [`f64`] and [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _0F_0\left[ \middle\| z\right] = \sum\_{n=0}^\infty {z^n \over n!},
/// $$
///
/// which is just the exponential function, $e^z$.
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function $_1F_0\left[a\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
#[inline]
pub fn hyp0f0<T: HypergeometricArg>(z: T) -> T {
    z.hyp0f0()
}

/// Hypergeometric function $_1F_0\left[a\middle\| z\right]$ for real or complex $z$
///
/// This function is currently not implemented in SciPy, but is straightforward to evaluate.
/// Both [`f64`] and [`num_complex::Complex<f64>`](num_complex::Complex) are accepted for `z`.
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _1F_0\left[a\middle\| z\right] = \sum\_{n=0}^\infty \rpow a n {z^n \over n!}
/// $$
///
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
/// It's also equal to $(1 - z)^{-a}$.
///
/// # See also
/// - [`hyp0f0`]: Hypergeometric function $_0F_0\left[\middle\| z\right]$
/// - [`hyp0f1`]: Confluent hypergeometric limit function, $_0F_1\left[b\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
#[inline]
pub fn hyp1f0<T: HypergeometricArg>(a: f64, z: T) -> T {
    z.hyp1f0(a)
}

/// Confluent hypergeometric limit function $_0F_1\left[b\middle\| z\right]$ for real or complex $z$
///
/// This is a translation of the [`scipy.special.hyp0f1`][hyp0f1] Cython implementation into Rust,
/// as in [scipy/xsf][xsf], which uses AMOS (instead of Cephes) for $J_v$.
///
/// [hyp0f1]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp0f1.html
/// [xsf]: https://github.com/scipy/xsf
///
/// # Notes
///
/// This function is defined as
///
/// $$
/// _0F_1\left[b\middle\| z\right] = \sum\_{n=0}^\infty {1\over \rpow b n} {z^n \over n!}.
/// $$
///
/// Here $\rpow \square n$ is the rising factorial; see [`pow_rising`](crate::pow_rising).
///
/// It's also the limit as $a \to \infty$ of $\hyp 1 1 a b z$, and satisfies the
/// differential equation $f\'\'(z) + b f\'(z) = f(z)$. See [^1] for more information.
///
/// # See also
/// - [`hyp1f0`]: Hypergeometric function, $_1F_0\left[a\middle\| z\right]$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
///
/// [^1]: Weisstein, Eric W. "Confluent Hypergeometric Limit Function." From MathWorld -- A Wolfram
/// Resource. <https://mathworld.wolfram.com/ConfluentHypergeometricLimitFunction.html>
///
#[inline]
pub fn hyp0f1<T: HypergeometricArg>(b: f64, z: T) -> T {
    z.hyp0f1(b)
}

#[cfg(test)]
mod tests {
    use crate::{hyp0f1, np_assert_allclose};
    use core::f64::consts::{E, PI};
    use num_complex::c64;

    #[test]
    fn test_hyp0f0_f64() {
        np_assert_allclose!([crate::hyp0f0(0.0)], [1.0], atol = 1e-15);
        np_assert_allclose!([crate::hyp0f0(1.0)], [E], atol = 1e-15);
        np_assert_allclose!([crate::hyp0f0(-1.0)], [1.0 / E], atol = 1e-15);
    }

    #[test]
    fn test_hyp0f0_c64() {
        np_assert_allclose!(
            [crate::hyp0f0(c64(0.0, 0.0))],
            [c64(1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!([crate::hyp0f0(c64(1.0, 0.0))], [c64(E, 0.0)], atol = 1e-15);
        np_assert_allclose!(
            [crate::hyp0f0(c64(0.0, PI))],
            [c64(-1.0, 0.0)],
            atol = 1e-15
        );
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_hyp1f0_f64() {
        assert_eq!(crate::hyp1f0(-1.0, -1.0), 2.0);
        assert_eq!(crate::hyp1f0(0.0, -1.0), 1.0);
        assert_eq!(crate::hyp1f0(1.0, -1.0), 0.5);
        assert_eq!(crate::hyp1f0(2.0, -1.0), 0.25);

        assert_eq!(crate::hyp1f0(-1.0, 0.0), 1.0);
        assert_eq!(crate::hyp1f0(0.0, 0.0), 1.0);
        assert_eq!(crate::hyp1f0(1.0, 0.0), 1.0);
        assert_eq!(crate::hyp1f0(2.0, 0.0), 1.0);

        assert_eq!(crate::hyp1f0(-1.0, 1.0), 0.0);
        assert_eq!(crate::hyp1f0(0.0, 1.0), 1.0);
        assert!(crate::hyp1f0(1.0, 1.0).is_nan());
        assert!(crate::hyp1f0(2.0, 1.0).is_nan());

        assert_eq!(crate::hyp1f0(-1.0, 2.0), -1.0);
        assert_eq!(crate::hyp1f0(0.0, 2.0), 1.0);
        assert_eq!(crate::hyp1f0(1.0, 2.0), -1.0);
        assert_eq!(crate::hyp1f0(2.0, 2.0), 1.0);
    }

    #[test]
    fn test_hyp1f0_c64() {
        assert_eq!(crate::hyp1f0(-1.0, c64(-1.0, 0.0)), c64(2.0, 0.0));
        assert_eq!(crate::hyp1f0(0.0, c64(-1.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(1.0, c64(-1.0, 0.0)), c64(0.5, 0.0));
        assert_eq!(crate::hyp1f0(2.0, c64(-1.0, 0.0)), c64(0.25, 0.0));

        assert_eq!(crate::hyp1f0(-1.0, c64(0.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(0.0, c64(0.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(1.0, c64(0.0, 0.0)), c64(1.0, 0.0));
        assert_eq!(crate::hyp1f0(2.0, c64(0.0, 0.0)), c64(1.0, 0.0));

        assert_eq!(crate::hyp1f0(-1.0, c64(1.0, 0.0)), c64(0.0, 0.0));
        assert_eq!(crate::hyp1f0(0.0, c64(1.0, 0.0)), c64(1.0, 0.0));
        assert!(crate::hyp1f0(1.0, c64(1.0, 0.0)).is_nan());
        assert!(crate::hyp1f0(2.0, c64(1.0, 0.0)).is_nan());

        np_assert_allclose!(
            [crate::hyp1f0(-1.0, c64(2.0, 0.0))],
            [c64(-1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!(
            [crate::hyp1f0(0.0, c64(2.0, 0.0))],
            [c64(1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!(
            [crate::hyp1f0(1.0, c64(2.0, 0.0))],
            [c64(-1.0, 0.0)],
            atol = 1e-15
        );
        np_assert_allclose!(
            [crate::hyp1f0(2.0, c64(2.0, 0.0))],
            [c64(1.0, 0.0)],
            atol = 1e-15
        );
    }

    /// Corresponds to `scipy.special.tests.test_basic.TestHyper.test_hyp0f1`
    /// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/tests/test_basic.py#L3144-L3168>
    #[test]
    fn test_hyp0f1() {
        // scalar input
        np_assert_allclose!([hyp0f1(2.5, 0.5)], [1.214_827_026_899_997], rtol = 1e-12);
        np_assert_allclose!([hyp0f1(2.5, 0.0)], [1.0], rtol = 1e-15);

        // float input, expected values match mpmath
        let x = [-1.5, -1.0, 0.0, 1.0, 1.5].map(|z| hyp0f1(3.0, z));
        let expected = [
            0.584_936_592_291_43,
            0.705_668_057_231_27,
            1.0,
            1.377_896_895_397_47,
            1.603_736_852_884_80,
        ];
        np_assert_allclose!(x, expected, rtol = 1e-12);

        // complex input
        let x = [-1.5, -1.0, 0.0, 1.0, 1.5].map(|re| hyp0f1(3.0, c64(re, 0.0)));
        np_assert_allclose!(x, expected.map(|v| c64(v, 0.0)), rtol = 1e-12);
    }

    /// Corresponds to `scipy.special.tests.test_basic.TestHyper.test_hyp0f1_gh5764`
    /// <https://github.com/scipy/scipy/blob/55053f8/scipy/special/tests/test_basic.py#L3170-L3176>
    /// which looks like a regression test for scipy/scipy#5764
    #[test]
    fn test_hyp0f1_gh5764() {
        // Just checks the point that failed; there's a more systematic test in test_mpmath
        let res = hyp0f1(0.8, c64(0.5, 0.5));
        // The expected value was generated using mpmath
        np_assert_allclose!(
            [res],
            [c64(1.613_971_977_644_111_5, 0.808_930_540_617_906_7)],
            atol = 1.5e-7
        );
    }

    #[test]
    fn test_hyp0f1_large_b() {
        let bz = [
            // large-v asymptotic expansion
            (300.0, 1.0),
            (-661.112_997_069_224_4, 0.062_385_007_293_382_27),
            // b near a non-positive integer
            (-9.000_000_000_000_306, -0.019_499_612_460_819_94),
            // over- or underflow of Gamma(b), |z|^((1-b)/2), or J_{b-1}(2 sqrt(-z))
            (200.0, -1.0),
            (233.063_303_394_596_46, -511.589_207_216_921_75),
            (873.979_261_659_223_8, -1_963.183_700_462_148_6),
            (-628.010_908_208_124_6, -2_515.795_365_646_08),
            (130_356_679_484_338.69, 229_453_612.840_189_75),
            (-127.999_999_999_999_96, -0.023),
        ];
        // computed with mpmath
        let expected = [
            1.003_338_876_548_527_4,
            0.999_905_640_851_462_9,
            1.002_169_266_580_656_6,
            0.995_012_417_311_756_4,
            0.110_196_664_254_785_33,
            0.105_488_627_754_553_3,
            55.638_844_554_215_72,
            1.000_001_760_199_904_5,
            1.000_179_703_771_905_4,
        ];
        np_assert_allclose!(bz.map(|(b, z)| hyp0f1(b, z)), expected, rtol = 1e-12);

        // near zeros of J_{b-1}(2 sqrt(-z)), where it rounds to 0 (instead of 5.2e-17 and 6.1e-15)
        assert!(hyp0f1(1.674_445_551_227_833_6, -2.866_969_344_179_029).abs() < 1e-14);
        assert!(hyp0f1(-3.854_390_408_608_444_4, -5.550_754_014_064_338).abs() < 1e-14);
    }

    #[test]
    fn test_hyp0f1_nan() {
        assert!(hyp0f1(f64::NAN, 0.0).is_nan());
        assert!(hyp0f1(1.0, c64(f64::NAN, 0.0)).is_nan());
        // pole at v = -inf
        assert!(hyp0f1(f64::NEG_INFINITY, 0.0).is_nan());
        assert!(hyp0f1(f64::NEG_INFINITY, c64(0.0, 0.0)).is_nan());
    }
}
