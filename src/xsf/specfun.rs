mod bessel;
mod erf;
mod gamma;
mod hyperg;
mod legendre;
mod mathieu;
mod numbers;
mod sphd;

pub(crate) use bessel::{jynb, jyzo, lamn, lamv, msta1, msta2, sphj, sphy};
pub(crate) use erf::cerzo;
pub(crate) use gamma::{cgama, gamma2, psi_spec};
pub(crate) use hyperg::{cchg, chgm, chgu};
pub(crate) use legendre::{lpmns, lpmv, lqmns};
pub(crate) use mathieu::{cva2, fcoef, mtu0, mtu12};
pub(crate) use numbers::{bernob, eulerb};
pub(crate) use sphd::{aswfa, rswfo, rswfp, segv};

/// `specfun::Status`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    Ok,
    NoMemory,
    #[allow(dead_code)]
    Other,
}

/// `new (std::nothrow) T[n]()`, i.e. `None` if the allocation fails
pub(crate) fn try_zeroed(n: usize) -> Option<Vec<f64>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).ok()?;
    v.resize(n, 0.0);
    Some(v)
}

use num_complex::Complex;

pub trait Hyp1F1Arg: crate::sealed::Sealed {
    fn hyp1f1(self, a: f64, b: f64) -> Self;
}

impl Hyp1F1Arg for f64 {
    #[inline]
    fn hyp1f1(self, a: f64, b: f64) -> Self {
        // `xsf::hyp1f1`
        chgm(self, a, b)
    }
}

impl Hyp1F1Arg for Complex<f64> {
    #[inline]
    fn hyp1f1(self, a: f64, b: f64) -> Self {
        // `xsf::hyp1f1`
        if a.is_nan() || b.is_nan() || self.re.is_nan() || self.im.is_nan() {
            return Complex::new(f64::NAN, f64::NAN);
        }
        let mut outz = cchg(a, b, self);
        #[allow(clippy::float_cmp)]
        if outz.re == 1e300 {
            outz.re = f64::INFINITY;
        }
        outz
    }
}

/// Bernoulli numbers $B_0, \ldots, B_{N-1}$
///
/// Corresponds to [`scipy.special.bernoulli`][scipy] in SciPy, and calls the FFI function
/// `xsf::specfun::bernob`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bernoulli.html
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
/// # Panics
/// - If `N` is so large that `N - 1` overflows `i32`.
///
/// # See also
/// - [`euler`]: Euler numbers E<sub>0</sub>, ..., E<sub>N-1</sub>
#[must_use]
#[inline]
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
        bernob((N - 1).try_into().unwrap(), &mut bn);
    }
    bn
}

/// Euler numbers $E_0, \ldots, E_{N-1}$
///
/// Corresponds to [`scipy.special.euler`][scipy] in SciPy, and calls the FFI function
/// `xsf::specfun::eulerb`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.euler.html
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
/// # Panics
/// - If `N` is so large that `N - 1` overflows `i32`.
///
/// # See also
/// - [`bernoulli`]: Bernoulli numbers B<sub>0</sub>, ..., B<sub>N-1</sub>
#[must_use]
#[inline]
pub fn euler<const N: usize>() -> [f64; N] {
    let mut en = [0.0; N];
    if N < 3 {
        // `eulerb` requires `N >= 3` (i.e. `n >= 2`), so handle small `N` manually
        if N >= 1 {
            en[0] = 1.0;
        }
    } else {
        eulerb((N - 1).try_into().unwrap(), &mut en);
    }
    en
}

/// Kummer's Confluent hypergeometric function $_1F_1$
///
/// Corresponds to [`scipy.special.hyp1f1`][scipy] in SciPy, and accepts both `f64` and
/// `num_complex::Complex<f64>` inputs for `z`.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp1f1.html
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
/// - [`hypu`]: Tricomi's confluent hypergeometric function $U(a,b,x)$
/// - [`hyp0f1`](crate::hyp0f1): Confluent hypergeometric limit function,
///   $_0F_1\left[b\middle\| z\right]$
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
#[inline]
pub fn hyp1f1<T: Hyp1F1Arg>(a: f64, b: f64, z: T) -> T {
    z.hyp1f1(a, b)
}

/// `xsf::hypu`
fn xsf_hypu(a: f64, b: f64, x: f64) -> f64 {
    let mut out;
    let mut md = 0; /* method code --- not returned */
    let mut isfer = 0;

    out = chgu(x, a, b, &mut md, &mut isfer);
    #[allow(clippy::float_cmp)]
    if out == 1e300 {
        out = f64::INFINITY;
    }
    if isfer != 0 {
        out = f64::NAN;
    }
    out
}

/// Tricomi's confluent hypergeometric function $U(a,b,x)$
///
/// Corresponds to [`scipy.special.hyperu`][scipy] in SciPy.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyperu.html
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
/// - [`hyp1f1`]: Kummer's confluent hypergeometric function $M(a,b,z)$
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/13.2#E6>
#[doc(alias = "hyperu")]
#[must_use]
#[allow(clippy::float_cmp)]
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
            crate::pow_rising(1.0 - b + a, -a)
        }
    } else if x < 1.0 && a > -0.25 && a < 0.3 && b == 1.0 {
        // DLMF 13.3.7. Fixes scipy/scipy#15650
        let a1 = a + 1.0;
        (x + a + a1) * xsf_hypu(a1, 1.0, x) - a1 * a1 * xsf_hypu(a + 2.0, 1.0, x)
    } else {
        xsf_hypu(a, b, x)
    }
}

/// Associated Legendre function for $|x| \leq 1$
///
/// Corresponds to [`scipy.special.lpmv`][scipy] in SciPy.
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lpmv.html
///
/// # See also
/// - [`assoc_legendre_p_all`](crate::assoc_legendre_p_all)
#[must_use]
#[inline]
pub fn pmv(m: i32, v: f64, x: f64) -> f64 {
    // `xsf::pmv`
    let mut out = lpmv(x, m, v);
    #[allow(clippy::float_cmp)]
    if out == 1e300 {
        out = f64::INFINITY;
    }
    #[allow(clippy::float_cmp)]
    if out == -1e300 {
        out = f64::NEG_INFINITY;
    }
    out
}

#[cfg(test)]
mod tests {
    use num_complex::c64;
    use num_traits::ToPrimitive;

    #[test]
    fn test_hypu() {
        // the table is called "hyperu" instead of "hypu"
        xsref::test("hyperu", "d_d_d-d", |x| crate::hypu(x[0], x[1], x[2]));
    }

    #[test]
    fn test_hypu_tiny_a() {
        // NaN in xsf, because of an undefined `int` conversion
        assert_eq!(crate::hypu(1e-17, 0.5, 100.0), 1.0);
    }

    #[test]
    fn test_hypu_huge_a() {
        // NaN without ~2^31 iterations, like in xsf
        assert!(crate::hypu(-2_147_483_649.0, 1.0, 1.0).is_nan());
        assert!(crate::hypu(1.0, f64::INFINITY, 1.0).is_nan());
    }

    #[test]
    fn test_hyp1f1() {
        xsref::test("hyp1f1", "d_d_cd-cd", |x| {
            crate::hyp1f1(x[0], x[1], c64(x[2], x[3]))
        });
    }

    #[test]
    fn test_pmv() {
        xsref::test("pmv", "d_d_d-d", |x| {
            crate::pmv(x[0].to_i32().unwrap(), x[1], x[2])
        });
    }

    #[test]
    fn test_pmv_large_m() {
        // NaN in xsf, which overflows `-m`
        assert!(crate::pmv(i32::MIN, 0.5, 0.5).is_nan());
        // no `int` overflow of `j * j` in `lpmv0`
        let _ = crate::pmv(50_000, 0.5, 0.5);
    }

    #[test]
    fn test_hyp1f1_huge_a() {
        // `(int)a` is INT_MIN like in xsf, instead of ~2^31 iterations
        let _ = crate::hyp1f1(3e9, 1.0, 1.0);
        let _ = crate::hyp1f1(3e9, 1.0, c64(1.0, 0.0));
    }
}
