//! Translated from <https://github.com/scipy/scipy/blob/4845755/scipy/special/_convex_analysis.pxd>

/// Elementwise function for computing the entropy
///
/// $$
/// \mathrm{entr}(x_i) = \begin{cases}
///     -x_i \ln(x_i) \& x_i > 0 \\\\
///     0             \& x_i = 0 \\\\
///     -\infty       \& \text{otherwise}
/// \end{cases}
/// $$
///
/// Translated from [`scipy.special.entr`][entr] into pure Rust.
///
/// [entr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.entr.html
///
/// # Notes
///
/// This function is concave.
///
/// The origin of this function is in convex programming [^1]. Given a discrete probability
/// distribution $p_1,\dotsc,p_n$, the definition of entropy in the context of
/// *information theory* is
///
/// $$
/// H\[ \\{ p_i \\} \] = \sum\_{i} \mathrm{entr}(p_i).
/// $$
///
/// [^1]: S. Boyd and L. Vandenberghe. *Convex optimization*. Cambridge University Press, 2004.
///       <https://doi.org/10.1017/CBO9780511804441>
///
/// # See also
/// - [`rel_entr`]
/// - [`kl_div`]
#[must_use]
#[inline]
pub fn entr(x_i: f64) -> f64 {
    if x_i.is_nan() {
        f64::NAN
    } else if x_i > 0.0 {
        -x_i * x_i.ln()
    } else if x_i == 0.0 {
        0.0
    } else {
        f64::NEG_INFINITY
    }
}

/// Elementwise function for computing relative entropy
///
/// Translated from [`scipy.special.rel_entr`][rel_entr] into pure Rust.
///
/// $$
/// \mathrm{rel\\_entr}(x_i,\\, y_i) = \begin{cases}
///     x_i \ln(x_i / y_i) \& x_i > 0 \wedge y_i > 0 \\\\
///     0                  \& x_i = 0 \wedge y_i \ge 0 \\\\
///     \infty             \& \text{otherwise}
/// \end{cases}
/// $$
///
/// [rel_entr]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html
///
/// # Notes
///
/// This function is jointly convex in $x$ and $y$.
///
/// The origin of this function is in convex programming [^1]. Given two discrete probability
/// distributions $p_1,\dotsc,p_n$ and $q_1,\dotsc,q_n$, the definition of relative entropy in the
/// context of *information theory* is [^2]
///
/// $$
/// H\left[ \\{ p_i \\} \\, \rvert \\, \\{ q_i \\} \right]
///     = \sum_i \mathrm{rel\\_entr}(p_i,\\, q_i).
/// $$
///
/// [^1]: S. Boyd and L. Vandenberghe. *Convex optimization*. Cambridge University Press, 2004.
///       <https://doi.org/10.1017/CBO9780511804441>
/// [^2]: Kullback-Leibler divergence, <https://wikipedia.org/wiki/KL-divergence>
///
/// # See also
/// - [`entr`]
/// - [`kl_div`]
#[must_use]
#[inline]
pub fn rel_entr(x_i: f64, y_i: f64) -> f64 {
    if x_i.is_nan() || y_i.is_nan() {
        f64::NAN
    } else if x_i <= 0.0 || y_i <= 0.0 {
        if x_i == 0.0 && 0.0 <= y_i {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        let ratio = x_i / y_i;
        if (0.5..2.0).contains(&ratio) {
            // When x_i and y_i are close, this is more accurate
            x_i * ((x_i - y_i) / y_i).ln_1p()
        } else if (f64::MIN_POSITIVE..f64::INFINITY).contains(&ratio) {
            // There are no underflow or overflow issues
            x_i * ratio.ln()
        } else {
            // x_i and y_i are so far apart that taking x_i / y_i results in either an underflow,
            // overflow, or subnormal number -- do the logarithm first
            x_i * (x_i.ln() - y_i.ln())
        }
    }
}

/// Elementwise function for computing the Kullback-Leibler divergence
///
/// $$
/// \mathrm{kl\\_div}(x_i,\\, y_i) = \begin{cases}
///     x_i \ln(x_i / y_i) - x_i + y_i  \& x_i > 0 \wedge y_i > 0 \\\\
///     y_i                             \& x_i = 0 \wedge y_i \ge 0 \\\\
///     \infty                          \& \text{otherwise}
/// \end{cases}
/// $$
///
/// Translated from [`scipy.special.kl_div`][kl_div] into pure Rust.
///
/// [kl_div]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kl_div.html
///
/// # Notes
///
/// This function is non-negative and is jointly convex in $x$ and $y$.
///
/// The origin of this function is in convex programming [^1]. This is why the function contains
/// the extra $-x_i + y_i$ terms over what might be expected from the Kullback-Leibler divergence.
/// For a version of the function without the extra terms, see [`rel_entr`].
///
/// [^1]: S. Boyd and L. Vandenberghe. *Convex optimization*. Cambridge University Press, 2004.
///       <https://doi.org/10.1017/CBO9780511804441>
///
/// # See also
/// - [`entr`]
/// - [`rel_entr`]
#[must_use]
#[inline]
pub fn kl_div(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        f64::NAN
    } else if 0.0 < x && 0.0 < y {
        // Unlike in SciPy (scipy/xsf#362), `rel_entr` is reused, which takes the logarithm first if
        // x / y under- or overflows, e.g. kl_div(1e-184, 1e258) was -inf instead of 1e258.
        //
        // For x / y in [0.5, 2), with x = y (1 + e), it's y ((1 + e) ln(1 + e) - e), written so
        // that nothing cancels, e.g. kl_div(1 + 1e-12, 1) had a relative error of 28.
        let e = (x - y) / y;
        if (-0.5..1.0).contains(&e) {
            y * (crate::log1pmx(e) + e * e.ln_1p())
        } else {
            rel_entr(x, y) - x + y
        }
    } else if x == 0.0 && 0.0 <= y {
        y
    } else {
        f64::INFINITY
    }
}

/// Huber Loss function
///
/// $$
/// L_\delta(r) = \begin{cases}
///     \infty                                          & \delta < 0  \\\\
///     \frac{1}{2} r^2                                 & \|r\| \le \delta \\\\
///     \delta \left(\|r\| - \frac{1}{2} \delta \right) & \text{otherwise}
/// \end{cases}
/// $$
///
/// Translated from [`scipy.special.huber`][huber] into pure Rust.
///
/// [huber]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.huber.html
///
/// # Notes
///
/// `huber` is useful as a loss function in robust statistics or machine learning to reduce the
/// influence of outliers as compared to the common squared error loss, residuals with a magnitude
/// higher than $\delta$ are not squared [^1].
///
/// Typically, $r$ represents residuals, the difference between a model prediction and data.
/// Then, for $\|r\| \leq \delta$, `huber` resembles the squared error and for $\|r\| > \delta$ the
/// absolute error. This way, the Huber loss often achieves a fast convergence in model fitting for
/// small residuals like the squared error loss function and still reduces the influence of
/// outliers ($\|r\| > \delta$) like the absolute error loss. As $\delta$ is the cutoff between
/// squared and absolute error regimes, it has to be tuned carefully for each problem. `huber` is
/// also convex, making it suitable for gradient based optimization.
///
/// [^1]: Peter Huber, "Robust Estimation of a Location Parameter".
///       1964. Annals of Statistics. 53 (1): 73 - 101
///
/// # See also
/// - [`pseudo_huber`]
#[must_use]
#[inline]
pub fn huber(delta: f64, r: f64) -> f64 {
    if delta < 0.0 {
        f64::INFINITY
    } else {
        let r_abs = r.abs();
        if r_abs <= delta {
            0.5 * r * r
        } else {
            delta * (r_abs - 0.5 * delta)
        }
    }
}

/// Pseudo-Huber Loss function
///
/// $$
/// \widetilde{L}_\delta(r) =
///     \delta^2 \left( \sqrt{ 1 + \left( \frac{r}{\delta} \right)^2 } - 1 \right)
/// $$
///
/// Translated from [`scipy.special.pseudo_huber`][phuber] into pure Rust.
///
/// [phuber]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pseudo_huber.html
///
/// # Notes
///
/// Like [`huber`], `pseudo_huber` often serves as a robust loss function in statistics or machine
/// learning to reduce the influence of outliers. Unlike [`huber`], `pseudo_huber` is smooth.
///
/// Typically, $r$ represents residuals, the difference between a model prediction and data.
/// Then, for $\|r\| \leq \delta$, `pseudo_huber` resembles the squared error and for
/// $\|r\| > \delta$ the absolute error. This way, the Pseudo-Huber loss often achieves a fast
/// convergence in model fitting for small residuals like the squared error loss function and still
/// reduces the influence of outliers ($\|r\| > \delta$) like the absolute error loss. As $\delta$
/// is the cutoff between squared and absolute error regimes, it has to be tuned carefully for each
/// problem. `pseudo_huber` is also convex, making it suitable for gradient based optimization
/// [^1] [^2].
///
/// [^1]: Hartley, Zisserman, "Multiple View Geometry in Computer Vision".
///       2003. Cambridge University Press. p. 619
/// [^2]: Charbonnier et al. "Deterministic edge-preserving regularization in computed imaging".
///       1997. IEEE Trans. Image Processing. 6 (2): 298 - 311
///
/// # See also
/// - [`huber`]
#[must_use]
#[inline]
pub fn pseudo_huber(delta: f64, r: f64) -> f64 {
    if delta < 0.0 {
        f64::INFINITY
    } else if delta == 0.0 || r == 0.0 {
        0.0
    } else {
        // Unlike in SciPy (scipy/xsf#362), u^2 (sqrt(1 + v^2) - 1) with u = delta and v = r / u is
        // not computed as u^2 expm1(log1p(v^2) / 2), where u^2 or v^2 over- or underflows if u or
        // |v| is huge or tiny, and the log-exp round trip loses precision if |v| is large. Instead,
        // the equivalent r^2 / (sqrt(1 + v^2) + 1) is used, scaled like `hypot`:
        //
        // - |r| <= u: |r| (|r| / (sqrt(1 + v^2) + 1))
        // - |r| > u: |r| (u / (sqrt(1 + w^2) + w)), with w = u / |r|
        //
        // E.g. pseudo_huber(1e-99, 1e211) was inf instead of 1e112.
        let r = r.abs();
        if r <= delta {
            let v = r / delta;
            r * (r / ((1.0 + v * v).sqrt() + 1.0))
        } else {
            let w = delta / r;
            r * (delta / ((1.0 + w * w).sqrt() + w))
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    fn get_test_xy() -> [(f64, f64); 25] {
        let z1 = [-1.0, -0.5, 0.0, 0.5, 1.0];
        let z2 = z1.map(|x| z1.map(|y| (x, y)));
        z2.as_flattened().try_into().unwrap()
    }

    fn map2<const N: usize>(func: impl Fn(f64, f64) -> f64, xy: [(f64, f64); N]) -> [f64; N] {
        xy.map(|(x, y)| func(x, y))
    }

    /// translated from `scipy.special.tests.test_basic.test_entr`
    #[test]
    fn test_entr() {
        fn xfunc(x: f64) -> f64 {
            if x < 0.0 {
                f64::NEG_INFINITY
            } else {
                -crate::xlogy(x, x)
            }
        }

        let z = [f64::NEG_INFINITY, -1.0, -0.5, 0.0, 0.5, 1.0, f64::INFINITY];
        let w = z.map(xfunc);
        np_assert_allclose!(z.map(crate::entr), w, rtol = 1e-13, atol = 1e-13);
    }

    /// translated from `scipy.special.tests.test_basic.test_kl_div`
    #[test]
    fn test_kl_div() {
        fn xfunc(x: f64, y: f64) -> f64 {
            if x < 0.0 || y < 0.0 || (y == 0.0 && x != 0.0) {
                // extension of natural domain to preserve convexity
                f64::INFINITY
            } else if (x.is_infinite() && x > 0.0) || (y.is_infinite() && y > 0.0) {
                // limits within the natural domain
                f64::INFINITY
            } else if x == 0.0 {
                y
            } else {
                crate::xlogy(x, x / y) - x + y
            }
        }

        let z = get_test_xy();
        let w = map2(xfunc, z);
        np_assert_allclose!(map2(crate::kl_div, z), w, rtol = 1e-13, atol = 1e-13);
    }

    #[test]
    fn test_kl_div_ratio_overflow() {
        // x / y underflows or overflows
        let inputs = [
            (1.126_513_958_653_338_2e-184, 1.521_698_452_194_790_4e258),
            (1e300, 1e-10),
        ];
        // mpmath, x * (log(x) - log(y)) - x + y
        let expected = [1.521_698_452_194_790_4e258, 7.128_013_788_281_542e302];
        np_assert_allclose!(map2(crate::kl_div, inputs), expected, rtol = 1e-15);
    }

    /// translated from `scipy.special.tests.test_basic.test_rel_entr`
    #[test]
    fn test_rel_entr() {
        fn xfunc(x: f64, y: f64) -> f64 {
            if x > 0.0 && y > 0.0 {
                crate::xlogy(x, x / y)
            } else if x == 0.0 && y >= 0.0 {
                0.0
            } else {
                f64::INFINITY
            }
        }

        let z = get_test_xy();
        let w = map2(xfunc, z);
        np_assert_allclose!(map2(crate::rel_entr, z), w, rtol = 1e-13, atol = 1e-13);
    }

    /// translated from `scipy.special.tests.test_basic.test_rel_entr_gh_20710_near_zero`
    #[test]
    fn test_rel_entr_gh_20710_near_zero() {
        // Check accuracy of inputs which are very close
        let inputs = [
            // x, y
            (0.945_665_771_343_000_1, 0.945_665_771_343_009_4),
            (0.480_660_985_647_915_15, 0.480_660_985_647_947_74),
            (0.786_048_657_854_401, 0.786_048_657_854_236_7),
        ];
        // Known values produced using `x * mpmath.log(x / y)` with dps=30
        let expected = [
            -9.325_873_406_851_269e-15,
            -3.258_504_577_274_724e-14,
            1.643_130_076_445_403_3e-13,
        ];
        np_assert_allclose!(map2(crate::rel_entr, inputs), expected, rtol = 1e-13);
    }

    /// translated from `scipy.special.tests.test_basic.test_rel_entr_gh_20710_overflow`
    #[test]
    fn test_rel_entr_gh_20710_overflow() {
        let inputs = [
            // x, y
            // Overflow
            (4.0, 2.22e-308),
            // Underflow
            (1e-200, 1e+200),
            // Subnormal
            (2.22e-308, 1e15),
        ];
        // Known values produced using `x * mpmath.log(x / y)` with dps=30
        let expected = [
            2_839.139_983_229_607,
            -9.210_340_371_976_183e-198,
            -1.649_321_200_807_447_5e-305,
        ];
        np_assert_allclose!(map2(crate::rel_entr, inputs), expected, rtol = 1e-13);
    }

    /// translated from `scipy.special.tests.test_basic.test_huber`
    #[test]
    #[allow(clippy::float_cmp)]
    fn test_huber() {
        assert_eq!(crate::huber(-1.0, 1.5), f64::INFINITY);
        np_assert_allclose!([crate::huber(2.0, 1.5)], [0.5 * 1.5_f64.powi(2)]);
        np_assert_allclose!([crate::huber(2.0, 2.5)], [2.0 * (2.5 - 0.5 * 2.0)]);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_pseudo_huber() {
        assert_eq!(crate::pseudo_huber(-1.0, 1.5), f64::INFINITY);
        np_assert_allclose!([crate::pseudo_huber(2.0, 1.5)], [1.0]);
        np_assert_allclose!([crate::pseudo_huber(2.0, 2.5)], [41.0_f64.sqrt() - 4.0]);
        // delta^2 underflows to 0, but r / delta doesn't
        assert_eq!(crate::pseudo_huber(1e-170, 1e-160), 0.0);
    }

    #[test]
    fn test_pseudo_huber_extreme() {
        // delta^2 or (r / delta)^2 over- or underflows, or |r / delta| is large
        let inputs = [
            (5.632_860_139_466_413e-99, 1.158_620_610_502_460_7e211),
            (1e200, 1.0),
            (1e160, -1e150),
            (1e-170, 1e-100),
            (1.0, 1e100),
            (f64::from_bits(1), 1e300),
            (f64::INFINITY, 1.0),
        ];
        // mpmath, r^2 / (sqrt(1 + (r / delta)^2) + 1)
        let expected = [
            6.526_347_853_663_551e112,
            0.5,
            5e299,
            1e-270,
            1e100,
            4.940_656_458_412_466e-24,
            0.5,
        ];
        np_assert_allclose!(map2(crate::pseudo_huber, inputs), expected, rtol = 1e-15);
    }
}
