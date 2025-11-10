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
#[inline]
pub fn kl_div(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        f64::NAN
    } else if 0.0 < x && 0.0 < y {
        x * (x / y).ln() - x + y
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
#[inline]
pub fn pseudo_huber(delta: f64, r: f64) -> f64 {
    if delta < 0.0 {
        f64::INFINITY
    } else if delta == 0.0 || r == 0.0 {
        0.0
    } else {
        let u2 = delta * delta;
        let v2 = r * r / u2;
        // The formula is u^2 * (sqrt(v^2 + 1) - 1), but to maintain precision with small v, we use
        // sqrt(1 + v^2) - 1
        //   =  exp(log(v^2 + 1) / 2) - 1
        //   =  expm1(log1p(v^2) / 2)
        u2 * (0.5 * v2.ln_1p()).exp_m1()
    }
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    fn get_test_xy() -> [(f64, f64); 25] {
        let z1 = [-1.0, -0.5, 0.0, 0.5, 1.0];
        let z2 = z1.map(|x| z1.map(|y| (x, y)));
        unsafe { core::mem::transmute::<[[(f64, f64); 5]; 5], [(f64, f64); 25]>(z2) }
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
        np_assert_allclose!(z.map(crate::entr), w, rtol = 1e-13, atol = 1e-13)
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
        np_assert_allclose!(map2(crate::kl_div, z), w, rtol = 1e-13, atol = 1e-13)
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
        np_assert_allclose!(map2(crate::rel_entr, z), w, rtol = 1e-13, atol = 1e-13)
    }

    /// translated from `scipy.special.tests.test_basic.test_rel_entr_gh_20710_near_zero`
    #[test]
    fn test_rel_entr_gh_20710_near_zero() {
        // Check accuracy of inputs which are very close
        let inputs = [
            // x, y
            (0.9456657713430001, 0.9456657713430094),
            (0.48066098564791515, 0.48066098564794774),
            (0.786048657854401, 0.7860486578542367),
        ];
        // Known values produced using `x * mpmath.log(x / y)` with dps=30
        let expected = [
            -9.325873406851269e-15,
            -3.258504577274724e-14,
            1.6431300764454033e-13,
        ];
        np_assert_allclose!(map2(crate::rel_entr, inputs), expected, rtol = 1e-13)
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
            2839.139983229607,
            -9.210340371976183e-198,
            -1.6493212008074475e-305,
        ];
        np_assert_allclose!(map2(crate::rel_entr, inputs), expected, rtol = 1e-13)
    }

    /// translated from `scipy.special.tests.test_basic.test_huber`
    #[test]
    fn test_huber() {
        assert_eq!(crate::huber(-1.0, 1.5), f64::INFINITY);
        np_assert_allclose!([crate::huber(2.0, 1.5)], [0.5 * 1.5_f64.powi(2)]);
        np_assert_allclose!([crate::huber(2.0, 2.5)], [2.0 * (2.5 - 0.5 * 2.0)]);
    }

    #[test]
    fn test_pseudo_huber() {
        assert_eq!(crate::pseudo_huber(-1.0, 1.5), f64::INFINITY);
        np_assert_allclose!([crate::pseudo_huber(2.0, 1.5)], [1.0]);
        np_assert_allclose!([crate::pseudo_huber(2.0, 2.5)], [41.0_f64.sqrt() - 4.0]);
    }
}
