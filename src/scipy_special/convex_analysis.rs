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
        // Unlike in SciPy (scipy/xsf#362), `rel_entr` is reused, and y ((1 + e) ln(1 + e) - e) is
        // used near x = y, e.g. kl_div(1e-184, 1e258) was -inf instead of 1e258.
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
        // Unlike in SciPy (scipy/xsf#362), this uses r^2 / (sqrt(1 + v^2) + 1) with v = r / delta,
        // scaled like `hypot`, instead of delta^2 expm1(log1p(v^2) / 2), which over- or underflows,
        // and loses precision for large |v|, e.g. pseudo_huber(1e-99, 1e211) was inf instead of
        // 1e112.
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
