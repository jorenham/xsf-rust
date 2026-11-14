//! Translated from scipy/special/_ellip_harm_2.pyx
//!
//! Unlike in SciPy, the integrals are computed after a change of variables that removes the
//! (weighted) endpoint singularities, so that [`quad`] is enough, instead of QUADPACK's `dqags`
//! and `dqaws`.

use core::f64::consts::{FRAC_PI_2, PI};

use super::ellip_harm::{Lame, sq_sub};
use super::quadpack::quad;

/// The relative tolerance of the integrals
///
/// Unlike in SciPy, not `1e-15`, which [`quad`] can't reach. Its error estimate is so pessimistic
/// that `1e-13` already gives full accuracy.
const EPSREL: f64 = 1e-13;

/// `(t, t^2 - h^2, t^2 - k^2)`, with the differences computed without cancellation
type Point = (f64, f64, f64);

/// Ellipsoidal harmonic function $F^p_n(s)$
///
/// Also known as the Lamé function of the second kind, a solution of the Lamé equation (see
/// [`ellip_harm`](crate::ellip_harm)), related to the one of the first kind by
///
/// $$
/// F^p_n(s) = (2n + 1) E^p_n(s)
/// \int_0^{1/s} \frac{du}{(E^p_n(1/u))^2 \sqrt{(1 - u^2 k^2)(1 - u^2 h^2)}}
/// $$
///
/// Corresponds to [`scipy.special.ellip_harm_2`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellip_harm_2.html
///
/// # Arguments
/// - `h2`: $h^2$
/// - `k2`: $k^2$, should be larger than `h2`
/// - `n`: degree
/// - `p`: order, `1 <= p <= 2n + 1`
/// - `s`: coordinate, `s >= k`
///
/// Returns NaN if `p` is invalid, or if the integral doesn't converge, e.g. at `s = k` for the
/// functions of type M and N (see [`ellip_harm`](crate::ellip_harm)).
///
/// # See also
/// - [`ellip_harm`](crate::ellip_harm): Ellipsoidal harmonic function $E^p_n(s)$
/// - [`ellip_normal`]: Ellipsoidal harmonic normalization constant $\gamma^p_n$
#[must_use]
pub fn ellip_harm_2(h2: f64, k2: f64, n: u32, p: u32, s: f64) -> f64 {
    let Some(lame) = Lame::new(h2, k2, n, p) else {
        return f64::NAN;
    };

    // `sqrt(k2) = k + dk`, and `s - sqrt(k2)`, where `s = k` is taken as `s >= sqrt(k2)`
    let k = k2.sqrt();
    let dk = -sq_sub(k, k2) / (2.0 * k);
    let s_k = (s - k) - dk;
    let s_k = if s == k { s_k.max(0.0) } else { s_k };
    // `s^2 - k2` and `s^2 - h2`, consistent with `1 / t^2 - k2` and `1 / t^2 - h2` below at
    // `t = 1/s`, which matters near `s = k`
    let s2_k2 = s_k * (s + k);
    let e_s = lame.eval_with(s, s2_k2 + (k2 - h2), s2_k2, 1.0, 1.0);
    if e_s.is_infinite() {
        // `F(s) ~ 1 / (s E(s))` underflows
        return 0.0;
    }

    // `_F_integrand` over `[0, 1/s]`, with `t = (1 - v^2) / sqrt(k2)` (up to the rounding of `k`),
    // so that its `1 / sqrt(1 - t sqrt(k2))` is `1 / v`, which cancels, also for `s = k`
    //
    // - `1 / t^2 - k2` and `1 / t^2 - h2` follow from `v`, without cancellation;
    // - `E^2(s)` is folded in, so that the integral doesn't underflow for large `s`.
    let integrand = |v: f64, t: f64| {
        let one_m_t2k2 = v * v * (1.0 + t * k);
        let one_m_t2h2 = one_m_t2k2 + t * t * (k2 - h2);
        let i = lame.eval_with(
            1.0 / t,
            one_m_t2h2 / (t * t),
            one_m_t2k2 / (t * t),
            1.0,
            1.0,
        );
        let r = e_s / i;
        2.0 / k * r * r / ((1.0 + t * k).sqrt() * one_m_t2h2.sqrt())
    };
    // over `v in [v1, 1]`, split at `v = 1/2`, so that `v ~ 0` (near `s = k`) and `w = 1 - v ~ 0`
    // (`t ~ 0` for large `s`) are both resolved
    let v1 = (s_k / s).sqrt();
    let (res_v, err_v) = if v1 < 0.5 {
        quad(|v| integrand(v, (1.0 - v * v) / k), v1, 0.5, 0.0, EPSREL)
    } else {
        (0.0, 0.0)
    };
    let w1 = k / s / (1.0 + v1);
    let (res_w, err_w) = quad(
        |w| integrand(1.0 - w, w * (2.0 - w) / k),
        0.0,
        w1.min(0.5),
        0.0,
        EPSREL,
    );
    let (res, err) = (res_v + res_w, err_v + err_w);
    if err > 1e-10 * res.abs() + 1e-290 {
        return f64::NAN;
    }
    res * (2.0 * f64::from(n) + 1.0) / e_s
}

/// Ellipsoidal harmonic normalization constant $\gamma^p_n$
///
/// $$
/// \gamma^p_n = 8 \int_0^h dx \int_h^k dy
/// \frac{(y^2 - x^2) (E^p_n(y) E^p_n(x))^2}{\sqrt{(k^2 - y^2)(y^2 - h^2)(h^2 - x^2)(k^2 - x^2)}}
/// $$
///
/// Corresponds to [`scipy.special.ellip_normal`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellip_normal.html
///
/// # Arguments
/// - `h2`: $h^2$
/// - `k2`: $k^2$, should be larger than `h2`
/// - `n`: degree
/// - `p`: order, `1 <= p <= 2n + 1`
///
/// Returns NaN if `p` is invalid, or if the integrals don't converge.
///
/// # See also
/// - [`ellip_harm`](crate::ellip_harm): Ellipsoidal harmonic function $E^p_n(s)$
/// - [`ellip_harm_2`]: Ellipsoidal harmonic function $F^p_n(s)$
#[must_use]
pub fn ellip_normal(h2: f64, k2: f64, n: u32, p: u32) -> f64 {
    let Some(lame) = Lame::new(h2, k2, n, p) else {
        return f64::NAN;
    };
    let (h, k) = (h2.sqrt(), k2.sqrt());
    let e2 = |(t, t2_h2, t2_k2): Point| {
        let i = lame.eval_with(t, t2_h2, t2_k2, 1.0, 1.0);
        i * i
    };

    // the integrals of `f(t(x))` and `t(x)^2 f(t(x))` over `[0, b]`
    let moments = |f: &dyn Fn(Point) -> f64, t: &dyn Fn(f64) -> Point, b: f64| {
        let m0 = quad(|x| f(t(x)), 0.0, b, 0.0, EPSREL);
        let m2 = quad(
            |x| {
                let t = t(x);
                t.0 * t.0 * f(t)
            },
            0.0,
            b,
            0.0,
            EPSREL,
        );
        (m0, m2)
    };

    // `_F_integrand1` and `_F_integrand2` over `[h, k]` with the weight `1 / sqrt((t - h)(k - t))`,
    // so with `t = h + (k - h) sin^2(theta / 2)` over `[0, pi]`
    let k_h = (k2 - h2) / (k + h);
    let ((res, err), (res1, err1)) = moments(
        &|p| e2(p) / ((p.0 + h) * (p.0 + k)).sqrt(),
        &|theta| {
            let (sin, cos) = (0.5 * theta).sin_cos();
            let (t_h, k_t) = (k_h * sin * sin, k_h * cos * cos);
            let t = h + t_h;
            (t, t_h * (t + h), -k_t * (k + t))
        },
        PI,
    );

    // `_F_integrand3` and `_F_integrand4` over `[0, h]` with the weight `1 / sqrt(h - t)`, so with
    // `t = h sin(phi)` over `[0, pi/2]`, which leaves `E^2(t) / sqrt(k^2 - t^2)`
    let ((res2, err2), (res3, err3)) = moments(
        &|p| e2(p) / (-p.2).sqrt(),
        &|phi| {
            let (sin, cos) = phi.sin_cos();
            let h2_t2 = h2 * cos * cos;
            (h * sin, -h2_t2, -((k2 - h2) + h2_t2))
        },
        FRAC_PI_2,
    );

    let error = 8.0 * (res2 * err1 + err2 * res1 + res * err3 + res3 * err);
    let result = 8.0 * (res1 * res2 - res * res3);
    if error > 10e-8 * result.abs() {
        return f64::NAN;
    }
    result
}
