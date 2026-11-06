//! Translated into pure Rust from `xsf/cephes/besselpoly.h` (scipy/xsf@09c93c7).
//!
//! This was not part of the original cephes library.
//!
//! Unlike in xsf, the parity of integer `nu` is computed with `nu % 2.0` instead of an `int`
//! cast, which is undefined behavior if `|nu| >= 2^31`.

use crate::xsf::cephes::gamma::gamma;
use crate::xsf::cephes::jv;
use crate::xsf::specfun::msta2;

const BESSELPOLY_EPS: f64 = 1.0e-17;

pub(crate) fn besselpoly(a: f64, lambda: f64, mut nu: f64) -> f64 {
    let mut factor = false;
    let mut sum = 0.0;

    /* Special handling for a = 0.0 */
    if a == 0.0 {
        if nu == 0.0 {
            return 1.0 / (lambda + 1.0);
        }
        return 0.0;
    }
    /* Special handling for negative and integer nu */
    if (nu < 0.0) && (nu.floor() == nu) {
        nu = -nu;
        factor = nu % 2.0 == 1.0;
    }
    // Unlike in xsf (scipy/xsf#452), use a Neumann series where the power series below cancels.
    // For a < 1.5 the power series is more accurate, and for a >= 500 the backward recurrence can
    // overflow.
    if (1.5..500.0).contains(&a)
        && (-1.0..a * a - 1.0).contains(&nu)
        && (-1.0..f64::INFINITY).contains(&lambda)
    {
        let sum = besselpoly_neumann(a, lambda, nu);
        return if factor { -sum } else { sum };
    }
    let mut sm = (nu * a.ln()).exp() / (gamma(nu + 1.0) * (lambda + nu + 1.0));
    let mut m = 0;
    loop {
        sum += sm;
        let sol = sm;
        let fm = f64::from(m);
        sm *= -a * a * (lambda + nu + 1.0 + 2.0 * fm)
            / ((nu + fm + 1.0) * (fm + 1.0) * (lambda + nu + 1.0 + 2.0 * fm + 2.0));
        m += 1;
        let relerr = ((sm - sol) / sm).abs();
        if !(relerr > BESSELPOLY_EPS && m < 1000) {
            break;
        }
    }
    if factor { -sum } else { sum }
}

/// The Neumann series [DLMF 10.22.10] with `x = 2a` and `mu = lambda`,
///
/// ```text
/// 1/(2a) sum_k (nu + 2k + 1) (c1)_k / (c0)_(k+1) J_(nu+2k+1)(2a),
/// ```
///
/// with `c0 = (lambda + nu + 1)/2` and `c1 = (nu - lambda + 1)/2`. Its terms don't grow, since
/// `|c1 + k| <= c0 + k + 1` for `lambda >= -1` and `nu > -1`. The J's come from Miller's backward
/// recurrence, fitted to `jv` at orders `nu + 1` and `nu + 2` by least squares.
///
/// [DLMF 10.22.10]: https://dlmf.nist.gov/10.22.E10
fn besselpoly_neumann(a: f64, lambda: f64, nu: f64) -> f64 {
    let b = 2.0 * a;
    let c0 = (lambda + nu + 1.0) / 2.0;
    let c1 = c0 - lambda;
    let k_max = (f64::from(msta2(b, (nu + 1.0).ceil() as i32, 15)) - nu) as u32 / 2;

    // `(f, g)` are `(J_mu(b), J_(mu+1)(b))` up to a common factor
    let (mut f, mut g) = (1.0, 0.0);
    let mut sum = 0.0;
    for k in (0..=k_max).rev() {
        let k = f64::from(k);
        let mu = nu + 1.0 + 2.0 * k;
        // two steps of J_(mu-1) = 2mu/b J_mu - J_(mu+1), down to order `mu`
        (f, g) = (2.0 * (mu + 2.0) / b * f - g, f);
        (f, g) = (2.0 * (mu + 1.0) / b * f - g, f);
        sum = mu * f + (c1 + k) / (c0 + k + 1.0) * sum;
    }
    sum * (jv(nu + 1.0, b) * f + jv(nu + 2.0, b) * g) / (f * f + g * g) / (b * c0)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_besselpoly_odd_nu_large() {
        // `static_cast<int>(nu)` is UB in xsf for `|nu| >= 2^31`; this `nu` is odd
        let r = super::besselpoly(0.5, 0.0, -2_147_483_649.0);
        assert!(r == 0.0 && r.is_sign_negative(), "{r}");
    }
}
