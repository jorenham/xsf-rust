//! Translated into pure Rust from `xsf/cephes/besselpoly.h` (scipy/xsf@09c93c7).
//!
//! This was not part of the original cephes library.
//!
//! Unlike in xsf, the parity of integer `nu` is computed with `nu % 2.0` instead of an `int`
//! cast, which is undefined behavior if `|nu| >= 2^31`.

use crate::xsf::cephes::gamma::{gamma, lgam_sgn};
use crate::xsf::cephes::pow_rgamma;

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
    // Unlike in xsf, `a^nu / Gamma(nu + 1)` comes from `pow_rgamma`, since `exp(nu ln(a))` loses
    // `|nu ln(a)| eps`
    let mut sm = match pow_rgamma(a, nu, 0, nu + 1.0, 0.0) {
        Some(t) => t / (lambda + nu + 1.0),
        None => (nu * a.ln()).exp() / (gamma(nu + 1.0) * (lambda + nu + 1.0)),
    };
    if !sm.is_normal() {
        // Unlike in xsf (scipy/xsf#538), the leading term is recomputed in log space if
        // `gamma(nu + 1)` or `a^nu` over- or underflows (e.g. `nu > 171`)
        let (lg, sign) = lgam_sgn(nu + 1.0);
        let d = lambda + nu + 1.0;
        sm = f64::from(sign) * (nu * a.ln() - lg - d.abs().ln()).exp().copysign(d);
    }
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
