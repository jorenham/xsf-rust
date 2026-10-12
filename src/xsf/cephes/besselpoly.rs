//! Translated into pure Rust from `xsf/cephes/besselpoly.h` (xsf v0.2.2), which was translated
//! into C++ by SciPy developers in 2024.
//!
//! This was not part of the original cephes library.
//!
//! Unlike in xsf, the parity of integer `nu` is computed with `nu % 2.0` instead of an `int`
//! cast, which is undefined behavior if `|nu| >= 2^31`.

use crate::xsf::cephes::gamma::gamma;

const BESSELPOLY_EPS: f64 = 1.0e-17;

#[allow(clippy::float_cmp)]
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

#[cfg(test)]
mod tests {
    #[test]
    fn test_besselpoly_odd_nu_large() {
        // `static_cast<int>(nu)` is UB in xsf for `|nu| >= 2^31`; this `nu` is odd
        let r = super::besselpoly(0.5, 0.0, -2_147_483_649.0);
        assert!(r == 0.0 && r.is_sign_negative(), "{r}");
    }
}
