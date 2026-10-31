//! Translated into pure Rust from `xsf/ellint_carlson/ellint_arith_aux.hh` (xsf main).
//!
//! ```text
//! Auxiliary floating-point manipulation utilities.
//! Ref:
//!
//! S. M. Rump, T. Ogita, S. Oishi: Accurate Floating-Point Summation.
//!      2005-11-13, Tech. Rep. 05.1, Fac. Inf. Commun. Sci, Hambg. Univ.
//!      Technol., also known as SIAM J. Sci. Comput. Volume 31, Issue 1, pp.
//!      189-224 (2008),
//!      http://www.ti3.tuhh.de/paper/rump/RuOgOi07I.pdf
//!
//! The main purpose of these templates is to implement Algorithm 6.1 (here as
//! acc_sum) which is required by compensated complex floating-point polynomial
//! evaluation.
//! ```

/// `std::numeric_limits<double>::epsilon() * 0.5`
const HALF_EPSILON: f64 = f64::EPSILON * 0.5;

/* Algorithm 3.6 */
#[inline]
fn next_power_two(p: f64) -> f64 {
    let q = p / HALF_EPSILON;
    let mut l = ((q + p) - q).abs();
    if l == 0.0 {
        l = p.abs();
    }
    l
}

/* Algorithm 3.2 */
#[inline]
fn extract_scalar(sigma: f64, p: &mut f64, q: &mut f64) {
    *q = (sigma + *p) - sigma;
    *p -= *q;
}

/* Algorithm 3.4 */
#[inline]
fn extract_vector<const LEN: usize>(
    sigma: f64,
    p: &mut [f64; LEN],
    mask: &mut [bool; LEN],
    tau: &mut f64,
) {
    let mut q = 0.0;
    *tau = 0.0;
    for i in 0..LEN {
        if mask[i] {
            extract_scalar(sigma, &mut p[i], &mut q);
            if p[i] == 0.0 {
                mask[i] = false;
            }
            *tau += q;
        }
    }
}

/* The functions make, none, and count implements boolean-array
 * initialization, test for all-false, and counting for number of true
 * values. */
#[inline]
pub(super) fn make<const LEN: usize>(buf: &[f64; LEN], mask: &mut [bool; LEN]) {
    for i in 0..LEN {
        mask[i] = buf[i] != 0.0;
    }
}

#[inline]
fn none<const LEN: usize>(mask: &[bool; LEN]) -> bool {
    !mask.contains(&true)
}

#[inline]
fn count<const LEN: usize>(mask: &[bool; LEN]) -> usize {
    mask.iter().filter(|&&m| m).count()
}

/* Search for maximal absolute value in masked array. */
#[inline]
fn vmax<const LEN: usize>(p: &[f64; LEN], mask: &[bool; LEN]) -> f64 {
    let mut v = 0.0;
    for i in 0..LEN {
        if mask[i] {
            // std::max(v, std::abs(p[i]))
            let a = p[i].abs();
            if v < a {
                v = a;
            }
        }
    }
    v
}

/* Summation for masked array */
#[inline]
fn masked_sum<const LEN: usize>(p: &[f64; LEN], mask: &[bool; LEN]) -> f64 {
    let mut v = 0.0;
    for i in 0..LEN {
        if mask[i] {
            v += p[i];
        }
    }
    v
}

/* Algorithm 6.1 */
#[allow(clippy::cast_precision_loss)]
pub(super) fn acc_sum<const LEN: usize>(p: &mut [f64; LEN], mask: &mut [bool; LEN]) -> f64 {
    if LEN == 0 || none(mask) {
        return 0.0;
    }

    let mu = vmax(p, mask);
    if mu == 0.0 {
        return mu;
    }

    let two_pm = next_power_two((count(mask) + 2) as f64);
    let mut sigma = two_pm * next_power_two(mu);
    let phi = HALF_EPSILON * two_pm;
    let factor = f64::EPSILON * two_pm * two_pm;

    let mut t = 0.0;
    loop {
        let mut tau = 0.0;
        extract_vector(sigma, p, mask, &mut tau);
        let tau1 = t + tau;
        if tau1.abs() >= factor * sigma || sigma <= f64::MIN_POSITIVE {
            /* Dekker "FastTwoSum" algorithm */
            let tau2 = tau - (tau1 - t);
            return tau1 + (tau2 + masked_sum(p, mask));
        }
        t = tau1;
        if t == 0.0 {
            return acc_sum(p, mask);
        }
        sigma *= phi;
    }
}
