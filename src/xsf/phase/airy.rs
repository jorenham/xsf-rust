use super::{INV_PIO2, Quadrant, RECIPROCAL_BITS, center_pi_4, decompose, read_quadrant_window};
use crate::utils::pow2;

// The fixed-point root has at most 1728 fractional bits and 27 integer bits.
const AIRY_ROOT_WORDS: usize = 29;

fn sqrt_scaled(radicand: u64, fractional_bits: usize) -> [u64; AIRY_ROOT_WORDS] {
    let prefix = u128::from(radicand) << 64;
    let first = prefix.isqrt();
    let mut root = [0_u64; AIRY_ROOT_WORDS];
    let mut remainder = [0_u64; AIRY_ROOT_WORDS];
    root[0] = first as u64;
    remainder[0] = (prefix - first * first) as u64;
    let mut negative = false;
    // Non-restoring square root. A negative r stores the true remainder minus 2q + 1.
    // Subtract the next trial, adjusting for that offset, then append the decision bit to q.
    for bit in 32..fractional_bits {
        // Allow 27 integer bits for q and four more for the signed trial.
        let active = (bit + 31).div_ceil(64);
        let (trial_factor, mut carry) = if negative { (4, 3) } else { (-4, -1) };
        let mut root_carry = 0;
        for (r, s) in remainder[..active].iter_mut().zip(&mut root) {
            let digit = 4 * i128::from(*r) + trial_factor * i128::from(*s) + carry;
            *r = digit as u64;
            carry = digit >> 64;
            let next = *s >> 63;
            *s = (*s << 1) | root_carry;
            root_carry = next;
        }
        negative = remainder[active - 1] >> 63 != 0;
        // Sign-extend the remainder in case the next step needs another limb.
        remainder[active] = if negative { u64::MAX } else { 0 };
        root[0] |= u64::from(!negative);
    }
    root
}

// Reduce (2/3) x^(3/2) + pi/4 in the Airy asymptotic expansion.
#[cold]
pub(in crate::xsf) fn reduce_airy_phase(x: f64) -> Option<(Quadrant, f64)> {
    // (x >= 1 keeps d + 53 + 195 > 0 below)
    if !(1.0..=f64::MAX).contains(&x) {
        return None;
    }
    let (m, e) = decompose(x);
    // x^(3/2) = m sqrt(m 2^(e mod 2)) 2^d
    let d = e + e.div_euclid(2);
    // The 53-bit mantissa times the square root error costs < 2^-195 quadrants.
    let frac_bits = ((d + 53 + 195) as usize).next_multiple_of(64);
    let mut root = sqrt_scaled(m << e.rem_euclid(2), frac_bits);
    // Reserve 27 integer bits for the root and 54 for twice the significand.
    let digits = &mut root[..(frac_bits + 81).div_ceil(64)];

    let mut carry = 0_u128;
    for v in &mut *digits {
        let product = u128::from(*v) * u128::from(2 * m) + carry;
        *v = product as u64;
        carry = product >> 64;
    }
    let mut rem = 0;
    for v in digits.iter_mut().rev() {
        let wide = (rem << 64) | u128::from(*v);
        *v = (wide / 3) as u64;
        rem = wide % 3;
    }

    let reciprocal = &INV_PIO2[..INV_PIO2.len() - 2];
    let mut product = [0_u64; AIRY_ROOT_WORDS + INV_PIO2.len() - 2];
    for (i, &digit) in digits.iter().enumerate() {
        let mut carry = 0;
        for (&constant, out) in reciprocal.iter().zip(&mut product[i..]) {
            let row = u128::from(digit) * u128::from(constant) + u128::from(*out) + carry;
            *out = row as u64;
            carry = row >> 64;
        }
        product[i + reciprocal.len()] = carry as u64;
    }

    let shift = (RECIPROCAL_BITS - d) as usize + frac_bits;
    let (word, bit) = (shift / 64, shift % 64);
    let words = product.get(word - 3..)?.first_chunk()?;
    let window = read_quadrant_window(words, bit);
    let [_, lo, mid, hi] = window;
    // Unlike in xsf, decline if the quadrant or the remainder sign cannot be resolved.
    // Truncation is downward and costs < two quanta of the 192-bit window.
    if hi == u64::MAX && mid == u64::MAX && lo >= u64::MAX - 2 {
        return None;
    }
    let (quadrant, angle) = center_pi_4(window);
    // Leave enough bits in the remainder to resolve its sign and magnitude.
    (angle.abs() >= pow2(-119)).then_some((quadrant, angle))
}
