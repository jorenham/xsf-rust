use core::f64::consts::FRAC_PI_2;

use crate::utils::pow2;

// The needed low limbs of floor(2^2305/pi), in little-endian order.
// Regenerate with Arb: floor(2^2305/pi) mod 2^(23*64).
const RECIPROCAL_BITS: i32 = 2304;
const INV_PIO2: [u64; 23] = [
    0xc4ad414d2c5d000c,
    0x15c614b59d19c3c2,
    0xfa6ed5772d30433b,
    0x87f121907c7c246a,
    0x9f3a1f35caf27f1d,
    0xc33d26ef6b1e5ef8,
    0x32c2de4f98327dbb,
    0xa5ff07053f7e33e8,
    0xddaf44d15719053e,
    0x8359c4768b961ca6,
    0x19c367cddce8092a,
    0x60e27bc08c6b47c4,
    0x06061556ca73a8c9,
    0x8dffd8804d732731,
    0x6599855f14a06840,
    0xa9e391615ee61b08,
    0xf0cfbc209af4361d,
    0x56033046fc7b6bab,
    0x6bfb5fb11f8d5d08,
    0x3d0739f78a5292ea,
    0x7527bac7ebe5f17b,
    0x4f463f669e5fea2d,
    0x6d367ecf27cb09b7,
];

// Payne–Hanek reduction of pi/4 + 2*x^2 when the doubled square overflows.
#[cold]
pub(super) fn sin_diagonal_phase(x: f64) -> f64 {
    debug_assert!(x.is_normal() && (2.0 * x * x).is_infinite());
    let bits = x.to_bits();
    let mantissa = u128::from((bits & ((1 << 52) - 1)) | (1 << 52));
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1075;
    let (quadrant, angle) = reduce_diagonal_phase(mantissa * mantissa, 2 * exponent + 1);
    match quadrant {
        0 => angle.sin(),
        1 => angle.cos(),
        2 => -angle.sin(),
        _ => -angle.cos(),
    }
}

fn reduce_diagonal_phase(mantissa: u128, exponent: i32) -> (u64, f64) {
    let shift = (RECIPROCAL_BITS - exponent) as usize;
    let word = shift / 64;
    let offset = shift % 64;
    // Two extra reciprocal words absorb all but 2^-214 quadrants of the omitted product.
    let mut carry = 0_u128;
    let [_, _, lo, mid, hi, top, next] = core::array::from_fn(|i| {
        let digit = u128::from(INV_PIO2[word - 5 + i]);
        let low = u128::from(mantissa as u64) * digit + u128::from(carry as u64);
        carry = (mantissa >> 64) * digit + (carry >> 64) + (low >> 64);
        low as u64
    });
    let join = |lo, hi| (((u128::from(hi) << 64) | u128::from(lo)) >> offset) as u64;
    let quadrant = ((join(top, next) & 3) + 1) & 3;
    let hi_bits = join(hi, top);
    // Adding pi/4 gives remainder f - 1/2 and advances the quadrant by one.
    // The exact distance from an axis is > 2^-117 quadrants (continued-fraction bound).
    let sign = if hi_bits >> 63 == 0 { u64::MAX } else { 0 };
    // One's complement differs from the negative magnitude by only 2^-192 quadrants.
    let lo = (join(lo, mid) ^ sign) as f64;
    let mid = (join(mid, hi) ^ sign) as f64;
    let hi = (hi_bits ^ (1 << 63) ^ sign) as f64;
    let frac = ((lo * pow2(-64) + mid) * pow2(-64) + hi) * pow2(-64);
    debug_assert!(frac > pow2(-118));
    let angle = frac * FRAC_PI_2;
    (quadrant, if sign == 0 { angle } else { -angle })
}
