use core::f64::consts::FRAC_PI_2;

use crate::utils::pow2;

mod airy;

pub(super) use airy::reduce_airy_phase;

#[derive(Clone, Copy)]
pub(super) enum Quadrant {
    Zero,
    QuarterTurn,
    HalfTurn,
    ThreeQuarterTurns,
}

impl Quadrant {
    pub(super) fn from_index(index: u64) -> Self {
        match index & 3 {
            0 => Self::Zero,
            1 => Self::QuarterTurn,
            2 => Self::HalfTurn,
            _ => Self::ThreeQuarterTurns,
        }
    }

    // Rotate a reduced sine and cosine, returning (cos, sin).
    pub(super) fn rotate(self, sin: f64, cos: f64) -> (f64, f64) {
        match self {
            Self::Zero => (cos, sin),
            Self::QuarterTurn => (-sin, cos),
            Self::HalfTurn => (-cos, -sin),
            Self::ThreeQuarterTurns => (sin, -cos),
        }
    }
}

// Little-endian limbs of floor(2^2305/pi), with two high zero words for the window.
const RECIPROCAL_BITS: i32 = 2304;
const INV_PIO2: [u64; 38] = [
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
    0xef2f118b5a0a6d1f,
    0x1ff897ffde05980f,
    0x9c845f8bbdf9283b,
    0x3991d639835339f4,
    0xe99c7026b45f7e41,
    0xe88235f52ebb4484,
    0xfe1deb1cb129a73e,
    0x06492eea09d1921c,
    0xb7246e3a424dd2e0,
    0xfe5163abdebbc561,
    0xdb6295993c439041,
    0xfc2757d1f534ddc0,
    0xa2f9836e4e441529,
    0,
    0,
];

fn decompose(x: f64) -> (u64, i32) {
    let bits = x.to_bits();
    let mantissa = (bits & ((1 << 52) - 1)) | (1 << 52);
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1075;
    (mantissa, exponent)
}

// In each binade from 2^55 to 2^82, the phase separation exceeds 8/(4*x^2).
pub(super) const MIN_DIAGONAL_TAIL_ARG: f64 = pow2(55);

// Payne–Hanek reduction of pi/4 + 2*x^2 without rounding the doubled square.
#[cold]
pub(super) fn sin_diagonal_tail_phase(x: f64) -> f64 {
    let (quadrant, angle) = reduce_diagonal_phase(x, 2.0);
    // Unlike in xsf, absorb the first erfc tail correction into the reduced angle.
    // The leading term alone has a relative remainder below 2^-57 from 2^82.
    let angle = if x.abs() < pow2(82) {
        angle - 0.25 / (x * x)
    } else {
        angle
    };
    match quadrant {
        Quadrant::Zero => angle.sin(),
        Quadrant::QuarterTurn => angle.cos(),
        Quadrant::HalfTurn => -angle.sin(),
        Quadrant::ThreeQuarterTurns => -angle.cos(),
    }
}

// The two low integer bits and the high 192 fractional bits of the product in quadrants.
fn quadrant_window(mantissa: u128, exponent: i32) -> [u64; 4] {
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
    read_quadrant_window(&[lo, mid, hi, top, next], offset)
}

// The first three words are below the quadrant, with a high word for alignment.
fn read_quadrant_window(words: &[u64; 5], offset: usize) -> [u64; 4] {
    let join =
        |i: usize| (((u128::from(words[i + 1]) << 64) | u128::from(words[i])) >> offset) as u64;
    [join(3) & 3, join(0), join(1), join(2)]
}

// Reduction of pi/4 + k*x^2 for k = 1 or 2.
pub(super) fn reduce_diagonal_phase(x: f64, k: f64) -> (Quadrant, f64) {
    debug_assert!(x.is_finite() && x.abs() >= MIN_DIAGONAL_TAIL_ARG && (k == 1.0 || k == 2.0));
    let (mantissa, exponent) = decompose(x);
    // The exact distance from an axis is > 2^-117 quadrants (continued-fraction bound).
    center_pi_4(quadrant_window(u128::from(mantissa).pow(2), 2 * exponent + i32::from(k == 2.0)))
}

fn center_pi_4([quadrant, lo, mid, hi]: [u64; 4]) -> (Quadrant, f64) {
    // Adding pi/4 gives remainder f - 1/2 and advances the quadrant by one.
    let sign = if hi >> 63 == 0 { u64::MAX } else { 0 };
    // One's complement differs from the negative magnitude by only 2^-192 quadrants.
    let lo = (lo ^ sign) as f64;
    let mid = (mid ^ sign) as f64;
    let hi = (hi ^ (1 << 63) ^ sign) as f64;
    let frac = ((lo * pow2(-64) + mid) * pow2(-64) + hi) * pow2(-64);
    let angle = frac * FRAC_PI_2;
    (Quadrant::from_index(quadrant + 1), if sign == 0 { angle } else { -angle })
}

// The erfc sign margin needs only 64 fractional bits. The wide product implies
// normal inputs and a shift in 361..=2356, within the reciprocal window above.
#[cold]
pub(super) fn reduce_product_phase(x: f64, y: f64) -> (Quadrant, f64) {
    debug_assert!(x.is_normal() && y.is_normal() && (2.0 * (x * y)).abs() > pow2(53));
    let (mx, ex) = decompose(x);
    let (my, ey) = decompose(y);
    let [q, _, _, frac] = quadrant_window(u128::from(mx) * u128::from(my), ex + ey + 1);
    let q = q + (frac >> 63);
    let r = (i64::from_ne_bytes(frac.to_ne_bytes()) as f64) * (FRAC_PI_2 * pow2(-64));
    if x.is_sign_negative() ^ y.is_sign_negative() {
        (Quadrant::from_index(q.wrapping_neg()), -r)
    } else {
        (Quadrant::from_index(q), r)
    }
}
