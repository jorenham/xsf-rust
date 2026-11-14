use core::f64::consts::PI;

/// ln(pi)
pub(crate) const LN_PI: f64 = 1.144_729_885_849_400_2;
/// ln(2 pi), as `TAU.ln()`
pub(crate) const LN_TAU: f64 = 1.837_877_066_409_345_3;
/// ln(sqrt(2 pi)), correctly rounded (1 ulp above `0.5 * LN_TAU`)
pub(crate) const LN_SQRT_TAU: f64 = 0.918_938_533_204_672_8;
/// sqrt(pi), as `PI.sqrt()` (1 ulp below the correctly rounded value)
pub(crate) const SQRT_PI: f64 = 1.772_453_850_905_515_9;
/// sqrt(2 pi), as `TAU.sqrt()` (1 ulp below the correctly rounded `SQRT_TAU`)
pub(crate) const SQRT_TAU_DOWN: f64 = 2.506_628_274_631_000_2;
/// sqrt(2 pi), correctly rounded
pub(crate) const SQRT_TAU: f64 = 2.506_628_274_631_000_7;
/// sqrt(f64::EPSILON)
pub(crate) const SQRT_EPSILON: f64 = 1.490_116_119_384_765_6e-8;
/// f64::EPSILON^2
pub(crate) const EPSILON_SQ: f64 = f64::EPSILON * f64::EPSILON;
/// pi^2 / 6
pub(crate) const FRAC_PI_SQ_6: f64 = PI * PI / 6.0;
/// Euler-Mascheroni constant
pub(crate) const EGAMMA: f64 = 0.577_215_664_901_532_9;
/// sqrt(1/pi)
pub(crate) const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;
/// sqrt(2/pi)
pub(crate) const SQRT_FRAC_2_PI: f64 = 0.797_884_560_802_865_4;
/// sqrt(3)
pub(crate) const SQRT_3: f64 = 1.732_050_807_568_877_2;

/// `2^k` for a normal result (`-1022 <= k <= 1023`)
pub(crate) const fn pow2(k: i32) -> f64 {
    assert!(-1022 <= k && k <= 1023);
    f64::from_bits(((1023 + k) as u64) << 52)
}

/// `m * e^e`, without spurious over- or underflow of `e^e`
///
/// This is `m * h * h` with `h = e^(e/2)`, which is a normal float for `|e| <= 1400 = 2 * 700`
/// (`ln(f64::MIN_POSITIVE) ~ -708.4`), and `e^(ln|m| + e)` beyond that.
pub(crate) fn mul_exp(m: f64, e: f64) -> f64 {
    if e.abs() > 1400.0 {
        return m.signum() * (m.abs().ln() + e).exp();
    }
    let h = (0.5 * e).exp();
    m * h * h
}

/// Whether `x` is zero or a negative integer, i.e. a pole of the Gamma function, as
/// `x <= 0 && x == std::floor(x)` in xsf (also true for `-inf`)
pub(crate) fn is_nonpos_int(x: f64) -> bool {
    x <= 0.0 && x == x.floor()
}

/// The rows of a row-major matrix with `cols` columns, with the elements converted to `U`
#[inline]
pub(crate) fn vec_to_vecvec<T: Copy + Into<U>, U>(vec: &[T], cols: usize) -> Vec<Vec<U>> {
    vec.chunks_exact(cols)
        .map(|row| row.iter().map(|&x| x.into()).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::{FRAC_2_PI, FRAC_2_SQRT_PI, LN_2, LN_10, SQRT_2, TAU};
    use std::hint::black_box;

    /// The constants must be bit-identical to the expressions they replace
    #[test]
    fn test_consts_bits() {
        let pi = black_box(PI);
        let tau = black_box(TAU);
        assert_eq!(LN_PI.to_bits(), pi.ln().to_bits());
        assert_eq!(LN_TAU.to_bits(), tau.ln().to_bits());
        assert_eq!(LN_TAU.to_bits(), (2.0 * pi).ln().to_bits());
        assert_eq!(LN_SQRT_TAU.to_bits(), (0.5 * LN_TAU).to_bits() + 1);
        assert_eq!(SQRT_PI.to_bits(), pi.sqrt().to_bits());
        assert_eq!(SQRT_TAU_DOWN.to_bits(), tau.sqrt().to_bits());
        assert_eq!(SQRT_TAU.to_bits(), SQRT_TAU_DOWN.to_bits() + 1);
        assert_eq!(SQRT_3.to_bits(), black_box(3.0_f64).sqrt().to_bits());
        assert_eq!(
            FRAC_PI_SQ_6.to_bits(),
            1.644_934_066_848_226_4_f64.to_bits()
        );
        assert_eq!(
            SQRT_EPSILON.to_bits(),
            black_box(f64::EPSILON).sqrt().to_bits()
        );
        assert_eq!(FRAC_1_SQRT_PI.to_bits(), (1.0 / pi.sqrt()).to_bits());
        assert_eq!(SQRT_FRAC_2_PI.to_bits(), (2.0 / pi).sqrt().to_bits());
        assert_eq!(FRAC_2_PI.to_bits(), (2.0 / pi).to_bits());
        assert_eq!(FRAC_2_SQRT_PI.to_bits(), (2.0 / pi.sqrt()).to_bits());
        assert_eq!(SQRT_2.to_bits(), black_box(2.0_f64).sqrt().to_bits());
        let (ln10, ln2) = (black_box(10.0_f64).ln(), black_box(2.0_f64).ln());
        assert_eq!((LN_10 / LN_2).to_bits(), (ln10 / ln2).to_bits());
        assert!((-1022..=1023).all(|k| pow2(k) == 2.0_f64.powi(k)));
    }
}
