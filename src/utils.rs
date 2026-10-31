use core::f64::consts::PI;

/// ln(pi)
pub(crate) const LN_PI: f64 = 1.144_729_885_849_400_2;
/// ln(2 pi), as `TAU.ln()`
pub(crate) const LN_TAU: f64 = 1.837_877_066_409_345_3;
/// ln(sqrt(2 pi)), correctly rounded (1 ulp above `0.5 * LN_TAU`)
pub(crate) const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;
/// sqrt(pi), as `PI.sqrt()` (1 ulp below the correctly rounded value)
pub(crate) const SQRT_PI: f64 = 1.772_453_850_905_515_9;
/// sqrt(2 pi), as `TAU.sqrt()` (1 ulp below the correctly rounded `SQRT2PI`)
pub(crate) const SQRT_TAU: f64 = 2.506_628_274_631_000_2;
/// sqrt(2 pi), correctly rounded
pub(crate) const SQRT2PI: f64 = 2.506_628_274_631_000_7;
/// sqrt(f64::EPSILON)
pub(crate) const SQRT_EPS: f64 = 1.490_116_119_384_765_6e-8;
/// pi^2 / 6
pub(crate) const PISQ_6: f64 = PI * PI / 6.0;

// cephes/const.h

/// 2**-53
pub(crate) const MACHEP: f64 = 1.110_223_024_625_156_5e-16;
/// log(DBL_MAX)
pub(crate) const MAXLOG: f64 = 709.782_712_893_384;
/// log(2**(-1075))
pub(crate) const MINLOG: f64 = -745.133_219_101_941_2;
/// Largest x such that Gamma(x) is finite
pub(crate) const MAXGAM: f64 = 171.624_376_956_302_7;
/// Euler-Mascheroni constant
pub(crate) const SCIPY_EULER: f64 = 0.577_215_664_901_532_9;
/// sqrt(1/pi)
pub(crate) const SQRT1OPI: f64 = 0.564_189_583_547_756_3;
/// sqrt(2/pi)
pub(crate) const SQRT2OPI: f64 = 0.797_884_560_802_865_4;
/// sqrt(3)
pub(crate) const SQRT3: f64 = 1.732_050_807_568_877_2;

#[inline]
pub(crate) fn vec_into<S: Into<T>, T>(xs: Vec<S>) -> Vec<T> {
    xs.into_iter().map(S::into).collect()
}

#[inline]
pub(crate) fn vec_to_vecvec<T: Clone>(
    vec: &[T],
    rows: usize,
    cols: usize,
    transpose: bool,
) -> Vec<Vec<T>> {
    if transpose {
        // Transpose: convert from (rows, cols) to (cols, rows)
        (0..cols)
            .map(|j| (0..rows).map(|i| vec[i * cols + j].clone()).collect())
            .collect()
    } else {
        // Normal: convert from flat vec to (rows, cols)
        (0..rows)
            .map(|i| vec[i * cols..(i + 1) * cols].to_vec())
            .collect()
    }
}

/// Keep the first `rows` rows and `cols` columns
#[inline]
pub(crate) fn vecvec_truncate<T>(mut vv: Vec<Vec<T>>, rows: usize, cols: usize) -> Vec<Vec<T>> {
    vv.truncate(rows);
    for v in &mut vv {
        v.truncate(cols);
    }
    vv
}

#[inline]
pub(crate) fn vec_into_vecvec<S: Into<T>, T: Clone>(
    vec: Vec<S>,
    rows: usize,
    cols: usize,
    transpose: bool,
) -> Vec<Vec<T>> {
    vec_to_vecvec(&vec_into(vec), rows, cols, transpose)
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
        assert_eq!(LN_SQRT_2PI.to_bits(), (0.5 * LN_TAU).to_bits() + 1);
        assert_eq!(SQRT_PI.to_bits(), pi.sqrt().to_bits());
        assert_eq!(SQRT_TAU.to_bits(), tau.sqrt().to_bits());
        assert_eq!(SQRT2PI.to_bits(), SQRT_TAU.to_bits() + 1);
        assert_eq!(SQRT3.to_bits(), black_box(3.0_f64).sqrt().to_bits());
        assert_eq!(PISQ_6.to_bits(), 1.644_934_066_848_226_4_f64.to_bits());
        assert_eq!(SQRT_EPS.to_bits(), black_box(f64::EPSILON).sqrt().to_bits());
        assert_eq!(MACHEP, f64::EPSILON / 2.0);
        assert_eq!(MAXLOG.to_bits(), black_box(f64::MAX).ln().to_bits());
        assert_eq!(SQRT1OPI.to_bits(), (1.0 / pi.sqrt()).to_bits());
        assert_eq!(SQRT2OPI.to_bits(), (2.0 / pi).sqrt().to_bits());
        assert_eq!(FRAC_2_PI.to_bits(), (2.0 / pi).to_bits());
        assert_eq!(FRAC_2_SQRT_PI.to_bits(), (2.0 / pi.sqrt()).to_bits());
        assert_eq!(SQRT_2.to_bits(), black_box(2.0_f64).sqrt().to_bits());
        let (ln10, ln2) = (black_box(10.0_f64).ln(), black_box(2.0_f64).ln());
        assert_eq!((LN_10 / LN_2).to_bits(), (ln10 / ln2).to_bits());
    }
}
