//! The machine constants of `cephes/const.h`

/// 2**-53
pub(crate) const MACHEP: f64 = 1.110_223_024_625_156_5e-16;
/// log(DBL_MAX)
pub(crate) const MAXLOG: f64 = 709.782_712_893_384;
/// log(2**(-1075))
pub(crate) const MINLOG: f64 = -745.133_219_101_941_2;
/// Largest x such that Gamma(x) is finite
pub(crate) const MAXGAM: f64 = 171.624_376_956_302_7;

#[cfg(test)]
mod tests {
    use super::*;
    use std::hint::black_box;

    #[test]
    fn test_consts_bits() {
        assert_eq!(MACHEP, f64::EPSILON / 2.0);
        assert_eq!(MAXLOG.to_bits(), black_box(f64::MAX).ln().to_bits());
    }
}
