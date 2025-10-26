use num_traits::{PrimInt, Unsigned};

/// Compute *n*!
///
/// Pure rust implementation of [`scipy.special.factorial(n, exact=True)`][factorial].
///
/// [factorial]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.factorial.html
///
/// # Maximum `n`
///
/// Anything larger than the maximum `n` for the given unsigned integer type will overflow.
///
/// | type     | `n <= _` |
/// |----------|---------:|
/// | [`u8`]   |        5 |
/// | [`u16`]  |        8 |
/// | [`u32`]  |       12 |
/// | [`u64`]  |       20 |
/// | [`u128`] |       34 |
///
/// # Types
/// - `N`: Any unsigned integer
///
/// # Examples
///
/// ```
/// assert_eq!([0_u32, 1, 2, 4, 8].map(xsf::factorial), [1, 1, 2, 24, 40_320]);
/// ```
///
/// # See also
/// - [`factorial_checked`] for the checked version
#[inline]
pub fn factorial<N: PrimInt + Unsigned>(n: N) -> N {
    let (mut val, mut i) = (N::one(), N::one());
    while i <= n {
        val = val * i;
        i = i + N::one();
    }
    val
}

/// Compute *n*! or return `None` on overflow
///
/// # Types
/// - `N`: Any unsigned integer
///
/// # Examples
///
/// ```
/// assert_eq!(xsf::factorial_checked(12_u32), Some(479_001_600));
/// assert_eq!(xsf::factorial_checked(13_u32), None);
/// ```
///
/// # See also
/// - [`factorial`] for the unchecked version
#[inline]
pub fn factorial_checked<N: PrimInt + Unsigned>(n: N) -> Option<N> {
    let (mut val, mut i) = (N::one(), N::one());
    while i <= n {
        val = val.checked_mul(&i)?;
        i = i + N::one();
    }
    Some(val)
}

#[cfg(test)]
mod tests {
    const FACT_5: u8 = 120;
    const FACT_8: u16 = 40_320;
    const FACT_12: u32 = 479_001_600;
    const FACT_20: u64 = 2_432_902_008_176_640_000;
    const FACT_34: u128 = 295_232_799_039_604_140_847_618_609_643_520_000_000;

    #[test]
    fn test_factorial_min() {
        assert_eq!(crate::factorial(0_u8), 1_u8);
        assert_eq!(crate::factorial(0_u16), 1_u16);
        assert_eq!(crate::factorial(0_u32), 1_u32);
        assert_eq!(crate::factorial(0_u64), 1_u64);
        assert_eq!(crate::factorial(0_u128), 1_u128);
    }

    #[test]
    fn test_factorial_max() {
        // anything larger will overflow
        assert_eq!(crate::factorial(5_u8), FACT_5);
        assert_eq!(crate::factorial(8_u16), FACT_8);
        assert_eq!(crate::factorial(12_u32), FACT_12);
        assert_eq!(crate::factorial(20_u64), FACT_20);
        assert_eq!(crate::factorial(34_u128), FACT_34);
    }

    #[test]
    fn test_factorial_checked_max() {
        assert_eq!(crate::factorial_checked(5_u8), Some(FACT_5));
        assert_eq!(crate::factorial_checked(8_u16), Some(FACT_8));
        assert_eq!(crate::factorial_checked(12_u32), Some(FACT_12));
        assert_eq!(crate::factorial_checked(20_u64), Some(FACT_20));
        assert_eq!(crate::factorial_checked(34_u128), Some(FACT_34));
    }

    #[test]
    fn test_factorial_checked_overflow() {
        assert_eq!(crate::factorial_checked(6_u8), None);
        assert_eq!(crate::factorial_checked(9_u16), None);
        assert_eq!(crate::factorial_checked(13_u32), None);
        assert_eq!(crate::factorial_checked(21_u64), None);
        assert_eq!(crate::factorial_checked(35_u128), None);
    }
}
