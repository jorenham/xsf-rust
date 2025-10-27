use num_traits::{PrimInt, Unsigned};

/// Factorial *n*!
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
/// # Examples
///
/// ```
/// assert_eq!([0_u32, 1, 2, 4, 8].map(xsf::factorial), [1, 1, 2, 24, 40_320]);
/// ```
///
/// # See also
/// - [`factorial_checked`]: factorial with overflow checking
/// - [`multifactorial`]: multifactorial *n*!<sub>(*k*)</sub>
#[inline]
pub fn factorial<N: PrimInt + Unsigned>(n: N) -> N {
    multifactorial(n, N::one())
}

/// Factorial *n*! with overflow checking
///
/// # Examples
///
/// ```
/// assert_eq!(xsf::factorial_checked(12_u32), Some(479_001_600));
/// assert_eq!(xsf::factorial_checked(13_u32), None);
/// ```
///
/// # See also
/// - [`factorial`]: unchecked factorial
/// - [`multifactorial_checked`]: multifactorial with overflow checking
#[inline]
pub fn factorial_checked<N: PrimInt + Unsigned>(n: N) -> Option<N> {
    multifactorial_checked(n, N::one())
}

/// Multifactorial *n*!<sub>(*k*)</sub> for positive *k*
///
/// Corresponds to [`scipy.special.factorialk(n, k, exact=True)`][factorialk].
///
/// [factorialk]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.factorialk.html
///
/// # Panics
/// - If *k* = 0 and *n* = 1, as 1<sup>∞</sup> is indeterminate.
/// - If *k* = 0 and *n* > 1, as n<sup>∞</sup> is divergent.
///
/// # See also
/// - [`multifactorial_checked`]: multifactorial with overflow checking
/// - [`factorial`]: factorial *n*!, the special case *k* = 1
#[doc(alias = "factorialk")]
#[doc(alias = "factorial2")]
#[inline]
pub fn multifactorial<N: PrimInt + Unsigned>(n: N, k: N) -> N {
    if k.is_zero() {
        if n.is_zero() {
            n
        } else if n.is_one() {
            panic!("1^∞ is indeterminate");
        } else {
            panic!("n^∞ is divergent for n>1");
        }
    } else if n <= k {
        if n.is_zero() { N::one() } else { n }
    } else {
        let (mut out, mut val) = (n, n);
        while val > k {
            val = val - k;
            out = out * val;
        }
        out
    }
}

/// Multifactorial *n*!<sub>(*k*)</sub> for positive *k* with overflow checking
///
/// Corresponds to [`scipy.special.factorialk(n, k, exact=True)`][factorialk].
///
/// [factorialk]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.factorialk.html
///
/// # See also
/// - [`multifactorial`]: unchecked multifactorial
/// - [`factorial_checked`]: factorial *n*! with overflow checking
#[doc(alias = "factorialk_checked")]
#[doc(alias = "factorial2_checked")]
#[inline]
pub fn multifactorial_checked<N: PrimInt + Unsigned>(n: N, k: N) -> Option<N> {
    if k.is_zero() {
        if n.is_zero() { Some(n) } else { None }
    } else if n <= k {
        Some(if n.is_zero() { N::one() } else { n })
    } else {
        let (mut out, mut val) = (n, n);
        while val > k {
            val = val - k;
            out = out.checked_mul(&val)?;
        }
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    const FACT1_5: u8 = 120;
    const FACT1_8: u16 = 40_320;
    const FACT1_12: u32 = 479_001_600;
    const FACT1_20: u64 = 2_432_902_008_176_640_000;
    const FACT1_34: u128 = 295_232_799_039_604_140_847_618_609_643_520_000_000;

    #[test]
    fn test_factorial_one() {
        for n in [0, 1] {
            assert_eq!(crate::factorial(n as u8), 1);
            assert_eq!(crate::factorial(n as u16), 1);
            assert_eq!(crate::factorial(n as u32), 1);
            assert_eq!(crate::factorial(n as u64), 1);
            assert_eq!(crate::factorial(n as u128), 1);
        }
    }

    #[test]
    fn test_factorial_max() {
        // anything larger will overflow
        assert_eq!(crate::factorial(5_u8), FACT1_5);
        assert_eq!(crate::factorial(8_u16), FACT1_8);
        assert_eq!(crate::factorial(12_u32), FACT1_12);
        assert_eq!(crate::factorial(20_u64), FACT1_20);
        assert_eq!(crate::factorial(34_u128), FACT1_34);
    }

    #[test]
    fn test_factorial_checked_max() {
        assert_eq!(crate::factorial_checked(5_u8), Some(FACT1_5));
        assert_eq!(crate::factorial_checked(8_u16), Some(FACT1_8));
        assert_eq!(crate::factorial_checked(12_u32), Some(FACT1_12));
        assert_eq!(crate::factorial_checked(20_u64), Some(FACT1_20));
        assert_eq!(crate::factorial_checked(34_u128), Some(FACT1_34));
    }

    #[test]
    fn test_factorial_checked_overflow() {
        assert_eq!(crate::factorial_checked(6_u8), None);
        assert_eq!(crate::factorial_checked(9_u16), None);
        assert_eq!(crate::factorial_checked(13_u32), None);
        assert_eq!(crate::factorial_checked(21_u64), None);
        assert_eq!(crate::factorial_checked(35_u128), None);
    }

    const FACT2_7: u8 = 105;
    const FACT2_12: u16 = 46_080;
    const FACT2_20: u32 = 3_715_891_200;
    const FACT2_33: u64 = 6_332_659_870_762_850_625;
    const FACT2_56: u128 = 81_842_841_814_930_553_085_241_614_925_824_000_000;

    #[test]
    fn test_multifactorial_k2_one() {
        for n in 0..=2 {
            let y = if n == 0 { 1 } else { n };
            assert_eq!(crate::multifactorial(n as u8, 2), y as u8, "{n}");
            assert_eq!(crate::multifactorial(n as u16, 2), y as u16, "{n}");
            assert_eq!(crate::multifactorial(n as u32, 2), y as u32, "{n}");
            assert_eq!(crate::multifactorial(n as u64, 2), y as u64, "{n}");
            assert_eq!(crate::multifactorial(n as u128, 2), y as u128, "{n}");
        }
    }

    #[test]
    fn test_multifactorial_k2_max() {
        assert_eq!(crate::multifactorial(7_u8, 2), FACT2_7);
        assert_eq!(crate::multifactorial(12_u16, 2), FACT2_12);
        assert_eq!(crate::multifactorial(20_u32, 2), FACT2_20);
        assert_eq!(crate::multifactorial(33_u64, 2), FACT2_33);
        assert_eq!(crate::multifactorial(56_u128, 2), FACT2_56);
    }

    #[test]
    fn test_multifactorial_k() {
        const N: [u64; 9] = [0, 1, 2, 3, 4, 5, 6, 7, 8];
        assert_eq!(
            N.map(|n| crate::multifactorial(n, 1)),
            [1, 1, 2, 6, 24, 120, 720, 5040, 40320],
        );
        assert_eq!(
            N.map(|n| crate::multifactorial(n, 2)),
            [1, 1, 2, 3, 8, 15, 48, 105, 384],
        );
        assert_eq!(
            N.map(|n| crate::multifactorial(n, 3)),
            [1, 1, 2, 3, 4, 10, 18, 28, 80],
        );
        assert_eq!(
            N.map(|n| crate::multifactorial(n, 4)),
            [1, 1, 2, 3, 4, 5, 12, 21, 32],
        );
    }
}
