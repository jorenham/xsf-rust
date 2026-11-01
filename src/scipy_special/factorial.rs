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
#[doc(alias = "factorialk", alias = "factorial2")]
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
#[doc(alias = "factorialk_checked", alias = "factorial2_checked")]
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
