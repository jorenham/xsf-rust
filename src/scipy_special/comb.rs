use num_traits::{PrimInt, Unsigned};

/// *k*-combinations of *n* things, <sub>*n*</sub>C<sub>*k*</sub>
///
/// This is often expressed as "*n* choose *k*".
///
/// Pure rust implementation of [`scipy.special.comb(n, k, exact=True)`][comb].
///
/// [comb]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.comb.html
///
/// # Types
/// - `N`: Any unsigned integer
///
/// # See also
/// - [`comb_rep`]: choosing with replacement
/// - [`perm`](crate::perm): *k*-permutations of *n*, <sub>*n*</sub>P<sub>*k*</sub>
/// - [`binom`](crate::binom): the binomial coefficient as a floating point function
#[inline]
pub fn comb<N: PrimInt + Unsigned>(n: N, k: N) -> N {
    if n < k {
        N::zero()
    } else {
        let m = n + N::one();

        // manual unrolling because fold and range require the `std::iter::Step` nightly trait
        let (mut val, mut j) = (N::one(), N::one());
        while j <= k.min(n - k) {
            val = val * (m - j) / j;
            j = j + N::one();
        }
        val
    }
}

/// *k*-combinations of *n* things with replacement
///
/// Because this (also) counts the number of multisets of size *k* drawn from an *n*-element set, it
/// is also known as the *multiset coefficient*, and sometimes expressed as "*n* multichoose *k*".
///
/// It is equivalent to [`comb(n + k - 1, k)`](comb).
///
/// Pure rust implementation of [`scipy.special.comb(n, k, exact=True, repetition=True)`][comb].
///
/// [comb]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.comb.html
///
/// # Types
/// - `N`: Any unsigned integer
///
/// # See also
/// - [`comb`]: *k*-combinations of *n* without replacement
/// - [`perm`](crate::perm): *k*-permutations of *n*, <sub>*n*</sub>P<sub>*k*</sub>
/// - [`binom`](crate::binom): the binomial coefficient as a floating point function
#[inline]
#[doc(alias = "multiset")]
#[doc(alias = "multichoose")]
pub fn comb_rep<N: PrimInt + Unsigned>(n: N, k: N) -> N {
    if n == N::zero() {
        if k == n { N::one() } else { N::zero() }
    } else {
        let m = n + k;

        // manual unrolling because fold and range require the `std::iter::Step` nightly trait
        let (mut val, mut j) = (N::one(), N::one());
        while j <= k.min(n - N::one()) {
            val = val * (m - j) / j;
            j = j + N::one();
        }
        val
    }
}
