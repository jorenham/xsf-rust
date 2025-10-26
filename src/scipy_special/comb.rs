use num_traits::{PrimInt, Unsigned};

/// The number of combinations of *n* things taken *k* at a time
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
/// - [`comb_rep`](crate::comb_rep): choosing with replacement
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

/// The number of combinations of *n* things taken *k* at a time with replacement
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
/// - [`comb`](comb): choosing without replacement
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

#[cfg(test)]
mod tests {
    const PASCAL_L8: [[usize; 8]; 8] = [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [1, 2, 1, 0, 0, 0, 0, 0],
        [1, 3, 3, 1, 0, 0, 0, 0],
        [1, 4, 6, 4, 1, 0, 0, 0],
        [1, 5, 10, 10, 5, 1, 0, 0],
        [1, 6, 15, 20, 15, 6, 1, 0],
        [1, 7, 21, 35, 35, 21, 7, 1],
    ];
    const PASCAL_S08: [[usize; 8]; 9] = [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1],
        [1, 2, 3, 4, 5, 6, 7, 8],
        [1, 3, 6, 10, 15, 21, 28, 36],
        [1, 4, 10, 20, 35, 56, 84, 120],
        [1, 5, 15, 35, 70, 126, 210, 330],
        [1, 6, 21, 56, 126, 252, 462, 792],
        [1, 7, 28, 84, 210, 462, 924, 1716],
        [1, 8, 36, 120, 330, 792, 1716, 3432],
    ];

    #[test]
    fn test_comb() {
        for (n, nc) in PASCAL_L8.iter().enumerate() {
            for (k, &nck) in nc.iter().enumerate() {
                assert_eq!(crate::comb(n, k), nck);
            }
        }
    }

    #[test]
    fn test_comb_rep() {
        for (n, nr) in PASCAL_S08.iter().enumerate() {
            for (k, &nrk) in nr.iter().enumerate() {
                assert_eq!(crate::comb_rep(n, k), nrk);
            }
        }
    }
}
