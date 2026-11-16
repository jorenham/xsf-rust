use core::iter;

use num_traits::{ConstOne, NumOps, Zero};

/// $k$-permutations of $n$ things, $_nP_k$
///
/// Permutations of $n$ things taken $k$ at a time, also known as "partial permutations".
/// It can be written in terms of the falling factorial, but unlike
/// [`xsf::pow_falling`](crate::pow_falling), this function computes the result exactly for integer
/// types.
///
/// Corresponds to [`scipy.special.perm(n, k, exact=True)`][perm].
///
/// [perm]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.perm.html
///
/// # Types
/// - `N`: Any numeric type supporting basic arithmetic operations
///
/// # Returns
/// The number of $k$-permutations of $n$, i.e. $n(n-1)(n-2) \cdots (n-k+1)$, which is 0 if $n$ is a
/// non-negative integer below $k$ (unless a partial product overflows)
///
/// # See also
/// - [`comb`](crate::comb): $k$-combinations of $n$ things, $_nC_k$
/// - [`factorial`](crate::factorial): factorial $n!$, the special case $k = n$
/// - [`pow_falling`](crate::pow_falling): falling factorial power function
#[inline]
pub fn perm<N: NumOps + ConstOne + Zero + Copy>(n: N, k: u32) -> N {
    // n, n - 1, ..., up to a zero factor, past which `j - 1` could underflow
    let factors = iter::successors(Some(n), |&j| (!j.is_zero()).then(|| j - N::ONE));
    // (`0..k` first, so that `zip` stops before computing another factor)
    (0..k).zip(factors).fold(N::ONE, |p, (_, j)| p * j)
}
