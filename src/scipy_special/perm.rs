use num_traits::{ConstOne, NumOps};

/// *k*-permutations of *n* things, <sub>*n*</sub>P<sub>*k*</sub>
///
/// Permutations of *n* things taken *k* at a time, also known as "partial permutations".
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
/// The number of *k*-permutations of *n*, i.e. *n(n-1)(n-2)...(n-k+1)*
///
/// # See also
/// - [`comb`](crate::comb): *k*-combinations of *n* things, <sub>*n*</sub>C<sub>*k*</sub>
/// - [`factorial`](crate::factorial): factorial *n*!, the special case *k*=*n*
/// - [`pow_falling`](crate::pow_falling): falling factorial power function
#[inline]
pub fn perm<N: NumOps + ConstOne + Copy>(n: N, k: u32) -> N {
    (0..k).fold((N::ONE, n), |(p, j), _| (p * j, j - N::ONE)).0
}
