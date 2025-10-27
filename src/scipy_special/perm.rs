use num_traits::{ConstOne, NumOps};

/// *k*-permutations of *n* things, <sub>*n*</sub>P<sub>*k*</sub>
///
/// Permutations of *n* things taken *k* at a time, also known as "partial permutations".
/// It can be written in terms of the falling factorial, but unlike [`pow_falling`](pow_falling),
/// this function computes the result exactly for integer types.
///
/// Corresponds to [`scipy.special.perm(n, k, exact=True)`][perm].
///
/// [perm]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.perm.html
///
/// # See also
/// - [`comb`](crate::comb): *k*-combinations of *n*, <sub>*n*</sub>C<sub>*k*</sub>
/// - [`factorial`](crate::factorial): factorial *n*!, the special case *k*=*n*
/// - [`pow_falling`](crate::pow_falling): falling factorial power function
#[inline]
pub fn perm<N: NumOps + ConstOne + Copy>(n: N, k: u32) -> N {
    (0..k).fold((N::ONE, n), |(p, j), _| (p * j, j - N::ONE)).0
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_perm_pos() {
        assert_eq!(crate::perm(10, 0), 1);
        assert_eq!(crate::perm(10, 1), 10);
        assert_eq!(crate::perm(10, 2), 90);
        assert_eq!(crate::perm(10, 3), 720);
        assert_eq!(crate::perm(10, 4), 5_040);
        assert_eq!(crate::perm(10, 5), 30_240);
    }

    #[test]
    fn test_perm_neg() {
        assert_eq!(crate::perm(-10, 0), 1);
        assert_eq!(crate::perm(-10, 1), -10);
        assert_eq!(crate::perm(-10, 2), 110);
        assert_eq!(crate::perm(-10, 3), -1_320);
        assert_eq!(crate::perm(-10, 4), 17_160);
        assert_eq!(crate::perm(-10, 5), -240_240);
    }
}
