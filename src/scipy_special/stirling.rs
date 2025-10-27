use num_traits::Unsigned;

/// Stirling number of the second kind *S(n,k)*
///
/// Stirling numbers of the second kind count the number of ways to partition a set with *n*
/// elements into *k* non-empty subsets.
///
/// This is a pure rust translation of [`scipy.special.stirling2(n, k, exact=True)`][stirling2].
///
/// [stirling2]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.stirling2.html
///
/// # Examples
///
/// ```
/// use xsf::stirling2;
///
/// assert_eq!(stirling2::<u64>(10, 3), 9_330);
/// assert_eq!(stirling2::<u64>(10, 0), 0);
/// assert_eq!(stirling2::<u64>(9, 3), 3_025);
///
/// // won't fit in u64
/// assert_eq!(stirling2::<u128>(42, 4), 805_932_309_912_663_709_372_025);
/// ```
///
/// # Types
/// - `T` (required): unsigned return type
///
/// # See also
/// - [`comb`](crate::comb): *k*-combinations of *n* things, <sub>*n*</sub>C<sub>*k*</sub>
/// - [`perm`](crate::perm): *k*-permutations of *n* things, <sub>*n*</sub>P<sub>*k*</sub>
#[inline]
pub fn stirling2<T: num_traits::PrimInt + Unsigned>(n: u32, k: u32) -> T {
    if n == k || (n > 0 && k == 1) {
        return T::one();
    } else if n < k || k == 0 {
        return T::zero();
    } else if k == 2 {
        // S(n, 2) = 2^(n-1) - 1
        return T::from((1 << (n - 1)) - 1).unwrap();
    } else if n == k + 1 {
        // S(n, n-1) = C(n, 2) = n(n-1)/2
        return T::from(n as u128 * (n as u128 - 1) / 2).unwrap();
    }

    let k = k as usize;
    let mut n_row = vec![T::one(); k];
    for i in 2..n as usize {
        for j in (1..i.min(k)).rev() {
            n_row[j] = T::from(j + 1).unwrap() * n_row[j] + n_row[j - 1];
        }
    }
    n_row[k - 1]
}

#[cfg(test)]
mod tests {
    const STIRLING2_TABLE: [[u64; 11]; 11] = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 3, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 7, 6, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 15, 25, 10, 1, 0, 0, 0, 0, 0],
        [0, 1, 31, 90, 65, 15, 1, 0, 0, 0, 0],
        [0, 1, 63, 301, 350, 140, 21, 1, 0, 0, 0],
        [0, 1, 127, 966, 1_701, 1_050, 266, 28, 1, 0, 0],
        [0, 1, 255, 3_025, 7_770, 6_951, 2_646, 462, 36, 1, 0],
        [0, 1, 511, 9_330, 34_105, 42_525, 22_827, 5_880, 750, 45, 1],
    ];

    #[test]
    fn test_stirling2() {
        for (n, row) in STIRLING2_TABLE.iter().enumerate() {
            for (k, &s0) in row.iter().enumerate() {
                let s: u64 = crate::stirling2(n as u32, k as u32);
                assert_eq!(s, s0, "S({}, {})", n, k);
            }
        }
    }
}
